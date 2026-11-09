// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 NeoFOAM authors

// Reads the same case files as OpenFOAM's icoUncoupledKinematicParcelFoam:
//   constant/g                         gravitational acceleration
//   constant/transportProperties       rhoInf, nu
//   constant/<cloud>Properties         cloud settings (supported subset: see the solver header)
//   constant/<positionsFile>           injection positions
//   0/U                                frozen carrier velocity

// ─────────────────────────────────────────────────────────────────────
// Fluid fields
// ─────────────────────────────────────────────────────────────────────

Info << "\nReading g" << endl;
const Foam::meshObjects::gravity& g = Foam::meshObjects::gravity::New(runTime);

Info << "Reading transportProperties\n" << endl;

Foam::IOdictionary transportProperties(Foam::IOobject(
    "transportProperties",
    runTime.constant(),
    mesh,
    Foam::IOobject::MUST_READ,
    Foam::IOobject::NO_WRITE
));

// rhoInf: carrier density used by the cloud (as icoUncoupledKinematicParcelFoam)
const Foam::scalar rhoInf = transportProperties.get<Foam::scalar>("rhoInf");
const Foam::scalar nu = transportProperties.get<Foam::scalar>("nu");

Info << "Reading field U (frozen carrier velocity)\n" << endl;
Foam::volVectorField
    ofU(Foam::IOobject(
            "U",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE,
            Foam::IOobject::NO_REGISTER
        ),
        mesh);

// ─── Carrier properties for the particle solver ───
const NeoFOAM::lagrangian::CarrierProperties carrier {
    rhoInf,    // fluid density [kg/m^3]
    rhoInf* nu // dynamic viscosity mu = rho * nu
};

// ─────────────────────────────────────────────────────────────────────
// Particle setup
// ─────────────────────────────────────────────────────────────────────

const Foam::word cloudName = args.getOrDefault<Foam::word>("cloud", "kinematicCloud");

Info << "Reading " << cloudName << "Properties\n" << endl;

Foam::IOdictionary cloudProperties(Foam::IOobject(
    cloudName + "Properties",
    runTime.constant(),
    mesh,
    Foam::IOobject::MUST_READ,
    Foam::IOobject::NO_WRITE
));

const Foam::dictionary& solutionDict = cloudProperties.subDict("solution");
const Foam::dictionary& subModels = cloudProperties.subDict("subModels");

// ─── Particle properties ───
const Foam::scalar rho0 = cloudProperties.subDict("constantProperties").get<Foam::scalar>("rho0");

// ─── Locate method (Q1, Q2, Q10): solution/locate, chosen like the drag law ───
//   findCell (default)  OpenFOAM mesh.findCell on the host, cell value of U
//   walk                NeoFOAM locate + interpolate on the executor (GPU-capable):
//                       interpolationSchemes/U cell or cellPoint, periodic (cyclic) patches
const Foam::word locateMode = solutionDict.getOrDefault<Foam::word>("locate", "findCell");
const auto locateMethod = NeoFOAM::lagrangian::locateMethodFromName(locateMode);
const bool walkLocate = (locateMethod == NeoFOAM::lagrangian::LocateMethod::walk);

// ─── Carrier interpolation ───
const Foam::word interpolationU =
    solutionDict.subDict("interpolationSchemes").getOrDefault<Foam::word>("U", "cell");
if (!walkLocate && interpolationU != "cell")
{
    WarningInFunction << "interpolationSchemes/U " << interpolationU
                      << " needs solution/locate walk: using cell (cell value of U)" << endl;
}

// ─── Time integration ───
const auto scheme = NeoFOAM::lagrangian::integrationSchemeFromName(
    solutionDict.subDict("integrationSchemes").get<Foam::word>("U")
);
if (scheme == NeoFOAM::lagrangian::IntegrationScheme::analytical)
{
    FatalErrorInFunction
        << "integrationSchemes/U analytical is not implemented yet (use Euler or RK2)"
        << Foam::exit(Foam::FatalError);
}

// ─── Forces: subModels/particleForces ───
//   sphereDrag           Putnam (as OpenFOAM)
//   SchillerNaumannDrag  Schiller-Naumann (NeoFOAM only)
//   stokesDrag           Stokes (NeoFOAM only)
//   gravity              gravity + buoyancy, g from constant/g
NeoFOAM::lagrangian::DragModel particleDrag;

NeoFOAM::lagrangian::ForceModels particleForces;
particleForces.gravity.g = NeoFOAM::convert(g.value());

for (const Foam::word& force : subModels.subDict("particleForces").toc())
{
    if (force == "sphereDrag" || force == "SchillerNaumannDrag" || force == "stokesDrag")
    {
        if (particleForces.useDrag)
        {
            FatalErrorInFunction << "particleForces: give only one drag model"
                                 << Foam::exit(Foam::FatalError);
        }
        particleDrag.type = (force == "sphereDrag") ? NeoFOAM::lagrangian::DragType::sphere
                          : (force == "SchillerNaumannDrag")
                              ? NeoFOAM::lagrangian::DragType::SchillerNaumann
                              : NeoFOAM::lagrangian::DragType::Stokes;
        particleForces.useDrag = true;
    }
    else if (force == "gravity")
    {
        particleForces.useGravity = true;
    }
    else
    {
        FatalErrorInFunction << "particleForces: " << force
                             << " is not supported (valid: sphereDrag, "
                             << "SchillerNaumannDrag, stokesDrag, gravity)"
                             << Foam::exit(Foam::FatalError);
    }
}
if (!particleForces.useDrag)
{
    FatalErrorInFunction << "particleForces: a drag model is required "
                         << "(sphereDrag, SchillerNaumannDrag or stokesDrag)"
                         << Foam::exit(Foam::FatalError);
}

// ─── Injection: subModels/injectionModels, manualInjection at t = 0 ───
std::vector<NeoN::Vec3> x0;
std::vector<NeoN::Vec3> U0;
std::vector<NeoN::scalar> d0;
std::vector<NeoN::scalar> rho0Vec;
std::vector<NeoN::label> cell0;

for (const Foam::entry& e : subModels.subDict("injectionModels"))
{
    const Foam::dictionary& inj = e.dict();
    const Foam::dictionary& sizeDict = inj.subDict("sizeDistribution");
    if (inj.get<Foam::word>("type") != "manualInjection"
        || sizeDict.get<Foam::word>("type") != "fixedValue" || inj.get<Foam::scalar>("SOI") != 0)
    {
        FatalErrorInFunction << "injector " << e.keyword() << ": only manualInjection with a "
                             << "fixedValue sizeDistribution and SOI 0 is supported"
                             << Foam::exit(Foam::FatalError);
    }
    const Foam::scalar d = sizeDict.subDict("fixedValueDistribution").get<Foam::scalar>("value");
    const Foam::vector Ui = inj.get<Foam::vector>("U0");

    Foam::vectorIOField positions(Foam::IOobject(
        inj.get<Foam::word>("positionsFile"),
        runTime.constant(),
        mesh,
        Foam::IOobject::MUST_READ,
        Foam::IOobject::NO_WRITE
    ));
    for (const Foam::vector& p : positions)
    {
        // Q1 at injection: the cell containing the particle
        const Foam::label celli = mesh.findCell(p);
        if (celli < 0)
        {
            FatalErrorInFunction << "injector " << e.keyword() << ": position " << p
                                 << " is outside the mesh" << Foam::exit(Foam::FatalError);
        }
        x0.push_back(NeoFOAM::convert(p));
        U0.push_back(NeoFOAM::convert(Ui));
        d0.push_back(d);
        rho0Vec.push_back(rho0);
        cell0.push_back(static_cast<NeoN::label>(celli));
    }
}

// ─── Build the particle cloud (memory of the executor: Serial, CPU or GPU) ───
NeoFOAM::lagrangian::ParticleCloud cloud(rt.exec, x0, U0, d0, rho0Vec, cell0);

// ─── Create the particle solver (with locate walk it also builds the tracking mesh
//     and the carrier field on the executor) ───
NeoFOAM::lagrangian::ParticleSolver particleSolver(
    particleDrag,
    particleForces,
    scheme,
    carrier,
    locateMethod,
    walkLocate ? NeoFOAM::lagrangian::interpolationFromName(interpolationU)
               : NeoFOAM::lagrangian::Interpolation::cell,
    rt.exec,
    mesh,
    ofU
);

// ─── Fluid velocity seen by each particle (updated every time step) ───
NeoN::Vector<NeoN::Vec3> particleUf(rt.exec, cloud.size(), NeoN::Vec3(0.0, 0.0, 0.0));

// ─── Frozen air: the cell values of U, copied once from OpenFOAM into NeoN ───
// Kept on the executor so the particle solver never reads the OpenFOAM field itself;
// the coupled program passes the velocity of the current time step here instead.
const auto Ufrozen = NeoFOAM::fromFoamField(rt.exec, ofU.primitiveField());

Info << "Particle cloud " << cloudName << ": " << Foam::label(cloud.size()) << " particle(s), rho0 "
     << rho0 << ", " << particleForces.size() << " force(s), integration "
     << solutionDict.subDict("integrationSchemes").get<Foam::word>("U") << ", locate " << locateMode
     << ", interpolation " << (walkLocate ? interpolationU : Foam::word("cell")) << nl << endl;

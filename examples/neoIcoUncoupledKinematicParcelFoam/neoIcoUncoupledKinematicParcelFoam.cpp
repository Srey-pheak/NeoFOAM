// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 NeoFOAM authors

/* neoIcoUncoupledKinematicParcelFoam
 *
 * Lagrangian particles in a frozen carrier flow (one-way coupling), the NeoFOAM
 * counterpart of OpenFOAM's icoUncoupledKinematicParcelFoam. Runs on the same
 * case directory. The particle physics is in the NeoFOAM library
 * (include/NeoFOAM/lagrangian/).
 *
 * Supported subset of <cloud>Properties:
 *   solution/integrationSchemes/U     Euler | RK2
 *   solution/interpolationSchemes/U   cell | cellPoint (cellPoint needs locate walk)
 *   solution/locate                   findCell (default) | walk
 *   constantProperties/rho0
 *   subModels/particleForces          sphereDrag | SchillerNaumannDrag | stokesDrag (one required),
 *                                     gravity (optional)
 *   subModels/injectionModels/*       manualInjection with fixedValue size, SOI 0
 *
 * Q1, Q2, Q10 (where is the particle, what is the air doing there, did it leave):
 *   locate findCell (default): OpenFOAM's mesh.findCell and the cell value of U, on the
 *     host: the particle data is copied to the host and back every time step.
 *   locate walk: NeoFOAM's locate (walk from cell to cell, periodic patches) and
 *     interpolate (cell or cellPoint), in a kernel on the executor: no copies.
 * A particle outside the mesh is deactivated. One particle step per time step
 * (maxCo is not used).
 *
 * Q3 - Q8 (forces and time integration) run in one kernel on the executor
 * (Serial, CPU or GPU).
 *
 * Options: -executor Serial|CPU|GPU, -cloud <name> (default kinematicCloud)
 *
 * Output: <case>/postProcessing/lagrangian/<cloud>/particles.csv
 *         columns t,id,x,y,z,ux,uy,uz,d, written at t = 0 and at every write time
 */


#include "NeoN/NeoN.hpp"
#include "NeoFOAM/NeoFOAM.hpp"
#include "NeoFOAM/lagrangian/CarrierField.hpp"
#include "NeoFOAM/lagrangian/CarrierProperties.hpp"
#include "NeoFOAM/lagrangian/IntegrationScheme.hpp"
#include "NeoFOAM/lagrangian/LocateMethod.hpp"
#include "NeoFOAM/lagrangian/ParticleCloud.hpp"
#include "NeoFOAM/lagrangian/ParticleOutput.hpp"
#include "NeoFOAM/lagrangian/ParticleSolver.hpp"
#include "NeoFOAM/lagrangian/drag/DragModel.hpp"
#include "NeoFOAM/lagrangian/force/ForceModel.hpp"
#include "fvCFD.H"
#include "vectorIOField.H"
#include <fstream>


using Foam::Info;
using Foam::endl;
using Foam::nl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/// Main program:
int main(int argc, char* argv[])
{
    Foam::argList::addOption("executor", "word", "NeoN executor type (Serial/CPU/GPU)");
    Foam::argList::addOption("cloud", "name", "Cloud name (default: kinematicCloud)");

// addCheckCaseOptions.H inside the NeoN init/finalize bracket, see neoIcoFoam
#include "addCheckCaseOptions.H"
#include "setRootCase.H"
    NeoN::initialize(argc, argv);
    {
// createTime.H inside the NeoN init/finalize bracket, see neoIcoFoam
#include "createTime.H"
        auto rt = NeoFOAM::createAdapterRunTime(runTime, args);
        auto& mesh = rt.mesh;

#include "createFields.H"

        // output
        const Foam::fileName outDir =
            runTime.globalPath() / "postProcessing" / "lagrangian" / cloudName;
        Foam::mkDir(outDir);
        std::ofstream csv(outDir / "particles.csv");
        NeoFOAM::lagrangian::writeParticlesHeader(csv);
        NeoFOAM::lagrangian::writeParticles(csv, runTime.value(), cloud);

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        Info << "Starting time loop\n" << endl;
        NeoN::label left = 0;
        while (runTime.loop())
        {
            rt.t = runTime.value();
            rt.dt = runTime.deltaTValue();
            Info << "Time = " << runTime.timeName() << nl << endl;

            // Find the cell the particle in and check what the air is doing there
            // (locate method chosen in the case: findCell or walk) from case directory
            left += particleSolver.locateAndInterpolate(cloud, particleUf, Ufrozen);


            // After locateAndInterpolate(), the particleUf vector contains the fluid velocity at
            // each particle's position, Advance the particles by one time step using the chosen
            // integration scheme (Euler or RK2) and the computed forces (drag, gravity, etc.).
            particleSolver.evolve(cloud, particleUf, rt.dt);


            // Write particle data to CSV file at specified write times, and report the number of
            // particles that have left the domain.
            if (runTime.writeTime())
            {
                NeoFOAM::lagrangian::writeParticles(csv, rt.t, cloud);
                if (left > 0)
                {
                    Info << "    particles that left the domain: " << Foam::label(left) << endl;
                }
            }

            runTime.printExecutionTime(Info);
        }

        Info << "End\n" << endl;
    }
    NeoN::finalize();

    return 0;
}

// ************************************************************************* //

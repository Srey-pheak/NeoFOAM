# neoIcoUncoupledKinematicParcelFoam

One-way coupled Lagrangian particle tracking solver for solid particles
in a frozen incompressible fluid field.

# Workflow of lagrangain particle tracking

1. Where is the particle is - cell localization
2. What are the air in the position of the particle is doing ?  - Field Interpolation
3. How fast is the particle is slippig - the relative velocity ( Uf-up)
4. what is reynold number of the particle  - Particle Reynolds number
5. What is Drag Correlation of the particle - Drag correlation
6. what forces acting on the particle  - Force Balance
7. How face the particle accelerate - Newton's 2nd law
8. Base on accelerate, where the particle go next ? - Time integration
9. Did the particle cross a cell? - Face-to-face tracking
10. Did the particle leave the domain? - Boundary handling


## Physics

- **Fluid**: Frozen (pre-computed by a separate solver)
- **Particles**: Solid, inert, spherical
- **Coupling**: One-way (particles feel fluid, fluid does not feel particles)
- **Forces**: Drag (sphereDrag = Putnam as OpenFOAM, SchillerNaumannDrag, stokesDrag), gravity + buoyancy
- **Integration**: Euler, RK2 (selectable); analytical planned

## Status

| Feature | Status |
|:---|:---|
| One-way coupling | ✅ Done |
| Drag + gravity | ✅ Done |
| Euler / RK2 integrators | ✅ Done |
| Cell localization (OpenFOAM findCell), cell-value U | ✅ MVP |
| cellPoint interpolation, face-to-face tracking, periodic boundaries | ⬜ Planned |
| Wall interaction | Out of the plans |
| Dispersion | ⬜ Planned |
| Two-way coupling | ⬜ Planned | - in neoIcocoupledKinematicParcelFoam

## Build

```bash
cmake --preset production
cmake --build build/production --target neoIcoUncoupledKinematicParcelFoam
cmake --install build/production

# PLIC

Verification of the PLIC volume-of-fluid library (`PLIC`), as a small solver would use it: an
`AmrCore` that owns the mesh and the ghost cells, with a prescribed velocity, and advances the
volume fraction with `PLIC::ComputeEdgeStateAndFluxes` in its own `MFIter` loop (`AdvectTest::Advance`
in `main.cpp`). A run exits with a non-zero status when one of its checks does not hold.

```
make -j                  # DIM = 2 by default; DIM=3 for inputs.3d.deform
mpirun -np 4 ./main2d.gnu.MPI.ex inputs.kernels
mpirun -np 4 ./main2d.gnu.MPI.ex inputs.2d.vortex
```

| Input                   | Case                                                            |
|-------------------------|-----------------------------------------------------------------|
| `inputs.kernels`        | Plane geometry, normals, plane fit, wall ghost cells, face density |
| `inputs.2d.vortex`      | Circle in the time-reversed single vortex, 3 levels             |
| `inputs.2d.zalesak`     | Zalesak's slotted disk, one revolution, 3 levels                |
| `inputs.2d.wall_vortex` | Circle cut by a wall, in the same flow: a free contact line     |
| `inputs.3d.deform`      | Sphere in LeVeque's deformation field, 3 levels                 |

The advection cases return to their initial shape at `t = prob.period`, so the shape error is the L1
difference from the initial volume fraction. `plic.reconstruction` selects the reconstruction
(default `myc_height`). Set `amr.max_level = 0` and `amr.n_cell` for a uniform grid. The `test.max_*`
parameters set the checks.

`PLICTestUtils.H` holds what the test solver does around the library and a real solver would do its
own way: the prescribed velocities, the time step, the refinement tags and the liquid volume.

## Expected results

`inputs.kernels` checks, each to round-off, that

- `volume_intercept` and `cut_volume` invert each other, including grid-aligned planes;
- the plane of each reconstruction holds the volume of its cell;
- the `myc_height` normal is exact for a planar interface (`youngs` and `myc` are not, and their
  largest normal errors are printed);
- `FitPlanes` on a level of several boxes fits the mixed cells, valid and in the first ghost layer,
  and zeroes the pure ones;
- `FillWallGhosts` continues a planar interface through walls, and keeps the ghost cells beyond
  `ext_dir` (inflow) and `reflect_even` (symmetry) faces as they were filled;
- `SharpFaceDensity` balances gravity normal to a planar interface: the curl of `rho_face g`
  vanishes.

In every advection case the liquid volume is conserved to round-off per step (about 1e-15) and F
leaves [0,1] by round-off only. The results do not depend on the box layout or the number of ranks.

| Input                   | L1 shape error |
|-------------------------|----------------|
| `inputs.2d.vortex`      | 1.39e-4        |
| `inputs.2d.zalesak`     | 4.49e-4        |
| `inputs.2d.wall_vortex` | 6.17e-4        |
| `inputs.3d.deform`      | 5.48e-4        |

L1 shape error on uniform grids (single vortex with `amr.fixed_dt = 0.5/n`):

| n     | Zalesak  | single vortex |
|-------|----------|---------------|
| 64²   | 2.96e-3  | 1.87e-3       |
| 128²  | 1.03e-3  | 4.68e-4       |
| 256²  | 4.68e-4  | 1.39e-4       |

Zalesak's disk is limited to about first order by the unresolved slot corners, as expected for PLIC
schemes. With the 3 levels of `inputs.2d.vortex` the error is that of the finest level.

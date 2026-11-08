# PLIC

Verification of the PLIC volume-of-fluid library (`PLIC`), as a small solver would use it: an
`AmrCore` that owns the mesh and the ghost cells, with a prescribed velocity, and advances the
volume fraction with `PLIC::ComputeEdgeStateAndFluxes` in its own `MFIter` loop (`AdvectTest::Advance`
in `test_advect.cpp`). A run exits with a non-zero status when one of its checks does not hold.

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

## Sources

| File                | Checks                                                     |
|---------------------|------------------------------------------------------------|
| `test_geometry.cpp` | volume <-> plane relations, normals, plane fit, `FitPlanes` |
| `test_walls.cpp`    | `FillWallGhosts`                                           |
| `test_density.cpp`  | `SharpFaceDensity`                                         |
| `test_advect.cpp`   | the advection cases: `ComputeEdgeStateAndFluxes` on an `AmrCore` |
| `test_util.H/.cpp`  | shared helpers: planes on a level, the largest error over regions of a MultiFab |

## Expected results

`inputs.kernels` checks, each to round-off, that

- `volume_intercept` and `cut_volume` invert each other, including grid-aligned planes, and planes
  close to a diagonal of a cell face (one tiny coefficient, the other two nearly equal);
- the plane of each reconstruction holds the volume of its cell;
- the `myc_height` normal is exact for a planar interface, in every orientation (`youngs` and
  `myc` are not, and their largest normal errors are printed), and in 3D over random orientations,
  among them one at 44 degrees with a tilt, where a height column centered on the cell would miss
  the crossing of a neighboring column;
- the `myc_height` normal is exact for a planar interface on cells of aspect ratio 1:4:2;
- with pure cells overshooting [0, 1] by round-off, the `myc_height` normal stays exact for a
  planar interface, and by 1e-3 every reconstruction still gives planes that hold the volume;
- `myc_height` falls back to the `myc` normal when a height column holds a liquid sliver above a
  thin gas gap, so that it crosses the interface more than once;
- the `myc` and `myc_height` normals keep the symmetries of spheres (swapping x and y, x and z,
  mirroring in x): where centered candidates of MYC tie, MYC takes their mean, and where components
  of the MYC normal tie, `myc_height` averages the heights along each;
- on films (normal to each axis, and diagonal), gaps, filaments and drops one cell thick, symmetric
  across them or with noise along them, the plane of each reconstruction lies across the structure
  (along the principal direction of the Hessian of `F`) and holds the volume of the cell, and the
  cells of a film agree on its orientation;
- on spheres of radius 0.6 to 6 cells, no cell is taken for such a structure, and every
  reconstruction keeps its normal;
- `FitPlanes` on a level of several boxes fits the mixed cells, valid and in the first ghost layer,
  and zeroes the pure ones;
- `FillWallGhosts` continues planar interfaces through walls on every face, in every ghost cell of
  the `nGrowF` layers, edges and corners included, also when `F` has more ghost layers than the
  search reads (15 for `myc_height`), and keeps the ghost cells beyond
  `ext_dir` (inflow) and `reflect_even` (symmetry) faces as they were filled;
- on a fine level that covers part of the domain, `FillWallGhosts` continues a planar interface
  through walls near its coarse/fine boundaries, where the search reaches past the ghost cells of
  `F`;
- `FillWallGhosts` does not read the boundary condition of a periodic direction, which may be left
  unset;
- `FillWallGhosts` continues a plane through a wall near its contact line with each reconstruction
  where its plane is exact: `youngs` at 45 degrees, `myc` and `myc_height` at 45, 68, 11 and 3
  degrees (at 11, `myc_height` falls back to the MYC plane; at 3, the interface runs inside the wall
  layer and both take its one-sided height function);
- `myc_height` continues random planes through a wall near their contact line, half of them near
  27 degrees, where the search must reach past the MYC range;
- `FillWallGhosts` gives the same ghost values, exactly, with one box, boxes of 16 cells and boxes
  of 8 cells, for a curved interface meeting walls on every face, and periodic in x with the
  interface across the periodic boundary, in every ghost cell of every box (a ghost cell that
  several boxes hold is compared in each);
- `FillWallGhosts` keeps the symmetries of a curved interface meeting walls on every face (swapping x
  and y, x and z, mirroring in x), at edges and corners too, with each reconstruction: no face or
  direction along a wall is preferred; boxes that share a wall ghost cell give it the same value.
- `SharpFaceDensity` balances gravity normal to a planar interface: the curl of `rho_face g`
  vanishes (in every component over random planes), its ghost faces are left as they are, and over
  several F it is the mean of their densities;
- with the phase densities as fields, uniform fields give the density of their values to the bit; on
  faces with no liquid the density is the mean of the gas density field of the two cells, and balances
  gravity for a gas density linear along it; across a planar interface normal to gravity such a gas
  density is balanced to |grad rho_g| h / 8 on the cut faces.

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

# PLIC

Verification of the PLIC volume-of-fluid library (`PLIC`). A run exits with a non-zero status when
one of its checks does not hold.

```
make -j                  # DIM = 2 by default
mpirun -np 4 ./main2d.gnu.MPI.ex inputs.kernels
```

| Input            | Case                                           |
|------------------|------------------------------------------------|
| `inputs.kernels` | Plane geometry, normals, plane fit             |

## Sources

| File                | Checks                                                      |
|---------------------|-------------------------------------------------------------|
| `test_geometry.cpp` | volume <-> plane relations, normals, plane fit, `FitPlanes` |
| `test_util.H/.cpp`  | shared helpers: counting the checks, round-off tolerance    |

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
  (along the principal direction of the Hessian of `F`) and holds the volume of the cell, and on
  spheres of radius 0.6 to 6 cells the normals are kept (not taken for such structures);
- `FitPlanes` on a level of several boxes fits the mixed cells, valid and in the first ghost layer,
  and zeroes the pure ones.

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

## Expected results

`inputs.kernels` checks, each to round-off, that

- `volume_intercept` and `cut_volume` invert each other, including grid-aligned planes, and planes
  close to a diagonal of a cell face (one tiny coefficient, the other two nearly equal);
- the plane of each reconstruction holds the volume of its cell;
- the `myc_height` normal is exact for a planar interface, in every orientation (`youngs` and
  `myc` are not, and their largest normal errors are printed), and in 3D over random orientations,
  among them one at 44 degrees with a tilt, where a height column centered on the cell would miss
  the crossing of a neighboring column;
- `myc_height` falls back to the `myc` normal when a height column holds a liquid sliver above a
  thin gas gap, so that it crosses the interface more than once;
- the `myc_height` normal keeps the symmetries of a sphere (swapping x and y, x and z, mirroring
  in x): where components of the MYC normal tie, it averages the heights along each;
- `FitPlanes` on a level of several boxes fits the mixed cells, valid and in the first ghost layer,
  and zeroes the pure ones.

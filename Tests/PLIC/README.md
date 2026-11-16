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
- the `myc_height` normal is exact for a planar interface (`youngs` and `myc` are not, and their
  largest normal errors are printed);
- `myc_height` falls back to the `myc` normal when a height column holds a liquid sliver above a
  thin gas gap, so that it crosses the interface more than once;
- `FitPlanes` on a level of several boxes fits the mixed cells, valid and in the first ghost layer,
  and zeroes the pure ones.

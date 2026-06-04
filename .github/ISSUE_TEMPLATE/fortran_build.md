---
name: Fortran build problem
about: The shared library fails to compile or load
title: "[FORTRAN] "
labels: fortran, build
assignees: niki8885
---

## Which library

- [ ] `libsentinel_validation`  (`validation.f90`)
- [ ] `libsentinel_indices`     (`indices_mod.f90`)
- [ ] `libsentinel_raster_ops`  (`raster_ops.f90`)
- [ ] `libsentinel_processing`  (`pansharpening.f90`)
- [ ] `libsentinel_filters`     (`filters.f90`)

## Stage of failure

- [ ] Compile-time error (`gfortran` command fails)
- [ ] Library not found at runtime (`FileNotFoundError`)
- [ ] Library loads but a symbol is missing (`AttributeError`)
- [ ] Wrong results / crash during computation

## Compile command used

```bash
# paste the exact gfortran command here
```

## Error output

<details>
<summary>Full error</summary>

```
paste full compiler or Python traceback here
```

</details>

## Environment

| | |
|---|---|
| OS | |
| gfortran version | <!-- gfortran --version --> |
| Python | |
| sentinel-processor | |
| Installation method | <!-- MSYS2 UCRT64 / apt / homebrew / other --> |

## Additional context

<!-- Anything else — custom install prefix, conda environment, cross-compilation, etc. -->

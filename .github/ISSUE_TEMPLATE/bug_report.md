---
name: Bug report
about: Something is broken or behaving unexpectedly
title: "[BUG] "
labels: bug
assignees: niki8885
---

## Describe the bug

<!-- A clear, concise description of what is wrong. -->

## Steps to reproduce

```python
# Minimal code that reproduces the issue
import sentinel_processor as sp

```

## Expected behaviour

<!-- What did you expect to happen? -->

## Actual behaviour

<!-- What actually happened? Include the full traceback if there is one. -->

<details>
<summary>Traceback</summary>

```
paste traceback here
```

</details>

## Environment

| | |
|---|---|
| OS | <!-- e.g. Windows 11 / Ubuntu 24.04 / macOS 14 --> |
| Python | <!-- e.g. 3.12.3 --> |
| sentinel-processor | <!-- pip show sentinel-processor --> |
| rasterio | |
| rioxarray | |
| numpy | |

## Fortran libraries (if relevant)

- [ ] `libsentinel_validation` compiled and present
- [ ] `libsentinel_indices` compiled and present
- [ ] `libsentinel_raster_ops` compiled and present
- [ ] `libsentinel_processing` compiled and present
- [ ] `libsentinel_filters` compiled and present

```
gfortran --version
```

## Additional context

<!-- Anything else that might be relevant — input file format, bbox size, etc. -->

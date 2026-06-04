# Contributing to sentinel-processor

Thank you for taking the time to contribute. This document explains how to get set up, what conventions to follow, and how to submit changes.

---

## Table of contents

- [Getting started](#getting-started)
- [Project structure](#project-structure)
- [Making changes](#making-changes)
- [Tests](#tests)
- [Fortran conventions](#fortran-conventions)
- [Submitting a pull request](#submitting-a-pull-request)
- [Reporting issues](#reporting-issues)

---

## Getting started

```bash
git clone https://github.com/niki8885/sentinel-processor
cd sentinel-processor
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev,netcdf]"
```

Build the Fortran libraries following the commands in [README.md](../README.md#fortran-libraries). Tests that exercise Fortran paths are skipped automatically if a library is not found.

---

## Project structure

```
sentinel_processor/
├── filters/        ← Fortran convolution filters
├── indices/        ← Fortran spectral index kernels
├── input/          ← Downloader and STAC search
├── processing/     ← Fortran pansharpening + raster ops
├── utils/          ← Shared data classes (SpectralBands, LocationSpec, …)
├── validation/     ← Fortran SCL / radiometry checks
└── visualisation/  ← Plotly visualisation

tests/              ← pytest suite — one file per module
docs/               ← Markdown reference docs
.github/
├── workflows/ci.yml
├── ISSUE_TEMPLATE/
├── PULL_REQUEST_TEMPLATE.md
├── CONTRIBUTING.md
├── CODE_OF_CONDUCT.md
└── SECURITY.md
```

---

## Making changes

1. Fork the repository and create a branch from `master`:

   ```bash
   git checkout -b feat/my-feature
   ```

2. **One concern per PR.** A new filter, a bug fix, or a documentation update — not all three at once.

3. **Python style:**
   - Type annotations on all public functions and methods.
   - Docstrings on all public classes and functions.
   - No strict formatter enforced; follow the existing style.

4. **Fortran style:**
   - Use `iso_c_binding` and `bind(C, name="...")` on every exported subroutine.
   - Pass arrays as flat 1-D Fortran-order buffers — let the Python bridge handle transposition.
   - Add a `try / except FileNotFoundError` NumPy fallback in every bridge function so the package works without a compiled library.
   - Register `restype` and `argtypes` on every function in `_get_lib()`.

5. **Documentation:** update the relevant `docs/*.md` file for any public API change.

---

## Tests

```bash
# Full suite
pytest tests/ -v

# Only fast tests (no file I/O)
pytest tests/ -v -m "not integration"

# Single module
pytest tests/test_filters.py -v

# With coverage
pytest tests/ --cov=sentinel_processor --cov-report=term-missing
```

**Adding tests:**

- Every new public function needs at least one test in `tests/test_<module>.py`.
- New Fortran routines need a round-trip test that verifies the math against a Python reference value.
- Use `pytest.approx` for floating-point comparisons.
- Use `tmp_path` for any test that writes files — never write to the project tree.

---

## Fortran conventions

| Rule | Reason |
|---|---|
| `real(c_double)` / `integer(c_int)` throughout | Exact match with ctypes `c_double` / `c_int` |
| Scalar arguments as `intent(in), value` | Avoids pointer-to-scalar confusion at the ABI boundary |
| Arrays as flat 1-D buffers | Keeps ctypes argtypes declarations simple |
| Allocate temporaries inside the subroutine | Keeps the bridge stateless and thread-safe |
| Add `+ 1.0d-12` to divisors | Prevents division-by-zero without masking — matches existing convention |

---

## Submitting a pull request

1. Push your branch and open a PR against `master`.
2. Fill in the [pull request template](PULL_REQUEST_TEMPLATE.md).
3. CI must be green — all test and lint jobs — before merging.
4. A maintainer will review and may request changes.
5. Squash or rebase before merging to keep the history clean.

---

## Reporting issues

Use the [issue templates](ISSUE_TEMPLATE/) — there are templates for bug reports, feature requests, and Fortran build problems.

---

**Maintainer:** Nikita Manaenkov — [nick.maanenkov@gmail.com](mailto:nick.maanenkov@gmail.com)
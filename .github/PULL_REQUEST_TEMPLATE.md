## Summary

<!-- One or two sentences describing what this PR does and why. -->

## Type of change

- [ ] Bug fix
- [ ] New feature
- [ ] Fortran kernel (new or modified `.f90`)
- [ ] Documentation
- [ ] Refactor / performance
- [ ] Tests only

## Related issues

<!-- Closes #NNN -->

---

## Changes

<!-- Bullet list of what changed. Include file names for non-obvious changes. -->

-
-

## Fortran changes (if applicable)

- [ ] New subroutine added to `.f90` with `bind(C, name="...")`
- [ ] `restype` and `argtypes` registered in `_get_lib()`
- [ ] NumPy fallback added in the bridge function (`try / except FileNotFoundError`)
- [ ] Build command added to `README.md` or the relevant `docs/*.md`

---

## Tests

- [ ] New tests added for the changed/added code
- [ ] All existing tests pass locally (`pytest tests/ -v`)
- [ ] Coverage is not reduced

<!-- Paste the relevant pytest output section here if helpful -->

---

## Documentation

- [ ] Docstrings updated for changed public functions
- [ ] Relevant `docs/*.md` file updated
- [ ] `README.md` updated (if the public API changed)

---

## Checklist

- [ ] Branch is up to date with `master`
- [ ] CI is green
- [ ] No debug prints, commented-out code, or `.pyc` files committed
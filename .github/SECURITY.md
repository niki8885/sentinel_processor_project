# Security Policy

## Supported versions

| Version | Supported |
|---|---|
| latest (`master`) | ✅ |
| older releases | ❌ — please upgrade |

## Reporting a vulnerability

**Do not open a public GitHub issue for security vulnerabilities.**

Please report security issues by email to:

**Nikita Manaenkov** — [nick.maanenkov@gmail.com](mailto:nick.maanenkov@gmail.com)

Include as much detail as possible:

- A description of the vulnerability and its potential impact
- Steps to reproduce or a minimal proof-of-concept
- Affected version(s)
- Any suggested remediation

You will receive an acknowledgement within **72 hours**. After the issue is confirmed and a fix is available, a coordinated disclosure will be made via a GitHub Security Advisory.

## Scope

This project is a data-processing library. The main attack surfaces relevant to security reports are:

- **Arbitrary file writes** — the downloader writes raster files to paths derived from STAC metadata. A malicious STAC server could potentially craft filenames that escape the intended output directory.
- **Fortran shared libraries** — the `.dll` / `.so` files are loaded via ctypes. Shipping pre-built binaries from untrusted sources is outside the scope of this project; always build from source.
- **Dependency vulnerabilities** — if you find a vulnerability in a direct dependency (`rasterio`, `rioxarray`, `pystac-client`, etc.) please report it to that project directly. You may also notify us so we can pin or update the dependency.

## Out of scope

- Denial-of-service through large or malformed raster files (expected use case involves untrusted data sources)
- Issues that require physical access to the machine
---
name: Feature request
about: Suggest a new index, filter, algorithm, or API improvement
title: "[FEAT] "
labels: enhancement
assignees: niki8885
---

## Is this related to a problem?

<!-- A clear description of what the problem or limitation is.
     Example: "I need to compute NDRE but it is not in the indices registry." -->

## Proposed solution

<!-- Describe the feature you would like to see. Be as specific as possible.
     For new indices: include the formula and the required Sentinel-2 bands.
     For new filters: describe the algorithm and expected parameters.
     For API changes: show what the new call would look like. -->

```python
# Example of what the API would look like
from sentinel_processor.indices.compute import compute_indices

results = compute_indices(scene, ["ndre"])
```

## Alternatives considered

<!-- Any alternative approaches you have thought about. -->

## Additional context

<!-- References, papers, links to similar implementations, etc. -->

## Are you willing to contribute this feature?

- [ ] Yes, I can open a PR
- [ ] I can help review a PR
- [ ] No, just suggesting

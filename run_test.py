from sentinel_processor.indices.compute import compute_indices

results = compute_indices(
    source=r"downloads/spectral/budapest_target_20260520T094746.nc",
    indices=["ndvi", "evi", "ndwi", "ndbi", "nbr"],
)
print(results)
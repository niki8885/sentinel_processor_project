from sentinel_processor.visualisation.plot import plot_band, plot_rgb, plot_grid, plot_mask

plot_band("data/spectral/budapest_test_20260526T095725.nc", band="nir", save_html="data/vis/nir.html").show()
plot_band("data/indices/indices_budapest_test_20260526T095725_ndvi.tif", save_html="data/vis/ndvi.html").show()

plot_rgb("data/visual/vis_budapest_test_20260526T095725.nc").show()
plot_rgb("data/spectral/budapest_test_20260526T095725.nc", "red", "green", "blue").show()
plot_rgb("data/spectral/budapest_test_20260526T095725.nc", "nir", "red", "green").show()

plot_grid([
    {"file": "data/spectral/budapest_test_20260526T095725.nc", "band": "red",   "label": "Red"},
    {"file": "data/spectral/budapest_test_20260526T095725.nc", "band": "nir",   "label": "NIR"},
    {"file": "data/spectral/budapest_test_20260526T095725.nc", "band": "swir16","label": "SWIR1"},
    {"file": "data/indices/indices_budapest_test_20260526T095725_ndvi.tif",    "label": "NDVI", "colorscale": "RdYlGn"},
], ncols=4, save_html="data/vis/grid.html").show()

plot_mask("data/technical/scl_budapest_test_20260526T095725.nc").show()
plot_mask("data/technical/scl_budapest_test_20260526T095725.nc", bad_classes=[8, 9, 10]).show()
fig, arr = plot_mask("data/technical/scl_budapest_test_20260526T095725.nc", return_mask=True)
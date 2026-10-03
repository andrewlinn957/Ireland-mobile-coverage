# Ireland Mobile Signal Map

Source backup of the national Ireland mobile coverage map.

- Live Site: https://ireland-mobile-signal-map.andrewlinn.chatgpt.site
- Target Site version: 6 (CLC+ clutter scenario; version 5 remains the rollback point)
- Static map interface and current coverage data: `dist/`
- National coverage calculation: `scripts/build-national-coverage.py`
- Model inputs: `model-inputs/`

The map raster uses 250 m projected spacing (about 150 m on the ground at Irish latitudes). The model uses planning-derived estimates of the highest/top-mounted antenna height where a structure height is linked to a ComReg site within 30 m; elsewhere it retains the nominal 30 m transmitter-height default. The estimates, status labels, method and source references are in `model-inputs/mobile-site-height-estimates.csv` and its accompanying notes. The target clutter scenario uses Sentinel-derived Copernicus CLC+ Backbone 2021 classes converted to representative heights using an explicit assumed crosswalk; CORINE 2018 is retained for CLC+ no-data/outside-area fallback. This is not RF-validated. Proposed/ambiguous records are a planning-informed scenario, not confirmation of installed equipment, and the map remains indicative and uncalibrated.

Use `scripts/build-national-coverage.py` for the configured national build and `scripts/run-national-250m-height-informed.py` for an isolated research recalculation. Research outputs and run logs are kept outside `dist/` until reviewed.

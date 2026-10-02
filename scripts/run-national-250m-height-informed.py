#!/usr/bin/env python3
"""Full Republic 250 m planning-height scenario; never writes to dist/."""
from __future__ import annotations

import importlib.util
import json
import multiprocessing
import sys
import time
from itertools import islice
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "scripts/build-national-coverage.py"
OUT = ROOT / "research/height-informed-2026-10-02/data"
GRID_SPACING_M = 250.0


def load_model():
    spec = importlib.util.spec_from_file_location("coverage_model_shadow", MODEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main():
    started_all = time.monotonic()
    m = load_model()
    OUT.mkdir(parents=True, exist_ok=True)
    m.OUTPUT_DIR = OUT
    m.GRID_M = GRID_SPACING_M
    sites = json.loads(m.SITES_PATH.read_text())
    site_heights, height_metadata = m.load_site_height_estimates(sites)
    _, polygons, bounds = m.read_boundary()
    grid = m.grid_for_boundary(polygons, bounds)
    land_cells = len(grid["rows"])
    print(f"National 250 m grid: {grid['width']} x {grid['height']}; {land_cells:,} in-boundary cells", flush=True)
    started_inputs = time.monotonic()
    terrain = m.load_elevation_tiles(grid)
    clutter = m.ClutterSampler()
    climate = m.load_climate()
    indexes = m.build_site_band_indexes(sites, site_heights)
    input_seconds = time.monotonic() - started_inputs
    networks = list(sites["networks"])
    keys = [(network, tech) for tech in m.TECH_BITS for network in ["all", *networks]]
    signals = {key: np.full((grid["height"], grid["width"]), np.nan, np.float32) for key in keys}
    for tech in m.TECH_BITS:
        signals[("all", tech)][grid["rows"], grid["columns"]] = -np.inf
        for network in networks:
            signals[(network, tech)][grid["rows"], grid["columns"]] = -np.inf
    m._MODEL_WORKER_STATE = (terrain, clutter, climate, indexes, networks)
    cells = ((int(row), int(col), float(lon), float(lat))
             for row, col, lon, lat in zip(grid["rows"], grid["columns"], grid["lons"], grid["lats"]))
    output_keys = [(network, tech) for tech in m.TECH_BITS for network in ["all", *networks]]
    pool_size = min(8, multiprocessing.cpu_count())
    started_model = time.monotonic()
    completed = 0
    with ProcessPoolExecutor(max_workers=pool_size, mp_context=multiprocessing.get_context("fork")) as pool:
        while batch := list(islice(cells, 8192)):
            for result in pool.map(m.model_grid_cell, batch, chunksize=32):
                row, col, *values = result
                for key, value in zip(output_keys, values):
                    signals[key][row, col] = value
                completed += 1
                if completed % 25_000 == 0 or completed == land_cells:
                    elapsed = time.monotonic() - started_model
                    rate = completed / elapsed if elapsed else 0
                    remaining = (land_cells - completed) / rate if rate else 0
                    print(f"modelled {completed:,}/{land_cells:,} cells; {elapsed/60:.1f} min elapsed; {remaining/60:.1f} min estimated remaining", flush=True)
    model_seconds = time.monotonic() - started_model

    for (network, tech), signal in signals.items():
        provider = "all" if network == "all" else m.NETWORK_SLUGS[network]
        coverage, sensitivity = m.output_rgba(signal, grid["mask"])
        coverage.save(OUT / f"coverage-{provider}-{tech}.png", optimize=True)
        sensitivity.save(OUT / f"sensitivity-{provider}-{tech}.png", optimize=True)
        print(f"wrote {provider}-{tech} PNGs", flush=True)

    m.write_metadata(sites, grid, land_cells, height_metadata)
    metadata_path = OUT / "coverage-metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["gridSpacingProjectedM"] = GRID_SPACING_M
    metadata["approximateGroundResolutionM"] = 150
    metadata["shadowBuild"] = True
    metadata["productionOutputModified"] = False
    metadata["heightInformedScenario"] = True
    metadata_path.write_text(json.dumps(metadata, separators=(",", ":")) + "\n")
    total_seconds = time.monotonic() - started_all
    result = {
        "gridSpacingProjectedM": GRID_SPACING_M,
        "approximateGroundResolutionM": 150,
        "width": grid["width"], "height": grid["height"], "landCells": land_cells,
        "siteRecords": len(sites["records"]), "workers": pool_size,
        "heightEstimates": height_metadata,
        "inputPreparationSeconds": round(input_seconds, 2),
        "modelSeconds": round(model_seconds, 2), "totalSeconds": round(total_seconds, 2),
        "productionOutputModified": False,
        "outputDirectory": str(OUT),
    }
    (OUT.parent / "run-summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"completed: {json.dumps(result)}", flush=True)


if __name__ == "__main__":
    main()

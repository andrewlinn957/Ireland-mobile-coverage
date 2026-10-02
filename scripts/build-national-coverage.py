#!/usr/bin/env python3
"""Build the national raster layers for the Ireland mobile signal map.

Uses the ComReg catalogue and the same Δ-Bullington implementation and
link-budget assumptions as the trail atlas. Site-specific planning-derived
antenna-height estimates are used where available; otherwise the nominal 30 m
transmitter-height default applies. Coverage is calculated on a Web Mercator
grid whose spacing is configured by GRID_M.
"""

from __future__ import annotations

import json
import math
import multiprocessing
import csv
import sys
import time
from itertools import islice
import urllib.error
import urllib.request
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from io import BytesIO
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
SITES_PATH = ROOT / "dist/data/mobile-sites.json"
SITE_HEIGHTS_PATH = ROOT / "model-inputs/mobile-site-height-estimates.csv"
BOUNDARY_PATH = ROOT / "model-inputs/ireland-boundary.geojson"
CLIMATE_PATH = ROOT / "model-inputs/mobile-climate.json"
CLUTTER_DIR = ROOT / "model-inputs/mobile-clutter"
OUTPUT_DIR = ROOT / "dist/data"
CACHE_DIR = Path("/tmp/ireland-mobile-terrain-z10")

EARTH_RADIUS_M = 6_371_000.0
WEB_MERCATOR_RADIUS_M = 6_378_137.0
GRID_M = 500.0
DEM_ZOOM = 10
PATH_SPACING_M = 160.0
TX_HEIGHT_M = 30.0
RX_HEIGHT_M = 1.5
LINK_BUDGET_ALLOWANCE_DB = 48.0
UNCERTAINTY_MARGIN_DB = 10.0
MAX_RANGE_M = {"4g": 42_000.0, "5g": 32_000.0}
TECH_BITS = {"4g": 1, "5g": 2}
NETWORK_SLUGS = {"Eir": "eir", "Three": "three", "Vodafone": "vodafone"}
SIGNAL_COLORS = {
    "high": (37, 105, 71, 230),
    "usable": (117, 164, 91, 230),
    "fringe": (227, 165, 72, 230),
    "gap": (200, 81, 70, 225),
    "unknown": (139, 147, 142, 220),
}
SIGNAL_THRESHOLDS = (-95.0, -105.0, -115.0)
CLUTTER_TILE_M = 163_840.0
CLUTTER_PIXEL_M = 160.0
CLUTTER_TILE_PX = 1_024
TERRAIN_URL = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"
_MODEL_WORKER_STATE = None


def log(message: str) -> None:
    print(message, flush=True)


def lon_lat_to_mercator(lon: float, lat: float) -> tuple[float, float]:
    lat = max(-85.05112878, min(85.05112878, lat))
    return (
        WEB_MERCATOR_RADIUS_M * math.radians(lon),
        WEB_MERCATOR_RADIUS_M * math.log(math.tan(math.pi / 4 + math.radians(lat) / 2)),
    )


def mercator_to_lon_lat(x: float, y: float) -> tuple[float, float]:
    return (
        math.degrees(x / WEB_MERCATOR_RADIUS_M),
        math.degrees(2 * math.atan(math.exp(y / WEB_MERCATOR_RADIUS_M)) - math.pi / 2),
    )


def mercator_to_tile_pixel(x: float, y: float) -> tuple[float, float]:
    n = 2**DEM_ZOOM
    px = (x / (2 * math.pi * WEB_MERCATOR_RADIUS_M) + 0.5) * n * 256
    py = (0.5 - y / (2 * math.pi * WEB_MERCATOR_RADIUS_M)) * n * 256
    return px, py


def geometry_rings(geometry: dict):
    if geometry["type"] == "Polygon":
        yield from geometry["coordinates"]
    elif geometry["type"] == "MultiPolygon":
        for polygon in geometry["coordinates"]:
            yield from polygon
    else:
        raise ValueError(f"Unsupported boundary geometry: {geometry['type']}")


def read_boundary() -> tuple[dict, list[list[list[float]]], tuple[float, float, float, float]]:
    boundary = json.loads(BOUNDARY_PATH.read_text())
    feature = next(f for f in boundary["features"] if f.get("properties", {}).get("name") == "Ireland")
    geometry = feature["geometry"]
    polygons = geometry["coordinates"] if geometry["type"] == "MultiPolygon" else [geometry["coordinates"]]
    points = [point for polygon in polygons for ring in polygon for point in ring]
    bounds = (
        min(point[0] for point in points), min(point[1] for point in points),
        max(point[0] for point in points), max(point[1] for point in points),
    )
    return boundary, polygons, bounds


def grid_for_boundary(polygons: list, bounds: tuple[float, float, float, float]):
    min_lon, min_lat, max_lon, max_lat = bounds
    projected = [lon_lat_to_mercator(lon, lat) for lon, lat in [
        (min_lon, min_lat), (min_lon, max_lat), (max_lon, min_lat), (max_lon, max_lat)
    ]]
    min_x = min(p[0] for p in projected)
    max_x = max(p[0] for p in projected)
    min_y = min(p[1] for p in projected)
    max_y = max(p[1] for p in projected)
    width = math.ceil((max_x - min_x) / GRID_M)
    height = math.ceil((max_y - min_y) / GRID_M)
    mask_image = Image.new("L", (width, height), 0)
    draw = ImageDraw.Draw(mask_image)
    for polygon in polygons:
        if not polygon:
            continue
        outer, *holes = polygon
        outer_xy = [
            ((lon_lat_to_mercator(lon, lat)[0] - min_x) / GRID_M,
             (max_y - lon_lat_to_mercator(lon, lat)[1]) / GRID_M)
            for lon, lat in outer
        ]
        draw.polygon(outer_xy, fill=1)
        for ring in holes:
            hole_xy = [
                ((lon_lat_to_mercator(lon, lat)[0] - min_x) / GRID_M,
                 (max_y - lon_lat_to_mercator(lon, lat)[1]) / GRID_M)
                for lon, lat in ring
            ]
            draw.polygon(hole_xy, fill=0)
    mask = np.asarray(mask_image, dtype=bool)
    rows, columns = np.nonzero(mask)
    xs = min_x + (columns + 0.5) * GRID_M
    ys = max_y - (rows + 0.5) * GRID_M
    lons = xs / WEB_MERCATOR_RADIUS_M * 180 / math.pi
    lats = (2 * np.arctan(np.exp(ys / WEB_MERCATOR_RADIUS_M)) - math.pi / 2) * 180 / math.pi
    image_coordinates = [
        list(mercator_to_lon_lat(min_x, max_y)),
        list(mercator_to_lon_lat(max_x, max_y)),
        list(mercator_to_lon_lat(max_x, min_y)),
        list(mercator_to_lon_lat(min_x, min_y)),
    ]
    grid = {
        "mask": mask,
        "rows": rows,
        "columns": columns,
        "lons": lons,
        "lats": lats,
        "width": width,
        "height": height,
        "min_x": min_x,
        "max_y": max_y,
        "image_coordinates": image_coordinates,
        "bounds_4326": [min_lon, min_lat, max_lon, max_lat],
    }
    return grid


def make_cell_sites_geojson(sites_data: dict) -> None:
    features = []
    for record_index, (network_index, lon, lat, bands) in enumerate(sites_data["records"]):
        frequencies = sorted({int(row[0]) for row in bands})
        technologies = []
        if any(int(row[2]) & TECH_BITS["4g"] for row in bands):
            technologies.append("4G")
        if any(int(row[2]) & TECH_BITS["5g"] for row in bands):
            technologies.append("5G")
        network = sites_data["networks"][network_index]
        band_rows = [
            {"frequencyMHz": int(freq), "maxEirpDbm": float(eirp),
             "technology": [name for name, bit in (("4G", 1), ("5G", 2)) if int(flags) & bit]}
            for freq, eirp, flags in bands
        ]
        features.append({
            "type": "Feature",
            "id": f"{network_index}-{record_index}",
            "properties": {
                "network": network,
                "networkIndex": int(network_index),
                "technologies": technologies,
                "frequenciesMHz": frequencies,
                "bandsJson": json.dumps(band_rows, separators=(",", ":")),
            },
            "geometry": {"type": "Point", "coordinates": [float(lon), float(lat)]},
        })
    target = OUTPUT_DIR / "cell-sites.geojson"
    target.write_text(json.dumps({"type": "FeatureCollection", "features": features}, separators=(",", ":")) + "\n")
    log(f"Wrote {len(features):,} ComReg site records ({target.stat().st_size / 1024 / 1024:.1f} MiB GeoJSON)")


def load_climate() -> dict:
    return json.loads(CLIMATE_PATH.read_text())


def climatology_at(lon: float, lat: float, climate: dict) -> tuple[float, float] | None:
    longitude = (lon % 360 + 360) % 360
    step = float(climate.get("gridResolutionDegrees") or 1.5)
    row_position = (float(climate["latitudes"][0]) - lat) / step
    column_position = (longitude - float(climate["longitudes"][0])) / step
    row0, column0 = math.floor(row_position), math.floor(column_position)
    row1, column1 = row0 + 1, column0 + 1
    if row0 < 0 or column0 < 0 or row1 >= len(climate["latitudes"]) or column1 >= len(climate["longitudes"]):
        return None
    row_fraction = row_position - row0
    column_fraction = column_position - column0

    def interpolate(grid: list[list[float]]) -> float:
        top = grid[row0][column0] * (1 - column_fraction) + grid[row0][column1] * column_fraction
        bottom = grid[row1][column0] * (1 - column_fraction) + grid[row1][column1] * column_fraction
        return top * (1 - row_fraction) + bottom * row_fraction

    return interpolate(climate["deltaN"]), interpolate(climate["n0"])


def fetch_elevation_tile(tile_x: int, tile_y: int) -> tuple[tuple[int, int], bytes | None]:
    key = (tile_x, tile_y)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = CACHE_DIR / f"{DEM_ZOOM}_{tile_x}_{tile_y}.png"
    if path.exists():
        return key, path.read_bytes()
    url = TERRAIN_URL.format(z=DEM_ZOOM, x=tile_x, y=tile_y)
    request = urllib.request.Request(url, headers={"User-Agent": "IrelandMobileSignalMap/1.0"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=25) as response:
                data = response.read()
            path.write_bytes(data)
            return key, data
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return key, None
            if attempt == 2:
                raise
        except Exception:
            if attempt == 2:
                raise
        time.sleep(0.5 * (attempt + 1))
    return key, None


def load_elevation_tiles(grid: dict) -> dict:
    min_lon, min_lat, max_lon, max_lat = grid["bounds_4326"]
    # 0.8 degrees safely exceeds the model's 42 km maximum path range.
    buffered = (min_lon - 0.8, min_lat - 0.8, max_lon + 0.8, max_lat + 0.8)
    xs, ys = [], []
    for lon in (buffered[0], buffered[2]):
        for lat in (buffered[1], buffered[3]):
            x, y = lon_lat_to_mercator(lon, lat)
            px, py = mercator_to_tile_pixel(x, y)
            xs.append(int(math.floor(px / 256)))
            ys.append(int(math.floor(py / 256)))
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    keys = [(x, y) for y in range(min_y, max_y + 1) for x in range(min_x, max_x + 1)]
    log(f"Loading {len(keys)} elevation tiles at zoom {DEM_ZOOM} (about 90 m ground pixels)")
    result = {}
    with ThreadPoolExecutor(max_workers=16) as pool:
        futures = [pool.submit(fetch_elevation_tile, x, y) for x, y in keys]
        for index, future in enumerate(as_completed(futures), start=1):
            key, data = future.result()
            if data:
                rgba = np.asarray(Image.open(BytesIO(data)).convert("RGBA"), dtype=np.uint8)
                elevation = (rgba[:, :, 0].astype(np.float32) * 256
                             + rgba[:, :, 1].astype(np.float32)
                             + rgba[:, :, 2].astype(np.float32) / 256 - 32768)
                valid = (rgba[:, :, 3] > 0) & (elevation > -1000) & (elevation < 10000)
                result[key] = (elevation, valid)
            if index % 40 == 0 or index == len(keys):
                log(f"  elevation tiles {index}/{len(keys)}")
    log(f"Loaded {len(result):,} elevation tiles")
    return result


def sample_elevation(lons: np.ndarray, lats: np.ndarray, tiles: dict) -> np.ndarray:
    x = WEB_MERCATOR_RADIUS_M * np.radians(lons)
    lat_radians = np.radians(np.clip(lats, -85.05112878, 85.05112878))
    y = WEB_MERCATOR_RADIUS_M * np.log(np.tan(np.pi / 4 + lat_radians / 2))
    px, py = mercator_to_tile_pixel_array(x, y)
    global_columns = np.floor(px).astype(np.int64)
    global_rows = np.floor(py).astype(np.int64)
    tile_x = np.floor_divide(global_columns, 256)
    tile_y = np.floor_divide(global_rows, 256)
    pixel_x = np.mod(global_columns, 256)
    pixel_y = np.mod(global_rows, 256)
    output = np.full(len(lons), np.nan, dtype=np.float32)
    pairs = np.column_stack((tile_x, tile_y))
    for tile_pair in np.unique(pairs, axis=0):
        key = (int(tile_pair[0]), int(tile_pair[1]))
        mask = (tile_x == key[0]) & (tile_y == key[1])
        tile = tiles.get(key)
        if tile is None:
            continue
        elevation, valid = tile
        rows, columns = pixel_y[mask], pixel_x[mask]
        tile_valid = valid[rows, columns]
        output_indices = np.flatnonzero(mask)
        output[output_indices[tile_valid]] = elevation[rows[tile_valid], columns[tile_valid]]
    return output


def mercator_to_tile_pixel_array(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n_pixels = (2**DEM_ZOOM) * 256
    px = (x / (2 * math.pi * WEB_MERCATOR_RADIUS_M) + 0.5) * n_pixels
    py = (0.5 - y / (2 * math.pi * WEB_MERCATOR_RADIUS_M)) * n_pixels
    return px, py


class ClutterSampler:
    def __init__(self):
        metadata = json.loads((CLUTTER_DIR / "metadata.json").read_text())
        self.tiles = set(metadata["tiles"])
        self.cache: dict[str, np.ndarray | None] = {}

    def _get_tile(self, filename: str):
        if filename in self.cache:
            return self.cache[filename]
        path = CLUTTER_DIR / filename
        if filename not in self.tiles or not path.exists():
            self.cache[filename] = None
            return None
        with Image.open(path) as image:
            values = np.asarray(image.convert("RGBA"), dtype=np.uint8)[:, :, 0]
        self.cache[filename] = values
        return values

    def sample(self, lons: np.ndarray, lats: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        lat_radians = np.radians(np.clip(lats, -85.05112878, 85.05112878))
        world_x = WEB_MERCATOR_RADIUS_M * np.radians(lons)
        world_y = WEB_MERCATOR_RADIUS_M * np.log(np.tan(np.pi / 4 + lat_radians / 2))
        tile_x = np.floor(world_x / CLUTTER_TILE_M).astype(np.int64)
        tile_y = np.floor(world_y / CLUTTER_TILE_M).astype(np.int64)
        pixel_x = np.floor((world_x - tile_x * CLUTTER_TILE_M) / CLUTTER_PIXEL_M).astype(np.int64)
        pixel_y = np.floor((tile_y * CLUTTER_TILE_M + CLUTTER_TILE_M - world_y) / CLUTTER_PIXEL_M).astype(np.int64)
        pixel_x = np.clip(pixel_x, 0, CLUTTER_TILE_PX - 1)
        pixel_y = np.clip(pixel_y, 0, CLUTTER_TILE_PX - 1)
        heights = np.full(len(lons), np.nan, dtype=np.float32)
        sea = np.zeros(len(lons), dtype=bool)
        pairs = np.column_stack((tile_x, tile_y))
        for pair in np.unique(pairs, axis=0):
            tx, ty = int(pair[0]), int(pair[1])
            mask = (tile_x == tx) & (tile_y == ty)
            values = self._get_tile(f"{tx}_{ty}.png")
            if values is None:
                continue
            codes = values[pixel_y[mask], pixel_x[mask]]
            indices = np.flatnonzero(mask)
            known = codes != 255
            heights[indices[known]] = np.where(codes[known] == 30, 0, codes[known])
            sea[indices[known]] = codes[known] == 30
        return heights, sea


def p1812_knife_edge_loss(v: float) -> float:
    if not math.isfinite(v) or v <= -0.78:
        return 0.0
    term = math.sqrt((v - 0.1) ** 2 + 1) + v - 0.1
    return 6.9 + 20 * math.log10(max(term, 1e-12))


def p1812_bullington_loss(profile: list[float], htc: float, hrc: float,
                          distance_km: float, wavelength_m: float,
                          effective_earth_km: float) -> float:
    curvature = 1 / effective_earth_km
    transmitter_slope = -math.inf
    receiver_slope = -math.inf
    samples = []
    for index in range(1, len(profile) - 1):
        di = distance_km * index / (len(profile) - 1)
        height = profile[index]
        curved_height = height + 500 * curvature * di * (distance_km - di)
        slope_t = (curved_height - htc) / di
        slope_r = (curved_height - hrc) / (distance_km - di)
        transmitter_slope = max(transmitter_slope, slope_t)
        receiver_slope = max(receiver_slope, slope_r)
        samples.append((di, height, curved_height))
    direct_slope = (hrc - htc) / distance_km
    knife_edge_v = -math.inf
    if transmitter_slope < direct_slope:
        for di, height, _ in samples:
            ray_height = (htc * (distance_km - di) + hrc * di) / distance_km
            v = (height + 500 * curvature * di * (distance_km - di) - ray_height) * math.sqrt(
                0.002 * distance_km / (wavelength_m * di * (distance_km - di))
            )
            knife_edge_v = max(knife_edge_v, v)
    else:
        denominator = transmitter_slope + receiver_slope
        if not math.isfinite(denominator) or abs(denominator) < 1e-12:
            return 0.0
        bullington_point_km = (hrc - htc + receiver_slope * distance_km) / denominator
        d_bp = max(0.001, min(distance_km - 0.001, bullington_point_km))
        ray_height = (htc * (distance_km - d_bp) + hrc * d_bp) / distance_km
        knife_edge_v = (htc + transmitter_slope * d_bp - ray_height) * math.sqrt(
            0.002 * distance_km / (wavelength_m * d_bp * (distance_km - d_bp))
        )
    uncorrected = p1812_knife_edge_loss(knife_edge_v)
    return uncorrected + (1 - math.exp(-uncorrected / 6)) * (10 + 0.02 * distance_km)


def p1812_first_term_spherical_loss(distance_km: float, frequency_ghz: float,
                                    effective_earth_km: float, hte: float, hre: float,
                                    sea_fraction: float) -> float:
    if not (distance_km > 0 and frequency_ghz > 0 and effective_earth_km > 0 and hte > 0 and hre > 0):
        return 0.0
    wavelength_m = 0.299792458 / frequency_ghz

    def first_term_for_surface(earth_km: float, epsilon: float, conductivity: float) -> float:
        kh = 0.036 * (earth_km * frequency_ghz) ** (-1 / 3) * (
            (epsilon - 1) ** 2 + (18 * conductivity / frequency_ghz) ** 2
        ) ** (-1 / 4)
        beta = 1.0
        x = 21.88 * beta * (frequency_ghz / earth_km**2) ** (1 / 3) * distance_km
        yt = 0.9575 * beta * (frequency_ghz**2 / earth_km) ** (1 / 3) * hte
        yr = 0.9575 * beta * (frequency_ghz**2 / earth_km) ** (1 / 3) * hre
        fx = (11 + 10 * math.log10(x) - 17.6 * x) if x >= 1.6 else (
            -20 * math.log10(max(x, 1e-12)) - 5.6488 * x**1.425
        )

        def g(height: float) -> float:
            b = beta * height
            value = (17.6 * math.sqrt(b - 1.1) - 5 * math.log10(b - 1.1) - 8) if b > 2 else (
                20 * math.log10(max(b + 0.1 * b**3, 1e-12))
            )
            return max(value, 2 + 20 * math.log10(kh))

        return -fx - g(yt) - g(yr)

    d_los = math.sqrt(2 * effective_earth_km) * (math.sqrt(0.001 * hte) + math.sqrt(0.001 * hre))
    if distance_km >= d_los:
        land = first_term_for_surface(effective_earth_km, 22, 0.003)
        sea = first_term_for_surface(effective_earth_km, 80, 5)
        return sea_fraction * sea + (1 - sea_fraction) * land
    sum_heights = hte + hre
    mc = 250 * distance_km**2 / (effective_earth_km * sum_heights)
    if not mc > 0:
        return 0.0
    c = (hte - hre) / sum_heights
    acos_argument = (3 * c / 2) * math.sqrt(3 * mc / (mc + 1) ** 3)
    b = 2 * math.sqrt((mc + 1) / (3 * mc)) * math.cos(
        math.pi / 3 + math.acos(max(-1, min(1, acos_argument))) / 3
    )
    d_se1 = distance_km / 2 * (1 + b)
    d_se2 = distance_km - d_se1
    if not (d_se1 > 0 and d_se2 > 0):
        return 0.0
    h_se = ((hte - 500 * d_se1**2 / effective_earth_km) * d_se2
            + (hre - 500 * d_se2**2 / effective_earth_km) * d_se1) / distance_km
    h_req = 17.456 * math.sqrt(d_se1 * d_se2 * wavelength_m / distance_km)
    if h_se > h_req or not h_req > 0:
        return 0.0
    modified_earth_km = 500 * (distance_km / (math.sqrt(hte) + math.sqrt(hre))) ** 2
    land = first_term_for_surface(modified_earth_km, 22, 0.003)
    sea = first_term_for_surface(modified_earth_km, 80, 5)
    diffraction = sea_fraction * sea + (1 - sea_fraction) * land
    if diffraction < 0:
        return 0.0
    return max(0.0, (1 - h_se / h_req) * diffraction)


def p1812_delta_bullington(terrain: np.ndarray, clutter: np.ndarray, distance_m: float,
                           frequency_mhz: float, delta_n: float,
                           sea_flags: np.ndarray, transmitter_height_m: float = TX_HEIGHT_M) -> float | None:
    if len(terrain) < 3 or len(terrain) != len(clutter) or not math.isfinite(distance_m) or distance_m <= 0:
        return None
    if not np.isfinite(terrain).all() or not np.isfinite(clutter).all():
        return None
    frequency_ghz = frequency_mhz / 1000
    wavelength_m = 0.299792458 / frequency_ghz
    effective_earth_km = 157 / (157 - delta_n) * 6371
    if not (frequency_ghz > 0 and effective_earth_km > 0):
        return None
    distance_km = distance_m / 1000
    htc = float(terrain[0]) + transmitter_height_m
    hrc = float(terrain[-1]) + RX_HEIGHT_M
    surface = terrain.astype(np.float64).copy()
    surface[1:-1] += clutter[1:-1]
    bullington_actual = p1812_bullington_loss(surface.tolist(), htc, hrc, distance_km, wavelength_m, effective_earth_km)
    v1, v2 = 0.0, 0.0
    count = len(terrain)
    for index in range(1, count):
        d0 = distance_km * (index - 1) / (count - 1)
        d1 = distance_km * index / (count - 1)
        delta_d = d1 - d0
        v1 += delta_d * (terrain[index] + terrain[index - 1])
        v2 += delta_d * (terrain[index] * (2 * d1 + d0) + terrain[index - 1] * (d1 + 2 * d0))
    hst = (2 * v1 * distance_km - v2) / distance_km**2
    hsr = (v2 - v1 * distance_km) / distance_km**2
    h_obs, alpha_obt, alpha_obr = -math.inf, -math.inf, -math.inf
    for index in range(1, count - 1):
        di = distance_km * index / (count - 1)
        obstruction = terrain[index] - (htc * (distance_km - di) + hrc * di) / distance_km
        h_obs = max(h_obs, obstruction)
        alpha_obt = max(alpha_obt, obstruction / di)
        alpha_obr = max(alpha_obr, obstruction / (distance_km - di))
    ratio = 0.5 if alpha_obt + alpha_obr == 0 else alpha_obt / (alpha_obt + alpha_obr)
    hst_prime = hst if h_obs <= 0 else hst - h_obs * ratio
    hsr_prime = hsr if h_obs <= 0 else hsr - h_obs * (1 - ratio)
    hstd = min(hst_prime, terrain[0])
    hsrd = min(hsr_prime, terrain[-1])
    htc_smooth = htc - hstd
    hrc_smooth = hrc - hsrd
    bullington_smooth = p1812_bullington_loss(
        np.zeros(count).tolist(), htc_smooth, hrc_smooth, distance_km, wavelength_m, effective_earth_km
    )
    sea_fraction = float(np.count_nonzero(sea_flags)) / len(sea_flags) if len(sea_flags) else 0.0
    spherical = p1812_first_term_spherical_loss(
        distance_km, frequency_ghz, effective_earth_km, htc_smooth, hrc_smooth, sea_fraction
    )
    return bullington_actual + max(spherical - bullington_smooth, 0.0)


def load_site_height_estimates(sites_data: dict) -> tuple[list[dict], dict]:
    """Load planning-derived top-antenna estimates aligned to ComReg record order."""
    networks = sites_data["networks"]
    records = sites_data["records"]
    with SITE_HEIGHTS_PATH.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != len(records):
        raise ValueError(f"Site-height input has {len(rows)} rows; expected {len(records)}")
    aligned = []
    status_counts: dict[str, int] = {}
    applied = 0
    for index, (row, site) in enumerate(zip(rows, records)):
        if int(row["recordIndex"]) != index:
            raise ValueError(f"Site-height record order mismatch at row {index}")
        if (row["network"] != networks[int(site[0])] or
                abs(float(row["longitude"]) - float(site[1])) > 1e-6 or
                abs(float(row["latitude"]) - float(site[2])) > 1e-6):
            raise ValueError(f"Site-height identity/coordinate mismatch at record {index}")
        raw_height = row["estimatedAntennaHeightM"].strip()
        height = TX_HEIGHT_M if not raw_height else float(raw_height)
        if not math.isfinite(height) or not 1.0 <= height <= 120.0:
            raise ValueError(f"Invalid transmitter height {height} at record {index}")
        status = row["estimateStatus"] or "default_30m"
        if raw_height:
            applied += 1
            status_counts[status] = status_counts.get(status, 0) + 1
        aligned.append({"height_m": height, "status": status})
    summary = {
        "source": "model-inputs/mobile-site-height-estimates.csv",
        "method": "For a planning match within 30 m, use stated mast/support structure height minus 1 m as an estimated highest/top-mounted antenna height; otherwise retain the 30 m default.",
        "defaultHeightM": TX_HEIGHT_M,
        "recordsTotal": len(records),
        "recordsWithEstimate": applied,
        "recordsUsingDefault": len(records) - applied,
        "heightStatusCounts": status_counts,
        "planningMatchMaximumDistanceM": 30,
        "identityCaveat": "Planning/site matches are largely automated coordinate/text candidates and are not individually verified; proposed and ambiguous records are retained with status labels.",
    }
    return aligned, summary


def build_site_band_indexes(sites_data: dict, site_heights: list[dict] | None = None) -> dict:
    indexes = {}
    for technology, bit in TECH_BITS.items():
        indexes[technology] = {}
        for network_index, network in enumerate(sites_data["networks"]):
            lons, lats, frequency, base_power, tx_heights = [], [], [], [], []
            for record_index, site in enumerate(sites_data["records"]):
                if site[0] != network_index:
                    continue
                best = None
                for freq_mhz, eirp_dbm, flags in site[3]:
                    if not (int(flags) & bit) or freq_mhz <= 0:
                        continue
                    constant = float(eirp_dbm) - 32.45 - 20 * math.log10(float(freq_mhz)) - LINK_BUDGET_ALLOWANCE_DB
                    if best is None or constant > best[0]:
                        best = (constant, float(freq_mhz))
                if best is not None:
                    lons.append(float(site[1]))
                    lats.append(float(site[2]))
                    frequency.append(best[1])
                    base_power.append(best[0])
                    tx_heights.append(float(site_heights[record_index]["height_m"]) if site_heights else TX_HEIGHT_M)
            indexes[technology][network] = {
                "lon": np.asarray(lons, dtype=np.float64),
                "lat": np.asarray(lats, dtype=np.float64),
                "frequency": np.asarray(frequency, dtype=np.float64),
                "base_power": np.asarray(base_power, dtype=np.float64),
                "tx_height": np.asarray(tx_heights, dtype=np.float64),
            }
            log(f"{network} {technology.upper()}: {len(lons):,} site-band candidates")
    return indexes


def haversine(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    lat1_r, lat2_r = math.radians(lat1), math.radians(lat2)
    dlat = lat2_r - lat1_r
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1_r) * math.cos(lat2_r) * math.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_M * math.asin(math.sqrt(max(0.0, min(1.0, a))))


def candidates_at(lon: float, lat: float, sites: dict, technology: str) -> list[dict]:
    radius = MAX_RANGE_M[technology]
    cos_lat = max(0.35, math.cos(math.radians(lat)))
    per_network = []
    for network, data in sites[technology].items():
        lat_difference = (data["lat"] - lat) * 111_320.0
        if not len(data["lat"]):
            per_network.append([])
            continue
        within_lat = np.abs(lat_difference) <= radius
        if not within_lat.any():
            per_network.append([])
            continue
        indices = np.flatnonzero(within_lat)
        lon_difference = (data["lon"][indices] - lon) * 111_320.0 * cos_lat
        approximate_distance = np.hypot(lat_difference[indices], lon_difference)
        nearby = indices[approximate_distance <= radius * 1.003]
        if not len(nearby):
            per_network.append([])
            continue
        distance_m = np.maximum(approximate_distance[np.isin(indices, nearby)], 1.0)
        nominal = data["base_power"][nearby] - 20 * np.log10(np.maximum(distance_m / 1000, 0.05))
        within_range = distance_m <= radius * 1.003
        nearby = nearby[within_range]
        nominal = nominal[within_range]
        if not len(nearby):
            per_network.append([])
            continue
        count = min(2, len(nearby))
        selected_positions = np.argpartition(nominal, -count)[-count:]
        selected_positions = selected_positions[np.argsort(nominal[selected_positions])[::-1]]
        selected = []
        top_nominal = float(nominal[selected_positions[0]])
        for position in selected_positions:
            if top_nominal - float(nominal[position]) > 18:
                continue
            site_index = int(nearby[position])
            site_lon, site_lat = float(data["lon"][site_index]), float(data["lat"][site_index])
            exact_distance = haversine(lon, lat, site_lon, site_lat)
            if exact_distance > radius:
                continue
            exact_nominal = float(data["base_power"][site_index]) - 20 * math.log10(max(exact_distance / 1000, 0.05))
            selected.append({
                "network": network,
                "lon": site_lon,
                "lat": site_lat,
                "distance_m": exact_distance,
                "nominal_dbm": exact_nominal,
                "frequency_mhz": float(data["frequency"][site_index]),
                "tx_height_m": float(data["tx_height"][site_index]),
            })
        per_network.append(selected)
    return [candidate for network_candidates in per_network for candidate in network_candidates]


def link_signal(candidate: dict, receiver_lon: float, receiver_lat: float,
                terrain_tiles: dict, clutter_sampler: ClutterSampler, climate: dict) -> float | None:
    distance_m = candidate["distance_m"]
    steps = max(2, min(320, math.ceil(distance_m / PATH_SPACING_M)))
    fractions = np.linspace(0, 1, steps + 1, dtype=np.float64)
    lons = candidate["lon"] + (receiver_lon - candidate["lon"]) * fractions
    lats = candidate["lat"] + (receiver_lat - candidate["lat"]) * fractions
    elevations = sample_elevation(lons, lats, terrain_tiles)
    clutter, sea = clutter_sampler.sample(lons, lats)
    if not np.isfinite(elevations).all() or not np.isfinite(clutter).all():
        return None
    midpoint_lon = (candidate["lon"] + receiver_lon) / 2
    midpoint_lat = (candidate["lat"] + receiver_lat) / 2
    climatology = climatology_at(midpoint_lon, midpoint_lat, climate)
    if climatology is None:
        return None
    diffraction_db = p1812_delta_bullington(
        elevations, clutter, distance_m, candidate["frequency_mhz"], climatology[0], sea,
        candidate.get("tx_height_m", TX_HEIGHT_M)
    )
    if diffraction_db is None or not math.isfinite(diffraction_db):
        return None
    return candidate["nominal_dbm"] - diffraction_db


def signal_band(signal_dbm: float) -> str:
    if not math.isfinite(signal_dbm):
        return "gap"
    if signal_dbm >= -95:
        return "high"
    if signal_dbm >= -105:
        return "usable"
    if signal_dbm >= -115:
        return "fringe"
    return "gap"


def output_rgba(signal: np.ndarray, mask: np.ndarray) -> tuple[Image.Image, Image.Image]:
    height, width = signal.shape
    coverage = np.zeros((height, width, 4), dtype=np.uint8)
    sensitivity = np.zeros((height, width, 4), dtype=np.uint8)
    finite = np.isfinite(signal)
    unknown = np.isnan(signal) & mask
    gap = np.isneginf(signal) & mask
    coverage[unknown] = SIGNAL_COLORS["unknown"]
    coverage[gap] = SIGNAL_COLORS["gap"]
    for class_name, predicate in (
        ("high", finite & (signal >= -95)),
        ("usable", finite & (signal >= -105) & (signal < -95)),
        ("fringe", finite & (signal >= -115) & (signal < -105)),
        ("gap", finite & (signal < -115)),
    ):
        coverage[predicate & mask] = SIGNAL_COLORS[class_name]
    near_cutoff = finite & mask & np.any(
        np.abs(signal[:, :, None] - np.asarray(SIGNAL_THRESHOLDS, dtype=np.float32)) <= UNCERTAINTY_MARGIN_DB,
        axis=2,
    )
    rows, columns = np.indices((height, width))
    stripes = near_cutoff & (((rows + columns) % 7) < 3)
    sensitivity[stripes] = (22, 34, 36, 145)
    return Image.fromarray(coverage, "RGBA"), Image.fromarray(sensitivity, "RGBA")


def model_grid_cell(cell: tuple[int, int, float, float]) -> tuple:
    """Model one independent grid cell in a forked worker process."""
    if _MODEL_WORKER_STATE is None:
        raise RuntimeError("Model worker state was not inherited")
    terrain_tiles, clutter_sampler, climate, site_indexes, networks = _MODEL_WORKER_STATE
    row, column, lon, lat = cell
    cell_signals = {}
    for technology in TECH_BITS:
        candidates = candidates_at(lon, lat, site_indexes, technology)
        by_network = {}
        for network in networks:
            best = -math.inf
            found_candidate = False
            had_valid_profile = False
            for candidate in candidates:
                if candidate["network"] != network:
                    continue
                found_candidate = True
                estimate = link_signal(candidate, lon, lat, terrain_tiles, clutter_sampler, climate)
                if estimate is None:
                    continue
                had_valid_profile = True
                best = max(best, estimate)
            if not found_candidate:
                value = -math.inf
            elif had_valid_profile:
                value = best
            else:
                value = math.nan
            by_network[network] = value
            cell_signals[(network, technology)] = value
        finite_values = [value for value in by_network.values() if math.isfinite(value)]
        if finite_values:
            aggregate = max(finite_values)
        elif any(math.isnan(value) for value in by_network.values()):
            aggregate = math.nan
        else:
            aggregate = -math.inf
        cell_signals[("all", technology)] = aggregate
    output_keys = [(network, technology) for technology in TECH_BITS for network in ["all", *networks]]
    return row, column, *(cell_signals[key] for key in output_keys)


def write_metadata(sites_data: dict, grid: dict, land_cells: int,
                   height_metadata: dict | None = None) -> None:
    info = {
        "title": "Ireland mobile signal estimate",
        "siteSource": sites_data["source"],
        "siteSourceUrl": sites_data["sourceUrl"],
        "siteCounts": sites_data["counts"],
        "siteRecordCount": len(sites_data["records"]),
        "networks": sites_data["networks"],
        "technologyBits": sites_data["technologyBits"],
        "gridSpacingProjectedM": GRID_M,
        "approximateGroundResolutionM": round(GRID_M * math.cos(math.radians(53.5))),
        "landGridCells": land_cells,
        "imageCoordinates": grid["image_coordinates"],
        "imageWidth": grid["width"],
        "imageHeight": grid["height"],
        "coverageFiles": [f"coverage-{provider}-{tech}.png" for provider in ["all", "eir", "three", "vodafone"] for tech in ["4g", "5g"]],
        "sensitivityFiles": [f"sensitivity-{provider}-{tech}.png" for provider in ["all", "eir", "three", "vodafone"] for tech in ["4g", "5g"]],
        "signalBandsDbm": {"high": ">= -95", "usable": "-105 to < -95", "fringe": "-115 to < -105", "gap": "< -115"},
        "sensitivityRule": f"Striped cells lie within {UNCERTAINTY_MARGIN_DB} dB of a signal-band cutoff; this is a threshold-sensitivity flag, not a confidence interval.",
        "model": {
            "propagation": "ITU-R P.1812-8 Δ-Bullington diffraction and spherical-Earth correction",
            "terrainSource": "Mapzen Terrain Tiles / AWS elevation-tiles-prod, zoom 10, about 90 m ground pixels",
            "clutterSource": "Copernicus / EEA CORINE Land Cover 2018, representative 160 m clutter categories",
            "climateSource": "ITU-R P.1812-8 median annual ΔN grid",
            "transmitterHeightM": TX_HEIGHT_M,
            "transmitterHeightDefaultM": TX_HEIGHT_M,
            "receiverHeightM": RX_HEIGHT_M,
            "linkBudgetAllowanceDb": LINK_BUDGET_ALLOWANCE_DB,
            "maxRangeM": MAX_RANGE_M,
            "pathProfileSpacingM": PATH_SPACING_M,
            "weatherAdjustment": False,
            "notes": [
                "Licensed maximum EIRP and band schedules are not measurements of active handset coverage.",
                "A planning-informed estimated highest antenna height is used only where a matched planning description states a support-structure height; the estimate is structure height minus 1 m. Other records retain the 30 m default. This is not verified current equipment; proposal and ambiguous statuses remain labelled in the input.",
                "The 48 dB link allowance is uncalibrated; map classes are indicative and do not have a statistical confidence level.",
                "The national image grid is for country-scale display; the underlying route model remains more detailed along selected trails."
            ]
        },
        "boundarySource": "Natural Earth 1:10m Admin 0 – Countries, v5.1.1, Ireland polygon",
        "boundarySourceUrl": "https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/",
    }
    if height_metadata is not None:
        info["transmitterHeightEstimates"] = height_metadata
    (OUTPUT_DIR / "coverage-metadata.json").write_text(json.dumps(info, separators=(",", ":")) + "\n")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    sites_data = json.loads(SITES_PATH.read_text())
    site_heights, height_metadata = load_site_height_estimates(sites_data)
    boundary, polygons, bounds = read_boundary()
    grid = grid_for_boundary(polygons, bounds)
    land_cells = len(grid["rows"])
    log(f"Republic grid: {grid['width']} × {grid['height']} pixels; {land_cells:,} in-country cells at {GRID_M:.0f} m projected spacing")
    make_cell_sites_geojson(sites_data)
    terrain_tiles = load_elevation_tiles(grid)
    clutter_sampler = ClutterSampler()
    climate = load_climate()
    site_indexes = build_site_band_indexes(sites_data, site_heights)
    networks = list(sites_data["networks"])
    signals = {
        (network, technology): np.full((grid["height"], grid["width"]), np.nan, dtype=np.float32)
        for network in ["all", *networks] for technology in TECH_BITS
    }
    for technology in TECH_BITS:
        signals[("all", technology)][grid["rows"], grid["columns"]] = -np.inf
        for network in networks:
            signals[(network, technology)][grid["rows"], grid["columns"]] = -np.inf

    global _MODEL_WORKER_STATE
    _MODEL_WORKER_STATE = (terrain_tiles, clutter_sampler, climate, site_indexes, networks)
    started = time.monotonic()
    rows, columns = grid["rows"], grid["columns"]
    cells = ((int(row), int(column), float(lon), float(lat))
             for row, column, lon, lat in zip(rows, columns, grid["lons"], grid["lats"]))
    output_keys = [(network, technology) for technology in TECH_BITS for network in ["all", *networks]]
    worker_count = min(8, multiprocessing.cpu_count())
    log(f"Modelling with {worker_count} worker processes")
    # Bound queued work and memory while keeping each worker busy between updates.
    with ProcessPoolExecutor(max_workers=worker_count, mp_context=multiprocessing.get_context("fork")) as pool:
        cell_index = 0
        while batch := list(islice(cells, 8_192)):
            for result in pool.map(model_grid_cell, batch, chunksize=32):
                row, column, *cell_values = result
                for key, value in zip(output_keys, cell_values):
                    signals[key][row, column] = value
                cell_index += 1
                if cell_index % 25_000 == 0 or cell_index == land_cells:
                    elapsed = time.monotonic() - started
                    log(f"  modelled {cell_index:,}/{land_cells:,} grid cells ({elapsed / 60:.1f} min)")

    for (network, technology), signal in signals.items():
        provider = "all" if network == "all" else NETWORK_SLUGS[network]
        coverage_image, sensitivity_image = output_rgba(signal, grid["mask"])
        coverage_path = OUTPUT_DIR / f"coverage-{provider}-{technology}.png"
        sensitivity_path = OUTPUT_DIR / f"sensitivity-{provider}-{technology}.png"
        coverage_image.save(coverage_path, optimize=True)
        sensitivity_image.save(sensitivity_path, optimize=True)
        log(f"Wrote {coverage_path.name} ({coverage_path.stat().st_size / 1024:.1f} KiB)")
    write_metadata(sites_data, grid, land_cells, height_metadata)
    log(f"Done in {(time.monotonic() - started) / 60:.1f} minutes")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        log(f"National coverage build failed: {error}")
        raise

#!/usr/bin/env python3
"""
Split an oversized Bambu/3MF project into printer-sized plates.

The p1s bed is roughly 256 mm square.  When a model (for example a 0.4 mm
thick advertising decal) is larger than the printable area this helper
projects the geometry to 2D, tiles it, and extrudes each tile back to 3D so
that every slice fits on a single build plate.

Typical usage:
    python split_3mf.py input.3mf --plate-width-mm 250 --plate-height-mm 250

Dependencies: trimesh, shapely, numpy
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import shapely.affinity
import shapely.geometry
import shapely.ops
import trimesh

# Small floating point tolerance for mm-sized models
EPS = 1e-6


@dataclass
class TileRegion:
    """Metadata describing a single tile extracted from the main footprint."""

    row: int
    col: int
    crop_bounds: Tuple[float, float, float, float]  # xmin, ymin, xmax, ymax in mm
    shape: shapely.geometry.base.BaseGeometry


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split an oversized 3MF/STL file into multiple printer plates by "
            "tiling its XY footprint."
        )
    )
    parser.add_argument(
        "input_file", type=Path, help="Path to the source 3MF or STL file."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Directory where the generated plates are stored "
            "(default: <input>_plates next to the source file)."
        ),
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="Optional prefix for generated files (default: input file stem).",
    )
    parser.add_argument(
        "--plate-width-mm",
        type=float,
        default=250.0,
        help="Printable width of a single bed in millimeters.",
    )
    parser.add_argument(
        "--plate-height-mm",
        type=float,
        default=250.0,
        help="Printable depth of a single bed in millimeters.",
    )
    parser.add_argument(
        "--overlap-mm",
        type=float,
        default=2.0,
        help=(
            "Desired overlap between neighbouring tiles in millimeters. "
            "Overlap helps with post-print trimming."
        ),
    )
    parser.add_argument(
        "--bed-margin-mm",
        type=float,
        default=3.0,
        help=(
            "Margin to keep between the recentered tile and the bed edge "
            "so skirts/brims have space."
        ),
    )
    parser.add_argument(
        "--min-feature-area-mm2",
        type=float,
        default=4.0,
        help=(
            "Ignore tiny scraps that are smaller than this area when "
            "creating tiles."
        ),
    )
    parser.add_argument(
        "--min-normal-z",
        type=float,
        default=0.6,
        help=(
            "Triangles whose face normal has a Z component smaller than this "
            "value are ignored when extracting the planar footprint. "
            "Lower it if the surface is slightly curved."
        ),
    )
    parser.add_argument(
        "--output-format",
        choices=["3mf", "stl"],
        default=None,
        help=(
            "Optional override for the generated file format. "
            "Defaults to the same extension as the input (3MF fallback)."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Calculate the tiling and write the manifest without exporting meshes.",
    )
    return parser.parse_args()


def load_mesh(in_file: Path) -> trimesh.Trimesh:
    loaded = trimesh.load(in_file, force="scene")
    if isinstance(loaded, trimesh.Scene):
        if not loaded.geometry:
            raise ValueError(f"No geometry found inside {in_file}")
        mesh = loaded.dump(concatenate=True)
    else:
        mesh = loaded
    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError(f"Unsupported mesh type {type(mesh)}")
    mesh.remove_unreferenced_vertices()
    return mesh


def mesh_thickness(mesh: trimesh.Trimesh) -> Tuple[float, float, float]:
    bounds = mesh.bounds
    z_min = float(bounds[0][2])
    z_max = float(bounds[1][2])
    thickness = z_max - z_min
    if thickness <= EPS:
        raise ValueError("Model thickness could not be determined (check orientation).")
    return thickness, z_min, z_max


def build_footprint(
    mesh: trimesh.Trimesh,
    min_normal_z: float,
    min_area: float,
) -> shapely.geometry.base.BaseGeometry:
    """Project the upward-facing surface triangles onto the XY plane."""
    triangles = mesh.triangles  # (n, 3, 3)
    normals = mesh.face_normals  # (n, 3)
    polys: List[shapely.geometry.Polygon] = []
    for tri, normal in zip(triangles, normals):
        if normal[2] < min_normal_z:
            continue
        poly = shapely.geometry.Polygon([(float(x), float(y)) for x, y in tri[:, :2]])
        if poly.area <= min_area:
            continue
        polys.append(poly)
    if not polys:
        raise ValueError(
            "No upward oriented faces were found. "
            "Try lowering --min-normal-z or check the model orientation."
        )
    footprint = shapely.ops.unary_union(polys)
    footprint = footprint.buffer(0)  # Clean self-intersections.
    if footprint.area <= EPS:
        raise ValueError("Failed to extract a valid 2D footprint from the mesh.")
    return footprint


def compute_grid_count(length: float, plate: float, stride: float) -> int:
    if length <= plate + EPS:
        return 1
    usable = max(0.0, length - plate)
    return int(math.ceil(usable / stride)) + 1


def generate_tiles(
    footprint: shapely.geometry.base.BaseGeometry,
    plate_w: float,
    plate_h: float,
    overlap: float,
    min_area: float,
) -> List[TileRegion]:
    if overlap >= min(plate_w, plate_h):
        raise ValueError("Overlap must be smaller than the plate dimensions.")
    stride_x = plate_w - overlap
    stride_y = plate_h - overlap
    if stride_x <= 0 or stride_y <= 0:
        raise ValueError("Overlap too large for the plate size.")

    minx, miny, maxx, maxy = footprint.bounds
    width = maxx - minx
    height = maxy - miny
    cols = compute_grid_count(width, plate_w, stride_x)
    rows = compute_grid_count(height, plate_h, stride_y)

    tiles: List[TileRegion] = []
    for row in range(rows):
        for col in range(cols):
            x0 = minx + col * stride_x
            x1 = x0 + plate_w
            y0 = miny + row * stride_y
            y1 = y0 + plate_h
            box = shapely.geometry.box(x0, y0, x1, y1)
            clipped = footprint.intersection(box)
            clipped = clipped.buffer(0)
            clipped = _polygonal_only(clipped)
            if clipped is None or clipped.area < min_area:
                continue
            tiles.append(TileRegion(row=row, col=col, crop_bounds=(x0, y0, x1, y1), shape=clipped))
    if not tiles:
        raise ValueError("No tiles were generated. Check your overlap and area settings.")
    return tiles


def _polygonal_only(
    geom: shapely.geometry.base.BaseGeometry,
) -> Optional[shapely.geometry.base.BaseGeometry]:
    if geom.is_empty:
        return None
    if isinstance(geom, shapely.geometry.Polygon):
        return geom
    if isinstance(geom, shapely.geometry.MultiPolygon):
        return geom
    if isinstance(geom, shapely.geometry.GeometryCollection):
        polys = [g for g in geom.geoms if isinstance(g, shapely.geometry.Polygon)]
        if not polys:
            return None
        if len(polys) == 1:
            return polys[0]
        return shapely.geometry.MultiPolygon(polys)
    return None


def recenter_shape(
    geom: shapely.geometry.base.BaseGeometry,
    plate_w: float,
    plate_h: float,
    margin: float,
) -> Tuple[shapely.geometry.base.BaseGeometry, Tuple[float, float], Tuple[float, float]]:
    if margin * 2 >= min(plate_w, plate_h):
        raise ValueError("Bed margin is too large for the configured plate size.")
    minx, miny, maxx, maxy = geom.bounds
    width = maxx - minx
    height = maxy - miny
    if width <= EPS or height <= EPS:
        raise ValueError("Degenerate tile encountered.")
    usable_w = plate_w - 2 * margin
    usable_h = plate_h - 2 * margin
    if width > usable_w + EPS or height > usable_h + EPS:
        raise ValueError(
            "Tile does not fit within the specified plate plus margin. "
            "Reduce overlap or margin."
        )
    # Move min corner to origin then shift into the usable area so the part is centered.
    shift_x = margin + (usable_w - width) / 2.0 - minx
    shift_y = margin + (usable_h - height) / 2.0 - miny
    shifted = shapely.affinity.translate(geom, xoff=shift_x, yoff=shift_y)
    return shifted, (minx, miny), (width, height)


def extrude_tile(
    geom: shapely.geometry.base.BaseGeometry,
    thickness: float,
    z_min: float,
    min_area: float,
) -> Optional[trimesh.Trimesh]:
    def _extrude_polygon(poly: shapely.geometry.Polygon) -> Optional[trimesh.Trimesh]:
        if poly.area <= min_area:
            return None
        mesh = trimesh.creation.extrude_polygon(poly, height=thickness)
        mesh.apply_translation([0.0, 0.0, z_min])
        return mesh

    if geom.is_empty:
        return None
    meshes: List[trimesh.Trimesh] = []
    if isinstance(geom, shapely.geometry.Polygon):
        mesh = _extrude_polygon(geom)
        if mesh is not None:
            meshes.append(mesh)
    elif isinstance(geom, shapely.geometry.MultiPolygon):
        for poly in geom.geoms:
            mesh = _extrude_polygon(poly)
            if mesh is not None:
                meshes.append(mesh)
    else:
        raise ValueError(f"Cannot extrude geometry type {geom.geom_type}")
    if not meshes:
        return None
    if len(meshes) == 1:
        return meshes[0]
    return trimesh.util.concatenate(meshes)


def ensure_output_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_output_format(in_file: Path, requested: Optional[str]) -> str:
    if requested:
        return requested.lower()
    suffix = in_file.suffix.lower().lstrip(".")
    if suffix in {"3mf", "stl"}:
        return suffix
    return "3mf"


def main() -> None:
    args = parse_args()
    in_file: Path = args.input_file
    if not in_file.exists():
        raise SystemExit(f"{in_file} does not exist.")
    output_format = resolve_output_format(in_file, args.output_format)
    prefix = args.prefix or in_file.stem
    out_dir = args.output_dir or (in_file.parent / f"{in_file.stem}_plates")
    ensure_output_dir(out_dir)

    mesh = load_mesh(in_file)
    thickness, z_min, z_max = mesh_thickness(mesh)
    footprint = build_footprint(mesh, args.min_normal_z, args.min_feature_area_mm2)
    tiles = generate_tiles(
        footprint=footprint,
        plate_w=args.plate_width_mm,
        plate_h=args.plate_height_mm,
        overlap=args.overlap_mm,
        min_area=args.min_feature_area_mm2,
    )

    manifest = {
        "source": str(in_file),
        "plate_size_mm": [args.plate_width_mm, args.plate_height_mm],
        "overlap_mm": args.overlap_mm,
        "bed_margin_mm": args.bed_margin_mm,
        "output_format": output_format,
        "thickness_mm": thickness,
        "z_range_mm": [z_min, z_max],
        "tiles": [],
    }

    for tile in tiles:
        centered, global_origin, actual_size = recenter_shape(
            tile.shape, args.plate_width_mm, args.plate_height_mm, args.bed_margin_mm
        )
        tile_mesh = extrude_tile(centered, thickness, z_min, args.min_feature_area_mm2)
        tile_name = f"{prefix}_r{tile.row + 1:02d}_c{tile.col + 1:02d}.{output_format}"
        tile_path = out_dir / tile_name
        if tile_mesh is not None and not args.dry_run:
            tile_mesh.export(tile_path)
        manifest["tiles"].append(
            {
                "row": tile.row,
                "col": tile.col,
                "file": tile_name,
                "global_bounds_mm": list(tile.crop_bounds),
                "global_origin_mm": list(global_origin),
                "local_size_mm": list(actual_size),
            }
        )
        print(
            f"Tile (row={tile.row + 1}, col={tile.col + 1}) -> "
            f"{tile_name} {'(dry run)' if args.dry_run else ''}"
        )

    manifest_path = out_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)
    print(f"Manifest written to {manifest_path}")
    if args.dry_run:
        print("Dry-run requested, no meshes were exported.")


if __name__ == "__main__":
    main()

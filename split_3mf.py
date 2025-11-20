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
import shutil
import tempfile
import uuid
import zipfile
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import List, Optional, Tuple
from xml.sax.saxutils import escape

import shapely.affinity
import shapely.geometry
import shapely.ops
import trimesh

# Small floating point tolerance for mm-sized models
EPS = 1e-6
DEFAULT_BAMBU_TEMPLATE = Path(__file__).with_name("bambu_template")
BAMBU_OBJECT_FILENAME = "Cube_1.model"
BAMBU_OBJECT_PATH = f"3D/Objects/{BAMBU_OBJECT_FILENAME}"
BAMBU_REL_PATH = "3D/_rels/3dmodel.model.rels"
BAMBU_DYNAMIC_FILES = {
    "3D/3dmodel.model",
    "Metadata/model_settings.config",
    BAMBU_REL_PATH,
}
BAMBU_DYNAMIC_PREFIXES = ("3D/Objects/",)
BAMBU_REQUIRED_TEMPLATE_FILES = {
    "[Content_Types].xml",
    "_rels/.rels",
    "3D/_rels/3dmodel.model.rels",
    "Metadata/project_settings.config",
    "Metadata/slice_info.config",
}


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
    parser.add_argument(
        "--bambu-template-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing static resources for the Bambu-style 3MF wrapper. "
            "Defaults to the bundled 'bambu_template' folder."
        ),
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


def _format_float(value: float) -> str:
    if abs(value) < 1e-9:
        value = 0.0
    text = f"{value:.6f}"
    text = text.rstrip("0").rstrip(".")
    if text in {"", "-"}:
        text = "0"
    if text == "-0":
        text = "0"
    return text


def _build_bambu_object_xml(mesh: trimesh.Trimesh) -> str:
    vertices = mesh.vertices
    faces = mesh.faces
    vert_lines = [
        f'     <vertex x="{_format_float(x)}" y="{_format_float(y)}" z="{_format_float(z)}"/>'
        for x, y, z in vertices
    ]
    tri_lines = [
        f'     <triangle v1="{int(v1)}" v2="{int(v2)}" v3="{int(v3)}"/>'
        for v1, v2, v3 in faces
    ]
    return "\n".join(
        [
            '<?xml version="1.0" encoding="UTF-8"?>',
            '<model unit="millimeter" xml:lang="en-US" '
            'xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" '
            'xmlns:slic3rpe="http://schemas.slic3r.org/3mf/2017/06" '
            'xmlns:p="http://schemas.microsoft.com/3dmanufacturing/production/2015/06" '
            'requiredextensions="p">',
            ' <metadata name="BambuStudio:3mfVersion">1</metadata>',
            " <resources>",
            '  <object id="1" type="model">',
            "   <mesh>",
            "    <vertices>",
            *vert_lines,
            "    </vertices>",
            "    <triangles>",
            *tri_lines,
            "    </triangles>",
            "   </mesh>",
            "  </object>",
            " </resources>",
            "</model>",
        ]
    )


def _build_bambu_main_model_xml(tile_label: str) -> str:
    today = date.today().isoformat()
    tile_safe = escape(tile_label)
    obj_uuid = uuid.uuid4()
    build_uuid = uuid.uuid4()
    item_uuid = uuid.uuid4()
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<model unit="millimeter" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" xmlns:slic3rpe="http://schemas.slic3r.org/3mf/2017/06" xmlns:p="http://schemas.microsoft.com/3dmanufacturing/production/2015/06" requiredextensions="p">
 <metadata name="Application">BambuStudio-PlateSplitter</metadata>
 <metadata name="BambuStudio:3mfVersion">1</metadata>
 <metadata name="CopyRight"></metadata>
 <metadata name="CreationDate">{today}</metadata>
 <metadata name="Description"></metadata>
 <metadata name="Designer"></metadata>
 <metadata name="DesignerCover"></metadata>
 <metadata name="DesignerUserId"></metadata>
 <metadata name="License"></metadata>
 <metadata name="ModificationDate">{today}</metadata>
 <metadata name="Origin"></metadata>
 <metadata name="Title">{tile_safe}</metadata>
 <resources>
  <object id="2" p:uuid="{obj_uuid}" type="model">
   <components>
    <component p:path="/{BAMBU_OBJECT_PATH}" objectid="1" transform="1 0 0 0 1 0 0 0 1 0 0 0"/>
   </components>
  </object>
 </resources>
 <build p:uuid="{build_uuid}">
  <item objectid="2" p:uuid="{item_uuid}" transform="1 0 0 0 1 0 0 0 1 0 0 0" printable="1"/>
 </build>
</model>
"""


def _build_bambu_model_settings_xml(tile_label: str) -> str:
    safe = escape(tile_label)
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<config>
  <object id="2">
    <metadata key="name" value="{safe}"/>
    <metadata key="extruder" value="1"/>
    <part id="1" subtype="normal_part">
      <metadata key="name" value="{safe}"/>
      <metadata key="matrix" value="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"/>
      <mesh_stat edges_fixed="0" degenerate_facets="0" facets_removed="0" facets_reversed="0" backwards_edges="0"/>
    </part>
  </object>
  <plate>
    <metadata key="plater_id" value="1"/>
    <metadata key="plater_name" value="{safe}"/>
    <metadata key="locked" value="false"/>
    <metadata key="thumbnail_file" value="Metadata/plate_1.png"/>
    <metadata key="top_file" value="Metadata/top_1.png"/>
    <metadata key="pick_file" value="Metadata/pick_1.png"/>
    <model_instance>
      <metadata key="object_id" value="2"/>
      <metadata key="instance_id" value="0"/>
      <metadata key="identify_id" value="0"/>
    </model_instance>
  </plate>
  <assemble>
   <assemble_item object_id="2" instance_id="0" transform="1 0 0 0 1 0 0 0 1 0 0 0" offset="0 0 0" />
  </assemble>
</config>
"""


def _build_bambu_relationships_xml() -> str:
    return """<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
 <Relationship Target="/%(target)s" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>
</Relationships>
""" % {
        "target": BAMBU_OBJECT_PATH
    }


def export_bambu_3mf(
    mesh: trimesh.Trimesh,
    out_path: Path,
    tile_label: str,
    template_dir: Path,
) -> None:
    if not template_dir.exists():
        raise FileNotFoundError(
            f"Bambu template directory '{template_dir}' was not found. "
            "Specify --bambu-template-dir or restore the bundled template."
        )
    template_dir = template_dir.resolve()
    static_files = []
    for path in template_dir.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(template_dir).as_posix()
        if rel in BAMBU_DYNAMIC_FILES or any(rel.startswith(prefix) for prefix in BAMBU_DYNAMIC_PREFIXES):
            continue
        static_files.append((path, rel))
    mesh_xml = _build_bambu_object_xml(mesh)
    main_model_xml = _build_bambu_main_model_xml(tile_label)
    model_settings_xml = _build_bambu_model_settings_xml(tile_label)
    with zipfile.ZipFile(out_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for src, rel in static_files:
            archive.write(src, arcname=rel)
        archive.writestr(BAMBU_OBJECT_PATH, mesh_xml)
        archive.writestr("3D/3dmodel.model", main_model_xml)
        archive.writestr("Metadata/model_settings.config", model_settings_xml)
        archive.writestr(BAMBU_REL_PATH, _build_bambu_relationships_xml())


def _extract_bambu_template_from_source(source: Path) -> Optional[Path]:
    if source.suffix.lower() != ".3mf":
        return None
    try:
        with zipfile.ZipFile(source, "r") as archive:
            names = set(archive.namelist())
            if not BAMBU_REQUIRED_TEMPLATE_FILES.issubset(names):
                return None
            temp_dir = Path(tempfile.mkdtemp(prefix="bambu_template_"))
            archive.extractall(temp_dir)
            return temp_dir
    except (zipfile.BadZipFile, FileNotFoundError):
        return None


def resolve_bambu_template_dir(
    requested: Optional[Path],
    input_file: Path,
) -> Tuple[Optional[Path], Optional[Path], str]:
    """
    Return (template_path, cleanup_path, descriptor).
    If cleanup_path is not None it will be deleted once processing is done.
    """
    if requested is not None:
        return requested, None, str(requested)
    extracted = _extract_bambu_template_from_source(input_file)
    if extracted is not None:
        return extracted, extracted, f"from:{input_file.name}"
    return DEFAULT_BAMBU_TEMPLATE, None, "default"


def main() -> None:
    args = parse_args()
    in_file: Path = args.input_file
    if not in_file.exists():
        raise SystemExit(f"{in_file} does not exist.")
    output_format = resolve_output_format(in_file, args.output_format)
    prefix = args.prefix or in_file.stem
    out_dir = args.output_dir or (in_file.parent / f"{in_file.stem}_plates")
    ensure_output_dir(out_dir)
    bambu_template_dir = args.bambu_template_dir or DEFAULT_BAMBU_TEMPLATE

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
        "bambu_template": None,
        "thickness_mm": thickness,
        "z_range_mm": [z_min, z_max],
        "tiles": [],
    }
    template_cleanup: Optional[Path] = None
    bambu_template_dir: Optional[Path] = None
    if output_format == "3mf":
        bambu_template_dir, template_cleanup, template_descriptor = resolve_bambu_template_dir(
            args.bambu_template_dir, in_file
        )
        manifest["bambu_template"] = template_descriptor

    for tile in tiles:
        centered, global_origin, actual_size = recenter_shape(
            tile.shape, args.plate_width_mm, args.plate_height_mm, args.bed_margin_mm
        )
        tile_mesh = extrude_tile(centered, thickness, z_min, args.min_feature_area_mm2)
        tile_name = f"{prefix}_r{tile.row + 1:02d}_c{tile.col + 1:02d}.{output_format}"
        tile_path = out_dir / tile_name
        if tile_mesh is not None and not args.dry_run:
            if output_format == "stl":
                tile_mesh.export(tile_path)
            else:
                if bambu_template_dir is None:
                    raise SystemExit(
                        "Bambu template directory could not be resolved. "
                        "Pass --bambu-template-dir explicitly."
                    )
                export_bambu_3mf(
                    mesh=tile_mesh,
                    out_path=tile_path,
                    tile_label=tile_name.rsplit(".", 1)[0],
                    template_dir=bambu_template_dir,
                )
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
    if template_cleanup is not None:
        shutil.rmtree(template_cleanup, ignore_errors=True)


if __name__ == "__main__":
    main()

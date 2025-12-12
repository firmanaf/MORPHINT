# -*- coding: utf-8 -*-
"""
***************************************************************************
* MORPHINT: Urban Morphology and Intensity Analyzer                        *
* Parcels: Voronoi / Grid / Momepy tessellation (optional tiling)          *
* Intensity: BCR (%), FAR, OSR (1 - BCR/100), Veg_Prop (0..1 from NDVI)    *
* Buildings: Height, Floors, and 6 morphometrics                           *
*                                                                          *
* Author: Firman Afrianto                                                  *
* PT. Sagamartha Ultima SkillShare 2025                                    *
***************************************************************************
"""

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterNumber,
    QgsProcessingParameterEnum,
    QgsProcessingParameterString,
    QgsProcessingParameterBoolean,
    QgsProcessingContext,
    QgsProcessingFeedback,
    QgsFeature,
    QgsFeatureSink,
    QgsFields,
    QgsField,
    QgsGeometry,
    QgsPointXY,
    QgsRectangle,
    QgsDistanceArea,
    QgsRasterLayer,
    QgsSpatialIndex,
    QgsVectorLayer,
    QgsVectorFileWriter,
    QgsProcessingException,
    QgsProcessingUtils,
    NULL
)
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry, QgsZonalStatistics
from qgis.PyQt.QtCore import QCoreApplication, QVariant

import processing
import math
import os
import tempfile


class MORPHINT(QgsProcessingAlgorithm):

    # Inputs
    INPUT_BUILDINGS = 'INPUT_BUILDINGS'
    INPUT_HEIGHT_RASTER = 'INPUT_HEIGHT_RASTER'
    INPUT_FLOOR_RASTER = 'INPUT_FLOOR_RASTER'
    INPUT_NDVI_RASTER = 'INPUT_NDVI_RASTER'
    INPUT_LIMIT_POLYGON = 'INPUT_LIMIT_POLYGON'
    UNIQUE_ID_FIELD = 'UNIQUE_ID_FIELD'

    # Floor
    FLOOR_CALC_METHOD = 'FLOOR_CALC_METHOD'
    METHOD_FROM_FLOOR_RASTER = 0
    METHOD_FROM_HEIGHT = 1
    FLOOR_METHODS = ['From floor raster (if provided)', 'Derive from height raster']
    FLOOR_HEIGHT_PARAM = 'FLOOR_HEIGHT_PARAM'

    # NDVI
    NDVI_THRESHOLD = 'NDVI_THRESHOLD'

    # Parcelling
    PARCELLING_METHOD = 'PARCELLING_METHOD'
    METHOD_VORONOI = 0
    METHOD_GRID = 1
    METHOD_MOMEPY = 2
    PARCELLING_OPTIONS = [
        'Voronoi by building points',
        'Grid tessellation (square)',
        'Momepy tessellation (by footprints)'
    ]
    GRID_CELL_SIZE = 'GRID_CELL_SIZE'

    # Elongation
    ELONGATION_MODE = 'ELONGATION_MODE'
    ELONG_AXIS_ALIGNED = 0
    ELONG_ORIENTED = 1
    ELONGATION_OPTIONS = ['Axis-aligned bounding box', 'Oriented minimum bounding box']

    # Momepy tiling
    MOMEPY_USE_TILING = 'MOMEPY_USE_TILING'
    TILE_SIZE = 'TILE_SIZE'
    TILE_OVERLAP = 'TILE_OVERLAP'

    # Clip parcels to boundary (Voronoi/Grid; Momepy already uses limit)
    CLIP_PARCELS_TO_LIMIT = 'CLIP_PARCELS_TO_LIMIT'

    # Outputs
    OUTPUT_BUILDINGS = 'OUTPUT_BUILDINGS'
    OUTPUT_PARCELS = 'OUTPUT_PARCELS'

    # Stable parcel id field
    PARC_UID = 'parc_uid'

    # Output fields (buildings)
    FIELD_HEIGHT = 'height_m'
    FIELD_FLOORS = 'floors'
    FIELD_PARCEL_AREA_HA = 'parcel_area_ha'
    FIELD_BLDG_AREA_M2 = 'parcel_bldg_area_m2'
    FIELD_BCR_PCT = 'bcr_pct'
    FIELD_FAR = 'far'
    FIELD_OSR = 'osr'
    FIELD_VEG_PROP = 'veg_prop'

    # Morphometrics
    FIELD_COMPACTNESS = 'compactness'
    FIELD_SHAPE_INDEX = 'shape_index'
    FIELD_FRACTAL_DIM = 'fractal_dim'
    FIELD_SOLIDITY = 'solidity'
    FIELD_CONVEXITY = 'convexity'
    FIELD_ELONGATION = 'elongation'

    def tr(self, s): return QCoreApplication.translate('Processing', s)
    def createInstance(self): return MORPHINT()
    def name(self): return 'morphint_morphology_intensity_analyzer'
    def displayName(self): return self.tr('MORPHINT (Urban Morphology and Intensity Analyzer)')
    def group(self): return self.tr('Urban and Regional Planning Analysis')
    def groupId(self): return 'urbanremotesensing'

    def shortHelpString(self):
        return self.tr(
            "<p><b>Created by Firman Afrianto and Maya Safira</b></p>"
            "<p>A comprehensive QGIS Processing tool to compute <b>parcel-based intensity metrics</b> and "
            "<b>building morphometrics</b> from building footprints, a height raster, an optional floor raster, "
            "and an optional NDVI raster. The tool generates parcels (Voronoi, grid, or Momepy tessellation), "
            "aggregates intensity per parcel, and attaches parcel metrics back to each building for urban planning analysis.</p>"

            "<p><b>Key outputs and how planners use them</b></p>"
            "<ul>"
            "<li><b>height_m</b> (meters): sampled building height. "
            "Useful for <i>urban form profiling</i>, skyline control, and identifying vertical growth hotspots.</li>"
            "<li><b>floors</b> (integer): floor count from raster or derived from height using a standard floor height. "
            "Used as a proxy for <i>development intensity</i> and estimating gross floor area when detailed BIM/cadastre is unavailable.</li>"
            "<li><b>parcel_area_ha</b>: parcel area (hectares) assigned to the building. "
            "Helpful for <i>parcel-level benchmarking</i> and comparing density patterns across neighborhoods.</li>"
            "<li><b>parcel_bldg_area_m2</b>: total building footprint area inside the parcel (m²). "
            "Supports <i>coverage auditing</i> and diagnosing overbuilt or underutilized parcels.</li>"
            "<li><b>bcr_pct</b> (Building Coverage Ratio, 0–100): percent of parcel area covered by buildings. "
            "Commonly used to assess <i>site coverage compliance</i>, imperviousness proxy, and compactness of land consumption.</li>"
            "<li><b>far</b> (Floor Area Ratio, ≥0): total floor area divided by parcel area. "
            "Directly supports <i>zoning intensity</i>, capacity estimation, and scenario testing for infill versus expansion.</li>"
            "<li><b>osr</b> (Open Space Ratio, 0–1): 1 − (BCR/100). "
            "Used to evaluate <i>private open space availability</i> and balance between built and unbuilt land on each parcel.</li>"
            "<li><b>veg_prop</b> (0–1, optional): mean of a binary NDVI mask (NDVI ≥ threshold) within each parcel. "
            "Supports <i>green coverage screening</i>, microclimate co-benefit assessment, and identifying parcels needing greening interventions.</li>"
            "</ul>"

            "<p><b>Building morphometrics (shape indices)</b></p>"
            "<ul>"
            "<li><b>compactness</b>: (4πA)/(P²), 0–1. Higher values indicate more compact footprints. "
            "Useful for diagnosing <i>fragmented or irregular building forms</i> often linked to inefficient land use.</li>"
            "<li><b>shape_index</b>: P / (2√(πA)). Higher values indicate more irregular shapes. "
            "Useful for identifying <i>complex building envelopes</i> that may correlate with energy, ventilation, or plot constraints.</li>"
            "<li><b>fractal_dim</b>: 2 × ln(P) / ln(A). Higher values imply more boundary complexity. "
            "Useful for urban morphology research and detecting <i>fine-grain, jagged patterns</i>.</li>"
            "<li><b>solidity</b>: A / A<sub>convex</sub>, 0–1. Lower values suggest concavity or internal voids. "
            "Useful for recognizing <i>courtyard-like forms</i> or complex footprints.</li>"
            "<li><b>convexity</b>: P / P<sub>convex</sub>, 0–1 (values closer to 1 are more convex). "
            "Useful for distinguishing <i>simple vs. indented footprints</i> at scale.</li>"
            "<li><b>elongation</b>: max(bbox side) / min(bbox side), ≥1. "
            "Used to detect <i>linear/slab-type buildings</i> and corridor-oriented forms (supports axis-aligned or oriented MBB).</li>"
            "</ul>"

            "<p><b>Intensity metrics (parcel level)</b></p>"
            "<ul>"
            "<li><b>BCR (bcr_pct)</b>: 0–100 percent of parcel area covered by buildings.</li>"
            "<li><b>FAR (far)</b>: total floor area divided by parcel area, ratio ≥ 0.</li>"
            "<li><b>OSR (osr)</b>: 1 − (BCR/100), 0–1.</li>"
            "<li><b>Vegetation proportion (veg_prop)</b>: fraction of parcel pixels with NDVI ≥ threshold (optional).</li>"
            "</ul>"

            "<p><b>Parcelization methods</b></p>"
            "<ul>"
            "<li><b>Voronoi</b>: partitions space based on nearest building point; suitable when no parcel/cadastre data exists.</li>"
            "<li><b>Grid</b>: square tessellation clipped to the limit polygon (if provided) or analysis extent; fast for citywide screening.</li>"
            "<li><b>Momepy tessellation</b>: footprint-driven cells (optional tiling for large datasets); best for detailed urban morphology.</li>"
            "</ul>"

            "<p><b>Practical planning applications</b></p>"
            "<ul>"
            "<li><b>Zoning and capacity estimation</b>: use FAR and floors to approximate development capacity and intensity gradients.</li>"
            "<li><b>Compliance screening</b>: compare BCR and FAR against local controls (coverage, FAR caps) for rapid auditing.</li>"
            "<li><b>Infill vs. sprawl diagnostics</b>: identify underbuilt parcels (low BCR/FAR) versus overbuilt parcels (high BCR/FAR).</li>"
            "<li><b>Urban form and design review</b>: use morphometrics to map typologies (compact blocks, elongated slabs, irregular footprints).</li>"
            "<li><b>Greening prioritization</b>: use OSR and veg_prop to target parcels with low open space and low vegetation.</li>"
            "</ul>"

            "<p><b>Notes and constraints</b></p>"
            "<ul>"
            "<li><b>Projected CRS (meters)</b> is strongly recommended for accurate area and perimeter measurements.</li>"
            "<li>If a floor raster is missing, floors are derived from height using <i>standard floor height</i>.</li>"
            "<li>NDVI may be stored as scaled integers (e.g., 0–10000). If so, adjust the threshold accordingly.</li>"
            "</ul>"
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_BUILDINGS, self.tr('Building polygons'), [QgsProcessing.TypeVectorPolygon]
        ))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.INPUT_HEIGHT_RASTER, self.tr('Height raster (meters)')
        ))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.INPUT_FLOOR_RASTER, self.tr('Floor count raster (optional)'), optional=True
        ))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.INPUT_NDVI_RASTER, self.tr('NDVI raster (optional)'), optional=True
        ))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_LIMIT_POLYGON, self.tr('Limit polygon (optional)'),
            [QgsProcessing.TypeVectorPolygon], optional=True
        ))
        self.addParameter(QgsProcessingParameterBoolean(
            self.CLIP_PARCELS_TO_LIMIT, self.tr('Clip parcels to limit polygon (if provided)'),
            defaultValue=True
        ))
        self.addParameter(QgsProcessingParameterString(
            self.UNIQUE_ID_FIELD, self.tr('Unique ID field in buildings for Momepy (blank = auto)'),
            defaultValue='', optional=True
        ))

        self.addParameter(QgsProcessingParameterEnum(
            self.FLOOR_CALC_METHOD, self.tr('Floor calculation method'),
            options=self.FLOOR_METHODS, defaultValue=self.METHOD_FROM_HEIGHT
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.FLOOR_HEIGHT_PARAM, self.tr('Standard floor height (m)'),
            QgsProcessingParameterNumber.Double, defaultValue=3.0, minValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.NDVI_THRESHOLD, self.tr('NDVI vegetation threshold'),
            QgsProcessingParameterNumber.Double, defaultValue=0.2, minValue=-1.0, maxValue=1.0
        ))

        self.addParameter(QgsProcessingParameterEnum(
            self.PARCELLING_METHOD, self.tr('Parcelization method'),
            options=self.PARCELLING_OPTIONS, defaultValue=self.METHOD_VORONOI
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.GRID_CELL_SIZE, self.tr('Grid cell size (CRS units)'),
            QgsProcessingParameterNumber.Double, defaultValue=50.0, minValue=1.0
        ))

        self.addParameter(QgsProcessingParameterEnum(
            self.ELONGATION_MODE, self.tr('Elongation method'),
            options=self.ELONGATION_OPTIONS, defaultValue=self.ELONG_AXIS_ALIGNED
        ))

        self.addParameter(QgsProcessingParameterBoolean(
            self.MOMEPY_USE_TILING, self.tr('Use tiling for Momepy (faster for big data)'), defaultValue=True
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.TILE_SIZE, self.tr('Tile size (meters)'), QgsProcessingParameterNumber.Double,
            defaultValue=2000.0, minValue=100.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.TILE_OVERLAP, self.tr('Tile overlap (meters)'), QgsProcessingParameterNumber.Double,
            defaultValue=200.0, minValue=0.0
        ))

        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_BUILDINGS, self.tr('Buildings with metrics')
        ))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_PARCELS, self.tr('Parcels with metrics')
        ))

    # ---------------- helpers ----------------
    @staticmethod
    def _unique_field_name(fields: QgsFields, base: str) -> str:
        name = base
        i = 2
        while fields.indexOf(name) != -1:
            name = f"{base}_{i}"
            i += 1
        return name

    @staticmethod
    def _find_field_case_insensitive(layer: QgsVectorLayer, wanted: str):
        w = wanted.lower()
        for f in layer.fields():
            if f.name().lower() == w:
                return f.name()
        return None

    @staticmethod
    def _as_layer(any_layer, context: QgsProcessingContext):
        if hasattr(any_layer, 'isValid'):
            return any_layer
        lyr = QgsProcessingUtils.mapLayerFromString(str(any_layer), context)
        if lyr is not None:
            return lyr
        return QgsVectorLayer(str(any_layer), "tmp", "ogr")

    @staticmethod
    def _ensure_parc_uid(parcels_layer: QgsVectorLayer, feedback: QgsProcessingFeedback) -> str:
        field_name = MORPHINT.PARC_UID
        if parcels_layer.fields().indexOf(field_name) != -1:
            return field_name

        parcels_layer.startEditing()
        parcels_layer.dataProvider().addAttributes([QgsField(field_name, QVariant.Int)])
        parcels_layer.updateFields()
        idx = parcels_layer.fields().indexOf(field_name)

        v = 1
        for f in parcels_layer.getFeatures():
            parcels_layer.changeAttributeValue(f.id(), idx, v)
            v += 1
        parcels_layer.commitChanges()
        feedback.pushInfo(f"Added parcel UID field: {field_name}")
        return field_name

    @staticmethod
    def _mbb_oriented_dims(geom: QgsGeometry, da: QgsDistanceArea):
        if not geom or geom.isNull():
            return None, None
        try:
            obb = geom.orientedMinimumBoundingBox()
            if not obb or obb.isNull():
                raise Exception("oriented MBB is null")
            ring = obb.asPolygon()
            if not ring or len(ring[0]) < 4:
                raise Exception("oriented MBB ring invalid")
            p0, p1, p2 = ring[0][0], ring[0][1], ring[0][2]
            w = da.measureLine(QgsPointXY(p0.x(), p0.y()), QgsPointXY(p1.x(), p1.y()))
            h = da.measureLine(QgsPointXY(p1.x(), p1.y()), QgsPointXY(p2.x(), p2.y()))
            if h > w:
                w, h = h, w
            return w, h
        except Exception:
            bb = geom.boundingBox()
            return bb.width(), bb.height()

    def _momepy_tessellation(self, bldg_layer, limit_layer, unique_id_name,
                             tile_size, overlap, use_tiling,
                             context: QgsProcessingContext, feedback: QgsProcessingFeedback):
        """
        Build tessellation using momepy.
        If use_tiling=True, process per tile with overlap and merge.
        Returns QgsVectorLayer of tessellation polygons.
        """
        try:
            import geopandas as gpd
            import momepy
        except Exception as e:
            raise QgsProcessingException(
                "Momepy mode requires geopandas/shapely/momepy installed in QGIS Python. "
                f"Import error: {e}"
            )

        tmp_bldg = QgsProcessingUtils.generateTempFilename('bldg_tmp.gpkg')
        opts_b = QgsVectorFileWriter.SaveVectorOptions()
        opts_b.driverName = 'GPKG'
        opts_b.fileEncoding = 'UTF-8'
        res_b, err_b = QgsVectorFileWriter.writeAsVectorFormatV2(
            bldg_layer, tmp_bldg, context.transformContext(), opts_b
        )
        if res_b != QgsVectorFileWriter.NoError:
            raise QgsProcessingException(f"Failed to export buildings to GPKG (code={res_b}): {err_b}")

        tmp_limit = None
        if limit_layer and limit_layer.isValid():
            tmp_limit = QgsProcessingUtils.generateTempFilename('limit_tmp.gpkg')
            opts_l = QgsVectorFileWriter.SaveVectorOptions()
            opts_l.driverName = 'GPKG'
            opts_l.fileEncoding = 'UTF-8'
            res_l, err_l = QgsVectorFileWriter.writeAsVectorFormatV2(
                limit_layer, tmp_limit, context.transformContext(), opts_l
            )
            if res_l != QgsVectorFileWriter.NoError:
                raise QgsProcessingException(f"Failed to export limit polygon to GPKG (code={res_l}): {err_l}")

        gdf = gpd.read_file(tmp_bldg)
        if gdf is None or gdf.empty:
            raise QgsProcessingException("Buildings GeoDataFrame is empty (Momepy).")

        # Determine unique id field
        if unique_id_name and unique_id_name in gdf.columns:
            uid_field = unique_id_name
        else:
            uid_field = '_uid_'
            gdf[uid_field] = range(1, len(gdf) + 1)

        limit_gdf = None
        if tmp_limit:
            limit_gdf = gpd.read_file(tmp_limit)
            if limit_gdf is not None and limit_gdf.empty:
                limit_gdf = None

        if not use_tiling:
            feedback.pushInfo(f"Momepy tessellation (no tiling) on {len(gdf)} buildings.")
            tess = momepy.Tessellation(gdf, unique_id=uid_field, limit=limit_gdf)
            cells = tess.tessellation
        else:
            # Use bounds from limit if available else buildings
            if limit_gdf is not None:
                xmin, ymin, xmax, ymax = limit_gdf.total_bounds
            else:
                xmin, ymin, xmax, ymax = gdf.total_bounds

            ts = float(tile_size)
            ov = float(overlap)

            tiles = []
            xi = xmin
            while xi < xmax:
                yi = ymin
                while yi < ymax:
                    tiles.append((xi - ov, yi - ov, xi + ts + ov, yi + ts + ov))
                    yi += ts
                xi += ts

            feedback.pushInfo(f"Momepy tiling: {len(tiles)} tiles, tile={ts}m, overlap={ov}m")

            import pandas as pd
            out_list = []

            for i, (tx0, ty0, tx1, ty1) in enumerate(tiles, 1):
                sub = gdf.cx[tx0:tx1, ty0:ty1]
                if sub is None or sub.empty:
                    continue

                lim_tile = None
                if limit_gdf is not None:
                    lim_sub = limit_gdf.cx[tx0:tx1, ty0:ty1]
                    if lim_sub is not None and not lim_sub.empty:
                        lim_tile = lim_sub

                try:
                    tess = momepy.Tessellation(sub, unique_id=uid_field, limit=lim_tile)
                    cells_i = tess.tessellation

                    # Keep only core area (reduce seam artifacts)
                    core_x0, core_y0 = tx0 + ov, ty0 + ov
                    core_x1, core_y1 = tx1 - ov, ty1 - ov
                    if core_x1 > core_x0 and core_y1 > core_y0:
                        cells_i = cells_i.clip(mask=(core_x0, core_y0, core_x1, core_y1))

                    out_list.append(cells_i)
                except Exception as e:
                    feedback.pushWarning(f"Momepy tessellation failed on tile {i}: {e}")

            if not out_list:
                raise QgsProcessingException("No tessellation produced by tiling. Try disabling tiling or increasing tile size.")

            cells = pd.concat(out_list, ignore_index=True)

        tmp_tess = QgsProcessingUtils.generateTempFilename('momepy_tess.gpkg')
        cells.to_file(tmp_tess, driver='GPKG')

        parcels_layer = QgsVectorLayer(tmp_tess, "parcels_momepy", "ogr")
        if not parcels_layer or not parcels_layer.isValid():
            raise QgsProcessingException("Failed to load Momepy tessellation layer.")
        return parcels_layer

    def _compute_ndvi_veg_prop(self, parcels_layer: QgsVectorLayer, ndvi_raster: QgsRasterLayer,
                               threshold: float, context, feedback):
        """
        QGIS 3.40 safe:
        - Build binary raster with QgsRasterCalculator
        - No reliance on QgsRasterCalculator.Ok constant
        - Compute zonal mean of binary raster -> vegetation proportion
        """
        try:
            out_fd, out_path = tempfile.mkstemp(suffix='.tif')
            os.close(out_fd)

            ndvi_ref = QgsRasterCalculatorEntry()
            ndvi_ref.ref = 'ndvi@1'
            ndvi_ref.raster = ndvi_raster
            ndvi_ref.bandNumber = 1
            entries = [ndvi_ref]

            ndvi_provider = ndvi_raster.dataProvider()
            ndvi_nodata = ndvi_provider.sourceNoDataValue(1)
            has_nodata = ndvi_nodata is not None

            expr = f"({ndvi_ref.ref} >= {threshold}) * 1"
            if has_nodata:
                expr = f"(({ndvi_ref.ref} >= {threshold}) * 1) * ({ndvi_ref.ref} != {ndvi_nodata})"

            calc = QgsRasterCalculator(
                expr, out_path, 'GTiff',
                ndvi_raster.extent(), ndvi_raster.crs(),
                ndvi_raster.width(), ndvi_raster.height(),
                entries
            )
            try:
                calc.setNoDataValue(-9999)
            except Exception:
                pass

            _ = calc.processCalculation(feedback)

            ndvi_bin = QgsRasterLayer(out_path, "ndvi_bin", "gdal")
            if not ndvi_bin or not ndvi_bin.isValid():
                feedback.pushWarning("NDVI binary raster is invalid. Veg_Prop will be NULL.")
                return None

            zs = QgsZonalStatistics(parcels_layer, ndvi_bin, 'veg_', 1, QgsZonalStatistics.Mean)
            ret = zs.calculateStatistics(feedback)
            if ret != 0:
                feedback.pushWarning(f"Zonal statistics failed (code={ret}). Veg_Prop will be NULL.")
                return None

            return self._find_field_case_insensitive(parcels_layer, 'veg_mean')

        except Exception as e:
            feedback.pushWarning(f"NDVI vegetation proportion failed: {e}")
            return None

    # ---------------- main ----------------
    def processAlgorithm(self, parameters, context: QgsProcessingContext, feedback: QgsProcessingFeedback):

        bldg_layer = self.parameterAsVectorLayer(parameters, self.INPUT_BUILDINGS, context)
        if not bldg_layer or not bldg_layer.isValid():
            raise QgsProcessingException("Invalid building layer.")

        height_raster = self.parameterAsRasterLayer(parameters, self.INPUT_HEIGHT_RASTER, context)
        if not height_raster or not height_raster.isValid():
            raise QgsProcessingException("Invalid height raster.")

        floor_raster = self.parameterAsRasterLayer(parameters, self.INPUT_FLOOR_RASTER, context)
        ndvi_raster = self.parameterAsRasterLayer(parameters, self.INPUT_NDVI_RASTER, context)

        limit_layer = self.parameterAsVectorLayer(parameters, self.INPUT_LIMIT_POLYGON, context)
        clip_to_limit = self.parameterAsBool(parameters, self.CLIP_PARCELS_TO_LIMIT, context)
        unique_id_name = self.parameterAsString(parameters, self.UNIQUE_ID_FIELD, context)

        floor_method = self.parameterAsEnum(parameters, self.FLOOR_CALC_METHOD, context)
        std_floor_h = self.parameterAsDouble(parameters, self.FLOOR_HEIGHT_PARAM, context)
        ndvi_thr = self.parameterAsDouble(parameters, self.NDVI_THRESHOLD, context)

        parc_method = self.parameterAsEnum(parameters, self.PARCELLING_METHOD, context)
        grid_size = self.parameterAsDouble(parameters, self.GRID_CELL_SIZE, context)

        elong_mode = self.parameterAsEnum(parameters, self.ELONGATION_MODE, context)

        use_tiling = self.parameterAsBool(parameters, self.MOMEPY_USE_TILING, context)
        tile_size = self.parameterAsDouble(parameters, self.TILE_SIZE, context)
        overlap = self.parameterAsDouble(parameters, self.TILE_OVERLAP, context)

        crs = bldg_layer.sourceCrs()
        if not crs.isValid():
            raise QgsProcessingException("Building CRS invalid.")
        if crs.isGeographic():
            feedback.pushWarning("Geographic CRS detected. Use projected CRS (meters) for reliable results.")

        da = QgsDistanceArea()
        da.setSourceCrs(crs, context.transformContext())
        da.setEllipsoid(crs.ellipsoidAcronym() or 'WGS84')

        # Floors fallback
        if floor_method == self.METHOD_FROM_FLOOR_RASTER and not floor_raster:
            feedback.pushInfo("Floor raster not provided. Falling back to derive floors from height raster.")
            floor_method = self.METHOD_FROM_HEIGHT

        if floor_method == self.METHOD_FROM_HEIGHT and std_floor_h <= 0:
            raise QgsProcessingException("Standard floor height must be positive.")

        # 1) Sample height and floors at building points
        feedback.pushInfo("Sampling height and floors at building points.")
        hprov = height_raster.dataProvider()
        fprov = floor_raster.dataProvider() if floor_raster else None

        building_metrics = {}
        processed_fids = []

        for feat in bldg_layer.getFeatures():
            geom = feat.geometry()
            if not geom or geom.isNull() or not geom.isGeosValid():
                continue

            p = geom.pointOnSurface().asPoint()
            pt = QgsPointXY(p.x(), p.y())

            h_val, ok_h = hprov.sample(p, 1)
            height = None
            if ok_h and h_val is not None:
                try:
                    height = float(h_val)
                    if height < 0:
                        height = 0.0
                except Exception:
                    height = None

            floors = None
            if floor_method == self.METHOD_FROM_FLOOR_RASTER and fprov:
                f_val, ok_f = fprov.sample(p, 1)
                if ok_f and f_val is not None:
                    try:
                        floors = max(1, int(math.ceil(float(f_val))))
                    except Exception:
                        floors = None
            elif floor_method == self.METHOD_FROM_HEIGHT and height is not None and height > 0:
                floors = max(1, int(math.ceil(height / std_floor_h)))

            building_metrics[feat.id()] = {'pt': pt, 'height': height, 'floors': floors}
            processed_fids.append(feat.id())

        if not building_metrics:
            raise QgsProcessingException("No valid buildings after sampling.")
        feedback.pushInfo(f"Processed buildings: {len(building_metrics)}")

        # 2) Create parcels
        parcels_layer = None

        if parc_method == self.METHOD_VORONOI:
            feedback.pushInfo("Creating Voronoi parcels from building points.")
            pts = QgsVectorLayer(f"Point?crs={crs.toWkt()}", "bldg_pts", "memory")
            pr = pts.dataProvider()
            pr.addAttributes([QgsField('bldg_fid', QVariant.Int)])
            pts.updateFields()

            pts.startEditing()
            for fid in processed_fids:
                f = QgsFeature(pts.fields())
                f.setGeometry(QgsGeometry.fromPointXY(building_metrics[fid]['pt']))
                f['bldg_fid'] = int(fid)
                pr.addFeature(f)
            pts.commitChanges()

            extent = limit_layer.extent() if (limit_layer and limit_layer.isValid()) else bldg_layer.extent()
            buffer_dist = max(10.0, (extent.width() + extent.height()) * 0.1) if extent and extent.width() > 0 else 100.0

            res_v = processing.run(
                'native:voronoipolygons',
                {'INPUT': pts, 'BUFFER': buffer_dist, 'OUTPUT': 'TEMPORARY_OUTPUT'},
                context=context, feedback=feedback, is_child_algorithm=True
            )
            parcels_layer = self._as_layer(res_v['OUTPUT'], context)

            if clip_to_limit and limit_layer and limit_layer.isValid():
                feedback.pushInfo("Clipping parcels to limit polygon.")
                res_clip = processing.run(
                    'native:clip',
                    {'INPUT': parcels_layer, 'OVERLAY': limit_layer, 'OUTPUT': 'TEMPORARY_OUTPUT'},
                    context=context, feedback=feedback, is_child_algorithm=True
                )
                parcels_layer = self._as_layer(res_clip['OUTPUT'], context)

            parcels_layer.setName("parcels_voronoi")

        elif parc_method == self.METHOD_GRID:
            feedback.pushInfo("Creating grid parcels.")
            ext = limit_layer.extent() if (limit_layer and limit_layer.isValid()) else bldg_layer.extent()

            res_g = processing.run(
                'native:creategrid',
                {
                    'TYPE': 2,
                    'EXTENT': ext,
                    'HSPACING': grid_size,
                    'VSPACING': grid_size,
                    'HOVERLAY': 0.0,
                    'VOVERLAY': 0.0,
                    'CRS': crs,
                    'OUTPUT': 'TEMPORARY_OUTPUT'
                },
                context=context, feedback=feedback, is_child_algorithm=True
            )
            parcels_layer = self._as_layer(res_g['OUTPUT'], context)

            if clip_to_limit and limit_layer and limit_layer.isValid():
                feedback.pushInfo("Clipping parcels to limit polygon.")
                res_clip = processing.run(
                    'native:clip',
                    {'INPUT': parcels_layer, 'OVERLAY': limit_layer, 'OUTPUT': 'TEMPORARY_OUTPUT'},
                    context=context, feedback=feedback, is_child_algorithm=True
                )
                parcels_layer = self._as_layer(res_clip['OUTPUT'], context)

            parcels_layer.setName("parcels_grid")

        else:
            feedback.pushInfo("Creating Momepy tessellation (tiling optional).")
            parcels_layer = self._momepy_tessellation(
                bldg_layer=bldg_layer,
                limit_layer=limit_layer,
                unique_id_name=unique_id_name,
                tile_size=tile_size,
                overlap=overlap,
                use_tiling=use_tiling,
                context=context,
                feedback=feedback
            )
            parcels_layer.setName("parcels_momepy")

        if not parcels_layer or not parcels_layer.isValid():
            raise QgsProcessingException("Failed to create parcels layer.")

        parc_uid_field = self._ensure_parc_uid(parcels_layer, feedback)

        # Spatial indexes (best effort)
        for lyr in (bldg_layer, parcels_layer):
            try:
                processing.run('native:createspatialindex', {'INPUT': lyr},
                               context=context, feedback=feedback, is_child_algorithm=True)
            except Exception:
                pass

        # 3) NDVI vegetation proportion (optional)
        veg_mean_field = None
        if ndvi_raster and ndvi_raster.isValid():
            feedback.pushInfo("Computing vegetation proportion from NDVI (binary mean).")
            veg_mean_field = self._compute_ndvi_veg_prop(parcels_layer, ndvi_raster, ndvi_thr, context, feedback)

        # 4) Parcel intensity aggregation by intersection with buildings
        feedback.pushInfo("Computing parcel intensity (BCR, FAR, OSR) by intersection aggregation.")
        b_index = QgsSpatialIndex(bldg_layer.getFeatures())
        parcel_metrics = {}  # keyed by parc_uid

        for p in parcels_layer.getFeatures():
            pg = p.geometry()
            if not pg or pg.isNull():
                continue

            area_m2 = da.measureArea(pg)
            if area_m2 <= 0:
                continue

            parc_uid = p[parc_uid_field]
            if parc_uid in (None, NULL):
                continue

            total_b_area = 0.0
            total_floor_area = 0.0

            cand = b_index.intersects(pg.boundingBox())
            for bf_id in cand:
                if bf_id not in building_metrics:
                    continue
                bf = bldg_layer.getFeature(bf_id)
                bg = bf.geometry()
                if not bg or bg.isNull():
                    continue
                if not pg.intersects(bg):
                    continue
                try:
                    inter = pg.intersection(bg)
                    if inter and not inter.isNull():
                        a = da.measureArea(inter)
                        if a > 0:
                            total_b_area += a
                            fl = building_metrics[bf_id].get('floors')
                            if fl and fl > 0:
                                total_floor_area += a * fl
                except Exception as e:
                    feedback.pushWarning(f"Intersection failed (parcel {parc_uid} vs bldg {bf_id}): {e}")

            bcr = (total_b_area / area_m2) * 100.0
            if bcr < 0:
                bcr = 0.0
            if bcr > 100:
                bcr = 100.0

            far = (total_floor_area / area_m2)
            osr = 1.0 - (bcr / 100.0)

            veg_prop = None
            if veg_mean_field:
                try:
                    veg_prop = float(p[veg_mean_field]) if p[veg_mean_field] not in (None, NULL) else None
                except Exception:
                    veg_prop = None

            parcel_metrics[parc_uid] = {
                'area_ha': area_m2 / 10000.0,
                'b_area_m2': total_b_area,
                'bcr_pct': round(bcr, 2),
                'far': round(far, 2),
                'osr': round(osr, 4),
                'veg_prop': round(veg_prop, 4) if veg_prop is not None else None
            }

        # 5) Output buildings fields
        feedback.pushInfo("Writing output buildings layer.")
        out_fields = QgsFields(bldg_layer.fields())

        f_height = self._unique_field_name(out_fields, self.FIELD_HEIGHT)
        f_floors = self._unique_field_name(out_fields, self.FIELD_FLOORS)
        f_area = self._unique_field_name(out_fields, self.FIELD_PARCEL_AREA_HA)
        f_ba = self._unique_field_name(out_fields, self.FIELD_BLDG_AREA_M2)
        f_bcr = self._unique_field_name(out_fields, self.FIELD_BCR_PCT)
        f_far = self._unique_field_name(out_fields, self.FIELD_FAR)
        f_osr = self._unique_field_name(out_fields, self.FIELD_OSR)
        f_veg = self._unique_field_name(out_fields, self.FIELD_VEG_PROP)

        f_comp = self._unique_field_name(out_fields, self.FIELD_COMPACTNESS)
        f_si = self._unique_field_name(out_fields, self.FIELD_SHAPE_INDEX)
        f_fd = self._unique_field_name(out_fields, self.FIELD_FRACTAL_DIM)
        f_sol = self._unique_field_name(out_fields, self.FIELD_SOLIDITY)
        f_conv = self._unique_field_name(out_fields, self.FIELD_CONVEXITY)
        f_elong = self._unique_field_name(out_fields, self.FIELD_ELONGATION)

        out_fields.append(QgsField(f_height, QVariant.Double, 'double', 12, 2))
        out_fields.append(QgsField(f_floors, QVariant.Int))
        out_fields.append(QgsField(f_area, QVariant.Double, 'double', 15, 4))
        out_fields.append(QgsField(f_ba, QVariant.Double, 'double', 15, 2))
        out_fields.append(QgsField(f_bcr, QVariant.Double, 'double', 7, 2))
        out_fields.append(QgsField(f_far, QVariant.Double, 'double', 10, 2))
        out_fields.append(QgsField(f_osr, QVariant.Double, 'double', 8, 4))
        out_fields.append(QgsField(f_veg, QVariant.Double, 'double', 8, 4))

        out_fields.append(QgsField(f_comp, QVariant.Double, 'double', 8, 4))
        out_fields.append(QgsField(f_si, QVariant.Double, 'double', 8, 4))
        out_fields.append(QgsField(f_fd, QVariant.Double, 'double', 8, 4))
        out_fields.append(QgsField(f_sol, QVariant.Double, 'double', 8, 4))
        out_fields.append(QgsField(f_conv, QVariant.Double, 'double', 8, 4))
        out_fields.append(QgsField(f_elong, QVariant.Double, 'double', 8, 4))

        sink_b, dest_b = self.parameterAsSink(
            parameters, self.OUTPUT_BUILDINGS, context,
            out_fields, bldg_layer.wkbType(), crs
        )
        if sink_b is None:
            raise QgsProcessingException(f"Could not create buildings output sink: {dest_b}")

        parc_index = QgsSpatialIndex(parcels_layer.getFeatures())

        for fid in processed_fids:
            src = bldg_layer.getFeature(fid)
            bm = building_metrics[fid]
            pt = bm['pt']

            found = None
            ids = parc_index.intersects(QgsRectangle(pt, pt))
            for pid in ids:
                pf = parcels_layer.getFeature(pid)
                pg = pf.geometry()
                if pg and not pg.isNull() and pg.contains(QgsGeometry.fromPointXY(pt)):
                    parc_uid = pf[parc_uid_field]
                    found = parcel_metrics.get(parc_uid)
                    break

            height_val = bm.get('height')
            floors_val = bm.get('floors')

            g = src.geometry()
            try:
                A = da.measureArea(g)
            except Exception:
                A = None
            try:
                P = da.measurePerimeter(g)
            except Exception:
                P = g.length() if g else None

            try:
                ch = g.convexHull() if g else None
                A_ch = da.measureArea(ch) if ch and not ch.isNull() else None
                P_ch = da.measurePerimeter(ch) if ch and not ch.isNull() else None
            except Exception:
                A_ch, P_ch = None, None

            if elong_mode == self.ELONG_ORIENTED:
                w, h = self._mbb_oriented_dims(g, da)
            else:
                bb = g.boundingBox() if g else None
                w = bb.width() if bb else None
                h = bb.height() if bb else None

            compactness = None
            if A and A > 0 and P and P > 0:
                compactness = (4.0 * math.pi * A) / (P * P)
                if compactness > 1.0:
                    compactness = 1.0
                compactness = round(compactness, 4)

            shape_index = None
            if A and A > 0 and P and P > 0:
                denom = 2.0 * math.sqrt(math.pi * A)
                if denom > 0:
                    shape_index = round(P / denom, 4)

            fractal_dim = None
            if A and A > 1.0 and P and P > 1.0:
                try:
                    fractal_dim = round(2.0 * (math.log(P) / math.log(A)), 4)
                except Exception:
                    fractal_dim = None

            solidity = None
            if A and A > 0 and A_ch and A_ch > 0:
                val = A / A_ch
                if 0.0 <= val <= 1.0:
                    solidity = round(val, 4)

            convexity = None
            if P and P > 0 and P_ch and P_ch > 0:
                convexity = round(P / P_ch, 4)

            elongation = None
            if w is not None and h is not None:
                mn = min(w, h)
                mx = max(w, h)
                if mn and mn > 0:
                    elongation = round(mx / mn, 4)

            base_attrs = list(src.attributes())
            if found is None:
                extra = [
                    height_val, floors_val,
                    None, None, None, None, None, None,
                    compactness, shape_index, fractal_dim, solidity, convexity, elongation
                ]
            else:
                extra = [
                    height_val, floors_val,
                    found.get('area_ha'),
                    found.get('b_area_m2'),
                    found.get('bcr_pct'),
                    found.get('far'),
                    found.get('osr'),
                    found.get('veg_prop'),
                    compactness, shape_index, fractal_dim, solidity, convexity, elongation
                ]

            nf = QgsFeature()
            nf.setGeometry(src.geometry())
            nf.setFields(out_fields)
            nf.setAttributes(base_attrs + extra)
            sink_b.addFeature(nf, QgsFeatureSink.FastInsert)

        del sink_b

        # Output parcels
        feedback.pushInfo("Writing output parcels layer.")
        p_fields = QgsFields(parcels_layer.fields())

        p_area = self._unique_field_name(p_fields, 'area_ha')
        p_ba = self._unique_field_name(p_fields, 'bldg_area_m2')
        p_bcr = self._unique_field_name(p_fields, 'bcr_pct')
        p_far = self._unique_field_name(p_fields, 'far')
        p_osr = self._unique_field_name(p_fields, 'osr')
        p_veg = self._unique_field_name(p_fields, 'veg_prop')

        p_fields.append(QgsField(p_area, QVariant.Double, 'double', 15, 4))
        p_fields.append(QgsField(p_ba, QVariant.Double, 'double', 15, 2))
        p_fields.append(QgsField(p_bcr, QVariant.Double, 'double', 7, 2))
        p_fields.append(QgsField(p_far, QVariant.Double, 'double', 10, 2))
        p_fields.append(QgsField(p_osr, QVariant.Double, 'double', 8, 4))
        p_fields.append(QgsField(p_veg, QVariant.Double, 'double', 8, 4))

        sink_p, dest_p = self.parameterAsSink(
            parameters, self.OUTPUT_PARCELS, context,
            p_fields, parcels_layer.wkbType(), crs
        )
        if sink_p is None:
            raise QgsProcessingException(f"Could not create parcels output sink: {dest_p}")

        for pf in parcels_layer.getFeatures():
            parc_uid = pf[parc_uid_field]
            m = parcel_metrics.get(parc_uid, None)

            outp = QgsFeature()
            outp.setGeometry(pf.geometry())
            outp.setFields(p_fields)

            attrs = list(pf.attributes())
            if m is None:
                attrs.extend([None, None, None, None, None, None])
            else:
                attrs.extend([
                    m.get('area_ha'),
                    m.get('b_area_m2'),
                    m.get('bcr_pct'),
                    m.get('far'),
                    m.get('osr'),
                    m.get('veg_prop')
                ])

            outp.setAttributes(attrs)
            sink_p.addFeature(outp, QgsFeatureSink.FastInsert)

        del sink_p

        feedback.pushInfo("Done.")
        return {self.OUTPUT_BUILDINGS: dest_b, self.OUTPUT_PARCELS: dest_p}

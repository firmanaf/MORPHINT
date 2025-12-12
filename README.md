# MORPHINT  
## Urban Morphology and Intensity Analyzer for QGIS

**MORPHINT** is a QGIS Processing plugin for analyzing **urban morphology and development intensity** using building footprints, height or floor information, and optional vegetation indicators.  
The plugin integrates **parcel-based density metrics** with **building morphometrics**, making it especially useful in urban areas where cadastral parcels are incomplete or unavailable.

**Authors**  
Firman Afrianto, Maya Safira  

---

## 🔍 Key Capabilities

MORPHINT provides an integrated workflow to compute:

### 1. Parcel-Based Intensity Metrics

Computed on generated parcels (Voronoi, Grid, or Momepy tessellation):

- **BCR – Building Coverage Ratio**  
  Percentage of parcel area covered by building footprints (0–100).  
  Used to evaluate horizontal land-use intensity and zoning compliance.

- **FAR – Floor Area Ratio**  
  Ratio between total building floor area and parcel area (≥ 0).  
  Represents vertical development intensity and land-use capacity.

- **OSR – Open Space Ratio**  
  Defined as `1 − (BCR / 100)` (0–1).  
  Indicates the proportion of unbuilt or open space.

- **Vegetation Proportion (optional)**  
  Fraction of vegetated pixels within each parcel derived from NDVI thresholding.  
  Useful for green space screening and environmental assessment.

---

### 2. Building Morphometric Indicators

Calculated directly from building geometries:

- **Compactness**  
  Measures how close a footprint is to a circular shape.  
  Relevant for urban efficiency and energy-related morphology studies.

- **Shape Index**  
  Indicates footprint irregularity and geometric complexity.

- **Fractal Dimension**  
  Quantifies boundary complexity of building footprints.

- **Solidity**  
  Ratio of footprint area to its convex hull area.  
  Indicates indentation and built-form cohesion.

- **Convexity**  
  Ratio between convex hull perimeter and actual perimeter.  
  Highlights concave or articulated building shapes.

- **Elongation**  
  Ratio of the longest to the shortest side of the bounding box.  
  Supports both axis-aligned and oriented minimum bounding box methods.

---

## 🧩 Parcelization Methods

MORPHINT supports three parcel generation strategies:

1. **Voronoi Tessellation**  
   Parcels generated from building centroids.  
   Suitable for density approximation where cadastral data are unavailable.

2. **Grid Tessellation**  
   Regular square grid with user-defined cell size.  
   Recommended for rapid city-wide screening and comparative studies.

3. **Momepy Tessellation**  
   Morphology-aware parcels generated from building footprints.  
   Optional tiling is available for large datasets to improve performance.

All parcel outputs can optionally be **clipped to a boundary polygon**.

---

## 🏙️ Applications in Urban Planning

MORPHINT is designed to support:

- Zoning and regulation evaluation (RTRW, RDTR)
- Density and land-use intensity analysis
- Urban infill versus sprawl diagnostics
- Morphological characterization of urban blocks
- Urban design and form-based planning
- Green space and open space assessment
- Teaching and studio-based urban analysis
- Research in urban morphology and spatial planning

---

## ⚙️ Requirements

- **QGIS** ≥ 3.40  
- Recommended CRS: **Projected coordinate system (meters)**

### Optional Python Dependencies (for Momepy mode)

- NumPy  
- GeoPandas  
- Shapely  
- Momepy  

---

## 📥 Installation

1. Clone or download this repository:
   ```bash
   git clone https://github.com/firmanaf/MORPHINT.git

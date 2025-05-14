# Crop Productivity Monitoring with MODIS and Google Earth Engine

This Python module provides a modular and scalable pipeline to assess agricultural productivity trends using MODIS EVI data and Google Earth Engine (GEE).

---

## 🚀 Features

- Load and clean MODIS Terra and Aqua EVI time series
- Apply cropland masks using DEA Crop Mask
- Compute seasonal composites (mean EVI)
- Derive z-scores against 10-year and 5-year baselines
- Export geospatial rasters to Google Cloud Storage
- Generate plots and statistical summaries
- Compute zonal statistics by administrative regions

---

## 📁 Project Structure

```
crop_productivity_project/
├── crop_productivity/            # Main module
│   ├── __init__.py
│   ├── analyzer.py               # Class-based pipeline
│   └── utils.py                  # Preprocessing, plotting, helpers
├── main.py                       # CLI entry point
├── cropProductivity.par          # Parameter configuration file
├── requirements.txt              # Python dependencies
├── setup.py                      # Installation script
├── README.md                     # This file
└── docs/                         # Documentation and images
```

---

## 🛠️ Installation

1. Clone or download the repository:
```bash
cd crop-prod
```

2. Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Install the module:
```bash
pip install -e .
```

---

## 📌 Configuration

Edit the `cropProductivity.par` file to specify:
- Date ranges (overall and per year)
- Growing season months
- ADM boundary paths
- Google Cloud bucket name
- Input/output directories

---

## 🧪 Usage

Run the CLI pipeline using:

```bash
python main.py --params cropProductivity.par --output GEE
```

Or after pip install:
```bash
cropproductivity --params cropProductivity.par --output GEE
```

---

## 📊 Output

- CSV: zonal statistics by ADM regions
- GeoTIFFs: EVI z-score rasters by year
- Plots: seasonal trends, histograms, anomaly summaries

---

## 🧩 Dependencies

- Google Earth Engine Python API
- pandas, geopandas, numpy
- rasterio, branca, folium
- matplotlib, seaborn
- gee-zonal, google-cloud-storage

---

## 📧 Contact

Created by **Pietro Milillo**. Contributions, issues, and feature requests welcome.

#### Data Availability Statement

Restrictions may apply to the data that support the findings of this study. Data received from the private sector through the Development Data Partnership are subject to the terms and conditions of the data license agreement and the “Official Use Only” data classification. These data are available upon request through the [Development Data Partnership](https://datapartnership.org/). Licensing and access information for all other datasets are included in the documentation.


## License

This projects is licensed under the [**Mozilla Public License**](https://opensource.org/license/mpl-2-0/) - see the [LICENSE](LICENSE) file for details.


## Code of Conduct

The World Bank Data Lab <span style="color:#3EACAD">template</span> used to create this project maintains a [Code of Conduct](docs/CODE_OF_CONDUCT.md) to ensure an inclusive and respectful environment for everyone. Please adhere to it in all interactions within our community.

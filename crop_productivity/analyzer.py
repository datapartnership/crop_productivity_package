"""
crop_productivity_analyzer.py

Class-based implementation of crop productivity analysis using MODIS EVI data,
Earth Engine, and geospatial analysis tools.
"""

import os
import ee
import json
import pandas as pd
import geopandas as gpd
from datetime import datetime
from gee_zonal import ZonalStats, gee_helpers

class CropProductivityAnalyzer:
    def __init__(self, param_file):
        """
        Constructor that loads parameters from a file and sets up initial state.
        
        Args:
            param_file (str): Path to the .par configuration file.
        """
        self.params = self.load_parameters(param_file)
        self.adm0 = None
        self.adm1 = None
        self.adm2 = None
        self.adm3 = None
        self.eth_geometry = None
        self.crop_mask = None
        self.mxd13q1 = None

    def load_parameters(self, path):
        """
        Reads key=value pairs from a parameter file.

        Args:
            path (str): Path to the .par file

        Returns:
            dict: dictionary of parameters
        """
        print(f"Loading parameters from {path}")
        params = {}
        with open(path, "r") as f:
            for line in f:
                if "=" in line and not line.strip().startswith("#"):
                    key, value = line.strip().split("=", 1)
                    key = key.strip()
                    value = value.strip().strip('"').strip("'")
                    try:
                        value = int(value)
                    except ValueError:
                        pass
                    params[key] = value
        print(params)
        return params
        
    def initialize_gee(self):
        """
        Authenticates and initializes the Google Earth Engine (GEE) environment.
        Includes a check to confirm successful initialization.
        """
        print("Initializing Google Earth Engine...")

        try:
            # Try to initialize without authenticating (in case already initialized)
            if not ee.data._initialized:
                ee.Authenticate()
                ee.Initialize()

            # Confirm successful initialization
            if ee.data._initialized:
                print("Google Earth Engine successfully initialized.")
            else:
                raise Exception("Google Earth Engine initialization failed.")

        except Exception as e:
            print(f"GEE Initialization error: {e}")
            raise
        
    def load_boundaries(self):
        """
        Loads ADM0 to ADM3 boundaries as GeoDataFrames.
        Also exports ADM0 as GeoJSON and extracts the GEE geometry.
        """
        print("Loading boundary shapefiles...")
        self.adm0 = gpd.read_file(self.params["ADM0"])
        self.adm1 = gpd.read_file(self.params["ADM1"])
        self.adm2 = gpd.read_file(self.params["ADM2"])
        self.adm3 = gpd.read_file(self.params["ADM3"])

        # Export ADM0 to a standard UTF-8 GeoJSON file
        geojson_path = "data/boundaries/ethiopia_adm0.geojson"
        self.adm0.to_file(geojson_path, driver='GeoJSON')

        # Load raw GeoJSON and convert to ee.Geometry
        with open(geojson_path, encoding='utf-8') as f:
            eth_geometry = json.load(f)
        self.eth_geometry = ee.Geometry(eth_geometry['features'][0]['geometry'])

        print("Boundaries loaded and geometry extracted.")

    def load_modis_collections(self):
        """
        Loads and merges Terra and Aqua MODIS EVI collections.
        Applies cloud and quality masking.
        """
        print("Loading MODIS EVI collections...")
        start_date = ee.Date(self.params["start_date"])
        end_date = ee.Date(self.params["end_date"])
        
        
        terra = ee.ImageCollection("MODIS/061/MOD13Q1") \
            .select(["EVI", "SummaryQA", "DetailedQA"]) \
            .filterDate(start_date, end_date)

        aqua = ee.ImageCollection("MODIS/061/MYD13Q1") \
            .select(["EVI", "SummaryQA", "DetailedQA"]) \
            .filterDate(start_date, end_date)

        def bitwiseExtract(value, fromBit, toBit=None):
            if toBit is None:
                toBit = fromBit
            maskSize = ee.Number(1).add(toBit).subtract(fromBit)
            mask = ee.Number(1).leftShift(maskSize).subtract(1)
            return value.rightShift(fromBit).bitwiseAnd(mask)

        def modisQA_mask(image):
            sqa = image.select("SummaryQA")
            dqa = image.select("DetailedQA")
            viQualityFlagsS = bitwiseExtract(sqa, 0, 1)
            viQualityFlagsD = bitwiseExtract(dqa, 0, 1)
            viSnowIceFlagsD = bitwiseExtract(dqa, 14)
            viShadowFlagsD = bitwiseExtract(dqa, 15)
            mask = (
                viQualityFlagsS.eq(0)
                .And(viQualityFlagsD.eq(0))
                .And(viQualityFlagsS.eq(1))
                .And(viQualityFlagsD.eq(1))
                .And(viSnowIceFlagsD.eq(0))
            )
            return image.updateMask(mask)

        mod13q1_QC = terra.map(modisQA_mask)
        myd13q1_QC = aqua.map(modisQA_mask)

        merged = mod13q1_QC.select("EVI").merge(myd13q1_QC.select("EVI"))
        self.mxd13q1 = merged.sort("system:time_start")  # ✅ This is the fix

        print("MODIS EVI collection loaded and cleaned.")
        return self.mxd13q1  # Optional: return for inspection

    def apply_crop_mask(self):
        """
        Applies DEA cropland mask to the cleaned MODIS EVI image collection.
        Stores the result in self.mxd13q1.
        """
        print("Applying DEA cropland mask...")
        bool_dict = {
            "0": "ocean",
            "1": "non_crop",
            "2": "crop_irrigated",
            "3": "crop_rainfed",
        }

        dea = ee.ImageCollection(
            "projects/sat-io/open-datasets/DEAF/CROPLAND-EXTENT/mask"
        ).mosaic()
        self.crop_mask = dea.select("b1").rename("crop").clip(self.eth_geometry)
        

        def cropmask(img):
            return img.updateMask(self.crop_mask).clip(self.eth_geometry)

        self.mxd13q1 = self.mxd13q1.map(cropmask)
        # Count number of images in the ImageCollection
        count = self.mxd13q1.size().getInfo()
        print(f"Number of images in the collection: {count}")
        
        print("Cropland mask applied to MODIS collection.")
        return self.crop_mask, self.mxd13q1

    def apply_esa_2020v100_cropland_mask(self):
        """
        Applies ESA WorldCover 10m v100 cropland mask (class 40) to MODIS EVI collection.

        Returns:
            crop_mask (ee.Image): ESA cropland mask (class 40) clipped to AOI.
            masked_collection (ee.ImageCollection): MODIS collection masked by ESA cropland extent.
        """
        print("Applying ESA v100 2020 WorldCover cropland mask (class 40)...")

        # Load ESA WorldCover 10m land cover map (2020)
        esa_lc = ee.Image("ESA/WorldCover/v100/2020")

        # Select only pixels classified as Cropland (class value 40)
        crop_class = 40
        crop_mask = esa_lc.select("Map").eq(crop_class).rename("crop")

        # Clip to country boundary
        self.crop_mask = crop_mask.clip(self.eth_geometry)

        # Mask each MODIS image with the ESA cropland mask
        def apply_mask(image):
            return image.updateMask(self.crop_mask).clip(self.eth_geometry)

        self.mxd13q1 = self.mxd13q1.map(apply_mask)

        # Get number of remaining images
        try:
            count = self.mxd13q1.size().getInfo()
            print(f"Number of MODIS images after ESA masking: {count}")
        except Exception as e:
            print("Could not get image count. Reason:", e)

        print("ESA v100 2020 cropland mask successfully applied.")
        return self.crop_mask, self.mxd13q1
    
    
    def generate_evi_composites(self):
        """
        Generates mean seasonal EVI composites for different years:
        - 10-year mean (EVI_MEAN)
        - Year 1 (EVI_Y1)
        - Year 2 (EVI_Y2)
        
        Returns:
            ee.Image: A stacked image of EVI composites
        """
        print("Generating EVI composites...")

        gssm = self.params['gssm']
        gsem = self.params['gsem']

        sdmegs = self.params['sdmegs']
        edmegs = self.params['edmegs']
        sdy1egs = self.params['sdy1egs']
        edy1egs = self.params['edy1egs']
        sdy2egs = self.params['sdy2egs']
        edy2egs = self.params['edy2egs']

        growing_season = ee.Filter.calendarRange(gssm, gsem, 'month')

        evi_mean = self.mxd13q1.filterDate(sdmegs, edmegs).filter(growing_season).mean().rename("EVI_MEAN")
        evi_y1 = self.mxd13q1.filterDate(sdy1egs, edy1egs).filter(growing_season).mean().rename("EVI_Y1")
        evi_y2 = self.mxd13q1.filterDate(sdy2egs, edy2egs).filter(growing_season).mean().rename("EVI_Y2")

        stacked = evi_mean.addBands(evi_y1).addBands(evi_y2)
        print("EVI composites generated and stacked.")
        return stacked

    def compute_zonal_statistics(self, image, feature_level="adm1", stat_type="mean", scale=10, output_name="evi_stats", output_dir="GEE_UTILS"):
        """
        Computes zonal statistics using the given image over specified boundary level.

        Args:
            image (ee.Image): Image to extract stats from
            feature_level (str): One of 'adm0', 'adm1', 'adm2', 'adm3'
            stat_type (str): Statistic type (e.g., 'mean', 'median', 'count')
            scale (int): Spatial resolution in meters
            output_name (str): Name of the output
            output_dir (str): Output directory

        Returns:
            pd.DataFrame: DataFrame of zonal stats
        """
        print(f"Computing {stat_type} zonal statistics at {feature_level.upper()} level...")

        level_map = {
            "adm0": self.adm0,
            "adm1": self.adm1,
            "adm2": self.adm2,
            "adm3": self.adm3,
        }

        target_features = gee_helpers.gpd_to_gee(level_map[feature_level])
        # Fix starts here

        zs = ZonalStats(
            ee_dataset=image,
            target_features=target_features,
            statistic_type=stat_type,
            scale=scale,
            output_name=output_name,
            output_dir=output_dir
        )

        zs.runZonalStats()

        df = pd.read_csv(f"{output_dir}/{output_name}.csv")
        print(f"Zonal statistics saved to {output_dir}/{output_name}.csv")
        return df



    def subtract_five_years(self, date_str):
        """
        Utility function to subtract 5 years from a date string (YYYY-MM-DD).
        """
        from dateutil.relativedelta import relativedelta
        date_obj = datetime.strptime(date_str, "%Y-%m-%d")
        return (date_obj - relativedelta(years=5)).strftime("%Y-%m-%d")

    def extract_years_between(self, start_date, end_date):
        """
        Returns a list of years between two date strings.
        """
        start = datetime.strptime(start_date, "%Y-%m-%d").year
        end = datetime.strptime(end_date, "%Y-%m-%d").year
        return list(range(start, end + 1))

    def calculate_evi_zscore(self, year, historical_median, historical_stddev):
        """
        Computes the EVI z-score for a given year against historical stats.
        """
        gssm = self.params['gssm']
        gsem = self.params['gsem']
        growing_season = ee.Filter.calendarRange(gssm, gsem, 'month')
        
        start_date = f"{year}-01-01"
        end_date = f"{year}-12-31"

        year_evi = self.mxd13q1.filterDate(start_date, end_date)                               .filter(growing_season)                               .median()                               .clip(self.eth_geometry)

        evi_zscore = year_evi.subtract(historical_median)                             .divide(historical_stddev)                             .rename("EVI_ZScore")

        return evi_zscore

    def compute_all_zscores(self):
        """
        Calculates EVI z-scores for all years between filterStartDate and filterEndDate.
        Returns dictionaries of z-scores for 10-year and 5-year references.
        """
        print("Calculating historical statistics and yearly z-scores...")

        sdmegs = self.params['sdmegs']
        edmegs = self.params['edmegs']
        short_edmegs = self.subtract_five_years(edmegs)
        growing_season = ee.Filter.calendarRange(self.params['gssm'], self.params['gsem'], 'month')

        hist_collection = self.mxd13q1.filterDate(sdmegs, edmegs).filter(growing_season).select("EVI")
        hist_collection_5yr = self.mxd13q1.filterDate(sdmegs, short_edmegs).filter(growing_season).select("EVI")

        median_10yr = hist_collection.median().clip(self.eth_geometry)
        median_5yr = hist_collection_5yr.median().clip(self.eth_geometry)
        stddev_10yr = hist_collection.reduce(ee.Reducer.stdDev()).clip(self.eth_geometry)
        stddev_5yr = hist_collection_5yr.reduce(ee.Reducer.stdDev()).clip(self.eth_geometry)

        start = self.params['filterStartDate']
        end = self.params['filterEndDate']
        years = self.extract_years_between(start, end)

        zscores_10yr = {}
        zscores_5yr = {}
        for year in years:
            zscores_10yr[year] = self.calculate_evi_zscore(year, median_10yr, stddev_10yr)
            zscores_5yr[year] = self.calculate_evi_zscore(year, median_5yr, stddev_5yr)

        print("Z-score computation complete.")
        return zscores_10yr, zscores_5yr

    def export_to_gcs(self, image, gcs_file_name):
        """
        Exports an image to Google Cloud Storage (GCS).

        Args:
            image (ee.Image): The image to export
            gcs_file_name (str): Name to use for the exported file
        """
        bucket = self.params["bucket_name"]
        print(f"Exporting {gcs_file_name} to GCS bucket {bucket}...")

        export_task = ee.batch.Export.image.toCloudStorage(
            image=image,
            description=gcs_file_name,
            bucket=bucket,
            fileNamePrefix=gcs_file_name,
            scale=250,
            region=self.eth_geometry,
            maxPixels=1e9
        )

        export_task.start()

        import time
        while export_task.status()['state'] in ['READY', 'RUNNING']:
            print(f"Task status: {export_task.status()['state']}")
            time.sleep(30)

        print(f"Task complete with state: {export_task.status()['state']}")

    def batch_export_zscores(self, zscore_dict):
        """
        Automatically export all EVI z-score images in the provided dictionary to GCS.

        Args:
            zscore_dict (dict): Dictionary of year -> ee.Image
        """
        for year, image in zscore_dict.items():
            file_name = f"evi_zscore_{year}"
            print(f"Scheduling export for {file_name}")
            self.export_to_gcs(image, file_name)
        print("All export tasks initiated.")

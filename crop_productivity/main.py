"""
main.py

Command-line entry point for running the crop productivity analysis pipeline.
"""

import argparse
from crop_productivity.analyzer import CropProductivityAnalyzer

def main():
    # Set up CLI argument parser
    parser = argparse.ArgumentParser(description="Crop Productivity Analysis using MODIS and GEE")
    parser.add_argument('--params', type=str, required=True, help="Path to the .par configuration file")
    parser.add_argument('--output', type=str, default="GEE_ETH", help="Directory to save zonal stats output")

    args = parser.parse_args()

    # Initialize the analyzer
    analyzer = CropProductivityAnalyzer(args.params)

    # Execute the analysis pipeline
    analyzer.initialize_gee()
    analyzer.load_boundaries()
    analyzer.load_modis_collections()
    # Inside main.py, after initializing analyzer and loading MODIS
    crop_mask_type = analyzer.params.get("cropland", "DEA").strip().upper()

    if crop_mask_type == "ESAV100":
        masked, mxd13q1 = analyzer.apply_esa_cropland_mask()
    elif crop_mask_type == "DEA":
        masked, mxd13q1 = analyzer.apply_crop_mask()
    else:
        raise ValueError(f"Unsupported cropland mask type: {crop_mask_type}")

    # Compute zonal statistics Number of crop pixels in each admin region
    #analyzer.compute_zonal_statistics(
    #    image=masked,
    #    feature_level="adm1",
    #    stat_type="count",
    #    output_name="adm1_crop_count",
    #    output_dir=args.output
    #)

    #mxd13q1_recent = mxd13q1.filterDate(filterStartDate, filterEndDate)

    # Convert your AOI from GeoDataFrame to GEE FeatureCollection
    #aoi = gee_helpers.gpd_to_gee(analyzer.adm1)

    #zs = ZonalStats(
    #    ee_dataset=mxd13q1_recent,
    #    target_features=aoi,
    #    statistic_type="median",
    #    scale=250,
    #    output_name="eth_phenology",
    #    output_dir="GEE_ETH"
    #)

    #zs.runZonalStats()
    
    # Generate seasonal composites
    evi_composites = analyzer.generate_evi_composites()

    # Compute zonal statistics
    #analyzer.compute_zonal_statistics(
    #    image=evi_composites,
    #    feature_level="adm1",
    #    stat_type="mean",
    #    output_name="evi_composite_stats",
    #    output_dir=args.output
    #)

    # Compute z-scores and export
    zscore_10yr, _ = analyzer.compute_all_zscores()
    analyzer.batch_export_zscores(zscore_10yr)

    # ➕ Compute cropland histogram stats
    from crop_productivity.utils import compute_crop_histogram_from_mask
    df_crop_hist = compute_crop_histogram_from_mask(analyzer.adm1, analyzer.crop_mask)
    print("\n🧮 Crop Histogram Summary:")
    print(df_crop_hist.head())
    output_csv = os.path.join(args.output, 'adm1_crop_histogram_stats.csv')
    os.makedirs(args.output, exist_ok=True)
    df_crop_hist.to_csv(output_csv, index=False)
    print(f"✔ Crop histogram saved to {output_csv}")

    print("Crop productivity analysis completed successfully.")

if __name__ == "__main__":
    main()

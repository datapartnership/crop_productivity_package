"""
utils.py

General-purpose utility functions for crop productivity analysis.
Includes time series processing, plotting, and date utilities.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from datetime import datetime
from dateutil.relativedelta import relativedelta
import seaborn as sns
import rasterio
import geopandas as gpd
import numpy as np
from rasterstats import zonal_stats

def preprocess_series(series):
    """
    Preprocess EVI time series using TIMESAT-style filtering.

    Steps:
    1. Remove outliers
    2. Interpolate missing values
    3. Apply Savitzky-Golay smoothing

    Args:
        series (pd.Series): Time series of EVI values

    Returns:
        np.ndarray: Smoothed EVI series
    """
    series = series.interpolate(limit_direction="both")
    median = series.rolling(window=5, center=True).median()
    std_dev = series.std()
    series = series.mask((series - median).abs() > 2 * std_dev)
    series = series.interpolate(limit_direction="both")
    return savgol_filter(series, window_length=5, polyorder=2)


def extract_sos_mos_eos(smoothed, dates, threshold=0.2):
    """
    Extract Start (SOS), Maximum (MOS), and End (EOS) of Season from EVI curve.

    Args:
        smoothed (np.ndarray): Smoothed EVI time series
        dates (pd.DatetimeIndex): Corresponding dates
        threshold (float): Threshold fraction of amplitude

    Returns:
        tuple: (SOS, MOS, EOS) as pd.Timestamp or None
    """
    max_val = smoothed.max()
    min_val = smoothed.min()
    amp = max_val - min_val
    sos = mos = eos = None
    for i in range(1, len(smoothed)):
        if sos is None and smoothed[i] > min_val + threshold * amp:
            sos = dates[i]
        if smoothed[i] == max_val:
            mos = dates[i]
        if sos is not None and smoothed[i] < min_val + threshold * amp:
            eos = dates[i]
            break
    return sos, mos, eos


def subtract_five_years(date_str):
    """
    Subtract five years from a date string (YYYY-MM-DD).

    Args:
        date_str (str): Input date

    Returns:
        str: New date string
    """
    date_obj = datetime.strptime(date_str, "%Y-%m-%d")
    return (date_obj - relativedelta(years=5)).strftime("%Y-%m-%d")


def extract_years_between(start_date, end_date):
    """
    Return a list of years between two date strings.

    Args:
        start_date (str): Start date
        end_date (str): End date

    Returns:
        list: List of years (int)
    """
    start = datetime.strptime(start_date, "%Y-%m-%d").year
    end = datetime.strptime(end_date, "%Y-%m-%d").year
    return list(range(start, end + 1))


def plot_zscore_histograms(file_paths, labels=None, colors=None, output_path=None):
    """
    Plot histograms for EVI Z-score raster files and compare distributions.

    Args:
        file_paths (list): List of .tif file paths
        labels (list): Optional labels for each file
        colors (list): Optional color scheme
        output_path (str): Optional file save path
    """
    if not file_paths:
        raise ValueError("At least one file must be specified.")

    if labels is None:
        labels = [f"File {i+1}" for i in range(len(file_paths))]
    if colors is None:
        colors = sns.color_palette("husl", len(file_paths))
    if len(labels) != len(file_paths) or len(colors) != len(file_paths):
        raise ValueError("Labels and colors must match number of files.")

    stats_list = []
    all_data = []

    for file_path in file_paths:
        with rasterio.open(file_path) as src:
            z = src.read(1)
            nodata = src.nodata
        if nodata is not None:
            z = z[z != nodata]
        z = z[(z >= -5) & (z <= 5)]

        stats = {
            'mean': np.mean(z),
            'median': np.median(z),
            'std': np.std(z),
            'sig_neg': np.sum(z < -1.96) / len(z) * 100,
            'sig_pos': np.sum(z > 1.96) / len(z) * 100
        }

        all_data.append(z)
        stats_list.append(stats)

    fig, ax = plt.subplots(figsize=(12, 6))
    for z, label, color in zip(all_data, labels, colors):
        sns.histplot(z, kde=True, bins=50, color=color, label=label, ax=ax, stat="density", alpha=0.6)

    ax.axvline(0, linestyle="--", color="black", label="Mean")
    ax.axvline(-1.96, linestyle="--", color="red", label="Below Normal")
    ax.axvline(1.96, linestyle="--", color="green", label="Above Normal")
    ax.set_title("EVI Z-Score Distributions")
    ax.set_xlabel("Z-Score")
    ax.set_ylabel("Density")
    ax.legend()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.tight_layout()
    plt.show()

    print("\nStatistical Summary:")
    for label, stats in zip(labels, stats_list):
        print(f"{label}: Mean={stats['mean']:.2f}, Median={stats['median']:.2f}, Std={stats['std']:.2f}, "
              f"< -1.96 = {stats['sig_neg']:.1f}%, > 1.96 = {stats['sig_pos']:.1f}%")



from gee_zonal import ZonalStats, gee_helpers
import ast

def compute_crop_histogram_from_mask(adm1_gdf, crop_mask_ee, scale=10):
    """
    Computes cropland/non-cropland pixel histograms and derives area stats by ADM1 region.

    Parameters:
        adm1_gdf (GeoDataFrame): ADM1 boundaries
        crop_mask_ee (ee.Image): DEA crop mask with binary classification
        scale (int): Spatial resolution in meters

    Returns:
        pd.DataFrame: Table with crop area stats and shares
    """
    bool_dict = {0: "non_crop", 1: "crop"}

    aoi_adm1 = gee_helpers.gpd_to_gee(adm1_gdf)

    zs = ZonalStats(
        ee_dataset=crop_mask_ee,
        target_features=aoi_adm1,
        statistic_type="histogram",
        scale=scale,
        output_name="adm1_crop_count_runtime",
        output_dir="GEE_ETH"
    )
    zs.runZonalStats()

    df = pd.read_csv("GEE_ETH/adm1_crop_count_runtime.csv")

    def load_dict(input_str):
        input_str = input_str.replace("null", "'null'")
        input_str = input_str.replace("=", ":")
        return ast.literal_eval(input_str)

    df["histogram"] = df["histogram"].apply(load_dict)

    histograms = df["histogram"]
    hist_dfs = [pd.DataFrame(index=[idx], data=h) for idx, h in histograms.items()]
    hist_df = pd.concat(hist_dfs).fillna(0)
    hist_df.rename(columns=bool_dict, inplace=True)

    df = df.join(hist_df)
    df["Crop Area ha."] = ((df["crop"] * 100) / 10000).apply(lambda x: f"{x:,.0f}")
    df["Crop Area Share (% of Region)"] = (df["crop"] / (df["crop"] + df["non_crop"])).apply(lambda x: f"{x:.2%}")
    df["Crop Area Share (% of Country)"] = (df["crop"] / df["crop"].sum()).apply(lambda x: f"{x:.2%}")
    df.rename(columns={"ADM1_EN": "Region"}, inplace=True)

    return df[["Region", "Crop Area ha.", "Crop Area Share (% of Country)", "Crop Area Share (% of Region)"]].reset_index(drop=True)


def aggregate_zscores(raster_path, adm1):
    """Aggregate pixel Z-scores to ADM1 regions"""
    # Load admin boundaries
    #adm1 = gpd.read_file(adm1_path)
    
    # Calculate zonal statistics
    stats = zonal_stats(
        adm1,
        raster_path,
        stats=['mean', 'median', 'min', 'max', 'std', 'count'],
        geojson_out=True
    )
    
    # Convert to GeoDataFrame
    result = gpd.GeoDataFrame.from_features(stats)
    
    # Add original attributes
    for col in ['ADM1_EN']:
        if col in adm1.columns:
            result[col] = adm1[col].values
    
    return result
    
    
    
def plot_zscore_maps(gdfs, stat='mean', titles=None, figsize_per_plot=(8, 4), main_title='EVI Z-Score'):
    """
    Plot aggregated Z-scores by region for multiple GeoDataFrames
    
    Parameters:
    gdfs: List of GeoDataFrames to plot
    stat: Single statistic to plot for each GDF (default 'mean')
    titles: Optional list of titles for each plot (default None)
    figsize_per_plot: Size of each individual plot
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    # If gdfs is not a list, convert it to one
    if not isinstance(gdfs, list):
        gdfs = [gdfs]
    
    # Calculate number of plots and grid dimensions
    n_plots = len(gdfs)
    
    # Calculate grid layout
    if n_plots == 1:
        nrows, ncols = 1, 1
    elif n_plots == 2:
        nrows, ncols = 1, 2
    elif n_plots <= 4:
        nrows, ncols = 2, 2
    elif n_plots <= 6:
        nrows, ncols = 2, 3
    else:
        # For more than 6, use a 3 column layout
        ncols = 3
        nrows = int(np.ceil(n_plots / ncols))
    
    # Create figure with appropriate size
    width = figsize_per_plot[0] * ncols
    height = figsize_per_plot[1] * nrows
    fig, axes = plt.subplots(nrows, ncols, figsize=(width, height))
    
    # Handle case where axes is not a 2D array
    if n_plots == 1:
        axes = np.array([[axes]])
    elif nrows == 1 or ncols == 1:
        axes = axes.reshape(nrows, ncols)
    
    # Generate default titles if none provided
    if titles is None:
        titles = [f'EVI Z-Score by Region {i+1}' for i in range(n_plots)]
    elif len(titles) < n_plots:
        # Extend titles list if needed
        titles.extend([f'EVI Z-Score by Region {i+1}' for i in range(len(titles), n_plots)])
    
    # Plot each GeoDataFrame
    for idx, gdf in enumerate(gdfs):
        row_idx = idx // ncols
        col_idx = idx % ncols
        ax = axes[row_idx, col_idx]
        
        if stat in gdf.columns:
            gdf.plot(column=stat, cmap='BrBG', linewidth=0.5,
                    edgecolor='0.5', legend=True, ax=ax, vmin=-3, vmax=3)
            ax.set_title(titles[idx])
        else:
            # If stat not found, show empty plot with message
            ax.text(0.5, 0.5, f'"{stat}" not found in GDF {idx+1}',
                   transform=ax.transAxes, ha='center', va='center', fontsize=12)
            ax.set_title(f'Missing: {stat} in {titles[idx]}')
        
        ax.set_axis_off()
    
    # Remove any unused subplots
    for idx in range(n_plots, nrows * ncols):
        row_idx = idx // ncols
        col_idx = idx % ncols
        fig.delaxes(axes[row_idx, col_idx])

    plt.suptitle(main_title, fontsize=16, weight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    return fig

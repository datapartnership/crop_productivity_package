from setuptools import setup, find_packages

setup(
    name='crop_productivity',
    version='0.1.0',
    description='Crop productivity analysis using MODIS and Google Earth Engine',
    author='Your Name',
    author_email='pmilillo@worldbank.org',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'earthengine-api',
        'geopandas',
        'pandas',
        'numpy',
        'rasterio',
        'folium',
        'matplotlib',
        'branca',
        'gee-zonal',
        'seaborn',
        'scipy',
        'google-cloud-storage'
    ],
    entry_points={
        'console_scripts': [
            'cropproductivity = crop_productivity.main:main',
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)

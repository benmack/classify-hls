import gdal
import geopandas as gpd
from pathlib import Path
import rasterio

def create_clc_raster_tile(path_single_tile_footprint, 
                           path_clc_raster, 
                           path_clc_raster_tile, overwrite=False):
    """Clip the CLC raster with vector footprint of a s2 tile."""
    if not Path(path_clc_raster_tile).exists() or overwrite:

        with rasterio.open(path_clc_raster) as clc_rio:
            pass


        tile_footprint = gpd.read_file(path_single_tile_footprint)
        tile_footprint = tile_footprint.to_crs(clc_rio.crs)
        tile_bounds = tile_footprint.bounds.values[0]

        from osgeo import gdal

        clc = gdal.Open(str(path_clc_raster))
        clc = gdal.Translate(str(path_clc_raster_tile), clc, 
                             projWin = [tile_bounds[0],
                                        tile_bounds[3],
                                        tile_bounds[2],
                                        tile_bounds[1]])

        clc = None

        assert Path(path_clc_raster_tile).exists()
        return 0    

create_clc_raster_tile(path_single_tile_footprint=snakemake.input.single_tile_footprint, 
                       path_clc_raster=snakemake.input.clc_raster, 
                       path_clc_raster_tile=snakemake.output[0], 
                       overwrite=True)

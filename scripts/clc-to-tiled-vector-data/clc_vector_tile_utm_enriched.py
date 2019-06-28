import geopandas as gpd
import pandas as pd
from pathlib import Path

def create_clc_vector_tile_utm_enriched(path_single_tile_footprint, 
                                        path_clc_vector_tile, 
                                        path_clc_legend,
                                        path_clc_vector_tile_utm, overwrite=False):
    """Reproject CLC vector tile to utm projection of the respective sentinel-2 tile, add polygon id, area and make subsets."""
    
    assert Path(path_single_tile_footprint).exists()
    assert Path(path_clc_vector_tile).exists()
    assert Path(path_clc_legend).exists()
    
    if not Path(path_clc_vector_tile_utm).exists() or overwrite:
        
        footprint = gpd.read_file(path_single_tile_footprint)
        epsg_utm = footprint["epsg"][0]
        tile = footprint["name"][0]

        clc_vector_tile = gpd.read_file(path_clc_vector_tile)
        clc_legend = pd.read_csv(path_clc_legend)

        clc_vector_tile = clc_vector_tile.rename({"DN": "cid_l3"}, axis=1)
        clc_vector_tile["pid"] = range(1, len(clc_vector_tile)+1)
        clc_vector_tile["area"] = clc_vector_tile.area
        clc_vector_tile["tile"] = tile

        clc_vector_tile = clc_vector_tile.merge(
            clc_legend[["cid_l1", "cid_l2", "cid_l3"]], on="cid_l3", how="left")

        clc_vector_tile_utm = clc_vector_tile.to_crs(epsg=epsg_utm)

        clc_vector_tile_utm.to_file(path_clc_vector_tile_utm, driver="GPKG")
    return 0

create_clc_vector_tile_utm_enriched(path_single_tile_footprint=snakemake.input[0], 
                                    path_clc_vector_tile=snakemake.input[1],
                                    path_clc_legend=snakemake.input[2], 
                                    path_clc_vector_tile_utm=snakemake.output[0], 
                                    overwrite=True)

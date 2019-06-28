import geopandas as gpd
from pathlib import Path

def create_clc_vector_tile_subset(path_clc_vector_tile_utm,
                                  max_area_threshold,
                                  path_clc_vector_tile_utm_subset,
                                  overwrite=False):
    """Reproject CLC vector tile to utm projection of the respective sentinel-2 tile, add polygon id and make subsets."""
    if not Path(path_clc_vector_tile_utm_subset).exists() or overwrite:
        gdf = gpd.read_file(path_clc_vector_tile_utm)
        gdf_sub = gdf.loc[gdf["area"] < max_area_threshold]
        gdf_sub.to_file(path_clc_vector_tile_utm_subset, driver="GPKG")
    return 0

create_clc_vector_tile_subset(path_clc_vector_tile_utm=snakemake.input[0],
                              max_area_threshold=float(snakemake.params[0]),
                              path_clc_vector_tile_utm_subset=snakemake.output[0],
                              overwrite=True)

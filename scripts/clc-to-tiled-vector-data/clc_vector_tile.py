from pathlib import Path
import subprocess

import eobox

def create_clc_vector_tile(path_clc_raster_tile, path_clc_vector_tile, overwrite=False):
    """Polygonize the CLC raster tile."""
    if not Path(path_clc_vector_tile).exists() or overwrite:
        path_gdal_polygonize_py = eobox.raster.gdalutils.POLYGONIZE_PATH
        cmd = f"python {path_gdal_polygonize_py} -b 1 {path_clc_raster_tile} {path_clc_vector_tile} -f 'GPKG' None DN"
        subprocess.check_call(cmd, shell=True)
    return 0

create_clc_vector_tile(path_clc_raster_tile=snakemake.input[0], 
                       path_clc_vector_tile=snakemake.output[0], 
                       overwrite=True)

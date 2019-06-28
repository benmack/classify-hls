from pathlib import Path
from nasa_hls import convert_hdf2tiffs


hdf_path = snakemake.input[0]
bands = snakemake.params.bands
dir__hls_tif = snakemake.output[0]

convert_hdf2tiffs(hdf_path=hdf_path, 
                  dstdir=Path(dir__hls_tif).parent.parent, 
                  bands=bands, 
                  max_cloud_coverage=100)

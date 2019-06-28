from pathlib import Path
import os

from src import configs
prjconf = configs.ProjectConfigParser()

TILES = prjconf.get("Params", "tiles").split(" ")
BANDS = ["Red", "NIR", "SWIR1", "SWIR2", "QA"]
# do we need that when it is also in the included Snakemake ?:
# CLC_MAX_AREA_THRESHOLD = prjconf.get("Params", "clc_max_area_threshold")


# snakemake -n -r
# snakemake -j 4
# snakemake --cores 8
# snakemake --dag data/interim/hls/{tile}/**/*__CLEAR.tif | dot -Tsvg > notebooks/01_raw/dag-graphs/dag_r2i_hls-prepare-l2-single-layer-files.svg
# snakemake --dag data/interim/hls/**/**/*__CLEAR.tif | dot -Tsvg > notebooks/01_raw/dag-graphs/dag_r2i_hls-prepare-l2-single-layer-files.svg


### clc-to-tiled-vector-data
# snakemake -s Snakefile_clc-to-tiled-vector-data --cores 8
# snakemake -s Snakefile_clc-to-tiled-vector-data --dag data/interim/clc/clc2018_*_subset_*.gpkg | dot -Tsvg > notebooks/01_raw/dag-graphs/dag_r2i_clc-to-tiled-vector-data.svg


cwd = os.getcwd()
expected_paths = []
for tile in TILES: #[TILES[0]]:
    print(tile)
    #expected_paths += [str(pth.relative_to(cwd)) for pth in list((prjconf.get_path("Raw", "hls") / tile).glob("*.hdf"))]
    sids_this_tile = [Path(pth).stem for pth in (prjconf.get_path("Raw", "hls") / tile).rglob("*.hdf")]
    absolut_paths_this_tile = [prjconf.get_path("Interim", "hls") / tile / sid / (sid + "__CLEAR.tif") for sid in sids_this_tile]
    expected_paths += [pth.relative_to(cwd) for pth in absolut_paths_this_tile]



rule all:
    input:
        expected_paths
        # expand("data/interim/hls/{tile}/**/*__CLEAR.tif", tile=TILES)
        # expand("data/interim/hls/{tile}/{sceneid}", tile=TILES, sceneid=SCENEIDS)

rule hls_create_clear_sky_band:
    input:
        "data/interim/hls/{tile}/{sceneid}/{sceneid}__QA.tif"
    output:
        "data/interim/hls/{tile}/{sceneid}/{sceneid}__CLEAR.tif"
    script:
        "scripts/hls-prepare-l2-single-layer-files/hls_create_clear_sky_band.py"

rule convert_hls_hdf_to_tif:
    input:
        "data/raw/hls/{tile}/{sceneid}.hdf"
    output:
        "data/interim/hls/{tile}/{sceneid}/{sceneid}__QA.tif"
        # directory("data/interim/hls/{tile}/{sceneid}")
        # WORKED - but in the script with Path(dir__hls_tif).parent.parent
        # expand("data/interim/hls/{tile}/{sceneid}/{sceneid}__{band}.tif", band=BANDS)
        # => No values given for wildcard 'tile'.
    params:
        bands = BANDS
    script:
        "scripts/hls-prepare-l2-single-layer-files/hls_convert_hdf_to_tif.py"

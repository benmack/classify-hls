import configparser
import geopandas as gpd
import json
import sys
from pathlib import Path
import pandas as pd

PROJECT_ROOT_DIR = rootdir = Path(sys.modules[__name__].__file__).parent.parent

class ProjectConfigParser(configparser.ConfigParser):
    def __init__(self, config_file=None):
        """Provides access to the project configuration file and more.
        
        Arguments:
            config_file {str} -- The configuration file, if ``None`` it is assumed to be the file \*project_configs.ini\* in the project root directory.  
        """

        configparser.ConfigParser.__init__(self, interpolation=configparser.ExtendedInterpolation())
        self._config_file = PROJECT_ROOT_DIR / "project_configs.ini"
        self.read(self._config_file)
    
    @property
    def config_file(self):
        """Get the file path of the configuration file."""
        return self._config_file

    def print_config_file_content(self):
        with open(self.config_file) as src:
            print(src.read())

    def get_path(self, section, option, tile=None):
        """Get an option from the configuration file and convert it to a Path object.
        
        Arguments:
            section {str} -- See ``configparser.ConfigParser.get()`` documentation.
            option {str} -- See ``configparser.ConfigParser.get()`` documentation.
            tile {str} -- If the path in the settings file contains a placeholder for the tile 
                (e.g. *./data/interim/s2_{tile}.tif*) and tile is given (e.g. *32UNU*) then the 
                function includes the tile in the path (e.g. returns *./data/interim/s2_32UNU.tif*)
        Returns:
            ``Path`` object -- ``pathlib.Path(self.get(section, option))``
        """
        pth = self.get(section, option)
        if tile:
            pth = pth.format(**{"tile": tile})
        return Path(pth)

    def get_clc_legend(self):
        """Get a look up table with the class IDs, names and colors of the three CLC levels."""
        clc_legend = pd.read_csv(self.get_path("Raw", "clc_legend"), delimiter=";").iloc[0:44, :]
        clc_legend.columns = ["l1_name", "l2_name", "l3_name", "grid_code", "rgb"]
        clc_legend_ids = clc_legend["l3_name"].str[:5].str.split(".", expand=True)
        clc_legend["l1_id"] = clc_legend_ids[0].astype("uint8")
        clc_legend["l2_id"] = (clc_legend_ids[0] + clc_legend_ids[1]).astype("uint8")
        clc_legend["l3_id"] = (clc_legend_ids[0] + clc_legend_ids[1] + clc_legend_ids[2]).astype("int")
        clc_legend["l1_name"] = clc_legend["l1_name"].str[3::]
        clc_legend["l2_name"] = clc_legend["l2_name"].str[4::]
        clc_legend["l3_name"] = clc_legend["l3_name"].str[6::]
        clc_legend = clc_legend.fillna(method="ffill")
        return clc_legend

    def get_clc_subset(self, tile="32UNU", area_threshold=500000):
        return self.get_path("Interim", "rootdir") / "clc" / f"clc2018_clc2018_v2018_20b2_raster100m_{tile}_epsg32623_area{area_threshold}.shp"

    def get_tilenames(self):
        return self.get("Params", "tiles").split(" ")

    def get_footprints(self, tile="ALL", epsg=4326):
        """Get tile footprint(s) as geopandas dataframe in the epsg specified or the original one if epsg=None."""
        if tile =="ALL":
            tilenames = self.get_tilenames()
        tile_footprints = []
        for tile in tilenames:
            path__tile_footprint = self.get_path("Raw", "tile_footprint", tile)
            if epsg is not None:
                tile_footprints.append(gpd.read_file(path__tile_footprint).to_crs(epsg=epsg))
            else:
                tile_footprints.append(gpd.read_file(path__tile_footprint))
        return pd.concat(tile_footprints)

    def get_scene_hdf(self, date, tile, product):
        """Get an option from the configuration file and convert it to a Path object.
        
        Arguments:
            date {str} -- e.g. 2018003
            tile {str} -- e.g. 32UNU
            product {str} -- L30 or S30
        Returns:
            ``Path`` object -- Scene path to the HDF file or directory with TIFFs.
        """
        pth = self.get("HLSL2", "scene_hdf")
        pth = pth.format(**{"date":date, "tile": tile, "product":product})
        return Path(pth)

    def get_scene_dir(self, date, tile, product):
        """Get an option from the configuration file and convert it to a Path object.
        
        Arguments:
            date {str} -- e.g. 2018003
            tile {str} -- e.g. 32UNU
            product {str} -- L30 or S30
        Returns:
            ``Path`` object -- Scene path to the HDF file or directory with TIFFs.
        """
        pth = self.get("HLSL2", "scene_dir_tif")
        pth = pth.format(**{"date":date, "tile": tile, "product":product})
        return Path(pth)

    def write_scene_collection(self, scenecoll, scenecoll_name, scenecoll_params, tile, exist_ok=False):
        sc_dir = self.get_path("Raw", "scene_colls")
        scenecoll_file = sc_dir / tile / ("df_" + scenecoll_name + ".csv")
        scenecoll_params_file = sc_dir / ("params_" + scenecoll_name + ".json")
        if scenecoll_params_file.exists():
            scenecoll_params_stored = self.read_scene_collection_params(scenecoll_name)
            if not scenecoll_params_stored == scenecoll_params:
                raise ValueError(f"The scene collection parameters already stored for {scenecoll_name} do not fit tothe given ones.")
        if scenecoll_file.exists() and exist_ok:
            print(f"{scenecoll_file} exists. It is NOT overwritten! Delete it first if you are sure to write it.")
        elif scenecoll_file.exists() and not exist_ok:
            raise ValueError(f"{scenecoll_file} exists. Delete it if you are sure to write it.")
        else:
            scenecoll_file.parent.mkdir(exist_ok=True, parents=True)
            scenecoll.to_csv(scenecoll_file, index=False)
            json.dump(scenecoll_params, open(scenecoll_params_file, "w"), indent=4)
            print(f" Scene collection written to {scenecoll_file}" )

    def read_scene_collection_params(self, scenecoll_name):
        sc_dir = self.get_path("Raw", "scene_colls")
        scenecoll_params_file = sc_dir / ("params_" + scenecoll_name + ".json")
        if not scenecoll_params_file.exists():
            raise ValueError(f"{scenecoll_file} doos not exists.")
        else:
            params = json.load(open(scenecoll_params_file,"r"))
            return params

    def read_scene_collection(self, scenecoll_name, tile):
        sc_dir = self.get_path("Raw", "scene_colls")
        scenecoll_file = sc_dir / tile / ("df_" + scenecoll_name + ".csv")
        scenecoll_params_file = sc_dir / ("params_" + scenecoll_name + ".json")
        if not scenecoll_file.exists():
            raise ValueError(f"{scenecoll_file} doos not exists.")
        else:
            df = pd.read_csv(scenecoll_file)
            df["date"] = pd.to_datetime(df["date"], format="%Y-%m-%d")
            return df

    def get_scene_collection_names(self):
        scenecolls = self.get_path("Raw", "scene_colls").glob("*.json")
        return [sc.stem.replace("params_", "") for sc in scenecolls]

    def get_layer_df_of_scene_collection(self, scenecoll_name, bands, tile):
        # print(prjconf.get_scene_collection_names())

        scoll = self.read_scene_collection(scenecoll_name, tile)

        dir_tiffs = self.get_path("Interim", "hls") / tile # TODO: use hls_layer: ${hls}/{tile}/{sid}/{sid}/{sid}__{band}.tif
        scoll["scenedir"] = [str(Path(dir_tiffs) / sceneid) for sceneid in scoll["sceneid"]]
        scoll.head()

        scoll_lays = []
        for i, row in scoll.iterrows():
            scoll_lays.append(
                pd.DataFrame({"sceneid":[row["sceneid"]]*len(bands),
                              "band": bands,
                              "uname": [f"{row['product']}_{row['date']}_{b}" for b in bands],
                              "product":[row["product"]]*len(bands),
                              "tile":[row["tile"]]*len(bands),
                              "date_Yj":[row["date_Yj"]]*len(bands),
                              "date":[row["date"]]*len(bands),
                              "cloud_cover":[row["cloud_cover"]]*len(bands),
                              "spatial_coverage":[row["spatial_coverage"]]*len(bands),
                              "path": [f"{row['scenedir']}/{row['sceneid']}__{b}.tif" for b in bands]
                              })
            )
        scoll_lays = pd.concat(scoll_lays).reset_index(drop=True)
        all_exists = True
        for pth in scoll_lays.path.values:
            if not Path(pth).exists():
                print(f"FILE DOES NOT EXIST: {pth}")
                all_exists = False
        if not all_exists:
            raise ValueError("Missing files for scenecoll {scenecoll_name} and bands {','.join(bands)}.")
        return scoll_lays

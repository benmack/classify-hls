import configparser
import sys
from pathlib import Path

class ProjectConfigParser(configparser.ConfigParser):
    def __init__(self, config_file=None):
        """Provides access to the project configuration file and more.
        
        Arguments:
            config_file {str} -- The configuration file, if ``None`` it is assumed to be the file \*project_configs.ini\* in the project root directory.  
        """

        configparser.ConfigParser.__init__(self, interpolation=configparser.ExtendedInterpolation())
        rootdir = sys.modules[__name__].__file__
        rootdir = Path(rootdir).parent.parent
        self._config_file = rootdir / "project_configs.ini"
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

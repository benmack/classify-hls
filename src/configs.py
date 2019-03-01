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

    def get_path(self, section, option):
        """Get an option from the configuration file and convert it to a Path object.
        
        Arguments:
            section {str} -- See ``configparser.ConfigParser.get()`` documentation.
            option {str} -- See ``configparser.ConfigParser.get()`` documentation.
        
        Returns:
            ``Path`` object -- ``pathlib.Path(self.get(section, option))``
        """
        return Path(self.get(section, option))

.. classify-hls documentation master file, created by
   sphinx-quickstart.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

classify-hls documentation!
==============================================

Scope of the project
--------------------

.. todo:: Short description of the project scope.


.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   getting-started
   notebooks/00_basics/00_project_config_parser

.. toctree::
   :maxdepth: 2
   :caption: Raw to Interim

   notebooks/01_raw/01a-r_bm_raw-data.ipynb
   notebooks/01_raw/01b-r_bm_scene-collections.ipynb
   notebooks/01_raw/02-r-eda_bm_raw-data-eda.ipynb
   notebooks/01_raw/03b-r2i-dag_clc-dag.ipynb
   notebooks/01_raw/04a-r2i-dev_bm_hls-dev.ipynb
   notebooks/01_raw/04b-r2i-dag_bm_hls-dag.ipynb
   
.. toctree::
   :maxdepth: 2
   :caption: Interim to Processed

   
   notebooks/02_interim/01-i-eda_bm_interim-data-eda.ipynb
   notebooks/02_interim/02a-i2p-dev_bm_aux-layers.ipynb
   notebooks/02_interim/03a-i2p-dev_bm_features-vts-dev.ipynb
   notebooks/02_interim/04a-i2p-dev_bm_extract-reference-data-dev.ipynb
    

.. toctree::
   :maxdepth: 2
   :caption: Processed
   
   notebooks/03_processed/01a-p-eda_bm_processed-data-eda.ipynb
   notebooks/03_processed/02a_bm_load-extracted-data-from-npy-files.ipynb
   

.. toctree::
   :maxdepth: 2
   :caption: Maybe useful

   notebooks/99_maybe-useful/polygon_to_footprint_border_distance.ipynb


: WARNING: document isn't included in any toctree
: WARNING: document isn't included in any toctree
: WARNING: document isn't included in any toctree


.. toctree::
   :maxdepth: 2
   :caption: src Package

   modules/src

.. toctree::
   :maxdepth: 2
   :caption: ToDo


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

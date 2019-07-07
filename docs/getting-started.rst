Project data and code structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This project borrows ideas from 
`Cookiecutter Data Science <https://drivendata.github.io/cookiecutter-data-science/>`_.

We follow the data structure::

    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.

Also, we try to follow the idea that Jupyter Notebooks are very effective for exploratory data analysis (EDA)
but not for reproducing an analysis.

All the data in the *raw* folder is described and partially generated in the notebook *Raw data*.
The whole project depends on this data and from here on we try to (...) 
run all following preocessing as a DAG (directed acyclic graph) managed with snakemake.

We have the following notebook types.

<step-type>_<user>_<description>.ipynb

with <step-type> being a number and tpye of the following::

    r       : raw data generation

    r-eda   : report-style EDA to see what we have got from *r* 

    r2i-dev : Developement notebook to prepare the following DAG

    r2i-dag : description of the DAG to be run

    i-eda   : report-style EDA to see what we have got by running the *r2i-dag*
    
    i2p     : ...

Typically, it might be a good idea to first create a *?-eda* notebook to see 

* what data is available since the last DAG run
* investigate the nature of the data to see what has to be done next

*?-eda* notebooks are a good place where can store some information about your 
data that you might to look up later from time to time. 

Then we can start from here, duplicate the notebook and work on the next DAG part. 


.. todo:: Installation / requirements to run the project.

============
Installation
============

First install miniforge following the instructions given here https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

User installation
---------------------

OpenAlea Adel may by installed simply on a conda environments:

::

    mamba create -n adel -c conda-forge -c openalea3 openalea.adel
    mamba activate adel

That creates a conda environment called *adel*, install in it *openalea.adel* with all the dependencies and
activate the environment. Then just open an Ipython session and enjoy.

If you want notebook support, run for example:

::

    mamba install jupyterlab

Developer installation
-------------------------

First fork the git repository (https://github.com/openalea/adel) and clone it locally see https://docs.github.com/en/get-started/quickstart/fork-a-repo.

Just run the following command:

::

    mamba create -f conda/environment.yml
    mamba activate adel
    pip install -e .[options]

This will first create a conda environment called *adel_dev* with the proper dependencies, then the environment will be activated,
and finally openalea.adel will be installed in development mode. As above to have notebook support run `mamba install jupyterlab`.
[options] is optional, and allows to install additional dependencies defined in the [project.optional-dependencies] section of your
pyproject.toml file (usually "dev", or "doc", ...)

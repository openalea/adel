# OpenAlea.Adel

**Authors** : C. Fournier, C. Pradal, B. Andrieu

**Contributors** :

-   M Abichou (wheat parameteristion, documentation, plantgen),
-   C Chambon (plantgen, documentation, povray, fit, io widgets)

**Institutes** : INRAe, CIRAD

**Status** : R and Python package (+ L-system)

**License** : Cecill-C

# About

[![Last version](https://anaconda.org/openalea3/openalea.adel/badges/version.svg)](https://anaconda.org/OpenAlea3/openalea.adel/files)
[![Documentation Status](https://readthedocs.org/projects/adel/badge/?version=latest)](https://adel.readthedocs.io/en/latest/?badge=latest)
[![Licence](https://anaconda.org/openalea3/openalea.adel/badges/license.svg)](https://cecill.info/licences/Licence_CeCILL_V2.1-en.html)
[![Platform](https://anaconda.org/openalea3/openalea.adel/badges/platforms.svg)](https://anaconda.org/openalea3/openalea.adel)
[![Downloads](https://anaconda.org/openalea3/openalea.adel/badges/downloads.svg)](https://anaconda.org/openalea3/openalea.adel)

## Description

OpenAlea.Adel (Architectural model of DEvelopment based on L-systems)
allows to simulate the 3D architectural development of the shoot of
gramineaous plant.

## Content

The package hosts generic data structure and simulation tools for
gramineaous plants(Fournier & Pradal, unpublished), the Adel-Maize
(Fournier & Andrieu, 1998), Adel-Wheat (Fournier et al. 2003) models,
together with the wheat parameterization model of Abichou et al. (2013)
and the plastic leaf model of Fournier & Pradal (2012)


## Installation

### Users

```bash
mamba env create -n adel -c openalea3 -c conda-forge openalea.adel 
``` 

### Developers

```bash
git clone 'https://github.com/openalea/adel.git'
cd caribu
mamba env create -n adel_dev -f ./conda/environment.yml
```

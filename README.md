# Analysis Tools for CP2K Simulations

This repository contains Python-based analysis tools to process results from density functional theory simulations. I've used these scripts mainly to process results from simulations with the [CP2K](https://cp2k.org) software, but you can use the scripts with your favorite simulation code as long as you preprocess the data appropriately.

The following scripts are included (see folder [src](https://github.com/nholmber/cp2k-analysis-tools/tree/master/src)):
* [`dos.py`](https://github.com/nholmber/cp2k-analysis-tools/blob/master/src/dos.py)
	* Analyze properties related to the system's density of states (DOS) and electrode potential.
* [`potential.py`](https://github.com/nholmber/cp2k-analysis-tools/blob/master/src/potential.py)
	* Tool to process one- or two-dimensional potential energy surfaces.

You can find Jupyter notebooks in the [examples](https://github.com/nholmber/cp2k-analysis-tools/tree/master/examples) directory, which show you how to use the scripts and highlight their features.

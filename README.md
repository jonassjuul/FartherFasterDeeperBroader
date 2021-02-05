# FartherFasterDeeperBroader

## What does this repository contain?
This repository contains code to reproduce figures and data from the manuscript

"Farther, faster, deeper, broader cascades -- or simply larger?"
by [Jonas L. Juul](https://people.cam.cornell.edu/jsj85/) and [Johan Ugander](https://web.stanford.edu/~jugander/).

To run these scripts, one should acquire the dataset of the spread of false and true news on Twitter. A link to these data can be found in the acknowledgment section of Vosoughi et al. [DOI: 10.1126/science.aap9559](https://science.sciencemag.org/content/359/6380/1146).

Due to space constraints, synthetic cascade data files has not been uploaded to this repository. Instead, we have shared scripts to simulate such data.

## Getting started.
The directory `code` contains the code for the paper. 
The jupyter notebook `main_clean.py` is the main script. This reproduces all figures and tables from the paper. To run this notebook, you will need to 

* Specify the path to the Vosoughi et al. data in the file `definitions_clean.py`. 
* Download data of simulated SIR and IC processes and save the files to the correct subfolders. Instructions are found in the files `download_data_for_this_folder.txt` found in the directories [`IC_data`](https://github.com/jonassjuul/FartherFasterDeeperBroader/tree/main/code/IC_data), [`IC_data/Network_simulations`](https://github.com/jonassjuul/FartherFasterDeeperBroader/tree/main/code/IC_data/Network_simulations), [`SIR_data`](https://github.com/jonassjuul/FartherFasterDeeperBroader/tree/main/code/SIR_data), and [`SIR_data/Network_simulations`](https://github.com/jonassjuul/FartherFasterDeeperBroader/tree/main/code/SIR_data/Network_simulations).
 * Instead of downloading these data, you can also simulate SIR and IC processes on networks by
   * Specifying edge-lists for the empirical networks.
   * Running each `go` file in the folders `IC_data` and `SIR_data` (and subfolders `Network_simulations`)
   * Running all Python files in the folders `IC_data` and `SIR_data` (and subfolders `Network_simulations`).

You can also write me an email to obtain the SIR and IC data files.

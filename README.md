# Title

This repository contains the dataset and code to reproduce results of the following conference paper:

## Abstract of the Paper


## Content of Code Package

This code package is structured as follows:

- `main.py`: Main script to run the simulations.
- `represent_results`: This script plots the figures from the paper.

See each file for further documentation.

# Associated dataset

This repository is associated with a dataset that contains the simulation of 5 distinct cell-free massive MIMO scenarios, where each scenario includes 20,000 setups. Each scenario considered different configurations that are indicated in the name of the variables as '**parameter_precoder**_**cooperation**' where:

* **parameter**: is the estimated variable which values can be; (1) **parameter** = 'R' for the spectral efficiency (SE) calculated as Pareto-optimal rate, (2) **parameter** = 'SE' for the SE refined version, and (3) **parameter** = 'time' for the execution time of that configuration in seconds.
* **precoder**: is the precoding scheme used; (1) **precoder** = 'MMSE' for the minimum mean-squared error (MMSE), (2) **precoder** = 'LTMMSE' for the local team MMSE (LTMMSE), (3) **precoder** = 'UTMMSE' for the unidirectional team MMSE (UTMMSE), and (4) **precoder** = 'LPMMSE' for the local partial MMSE (LPMMSE).
* **cooperation**: is the APs cooperation schemes; (1) **cooperation** = 'all' for all the APs serving all UEs, and (2) **cooperation** = 'DCC' for dynamic cooperation cluster (DCC) scheme.

The dataset is available at

**NOTE:** The downloaded files should be placed in a directory named `results/`.

# Acknowledgments

This work is supported by the Spanish ministry of economic affairs and digital transformation and the European Union - NextGenerationEU [UNICO I+D 6G/INSIGNIA] (TSI-064200-2022-006).

# License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

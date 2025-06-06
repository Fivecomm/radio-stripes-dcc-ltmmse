# Efficient DCC-Enhanced Precoder for Radio Stripes in Cell-Free Massive MIMO

This repository contains the dataset and code to reproduce results of the following conference paper:

G. García-Barrios and M. Fuentes, "Efficient DCC-Enhanced Precoder for Radio Stripes in Cell-Free Massive MIMO,"  *2025 IEEE Conference on Standards for Communications and Networking (CSCN)*, Bologna, Italy, 2025. [Submitted]

## Abstract of the Paper

Future wireless networks will require architectures that can deliver uniformly high service quality, especially in dense deployments and dynamic environments. Cell-free massive multiple-input multiple-output (MIMO) systems have emerged as a promising solution by distributing many access points (APs) across a wide area to jointly serve users without the constraints of traditional cell boundaries. However, their practical deployment poses challenges in terms of cost, scalability, and coordination. To address these issues, the radio stripes architecture has been proposed, enabling low-complexity and scalable implementations through the serial interconnection of simple APs along a stripe. In this work, we address a key challenge in radio stripes: the design of efficient distributed precoding algorithms that can operate under the constraints of limited processing and communication capabilities. We propose a novel solution by applying the dynamic cooperation cluster (DCC) framework to the local team minimum mean-squared error (LTMMSE) precoder, resulting in a highly efficient and scalable design. Our results show that the DCC-enhanced LTMMSE achieves nearly the same spectral efficiency as the centralized baseline while reducing computation time by almost an order of magnitude. Moreover, it outperforms the only other known DCC-based method, the local partial MMSE (LP-MMSE), in both performance and computational complexity. A detailed computational cost analysis further reveals how the AP-to-UE ratio critically impacts system efficiency. These findings demonstrate the practical potential of DCC-enhanced LTMMSE for enabling real-world radio stripe deployments in future cell-free massive MIMO networks.


## Content of Code Package

This code package is structured as follows:

- `main.py`: Main script to run the simulations.
- `represent_results`: This script plots the figures from the paper.

See each file for further documentation.

# Associated dataset

This repository is associated with a dataset that contains the simulation of 4 distinct cell-free massive MIMO scenarios, where each scenario includes 20,000 setups. Each scenario considered different configurations that are indicated in the name of the variables as '**parameter_precoder**_**cooperation**' where:

* **parameter**: is the estimated variable which values can be; (1) **parameter** = 'R' for the spectral efficiency (SE) calculated as Pareto-optimal rate, (2) **parameter** = 'SE' for the SE refined version, and (3) **parameter** = 'time' for the execution time of that configuration in seconds.
* **precoder**: is the precoding scheme used; (1) **precoder** = 'MMSE' for the minimum mean-squared error (MMSE), (2) **precoder** = 'LTMMSE' for the local team MMSE (LTMMSE), (3) **precoder** = 'UTMMSE' for the unidirectional team MMSE (UTMMSE), and (4) **precoder** = 'LPMMSE' for the local partial MMSE (LPMMSE).
* **cooperation**: is the APs cooperation schemes; (1) **cooperation** = 'all' for all the APs serving all UEs, and (2) **cooperation** = 'DCC' for dynamic cooperation cluster (DCC) scheme.

The dataset is available at [https://zenodo.org/records/15598597](https://zenodo.org/records/15598597).

**NOTE:** The downloaded files should be placed in a directory named `results/`.

# Acknowledgments

This work is supported by the Spanish ministry of economic affairs and digital transformation and the European Union - NextGenerationEU [UNICO I+D 6G/INSIGNIA] (TSI-064200-2022-006).

# License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

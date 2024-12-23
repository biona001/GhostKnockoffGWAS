# Downloads page

Here is the main downloads page. New software and pre-processed knockoff data will be released here.

## Software

| Operating System | v0.2.3 (Nov 7th, 2024)      |
| :---             |       :----:               |
| Linux 64-bit     | [Download](https://github.com/biona001/GhostKnockoffGWAS/releases/tag/v0.2.3)       |

After unzipping, the executable will be located inside `bin/GhostKnockoffGWAS`. We recommend adding the folder containing the `GhostKnockoffGWAS` executable to `PATH` for easier access.

## Pre-processed LD files

We welcome `solveblock` users to upload its output to the cloud and share the download link with us. After checking for its quality, we will include it in the following table.

| Population              | Link        | Number of SNPs | Description   |  Citation  |
| :---                    |    :----:   |      :---:     |   :---:     |   :---:    |
| EUR (Europeans)         | [download](https://zenodo.org/records/10433663) (7.5GB)  |650826 | See **Note 1** |  [paper](https://www.biorxiv.org/content/10.1101/2024.02.28.582621v2)   |
| ASN (East Asians)       | TBD        |       |  |
| AFR (Africans)          | TBD        |       |  |
| AMR (Admixed Americans) | TBD        |       |  |  |

+ **Note 1**: This file contain pre-processed LD files generated from the typed SNPs of the EUR cohort from the Pan-UKB panel. The quasi-independent regions were obtained by directly adapting [the output of ldetect](https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/)

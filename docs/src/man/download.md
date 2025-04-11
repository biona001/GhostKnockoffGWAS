# Downloads page

Here is the main downloads page. New software and pre-processed knockoff data will be released here.

## Software

| Operating System | v0.2.3 (Nov 7th, 2024)      |
| :---             |       :----:               |
| Linux 64-bit     | [Download](https://github.com/biona001/GhostKnockoffGWAS/releases/tag/v0.2.3)       |

After unzipping, the executable will be located inside `bin/GhostKnockoffGWAS`. We recommend adding the folder containing the `GhostKnockoffGWAS` executable to `PATH` for easier access.

## Pre-processed LD files

We welcome `solveblock` users to upload its output to the cloud and share the download link with us. 

| Population              | Link        | Number of samples | Number of blocks   |  How were blocks defined | HG Build  |  Citation  |
| :---                    |    :----:   |      :---:     |   :---:     |   :---:    |  :---:    |  :---:    |
| **EUR (British in UK-Biobank)     | [download](https://zenodo.org/records/15191305)  | 306604 | 636 |  Computed with [snp_ldsplit](https://privefl.github.io/bigsnpr/reference/snp_ldsplit.html)  | 19 |  |  
| **EUR (Europeans from Pan-UKB)    | [download](https://zenodo.org/records/10433663)  | See * |  1703 | Adapted from [ldetect](https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/) | 19 and 38  |  [paper](https://www.biorxiv.org/content/10.1101/2024.02.28.582621v2)   |  
| IND (Indians in UK-Biobank)      | [download](https://zenodo.org/records/15191862)  | 5951 | 615 | Computed with [snp_ldsplit](https://privefl.github.io/bigsnpr/reference/snp_ldsplit.html)  | 19 |  | 
| CRB (Caribbeans in UK-Biobank)   | [download](https://zenodo.org/records/15192021)  | 4517 | 489 | Computed with [snp_ldsplit](https://privefl.github.io/bigsnpr/reference/snp_ldsplit.html)  | 19 |  | 
| AFR (Africans in UK-Biobank)     | [download](https://zenodo.org/records/15198591)  |  3394 | 513 | Computed with [snp_ldsplit](https://privefl.github.io/bigsnpr/reference/snp_ldsplit.html)  | 19 |  | 
| CHN (Chinese in UK-Biobank)      | [download](https://zenodo.org/records/15198714)  |  1574 | 505 | Computed with [snp_ldsplit](https://privefl.github.io/bigsnpr/reference/snp_ldsplit.html)  | 19 |  | 

\*: This file contain pre-processed LD files generated from the typed SNPs of the EUR cohort from the Pan-UKB panel. 

\*\*: The difference between the 2 EUR panels mainly resides in how the quasi-independent blocks are computed. All else being equal, the first entry, **EUR (British in UK-Biobank)**, should be preferred. Their main difference is that the in the first, we compute the independent blocks using snp_ldsplit directly on the UK-Biobank British samples. The second one used quasi-independent blocks from LDetect computed on the 1000 genomes project. These block boundaries should be less accurate due to significantly fewer samples. 

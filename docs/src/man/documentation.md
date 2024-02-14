
# Command-line documentation and usage of GhostKnockoffGWAS

After downloading the software, e.g. `app_linux_x86.tar.gz`, unzip the file in linux command line via:
```shell
tar -xvzf app_linux_x86.tar.gz
```
This should create a folder called `app_linux_x86` in the current directory. The executable is located inside `app_linux_x86/bin/GhostKnockoffGWAS`. 

!!! warning

    Do NOT modify the content in unzipped folder. 


## Usage

Simple run

```
GhostKnockoffGWAS --zfile example_zfile.txt --LD-files EUR --N 506200 --genome-build 38 --out example_output
```

## Required inputs

| Option name              | Argument        | Default | Description   |
| :---                    |    :----:   |      :---:     |   :---     |
| `--zfile`        | String | NA | Input file containing Z-scores |
| `--LD-files`     | String | NA | Input directory to the pre-processed LD files |
| `--N`            | Int    | NA | Sample size for target (original) study |
| `--genome-build` | Int    | NA | The human genome build used for SNP positions in `zfile` (this value must be 19 or 38) |
| `--out`          | String | NA | Output file name (without extensions) |

## Optional inputs


| Option name              | Argument        | Default | Description   |
| :---                    |    :----:   |      :---:     |   :---     |
| `--CHR`        | Int    | Whichever column in `zfile` with header `CHR` | The column in `zfile` that will be read as chromosome number (note this must be an integer, e.g. chr22, X, chrX, ...etc are NOT acceptable) |
| `--POS`        | Int    | Whichever column in `zfile` with header `POS` | The column in `zfile` that will be read as SNP position |
| `--REF`        | Int    | Whichever column in `zfile` with header `REF` | The column in `zfile` that will be read as REF (non-effectiv) allele |
| `--ALT`        | Int    | Whichever column in `zfile` with header `ALT` | The column in `zfile` that will be read as ALT (effective allele) |
| `--Z`          | Int    | Whichever column in `zfile` with header `Z` | The column in `zfile` that will be read as Z-scores |
| `--seed`       | Int    | 2023 | Sets the random seed |
| `--verbose`    | Bool   | `true` | Whether to print intermediate messages |
| `--random-shuffle` | Bool | `true` | Whether to randomly permute the order of Z-scores and their knockoffs to adjust for potential ordering bias. |
| `--skip-shrinkage-check` | Bool | `false` | Whether to allow Knockoff analysis to proceed even with large (>0.25) LD shrinkages |


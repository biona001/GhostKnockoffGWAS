
# Command-line documentation and usage of GhostKnockoffGWAS

## Usage

Simple run

```
GhostKnockoffGWAS --zfile example_zfile.txt --LD-files EUR --N 506200 --genome-build 38 --out example_output
```

## Required inputs

| Option name              | Argument        | Description   |
| :---                     |    :----:       |   :---        |
| `--zfile`        | String | Input file containing Z-scores as well as CHR/POS/REF/ALT. See [Acceptable Z-score files](https://biona001.github.io/GhostKnockoffGWAS/dev/man/zfile) for detailed requirement on this file. |
| `--LD-files`     | String | Input directory to the pre-processed LD files. Most users downloads this from the [Downloads Page](https://biona001.github.io/GhostKnockoffGWAS/dev/man/download) |
| `--N`            | Int    | Sample size for target (original) study |
| `--genome-build` | Int    | The human genome build used for SNP positions in `zfile` (this value must be 19 or 38) |
| `--out`          | String | Output file name (without extensions) |

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


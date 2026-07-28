# Installation

The current package layout has `DESCRIPTION` at the repository root.

## GitHub installation

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("marianaboroni/Variants_analysis")
```

Then verify:

```r
library(tumoronly)
tumoronly_default_config()
```

## Local development installation

```bash
git clone https://github.com/marianaboroni/Variants_analysis.git
cd Variants_analysis
R CMD INSTALL .
Rscript -e 'library(tumoronly); check_installation()'
```

## Source-tree CLI

When developing from a checkout, use:

```bash
Rscript exec/tumoronly doctor
Rscript exec/tumoronly init --output config.yaml
```

The CLI loads local `R/*.R` files before falling back to an installed package,
which keeps development runs tied to the checked-out source.

## Conda

Create an environment from `environment.yml`:

```bash
conda env create -f environment.yml
conda activate tumoronly
R CMD INSTALL .
Rscript exec/tumoronly doctor
```

## Docker

Build and run:

```bash
docker build -t tumoronly:dev .
docker run --rm -v "$PWD":/work -w /work tumoronly:dev Rscript exec/tumoronly doctor
```

## Required R packages

Core:

- `data.table`
- `yaml`
- `jsonlite`
- `digest`

Recommended:

- `ggplot2`
- `rmarkdown`
- `knitr`
- `maftools`

Optional specialized packages:

- `GenomicRanges`, `IRanges`, `rtracklayer`, `Biostrings` for reference and
  liftover workflows.
- `httr2` for optional OncoKB post-hoc annotation.

Core tumoronly analysis does not require credentials.


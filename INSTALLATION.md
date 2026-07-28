# Installation

The package root is the repository root.

## GitHub

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("marianaboroni/Variants_analysis")
library(tumoronly)
check_installation()
```

The installed CLI script can be located with:

```bash
TO_CLI="$(Rscript -e 'cat(system.file("exec", "tumoronly", package = "tumoronly"))')"
Rscript "$TO_CLI" init --output /tmp/tumoronly_config.yml --force
```

The CLI currently expects invocation through `Rscript "$TO_CLI"` after package
installation. In a source checkout, use `Rscript exec/tumoronly`.

## Local Development

```bash
git clone https://github.com/marianaboroni/Variants_analysis.git
cd Variants_analysis
Rscript -e 'install.packages(c("data.table","yaml","jsonlite","digest"), repos="https://cloud.r-project.org")'
R CMD INSTALL .
Rscript -e 'library(tumoronly); check_installation()'
```

Recommended for HTML reports and figures:

```bash
Rscript -e 'install.packages(c("ggplot2","rmarkdown","knitr"), repos="https://cloud.r-project.org")'
```

## Conda

```bash
conda env create -f environment.yml
conda activate tumoronly
R CMD INSTALL .
Rscript -e 'library(tumoronly); check_installation()'
```

Run the demo from the checkout:

```bash
Rscript exec/tumoronly validate --input data/demo/variants.tsv --config config/demo_v2.yml
```

## Docker

```bash
docker build -t tumoronly:dev .
docker run --rm -v "$PWD":/work -w /work tumoronly:dev R CMD INSTALL .
docker run --rm -v "$PWD":/work -w /work tumoronly:dev Rscript exec/tumoronly validate --input data/demo/variants.tsv --config config/demo_v2.yml
```

## Dependencies

Core dependencies:

- `data.table`
- `yaml`
- `jsonlite`
- `digest`

Recommended dependencies:

- `ggplot2`
- `rmarkdown`
- `knitr`
- `maftools`

Optional specialized dependencies:

- `GenomicRanges`, `IRanges`, `rtracklayer`, `Biostrings` for reference-backed
  and liftover-related workflows.
- `httr2` for optional post-hoc annotation workflows.

The core v2 analysis does not require credentials.

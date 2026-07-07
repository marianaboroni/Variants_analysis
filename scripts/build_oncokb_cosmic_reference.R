#!/usr/bin/env Rscript
# DEPRECATED. This script mixed four unrelated concerns (COSMIC, OncoKB, driver
# genes, hotspots). They are now decoupled:
#   * COSMIC       -> tumoronly prepare-cosmic --config <config.yml>  (idempotent,
#                     lifted over ONCE, canonical-key, checksummed manifest)
#   * OncoKB       -> tumoronly annotate-oncokb --run-dir <run>       (post-hoc API)
#   * driver genes -> read directly from driver_resources.driver_genes at run time
#   * hotspots     -> read directly from driver_resources.hotspots at run time

stop(paste0(
  "scripts/build_oncokb_cosmic_reference.R is deprecated.\n",
  "Use: tumoronly prepare-cosmic --config <config.yml> to build the local COSMIC DB.\n",
  "OncoKB is now a separate step: tumoronly annotate-oncokb --run-dir <run>."), call. = FALSE)

test_that("export_review_template pre-fills current labels and blank review cols", {
  v <- data.frame(sample_id = "S1", variant_id = "1:100:A:T", gene = "TP53",
    consequence = "missense_variant", STATUS = "PASS", PRIORITY = "HIGH",
    AUTHENTICITY = "LIKELY_TRUE", VARIANT_REVIEW_SCORE = 0.3, stringsAsFactors = FALSE)
  p <- tempfile(fileext = ".tsv"); export_review_template(v, p)
  t <- read_variants(p, "\t")
  expect_equal(t$review_label, "NOT_REVIEWED")
  expect_true(all(c("reviewer", "review_reason", "validation_method") %in% names(t)))
})

test_that("import_variant_reviews is append-only and never converts PASS/FAIL to labels", {
  db <- tempfile("ev")
  mk_reviews <- function(labels) {
    d <- data.frame(sample_id = "S1", variant_id = sprintf("1:%d:A:T", seq_along(labels)),
      review_label = labels, reviewer = "curator1", review_reason = "orthogonal validation",
      stringsAsFactors = FALSE)
    f <- tempfile(fileext = ".tsv"); write_tsv(d, f); f
  }
  r1 <- suppressMessages(import_variant_reviews(mk_reviews(c("TRUE_POSITIVE", "FALSE_POSITIVE", "NOT_REVIEWED")), db_dir = db))
  expect_equal(r1$n_appended, 2)   # NOT_REVIEWED rows are not ingested
  hist1 <- read_variants(file.path(db, "review_history.tsv.gz"), "\t")
  r2 <- suppressMessages(import_variant_reviews(mk_reviews(c("UNCERTAIN")), db_dir = db))
  hist2 <- read_variants(file.path(db, "review_history.tsv.gz"), "\t")
  expect_gt(nrow(hist2), nrow(hist1))   # append-only history grows
  expect_true(file.exists(file.path(db, "feature_schema.json")))
  expect_true(file.exists(file.path(db, "label_dictionary.yml")))
})

test_that("TRUE_POSITIVE without a reviewer is rejected; invalid labels are rejected", {
  db <- tempfile("ev")
  d1 <- data.frame(variant_id = "1:1:A:T", review_label = "TRUE_POSITIVE", stringsAsFactors = FALSE)
  f1 <- tempfile(fileext = ".tsv"); write_tsv(d1, f1)
  expect_error(import_variant_reviews(f1, db_dir = db), "reviewer", ignore.case = TRUE)
  d2 <- data.frame(variant_id = "1:1:A:T", review_label = "SOMATIC", reviewer = "x", stringsAsFactors = FALSE)
  f2 <- tempfile(fileext = ".tsv"); write_tsv(d2, f2)
  expect_error(import_variant_reviews(f2, db_dir = db), "Invalid review_label", ignore.case = TRUE)
})

test_that("ML scaffold: build training set, train candidate (not active), gated activation", {
  db <- tempfile("ev"); ensure_variant_evidence_db(db)
  set.seed(1); n <- 60
  samples <- paste0("S", rep(1:6, each = 10))
  y <- rep(c("TRUE_POSITIVE", "FALSE_POSITIVE"), length.out = n)
  # a run_dir with features
  rd <- tempfile("run"); dir.create(file.path(rd, "tables"), recursive = TRUE)
  va <- data.frame(sample_id = samples, variant_id = sprintf("1:%d:A:T", 1:n),
    dp = ifelse(y == "TRUE_POSITIVE", 80, 20) + rnorm(n, 0, 5),
    alt_count = ifelse(y == "TRUE_POSITIVE", 38, 4) + rnorm(n, 0, 2),
    vaf = ifelse(y == "TRUE_POSITIVE", 0.45, 0.05), tlod = ifelse(y == "TRUE_POSITIVE", 60, 5),
    mbq = 35, mmq = 60, strand_artifact = ifelse(y == "TRUE_POSITIVE", 0.1, 0.7),
    orientation_bias = 0.1, pon_flag = y == "FALSE_POSITIVE",
    clustered_events = FALSE, weak_evidence = y == "FALSE_POSITIVE", stringsAsFactors = FALSE)
  write_tsv(va, file.path(rd, "tables", "variants_all.tsv.gz"))
  rev <- data.frame(sample_id = samples, variant_id = va$variant_id, review_label = y,
    reviewer = "curator", stringsAsFactors = FALSE)
  rf <- tempfile(fileext = ".tsv"); write_tsv(rev, rf)
  suppressMessages(import_variant_reviews(rf, db_dir = db))

  st <- suppressMessages(ml_training_status(db))
  expect_true(st$ready)
  suppressMessages(ml_build_training_set(db, rd, "P_TRUE_VARIANT"))
  res <- suppressMessages(ml_train_model(db, "P_TRUE_VARIANT"))
  expect_true(!is.na(res$metrics$roc_auc))
  # no ACTIVE model yet
  expect_false(file.exists(file.path(db, "model_registry", "ACTIVE")))
  # activation requires approval
  expect_error(ml_activate_model(db, res$model_id, approve = FALSE), "approval", ignore.case = TRUE)
  suppressMessages(ml_activate_model(db, res$model_id, approve = TRUE))
  expect_true(file.exists(file.path(db, "model_registry", "ACTIVE")))

  cfg <- tumoronly_default_config()
  cfg$ml$db_dir <- db
  va$final_class <- ifelse(y == "TRUE_POSITIVE", "probable_somatic", "likely_artifact")
  va$filter_status <- ifelse(y == "TRUE_POSITIVE", "PASS", "FAIL")
  pred <- apply_ml_predictions(va, cfg)
  expect_equal(pred$status$status, "evaluated")
  expect_true(any(!is.na(pred$variants$ml_true_positive_probability)))
  expect_equal(pred$variants$final_class, va$final_class)
  expect_equal(pred$variants$filter_status, va$filter_status)
})

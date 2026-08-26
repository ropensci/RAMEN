test_that("selectVariables output structure is correct", {
  expect_true(is.data.frame(selected_variables_test))
  expect_equal(ncol(selected_variables_test), 3)
  expect_equal(nrow(selected_variables_test), length(VML_cis_snps_test))
  expect_true(all(
    c("VML_index", "selected_genot", "selected_env") %in%
      colnames(selected_variables_test)
  ))
})

# Test that errors happen when expected
test_that("selectVariables throws errors when expected", {
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = "a",
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the input VML_wSNPs belongs to the GRanges class.",
    fixed = TRUE
  )
  # Test error when there are ID mismatches
  test_genot <- RAMEN::test_genotype_matrix
  colnames(test_genot) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = test_genot,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "The objects rownames(summarized_methyl_VML) and colnames(genotype_matrix) must match.",
    fixed = TRUE
  )
  # Test error when there is argument mismatch with the environmental_matrix
  test_env <- RAMEN::test_environmental_matrix
  rownames(test_env) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "The objects rownames(summarized_methyl_VML) and rownames(environmental_matrix) must match.",
    fixed = TRUE
  )
  # Test error when there is argument mismatch with the covariates
  test_cov <- RAMEN::test_covariates
  rownames(test_cov) <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = test_cov,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "The objects rownames(summarized_methyl_VML) and rownames(covariates) must match.",
    fixed = TRUE
  )

  # Test that matrix arguments throw errors if input is not a matrix
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = as.data.frame(RAMEN::test_environmental_matrix),
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the input environmental_matrix belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = as.data.frame(RAMEN::test_genotype_matrix),
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the input genotype_matrix belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = as.data.frame(RAMEN::test_covariates),
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the input covariates belongs to the matrix class.",
    fixed = TRUE
  )
  # Test missing columns in VML_df
  VML_no_SNP = VML_cis_snps_test
  VML_no_SNP$SNP <- NULL
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_no_SNP,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "The object S4Vectors::mcols(VML_wSNPs) does not have the required columns: VML_index, SNP .",
    fixed = TRUE
  )
  # Test missing values in genotype matrix
  # Introduce NA values
  test_genot_na <- RAMEN::test_genotype_matrix
  test_genot_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = test_genot_na,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the object genotype_matrix contains only finite numeric values (i.e., no NA, NaN or Inf)",
    fixed = TRUE
  )
  # Test missing values in environmental matrix
  # Introduce Inf value
  test_env_na <- RAMEN::test_environmental_matrix
  test_env_na[1, 1] <- Inf
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = test_env_na,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the object environmental_matrix contains only finite numeric values (i.e., no NA, NaN or Inf)",
    fixed = TRUE
  )
  # Test missing values in covariates matrix
  # Introduce NA values
  test_cov_na <- RAMEN::test_covariates
  test_cov_na[1, 1] <- NaN
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = test_cov_na,
      summarized_methyl_VML = summarized_methyl_VML_test
    ),
    "Please make sure the object covariates contains only finite numeric values (i.e., no NA, NaN or Inf)",
    fixed = TRUE
  )
  # Test missing values in summarized methylation VML
  test_summeth_na <- summarized_methyl_VML_test
  test_summeth_na[1, 1] <- NA
  expect_error(
    RAMEN::selectVariables(
      VML_wSNPs = VML_cis_snps_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      summarized_methyl_VML = test_summeth_na
    ),
    "Please make sure the object summarized_methyl_VML contains only finite numeric values (i.e., no NA, NaN or Inf)",
    fixed = TRUE
  )
})

# Test that a VML with a single candidate SNP does not throw an error
test_that("selectVariables does not error when a VML has only one candidate SNP", {
  # glmnet::cv.glmnet() requires a predictor matrix with 2+ columns. VML1 in
  # the test fixture has exactly one candidate SNP, so running it without
  # covariates (which would otherwise pad genot_VMLi to 2+ columns) used to
  # throw "x should be a matrix with 2 or more columns".
  VML_one_snp <- VML_cis_snps_test[VML_cis_snps_test$VML_index == "VML1"]
  expect_equal(length(VML_one_snp$SNP[[1]]), 1)

  result <- expect_no_error(
    suppressWarnings(RAMEN::selectVariables(
      VML_wSNPs = VML_one_snp,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = NULL,
      summarized_methyl_VML = summarized_methyl_VML_test,
      seed = 1
    ))
  )

  # The lone candidate SNP has nothing to be selected against, so it should
  # be kept as-is
  expect_true(
    VML_one_snp$SNP[[1]] %in% result$selected_genot[[1]]
  )
})

# Check that the LASSO stage is wired up correctly, by running the same fits
# outside of selectVariables() and comparing what they select.
# The reference fits are wrapped in a one-task %dorng% loop started from the
# same seed, so that they draw the same random numbers, and therefore build the
# same cross-validation folds, as the fits inside selectVariables(). Comparing
# against a plain cv.glmnet() call would compare different fold assignments.
test_that("the variables LASSO selects are the ones selectVariables reports", {
  foreach::registerDoSEQ()
  set.seed(42)
  n <- 100
  ids <- paste0("S", seq_len(n))
  genot <- matrix(rbinom(6 * n, 2, 0.5), nrow = 6,
                  dimnames = list(paste0("g", 1:6), ids))
  env <- matrix(rnorm(3 * n), ncol = 3,
                dimnames = list(ids, paste0("e", 1:3)))
  methyl <- matrix(2 * genot["g1", ] + env[, "e1"] + rnorm(n, sd = 0.5),
                   ncol = 1, dimnames = list(ids, "VML1"))
  vml <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = 1000, width = 100),
    VML_index = "VML1", SNP = I(list(paste0("g", 1:6)))
  )

  result <- suppressMessages(RAMEN::selectVariables(
    VML_wSNPs = vml,
    genotype_matrix = genot,
    environmental_matrix = env,
    covariates = NULL,
    summarized_methyl_VML = methyl,
    seed = 1
  ))

  # The variables kept by one fit, dropping the intercept
  selected_by <- function(fit) {
    coefficients <- stats::coef(fit, s = "lambda.min")
    names(coefficients[abs(coefficients[, 1]) > 0, ])[-1]
  }
  x_genot <- t(genot[paste0("g", 1:6), ids, drop = FALSE])
  x_env <- env[ids, , drop = FALSE]
  x_joint <- cbind(x_genot, x_env)

  k <- NULL # To avoid R CMD check notes
  `%dorng%` <- doRNG::`%dorng%`
  set.seed(1)
  reference <- foreach::foreach(k = 1) %dorng% {
    list(
      genot = selected_by(glmnet::cv.glmnet(
        x = x_genot, y = methyl, alpha = 1, nfolds = 5,
        penalty.factor = rep(1, ncol(x_genot)))),
      env = selected_by(glmnet::cv.glmnet(
        x = x_env, y = methyl, alpha = 1, nfolds = 5,
        penalty.factor = rep(1, ncol(x_env)))),
      joint = selected_by(glmnet::cv.glmnet(
        x = x_joint, y = methyl, alpha = 1, nfolds = 5,
        penalty.factor = rep(1, ncol(x_joint))))
    )
  }
  reference <- reference[[1]]

  # LASSO has to actually select something, or the comparison below is empty
  expect_gt(length(reference$genot), 0)
  expect_gt(length(reference$env), 0)
  # The SNP and the environmental variable the VML was simulated from are among
  # the variables LASSO keeps
  expect_true("g1" %in% reference$genot)
  expect_true("e1" %in% reference$env)

  # selectVariables() reports exactly the union the two relevant fits produce
  expect_setequal(
    result$selected_genot[[1]],
    setdiff(unique(c(reference$genot, reference$joint)), colnames(x_env))
  )
  expect_setequal(
    result$selected_env[[1]],
    setdiff(unique(c(reference$env, reference$joint)), colnames(x_genot))
  )
})

test_that("selectVariables runs on the genotype alone when environmental_matrix is NULL", {
  expect_no_error(suppressWarnings(RAMEN::selectVariables(VML_wSNPs = VML_cis_snps_test[1:3],
                                         genotype_matrix = RAMEN::test_genotype_matrix,
                                         environmental_matrix = NULL,
                                         covariates = RAMEN::test_covariates,
                                         summarized_methyl_VML = summarized_methyl_VML_test)))
  }
  )

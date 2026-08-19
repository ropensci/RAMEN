test_that("lmGE output structure is correct", {
  foreach::registerDoSEQ()
  lmge_res <- RAMEN::lmGE(
    selected_variables = selected_variables_test[1:5, ],
    summarized_methyl_VML = summarized_methyl_VML_test,
    genotype_matrix = RAMEN::test_genotype_matrix,
    environmental_matrix = RAMEN::test_environmental_matrix,
    covariates = RAMEN::test_covariates,
    model_selection = "AIC"
  )
  expect_true(is.data.frame(lmge_res))
  expect_equal(ncol(lmge_res), 13)
  expect_equal(nrow(lmge_res), 5)
})

test_that("lmGE works with BIC", {
  foreach::registerDoSEQ()
  lmge_res <- RAMEN::lmGE(
    selected_variables = selected_variables_test[1:5, ],
    summarized_methyl_VML = summarized_methyl_VML_test,
    genotype_matrix = RAMEN::test_genotype_matrix,
    environmental_matrix = RAMEN::test_environmental_matrix,
    covariates = RAMEN::test_covariates,
    model_selection = "BIC"
  )
  expect_true(is.data.frame(lmge_res))
  expect_equal(ncol(lmge_res), 13)
  expect_equal(nrow(lmge_res), 5)
})

test_that("lmGE works when there are no covariates", {
  foreach::registerDoSEQ()
  lmge_res <- RAMEN::lmGE(
    selected_variables = selected_variables_test[1:5, ],
    summarized_methyl_VML = summarized_methyl_VML_test,
    genotype_matrix = RAMEN::test_genotype_matrix,
    environmental_matrix = RAMEN::test_environmental_matrix,
    covariates = NULL,
    model_selection = "AIC"
  )
  expect_true(is.data.frame(lmge_res))
  expect_equal(nrow(lmge_res), 5)
  expect_false(any(is.na(lmge_res$tot_r_squared)))

  # VML4's winning model only involves one variable (E). With no covariates
  # to pad out the model, this used to break relaimpo::calc.relimp.lm(),
  # which requires a model with 2+ regressors
  vml_e <- lmge_res[lmge_res$model_group == "E", ][1,]
  expect_false(is.na(vml_e$e_r_squared))
  expect_true(is.na(vml_e$g_r_squared))
})

test_that("lmGE throws errors when expected", {
  expect_error(
    RAMEN::lmGE(
      selected_variables = "a",
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = "AIC"
    ),
    "Please make sure the input selected_variables belongs to the data.frame class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = "summarized_methyl_VML_test",
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = "AIC"
    ),
    "Please make sure the input summarized_methyl_VML belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = "test_genotype_matrix",
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = "AIC"
    ),
    "Please make sure the input genotype_matrix belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = "test_environmental_matrix",
      covariates = RAMEN::test_covariates,
      model_selection = "AIC"
    ),
    "Please make sure the input environmental_matrix belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = "test_covariates",
      model_selection = "AIC"
    ),
    "Please make sure the input covariates belongs to the matrix class.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = "a"
    ),
    "Please make sure the input model_selection is one of the following options: AIC, BIC .",
    fixed = TRUE
  )
  expect_error(
    RAMEN::lmGE(
      selected_variables = selected_variables_test[1:5, ],
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = RAMEN::test_covariates,
      model_selection = 1
    ),
    "Please make sure the input model_selection belongs to the character class.",
    fixed = TRUE
  )
})

test_that("lmGE handles non-syntactic covariate names", {
  cov_nonsyntactic <- RAMEN::test_covariates
  colnames(cov_nonsyntactic)[1] <- "cell type"
  expect_false(make.names(colnames(cov_nonsyntactic)[1]) ==
                 colnames(cov_nonsyntactic)[1])

  result <- expect_no_error(
    suppressWarnings(RAMEN::lmGE(
      selected_variables = selected_variables_test,
      summarized_methyl_VML = summarized_methyl_VML_test,
      genotype_matrix = RAMEN::test_genotype_matrix,
      environmental_matrix = RAMEN::test_environmental_matrix,
      covariates = cov_nonsyntactic,
      model_selection = "AIC"
    ))
  )
})

# A small, fully specified example whose correct answer is known by
# construction: VML_A's methylation is simulated from the SNP g1 and VML_B's
# from the environmental variable e1, so the winning model must be G for the
# first and E for the second. lmGE() matches every input by name rather than by
# position, so reordering the samples or the VML must not change the result.
test_that("lmGE matches inputs by name, not by position", {
  set.seed(42)
  n <- 60
  ids <- paste0("S", seq_len(n))
  g1 <- rbinom(n, 2, 0.5)
  g2 <- rbinom(n, 2, 0.5)
  e1 <- rnorm(n)
  e2 <- rnorm(n)
  c1 <- rnorm(n)
  meth_A <- 0.80 * g1 + 0.05 * c1 + rnorm(n, sd = 0.10) # driven by G
  meth_B <- 0.80 * e1 + 0.05 * c1 + rnorm(n, sd = 0.10) # driven by E

  genot <- rbind(g1 = g1, g2 = g2)
  colnames(genot) <- ids
  env <- cbind(e1 = e1, e2 = e2)
  rownames(env) <- ids
  covs <- cbind(c1 = c1)
  rownames(covs) <- ids
  methyl <- cbind(VML_A = meth_A, VML_B = meth_B)
  rownames(methyl) <- ids

  sel <- data.frame(VML_index = c("VML_A", "VML_B"))
  sel$selected_genot <- list(c("g1", "g2"), c("g1", "g2"))
  sel$selected_env <- list(c("e1", "e2"), c("e1", "e2"))

  # Return the two VML in a fixed order so the assertions do not depend on the
  # order in which lmGE() happens to emit its rows
  fit <- function(methyl, genot, env, covs) {
    res <- suppressWarnings(RAMEN::lmGE(
      selected_variables = sel,
      summarized_methyl_VML = methyl,
      genotype_matrix = genot,
      environmental_matrix = env,
      covariates = covs,
      model_selection = "AIC"
    ))
    res[match(c("VML_A", "VML_B"), res$VML_index), ]
  }

  # The winning model recovers how each VML was simulated
  base <- fit(methyl, genot, env, covs)
  expect_equal(base$model_group, c("G", "E"))
  expect_equal(base$variables[[1]], "g1")
  expect_equal(base$variables[[2]], "e1")

  # Shuffling the VML columns, and the sample order of every other input,
  # must leave the result unchanged
  set.seed(3)
  shuffled <- fit(
    methyl[sample(n), c(2, 1), drop = FALSE],
    genot[, sample(n), drop = FALSE],
    env[sample(n), , drop = FALSE],
    covs[sample(n), , drop = FALSE]
  )
  expect_equal(shuffled$model_group, c("G", "E"))
  expect_equal(shuffled$variables[[1]], "g1")
  expect_equal(shuffled$variables[[2]], "e1")
  expect_equal(shuffled$tot_r_squared, base$tot_r_squared)

  # Sensitivity check: swapping the two VML labels points each VML_index at the
  # other one's methylation, which must flip the winning models. Without this,
  # the assertions above could still pass if VML were matched positionally.
  mislabelled <- methyl
  colnames(mislabelled) <- c("VML_B", "VML_A")
  expect_equal(fit(mislabelled, genot, env, covs)$model_group, c("E", "G"))
})

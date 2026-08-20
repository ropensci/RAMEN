#' Fit linear G, E, G+E and GxE models and select the winning model
#'
#' For a set of Variable Methylated Loci (VML), this function fits a set of
#' genotype (G), environment (E), pairwise additive (G + E) or pairwise
#' interaction (G x E) models, one variable at a time, and selects the best
#' fitting one. Additional information for each winning model is provided,
#' such as its R2, its R2 increase comparing it to a basal model (i.e., a model
#' only fitted with the concomitant variables), the delta AIC/BIC to the next
#' best model from a different category, and the explained variance decomposed
#' for the G, E and GxE components (when applicable). If a VML has no variables
#' selected in the selected_variables object, it will be returned with "B"
#' (basal) as the best model (interpreted as no G or E associated effect).
#' For guidance on interpretation, please build and read the package's vignette.
#'
#' This function supports parallel computing for increased speed. To do so, you
#' have to set the parallel backend in your R session before running the
#' function (e.g., *doParallel::registerDoParallel(4)*)). After that, the
#' function can be run as usual. It is recommended to also set
#' options(future.globals.maxSize= +Inf).
#'
#' For each VML, this function computes a set of models using the variables
#' indicated in the selected_variables object. From the indicated G and E
#' variables, lmGE() fits four groups of models:
#'  - G: Genetics model - fitted one SNP at a time.
#'  - E: Environmental model - fitted one environmental variable at a time.
#'  - G+E: Additive model - fitted for each pairwise combination of G and E
#'  variables indicated in selected_variables.
#'  - GxE: Interaction model - fitted for each pairwise combination of G and E
#'  variables indicated in selected_variables.
#'
#'  These models are fit only if the VML has G or E variables in the
#'  selected_variables object. If a VML does not have neither G nor E variables,
#'   that VML will be ignored and will be returned in the output object with
#'   "B" (baseline) as the best explanatory model.
#'
#' **Model selection**
#'
#' Following the model fitting stage, the best model **per group** is selected
#' using Akaike Information Criterion (AIC) or Bayesian Information Criterion
#' (BIC). Both of these metrics are statistical approaches to select the best
#' model in the same data set, and they have strengths and limitations that make
#'  them excel in different situations. We recommend using AIC because BIC
#'  assumes that the true model is in the set of compared models. Since this
#'  function fits models with individual variables, and we assume that DNAme
#'  variability is more likely to be influenced by more than one single
#'  SNP/environmental exposure at a time, we hypothesize that in most cases,
#'  the true model will not be in the set of compared models. Also, AIC excels
#'  in situations where all models in the model space are "incomplete", and AIC
#'  is preferentially used in cases where the true underlying function is
#'  unknown and our selected model could belong to a very large class of
#'  functions where the relationship could be pretty complex. It is worth
#'  mentioning however that, both metrics tend to pick the same model in a large
#'   number of scenarios. We suggest the users to read Arijit
#'   Chakrabarti & Jayanta K. Ghosh, 2011 for further information about the
#'   difference between these metrics.
#'
#' After selecting the best model per group (G,E,G+E pr GxE), the model with the
#'  lowest AIC or BIC is declared as the winning model. The delta AIC/BIC and
#'  difference of R2 is computed relative to the model with the second lowest
#'  AIC/BIC (i.e., the best model from a different group to the winning one),
#'   and reported in the final object.
#'
#' **Analysis of variance and variance decomposition**
#'
#' Finally, the variance is decomposed and the relative R2 contribution of each
#' of the variables of interest (G, E and GxE) is reported. This decomposition
#' is done using the relaimpo R package, using the Lindeman, Merenda and Gold
#' (lmg) method, which is based on the heuristic approach of averaging the
#' relative R contribution of each variable over all input orders in the linear
#' model. The estimation of the partitioned R2 of each factor in the models
#' was conducted keeping the covariates always in the model as first entry
#' (i.e., the variables specified in covariates did not change order). For
#' further information, we suggest the users to read the documentation and
#' publication of the relaimpo R package (Grömping, 2006).
#'
#' @param selected_variables A data frame obtained with *RAMEN::selectVariables()*.
#' This data frame must contain three columns: 'VML_index' with characters of
#' an unique ID of each VML; ´selected_genot' and 'selected_env' with the SNPs
#' and environmental variables, respectively, that will be used for fitting the
#' genotype (G), environment (E), additive (G + E) or interaction (G x E)
#' models. The columns 'selected_env' and 'selected_genot' must contain lists
#' as elements; VML with no environmental or genotype selected variables
#' must contain an empty list (i.e., list(NULL), list(NA), list("") or
#' list(character(0)) ).
#' @param model_selection Which metric to use to select the best model for each
#' VML. Supported options are "AIC" or BIC".
#' @inheritParams selectVariables
#' @return A data frame with the following columns:
#'  - VML_index: The unique ID of the VML
#'  - model_group: The group to which the winning model belongs to (i.e., G, E,
#'  G+E or GxE)
#'  - variables: The variable(s) that are present in the winning model
#'  (excluding the covariates, which are included in all the models)
#'  - tot_r_squared: R squared of the winning model
#'  - g_r_squared: Estimated R2 allocated to the G in the winning model, if
#'  applicable.
#'  - e_r_squared: Estimated R2 allocated to the E in the winning model, if
#'  applicable.
#'  - gxe_r_squared: Estimated R2 allocated to the interaction in the winning
#'  model (GxE), if applicable.
#'  - AIC/BIC: AIC or BIC metric from the best model in each VML (depending on
#'  the option specified in the argument model_selection).
#'  - second_winner: The second group that possesses the next best model after
#'  the winning one (i.e., G, E, G+E or GxE). This column may have NA if the
#'  variables in selected_variables correspond only to one group (G or E), so
#'  that there is no other model groups to compare to.
#'  - delta_aic/delta_bic: The difference of AIC or BIC value (depending on the
#'  option specified in the argument model_selection) of the winning model and
#'  the best model from the second_winner group (i.e., G, E, G+E or GxE). This
#'  column may have NA if the variables in selected_variables correspond only to
#'   one group (G or E), so that there is no other groups to compare to.
#'  - delta_r_squared: The R2 of the winning model - R2 of the second winner
#'  model. This column may have NA if the variables in selected_variables
#'  correspond only to one group (G or E), so that there is no other groups to
#'  compare to.
#'  - basal_AIC/basal_BIC: AIC or BIC of the basal model (i.e., model fitted
#'  only with the concomitant variables specified in the *covariates* argument)
#'  - basal_rsquared: The R2 of the basal model (i.e., model fitted only with
#'  the concomitant variables specified in the *covariates* argument)
#' @importFrom foreach %dopar%
#' @export
#' @examples
#' # Evaluate sequentially
#' foreach::registerDoSEQ()
#' ## Find VML in test data
#' VML <- RAMEN::findVML(
#'   methylation_data = RAMEN::test_methylation_data,
#'   array_manifest = "IlluminaHumanMethylationEPICv1",
#'   cor_threshold = 0,
#'   var_method = "variance",
#'   var_distribution = "ultrastable",
#'   var_threshold_percentile = 0.99,
#'   max_distance = 1000
#' )
#' ## Find cis SNPs around VML
#' VML_with_cis_snps <- RAMEN::findCisSNPs(
#'   # Use only 5 for demonstration purposes
#'   VML = VML$VML [1:5, ],
#'   genotype_information = RAMEN::test_genotype_information,
#'   distance = 1e6
#' )
#'
#' ## Summarize methylation levels in VML
#' summarized_methyl_VML <- RAMEN::summarizeVML(
#'   methylation_data = RAMEN::test_methylation_data,
#'   VML = VML_with_cis_snps
#' )
#'
#' ## Select relevant genotype and environmental variables
#' selected_vars <- RAMEN::selectVariables(
#'   VML_wSNPs = VML_with_cis_snps,
#'   genotype_matrix = RAMEN::test_genotype_matrix,
#'   environmental_matrix = RAMEN::test_environmental_matrix,
#'   covariates = RAMEN::test_covariates,
#'   summarized_methyl_VML = summarized_methyl_VML,
#'   seed = 1
#' )
#'
#' ## Fit G, E, G+E and GxE models and select the winning one
#' lmge_res <- RAMEN::lmGE(
#'   selected_variables = selected_vars,
#'   summarized_methyl_VML = summarized_methyl_VML,
#'   genotype_matrix = RAMEN::test_genotype_matrix,
#'   environmental_matrix = RAMEN::test_environmental_matrix,
#'   covariates = RAMEN::test_covariates,
#'   model_selection = "AIC"
#' )
#'
lmGE <- function(selected_variables,
                 summarized_methyl_VML,
                 genotype_matrix,
                 environmental_matrix,
                 covariates = NULL,
                 model_selection = "AIC") {
  #### Binding of variables used within the tidyverse framework ####
  selected_env <- selected_genot <- i <- SNP <- env <- model_group <- NULL
  AIC <- BIC <- tot_r_squared <- VML_index <- variables <- g_r_squared <- NULL
  e_r_squared <- gxe_r_squared <- second_winner <- delta_aic <- NULL
  delta_r_squared <- basal_AIC <- basal_rsquared <- delta_bic <- basal_BIC <- NULL

  #### Check arguments ####
  # Input has the right structure
  argument_check(selected_variables, "data.frame")
  columns_exist(selected_variables, c("VML_index",
                                      "selected_genot",
                                      "selected_env"))
  argument_check(selected_variables$selected_env, "list")
  argument_check(selected_variables$selected_genot, "list")
  argument_check(selected_variables$VML_index, "character")
  argument_check(summarized_methyl_VML, "matrix")
  argument_check(genotype_matrix, "matrix")
  argument_check(environmental_matrix, "matrix")
  if (!is.null(covariates)) argument_check(covariates, "matrix")
  argument_char_options(object = model_selection, options = c("AIC", "BIC"))
  # All objects have matching IDs
  vectors_match(rownames(summarized_methyl_VML), colnames(genotype_matrix),
               object_1_name = "rownames(summarized_methyl_VML)",
               object_2_name = "colnames(genotype_matrix)")
  vectors_match(rownames(summarized_methyl_VML), rownames(environmental_matrix),
               object_1_name = "rownames(summarized_methyl_VML)",
               object_2_name = "rownames(environmental_matrix)")
  if (!is.null(covariates)) {
    vectors_match(rownames(summarized_methyl_VML), rownames(covariates),
                 object_1_name = "rownames(summarized_methyl_VML)",
                 object_2_name = "rownames(covariates)")
  }
  # Matrices have only finite numeric values
  finite_numeric_check(genotype_matrix)
  finite_numeric_check(environmental_matrix)
  finite_numeric_check(summarized_methyl_VML)
  if (!is.null(covariates)) finite_numeric_check(covariates)

  #### Select best winning model ####
  # Filter VML that have no selected G and no selected E
  no_vars_VML <- selected_variables |>
    dplyr::filter((selected_env %in% empty_lists &
      selected_genot %in% empty_lists))
  selected_variables <- selected_variables |>
    dplyr::filter(!(selected_env %in% empty_lists &
      selected_genot %in% empty_lists))
  # Create vectors that will be passed to the foreach loop
  sample_ids <- rownames(summarized_methyl_VML) # Get samples in DNAme object
  # Columns of selected_variables, accessed positionally
  sv_ids <- selected_variables$VML_index
  sv_genot <- selected_variables$selected_genot
  sv_env <- selected_variables$selected_env
  # Get the position of each VML in DNAme object taking select_variables as
  # reference
  methyl_col <- match(sv_ids, colnames(summarized_methyl_VML))
  if (anyNA(methyl_col)) { # Check that all IDs are present
    missing_ids <- sv_ids[is.na(methyl_col)]
    stop(paste("Please make sure every VML_index in selected_variables has a",
               "matching column in summarized_methyl_VML. Missing:",
               paste(missing_ids[seq_len(min(5, length(missing_ids)))],
                     collapse = ", ")))
  }
  # Get the positions of each individual in environmental matrix taking DNAme
  # as reference
  env_row <- match(sample_ids, rownames(environmental_matrix))
  # Identify VML with no genotype variables selected
  genot_empty <- vapply(seq_along(sv_genot),
                        function(j) sv_genot[j] %in% empty_lists,
                        logical(1))
  # Get an unlisted vector of all SNPs in the selected_variables object in order
  flat_snps <- unlist(sv_genot[!genot_empty], use.names = FALSE)
  used_snps <- unique(flat_snps) # get unique SNPs used across all VML
  # Get the positons of SNPs used in genotype matrix
  used_rows <- match(used_snps, rownames(genotype_matrix))
  if (anyNA(used_rows)) { # If a used SNP is not present in genotype_matrix throw error
    missing_snps <- used_snps[is.na(used_rows)]
    stop(paste("Please make sure every SNP in selected_variables$selected_genot",
               "is present in the row names of genotype_matrix. Missing:",
               paste(missing_snps[seq_len(min(5, length(missing_snps)))],
                     collapse = ", ")))
  }
  # Create a smaller matrix with used genotypes when worth it
  if (length(used_snps) < 0.7 * nrow(genotype_matrix)) {
    genotype_matrix <- genotype_matrix[used_rows, , drop = FALSE]
  }
  # Get the positions of each individual in genotype_matrix taking DNAme
  # as reference
  genot_col <- match(sample_ids, colnames(genotype_matrix))
  # Get the positions of each selected SNP in genotype matrix in VML order
  flat_rows <- match(flat_snps, rownames(genotype_matrix))
  # Separate the order per VML
  sv_genot_rows <- vector("list", length(sv_genot))
  sv_genot_rows[!genot_empty] <- split(
    flat_rows,
    rep(seq_len(sum(!genot_empty)), lengths(sv_genot[!genot_empty]))
  )

  # Create a covariate formula, which is shared across all VML
  if (is.null(covariates)) {
    covariates_i <- NULL
    basal_model_formula <- "1"
  } else {
    covariates_i <- covariates[sample_ids, , drop = FALSE]
    basal_model_formula <- colnames(covariates) |>
      make.names() |>
      paste(collapse = " + ")
  }

  # Select the winning model
  winning_models <- foreach::foreach(i = seq_len(nrow(selected_variables)),
                                     .combine = "rbind",
                                     .export = c("empty_lists", "sv_ids",
                                                 "sv_genot", "sv_env",
                                                 "sv_genot_rows", "methyl_col",
                                                 "env_row", "genot_col",
                                                 "covariates_i",
                                                 "basal_model_formula")) %dopar% { # For every VML
    #### Prepare data sets ####
    # Create the data frame with all the information for each VML.
    # Single-bracket indexing keeps each element wrapped in a length-1 list,
    # matching the structure of selected_variables[i, ]$selected_genot
    VML_index_i <- sv_ids[i]
    selected_genot_i <- sv_genot[i]
    selected_env_i <- sv_env[i]
    summ_vml_i <- summarized_methyl_VML[, methyl_col[i], drop = FALSE]
    colnames(summ_vml_i) <- "DNAme"
    if (!selected_env_i %in% empty_lists) {
      env_i <- environmental_matrix[env_row,
                                    unlist(selected_env_i),
                                    drop = FALSE]
    } else {
      env_i <- NULL
    }
    if (!selected_genot_i %in% empty_lists) {
      genot_i <- genotype_matrix[sv_genot_rows[[i]],
                                 genot_col,
                                 drop = FALSE] |>
        t()
      } else {
        genot_i <- NULL
      }
    full_data_vml_i <- cbind(summ_vml_i, env_i, genot_i, covariates_i)
    colnames(full_data_vml_i) <- make.names(colnames(full_data_vml_i))
    # Converted once and reused across every lm() call below, instead of
    # re-converting the same matrix on every single model fit
    full_data_vml_i_df <- as.data.frame(full_data_vml_i)
    # Whether this VML has G and E variables to fit
    genot_i_has_vars <- !selected_genot_i %in% empty_lists
    env_i_has_vars <- !selected_env_i %in% empty_lists
    snps_i <- if (genot_i_has_vars) unlist(selected_genot_i) else character(0)
    envs_i <- if (env_i_has_vars) unlist(selected_env_i) else character(0)

    #### Fit G Models ####
    ## Fit models involving G if G has selected variables
    if (genot_i_has_vars) {
      # One slot per model that will be fitted: for each SNP, the G model
      # followed by the GxE and G+E models of every environmental variable.
      # The slots are filled in place and bound in a single rbind below
      models_g_involving <- vector("list",
                                   length(snps_i) * (1 + 2 * length(envs_i)))
      pos <- 0L
      for (SNP in snps_i) { # For each SNP
        ### Fit G models
        model_g <- stats::lm(data = full_data_vml_i_df,
                             formula = paste("DNAme ~",
                                             make.names(SNP),
                                             "+",
                                             basal_model_formula))

        # Create data frame structure for the results
        model_g_df <- data.frame(model_group = "G")
        model_g_df$variables <- list(SNP)
        if (model_selection == "AIC") {
          model_g_df$AIC <- stats::AIC(model_g)
        } else if (model_selection == "BIC") {
          model_g_df$BIC <- stats::BIC(model_g)
        }
        model_g_df$tot_r_squared <- summary(model_g)$r.squared
        pos <- pos + 1L
        models_g_involving[[pos]] <- model_g_df
        #### Fit G+E and GxE models ####
        ### Fit GxE and G+E models if E is not empty
        for (env in envs_i) { # For every env var
          # Fit G + E
          model_ge <- stats::lm(data = full_data_vml_i_df,
                                formula = paste("DNAme ~",
                                                make.names(SNP),
                                                "+",
                                                make.names(env),
                                                "+",
                                                basal_model_formula)
                                )

          # Create data frame structure for the results
          model_ge_df <- data.frame(model_group = "G+E")
          model_ge_df$variables <- list(c(SNP, env))
          if (model_selection == "AIC") {
            model_ge_df$AIC <- stats::AIC(model_ge)
          } else if (model_selection == "BIC") {
            model_ge_df$BIC <- stats::BIC(model_ge)
          }
          model_ge_df$tot_r_squared <- summary(model_ge)$r.squared
          # Fit GxE
          model_gxe <- stats::lm(data = full_data_vml_i_df,
                                 formula = paste0("DNAme ~ ",
                                                 make.names(SNP),
                                                 " + ",
                                                 make.names(env),
                                                 " + ",
                                                 make.names(SNP),
                                                 "*", make.names(env),
                                                 " + ",
                                                 basal_model_formula)
                                 )
          # Create data frame structure for the results
          model_gxe_df <- data.frame(model_group = "GxE")
          model_gxe_df$variables <- list(c(SNP, env))
          if (model_selection == "AIC") {
            model_gxe_df$AIC <- stats::AIC(model_gxe)
          } else if (model_selection == "BIC") {
            model_gxe_df$BIC <- stats::BIC(model_gxe)
          }
          model_gxe_df$tot_r_squared <- summary(model_gxe)$r.squared
          # GxE is stored before G+E, keeping the order of the joint models
          pos <- pos + 1L
          models_g_involving[[pos]] <- model_gxe_df
          pos <- pos + 1L
          models_g_involving[[pos]] <- model_ge_df
        }
      }
      models_g_involving_df <- do.call(rbind, models_g_involving)
    } else {
      models_g_involving_df <- NULL
    }

    #### Fit E models ####
    # Only if E is not empty
    if (env_i_has_vars) { # For each env var
      models_e <- vector("list", length(envs_i))
      for (pos_e in seq_along(envs_i)) { # For every env var
        env <- envs_i[pos_e]
        # Fit E models
        model_e <- stats::lm(data = full_data_vml_i_df,
                             formula = paste("DNAme ~",
                                             make.names(env),
                                             "+",
                                             basal_model_formula)
                             )

        # Create data frame structure for the results
        model_e_df <- data.frame(model_group = "E")
        model_e_df$variables <- list(c(env))
        if (model_selection == "AIC") {
          model_e_df$AIC <- stats::AIC(model_e)
        } else if (model_selection == "BIC") {
          model_e_df$BIC <- stats::BIC(model_e)
        }
        model_e_df$tot_r_squared <- summary(model_e)$r.squared
        models_e[[pos_e]] <- model_e_df
      }
      models_e_df <- do.call(rbind, models_e)
    } else {
      models_e_df <- NULL
    }

    # Create object with the metrics for all the fitted models
    all_models_VML_i <- rbind(models_g_involving_df, models_e_df)

    # Select the best model per category (G,E,GxE,G+E) and compute its
    # delta AIC/BIC
    if (model_selection == "AIC") {
      best_models_VML_i <- all_models_VML_i |>
        dplyr::group_by(model_group) |>
        dplyr::filter(AIC == min(AIC)) |>
        # In case there are more than one model per group with the exact same
        # AIC, pick the first one
        dplyr::slice(1) |>
        dplyr::arrange(AIC, dplyr::desc(tot_r_squared)) |>
        dplyr::ungroup() |>
        dplyr::mutate(delta_aic = abs(AIC - dplyr::lead(AIC)))
    } else if (model_selection == "BIC") {
      best_models_VML_i <- all_models_VML_i |>
        dplyr::group_by(model_group) |>
        dplyr::filter(BIC == min(BIC)) |>
        # In case there are more than one model per group with the exact same
        # AIC, pick the first one
        dplyr::slice(1) |>
        dplyr::arrange(BIC, dplyr::desc(tot_r_squared)) |>
        dplyr::ungroup() |>
        dplyr::mutate(delta_bic = abs(BIC - dplyr::lead(BIC)))
    }


    # Create the final object that will be returned
    if (model_selection == "AIC") {
      winning_model_VML_i <- best_models_VML_i |>
        dplyr::filter(AIC == min(AIC)) |>
        # In case there is more than one model with the exact same AIC from
        # different groups, pick the one with the highest tot_r_squared
        dplyr::slice(1) |>
        dplyr::mutate(
          second_winner = best_models_VML_i$model_group[2],
          delta_r_squared = best_models_VML_i$tot_r_squared[1] - best_models_VML_i$tot_r_squared[2]
        )
    } else if (model_selection == "BIC") {
      winning_model_VML_i <- best_models_VML_i |>
        dplyr::filter(BIC == min(BIC)) |>
        # In case there is more than one model with the exact same AIC from
        # different groups, pick the one with the highest tot_r_squared
        dplyr::slice(1) |>
        dplyr::mutate(
          second_winner = best_models_VML_i$model_group[2],
          delta_r_squared = best_models_VML_i$tot_r_squared[1] - best_models_VML_i$tot_r_squared[2]
        )
    }

    # Test the winning model against the basal one and decompose variance for
    # the G, E and GxE components
    model_basal <- stats::lm(data = full_data_vml_i_df,
                             formula = paste("DNAme ~",
                                             basal_model_formula)
                             ) # set the basal model for comparing the rest
    if (model_selection == "AIC") {
      winning_model_VML_i$basal_AIC <- stats::AIC(model_basal)
    } else if (model_selection == "BIC") {
      winning_model_VML_i$basal_BIC <- stats::BIC(model_basal)
    }
    winning_model_VML_i$basal_rsquared <- summary(model_basal)$r.squared
    if (winning_model_VML_i$model_group == "G") {
      winning_lm <- stats::lm(data = full_data_vml_i_df,
                              formula = paste("DNAme ~",
                                              make.names(unlist(winning_model_VML_i$variables)),
                                              "+",
                                              basal_model_formula)
                              )
      if (is.null(covariates)) {
        # With no covariates, this model has a single predictor (the SNP), so
        # its unique contribution to R^2 is simply the model's own R^2.
        # relaimpo::calc.relimp.lm() requires 2+ regressors and errors out
        # otherwise.
        winning_model_VML_i$g_r_squared <- summary(winning_lm)$r.squared
      } else {
        r_decomp <- relaimpo::calc.relimp.lm(
          object = winning_lm,
          rela = FALSE,
          type = "last"
        ) # This would be the equivalent to using lmg and setting always = covariates.
        winning_model_VML_i$g_r_squared <- r_decomp$last[make.names(unlist(winning_model_VML_i$variables))[1]]
      }
      winning_model_VML_i$e_r_squared <- NA_real_
      winning_model_VML_i$gxe_r_squared <- NA_real_
    } else if (winning_model_VML_i$model_group == "E") {
      winning_lm <- stats::lm(data = full_data_vml_i_df,
                              formula = paste("DNAme ~",
                                              make.names(unlist(winning_model_VML_i$variables))[1],
                                              "+",
                                              basal_model_formula)
                              )
      if (is.null(covariates)) {
        # With no covariates, this model has a single predictor (the env.
        # variable), so its unique contribution to R^2 is simply the model's
        # own R^2. relaimpo::calc.relimp.lm() requires 2+ regressors and
        # errors out otherwise.
        winning_model_VML_i$e_r_squared <- summary(winning_lm)$r.squared
      } else {
        r_decomp <- relaimpo::calc.relimp.lm(
          object = winning_lm,
          rela = FALSE,
          type = "last"
        ) # This would be the equivalent to using lmg and setting always = covariates.
        winning_model_VML_i$e_r_squared <- r_decomp$last[make.names(unlist(winning_model_VML_i$variables))[1]]
      }
      winning_model_VML_i$g_r_squared <- NA_real_
      winning_model_VML_i$gxe_r_squared <- NA_real_
    } else if (winning_model_VML_i$model_group == "G+E") {
      winning_lm <- stats::lm(data = full_data_vml_i_df,
                              formula = paste("DNAme ~",
                                              make.names(unlist(winning_model_VML_i$variables))[1],
                                              "+",
                                              make.names(unlist(winning_model_VML_i$variables))[2],
                                              "+",
                                              basal_model_formula)
                              )
      r_decomp <- relaimpo::calc.relimp.lm(
        object = winning_lm,
        rela = FALSE,
        type = "lmg",
        # The columns of full_data_vml_i_df were passed through make.names()
        # above, so the covariate names must be too in order to refer to terms
        # that actually exist in the model
        always = make.names(colnames(covariates_i))
      )
      winning_model_VML_i$g_r_squared <- r_decomp$lmg[make.names(unlist(winning_model_VML_i$variables))[1]]
      winning_model_VML_i$e_r_squared <- r_decomp$lmg[make.names(unlist(winning_model_VML_i$variables))[2]]
      winning_model_VML_i$gxe_r_squared <- NA_real_
    } else if (winning_model_VML_i$model_group == "GxE") {
      winning_lm <- stats::lm(data = full_data_vml_i_df,
                              formula = paste0("DNAme ~ ",
                                              make.names(unlist(winning_model_VML_i$variables))[1],
                                              " + ",
                                              make.names(unlist(winning_model_VML_i$variables))[2],
                                              " + ",
                                              make.names(unlist(winning_model_VML_i$variables))[1],
                                              "*",
                                              make.names(unlist(winning_model_VML_i$variables))[2],
                                              " + ",
                                              basal_model_formula)
                              )
      r_decomp <- relaimpo::calc.relimp.lm(
        object = winning_lm,
        rela = FALSE,
        type = "lmg",
        # The columns of full_data_vml_i_df were passed through make.names()
        # above, so the covariate names must be too in order to refer to terms
        # that actually exist in the model
        always = make.names(colnames(covariates_i))
      )
      # This slightly underestimates the relative importance compared to not
      # using the covariates as the basal model, but in the interaction option
      # the computational time is greatly increased if the relative
      # contribution of all the other covariates is also estimated (which we
      # dont look at anyways). So, because of the high dimensional nature of
      # this package, this option will be used.
      winning_model_VML_i$g_r_squared <- r_decomp$lmg[make.names(unlist(winning_model_VML_i$variables))[1]]
      winning_model_VML_i$e_r_squared <- r_decomp$lmg[make.names(unlist(winning_model_VML_i$variables))[2]]
      winning_model_VML_i$gxe_r_squared <- r_decomp$lmg[paste0(
        make.names(unlist(winning_model_VML_i$variables))[1],
        ":",
        make.names(unlist(winning_model_VML_i$variables))[2])
        ]
    }

    winning_model_VML_i$VML_index <- VML_index_i
    # Return final object
    winning_model_VML_i
  }

  # Rearrange columns
  if (model_selection == "AIC") {
    winning_models <- winning_models |>
      dplyr::select(VML_index, model_group, variables, tot_r_squared,
                    g_r_squared, e_r_squared, gxe_r_squared, AIC,
                    second_winner, delta_aic, delta_r_squared, basal_AIC,
                    basal_rsquared)
  } else if (model_selection == "BIC") {
    winning_models <- winning_models |>
      dplyr::select(VML_index, model_group, variables, tot_r_squared,
                    g_r_squared, e_r_squared, gxe_r_squared, BIC, second_winner,
                    delta_bic, delta_r_squared, basal_BIC, basal_rsquared)
  }

  if (model_selection == "AIC") {
    return(winning_models |>
             # Attach VML with no variables selected in selectVariables()
             rbind(no_vars_VML |>
                     # remove empty columns
                     dplyr::select(-selected_genot, -selected_env) |>
                     dplyr::mutate(
                       model_group = "B",
                       variables = list(NA_character_),
                       tot_r_squared = NA_real_,
                       g_r_squared = NA_real_,
                       e_r_squared = NA_real_,
                       gxe_r_squared = NA_real_,
                       AIC = NA_real_,
                       second_winner = NA_character_,
                       delta_aic = NA_real_,
                       delta_r_squared = NA_real_,
                       basal_AIC = NA_real_,
                       basal_rsquared = NA_real_
                       )))
  }
  if (model_selection == "BIC") {
    return(winning_models |>
             # Attach VML with no variables selected in selectVariables()
             rbind(no_vars_VML |>
                     # remove empty columns
                     dplyr::select(-selected_genot, -selected_env) |>
                     dplyr::mutate(
                       model_group = "B",
                       variables = list(NA_character_),
                       tot_r_squared = NA_real_,
                       g_r_squared = NA_real_,
                       e_r_squared = NA_real_,
                       gxe_r_squared = NA_real_,
                       BIC = NA_real_,
                       second_winner = NA_character_,
                       delta_bic = NA_real_,
                       delta_r_squared = NA_real_,
                       basal_BIC = NA_real_,
                       basal_rsquared = NA_real_
                     )))
  }
}

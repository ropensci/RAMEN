test_that("medCorVMR throws errors when expected", {
  expect_error(
    RAMEN::medCorVMR(
      VML = "a",
      methylation_data = test_methylation_data
    ),
    "Please make sure the input VML belongs to the GRanges class."
  )
  VML_no_probes <- VML_test$VML
  S4Vectors::mcols(VML_no_probes)$probes <- NULL
  expect_error(
    RAMEN::medCorVMR(
      VML = VML_no_probes,
      methylation_data = test_methylation_data
    ),
    "Please make sure the VML object has the 'probes' column.",
    fixed = TRUE
  )
  expect_error(
    RAMEN::medCorVMR(
      VML = VML_test$VML,
      methylation_data = "a"
    ),
    "Please make sure the input methylation_data belongs to the data.frame class."
  )
})

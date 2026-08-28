argument_check <- function(object, data_type, extra_msg = NULL,
                           obj_name = deparse(substitute(object))) {
  correct_type <- switch(data_type,
                 list = is.list(object),
                 data.frame = is.data.frame(object),
                 character = is.character(object),
                 matrix = is.matrix(object),
                 numeric = is.numeric(object),
                 GRanges = methods::is(object, "GRanges")
                 )
  if (!correct_type) {
    msg <- paste0("Please make sure the input ", obj_name, " belongs to the ",
                  data_type, " class.")
    if (!is.null(extra_msg)) msg <- paste0(msg, " ", extra_msg)
    stop(msg)
  }

  is_empty <- switch(data_type,
                 list = length(object) == 0,
                 data.frame = nrow(object) == 0 || ncol(object) == 0,
                 character = length(object) == 0,
                 matrix = nrow(object) == 0 || ncol(object) == 0,
                 numeric = length(object) == 0,
                 GRanges = length(object) == 0
                 )
  if (is_empty) {
    msg <- paste0("Please make sure the input ", obj_name, " is not empty.")
    if (!is.null(extra_msg)) msg <- paste0(msg, " ", extra_msg)
    stop(msg)
  }
}

argument_char_options <- function(object, options, extra_msg = NULL) {
  obj_name <- deparse(substitute(object))
  argument_check(object, "character", obj_name = obj_name)
  if (!length(object) == 1) {
    stop("Please make sure the input ", obj_name,
         " is a character object of length 1")
  }
  if (!object %in% options) {
    msg <- paste0("Please make sure the input ", obj_name,
                  " is one of the following options: ",
                  paste(options, collapse = ", "), " .")
    if (!is.null(extra_msg)) msg <- paste0(msg, " ", extra_msg)
    stop(msg)
  }
}

columns_exist <- function(data.frame, columns, extra_msg = NULL) {
  if (!all(columns %in% colnames(data.frame))) {
    msg <- paste0("The object ", deparse(substitute(data.frame)),
                  " does not have the required columns: ",
                  paste(columns, collapse = ", "), " .")
    if (!is.null(extra_msg)) msg <- paste0(msg, " ", extra_msg)
    stop(msg)
  }
}

vectors_match <- function(object_1, object_2,
                          object_1_name = "object_1",
                          object_2_name = "object_2",
                          extra_msg = NULL) {
  if (!setequal(object_1, object_2)) {
    msg <- paste0("The objects ", object_1_name, " and ", object_2_name,
                  " must match.")
    if (!is.null(extra_msg)) msg <- paste0(msg, " ", extra_msg)
    stop(msg)
  }
}

finite_numeric_check <- function(object, extra_msg = NULL) {
  if (
    !is.numeric(object) ||
    anyNA(object) ||
    # is.infinite() builds a logical vector as large as the object it is given,
    # which is a substantial allocation for a genotype matrix. Integers cannot
    # hold Inf, and anyNA() above already covers NA_integer_, so the check is
    # only needed for doubles. Written as part of the same || chain so that it
    # is still only reached when the cheaper checks pass.
    (!is.integer(object) && any(is.infinite(object)))
  ) {
    msg <- paste0("Please make sure the object ", deparse(substitute(object)),
                  " contains only finite numeric values ",
                  "(i.e., no NA, NaN or Inf).")
    if (!is.null(extra_msg)) msg <- paste0(msg, " ", extra_msg)
    stop(msg)
  }
}

empty_lists <- c(list(NULL), list(""), list(NA), list(character(0)))

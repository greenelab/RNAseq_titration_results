check_options <- function(opt) {
  # this function checks standardized command line options given to scripts
  # options ending with "_input" are checked to see if the input file exists
  # options ending with "_output" are checked to see if the output directory exists
  #   and if the output file already exists, will it be overwritten or not
  # all messages and errors are reported
  # if there are any errors, the script stops

  my_errors <- list()
  my_messages <- list()

  for(option in names(opt)){

    if (is.na(opt[[option]])) { # all required options should default to NA_character_ or NA_integer_
      my_errors[[option]] <- stringr::str_c("\nOption given for --", option,
                                            " is missing and must be specified.")
    } else if (option == "cancer_type") {
      if (!(opt[[option]] %in% c("BRCA", "GBM"))) { # cancer type must be BRCA or GBM
        my_errors[[option]] <- stringr::str_c("\nCancer type given for --", option,
                                             " (", opt[[option]], ") ",
                                             " must be BRCA or GBM.")
      }
    } else if (option == "predictor") {
      if (!(opt[[option]] %in% c("subtype", "TP53", "PIK3CA"))) { # predictor must be subtype or TP53 or PIK3CA
        my_errors[[option]] <- stringr::str_c("\nPredictor given for --", option,
                                              " (", opt[[option]], ") ",
                                              " must be subtype, TP53, or PIK3CA.")
      }
    }  else if (option == "subtype_vs_subtype") {
      two_subtypes <- as.vector(stringr::str_split(opt[[option]], pattern = ",", simplify = TRUE))
      if (length(two_subtypes) != 2) {
        my_errors[[option]] <- stringr::str_c("\nSubtypes given for --", option,
                                             " (", opt[[option]], ") ",
                                             " must have (only) two comma-separated subtypes.")
      }
      
    } else if (stringr::str_ends(option, "_input")) { # option related to inputs
      if (!file.exists(opt[[option]])) {
        my_errors[[option]] <- stringr::str_c("\nInput file given for --", option,
                                              " (", opt[[option]], ") ",
                                              "does not exist.")
      }
    } else if (stringr::str_ends(option, "_output")) { # option related to outputs
      if (file.exists(opt[[option]])) { # if output file already exists
        if (opt$overwrite) { # overwrite is TRUE if given
          my_messages[[option]] <- stringr::str_c("\nOutput file given for --", option,
                                                  " (", opt[[option]], ") ",
                                                  "already exists and will be overwritten (--overwrite is set).")
        } else { # overwrite defaults to FALSE unless given
          my_errors[[option]] <- stringr::str_c("\nOutput file given for --", option,
                                                " (", opt[[option]], ") ",
                                                "already exists and will not be overwritten (use --overwrite).")
        }
      } else if (!dir.exists(dirname(opt[[option]]))) { # if output directory does not exist
        my_errors[[option]] <- stringr::str_c("\nOutput directory given for --", option,
                                              " (", dirname(opt[[option]]), ") ",
                                              "does not exist.")
      }
    } else if (option == "ncores") {
      if (!is.integer(opt[[option]]) | opt[[option]] < 1) {
        my_errors[[option]] <- stringr::str_c("\nNumber of cores given for --", option,
                                              " must be a positive integer.")
      }
    }
  }

  if (length(my_messages) > 0) {
    message("  Messages:", my_messages, "\n")
  }
  if (length(my_errors) > 0) {
    message("  Errors:", my_errors, "\n")
    stop()
  }
}

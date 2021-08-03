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

    if (is.null(opt[[option]])) { # all required options should default to NULL
      my_errors[[option]] <- stringr::str_c("\nOption given for --", option,
                                            " is missing and must be specified.")
    } else if (option == "cancer_type") {
      
      if (!(opt[[option]] %in% c("BRCA", "GBM"))) { # cancer type must be BRCA or GBM
        my_errors[[option]] <- strinr::str_c("\nCancer type given for --", option,
                                             " (", opt[[option]], ") ",
                                             " must be BRCA or GBM.")  
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

# S. Foltz Oct 2021
# The purpose of this script is to visualize normalized gene expression
# and compare values from matched microarray and RNA-seq samples
# USAGE: Rscript visualize_expression.R --cancer_type --predictor --null_model --seed

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NA_character_,
                        help = "Predictor used"),
  optparse::make_option("--null_model",
                        action = "store_true",
                        default = FALSE,
                        help = "Refer to models with permuted dependent variable (within subtype if predictor is a gene)"),
  optparse::make_option("--seed",
                        default = 1234,
                        help = "Set a seed to ensure reproducible results when subsampling genes")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(library(tidyverse))
source(here::here("util", "color_blind_friendly_palette.R"))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
null_model <- opt$null_model
file_identifier <- ifelse(null_model,
                          str_c(cancer_type, predictor, "null", sep = "_"),
                          str_c(cancer_type, predictor, sep = "_"))

# set seed
set.seed(opt$seed)

# define directories
plot.dir <- here::here("plots")
norm.dir <- here::here("normalized_data")
viz.dir <- file.path(plot.dir, "visualize_expression")

# define input files
normalized_test_data_filename <- list.files(norm.dir,
                                            pattern = str_c(file_identifier,
                                                            "_array_seq_test_data_normalized_list_"),
                                            full.names = TRUE)[1]

normalized_test_data <- read_rds(normalized_test_data_filename)

### functions ------------------------------------------------------------------

plot_matched_expression <- function(array_values, seq_values,
                                    method_title, plot_type,
                                    output_directory, filename_lead) {
  
  # This function creates a plot for expression values from matched array and RNA-seq samples
  # The function can produce a plot with points (alpha = 0.1) or a hex grid to show density
  # Inputs:
  #   array_values = a vector of array values
  #   seq_values = vector of seq values, matched to array values
  #   method_title = something informative that will define the plot title and output filename
  #   plot_type = either 'point' or 'hex' depending on the desired plot type
  #   output_directory = output directory of PDF
  #   filename_lead = start of the output filename
  # Outputs:
  #   a PDF of the plot is saved to output_directory
  
  this_plot <- ggplot(mapping = aes(x = array_values,
                                      y = seq_values))
    
  if (plot_type == "point") {
    this_plot <- this_plot +
      geom_point(alpha = 0.1,
                 shape = 16)
  } else if (plot_type == "hex") {
    this_plot <- this_plot +
      geom_hex()
  } else {
    stop("Plot type must be 'point' or 'hex'.")
  }
  
  this_plot <- this_plot + 
    geom_abline(lty = 2, # dashed red x-y line 
                color = "red") +
    geom_smooth(method = "gam", # fit a curve to the data
                formula = y ~ s(x, bs = "cs")) + # loess no good for large n
    labs(x = "Microarray expression values",
         y = "RNA-seq expression values",
         title = method_title) +
    theme_minimal()
  
  if (method_title != "UN") {
    this_plot <- this_plot +
      coord_fixed()
  }
  
  ggsave(plot = this_plot,
         filename = file.path(output_directory,
                              str_c(filename_lead,
                                    method_title,
                                    plot_type,
                                    "pdf",
                                    sep = ".")),
         height = 7.25,
         width = 7.25)
}

#### plot matched comparison of matched microarray and RNA-seq -----------------

gene_rows_included <- sort(sample(1:nrow(normalized_test_data$array$log),
                                  size = 1000, # select 1000 random genes
                                  replace = FALSE))

norm_methods <- names(normalized_test_data$seq) # get all normalization methods

for (nm in norm_methods) {
  
  if (nm %in% c("seurat")) next
  
  if (nm == "tdm") {
    # array has no TDM (it is already log)
    array_values <- as.vector(as.matrix(normalized_test_data$array[["log"]][gene_rows_included, -1]))
    for (pct_rna_seq in as.character(seq(0, 90, 10))) { # NULL at 100% RNA-seq
      # only seq varies across %RNA-seq
      seq_values <- as.vector(as.matrix(normalized_test_data$seq[[nm]][[pct_rna_seq]][gene_rows_included, -1]))
      method_title <- str_c(str_to_upper(nm), pct_rna_seq, sep = "_")
      
      plot_matched_expression(array_values, seq_values,
                              method_title, plot_type = "hex",
                              viz.dir, file_identifier)
      
    } 
  } else if (nm %in% c("qn", "qn-z")) {
    array_values <- as.vector(as.matrix(normalized_test_data$array[[nm]][gene_rows_included, -1]))
    for (pct_rna_seq in as.character(seq(0, 100, 10))) {
      # only seq varies across %RNA-seq
      seq_values <- as.vector(as.matrix(normalized_test_data$seq[[nm]][[pct_rna_seq]][gene_rows_included, -1]))
      method_title <- str_c(str_to_upper(nm), pct_rna_seq, sep = "_")
      
      plot_matched_expression(array_values, seq_values,
                              method_title, plot_type = "hex",
                              viz.dir, file_identifier)
      
    }
  } else { # test data for normalization methods that do not vary with RNA-seq % in training data
    array_values <- as.vector(as.matrix(normalized_test_data$array[[nm]][gene_rows_included, -1]))
    seq_values <- as.vector(as.matrix(normalized_test_data$seq[[nm]][gene_rows_included, -1]))
    method_title <- str_to_upper(nm)
    
    plot_matched_expression(array_values, seq_values,
                            method_title, plot_type = "hex",
                            viz.dir, file_identifier)
    
  }
}

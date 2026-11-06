#' Encode DNA sequences as nucleotide-position factors
#'
#' Converts equal-length DNA strings into a data frame, with one column
#' for each nucleotide position. Each column is a categorical factor.
#'
#' @param dna_strings A character vector of equal-length DNA strings
#'   containing the nucleotides A, T, C and G.
#'
#' @return A data frame with one row per DNA string and one factor column
#'   per nucleotide position. Columns are named nt_pos1, nt_pos2, and so on,
#'   and each factor has the levels A, T, C and G.
#'
#' @examples
#' # A runnable example will be added in Step 5.
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A probabilities and statuses for multiple RNA sites
#'
#' Uses a fitted random forest model to predict m6A probabilities for
#' multiple RNA sites and classifies each site using a probability threshold.
#'
#' @param ml_fit A fitted randomForest model for m6A prediction.
#' @param feature_df A data frame containing the columns gc_content,
#'   RNA_type, RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation and DNA_5mer. The length and distance
#'   variables must already be log2-transformed.
#' @param positive_threshold A numeric probability threshold, defaulting
#'   to 0.5. Probabilities strictly greater than this value are Positive;
#'   probabilities equal to or below it are Negative.
#'
#' @return The input data frame with two additional columns:
#'   predicted_m6A_prob, a numeric predicted probability, and
#'   predicted_m6A_status, a character string of either Positive or Negative.
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, example_df[1:3, ], positive_threshold = 0.6)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict m6A probability and status for a single RNA site
#'
#' Uses a fitted random forest model to predict the m6A probability of
#' one RNA site and classifies it using a probability threshold.
#'
#' @param ml_fit A fitted randomForest model for m6A prediction.
#' @param gc_content Numeric GC content of the RNA site.
#' @param RNA_type A character string: mRNA, lincRNA, lncRNA or pseudogene.
#' @param RNA_region A character string: CDS, intron, 3'UTR or 5'UTR.
#' @param exon_length Numeric log2-transformed exon length.
#' @param distance_to_junction Numeric log2-transformed distance to junction.
#' @param evolutionary_conservation Numeric evolutionary conservation score.
#' @param DNA_5mer A character string containing the five-nucleotide DNA sequence.
#' @param positive_threshold A numeric probability threshold, defaulting
#'   to 0.5. Probabilities strictly greater than this value are Positive;
#'   probabilities equal to or below it are Negative.
#'
#' @return A named character vector of length two. The element
#'   predicted_m6A_prob contains the predicted probability converted to
#'   text, and predicted_m6A_status contains either Positive or Negative.
#'
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8, DNA_5mer = "GGACA",
#'                   positive_threshold = 0.5)
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}

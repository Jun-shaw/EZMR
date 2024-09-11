#' @name result_MR
#' @title Organize and export the MR results.
#' @description Organize and export the results of the Mendelian Randomization (MR) analysis in a clear and structured manner. This involves collecting the key findings, such as effect estimates, confidence intervals, and p-values, and arranging them in a well-labeled table for easy interpretation.
#'
#' @param folder.path Directory for storing MR results files.
#' @param exposure.cat The category of the exposure factor.
#' @param outcome.name The name of the outcome factor.
#'
#' @author Junxiao Shen
#' @examples
#' result_MR (folder.path='MR_result/',
#'            exposure.cat='Biomarkers'
#'            outcome.name='PCA')
#'
#' @export
#'
result_MR <- function (folder.path='',
                       exposure.cat='',
                       outcome.name='')
{
  # Get the names of all .csv files in the specified folder
  file_list <- list.files(path = folder.path, pattern = "\\.csv$", full.names = TRUE)

  # Display the list of file names
  print(file_list)

  # Create an empty data frame to store results
  result_df <- data.frame(matrix(ncol = 7, nrow = 0))

  print('\n=> Looping through files and building data rows--------')
  for (file in file_list) {
    # Read the data from the file
    data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)

    # Extract the required data
    biomarkers <- data[data$method == "Inverse variance weighted", "id.exposure"]
    IVW_OR <- data[data$method == "Inverse variance weighted", "or"]
    MR_Egger_OR <- data[data$method == "MR Egger", "or"]
    Weighted_median_OR <- data[data$method == "Penalised weighted median", "or"]
    IVW_P <- data[data$method == "Inverse variance weighted", "pval"]
    MR_Egger_P <- data[data$method == "MR Egger", "pval"]
    Weighted_median_P <- data[data$method == "Penalised weighted median", "pval"]

    # Create a data row for the current file
    row_data <- c(biomarkers, IVW_OR, MR_Egger_OR, Weighted_median_OR, IVW_P, MR_Egger_P, Weighted_median_P)

    # Add the data row to the result data frame
    result_df <- rbind(result_df, row_data)
  }

  # Set the column names for the result data frame
  colnames(result_df) <- c(exposure.cat, "IVW_OR", "MR_Egger_OR", "Weighted_median_OR", "IVW_P", 'MR_Egger_P', 'Weighted_median_P')

  print('\n=> Printing the result data frame')
  # Write the result data frame to a CSV file
  write.csv(result_df, paste(outcome.name, "_all_result.csv"), row.names = FALSE)
}

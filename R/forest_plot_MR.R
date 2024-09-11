#' @name forest_plot_MR
#' @title Create a forest plot for the MR results.
#' @description Construct a forest plot to depict the results of the Mendelian Randomization (MR) analysis. This graphical representation should encapsulate the effect estimates along with their corresponding confidence intervals for each of the key variables. By elegantly illustrating both the magnitude and precision of the effects, the forest plot will provide a clear visual interpretation of the findings, allowing for an immediate grasp of the relationships under investigation.
#'
#' @import grid
#' @import forestploter
#' @import forestplot
#'
#' @param file The path to the input file.
#' @param column_num The total number of columns in the input file.
#' @param shape_col The column that defines the shape of the points in a plot
#' @param line_col The column used to define the color or type of lines in a plot
#' @param text_col The column that contains text labels for the points in the plot
#' @param P_threshold The p-value threshold from Mendelian Randomization (MR) results, indicating whether there is a causal relationship between the exposure and outcome factors. A lower p-value suggests stronger evidence for causality.
#' @param OR_range A range for odds ratios (OR) to filter results based on their effect size, which is a vector indicating the lower and upper bounds.
#' @param cutoff A cutoff value for the odds ratio (OR)
#'
#' @author Junxiao Shen
#' @examples
#' forest_plot_MR (file='data/pca_result.csv',
#'                 column.num=10,
#'                 shape.col='#4575b4',
#'                 line.col='#CE5C69',
#'                 text.ol='#4575b4',
#'                 P.threshold=0.05,
#'                 OR.range = c(0,2),
#'                 cutoff = c(0,1,2))
#'
#' @export
#'
forest_plot_MR <- function (file='',
                          column.num=10,
                          shape.col='#4575b4',
                          line.col='#CE5C69',
                          text.ol='#4575b4',
                          P.threshold=0.05,
                          OR.range = c(0,2),
                          cutoff = c(0,1,2))
{
  # Reading data
dt <- read.csv(file, header = TRUE, fileEncoding = 'GBK')
dt <- dt[, 1:column.num]  # (3, 5)
dt$` ` <- paste(rep(" ", 20), collapse = " ")
dt$'OR(95%CI)'<-ifelse(is.na(dt$or),"",
                       sprintf('%.2f(%.2f to %.2f)',
                               dt$or,dt$or_lci95,dt$or_uci95))
dt[is.na(dt)] <- " "

# Plotting
tm <- forest_theme(base_size = 10, # Base size
                   # Shape, line type, color, and width of the confidence interval points
                   ci_pch = 20,
                   ci_col = shape.col, # #762a83
                   ci_lty = 1,
                   ci_lwd = 2.3,
                   ci_Theight = 0.2, # Short vertical lines at both ends of the confidence interval

                   # Reference line width, shape, and color
                   refline_lwd = 1.5,
                   refline_lty = "dashed",
                   refline_col = line.col,

                   # Fill color and border color for the summary diamond
                   summary_fill = shape.ol,
                   summary_col = shape.col,

                   # Footnote size, font style, and color
                   footnote_cex = 1.1,
                   footnote_fontface = "italic",
                   footnote_col = text.col)

p <- forest(dt[, c(1:4, 11:12, 8:10)], # Selecting data from columns 1 to 3 and 8 to 9 as elements of the forest plot
            est = dt$or, # HR
            lower = dt$or_lci95, # Lower limit of the confidence interval
            upper = dt$or_uci95, # Upper limit of the confidence interval
            sizes = 0.6, # Size of point estimate boxes, mapped to the standard error
            ci_column = 5, # Column for displaying the confidence interval
            ref_line = 1,
            xlim = OR.range,
            ticks_at = cutoff,
            arrow_lab = c('Protective factor', 'Risk factor'),
            footnote = paste0('P<', P.threshold, ' was considered statistically significant'),
            theme = tm) # Theme
print(p)
}

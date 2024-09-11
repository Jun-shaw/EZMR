#' @name visualize_MR
#' @title The heatmap illustrates the results of the MR analysis.
#' @description The heatmap visually represents the outcomes of the Mendelian Randomization (MR) analysis, providing a comprehensive overview of the relationships and effects observed in the study. Through this colorful graphical display, one can easily discern patterns, correlations, and variations among the different variables examined, enhancing the clarity and interpretability of the results.
#'
#' @import ComplexHeatmap
#' @import data.table
#' @import tidyverse
#' @import circlize
#' @import dplyr
#' @import RColorBrewer
#'
#' @param file Filename of the results after MR organization.
#' @param outcome.num The number of the outcome factor.
#' @param outcome.name The name of the outcome factor.
#' @param color1 The gradient colors of the first tier of the heat map.
#' @param color2 The gradient colors of the second tier of the heat map.
#' @param color3 The gradient colors of the third tier of the heat map.
#' @param color4 The gradient colors of the forth tier of the heat map.
#' @param range Color range of the heatmap.
#' @param threshold1 The primary criterion for significance values.
#' @param threshold2 The secondary criterion for significance values.
#' @param category Logical variables that ascertain the presence of diverse categories of exposure factors.
#' @param category.colors Colors of various classifications.
#'
#' @author Junxiao Shen
#' @examples
#' visualize_MR (file='data/all_result.csv',
#'              outcome.num=3,
#'              outcome.name=c('A','B','C'),
#'              color1=c("#003399","#ffffff","#cc0033"),
#'              color2=c( "#663366", "#ffffff","#ff9933"),
#'              color3=c( "#cccc00", "#ffffff", "#66cc99"),
#'              color4=NA,
#'              range=c(0.5, 1, 2),
#'              threshold1=0.05,
#'              threshold2=0.01,
#'              category=T,
#'              category.colors=c('#CC88B0', '#DBE0ED', '#87B5B2', '#F4CEB4', '#F1DFA4', '#998DB7'))
#'
#' @export
#'
visualize_MR <-function(file,
                        outcome.num=3,
                        outcome.name=c('A','B','C'),
                        color1=c("#003399","#ffffff","#cc0033"),
                        color2=c( "#663366", "#ffffff","#ff9933"),
                        color3=c( "#cccc00", "#ffffff", "#66cc99"),
                        color4=NA,
                        range=c(0.5, 1, 2),
                        threshold1=0.05,
                        threshold2=0.01,
                        category=T,
                        category.colors=c('#CC88B0', '#DBE0ED', '#87B5B2', '#F4CEB4', '#F1DFA4', '#998DB7')
                        )
{
  print('\n=> Processing data...')
  data = fread(file)
  df = data %>% as.data.frame()  # Transform data into a data frame
  rownames(df) = df[,1]
  df = df[, -1]
  for (i in 1:outcome.num) {
    cir <- df[((i - 1) * 6 + 1):((i - 1) * 6 + 3)]
    star <- df[((i - 1) * 6 + 4):(i * 6)]
    assign(paste("cir", i, sep = ""), cir)
    assign(paste("star", i, sep = ""), star)
  }
  cir_group <- as.matrix(df[, ncol(df)])

  get_stars <- function(value) {
    if (value < threshold2) {
      stars <- "**"
    } else if (value < threshold1) {
      stars <- "*"
    } else {
      stars <- ""
    }
    return(stars)
  }
  # If the matrix data is grouped, the split parameter can specify the categorical variable
  if (category == TRUE) {
    ann_row = data.frame(pathway = df$category)  # Annotate rows for subsequent heatmap splitting
    row.names(ann_row) = rownames(cir1)
    ann_row <- as.matrix(ann_row)  # Required as matrix in the circlize function
  }

  print('\n=> Starting to generate heatmap.......')

  circos.clear()
  # Define the gradient colors for the heatmap:
  mycol1 = colorRamp2(range, color1)  # Set the legend color
  if (outcome.num > 1) { mycol2 = colorRamp2(range, color2) }
  if (outcome.num > 2) { mycol3 = colorRamp2(range, color3) }
  if (outcome.num > 3) { mycol4 = colorRamp2(range, color4) }

  mycol <- list(mycol1)
  if (outcome.num > 1) { mycol <- list(mycol1, mycol2) }
  if (outcome.num > 2) { mycol <- list(mycol1, mycol2, mycol3) }
  if (outcome.num > 3) { mycol <- list(mycol1, mycol2, mycol3, mycol4) }

  circos.par(gap.after = c(2, 2, 2, 2, 2, 90))  # Adjusts the spacing between the ends of the rings; larger numbers mean wider gaps

  circos.heatmap(cir1, col = mycol1,

                 split = ann_row,  # Split heatmap using row annotation

                 rownames.col = "black",

                 show.sector.labels = FALSE,

                 track.height = 0.15,  # Track height; larger values yield thicker rings

                 rownames.side = "outside",

                 rownames.cex = 0.8,  # Font size

                 rownames.font = 1,  # Font weight

                 bg.border = "black",  # Background border color

                 dend.side = "inside",  # Control the direction of the row clustering tree; 'inside' for showing inside the ring

                 cluster = FALSE,  # 'TRUE' for row clustering, 'FALSE' for no clustering display

                 dend.track.height = 0.2,  # Adjust the height of the row clustering tree

                 dend.callback = function(dend, m, si) {  # Callback function for clustering tree reordering or color addition
                   color_branches(dend, k = 10, col = 1:10)  # Change clustering tree color to black
                 }
  )

  if (outcome.num > 1) {
    circos.heatmap(cir2,

                   col = mycol2,

                   split = ann_row, track.height = 0.15,

                   bg.border = "black",  # Background border color

                   rownames.cex = 0.3)  # Add second heatmap
  }
  if (outcome.num > 2) {
    circos.heatmap(cir3,

                   col = mycol3,

                   split = ann_row, track.height = 0.15,

                   bg.border = "black",  # Background border color

                   rownames.cex = 0.3)
  }
  if (outcome.num > 3) {
    circos.heatmap(cir4,

                   col = mycol4,

                   split = ann_row, track.height = 0.15,

                   bg.border = "black",  # Background border color

                   rownames.cex = 0.3)
  }
  if (category == TRUE) {
    # Create a pairing of group classification information and colors
    group_names <- unique(df$category)
    group_colors <- category.colors

    # Create a named vector to pair these
    row_colors <- setNames(group_colors, group_names)
    circos.heatmap(cir_group,

                   col = row_colors,

                   split = ann_row,

                   track.height = 0.05,

                   bg.border = "black",  # Background border color

                   rownames.cex = 0.3,

                   show.sector.labels = FALSE)
  }

  # Add column names

  # First circular column names

  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 6) {  # The last sector
      cn = colnames(cir1)
      n = length(cn)

      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),  # x coordinates
                  (1:n) * 1.2 + 5,  # Adjust y coordinates: row spacing + distance from the center
                  cn, cex = 1, adj = c(0, 0), facing = "inside")
    }
  }, bg.border = NA)

  print('\n=> Adding legend......')
  library(circlize)
  library(gridBase)

  for (i in 1:outcome.num) {
    lg_Exp = Legend(title = outcome[i], col_fun = mycol[[i]], direction = c("vertical"))
    assign(paste("lg_Exp", i, sep = ""), lg_Exp)
  }

  # Create legend

  circle_size = unit(0.6, "snpc")

  h = dev.size()

  lgd_list = packLegend(lg_Exp1, lg_Exp2, lg_Exp3, max_height = unit(2 * h, "inch"))
  if (outcome.num > 1) { lgd_list = packLegend(lg_Exp1, lg_Exp2, max_height = unit(2 * h, "inch")) }
  if (outcome.num > 2) { lgd_list = packLegend(lg_Exp1, lg_Exp2, lg_Exp3, max_height = unit(2 * h, "inch")) }
  if (outcome.num > 3) { lgd_list = packLegend(lg_Exp1, lg_Exp2, lg_Exp3, lg_Exp4, max_height = unit(2 * h, "inch")) }
  draw(lgd_list, x = circle_size, just = "midle")

}


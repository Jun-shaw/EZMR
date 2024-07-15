#' Calculating the area under the curve after developing the category predictive model
#'
#' @param res.by.ML.Dev.Pred.Category.Sig  Output of function ML.Dev.Pred.Category.Sig
#' @param cohort.for.cal A data frame with the 'ID' and 'Var' as the first two columns. Starting in the fourth column are the variables that contain variables of the model you want to build. The second column 'Var' only contains 'Y' or 'N'.
#'
#' @return A data frame containing the AUC of each predictive model.
#' @export
#'
#' @examples
#'
forest_plot_MR <- function (file='',
                          column_num=10,
                          shape_col='#4575b4',
                          line_col='#CE5C69',
                          text_col='#4575b4',
                          P_threshold=0.05,
                          OR_range = c(0,2),
                          cutoff = c(0,1,2)
                          )
{
#布置环境
library(grid)
library(forestploter)
library(forestplot)
#读取数据
dt <- read.csv("data/前列腺癌result.csv",header = T,fileEncoding = 'GBK')
dt<-dt[,1:column_num]#(3,5)
dt$` ` <- paste(rep(" ", 20), collapse = " ")
dt$'OR(95%CI)'<-ifelse(is.na(dt$or),"",
                       sprintf('%.2f(%.2f to %.2f)',
                               dt$or,dt$or_lci95,dt$or_uci95))
dt[is.na(dt)] <- " "
#绘图
tm <- forest_theme(base_size = 10, # 基础大小
                   # 可信区间点的形状，线型、颜色、宽度
                   ci_pch = 20,
                   ci_col = shape_col, # #762a83
                   ci_lty = 1,
                   ci_lwd = 2.3,
                   ci_Theight = 0.2, # 可信区间两端加短竖线

                   # 参考线宽度、形状、颜色
                   refline_lwd = 1.5,
                   refline_lty = "dashed",
                   refline_col = line_col,

                   # 汇总菱形的填充色和边框色
                   summary_fill = shape_col,
                   summary_col = shape_col,

                   # 脚注大小、字体、颜色
                   footnote_cex = 1.1,
                   footnote_fontface = "italic",
                   footnote_col = text_col)
p <- forest(dt[,c(1:4, 11:12,8:10)],#选择数据1-3列和8-9列作为森林图元素
            est = dt$or,#HR
            lower = dt$or_lci95,#可信区间下限
            upper = dt$or_uci95,#可信区间上限
            sizes = 0.6,#点估计框大小，用标准误映射
            ci_column = 5,#可信区间在第几列展示
            ref_line = 1,
            xlim = OR_range,
            ticks_at = cutoff,
            arrow_lab = c('protective factor','risk factor'),
            footnote = paste0('P<',P_threshold, 'was considered statistically significant'),
            theme = tm)#主题
print(p)
}

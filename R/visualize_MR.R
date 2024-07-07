
#多轨热图绘制#
library(ComplexHeatmap)
library(data.table)
library(tidyverse)
library(circlize)
library(dplyr)
library(RColorBrewer)

visualize_MR <-function(filename,
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
                        group_colors=c('#CC88B0', '#DBE0ED', '#87B5B2', '#F4CEB4', '#F1DFA4', '#998DB7')
                        )
{
#假设有个热图的矩阵数据（这里仅为一组重复两次以作示范）

print('处理数据......')
data=fread(filename)
df=data %>% as.data.frame()#数据转化为数据框
rownames(df)=df[,1]
df=df[,-1]
for (i in 1:outcome.num) {
  cir <- df[((i-1)*6+1):((i-1)*6+3)]
  star <- df[((i-1)*6+4):(i*6)]
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
#如果矩阵数据分组，可用split参数来指定分类变量
if(category==T){
ann_row =data.frame(pathway=df$category)#对行进行注释，用于后续的热图分裂

row.names(ann_row) = rownames(cir1)

ann_row <- as.matrix(ann_row)}#在circlize函数中，需要为matrix

print('开始生成热图.......')

circos.clear()
#定义热图颜色梯度：
mycol1=colorRamp2(range,color1)#设置legend颜色，
if (outcome.num>1){mycol2 = colorRamp2(range,color2)}
if (outcome.num>2){mycol3 = colorRamp2(range,color3)}
if (outcome.num>3){mycol4 = colorRamp2(range,color4)}

mycol <- list(mycol1)
if (outcome.num>1){mycol <- list(mycol1,mycol2)}
if (outcome.num>2){mycol <- list(mycol1,mycol2,mycol3)}
if (outcome.num>3){mycol <- list(mycol1,mycol2,mycol3,mycol4)}

circos.par(gap.after=c(2,2,2,2,2,90))#circos.par()调整圆环首尾间的距离，数值越大，距离越宽#让分裂的一个口大一点，可以添加行信息

circos.heatmap(cir1,col=mycol1,

               split=ann_row, #用行注释分裂热图

               rownames.col="black",

               show.sector.labels = F,

               track.height = 0.15, #轨道的高度，数值越大圆环越粗

               rownames.side="outside",

               rownames.cex=0.8,#字体大小

               rownames.font=1,#字体粗细

               bg.border="black", #背景边缘颜色

               dend.side="inside",#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈

               cluster=F,#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类

               dend.track.height=0.2,#调整行聚类树的高度

               dend.callback=function(dend,m,si) { #dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称

                 color_branches(dend,k=10,col=1:10) #color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色

               }

)
if (outcome.num>1){
circos.heatmap(cir2,

               col = mycol2,

               split=ann_row,track.height = 0.15,

               bg.border="black", #背景边缘颜色

               rownames.cex=0.3)#加入第二个热图
}
if (outcome.num>2){
circos.heatmap(cir3,

               col = mycol3,

               split=ann_row,track.height = 0.15,

               bg.border="black", #背景边缘颜色

               rownames.cex=0.3,
)
}
if (outcome.num>3){
  circos.heatmap(cir4,

                 col = mycol4,

                 split=ann_row,track.height = 0.15,

                 bg.border="black", #背景边缘颜色

                 rownames.cex=0.3,
  )
}
if(category==T){
# 创建行名分类信息与颜色的配对
group_names <-  unique(df$category)
group_colors <- group_colors

# 创建命名向量将这两者配对
row_colors <- setNames(group_colors, group_names)
circos.heatmap(cir_group,

               col = row_colors,

               split=ann_row,

               track.height = 0.05,

               bg.border="black", #背景边缘颜色

               rownames.cex=0.3,

               show.sector.labels = F,)
}

#添加列名#

#第一个环形列名

circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){

  if(CELL_META$sector.numeric.index==6){# the last sector

    cn=colnames(cir1)

    n=length(cn)

    circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(1,"mm"),#x坐标

                (1:n)*1.2+5,#调整y坐标,行距+距离中心距(1:n)*1.2+5,

                cn,cex=1,adj=c(0,0),facing="inside")

  }

},bg.border=NA)

print('添加图例......')
library(circlize)
library(gridBase)

for (i in 1:outcome.num){
lg_Exp=Legend(title=outcome[i],col_fun=mycol[[i]],direction = c("vertical"))
assign(paste("lg_Exp", i, sep = ""), lg_Exp)
}

# 创建图例

circle_size= unit(0.6,"snpc")

h= dev.size()

lgd_list= packLegend(lg_Exp1,lg_Exp2,lg_Exp3, max_height = unit(2*h,"inch"))
if (outcome.num>1){lgd_list= packLegend(lg_Exp1,lg_Exp2, max_height = unit(2*h,"inch"))}
if (outcome.num>2){lgd_list= packLegend(lg_Exp1,lg_Exp2,lg_Exp3, max_height = unit(2*h,"inch"))}
if (outcome.num>3){lgd_list= packLegend(lg_Exp1,lg_Exp2,lg_Exp3,lg_Exp4, max_height = unit(2*h,"inch"))}
draw(lgd_list, x = circle_size, just ="midle")

}




result_MR <- function (folder.path='',
                       exposure.name=''
                       ){
# 指定本地文件夹路径
folder_path <- folder.path

# 获取文件夹中所有.gz文件的文件名
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# 显示文件名列表
print(file_list)


# 创建一个空的数据框用于存储结果
result_df <- data.frame(matrix(ncol = 7, nrow = 0))

print('循环读取文件并建立数据行--------')
for (file in file_list) {
  # 读取文件数据
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)

  # 提取所需的数据
  biomarkers <- data[data$method == "Inverse variance weighted", "id.exposure"]
  IVW_OR <- data[data$method == "Inverse variance weighted", "or"]
  MR_Egger_OR <- data[data$method == "MR Egger", "or"]
  Weighted_median_OR <- data[data$method == "Penalised weighted median", "or"]
  IVW_P <- data[data$method == "Inverse variance weighted", "pval"]
  MR_Egger_P<- data[data$method == "MR Egger", "pval"]
  Weighted_median_P<- data[data$method == "Penalised weighted median", "pval"]
  # 创建当前文件的数据行
  row_data <- c(biomarkers, IVW_OR, MR_Egger_OR, Weighted_median_OR, IVW_P,MR_Egger_P,Weighted_median_P)

  # 将数据行添加到结果数据框中
  result_df <- rbind(result_df, row_data)
}
colnames(result_df) <- c(exposure.name, "IVW_OR", "MR_Egger_OR", "Weighted_median_OR",  "IVW_P",'MR_Egger_P','Weighted_median_P')
print('打印结果数据框')
write.csv(result_df, paste(exposure.name,"_result.csv"), row.names=F)
}

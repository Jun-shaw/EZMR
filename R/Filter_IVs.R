library(MendelianRandomization)
library(TwoSampleMR)
filter_IVs <- function (folder.path='',
                        file.type='',
                        SNP='',
                        beta='',
                        se='',
                        effect='',
                        other='',
                        P='',
                        EAF='',
                        sample.num=50000,
                        Pfilter=5e-08,
                        kb=10000,
                        r2=0.001,
                        Ffilter=10)
{
  # 指定本地文件夹路径
  folder_path <- file.path

  # 获取文件夹中所有.gz文件的文件名
  file_list <- list.files(path = folder_path, pattern = paste0('\\',file.type,'$'), full.names = TRUE)

  # 显示文件名列表
  print(file_list)

  for (file in file_list) {


    exposure_dat <- read_exposure_data(filename = file,
                                       sep = ",",
                                       snp_col = SNP,
                                       beta_col = beta,
                                       se_col = se,
                                       effect_allele_col = effect,
                                       other_allele_col = other,
                                       pval_col = P,
                                       eaf_col = EAF,
                                       samplesize_col=sample.num,
                                       clump = FALSE)
    exposure_dat_filtered<-subset(exposure_dat, pval.exposure<Pfilter)
    options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
    exposure_dat_clumped <- clump_data(exposure_dat_filtered, clump_kb = kb, clump_r2 = r2)
    dir.create('去除连锁不平衡')
    write.csv(exposure_dat_clumped, file = paste0("去除连锁不平衡/", gsub(file.type, "", file), "_LD.csv"), row.names = FALSE)

    #去除弱工具变量
    Ffilter = 10
    dat<-read.csv(paste0("去除连锁不平衡/", gsub(file.type, "", file), "_LD.csv"), header=T, sep=",", check.names=F)
    N=dat[1,"sample.num"]
    dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
    dat=transform(dat,F=(N-2)*R2/(1-R2))
    outTab=dat[dat$F>Ffilter,]
    dir.create('工具变量')
    write.csv(outTab, file = paste0("工具变量/", gsub(file.type, "", file), "_F.csv"), row.names = FALSE)
  }

}

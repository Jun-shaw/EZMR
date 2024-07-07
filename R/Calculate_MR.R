library(MendelianRandomization)
library(TwoSampleMR)
Calculate_MR <- function (folder.path='',
                          outcome.id='',
                          outcome.name='',
                          valid=T,
                          valid_threshold=0.05)
{
folder_path <- folder.path
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
for (file in file_list) {

  exposure_dat<-read_exposure_data(filename = file,
                                   sep = ",",
                                   snp_col = "SNP",
                                   beta_col = "beta.exposure",
                                   se_col = "se.exposure",
                                   effect_allele_col = "effect_allele.exposure",
                                   pval_col = "pval.exposure",
                                   other_allele_col = "other_allele.exposure",
                                   eaf_col = "eaf.exposure",
                                   clump = F)



  outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP,
                                    outcomes=outcome.id, #finn-b-N14_PROSTHYPERPLA #finn-b-N14_PROSTATITIS #finn-b-C3_PROSTATE
                                    proxies = FALSE,
                                    maf_threshold = 0.01,
                                    access_token = NULL)
  #write.csv(outcome_dat, file="outcome_dat.csv", row.names=F)

  dat <- harmonise_data(exposure_dat = exposure_dat,outcome_dat = outcome_dat)

  mrResult=mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_penalised_weighted_median"))
  mrTab=generate_odds_ratios(mrResult)
  #mrTab$id.exposure <- gsub("工具变量/(.*?)\\.imp.*", "\\1", file)
  dir_name <- paste0(outcome.id,'结果')
  dir.create(dir_name)
  write.csv(mrTab, file=paste0(dir_name,"/", gsub("工具变量/", "", file), "_result.csv"), row.names=F)
  }
if (valid ==T) {
  use_package("MRPRESSO")
  library(MRPRESSO)
  mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure",
            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = data, NbDistribution = 1000,
            SignifThreshold = threshold)


  dir.create(outcome.name)
  #异质性检验
  heterTab=mr_heterogeneity(dat)
  write.csv(heterTab, file=paste0(outcome.name,"/", gsub("工具变量/", "", file), "_heterTab.csv"), row.names=F)
  #多效性检验
  pleioTab=mr_pleiotropy_test(dat)
  write.csv(pleioTab, file=paste0(outcome.name,"/", gsub("工具变量/", "", file), "_pleioTab.csv"), row.names=F)
  }

}


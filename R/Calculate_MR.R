#' @name calculate_MR
#' @title Calculate mendelian randomization.
#' @description Numerous MR and sensitivity analyses are employed to validate causal inferences.
#'
#' @import MendelianRandomization
#' @import TwoSampleMR
#' @import MRPRESSO
#'
#' @param folder.path This parameter typically refers to the path to the folder where data files are stored.
#' @param exposure.id Acquire the GWAS project ID from the website https://gwas.mrcieu.ac.uk/.
#' @param Pfilter P-value filtering threshold, used to select significant SNPs. Only SNPs with P-values below this threshold are considered significantly associated.
#' @param kb Kilobases, which may define the window size for analysis in the genome, often used in assessing associations of nearby SNPs.
#' @param r2 This represents the degree of linkage disequilibrium (LD) between SNPs, typically used to evaluate the correlation between SNPs. An r² value closer to 1 indicates a stronger association between two SNPs.
#' @param outcome.id Acquire the GWAS project ID from the website https://gwas.mrcieu.ac.uk/.
#' @param outcome.name The name of the outcome factor.
#' @param valid This parameter is a logical variable that determines whether to conduct sensitivity validation.
#' @param valid.threshold  This parameter sets the significance threshold for sensitivity analysis.
#'
#'#' @references
#' Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases. Nat Genet. 2018 May;50(5):693-698. doi: 10.1038/s41588-018-0099-7. Epub 2018 Apr 23. Erratum in: Nat Genet. 2018 Aug;50(8):1196. doi: 10.1038/s41588-018-0164-2. PMID: 29686387; PMCID: PMC6083837.
#'
#' @author Junxiao Shen
#' @examples
#' filter_IVs (exposure.id=NA, # Filter criteria required when using open GWAS data
#'             folder.path='Instrumental_Variables/',
#'             Pfilter=5e-08,
#'             kb=10000,
#'             r2=0.001,
#'             outcome.id='',
#'             outcome.name='',
#'             valid=T,
#'             valid.threshold=0.05)
#'
#' @export
#'
  calculate_MR <- function(exposure.id = NA, # Filter criteria required when using open GWAS data
                           folder.path = '',
                           Pfilter = 5e-08,
                           kb = 10000,
                           r2 = 0.001,
                           outcome.id = '',
                           outcome.name = '',
                           valid = TRUE,
                           valid.threshold = 0.05)
  {
    if (!is.na(exposure.id)) {
      for (i in exposure.id) {
        extract_instruments(outcomes = i, p < Pfilter, clump = TRUE, r2 = r2, kb = kb, access_token = NULL)
        outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP,
                                            outcomes = outcome.id, # Example: finn-b-N14_PROSTHYPERPLA #finn-b-N14_PROSTATITIS #finn-b-C3_PROSTATE
                                            proxies = FALSE,
                                            maf_threshold = 0.01,
                                            access_token = NULL)
        # write.csv(outcome_dat, file = "outcome_dat.csv", row.names = FALSE)

        dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)

        mrResult = mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_penalised_weighted_median"))
        mrTab = generate_odds_ratios(mrResult)
        # mrTab$id.exposure <- gsub("instrument/(.*?)\\.imp.*", "\\1", file)
        res.folder <- paste0(outcome.name, '_results')
        dir.create(res.folder)
        write.csv(mrTab, file = paste0(res.folder, "/", gsub("instrument/", "", file), "_result.csv"), row.names = FALSE)
      }
    } else {
      file_list <- list.files(path = folder.path, pattern = "\\.csv$", full.names = TRUE)

      for (file in file_list) {
        exposure_dat <- read_exposure_data(filename = file,
                                           sep = ",",
                                           snp_col = "SNP",
                                           beta_col = "beta.exposure",
                                           se_col = "se.exposure",
                                           effect_allele_col = "effect_allele.exposure",
                                           pval_col = "pval.exposure",
                                           other_allele_col = "other_allele.exposure",
                                           eaf_col = "eaf.exposure",
                                           clump = FALSE)
        outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP,
                                            outcomes = outcome.id, # Example: finn-b-N14_PROSTHYPERPLA #finn-b-N14_PROSTATITIS #finn-b-C3_PROSTATE
                                            proxies = FALSE,
                                            maf_threshold = 0.01,
                                            access_token = NULL)
        # write.csv(outcome_dat, file = "outcome_dat.csv", row.names = FALSE)

        dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)

        mrResult = mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_penalised_weighted_median"))
        mrTab = generate_odds_ratios(mrResult)
        # mrTab$id.exposure <- gsub("Instrumental_Variables/(.*?)\\.imp.*", "\\1", file)
        res.folder <- paste0(outcome.name, '_results')
        dir.create(res.folder)
        write.csv(mrTab, file = paste0(res.folder, "/", gsub("Instrumental_Variables/", "", file), "_result.csv"), row.names = FALSE)
      }
    }

    if (valid == TRUE) {
      mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
                NbDistribution = 1000, SignifThreshold = valid.threshold)

      val.folder <- paste0(outcome.name, '_valid')
      dir.create(val.folder)
      # Heterogeneity test
      heterTab = mr_heterogeneity(dat)
      write.csv(heterTab, file = paste0(val.folder, "/", gsub("Instrumental_Variables/", "", file), "_heterTab.csv"), row.names = FALSE)
      # Pleiotropy test
      pleioTab = mr_pleiotropy_test(dat)
      write.csv(pleioTab, file = paste0(val.folder, "/", gsub("Instrumental_Variables/", "", file), "_pleioTab.csv"), row.names = FALSE)
    }
  }



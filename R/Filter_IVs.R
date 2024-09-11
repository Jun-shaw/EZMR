#' @name filter_IVs
#' @title Refining instrumental variables for exposure factors.
#' @description Enhancing the instrumental variables pertinent to exposure factors is essential for ensuring robust and reliable analysis. This refinement process entails carefully selecting and fine-tuning the variables that serve as the instruments, thereby bolstering their validity and effectiveness in elucidating the causal relationships within the study. Such meticulous attention to detail will ultimately yield more accurate and insightful conclusions concerning the impact of exposure factors on the outcomes of interest.
#'
#' @import MendelianRandomization
#' @import TwoSampleMR
#'
#' @param folder.path This parameter typically refers to the path to the folder where data files are stored.
#' @param file.type This indicates the format of the input files (e.g., CSV, gz), which is necessary for correctly reading the data.
#' @param SNP This refers to Single Nucleotide Polymorphisms, which are variations at a single nucleotide position in the genome.
#' @param beta This parameter usually represents the effect size of the SNP on the trait, indicating the magnitude of the SNP's influence on the phenotype. It is often the regression coefficient in a linear regression model.
#' @param se Standard Error, which measures the precision of the beta estimate. A smaller standard error indicates a more reliable estimate of the effect size.
#' @param effect This parameter indicates the allele (variant) that is associated with an increase in the risk of the trait or disease being studied. It is the allele whose presence is linked to a higher value of the outcome variable (e.g., increased disease risk or trait value). .
#' @param other This parameter represents the allele that is considered the alternative or non-risk allele in comparison to the effect_allele.exposure. While the effect allele is associated with an increase in the trait or disease risk, the other allele is typically associated with a lower risk or a different phenotype.
#' @param P The P-value, used to assess the significance of the association between the SNP and the trait. Smaller P-values indicate stronger statistical significance.
#' @param sample.num Sample size, which refers to the number of individuals used in the GWAS analysis. Larger sample sizes generally increase the statistical power of the results.
#' @param Pfilter P-value filtering threshold, used to select significant SNPs. Only SNPs with P-values below this threshold are considered significantly associated.
#' @param kb Kilobases, which may define the window size for analysis in the genome, often used in assessing associations of nearby SNPs.
#' @param r2 This represents the degree of linkage disequilibrium (LD) between SNPs, typically used to evaluate the correlation between SNPs. An r² value closer to 1 indicates a stronger association between two SNPs.
#' @param Ffilter F-statistic filtering, often used to assess whether SNPs have sufficient variability in a specific sample. SNPs with higher F-values are generally considered more reliable.
#'
#'#' @references
#' Hemani G, Zheng J, Elsworth B, Wade KH, Baird D, Haberland V, Laurin C, Burgess S, Bowden J, Langdon R, Tan VY, Yarmolinsky J, Shihab HA, Timpson NJ, Evans DM, Relton C, Martin RM, Davey Smith G, Gaunt TR, Haycock PC, The MR-Base Collaboration. The MR-Base platform supports systematic causal inference across the human phenome. eLife 2018;7:e34408. doi: 10.7554/eLife.34408
#'
#' @author Junxiao Shen
#' @examples
#' filter_IVs (folder.path='document/',
#'             file.type='.gz',
#'             SNP='SNP',
#'             beta='beta.exposure',
#'             se='se.exposure',
#'             effect='effect_allele.exposure',
#'             other='other_allele.exposure',
#'             P='pval.exposure',
#'             EAF='eaf.exposure',
#'             sample.num=40000,
#'             Pfilter=5e-08,
#'             kb=10000,
#'             r2=0.001,
#'             Ffilter=10)
#'
#' @export
#'
filter_IVs <- function (folder.path='',
                        file.type='',
                        SNP='',
                        beta='',
                        se='',
                        effect='',
                        other='',
                        P='',
                        EAF='',
                        sample.num=40000,
                        Pfilter=5e-08,
                        kb=10000,
                        r2=0.001,
                        Ffilter=10)
{
  # Retrieve the filenames of all .gz files in the folder
  file_list <- list.files(path = folder_path, pattern = paste0('\\', file.type), full.names = TRUE)

  # Display the list of filenames
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
                                       samplesize_col = sample.num,
                                       clump = FALSE)

    exposure_dat_filtered <- subset(exposure_dat, pval.exposure < Pfilter)

    options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
    exposure_dat_clumped <- clump_data(exposure_dat_filtered, clump_kb = kb, clump_r2 = r2)

    dir.create('Removed_Linkage_Disequilibrium')
    write.csv(exposure_dat_clumped, file = paste0("Removed_Linkage_Disequilibrium/", gsub(file.type, "", file), "_LD.csv"), row.names = FALSE)

    # Eliminate weak instrumental variables
    Ffilter <- 10
    dat <- read.csv(paste0("Removed_Linkage_Disequilibrium/", gsub(file.type, "", file), "_LD.csv"), header = TRUE, sep = ",", check.names = FALSE)

    N <- dat[1, "sample.num"]
    dat <- transform(dat, R2 = 2 * ((beta.exposure)^2) * eaf.exposure * (1 - eaf.exposure))
    dat <- transform(dat, F = (N - 2) * R2 / (1 - R2))

    outTab <- dat[dat$F > Ffilter, ]
    dir.create('Instrumental_Variables')
    write.csv(outTab, file = paste0("Instrumental_Variables/", gsub(file.type, "", file), "_F.csv"), row.names = FALSE)
  }
}

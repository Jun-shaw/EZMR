
library(devtools)
library(roxygen2)
library(EZMR)
usethis::use_r("Filter_IVs")
usethis::create_package("adasdas")
devtools::load_all()

filter_IVs(folder.path = ,
            )

Calculate_MR(folder.path = '',outcome.id = '',outcome.name = '',valid = F, valid_threshold = 0.05)

Result_MR(folder.path = '',exposure.name='')

visualize_MR(filename = 'data/前列腺炎增癌result.csv',
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

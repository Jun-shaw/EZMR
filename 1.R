library(devtools)
library(roxygen2)
library(EZMR)
install_github("Jun-shaw/EZMR")
usethis::use_r("Forest_plot_MR")
devtools::load_all()
filter_IVs(folder.path='',
           file.type='',
           SNP='',
           beta='',
           se='',
           effect='',
           other='',
           P='',
           EAF='',
           sample.num=,
           Pfilter=5e-08,
           kb=10000,
           r2=0.001,
           Ffilter=10)

Calculate_MR(folder.path='',
             exposure.id=NA,#如果使用open GWAS数据则需要下列输入筛选条件
             Pfilter=5e-08,
             kb=10000,
             r2=0.001,
             outcome.id='',
             outcome.name='',
             valid=T,
             valid_threshold=0.05)

Result_MR(folder.path = '',exposure.name='')

Visualize_MR(filename = 'data/前列腺炎增癌result.csv',
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
              category_colors=c('#CC88B0', '#DBE0ED', '#87B5B2', '#F4CEB4', '#F1DFA4', '#998DB7')
             )

Forest_plot_MR(file = 'data/前列腺癌result.csv',
                column_num=10,
                shape_col='#4575b4',
                line_col='#CE5C69',
                text_col='#4575b4',
                P_threshold=0.05,
                OR_range = c(0,2),
                cutoff = c(0,1,2)
               )
devtools::check()

PhyloseqDATA
================
KentoLogan
11/20/2020

``` r
BiocManager::install("phyloseq")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.2 (2020-06-22)

    ## Installing package(s) 'phyloseq'

    ## package 'phyloseq' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\arakn\AppData\Local\Temp\RtmpyEA7Yr\downloaded_packages

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, MASS, mgcv, nlme, survival

    ## Old packages: 'cli', 'colorspace', 'digest', 'lubridate', 'magrittr', 'pillar',
    ##   'rlang', 'rprojroot', 'tibble', 'vctrs', 'xfun'

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.2
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.4.0     v forcats 0.5.0

    ## Warning: package 'readr' was built under R version 4.0.3

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(phyloseq)
library(RColorBrewer)
```

    ## Warning: package 'RColorBrewer' was built under R version 4.0.3

\#IMPORT DATA

``` r
count.tab<- read_rds("~/GitHub_Files/144l_students/Output_Data/seqtab-nochimtaxa.rds")
#TableOContents for each seq in each samp
tax.tab<- read_rds("~/GitHub_Files/144l_students/Output_Data/taxa.rds")
#Table matches ASV to seq
sample.tab<- read_rds("~/GitHub_Files/144l_students/Output_Data/DOC_BGE.rds") %>% 
  drop_na(DNA_SampleID) %>% 
  column_to_rownames("DNA_SampleID")
```

\#PHYLOSEQ OBJECT: Merge the Three Objects

``` r
OTU= otu_table(count.tab, taxa_are_rows = T)
TAX= tax_table(tax.tab)
SAM= sample_data(sample.tab)
ps= phyloseq(OTU, TAX, SAM)
```

\#FILTER SEQ: Filter Out Chloroplasts and Mitochondria

``` r
sub_ps<- ps %>% 
  subset_taxa(Family!= "mitochondria" & Order!="Chloroplast")
```

\#SAMPLE SUMMARY

``` r
sample_sum_df<- data.frame(sum= sample_sums(sub_ps))
ggplot(sample_sum_df, aes(x= sum))+
  geom_histogram(color= "black", fill= "#377EB8", binwidth = 1000)+
  ggtitle("Distribution of Sample Sequencing Depth")+
  xlab("Read Counts")+
  theme(axis.title.y = element_blank())+
  theme_bw()
```

![](PHYLOSEQ_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

\#BETA SAMPLING \#\#Subsample

``` r
ps_min<- rarefy_even_depth(sub_ps, sample.size = min(sample_sums(sub_ps)))
```

    ## You set `rngseed` to FALSE. Make sure you've set & recorded
    ##  the random seed of your session for reproducibility.
    ## See `?set.seed`

    ## ...

    ## 139OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
mean(sample_sums(sub_ps))
```

    ## [1] 32754.79

``` r
mean(sample_sums(ps_min))
```

    ## [1] 2487

\#\#NMDS

``` r
set.seed(1)
nmds<- ordinate(sub_ps, method = "NMDS", distance = "bray")
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.07118227 
    ## Run 1 stress 0.07308537 
    ## Run 2 stress 0.07146138 
    ## ... Procrustes: rmse 0.006885389  max resid 0.02456195 
    ## Run 3 stress 0.07118227 
    ## ... Procrustes: rmse 3.270612e-05  max resid 0.0001269942 
    ## ... Similar to previous best
    ## Run 4 stress 0.1136674 
    ## Run 5 stress 0.07118229 
    ## ... Procrustes: rmse 3.161281e-05  max resid 8.46401e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.08046092 
    ## Run 7 stress 0.1482789 
    ## Run 8 stress 0.07307733 
    ## Run 9 stress 0.1480984 
    ## Run 10 stress 0.1529913 
    ## Run 11 stress 0.1135783 
    ## Run 12 stress 0.07159431 
    ## ... Procrustes: rmse 0.008425193  max resid 0.02500106 
    ## Run 13 stress 0.0711823 
    ## ... Procrustes: rmse 3.81204e-05  max resid 0.0001254173 
    ## ... Similar to previous best
    ## Run 14 stress 0.08028487 
    ## Run 15 stress 0.0715943 
    ## ... Procrustes: rmse 0.008424011  max resid 0.02499544 
    ## Run 16 stress 0.07118227 
    ## ... Procrustes: rmse 2.216257e-05  max resid 7.701125e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.08123803 
    ## Run 18 stress 0.07118226 
    ## ... New best solution
    ## ... Procrustes: rmse 2.117723e-05  max resid 8.032869e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.08036031 
    ## Run 20 stress 0.07118226 
    ## ... New best solution
    ## ... Procrustes: rmse 1.07449e-05  max resid 3.979448e-05 
    ## ... Similar to previous best
    ## *** Solution reached

``` r
set.seed(1)
nmds_min<- ordinate(ps_min, method = "NMDS", distance = "bray")
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.0906512 
    ## Run 1 stress 0.1030418 
    ## Run 2 stress 0.09171425 
    ## Run 3 stress 0.09078387 
    ## ... Procrustes: rmse 0.01338821  max resid 0.04115407 
    ## Run 4 stress 0.09005025 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01855108  max resid 0.06416527 
    ## Run 5 stress 0.1615297 
    ## Run 6 stress 0.09005025 
    ## ... New best solution
    ## ... Procrustes: rmse 2.50495e-06  max resid 9.365734e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.1706283 
    ## Run 8 stress 0.09044209 
    ## ... Procrustes: rmse 0.01047307  max resid 0.03609418 
    ## Run 9 stress 0.1699919 
    ## Run 10 stress 0.1340213 
    ## Run 11 stress 0.1352634 
    ## Run 12 stress 0.09171467 
    ## Run 13 stress 0.09005025 
    ## ... Procrustes: rmse 6.332141e-06  max resid 1.865889e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.1669637 
    ## Run 15 stress 0.09044209 
    ## ... Procrustes: rmse 0.01047281  max resid 0.03609364 
    ## Run 16 stress 0.09069996 
    ## Run 17 stress 0.1031393 
    ## Run 18 stress 0.1335853 
    ## Run 19 stress 0.1024304 
    ## Run 20 stress 0.1031393 
    ## *** Solution reached

``` r
levels<-c("Control", "Ash Leachate", "Mud Leachate", "Glucose_Nitrate_Phosphate")

nmds.plot<- plot_ordination(sub_ps, nmds, title = "NMDs")+
  geom_point(aes(fill= days, shape= factor(Treatment, levels = levels)), alpha= 0.6, stroke= 2, size= 4)+
  scale_shape_manual(values = c(21,22,23))+
  scale_fill_gradient(low = "#0db5e6", high= "#d31f2a")+
  theme_bw()

nmds.plot$layers<- nmds.plot$layers[1]
nmds.plot+
  facet_grid()+
  guides(fill= guide_colorbar(title = "Days"), shape= guide_legend(title = "Treatment"))
```

![](PHYLOSEQ_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
levels<-c("Control", "Ash Leachate", "Mud Leachate", "Glucose_Nitrate_Phosphate")

nmds_min.plot<- plot_ordination(sub_ps, nmds_min, title = "NMDs")+
  geom_point(aes(fill= days, shape= factor(Treatment, levels = levels)), alpha= 0.6, stroke= 2, size= 4)+
  scale_shape_manual(values = c(21,22,23))+
  scale_fill_gradient(low = "#0db5e6", high= "#d31f2a")+
  theme_bw()
  

nmds_min.plot$layers<- nmds_min.plot$layers[1]
nmds_min.plot+
  facet_grid()+
  guides(fill= guide_colorbar(title = "Days"), shape= guide_legend(title = "Treatment"))
```

![](PHYLOSEQ_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

ALPHA DIVERSITY

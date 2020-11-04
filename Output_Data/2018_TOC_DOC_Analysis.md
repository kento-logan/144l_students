2018\_TOC\_Data
================
KentoLogan
11/3/2020

\#Load Libraries

``` r
library(tidyverse)
```

    ## -- Attaching packages ------------------------------------------------------------------ tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.2
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts --------------------------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(readxl)
library(lubridate)
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

\#Load Data

``` r
excel_sheets("~/GitHub_Files/144l_students/Input_Data/week4/144L_2018_Exp_TOC.xlsx")
```

    ## [1] "Metadata" "Data"

``` r
metadata<- read_excel("~/GitHub_Files/144l_students/Input_Data/week4/144L_2018_Exp_TOC.xlsx", sheet = "Metadata")
data<- read_excel("~/GitHub_Files/144l_students/Input_Data/week4/144L_2018_Exp_TOC.xlsx", sheet = "Data")
```

\#Prep Data

``` r
joined<- left_join(metadata, data)
```

    ## Joining, by = c("Bottle", "Timepoint")

``` r
toc<- joined %>% 
  mutate(Datetime= ymd_hm(Datetime)) %>% 
  group_by(Bottle) %>% 
  mutate(interv= interval(first(Datetime), Datetime),
         hours= interv/3600,
         days= hours/24) %>% 
  ungroup() %>% 
  select(Bottle:Datetime, hours, days, everything())

subset<- toc %>% 
  select(Bottle, Datetime, interv)
```

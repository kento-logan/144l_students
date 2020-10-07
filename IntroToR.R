##### Intro to R: CalFire#####
#Kento Logan#
#10/6/20#

library(tidyverse)
library(readxl)

##### LOAD DATA SET#####
excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")


calfire.metadata<- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet= "Metadata")
calfire.data<- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx",sheet = 2)

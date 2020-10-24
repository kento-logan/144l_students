####THOMAS FIRE PROGRESSION####
#10/24/20 KentoLogan

####Packages####
library(tidyverse)
library(dplyr)
library(readxl)
install.packages("praise")
library(praise)
library(lubridate)

#####Data Sets#####
excel_sheets("GitHub_Files/144l_students/Input_Data/week1/Thomas_Fire_Progression.xlsx")

data<- read_excel("GitHub_Files/144l_students/Input_Data/week1/Thomas_Fire_Progression.xlsx", sheet = "Data")
meta<- read_excel("GitHub_Files/144l_students/Input_Data/week1/Thomas_Fire_Progression.xlsx", sheet = "Metadata")

#####Data Org#####
thomas_fire<- data %>%
  select(Date:PM25) %>% 
  arrange(desc(Date)) %>% 
  filter(Acres_Burned >= 30000)

ggplot(thomas_fire,aes(x = Date, y = Acres_Burned))+ geom_point()+geom_line()+
  ggtitle("Kento Logan \nThomas Fire Data \n 12/5/17-1/12/18")




##### Intro to R: CalFire#####
#Kento Logan#
#10/6/20#

library(tidyverse)
library(readxl)

##### LOAD DATA SET#####
excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")

calfire.metadata<- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet= "Metadata")
calfire.data<- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx",sheet = 2)

--------------------------------------------------------------------------------
##### Initial DATA Exploration#####
names(calfire.data) #shows column names
dim(calfire.data)
#Rows/obs, columns
class(calfire.data) #data type/class
head(calfire.data) #shows first 6 lines
tail(calfire.data)
#shows last 6 lines

#single columns can be referred to with a $ 
county <- calfire.data$County_Unit
county

names(calfire.data)
max_acres <-  max(calfire.data$Total_Acres_Burned, na.rm = T)
 #NA example:
max(calfire.data$Structures_Destroyed)
max(calfire.data$Structures_Destroyed, na.rm = T)
 -------------------------------------------------------------------------------
#####Basic DATA Wrangling#####
df1 <- select(calfire.data, County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities )
view(df1)

df2 <- filter(df1, County_Unit%in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name == "THOMAS")
view(df2)

df3 <- arrange(df2, desc(Start_Date), Total_Acres_Burned)
view(df3)

df4 <- mutate_at(df3, vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0)
view(df4)

df5 <- mutate(df4, Fatalities = Fire_Fatalities + Civil_Fatalities)
view(df5)
--------------------------------------------------------------------------------
#####THE WORLD: TIME#####
install.packages("lubridate")
library(lubridate)

library(tidyverse)
df6 <- mutate(df5, interv = interval(Start_Date, Controlled_Date), dur = as.duration(interv), days = as.numeric(dur, "days"))
--------------------------------------------------------------------------------
#15Lines to do 5 dataframes. seems pretty inefficient. A better way is piping
#####Intro to Piping#####

#(control+shift+m) = Pipe Operator

socal.fires<- calfire.data %>%
  select(County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities ) %>%
  filter(County_Unit%in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name == "THOMAS") %>%
  arrange(desc(Start_Date), Total_Acres_Burned) %>% 
  mutate_at(vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0) %>% 
  mutate(Fatalities = Fire_Fatalities + Civil_Fatalities) %>% 
  mutate(Struct_Impact = Structures_Damaged + Structures_Destroyed) %>% 
  mutate(interv = interval(Start_Date, Controlled_Date), dur = as.duration(interv), days = as.numeric(dur, "days"))
  ------------------------------------------------------------------------------
  #####Intro to GGPLOT#####

ggplot(socal.fires, aes(x= Start_Date, y= Total_Acres_Burned)) + geom_point(aes(color = County_Unit))+ ggtitle("CA South Coast Major Fires \n2014-2018")+ labs(x= "", y= "Acres Burned", color = "County")+ theme_bw()+ facet_grid(rows = "County_Unit", scales = "free")

incidents<- socal.fires %>% 
  rename(county = County_Unit)

plot.data<- socal.fires %>% 
  rename(county=County_Unit, acres=Total_Acres_Burned,start=Start_Date, end=Controlled_Date) %>%
  mutate(year=year(start), county=ifelse(county=="VENTURA/SANTA BARBARA", "VENTURA", county))

incidents<-plot.data %>% 
  group_by(county, year) %>% 
 tally() 
  ungroup()

  all_incidents.plot<- all_incidents %>%
    ggplot(aes(x=year, y=n)) + 
    geom_point(color= "blue") + 
    geom_line(color = "blue") + 
    labs(title="CA South Coast Major Fire Incidents \n 2014-2018", x="", y= "Incidents")+ theme_bw()
 
  all_incidents<- plot.data %>% 
    group_by(year) %>% 
    tally() %>% 
    ungroup()
  
  #####SAVE DATA AND PLOTS#####
  #cntrl s 
  saveRDS(socal.fires, file = "Output_Data/socalfires.data.rds")
  write.csv(socal.fires, "Output_Data/socal.fies.data.csv")

  ggsave(filename="Fire_Incidents", all_incidents.plot, device ="jpeg", "Output_Data/")  
  
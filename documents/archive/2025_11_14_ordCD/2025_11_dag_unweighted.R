rm(list=ls())     # Remove everything from environment

pacman::p_load(
  OrdCD, 
  igraph,
  haven,      # importing datasets 
  here,       # using relative file paths affilated w/ project
  dplyr,      # selecting and renaming variables
  magrittr, 
  openxlsx,    #export to excel
  
  tidyverse,
  gtsummary,   # Creating tables
  janitor,
  flextable,
  skimr,
  rio, 
  kableExtra,
  
  
  gtsummary,  # summary statistics and tests 
  tidyverse,  # piping, graphics, labeling, etc
  rstatix,    # summary statistics and statistical tests
  janitor,     #displaying tabular data
  
  ggplot2,
  ggthemes,
  ggsci, 
  
  survey,
  svyVGAM
  
)

##### 1: read in survey data and format file #####
dat<-haven::read_sas('R:\\Hausmann_IRBNET1628373_SDOHSHEP\\Data\\YRV_2023\\R_proj\\YRV_analysis_2023_DF\\SHEP Paper 2\\raw_data\\yrv_fy23_pilot.sas7bdat' # data file
                     ,'R:\\Hausmann_IRBNET1628373_SDOHSHEP\\Data\\YRV_2023\\R_proj\\YRV_analysis_2023_DF\\SHEP Paper 2\\raw_data\\formats_th.sas7bdat' # format file
) 


##### 2: create new vars #####

# a. Jen's combined race-ethnicity variable
race_eth <- read.csv('R:\\Hausmann_IRBNET1628373_SDOHSHEP\\Data\\YRV_2023\\SPSS\\race and ethnicity combined variable csv.csv')
dat$SURVID <- as.numeric(dat$SURVID)
dat <- left_join(dat, race_eth, by = c('SURVID'))

dat<-dat %>%
  ##Age
  mutate(age_v1 = case_when(    #combine first two groups
    Age_T1 == 1 ~   1
    ,Age_T1 == 2 ~ 1
    ,Age_T1 == 3 ~ 2
    ,Age_T1 == 4 ~ 3
    ,Age_T1 == 5 ~ 4
    ,Age_T1 == 6 ~ 5
    ,Age_T1 == 7 ~ 6
    ,TRUE        ~ NA_real_),
    age_v2 = case_when(    #combine first three groups
      Age_T1 == 1 ~ 1
      ,Age_T1 == 2 ~ 1
      ,Age_T1 == 3 ~ 1
      ,Age_T1 == 4 ~ 2
      ,Age_T1 == 5 ~ 3
      ,Age_T1 == 6 ~ 4
      ,Age_T1 == 7 ~ 5
      ,TRUE        ~ NA_real_),
    #gender
    gender_v1= case_when(    
      SOGI1_T1==1 ~ 1
      ,SOGI1_T1==2 ~ 2
      ,SOGI1_T1>2  ~ 3
      ,TRUE        ~ NA_real_),
    
    #sexual orientation     
    sexor_v1 = case_when(
      SOGI2_T1==1    ~ "Straight"  #Heterosexual
      ,SOGI2_T1==2
      |SOGI2_T1==3
      |SOGI2_T1==4
      |SOGI2_T1==5 
      |SOGI2_T1==6 ~"LGB+" #LGB, other, not sure
      ,TRUE        ~ NA),
    #sex      
    
    sex_admin_num = if_else( #creating numeric sex
      SEX== "F", 1,2),
    
    BirthCert_T1= ifelse(   #converting NA to 0
      is.na(BirthCert_T1),0,BirthCert_T1),   
    
    sex_complete= case_when( #creating new sex var
      BirthCert_T1==1  ~ 1  #Female
      ,BirthCert_T1==2 ~ 2  #Male
      ,BirthCert_T1== 0 & sex_admin_num==1 ~1
      ,BirthCert_T1==0 & sex_admin_num ==2 ~2),  
    
    BirthCert_T1= na_if(BirthCert_T1,0), #convert 0 to NA
    #race and ethnicity 
    ##hispanic
    hispanic= case_when(
      RaceEthcombined==7  ~1  #Hispanic self report
      ,RaceEthcombined==0     #missing self report, use admin
      &Ethnicity=="Hispanic or Latino" ~ 1),
    ##nh black
    race_black = case_when(
      RaceEthcombined==4 ~ 1  #nh black self report
      ,RaceEthcombined==8     #race missing not hispanic, use admin
      & Race=="Black or African American" ~1 
      ,RaceEthcombined==0     #race/eth missing use admin
      & Race=="Black or African American"
      & Ethnicity!="Hispanic or Latino"   ~1
      ,RaceEthcombined==1 #multi use if not black and white
      &Race_AfricanAmerican_T1==3
      &Race_White_T1!=5 ~1),
    ##nh white
    race_white= case_when(
      RaceEthcombined==2 ~ 1 #nh white self report
      ,RaceEthcombined==8    #race missing not hispanic, use admin
      & Race=="White" ~1 
      ,RaceEthcombined==0   #race/eth missing use admin
      & Race=="White"
      & Ethnicity!="Hispanic or Latino" ~1
      ,RaceEthcombined==1   #multi use if not black and white
      &Race_White_T1==5
      &Race_AfricanAmerican_T1!=3 ~1),
    ##multi/other
    race_multiother= case_when(
      RaceEthcombined==3     #self report as PI,Asian, AI
      |RaceEthcombined==5
      |RaceEthcombined==6 ~1
      ,RaceEthcombined==1   #self report as both black and white
      &Race_White_T1==5
      &Race_AfricanAmerican_T1==3 ~1
      ,RaceEthcombined==8   #one self report as unknown but admin hispanic
      &Ethnicity=="Hispanic or Latino" ~1
      ,RaceEthcombined==1   #this case covers two people in dat_comp
      &Race_PacificIslander_T1==4 ~1),
    ##combine into one category
    raceEthcombined_v2= case_when(
      hispanic==1 ~ 1
      ,race_white==1 ~2
      ,race_black==1 ~3
      ,race_multiother==1 ~4
      ,TRUE         ~ NA_real_),
    
    RaceEthcombined= na_if(RaceEthcombined,0), #converting NA to 0
    #combining new race and sex vars
    race_sex= case_when(
      sex_complete==1 
      & raceEthcombined_v2==1 ~1 #female hispanic
      ,sex_complete==1
      & raceEthcombined_v2==2 ~2 #female nh white
      ,sex_complete==1
      & raceEthcombined_v2==3 ~3 #female nh black
      ,sex_complete==1
      & raceEthcombined_v2==4 ~4 #female multi/other
      ,sex_complete==2 
      & raceEthcombined_v2==1 ~5 #male hispanic
      ,sex_complete==2 
      & raceEthcombined_v2==2 ~6 #male white
      ,sex_complete==2 
      & raceEthcombined_v2==3 ~7 #male black
      ,sex_complete==2 
      & raceEthcombined_v2==4 ~8), #male multi/other
    #removing female/male multi/other for new race_sex
    race_sex_v2=case_when(
      race_sex==1   ~"Hispanic women" #female hispanic
      ,race_sex==2  ~"White women" #female nh white
      ,race_sex==3  ~"Black women" #female nh black
      ,race_sex==5  ~"Hispanic men" #male hispanic
      ,race_sex==6  ~"White men" #male white
      ,race_sex==7  ~"Black men" #male black
      ,TRUE ~ NA),
    #education
    edu_v1= case_when(
      School_T1== 1
      |School_T1==2
      |School_T1==3  ~1 #high school or less
      ,School_T1==4    ~2 #some college or 2-year
      ,School_T1==5    ~3 #4 year college
      ,School_T1==6    ~4 # more than 4 year college
      ,TRUE            ~ NA_real_),
    
    #age_v3 added 3.12.24
    #complete age_v1
    age_v3 = case_when(
      is.na(age_v1)
      & Age_cat_th==2   ~1
      ,is.na(age_v1)
      & Age_cat_th==3   ~2
      ,is.na(age_v1)
      & Age_cat_th==4   ~3
      ,is.na(age_v1)
      & Age_cat_th==5   ~4
      ,is.na(age_v1)
      & Age_cat_th==6   ~5
      ,is.na(age_v1)
      & Age_cat_th==7   ~6
      ,TRUE              ~age_v1),
    #age_v4 added 3.12.24
    #complete age_v2
    age_v4 = case_when(
      age_v3 ==1 ~ 1
      ,age_v3==2 ~ 1
      ,age_v3==3 ~ 2
      ,age_v3==4 ~ 3
      ,age_v3==5 ~ 4
      ,age_v3==6 ~ 5
      ,TRUE      ~ NA_real_
      
    ))


##### 3: read in data dictonaries ####
dat_dict <- readxl::read_xls("R:\\Hausmann_IRBNET1628373_SDOHSHEP\\Data\\YRV_2023\\R_proj\\YRV_analysis_2023_DF\\SHEP Paper 2\\data dictionaries\\data_dictionary_fy23t05_th_forSDOH.xls")

vars <- dat_dict |> 
  filter(str_detect(`Version 1 - Short`, "27|28|29|30|31|32|33|34|35|36|37|38|39|40|42|43|44|45|46|47|48|49")) |> 
  mutate(Name = case_when(is.na(TH_name) ~ Name, 
                          .default = TH_name)) |> 
  pull(Name)

dat_f<-dat%>%
  mutate_at(vars,factor)

#narrow down columns to those of interest
dat_f<-dat_f%>%
  select(all_of(vars)
         ,SURVID
         ,respond_source
         ,sta5a
         ,sta3n
         ,visn
         ,WEIGHT
         ,COMP_TH, 
         ,age_v4 ## uses self report and admin to complete missing from age_v2 (age category: combines first three levels of Age_T1)
         ,sexor_v1 ## sexual orientation: combines SOGI2_T1 levels 2,3,4,5,6
         ,race_sex_v2 ## combines raceEthcombined_v2 and sex_complete
  )

##### 4: create analytic dataset ####
dat_comp<-dat_f%>%
  filter(COMP_TH==1)  


need_vars <-c('HelpBasics_T1'
              ,'HelpCaregive_T1'
              ,'HelpChildcare_T1'
              ,'HelpWork_T1'
              ,'HelpFood_T1'
              ,'HelpHousing_T1'
              ,'HelpTransport_T1'
              ,'HelpInternet_T1'
              ,'HelpIsolation_T1'
              ,'HelpLoneliness_T1'
              ,'HelpDiscrim_T1'
              ,'HelpLegal_T1'
              ,'HelpEduc_T1'
              ,'SURVID'
)

risk_vars <- c("PayBasics_T1"               
               ,"GetCaregiving_T1"           
               ,"GetChildcare_T1"
               ,"FindWork_T1"                
               ,"LackFood_T1"
               ,"HaveHousing_T1"             
               ,"LackTransport_T1"
               ,"AccessInternet_T1"          
               ,"LackSocial_T1" 
               ,"FeelLonely_T1"              
               ,"DiscrimREC_T1" 
               ,"DiscrimGend_T1"             
               ,"DiscrimAge_T1"
               ,"DiscrimOther_T1"     
               ,"LegalNeed_T1"               
               ,"LackEduc_T1"
               ,'SURVID'
)

sat_vars <- c("WBS_FullySat_T1"
              ,"WBS_Involved_T1"     
              ,"WBS_Function_T1"
              ,'SURVID'
)

health_vars <- c("RateHealth_T1"
                 ,"RateMentalHealth_T1"
                 ,'SURVID'
)



#indicator for missing grid qs

## needs 
outcomes.need<-dat_comp%>%
  select(all_of(need_vars))

outcomes.need$miss_sum<-apply(X = is.na(outcomes.need)
                              , MARGIN = 1
                              , FUN = sum)

outcomes.risk<-dat_comp%>%
  select(all_of(risk_vars))

outcomes.risk$miss_sum<-apply(X = is.na(outcomes.risk)
                              , MARGIN = 1
                              , FUN = sum)


## satisfaction
outcomes.sat<-dat_comp%>%
  select(all_of(sat_vars))

outcomes.sat$miss_sum<-apply(X = is.na(outcomes.sat)
                             , MARGIN = 1
                             , FUN = sum)
## health
outcomes.health<-dat_comp%>%
  select(all_of(health_vars))

outcomes.health$miss_sum<-apply(X = is.na(outcomes.health)
                                , MARGIN = 1
                                , FUN = sum)


#indicator for missing all grid qs

## needs
outcomes.need<-outcomes.need%>%
  mutate(miss_all_need= if_else(
    miss_sum==13,1,0         
  ))

miss_grid_need<-outcomes.need%>%
  select(SURVID ,miss_all_need)


## risks
outcomes.risk<-outcomes.risk%>%
  mutate(miss_all_risk= if_else(
    miss_sum==13,1,0         
  ))

miss_grid_risk<-outcomes.risk%>%
  select(SURVID ,miss_all_risk)

## satisfaction
outcomes.sat<-outcomes.sat%>%
  mutate(miss_all_sat= if_else(
    miss_sum==3,1,0         
  ))

miss_grid_sat<-outcomes.sat%>%
  select(SURVID ,miss_all_sat)

## health
outcomes.health<-outcomes.health%>%
  mutate(miss_all_health= if_else(
    miss_sum==2,1,0         
  ))

miss_grid_health<-outcomes.health%>%
  select(SURVID ,miss_all_health)




#final dataset joined with missing all grid

dat_fin<-dat_comp%>%
  left_join(miss_grid_need,by =c("SURVID"="SURVID")) |> 
  left_join(miss_grid_risk,by =c("SURVID"="SURVID")) |> 
  left_join(miss_grid_sat,by =c("SURVID"="SURVID")) |> 
  left_join(miss_grid_health,by =c("SURVID"="SURVID")) 

test<-dat_fin%>% tabyl(miss_all_need)

#create final dataset for analysis with analysis indicator
#analysis = remove multi racial and all missing need, all missing health, and all missing satisfaction

dat_fin<-dat_fin%>%
  mutate(analysis=case_when(
    miss_all_need ==1 ~0
    ,miss_all_risk ==1 ~0
    ,miss_all_health ==1 ~0
    ,miss_all_sat ==1 ~0
    , TRUE ~1
  ))%>%
  select(-contains("miss_")) 


dat_fin <- dat_fin |> 
  rowwise() |> 
  mutate(WBS_total = sum(10*(as.numeric(WBS_FullySat_T1) - 1),
                         10*(as.numeric(WBS_Involved_T1) - 1), 
                         10*(as.numeric(WBS_Function_T1) - 1), na.rm = TRUE )/3) |> 
  ungroup() |> 
  mutate(RateMentalHealth_T1 = factor(RateMentalHealth_T1, 
                                      levels = c(1, 2, 3, 4, 5), 
                                      labels = c("4", "3", "2", "1", "0")))


rm(list = ls()[!ls() %in% c("dat_fin")])


####

needs <- dat_fin |> dplyr::select(contains(c("Help", "SURVID"))) |> drop_na()
names(needs) <- c("basics", 
                  "caregiving", "childcare", "work", "food", "housing", 
                  "transport", "internet", "isolation", "loneliness", "discrimination", 
                  "legal", "education", "survid")

outcomes <- dat_fin |> 
  rowwise() |> 
  mutate(wellbeing_outcome = sum(10*(as.numeric(WBS_FullySat_T1) - 1),
                                 10*(as.numeric(WBS_Involved_T1) - 1), 
                                 10*(as.numeric(WBS_Function_T1) - 1), na.rm = TRUE )/3) |> 
  ungroup() |> 
  mutate(mentalHealth_outcome = RateMentalHealth_T1) |> 
  select(wellbeing_outcome, mentalHealth_outcome, SURVID) |> drop_na()

outcomes$wellbeing_outcome <- cut(outcomes$wellbeing_outcome, 
                                  breaks = c(0, 25, 50, 75, 100), 
                                  labels = c("0", "1", "2", "3"))
colnames(outcomes) <- c("wellbeing", "mentalhealth", "survid")

data <- inner_join(needs, outcomes, by = "survid") |> dplyr::select(-survid)

PC = bnlearn::pc.stable(data.frame(data), test = "mi-sh", alpha = 0.05)
OrdCD(data.frame(data))


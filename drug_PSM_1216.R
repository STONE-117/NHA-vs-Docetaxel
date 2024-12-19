#session info-------
sessionInfo()
library(readxl)
library(stringr)
library(survival)
library(survminer)
library(dplyr)
library(lubridate)
library(rio)
library(tidyr)
library(parallel)

file_list <- list.files(pattern = "*.xlsx", full.names = TRUE)
data_frames <- lapply(file_list, read_excel)
drug<-bind_rows(data_frames)
colnames(drug)<-c("reference_key","dispensing_date","dispensing_institution\3","dispensing_specialty\4","dispensing_institution\5",
                  "dispensing_specialty\6","prescription_start","prescription_end","multiple_dosage_indicator","drug_item_code",
                  "drug_name","route","drug_strength","dosage","dosage_unit",
                  "drug_frequency","dispensing_duration","dispensing_duration_unit","quantity","base_unit",
                  "action_status","therapeutic_classification","type_of_patient")
export(drug,"drug.RDS")
drug$drug_name<-str_to_upper(drug$drug_name)
drug<-drug[!is.na(drug$drug_name),]

#cohort identification-------
Dx_pc<-Dx[str_which(Dx$all_diagnosis_code_ICD9,"^185"),]
Dx_list<-aggregate(Dx_pc$reference_date~Dx_pc$reference_key,Dx_pc,min)
colnames(Dx_list)<-c("reference_key","diagnosis_date")
Dx_list<-Dx_list[Dx_list$reference_key%in%Dx_pc[Dx_pc$reference_date>="2001-1-1"&Dx_pc$reference_date<"2022-1-1",]$reference_key,]
Dx_list<-left_join(Dx_list,PC_dx[,1:2],by="reference_key")
Dx_list<-Dx_list[!duplicated(Dx_list$reference_key),]
Dx_list$age_at_diagnosis<- as.integer(year(Dx_list$diagnosis_date) - year(Dx_list$`Date of Birth (yyyy-mm-dd)`))
Dx_list<-Dx_list[Dx_list$age_at_diagnosis>=40,]
no_pc<-Dx[!str_detect(Dx$all_diagnosis_code_ICD9,"^185"),]
no_pc<-merge(no_pc,Dx_list,by="reference_key",all.x = F,all.y = F)
no_pc<-no_pc[no_pc$reference_date<no_pc$diagnosis_date,]
Dx_list<-Dx_list[!Dx_list$reference_key%in%no_pc[str_which(no_pc$all_diagnosis_code_ICD9,"^1[4-9][0-9]|^20[0-9]|^23[0-9]"),]$reference_key,]
drug$drug_name<-str_to_upper(drug$drug_name)
colnames(Dx_list)[3]<-"DoB"
drug_cc<-drug[drug$reference_key%in%Dx_list$reference_key,]
abiraterone_cc<-drug_cc[str_which(drug_cc$drug_name,"ABIRATERONE"),]
abiraterone_cc<-aggregate(dispensing_date~reference_key,abiraterone_cc,min)
colnames(abiraterone_cc)[2]<-"ABI_first_dispensing"
enzalutamide_cc<-drug_cc[str_which(drug_cc$drug_name,"ENZALUTAMIDE"),]
enzalutamide_cc<-aggregate(dispensing_date~reference_key,enzalutamide_cc,min)
colnames(enzalutamide_cc)[2]<-"ENZA_first_dispensing"
docetaxel_cc<-drug_cc[str_which(drug_cc$drug_name,"DOCETAXEL"),]
docetaxel_cc<-aggregate(dispensing_date~reference_key,docetaxel_cc,min)
colnames(docetaxel_cc)[2]<-"DOCE_first_dispensing"
#first exposure
adjunct_cc<-abiraterone_cc %>% 
  full_join(enzalutamide_cc,by="reference_key") %>% 
  full_join(docetaxel_cc,by="reference_key")
adjunct_cc_no_na<-adjunct_cc
adjunct_cc_no_na[is.na(adjunct_cc_no_na)]<-as.Date("2023-12-31",tz="UTC")
DOCE_cc<-adjunct_cc_no_na %>% filter(DOCE_first_dispensing<ABI_first_dispensing&DOCE_first_dispensing<ENZA_first_dispensing)
NHA_cc<-adjunct_cc_no_na %>% filter(DOCE_first_dispensing>ABI_first_dispensing|DOCE_first_dispensing>ENZA_first_dispensing)
ABI_cc<-adjunct_cc_no_na %>% filter(DOCE_first_dispensing>ABI_first_dispensing&ENZA_first_dispensing>ABI_first_dispensing)
ENZA_cc<-adjunct_cc_no_na %>% filter(DOCE_first_dispensing>ENZA_first_dispensing&ENZA_first_dispensing<ABI_first_dispensing)
docetaxel_rk<-data.frame(therapy="Docetaxel",
                         reference_key=DOCE_cc$reference_key,
                         dispensing_date=DOCE_cc$DOCE_first_dispensing,
                         switch.1=DOCE_cc$ABI_first_dispensing,
                         switch.2=DOCE_cc$ENZA_first_dispensing)
abiraterone_rk<-data.frame(therapy="Abiraterone",
                           reference_key=ABI_cc$reference_key,
                           dispensing_date=ABI_cc$ABI_first_dispensing,
                           switch.1=ABI_cc$ENZA_first_dispensing,
                           switch.2=ABI_cc$DOCE_first_dispensing)
enzalutamide_rk<-data.frame(therapy="Enzalutamide",
                            reference_key=ENZA_cc$reference_key,
                            dispensing_date=ENZA_cc$ENZA_first_dispensing,
                            switch.1=ENZA_cc$ABI_first_dispensing,
                            switch.2=ENZA_cc$DOCE_first_dispensing)
adjunct_rk<-rbind(docetaxel_rk,abiraterone_rk,enzalutamide_rk)
adjunct_rk<-merge(adjunct_rk,PC_dx[,1:2],by="reference_key",all.y = FALSE)
adjunct_rk<-adjunct_rk[!duplicated(adjunct_rk),]
colnames(adjunct_rk)[3]<-"index_date"
adjunct_rk<-merge(adjunct_rk,PC_dx[,c(1,6)],by="reference_key",all.y = F)
adjunct_rk<-adjunct_rk[!duplicated(adjunct_rk),]
adjunct_rk$death<-ifelse(is.na(adjunct_rk$date_of_registered_death),0,1)
endtime<-as.Date("2022-12-31")
adjunct_rk$index_date<-as.Date(adjunct_rk$index_date)
adjunct_rk$switch.1<-as.Date(adjunct_rk$switch.1)
adjunct_rk$switch.2<-as.Date(adjunct_rk$switch.2)
adjunct_rk$date_of_registered_death<-as.Date(adjunct_rk$date_of_registered_death)
adjunct_rk$survtime<-as.numeric(pmin(adjunct_rk$switch.1,adjunct_rk$switch.2,adjunct_rk$date_of_registered_death,endtime,na.rm = T)-adjunct_rk$index_date)
#baseline
drug_baseline<-left_join(adjunct_rk,drug,by="reference_key")
drug_baseline$index_date<-as.Date(drug_baseline$index_date)
drug_baseline$dispensing_date<-as.Date(drug_baseline$dispensing_date)
drug_baseline<-drug_baseline %>% filter(dispensing_date<index_date&index_date-dispensing_date<=365)
Dx_baseline<-left_join(adjunct_rk,Dx,by="reference_key")
Dx_baseline$index_date<-as.Date(Dx_baseline$index_date)
Dx_baseline$reference_date<-as.Date(Dx_baseline$reference_date)
Dx_baseline<-Dx_baseline %>% filter(reference_date<=index_date)
Px_baseline<-left_join(adjunct_rk,px,by="reference_key")
Px_baseline$index_date<-as.Date(Px_baseline$index_date)
Px_baseline$admission_date<-as.Date(Px_baseline$admission_date)
Px_baseline<-Px_baseline %>% filter(admission_date<=index_date)
PSA_baseline<-left_join(adjunct_rk,test_PSA,by="reference_key")
PSA_baseline$index_date<-as.Date(PSA_baseline$index_date)
PSA_baseline<-PSA_baseline %>% filter(reference_date<=index_date&index_date-reference_date<=31) %>% 
  arrange(reference_key,desc(reference_date_posix) )%>% 
  group_by(reference_key) %>% 
  filter(row_number()==1) %>% 
  ungroup()

#CR-----
tstst <- readRDS("~/Desktop/PSA_screening/tstst.RDS")
PSA_raw <- readRDS("~/Desktop/PSA_screening/PSA_raw.RDS")
ca.date<-tstst %>% 
  filter(res_tstst_fl==TRUE) %>% 
  group_by(Reference_Key) %>% 
  summarise(castration.time=min(date))
CRPC <- PSA_raw %>% left_join(ca.date,by="Reference_Key") %>% 
  group_by(Reference_Key) %>%
  filter(date>=castration.time) %>%
  arrange(Reference_Key,date) %>%
  group_by(Reference_Key) %>%
  mutate(lag_result = lag(res),
         lag_date = lag(date),
         pct_change = (res - lag_result) / lag_result * 100,
         date_diff = as.numeric(difftime(date, lag_date, units = "weeks")),
         is_increased = ifelse(pct_change >= 50 & date_diff >= 1, 1, 0)) %>%
  mutate(is_increased = replace_na(is_increased, 0)) %>%
  mutate(three_consecutive_increases = is_increased==1&lag(is_increased)==1&lag(is_increased,2)==1) %>%
  ungroup() %>% 
  filter(three_consecutive_increases & res > 2)
CRPC.date <- CRPC %>%
  group_by(Reference_Key) %>%
  summarise(third_rise_date = min(date[res > 2 & is_increased])) %>%
  ungroup()
colnames(CRPC.date)[1:2]<-c("reference_key","crpc_date")
#metastasis----
Dx.afterPC<-Dx_list %>% 
  dplyr::select(reference_key,diagnosis_date) %>%
  left_join(Dx %>% dplyr::select(reference_key,reference_date,all_diagnosis_code_ICD9),by="reference_key") %>%
  filter(diagnosis_date<reference_date)
m.date<-Dx.afterPC %>%
  filter(str_detect(all_diagnosis_code_ICD9,"^19[6-8]")) %>%
  group_by(reference_key) %>%
  summarise(metastasis_date=min(diagnosis_date)) %>% 
  mutate(metastasis_date=as.Date(metastasis_date))
adjunct_new<-adjunct_rk %>% 
  dplyr::select(reference_key,therapy,index_date) %>%
  left_join(m.date,by="reference_key") %>%
  left_join(CRPC.date,by="reference_key") %>% 
  mutate(metastasis_date=as.Date(ifelse(is.na(metastasis_date),as.Date("2023-12-31"),metastasis_date)),
         crpc_date=as.Date(ifelse(is.na(crpc_date),as.Date("2023-12-31"),crpc_date)))
mhspc.cohort<-adjunct_new %>% 
  filter(index_date>=metastasis_date&index_date<crpc_date)
table(mhspc.cohort$therapy)
nmcrpc.cohort<-adjunct_new %>% 
  filter(index_date>=crpc_date&index_date<metastasis_date)
mcrpc.cohort<-adjunct_new %>% 
  filter(index_date>=crpc_date&index_date>=metastasis_date)


#adjunct table------
adjunct<-merge(adjunct_rk[,-c(4:6)],PSA_baseline[,c(1,19)],by="reference_key",all.x = T) %>% 
  left_join(Dx_list[,c(1:3)],by="reference_key") %>% 
  mutate(DoB=as.Date(DoB),
         diagnosis_date=as.Date(diagnosis_date)) %>%
  mutate(age_at_diagnosis=as.numeric(diagnosis_date-DoB)/365.25,,
         age_at_dispensing=as.numeric(index_date-DoB)/365.25,
         disease_durtion=as.numeric(index_date-diagnosis_date)) %>% 
  mutate(dx_IHD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^41[0-3]|^414.0|^414.8|^414.9|^429.7|^V45.81|^V45.82"),]$reference_key,1,0)) %>% 
  mutate(dx_CAD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^43[0-8]"),]$reference_key,1,0)) %>% 
  mutate(dx_PVD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^44[0-3]"),]$reference_key,1,0)) %>% 
  mutate(px_CABG=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.1|^36.1[0-7]|^36.19|^36.2"),]$reference_key,1,0)) %>% 
  mutate(px_PCI=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.01|^36.02|^36.0[5-7]"),]$reference_key,1,0)) %>% 
  mutate(px_coronary=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.0[0|3|4|9]|^36.2|^36.3|^36.91|^36.99"),]$reference_key,1,0)) %>% 
  mutate(dx_HF=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^398.91|^402.01|^402.11|^402.91|^404.0[1|3]|^404.1[1|3]|^404.9[1|3]|^428|^V42.1"),]$reference_key,1,0)) %>% 
  mutate(dx_liver=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^070\\.2[2-3]|^070\\.3[2-3]|^070\\.[4-5]4|^070\\.[6|9]|^57[0-1]|^573\\.[3|4|8|9]|^V42\\.7|^456\\.[0-2]|^572\\.[2-4|8]"),]$reference_key,1,0)) %>% 
  mutate(dx_HT=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^40[1-5]"),]$reference_key,1,0)) %>% 
  mutate(dx_DM=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^249|^250|^357.2|^362.0|^366.41|^648.0"),]$reference_key,1,0)) %>% 
  mutate(dx_lipid=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^272.[0-4]"),]$reference_key,1,0)) %>% 
  mutate(dx_obesity=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^278"),]$reference_key,1,0)) %>% 
  mutate(dx_AF=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^427.31"),]$reference_key,1,0)) %>% 
  mutate(dx_arrhythmia=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^42[6-7]"),]$reference_key,1,0)) %>% 
  mutate(dx_cardiomypathy=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^425"),]$reference_key,1,0)) %>% 
  mutate(dx_COPD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^49[1|2|6]"),]$reference_key,1,0)) %>% 
  mutate(dx_asthma=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^493"),]$reference_key,1,0)) %>% 
  mutate(dx_alcohol=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^291|^303|^305.0|^357.5|^425.5|^571.[0-3]|^980.[8|9]|^V11.3"),]$reference_key,1,0)) %>% 
  mutate(dx_smoke=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^305.1|^V15.82"),]$reference_key,1,0)) %>% 
  mutate(dx_depression=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^296.2|^296.3|^300.4|^311"),]$reference_key,1,0)) %>% 
  mutate(dx_hemiplagia=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^342"),]$reference_key,1,0)) %>% 
  mutate(dx_rheumatic=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^710|^714|^720|^725"),]$reference_key,1,0)) %>% 
  mutate(dx_dementia=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^290|^294.[0-2]|^331"),]$reference_key,1,0)) %>% 
  mutate(dx_renal=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^403|^58[2|3|5|6]|^588.0|^V42.0|^V45.1|^V56"),]$reference_key,1,0)) %>% 
  mutate(dx_ulcer=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^53[1-4]"),]$reference_key,1,0)) %>% 
  mutate(dx_AIDS=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^042|^079.53"),]$reference_key,1,0)) %>% 
  mutate(Rx_GnRHa=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("Leuprolide|LEUPRORELIN|Goserelin|TRIPTORELIN|OCTREOTIDE", ignore_case = T)),]$reference_key,1,0)) %>%
  mutate(Rx_GnRHat=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("Degarelix", ignore_case = T)),]$reference_key,1,0)) %>%
  mutate(Rx_ARB=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("BICALUTAMIDE|FLUTAMIDE", ignore_case = T)),]$reference_key,1,0)) %>% 
  mutate(Rx_cv_chemo=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("adriamycin|epirubicin|Trastuzumab", ignore_case = T)),]$reference_key,1,0)) %>% 
  mutate(Rx_chemo=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^8.[1|2]"),]$reference_key,1,0)) %>% 
  mutate(Rx_chemo=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^V58.1|^V66.2|^V67.2"),]$reference_key,1,Rx_chemo)) %>% 
  mutate(Rx_chemo=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^99.25"),]$reference_key,1,Rx_chemo)) %>% 
  mutate(Rx_radiation=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^V58.0|^V66.1|^V67.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_radiation=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^92.[2-3][0-9]"),]$reference_key,1,Rx_radiation)) %>%
  mutate(px_coronary=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.0[3|4|9]|^36.0|^36.[2|3]|^36.9[1|9]"),]$reference_key,1,0)) %>% 
  mutate(Rx_h2=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^1.3.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_glycoside=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.1.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_diuretic=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.2.1|^2.2.3|^2.2.4"),]$reference_key,1,0)) %>% 
  mutate(Rx_AA=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.3.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_beta=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.4"),]$reference_key,1,0)) %>% 
  mutate(Rx_vasodilator=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_aHT=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_alpha=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.4"),]$reference_key,1,0)) %>% 
  mutate(Rx_RAS=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.5.[1|2|3]"),]$reference_key,1,0)) %>% 
  mutate(Rx_CCB=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.6.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_anticoagulant=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.8.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_antiplatelet=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.9"),]$reference_key,1,0)) %>% 
  mutate(Rx_LR=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.12"),]$reference_key,1,0)) %>% 
  mutate(Rx_NRT=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^4.10"),]$reference_key,1,0)) %>% 
  mutate(Rx_dementia=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^4.11"),]$reference_key,1,0)) %>% 
  mutate(Rx_anticdiabetic=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^6.1.1.[1|2]|^6.1.2.[1|2|3]"),]$reference_key,1,0)) %>% 
  mutate(Rx_corticosteroids=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^6.3.[1|2]"),]$reference_key,1,0)) %>% 
  mutate(Rx_NSAID=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^10.1.1"),]$reference_key,1,0)) %>% 
  mutate(px_orchidectomy=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^62.4"),]$reference_key,1,0)) %>% 
  mutate(px_prostatectomy=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^60\\.[2-6]"),]$reference_key,1,0))

#remove NA PSA in adjunct
adjunct_nona<-adjunct[!is.na(adjunct$unified_results),] %>% 
  filter(index_date>="2013-1-1")


#Cohort Characteristics--------
cohort_character<-data.frame(treatment=c("doce","enza","abi","nha"),
                             number=c(length(adjunct_nona[adjunct_nona$therapy=="Docetaxel",]$reference_key),
                                      length(adjunct_nona[adjunct_nona$therapy=="Enzalutamide",]$reference_key),
                                      length(adjunct_nona[adjunct_nona$therapy=="Abiraterone",]$reference_key),
                                      length(adjunct_nona[adjunct_nona$therapy!="Docetaxel",]$reference_key)),
                             death=c(length(adjunct_nona[adjunct_nona$therapy=="Docetaxel"&adjunct_nona$death==1,]$reference_key),
                                     length(adjunct_nona[adjunct_nona$therapy=="Enzalutamide"&adjunct_nona$death==1,]$reference_key),
                                     length(adjunct_nona[adjunct_nona$therapy=="Abiraterone"&adjunct_nona$death==1,]$reference_key),
                                     length(adjunct_nona[adjunct_nona$therapy!="Docetaxel"&adjunct_nona$death==1,]$reference_key)),
                             mean_diagnosis_age=c(mean(adjunct_nona[adjunct_nona$therapy=="Docetaxel",]$age_at_diagnosis),
                                                  mean(adjunct_nona[adjunct_nona$therapy=="Enzalutamide",]$age_at_diagnosis),
                                                  mean(adjunct_nona[adjunct_nona$therapy=="Abiraterone",]$age_at_diagnosis),
                                                  mean(adjunct_nona[adjunct_nona$therapy!="Docetaxel",]$age_at_diagnosis)),
                             diagnosis_age_SD=c(sd(adjunct_nona[adjunct_nona$therapy=="Docetaxel",]$age_at_diagnosis),
                                                sd(adjunct_nona[adjunct_nona$therapy=="Enzalutamide",]$age_at_diagnosis),
                                                sd(adjunct_nona[adjunct_nona$therapy=="Abiraterone",]$age_at_diagnosis),
                                                sd(adjunct_nona[adjunct_nona$therapy!="Docetaxel",]$age_at_diagnosis)),
                             mean_dispensing_age=c(mean(adjunct_nona[adjunct_nona$therapy=="Docetaxel",]$age_at_dispensing),
                                                  mean(adjunct_nona[adjunct_nona$therapy=="Enzalutamide",]$age_at_dispensing),
                                                  mean(adjunct_nona[adjunct_nona$therapy=="Abiraterone",]$age_at_dispensing),
                                                  mean(adjunct_nona[adjunct_nona$therapy!="Docetaxel",]$age_at_dispensing)),
                             dispensing_age_SD=c(sd(adjunct_nona[adjunct_nona$therapy=="Docetaxel",]$age_at_dispensing),
                                                sd(adjunct_nona[adjunct_nona$therapy=="Enzalutamide",]$age_at_dispensing),
                                                sd(adjunct_nona[adjunct_nona$therapy=="Abiraterone",]$age_at_dispensing),
                                                sd(adjunct_nona[adjunct_nona$therapy!="Docetaxel",]$age_at_dispensing)),
                             median_followup=c(median(adjunct_nona[adjunct_nona$therapy=="Docetaxel",]$survtime),
                                               median(adjunct_nona[adjunct_nona$therapy=="Enzalutamide",]$survtime),
                                               median(adjunct_nona[adjunct_nona$therapy=="Abiraterone",]$survtime),
                                               median(adjunct_nona[adjunct_nona$therapy!="Docetaxel",]$survtime)),
                             followup_IQR=c(IQR(adjunct_nona[adjunct_nona$therapy=="Docetaxel",]$survtime),
                                            IQR(adjunct_nona[adjunct_nona$therapy=="Enzalutamide",]$survtime),
                                            IQR(adjunct_nona[adjunct_nona$therapy=="Abiraterone",]$survtime),
                                            IQR(adjunct_nona[adjunct_nona$therapy!="Docetaxel",]$survtime))) %>% 
  mutate(mean_diagnosis_age=paste0(round(mean_diagnosis_age, 2),"(",round(diagnosis_age_SD,2),")"),
         mean_dispensing_age=paste0(round(mean_dispensing_age, 2),"(",round(dispensing_age_SD,2),")"),
         median_followup=paste0(round(median_followup/30.5, 2),"(",round(followup_IQR/30.5,2),")")) %>% 
  select(-followup_IQR,-diagnosis_age_SD,-dispensing_age_SD)
export(cohort_character,"cohort_character.xlsx")

#HRU------
##Inpatient------
inpatient$admission_date<-as.Date(inpatient$admission_date)
inpatient_pc<-adjunct_rk %>% 
  left_join(inpatient,by="reference_key") %>% 
  filter(admission_date>=index_date) %>% 
  filter(admission_date<=pmin(switch.1,switch.2,date_of_registered_death,endtime,na.rm = T))
inpatient_num<-aggregate(inpatient_pc$number_of_episode~reference_key,
                         data = inpatient_pc,sum)
colnames(inpatient_num)<-c("reference_key","IP_num")
##AE------
AE$attendance_date<-as.Date(AE$attendance_date)
AE_pc<-adjunct_rk %>% 
  left_join(AE,by="reference_key") %>% 
  filter(attendance_date>=index_date) %>% 
  filter(attendance_date<=pmin(switch.1,switch.2,date_of_registered_death,endtime,na.rm = T))
AnE_num<-aggregate(AE_pc$attendance_date~reference_key, data = AE_pc, length)
colnames(AnE_num)[2]<-"AE_num"
AnE_num$AE_cost<-AnE_num$AE_num*1230

##outpatient-----
outpatient$appointment_date<-as.Date(outpatient$appointment_date)
outpatient_pc<-adjunct_rk %>% 
  left_join(outpatient,by="reference_key") %>% 
  filter(appointment_date>=index_date) %>% 
  filter(appointment_date<=pmin(switch.1,switch.2,date_of_registered_death,endtime,na.rm = T))
sopc_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^SOP|^NUR|^FM|^AHOP|^AHTC")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(sopc_cost=1190*service_group_EIS) %>% 
  rename(sopc_n=service_group_EIS)
inj_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^INJ|^GOP Injection")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(inj_cost=100*service_group_EIS) %>% 
  rename(inj_n=service_group_EIS)
gop_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^GOP$|^GOP IPC$|^NAHC")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(gop_cost=445*service_group_EIS)%>% 
  rename(gop_n=service_group_EIS)
ah_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^AHC")|str_detect(service_type_code_EIS,"^ICMA")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(ah_cost=1730*service_group_EIS)%>% 
  rename(ah_n=service_group_EIS)
cg_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^CGAT")|str_detect(service_type_code_EIS,"^ICMD|^ICMN")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(cg_cost=535*service_group_EIS)%>% 
  rename(cg_n=service_group_EIS)
gdh_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^GDH")|str_detect(service_type_code_EIS,"^AHGD")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(gdh_cost=1960*service_group_EIS)%>% 
  rename(gdh_n=service_group_EIS)
rdh_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^Day Rehabilitation$")|str_detect(service_type_code_EIS,"^AHRD")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(rdh_cost=1320*service_group_EIS)%>% 
  rename(rdh_n=service_group_EIS)
cp_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^PG Outreach$|^Community Psychiatric Services$")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(cp_cost=1550*service_group_EIS)%>% 
  rename(cp_n=service_group_EIS)
outpatient_num<-Reduce(function(x,y)full_join(x,y,by="reference_key"),
                       list(sopc_num,inj_num,gop_num,ah_num,cg_num,gdh_num,rdh_num,cp_num))
outpatient_num[is.na(outpatient_num)]<-0
outpatient_num<-outpatient_num %>% 
  mutate(OP_num=sopc_n+inj_n+gop_n+ah_n+cg_n+gdh_n+rdh_n+cp_n) %>% 
  mutate(OP_cost=sopc_cost+inj_cost+gop_cost+ah_cost+cg_cost+gdh_cost+rdh_cost+cp_cost)
##inpatient LOS----
inpatient_los<-aggregate(inpatient_pc$LOS_acute_general~reference_key,
                         data = inpatient_pc,sum)
colnames(inpatient_los)<-c("reference_key","IP_los")
inpatient_cost<-aggregate(cbind(inpatient_pc$LOS_acute_general_acute,
                                inpatient_pc$LOS_acute_general_high_dependency,
                                inpatient_pc$LOS_acute_general_intensive_care,
                                inpatient_pc$LOS_psychiatry,
                                inpatient_pc$LOS_convalescent_rehabilitation_infirmary)~reference_key,
                          data = inpatient_pc,sum)
colnames(inpatient_cost)<-c("reference_key","LOS_actue","LOS_hd","LOS_icu","lOS_psy","LOS_conv")
inpatient_cost<-inpatient_cost %>% 
  mutate(ip_cost=5100*(LOS_actue+LOS_conv)+2340*lOS_psy+24400*LOS_icu+13650*LOS_hd)

adjunct_nona<-adjunct_nona %>% 
  left_join(inpatient_num,by="reference_key") %>% 
  left_join(AnE_num,by="reference_key") %>% 
  left_join(outpatient_num[,c(1,18,19)],by="reference_key") %>% 
  left_join(inpatient_los,by="reference_key") %>% 
  left_join(inpatient_cost[,c(1,7)],by="reference_key") %>% 
  mutate(IP_yearlynum=IP_num/survtime*365.25,
         AE_yearlynum=AE_num/survtime*365.25,
         AE_yearlycost=AE_cost/survtime*365.25,
         OP_yearlynum=OP_num/survtime*365.25,
         OP_yearlycost=OP_cost/survtime*365.25,
         IP_yearlylos=IP_los/survtime*365.25,
         IP_yearlycost=ip_cost/survtime*365.25)
adjunct_nona[,65:78]<-adjunct_nona[,65:78] %>% 
  mutate(across(everything(), ~ replace_na(., 0)))


#HRU after six months------
##Inpatient------
inpatient_pc<-adjunct_rk %>% 
  left_join(inpatient,by="reference_key") %>% 
  filter(admission_date>=index_date+180) %>% 
  filter(admission_date<=pmin(switch.1,switch.2,date_of_registered_death,endtime,na.rm = T))
inpatient_num_sup<-aggregate(inpatient_pc$number_of_episode~reference_key,
                         data = inpatient_pc,sum)
colnames(inpatient_num_sup)<-c("reference_key","IP_num_sup")
##AE------
AE_pc<-adjunct_rk %>% 
  left_join(AE,by="reference_key") %>% 
  filter(attendance_date>=index_date+180) %>% 
  filter(attendance_date<=pmin(switch.1,switch.2,date_of_registered_death,endtime,na.rm = T))
AnE_num_sup<-aggregate(AE_pc$attendance_date~reference_key, data = AE_pc, length)
colnames(AnE_num_sup)[2]<-"AE_num_sup"
AnE_num_sup$AE_cost_sup<-AnE_num_sup$AE_num*1230

##outpatient-----
outpatient_pc<-adjunct_rk %>% 
  left_join(outpatient,by="reference_key") %>% 
  filter(appointment_date>=index_date+180) %>% 
  filter(appointment_date<=pmin(switch.1,switch.2,date_of_registered_death,endtime,na.rm = T))
sopc_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^SOP|^NUR|^FM|^AHOP|^AHTC")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(sopc_cost=1190*service_group_EIS) %>% 
  rename(sopc_n=service_group_EIS)
inj_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^INJ|^GOP Injection")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(inj_cost=100*service_group_EIS) %>% 
  rename(inj_n=service_group_EIS)
gop_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^GOP$|^GOP IPC$|^NAHC")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(gop_cost=445*service_group_EIS)%>% 
  rename(gop_n=service_group_EIS)
ah_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^AHC")|str_detect(service_type_code_EIS,"^ICMA")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(ah_cost=1730*service_group_EIS)%>% 
  rename(ah_n=service_group_EIS)
cg_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^CGAT")|str_detect(service_type_code_EIS,"^ICMD|^ICMN")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(cg_cost=535*service_group_EIS)%>% 
  rename(cg_n=service_group_EIS)
gdh_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^GDH")|str_detect(service_type_code_EIS,"^AHGD")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(gdh_cost=1960*service_group_EIS)%>% 
  rename(gdh_n=service_group_EIS)
rdh_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^Day Rehabilitation$")|str_detect(service_type_code_EIS,"^AHRD")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(rdh_cost=1320*service_group_EIS)%>% 
  rename(rdh_n=service_group_EIS)
cp_num<-outpatient_pc %>% 
  filter(str_detect(service_group_EIS,"^PG Outreach$|^Community Psychiatric Services$")) %>% 
  aggregate(service_group_EIS~reference_key, length) %>% 
  mutate(cp_cost=1550*service_group_EIS)%>% 
  rename(cp_n=service_group_EIS)
outpatient_num_sup<-Reduce(function(x,y)full_join(x,y,by="reference_key"),
                       list(sopc_num,inj_num,gop_num,ah_num,cg_num,gdh_num,rdh_num,cp_num))
outpatient_num_sup[is.na(outpatient_num_sup)]<-0
outpatient_num_sup<-outpatient_num_sup %>% 
  mutate(OP_num_sup=sopc_n+inj_n+gop_n+ah_n+cg_n+gdh_n+rdh_n+cp_n) %>% 
  mutate(OP_cost_sup=sopc_cost+inj_cost+gop_cost+ah_cost+cg_cost+gdh_cost+rdh_cost+cp_cost)
##inpatient LOS----
inpatient_los_sup<-aggregate(inpatient_pc$LOS_acute_general~reference_key,
                         data = inpatient_pc,sum)
colnames(inpatient_los_sup)<-c("reference_key","IP_los_sup")
inpatient_cost_sup<-aggregate(cbind(inpatient_pc$LOS_acute_general_acute,
                                inpatient_pc$LOS_acute_general_high_dependency,
                                inpatient_pc$LOS_acute_general_intensive_care,
                                inpatient_pc$LOS_psychiatry,
                                inpatient_pc$LOS_convalescent_rehabilitation_infirmary)~reference_key,
                          data = inpatient_pc,sum)
colnames(inpatient_cost_sup)<-c("reference_key","LOS_actue","LOS_hd","LOS_icu","lOS_psy","LOS_conv")
inpatient_cost_sup<-inpatient_cost_sup %>% 
  mutate(ip_cost_sup=5100*(LOS_actue+LOS_conv)+2340*lOS_psy+24400*LOS_icu+13650*LOS_hd)

adjunct_nona<-adjunct_nona %>% 
  left_join(inpatient_num_sup,by="reference_key") %>% 
  left_join(AnE_num_sup,by="reference_key") %>% 
  left_join(outpatient_num_sup[,c(1,18,19)],by="reference_key") %>% 
  left_join(inpatient_los_sup,by="reference_key") %>% 
  left_join(inpatient_cost_sup[,c(1,7)],by="reference_key") %>% 
  mutate(IP_yearlynum_sup=IP_num_sup/survtime*365.25,
         AE_yearlynum_sup=AE_num_sup/survtime*365.25,
         AE_yearlycost_sup=AE_cost_sup/survtime*365.25,
         OP_yearlynum_sup=OP_num_sup/survtime*365.25,
         OP_yearlycost_sup=OP_cost_sup/survtime*365.25,
         IP_yearlylos_sup=IP_los_sup/survtime*365.25,
         IP_yearlycost_sup=ip_cost_sup/survtime*365.25)
adjunct_nona[,79:92]<-adjunct_nona[,79:92] %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

#PSM------
library(MatchIt)
library(cobalt)
library(scales)
library(rstpm2)

##PS matching NHA&T------
adjunct_NHA<-adjunct_nona %>% 
  mutate(treat=ifelse(therapy=="Docetaxel",0,1)) 
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_NHA,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
love.plot(m.out, threshold = .1, which = "both")
#baseline table
#before match
smd <- function(treated, untreated) {
  num <- mean(treated) - mean(untreated)
  denom <- sqrt((var(treated) + var(untreated)) / 2)
  return(num / denom)
}
summary_table_before <- data.frame(
  covariate = character(),
  treat_group = numeric(),
  treat_sd.per = numeric(),
  untreat_group = numeric(),
  untreat_sd.per = numeric(),
  SMD = numeric(),
  stringsAsFactors = FALSE
)
for (i in c(11,7)) {
  treat_group_value <- mean(as.numeric(adjunct_NHA[adjunct_NHA$treat == 1,][,i]))
  treat_group_sd<- sd(as.numeric(adjunct_NHA[adjunct_NHA$treat == 1,][,i]))
  untreat_group_value <- mean(as.numeric(adjunct_NHA[adjunct_NHA$treat == 0,][,i]))
  untreat_group_sd<- sd(as.numeric(adjunct_NHA[adjunct_NHA$treat == 0,][,i]))
  smd_value <- smd(adjunct_NHA[adjunct_NHA$treat == 1,][,i], adjunct_NHA[adjunct_NHA$treat == 0,][,i])
  result_table<-data.frame(
    covariate = colnames(adjunct_NHA)[i],
    treat_group = round(as.numeric(treat_group_value),digits = 3),
    treat_sd.per = round(as.numeric(treat_group_sd),digits = 3),
    untreat_group = round(as.numeric(untreat_group_value),digits = 3),
    untreat_sd.per = round(as.numeric(untreat_group_sd),digits = 3),
    SMD = round(as.numeric(smd_value),digit = 3))
  summary_table_before<- rbind(summary_table_before,result_table)
}
for (i in c(64:13)) {
  treat_group_value <- sum(adjunct_NHA[adjunct_NHA$treat == 1,][,i])
  treat_group_per <- treat_group_value/length(adjunct_NHA[adjunct_NHA$treat == 1,][,i])
  untreat_group_value <- sum(adjunct_NHA[adjunct_NHA$treat == 0,][,i])
  untreat_group_per <- untreat_group_value/length(adjunct_NHA[adjunct_NHA$treat == 0,][,i])
  smd_value <- smd(adjunct_NHA[adjunct_NHA$treat == 1,][,i], adjunct_NHA[adjunct_NHA$treat == 0,][,i])
  result_table<-data.frame(
    covariate = colnames(adjunct_NHA)[i],
    treat_group = treat_group_value,
    treat_sd.per = round(as.numeric(treat_group_per),digits = 3)*100,
    untreat_group = untreat_group_value,
    untreat_sd.per = round(as.numeric(untreat_group_per),digits = 3)*100,
    SMD = round(as.numeric(smd_value),digits = 3))
  summary_table_before<- rbind(summary_table_before,result_table)
}
summary_table_before<-summary_table_before %>% 
  mutate(treat.combine=paste(treat_group,"(",
                             round(treat_sd.per,digits = 1),"%)"),
         untreat.combine=paste(untreat_group,"(",
                               round(untreat_sd.per,digits = 1),"%)")) %>% 
  select(covariate,untreat.combine,treat.combine,SMD)

#after match
summary_table_after <- data.frame(
  covariate = character(),
  treat_group = numeric(),
  treat_sd.per = numeric(),
  untreat_group = numeric(),
  untreat_sd.per = numeric(),
  SMD = numeric(),
  stringsAsFactors = FALSE
)
for (i in c(11,7)) {
  treat_group_value <- mean(as.numeric(matched_data[matched_data$treat == 1,][,i]))
  treat_group_sd<- sd(as.numeric(matched_data[matched_data$treat == 1,][,i]))
  untreat_group_value <- mean(as.numeric(matched_data[matched_data$treat == 0,][,i]))
  untreat_group_sd<- sd(as.numeric(matched_data[matched_data$treat == 0,][,i]))
  smd_value <- smd(matched_data[matched_data$treat == 1,][,i], matched_data[matched_data$treat == 0,][,i])
  result_table<-data.frame(
    covariate = colnames(matched_data)[i],
    treat_group = round(as.numeric(treat_group_value),digits = 3),
    treat_sd.per = round(as.numeric(treat_group_sd),digits = 3),
    untreat_group = round(as.numeric(untreat_group_value),digits = 3),
    untreat_sd.per = round(as.numeric(untreat_group_sd),digits = 3),
    SMD = round(as.numeric(smd_value),digit = 3))
  summary_table_after<- rbind(summary_table_after,result_table)
}
for (i in c(64:13)) {
  treat_group_value <- sum(matched_data[matched_data$treat == 1,][,i])
  treat_group_per <- treat_group_value/length(matched_data[matched_data$treat == 1,][,i])
  untreat_group_value <- sum(matched_data[matched_data$treat == 0,][,i])
  untreat_group_per <- untreat_group_value/length(matched_data[matched_data$treat == 0,][,i])
  smd_value <- smd(matched_data[matched_data$treat == 1,][,i], matched_data[matched_data$treat == 0,][,i])
  result_table<-data.frame(
    covariate = colnames(matched_data)[i],
    treat_group = treat_group_value,
    treat_sd.per = round(as.numeric(treat_group_per),digits = 3)*100,
    untreat_group = untreat_group_value,
    untreat_sd.per = round(as.numeric(untreat_group_per),digits = 3)*100,
    SMD = round(as.numeric(smd_value),digits = 3))
  summary_table_after<- rbind(summary_table_after,result_table)
}
summary_table_after<-summary_table_after %>% 
  mutate(treat.combine=paste(treat_group,"(",
                             round(treat_sd.per,digits = 1),"%)"),
         untreat.combine=paste(untreat_group,"(",
                               round(untreat_sd.per,digits = 1),"%)")) %>% 
  select(covariate,untreat.combine,treat.combine,SMD)

export(summary_table_before,"baseline_table_before.xlsx")
export(summary_table_after,"baseline_table_after.xlsx")

#survival analysis
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("OS_nha.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)", ylab="Overall survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","NHA","KM Docetaxel","KM NHA"),lty=1:2,lwd=3, col=c("black","black","steelblue","steelblue"))
title(main="NHA vs Docetaxel")
text(1600,0.5,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
dev.off()
#HR
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#"0.998 (95%CI 0.891 to 1.114)"

fit<-survfit(Surv(survtime, death)~1,data=matched_data)
fit #median 608 days
#5year survival probability
summary_fit <- summary(fit, times = 1826.25)
survival_at_max_time <- summary_fit$surv[summary_fit$time == 1826.25]
surv_ci_lower <- summary_fit$lower[summary_fit$time == 365.25*5]
surv_ci_upper <- summary_fit$upper[summary_fit$time == 365.25*5]
print(paste("5-year Survival probability:", round(survival_at_max_time,2),"(95%CI",round(surv_ci_lower,2),"to",round(surv_ci_upper,2),")"))
#5-year Survival probability: 0.18 (95%CI 0.14 to 0.23 )

# Summarize the fit
max(matched_data$survtime) #3533
summary_fit <- summary(fit, times = max(matched_data$survtime))  # Use the maximum time in your data as a reference
# The survival probability at the end of the study
survival_at_max_time <- summary_fit$surv
surv_ci_lower <- summary_fit$lower[summary_fit$time == max(matched_data$survtime)]
surv_ci_upper <- summary_fit$upper[summary_fit$time == max(matched_data$survtime)]
print(paste("Survival probability at end:", round(survival_at_max_time,2),"(95%CI",round(surv_ci_lower,2),"to",round(surv_ci_upper,2),")"))
#Survival probability at end: 0.16 (95%CI 0.12 to 0.21 )


matched_data$treat<-as.factor(matched_data$treat)
treatment_labels<-c("0"="Docetaxel","1"="NHA")
matched_data_NHA<-matched_data%>% 
  select(reference_key,IP_yearlynum,AE_yearlynum,OP_yearlynum,IP_yearlylos,treat) %>% 
  mutate(source="Overall")
#IP
wilcox.test(matched_data$IP_yearlynum[matched_data$treat==1],matched_data$IP_yearlynum[matched_data$treat==0])

#IP LOS
wilcox.test(matched_data$IP_yearlylos[matched_data$treat==1],matched_data$IP_yearlylos[matched_data$treat==0])

#AE
wilcox.test(matched_data$AE_yearlynum[matched_data$treat==1],matched_data$AE_yearlynum[matched_data$treat==0])

#OP
wilcox.test(matched_data$OP_yearlynum[matched_data$treat==1],matched_data$OP_yearlynum[matched_data$treat==0],conf.int = T)


matched_data_NHA_sup<-matched_data%>% 
  select(reference_key,IP_yearlynum_sup,AE_yearlynum_sup,OP_yearlynum_sup,IP_yearlylos_sup,treat) %>% 
  mutate(source="Overall")
#IP
wilcox.test(matched_data$IP_yearlynum[matched_data$treat==1],matched_data$IP_yearlynum[matched_data$treat==0])
median(matched_data$IP_yearlynum[matched_data$treat==1])
median(matched_data$IP_yearlynum[matched_data$treat==0])
#IP LOS
wilcox.test(matched_data$IP_yearlylos[matched_data$treat==1],matched_data$IP_yearlylos[matched_data$treat==0])

#AE
wilcox.test(matched_data$AE_yearlynum[matched_data$treat==1],matched_data$AE_yearlynum[matched_data$treat==0])

#OP
wilcox.test(matched_data$OP_yearlynum[matched_data$treat==1],matched_data$OP_yearlynum[matched_data$treat==0],conf.int = T)


##Mean cost difference-----

difference<-0
for (i in 1:length(unique(matched_data$subclass))) {
  treat<-filter(matched_data,subclass==i&treat==1)
  untreat<-filter(matched_data,subclass==i&treat==0)
  gap<-treat$AE_yearlycost+treat$OP_yearlycost+treat$IP_yearlycost-untreat$AE_yearlycost-untreat$OP_yearlycost-untreat$IP_yearlycost
  difference<-difference+gap
}
mean_cost_difference<-difference/length(unique(matched_data$subclass))
#-18292.11

##PS matching A&T------
adjunct_AT<-rbind(adjunct_nona[adjunct_nona$therapy=="Abiraterone",],
                  adjunct_nona[adjunct_nona$therapy=="Docetaxel",])
adjunct_AT$treat<-ifelse(adjunct_AT$therapy=="Abiraterone",1,0)
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_AT,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
love.plot(m.out, threshold = .1, which = "both")
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("OS_abi.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)", ylab="Overall survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","Abiraterone","KM Docetaxel","KM Abiraterone"),lty=1:2,lwd = 3, col=c("black","black","steelblue","steelblue"))
title(main="Abiraterone vs Docetaxel")
text(1800,0.5,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))

dev.off()
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#0.996 (95%CI 0.874 to 1.129)
matched_data$treat<-as.factor(matched_data$treat)
treatment_labels<-c("0"="Docetaxel","1"="Abiraterone")
matched_data_ABI<-matched_data%>% 
  select(reference_key,IP_yearlynum,AE_yearlynum,OP_yearlynum,IP_yearlylos,treat) %>% 
  mutate(source="Abiraterone subgroup")
#IP
wilcox.test(matched_data$IP_yearlynum[matched_data$treat==1],matched_data$IP_yearlynum[matched_data$treat==0])

#IP LOS
wilcox.test(matched_data$IP_yearlylos[matched_data$treat==1],matched_data$IP_yearlylos[matched_data$treat==0])

#AE
wilcox.test(matched_data$AE_yearlynum[matched_data$treat==1],matched_data$AE_yearlynum[matched_data$treat==0])
#0.002
#OP
wilcox.test(matched_data$OP_yearlynum[matched_data$treat==1],matched_data$OP_yearlynum[matched_data$treat==0])



##PS matching E&T------
adjunct_ET<-rbind(adjunct_nona[adjunct_nona$therapy=="Enzalutamide",],
                  adjunct_nona[adjunct_nona$therapy=="Docetaxel",])
adjunct_ET$treat<-ifelse(adjunct_ET$therapy=="Enzalutamide",1,0)
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_ET,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
love.plot(m.out, threshold = .1, which = "both")
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("OS_enz.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)", ylab="Overall survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","Enzalutamide","KM Docetaxel","KM Enzalutamide"),lty=1:2, lwd=4,col=c("black","black","steelblue","steelblue"))
title(main="Enzalutamide vs Docetaxel")
text(1200,0.55,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
dev.off()
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#0.925 (95%CI 0.782 to 1.084)

matched_data$treat<-as.factor(matched_data$treat)
treatment_labels<-c("0"="Docetaxel","1"="Enzalutamide")
matched_data_ENZ<-matched_data%>% 
  select(reference_key,IP_yearlynum,AE_yearlynum,OP_yearlynum,IP_yearlylos,treat) %>% 
  mutate(source="Enzalutamide subgroup")
#IP
wilcox.test(matched_data$IP_yearlynum[matched_data$treat==1],matched_data$IP_yearlynum[matched_data$treat==0])

#IP LOS
wilcox.test(matched_data$IP_yearlylos[matched_data$treat==1],matched_data$IP_yearlylos[matched_data$treat==0])

#AE
wilcox.test(matched_data$AE_yearlynum[matched_data$treat==1],matched_data$AE_yearlynum[matched_data$treat==0])

#OP
wilcox.test(matched_data$OP_yearlynum[matched_data$treat==1],matched_data$OP_yearlynum[matched_data$treat==0])




#dispensing age>=75-----
adjunct_NHA_75<-adjunct_NHA[adjunct_NHA$age_at_dispensing>=75,]
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_NHA_75,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
love.plot(m.out, threshold = .1, which = "both")
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("OS_old.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)", ylab="Overall survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","NHA","KM Docetaxel","KM NHA"),lty=1:2, lwd = 3 ,col=c("black","black","steelblue","steelblue"))
text(1300,0.5,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
title(main="NHA vs Docetaxel (age>=75)")
dev.off()
#HR
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#0.780 (95%CI 0.585 to 1.014)

matched_data$treat<-as.factor(matched_data$treat)
treatment_labels<-c("0"="Docetaxel","1"="NHA")
matched_data_75<-matched_data%>% 
  select(reference_key,IP_yearlynum,AE_yearlynum,OP_yearlynum,IP_yearlylos,treat) %>% 
  mutate(source="Age>=75 subgroup")
#IP
wilcox.test(matched_data$IP_yearlynum[matched_data$treat==1],matched_data$IP_yearlynum[matched_data$treat==0])

#IP LOS
wilcox.test(matched_data$IP_yearlylos[matched_data$treat==1],matched_data$IP_yearlylos[matched_data$treat==0])

#AE
wilcox.test(matched_data$AE_yearlynum[matched_data$treat==1],matched_data$AE_yearlynum[matched_data$treat==0])
#*
#OP
wilcox.test(matched_data$OP_yearlynum[matched_data$treat==1],matched_data$OP_yearlynum[matched_data$treat==0])

#dispensing age<75-----
adjunct_NHA_74<-adjunct_NHA[adjunct_NHA$age_at_dispensing<75,]
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_NHA_74,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
love.plot(m.out, threshold = .1, which = "both")
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("OS_young.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)", ylab="Overall survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","NHA","KM Docetaxel","KM NHA"),lty=1:2, lwd = 3,col=c("black","black","steelblue","steelblue"))
text(1800,0.5,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
title(main="NHA vs Docetaxel (age<75)")
dev.off()
#HR
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#0.936 (95%CI 0.820 to 1.064)

matched_data$treat<-as.factor(matched_data$treat)
treatment_labels<-c("0"="Docetaxel","1"="NHA")
matched_data_74<-matched_data %>% 
  select(reference_key,IP_yearlynum,AE_yearlynum,OP_yearlynum,IP_yearlylos,treat) %>% 
  mutate(source="Age<75 subgroup")
#IP
wilcox.test(matched_data$IP_yearlynum[matched_data$treat==1],matched_data$IP_yearlynum[matched_data$treat==0])

#IP LOS
wilcox.test(matched_data$IP_yearlylos[matched_data$treat==1],matched_data$IP_yearlylos[matched_data$treat==0])

#AE
wilcox.test(matched_data$AE_yearlynum[matched_data$treat==1],matched_data$AE_yearlynum[matched_data$treat==0])

#OP
wilcox.test(matched_data$OP_yearlynum[matched_data$treat==1],matched_data$OP_yearlynum[matched_data$treat==0])


#HRU merge plot-----
treatment_labels<-c("0"="Docetaxel","1"="NHA")
match_combine<-bind_rows(matched_data_NHA,matched_data_ABI,matched_data_ENZ,matched_data_75,matched_data_74)
match_combine$source<-factor(match_combine$source,levels=unique(match_combine$source))
#IP
sign_labels<-data.frame(source=unique(match_combine$source),
                        x=1.5,
                        y=100,
                        labels=c("**"))
ggplot(match_combine,aes(x=treat,y=IP_yearlynum))+
  geom_boxplot(aes(fill=treat))+
  facet_wrap(~ source, ncol=1) +
  coord_flip()+
  scale_fill_manual(name = "Treatment Group", 
                    values = c("0" = "blue", "1" = "red"), 
                    labels = treatment_labels) +
  ylim(c(NA,100))+
  labs(x="Treatment", y="Number of Admission (per year)", title="Inpatient Yearly Admission") +
  theme_cleveland()+
  theme(
    axis.text.y  = element_blank(),  # Remove x-axis text
    axis.ticks.y = element_blank()  # Optionally remove x-axis ticks as well
  )+
  geom_text(data = sign_labels, aes(x = x, y = y, label = labels), vjust = 0.5,label.size=0.25)

#IP los
sign_labels<-data.frame(source=unique(match_combine$source),
                        x=1.5,
                        y=300,
                        labels=c("**","**","**","**","**"))
ggplot(match_combine,aes(x=treat,y=IP_yearlylos))+
  geom_boxplot(aes(fill=treat))+
  facet_wrap(~ source, ncol=1) +
  coord_flip()+
  scale_fill_manual(name = "Treatment Group", 
                    values = c("0" = "blue", "1" = "red"), 
                    labels = treatment_labels) +
  ylim(c(NA,300))+
  labs(x="Treatment", y="LOS (days)", title="Inpatient Yearly Length of Stay") +
  theme_cleveland()+
  theme(
    axis.text.y  = element_blank(),  # Remove x-axis text
    axis.ticks.y = element_blank()  # Optionally remove x-axis ticks as well
  )+
  geom_text(data = sign_labels, aes(x = x, y = y, label = labels), vjust = 0.7,label.size=0.25)
#AE
sign_labels<-data.frame(source=unique(match_combine$source),
                        x=1.5,
                        y=75,
                        labels=c("**","**","**","*","**"))
ggplot(match_combine,aes(x=treat,y=AE_yearlynum))+
  geom_boxplot(aes(fill=treat))+
  facet_wrap(~ source, ncol=1) +
  coord_flip()+
  scale_fill_manual(name = "Treatment Group", 
                    values = c("0" = "blue", "1" = "red"), 
                    labels = treatment_labels) +
  ylim(c(NA,75))+
  labs(x="Treatment", y="Number of Admission (per year)", title="A&E Yearly Admission") +
  theme_cleveland()+
  theme(
    axis.text.y  = element_blank(),  # Remove x-axis text
    axis.ticks.y = element_blank()  # Optionally remove x-axis ticks as well
  )+
  geom_text(data = sign_labels, aes(x = x, y = y, label = labels), vjust = 0.5,label.size=0.25)
#OP
sign_labels<-data.frame(source=unique(match_combine$source),
                        x=1.5,
                        y=150,
                        labels=c("**"))
ggplot(match_combine,aes(x=treat,y=OP_yearlynum))+
  geom_boxplot(aes(fill=treat))+
  facet_wrap(~ source, ncol=1) +
  coord_flip()+
  scale_fill_manual(name = "Treatment Group", 
                    values = c("0" = "blue", "1" = "red"), 
                    labels = treatment_labels) +
  ylim(c(NA,150))+
  labs(x="Treatment", y="Number of Admission (per year)", title="Outpatient Yearly Admission") +
  theme_cleveland()+
  theme(
    axis.text.y  = element_blank(),  # Remove x-axis text
    axis.ticks.y = element_blank()  # Optionally remove x-axis ticks as well
  )+
  geom_text(data = sign_labels, aes(x = x, y = y, label = labels), vjust = 0.5,label.size=0.25)

#overall median
overall_nha_ip <- sprintf("%.3f (%.3f to %.3f)", 
                       median(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$IP_yearlynum), 
                       quantile(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                       quantile(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ip <- sprintf("%.3f (%.3f to %.3f)", 
                        median(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$IP_yearlynum), 
                        quantile(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                        quantile(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_op <- sprintf("%.3f (%.3f to %.3f)", 
                       median(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$OP_yearlynum), 
                       quantile(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                       quantile(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_op <- sprintf("%.3f (%.3f to %.3f)", 
                        median(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$OP_yearlynum), 
                        quantile(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                        quantile(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$IP_yearlylos), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_doce_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$IP_yearlylos), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_nha_ae <- sprintf("%.3f (%.3f to %.3f)", 
                             median(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$AE_yearlynum), 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Overall",]$AE_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ae <- sprintf("%.3f (%.3f to %.3f)", 
                              median(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$AE_yearlynum), 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Overall",]$AE_yearlynum, prob=c(0.25,0.75))[2])
#abiraterone median
overall_nha_ip <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$IP_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ip <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$IP_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_op <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$OP_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_op <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$OP_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                             median(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$IP_yearlylos), 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_doce_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                              median(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$IP_yearlylos), 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_nha_ae <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$AE_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Abiraterone subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ae <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$AE_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Abiraterone subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[2])
#enzaluatmide median
overall_nha_ip <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$IP_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ip <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$IP_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_op <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$OP_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_op <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$OP_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                             median(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$IP_yearlylos), 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_doce_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                              median(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$IP_yearlylos), 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_nha_ae <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$AE_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Enzalutamide subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ae <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$AE_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Enzalutamide subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[2])
#Age>=75 subgroup median
overall_nha_ip <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$IP_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ip <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$IP_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_op <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$OP_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_op <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$OP_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                             median(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$IP_yearlylos), 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_doce_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                              median(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$IP_yearlylos), 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_nha_ae <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$AE_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age>=75 subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ae <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$AE_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age>=75 subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[2])
#Age<75 subgroup
overall_nha_ip <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$IP_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ip <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$IP_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$IP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_op <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$OP_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_op <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$OP_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$OP_yearlynum, prob=c(0.25,0.75))[2])
overall_nha_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                             median(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$IP_yearlylos), 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                             quantile(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_doce_iplos <- sprintf("%.3f (%.3f to %.3f)", 
                              median(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$IP_yearlylos), 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[1], 
                              quantile(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$IP_yearlylos, prob=c(0.25,0.75))[2])
overall_nha_ae <- sprintf("%.3f (%.3f to %.3f)", 
                          median(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$AE_yearlynum), 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                          quantile(match_combine[match_combine$treat==1&match_combine$source=="Age<75 subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[2])
overall_doce_ae <- sprintf("%.3f (%.3f to %.3f)", 
                           median(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$AE_yearlynum), 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[1], 
                           quantile(match_combine[match_combine$treat==0&match_combine$source=="Age<75 subgroup",]$AE_yearlynum, prob=c(0.25,0.75))[2])
print(overall_nha_ip) 
print(overall_doce_ip)
print(overall_nha_op)
print(overall_doce_op)
print(overall_nha_iplos)
print(overall_doce_iplos)
print(overall_nha_ae)
print(overall_doce_ae)



#Sensitive analysis------
##Remove survtime<30----
data<-list(main=adjunct_NHA,abi=adjunct_AT,enz=adjunct_ET,old=adjunct_NHA_75,young=adjunct_NHA_74)
shortremove_results<-data.frame(group=NA,rm_30=NA)
for (i in names(data)) {
  adjunct_NHA_30<-data[[i]]%>% 
    filter(survtime>30)
  m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
                 data = adjunct_NHA_30,method = "nearest",caliper = 0.1)
  matched_data <- match.data(m.out)
  model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
  nha_30<-sprintf("%.3f (95%%CI %.3f to %.3f)",
                  exp(coef(model)["treat"]),
                  exp(confint(model)["treat", 1]),
                  exp(confint(model)["treat", 2]))
  shortremove_results<-rbind(shortremove_results,
                      data.frame(group=i,
                                 rm_30=sprintf("%.3f (95%%CI %.3f to %.3f)",
                                               exp(coef(model)["treat"]),
                                               exp(confint(model)["treat", 1]),
                                               exp(confint(model)["treat", 2]))))
}
shortremove_results<-shortremove_results[-1,]
##IPTW------
library(ipw)
iptw_results<-data.frame(group=NA,IPTW=NA)
for (i in names(data)) {
  adjunct_name<-data[[i]]
  colnames(adjunct_name)[11:64]<-paste0("cov_", 1:54)
  w1<-ipwpoint(
    exposure = treat,
    family = "binomial",
    link = "logit",
    numerator = ~1,
    denominator = ~cov_1+cov_3+cov_4+cov_5+cov_6+cov_7+cov_8+cov_9+cov_10+cov_11+cov_12+cov_13+cov_14+cov_15+cov_16+cov_17+cov_18+cov_19+cov_20+cov_21+cov_22+cov_23+cov_24+cov_25+cov_26+cov_27+cov_28+cov_29+cov_30+cov_31+cov_32+cov_33+cov_34+cov_35+cov_36+cov_37+cov_38+cov_39+cov_40+cov_41+cov_42+cov_43+cov_44+cov_45+cov_46+cov_47+cov_48+cov_49+cov_50+cov_51+cov_52+cov_53+cov_54,
    data = adjunct_name
  )
  adjunct_name$w1<-w1$ipw.weights
  survival_model<-coxph(Surv(survtime,death)~treat,data=adjunct_name,weights=w1)
  iptw_results<-rbind(iptw_results,
                      data.frame(group=i,
                                 IPTW=sprintf("%.3f (95%%CI %.3f to %.3f)",
                                              exp(coef(survival_model)["treat"]),
                                              exp(confint(survival_model)["treat", 1]),
                                              exp(confint(survival_model)["treat", 2]))))
}
iptw_results<-iptw_results[-1,]


##IPCW------
ipcw_results<-data.frame(group=NA,IPCW=NA)
for (i in names(data)) {
  m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
                 data = data[[i]],method = "nearest",caliper = 0.1)
  matched_data <- match.data(m.out)
  adjunct_IPCW<-adjunct_rk %>% 
    mutate(censor=ifelse(date_of_registered_death==pmin(switch.1,switch.2,date_of_registered_death,endtime,na.rm = T),0,1)) %>% 
    mutate(censor=ifelse(is.na(date_of_registered_death),1,censor)) %>% 
    right_join(matched_data[,c(1,7:79)],by="reference_key")
  c0<-coxph(Surv(survtime,censor)~1,data = adjunct_IPCW)
  c0fit<-summary(survfit(c0),times = adjunct_IPCW$survtime)
  adjunct_IPCW$k0ti<-c0fit$surv
  cz<-coxph(Surv(survtime,censor)~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_AF+dx_arrhythmia+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_renal+dx_ulcer+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
            data=adjunct_IPCW)
  process_row <- function(i) {
    datai <-adjunct_IPCW[i, ]
    czfit <- summary(survfit(cz, newdata = datai), times = datai$survtime)
    return(czfit$surv)
  }
  results <- mclapply(1:nrow(adjunct_IPCW), process_row, mc.cores = 10)
  adjunct_IPCW$kzti <- unlist(results)
  adjunct_IPCW<-adjunct_IPCW %>% 
    mutate(ipcw=1/kzti,
           ipcw_s=k0ti/kzti) %>% 
    mutate(ipcw_s_upper=quantile(ipcw_s, 0.95),
           ipcw_s_lower=quantile(ipcw_s, 0.05)) %>%
    mutate(ipcw_s=ifelse(ipcw_s>=ipcw_s_upper,ipcw_s_upper,ipcw_s),
           ipcw_s=ifelse(ipcw_s<=ipcw_s_lower,ipcw_s_lower,ipcw_s))
  survival_model<-coxph(Surv(survtime,death)~treat,data=adjunct_IPCW,weights=ipcw_s)
  ipcw_results<-rbind(ipcw_results,
                      data.frame(group=i,
                                 IPCW=sprintf("%.3f (95%%CI %.3f to %.3f)",
                                              exp(coef(survival_model)["treat"]),
                                              exp(confint(survival_model)["treat", 1]),
                                              exp(confint(survival_model)["treat", 2]))))
  
}
ipcw_results<-ipcw_results[-1,]
sensitive_results<-shortremove_results %>% 
  left_join(iptw_results,by="group") %>% 
  left_join(ipcw_results,by="group")
export(sensitive_results,"sensitive_results.xlsx")
##mhspc------
mhspc.table<- mhspc.cohort %>% 
  left_join(adjunct_rk %>% select(reference_key,date_of_registered_death,death,survtime),by="reference_key") %>%
  left_join(PSA_baseline[,c(1,19)],by="reference_key") %>% 
  left_join(Dx_list[,c(1:3)],by="reference_key") %>% 
  mutate(DoB=as.Date(DoB),
         diagnosis_date=as.Date(diagnosis_date)) %>%
  mutate(age_at_diagnosis=as.numeric(diagnosis_date-DoB)/365.25,,
         age_at_dispensing=as.numeric(index_date-DoB)/365.25,
         disease_durtion=as.numeric(index_date-diagnosis_date)) %>% 
  mutate(dx_IHD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^41[0-3]|^414.0|^414.8|^414.9|^429.7|^V45.81|^V45.82"),]$reference_key,1,0)) %>% 
  mutate(dx_CAD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^43[0-8]"),]$reference_key,1,0)) %>% 
  mutate(dx_PVD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^44[0-3]"),]$reference_key,1,0)) %>% 
  mutate(px_CABG=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.1|^36.1[0-7]|^36.19|^36.2"),]$reference_key,1,0)) %>% 
  mutate(px_PCI=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.01|^36.02|^36.0[5-7]"),]$reference_key,1,0)) %>% 
  mutate(px_coronary=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.0[0|3|4|9]|^36.2|^36.3|^36.91|^36.99"),]$reference_key,1,0)) %>% 
  mutate(dx_HF=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^398.91|^402.01|^402.11|^402.91|^404.0[1|3]|^404.1[1|3]|^404.9[1|3]|^428|^V42.1"),]$reference_key,1,0)) %>% 
  mutate(dx_liver=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^070\\.2[2-3]|^070\\.3[2-3]|^070\\.[4-5]4|^070\\.[6|9]|^57[0-1]|^573\\.[3|4|8|9]|^V42\\.7|^456\\.[0-2]|^572\\.[2-4|8]"),]$reference_key,1,0)) %>% 
  mutate(dx_HT=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^40[1-5]"),]$reference_key,1,0)) %>% 
  mutate(dx_DM=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^249|^250|^357.2|^362.0|^366.41|^648.0"),]$reference_key,1,0)) %>% 
  mutate(dx_lipid=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^272.[0-4]"),]$reference_key,1,0)) %>% 
  mutate(dx_obesity=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^278"),]$reference_key,1,0)) %>% 
  mutate(dx_AF=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^427.31"),]$reference_key,1,0)) %>% 
  mutate(dx_arrhythmia=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^42[6-7]"),]$reference_key,1,0)) %>% 
  mutate(dx_cardiomypathy=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^425"),]$reference_key,1,0)) %>% 
  mutate(dx_COPD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^49[1|2|6]"),]$reference_key,1,0)) %>% 
  mutate(dx_asthma=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^493"),]$reference_key,1,0)) %>% 
  mutate(dx_alcohol=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^291|^303|^305.0|^357.5|^425.5|^571.[0-3]|^980.[8|9]|^V11.3"),]$reference_key,1,0)) %>% 
  mutate(dx_smoke=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^305.1|^V15.82"),]$reference_key,1,0)) %>% 
  mutate(dx_depression=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^296.2|^296.3|^300.4|^311"),]$reference_key,1,0)) %>% 
  mutate(dx_hemiplagia=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^342"),]$reference_key,1,0)) %>% 
  mutate(dx_rheumatic=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^710|^714|^720|^725"),]$reference_key,1,0)) %>% 
  mutate(dx_dementia=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^290|^294.[0-2]|^331"),]$reference_key,1,0)) %>% 
  mutate(dx_renal=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^403|^58[2|3|5|6]|^588.0|^V42.0|^V45.1|^V56"),]$reference_key,1,0)) %>% 
  mutate(dx_ulcer=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^53[1-4]"),]$reference_key,1,0)) %>% 
  mutate(dx_AIDS=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^042|^079.53"),]$reference_key,1,0)) %>% 
  mutate(Rx_GnRHa=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("Leuprolide|LEUPRORELIN|Goserelin|TRIPTORELIN|OCTREOTIDE", ignore_case = T)),]$reference_key,1,0)) %>%
  mutate(Rx_GnRHat=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("Degarelix", ignore_case = T)),]$reference_key,1,0)) %>%
  mutate(Rx_ARB=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("BICALUTAMIDE|FLUTAMIDE", ignore_case = T)),]$reference_key,1,0)) %>% 
  mutate(Rx_cv_chemo=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("adriamycin|epirubicin|Trastuzumab", ignore_case = T)),]$reference_key,1,0)) %>% 
  mutate(Rx_chemo=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^8.[1|2]"),]$reference_key,1,0)) %>% 
  mutate(Rx_chemo=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^V58.1|^V66.2|^V67.2"),]$reference_key,1,Rx_chemo)) %>% 
  mutate(Rx_chemo=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^99.25"),]$reference_key,1,Rx_chemo)) %>% 
  mutate(Rx_radiation=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^V58.0|^V66.1|^V67.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_radiation=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^92.[2-3][0-9]"),]$reference_key,1,Rx_radiation)) %>%
  mutate(px_coronary=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.0[3|4|9]|^36.0|^36.[2|3]|^36.9[1|9]"),]$reference_key,1,0)) %>% 
  mutate(Rx_h2=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^1.3.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_glycoside=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.1.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_diuretic=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.2.1|^2.2.3|^2.2.4"),]$reference_key,1,0)) %>% 
  mutate(Rx_AA=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.3.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_beta=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.4"),]$reference_key,1,0)) %>% 
  mutate(Rx_vasodilator=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_aHT=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_alpha=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.4"),]$reference_key,1,0)) %>% 
  mutate(Rx_RAS=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.5.[1|2|3]"),]$reference_key,1,0)) %>% 
  mutate(Rx_CCB=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.6.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_anticoagulant=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.8.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_antiplatelet=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.9"),]$reference_key,1,0)) %>% 
  mutate(Rx_LR=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.12"),]$reference_key,1,0)) %>% 
  mutate(Rx_NRT=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^4.10"),]$reference_key,1,0)) %>% 
  mutate(Rx_dementia=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^4.11"),]$reference_key,1,0)) %>% 
  mutate(Rx_anticdiabetic=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^6.1.1.[1|2]|^6.1.2.[1|2|3]"),]$reference_key,1,0)) %>% 
  mutate(Rx_corticosteroids=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^6.3.[1|2]"),]$reference_key,1,0)) %>% 
  mutate(Rx_NSAID=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^10.1.1"),]$reference_key,1,0)) %>% 
  mutate(px_orchidectomy=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^62.4"),]$reference_key,1,0)) %>% 
  mutate(px_prostatectomy=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^60\\.[2-6]"),]$reference_key,1,0))

mhspc.table<-mhspc.table[!is.na(mhspc.table$unified_results),] %>% 
  filter(index_date>="2013-1-1")
mhspc.table<-mhspc.table %>% 
  left_join(inpatient_num,by="reference_key") %>% 
  left_join(AnE_num,by="reference_key") %>% 
  left_join(outpatient_num[,c(1,18,19)],by="reference_key") %>% 
  left_join(inpatient_los,by="reference_key") %>% 
  left_join(inpatient_cost[,c(1,7)],by="reference_key") %>% 
  mutate(IP_yearlynum=IP_num/survtime*365.25,
         AE_yearlynum=AE_num/survtime*365.25,
         AE_yearlycost=AE_cost/survtime*365.25,
         OP_yearlynum=OP_num/survtime*365.25,
         OP_yearlycost=OP_cost/survtime*365.25,
         IP_yearlylos=IP_los/survtime*365.25,
         IP_yearlycost=ip_cost/survtime*365.25)
mhspc.table[,67:80]<-mhspc.table[,67:80] %>% 
  mutate(across(everything(), ~ replace_na(., 0)))
mhspc_NHA<-mhspc.table %>% 
  mutate(treat=ifelse(therapy=="Docetaxel",0,1)) 
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = mhspc_NHA,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("OS_mhspc_nha.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","NHA","KM Docetaxel","KM NHA"),lty=1:2,lwd=3, col=c("black","black","steelblue","steelblue"))
title(main="FPS for NHA vs Docetaxel in mHSPC")
text(1300,0.5,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
dev.off()
#HR
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#1.181 (95%CI 1.028 to 1.349)

fit<-survfit(Surv(survtime, death)~1,data=matched_data)
fit #median 468 days
#5year survival probability
summary_fit <- summary(fit, times = 1826.25)
survival_at_max_time <- summary_fit$surv[summary_fit$time == 1826.25]
surv_ci_lower <- summary_fit$lower[summary_fit$time == 365.25*5]
surv_ci_upper <- summary_fit$upper[summary_fit$time == 365.25*5]
print(paste("5-year Survival probability:", round(survival_at_max_time,2),"(95%CI",round(surv_ci_lower,2),"to",round(surv_ci_upper,2),")"))
#5-year Survival probability: 0.13 (95%CI 0.09 to 0.18 )

# Summarize the fit
max(matched_data$survtime) #3062
summary_fit <- summary(fit, times = max(matched_data$survtime))  # Use the maximum time in your data as a reference
# The survival probability at the end of the study
survival_at_max_time <- summary_fit$surv
surv_ci_lower <- summary_fit$lower[summary_fit$time == max(matched_data$survtime)]
surv_ci_upper <- summary_fit$upper[summary_fit$time == max(matched_data$survtime)]
print(paste("Survival probability at end:", round(survival_at_max_time,2),"(95%CI",round(surv_ci_lower,2),"to",round(surv_ci_upper,2),")"))
#Survival probability at end: 0.12 (95%CI 0.08 to 0.17 )


matched_data$treat<-as.factor(matched_data$treat)
treatment_labels<-c("0"="Docetaxel","1"="NHA")
matched_data_NHA<-matched_data%>% 
  select(reference_key,IP_yearlynum,AE_yearlynum,OP_yearlynum,IP_yearlylos,treat) %>% 
  mutate(source="Overall")
#IP
wilcox.test(matched_data$IP_yearlynum[matched_data$treat==1],matched_data$IP_yearlynum[matched_data$treat==0])

#IP LOS
wilcox.test(matched_data$IP_yearlylos[matched_data$treat==1],matched_data$IP_yearlylos[matched_data$treat==0])

#AE
wilcox.test(matched_data$AE_yearlynum[matched_data$treat==1],matched_data$AE_yearlynum[matched_data$treat==0])

#OP
wilcox.test(matched_data$OP_yearlynum[matched_data$treat==1],matched_data$OP_yearlynum[matched_data$treat==0],conf.int = T)


#secondary outcome: PSA progression--------
PSA_progression<-left_join(adjunct_rk,test_PSA,by="reference_key")
PSA_progression$index_date<-as.Date(PSA_progression$index_date)
PSA_progression_rk<-PSA_progression %>% 
  group_by(reference_key) %>% 
  filter(reference_date>index_date) %>% 
  summarise(test_count=n()) %>% 
  filter(test_count>=2)
PSA_progression<-PSA_progression %>% 
  group_by(reference_key) %>% 
  filter(reference_date>index_date) %>% 
  filter(reference_key%in%PSA_progression_rk$reference_key) %>% 
  rbind(PSA_baseline)
PSA_progression<-PSA_progression[!duplicated(PSA_progression),]
PSA_progression<-PSA_progression[!is.na(PSA_progression$unified_results),]
#PSA not decreased
find_psa_progression_1<-function(df){
  df<-df %>% arrange(reference_date)
  if(nrow(df)==1|| all(diff(df$unified_results) == 0)){
    return(NA)
  }
  if(any(diff(df$unified_results)<0)){
    return(NA)
  }
  baseline<-df$unified_results[1]
  threshold<-baseline*1.25
  ab_threshold<-baseline+5
  for (i in seq_along(df$unified_results)) {
    if (df$unified_results[i]>=threshold && df$unified_results[i]>=ab_threshold){
      return(df$reference_date[i])
    }
  }
  return(NA)
}
PSA_progression_1<-PSA_progression %>% 
  group_by(reference_key) %>% 
  summarize(progression_time=find_psa_progression_1(cur_data()))%>% 
  filter(!is.na(progression_time)) %>% 
  select(reference_key,progression_time)

#PSA decreased <=50%
find_psa_progression_2<-function(df){
  df<-df %>% arrange(reference_date)
  response_criteria<-0.5*df$unified_results[1]
  if(nrow(df)!=1 && min(df$unified_results) >= response_criteria && min(df$unified_results) < df$unified_results[1]){
    nadir_index<-which.min(df$unified_results)
    if(length(nadir_index)==0){
      return(NA)
    }
    nadir<-df$unified_results[nadir_index]
    progression_time<-NA
    for (i in nadir_index:nrow(df)) {
      threshold<-nadir*1.25
      ab_threshold<-nadir+5
      if(df$unified_results[i]>=threshold &&
         df$unified_results[i]>=ab_threshold){
        progression_time<-df$reference_date[i]
        break
      }
    }
    return(progression_time)
  }
  else{return(NA)}
}
PSA_progression_2<-PSA_progression %>% 
  group_by(reference_key) %>% 
  summarize(progression_time=find_psa_progression_2(cur_data())) %>% 
  filter(!is.na(progression_time)) %>% 
  select(reference_key,progression_time)

#PSA decreased >50%
find_psa_progression_3<-function(df){
  df<-df %>% arrange(reference_date)
  response_criteria<-0.5*df$unified_results[1]
  if(nrow(df)!=1 && min(df$unified_results) < response_criteria){
    nadir_index<-which.min(df$unified_results)
    if(length(nadir_index)==0){
      return(NA)
    }
    nadir<-df$unified_results[nadir_index]
    progression_time<-NA
    for (i in nadir_index:nrow(df)) {
      threshold<-nadir*1.5
      ab_threshold<-nadir+5
      if(df$unified_results[i]>=threshold &&
         df$unified_results[i]>=ab_threshold){
        progression_time<-df$reference_date[i]
        break
      }
    }
    return(progression_time)
  }
  else{return(NA)}
}
PSA_progression_3<-PSA_progression %>% 
  group_by(reference_key) %>% 
  summarize(progression_time=find_psa_progression_3(cur_data())) %>% 
  filter(!is.na(progression_time)) %>% 
  select(reference_key,progression_time)

#PSA response
find_psa_response <- function(data) {
  data <- data %>% arrange(reference_date)
  baseline_psa <- data$unified_results[1]
  response_threshold <- baseline_psa * 0.50
  potential_responses <- which(data$unified_results <= response_threshold)
  for (index in potential_responses) {
    confirmatory_tests <- data$reference_date > (data$reference_date[index] + 28) & data$unified_results <= response_threshold
    if (any(confirmatory_tests)) {
      confirmed_date <- min(data$reference_date[confirmatory_tests])
      return("1")
    }
  }
  return(NA)
}
psa_responses <- PSA_progression %>%
  group_by(reference_key) %>%
  arrange(reference_date,.by_group = TRUE) %>% 
  filter(unified_results[1]!=0) %>% 
  do(response = find_psa_response(.))
psa_responses <- psa_responses %>%
  filter(!is.na(response))

PSA_progression_date<-rbind(PSA_progression_1,PSA_progression_2,PSA_progression_3)
PSA_progression_date<-aggregate(PSA_progression_date$progression_time~PSA_progression_date$reference_key,PSA_progression_date,min)
colnames(PSA_progression_date)<-c("reference_key","PSA_progression")
PSA_progression_date<-PSA_progression_date %>% 
  right_join(PSA_progression_rk[,1],by="reference_key")
adjunct_psa<-adjunct_rk %>% 
  inner_join(PSA_progression_date,by="reference_key") %>% 
  left_join(psa_responses,by="reference_key") %>% 
  left_join(PSA_baseline[,c(1,19)],by="reference_key") %>% 
  left_join(Dx_list[,c(1,4)],by="reference_key") %>%
  mutate(survtime.psa=as.numeric(pmin(switch.1,switch.2,date_of_registered_death,endtime,PSA_progression,na.rm = T)-index_date)) %>% 
  mutate(death.psa=ifelse(is.na(date_of_registered_death),0,1)) %>% 
  mutate(death.psa=ifelse(is.na(PSA_progression),death.psa,1)) %>% 
  left_join(Dx_list[,c(1:3)],by="reference_key") %>% 
  mutate(DoB=as.Date(DoB),
         diagnosis_date=as.Date(diagnosis_date)) %>%
  mutate(age_at_diagnosis=as.numeric(diagnosis_date-DoB)/365.25,,
         age_at_dispensing=as.numeric(index_date-DoB)/365.25,
         disease_durtion=as.numeric(index_date-diagnosis_date)) %>% 
  mutate(dx_IHD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^41[0-3]|^414.0|^414.8|^414.9|^429.7|^V45.81|^V45.82"),]$reference_key,1,0)) %>% 
  mutate(dx_CAD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^43[0-8]"),]$reference_key,1,0)) %>% 
  mutate(dx_PVD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^44[0-3]"),]$reference_key,1,0)) %>% 
  mutate(px_CABG=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.1|^36.1[0-7]|^36.19|^36.2"),]$reference_key,1,0)) %>% 
  mutate(px_PCI=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.01|^36.02|^36.0[5-7]"),]$reference_key,1,0)) %>% 
  mutate(px_coronary=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.0[0|3|4|9]|^36.2|^36.3|^36.91|^36.99"),]$reference_key,1,0)) %>% 
  mutate(dx_HF=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^398.91|^402.01|^402.11|^402.91|^404.0[1|3]|^404.1[1|3]|^404.9[1|3]|^428|^V42.1"),]$reference_key,1,0)) %>% 
  mutate(dx_liver=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^070\\.2[2-3]|^070\\.3[2-3]|^070\\.[4-5]4|^070\\.[6|9]|^57[0-1]|^573\\.[3|4|8|9]|^V42\\.7|^456\\.[0-2]|^572\\.[2-4|8]"),]$reference_key,1,0)) %>% 
  mutate(dx_HT=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^40[1-5]"),]$reference_key,1,0)) %>% 
  mutate(dx_DM=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^249|^250|^357.2|^362.0|^366.41|^648.0"),]$reference_key,1,0)) %>% 
  mutate(dx_lipid=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^272.[0-4]"),]$reference_key,1,0)) %>% 
  mutate(dx_obesity=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^278"),]$reference_key,1,0)) %>% 
  mutate(dx_AF=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^427.31"),]$reference_key,1,0)) %>% 
  mutate(dx_arrhythmia=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^42[6-7]"),]$reference_key,1,0)) %>% 
  mutate(dx_cardiomypathy=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^425"),]$reference_key,1,0)) %>% 
  mutate(dx_COPD=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^49[1|2|6]"),]$reference_key,1,0)) %>% 
  mutate(dx_asthma=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^493"),]$reference_key,1,0)) %>% 
  mutate(dx_alcohol=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^291|^303|^305.0|^357.5|^425.5|^571.[0-3]|^980.[8|9]|^V11.3"),]$reference_key,1,0)) %>% 
  mutate(dx_smoke=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^305.1|^V15.82"),]$reference_key,1,0)) %>% 
  mutate(dx_depression=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^296.2|^296.3|^300.4|^311"),]$reference_key,1,0)) %>% 
  mutate(dx_hemiplagia=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^342"),]$reference_key,1,0)) %>% 
  mutate(dx_rheumatic=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^710|^714|^720|^725"),]$reference_key,1,0)) %>% 
  mutate(dx_dementia=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^290|^294.[0-2]|^331"),]$reference_key,1,0)) %>% 
  mutate(dx_renal=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^403|^58[2|3|5|6]|^588.0|^V42.0|^V45.1|^V56"),]$reference_key,1,0)) %>% 
  mutate(dx_ulcer=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^53[1-4]"),]$reference_key,1,0)) %>% 
  mutate(dx_AIDS=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^042|^079.53"),]$reference_key,1,0)) %>% 
  mutate(Rx_GnRHa=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("Leuprolide|LEUPRORELIN|Goserelin|TRIPTORELIN|OCTREOTIDE", ignore_case = T)),]$reference_key,1,0)) %>%
  mutate(Rx_GnRHat=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("Degarelix", ignore_case = T)),]$reference_key,1,0)) %>%
  mutate(Rx_ARB=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("BICALUTAMIDE|FLUTAMIDE", ignore_case = T)),]$reference_key,1,0)) %>% 
  mutate(Rx_cv_chemo=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$drug_name,regex("adriamycin|epirubicin|Trastuzumab", ignore_case = T)),]$reference_key,1,0)) %>% 
  mutate(Rx_chemo=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^8.[1|2]"),]$reference_key,1,0)) %>% 
  mutate(Rx_chemo=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^V58.1|^V66.2|^V67.2"),]$reference_key,1,Rx_chemo)) %>% 
  mutate(Rx_chemo=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^99.25"),]$reference_key,1,Rx_chemo)) %>% 
  mutate(Rx_radiation=ifelse(reference_key%in%Dx_baseline[str_which(Dx_baseline$all_diagnosis_code_ICD9,"^V58.0|^V66.1|^V67.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_radiation=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^92.[2-3][0-9]"),]$reference_key,1,Rx_radiation)) %>%
  mutate(px_coronary=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^36.0[3|4|9]|^36.0|^36.[2|3]|^36.9[1|9]"),]$reference_key,1,0)) %>% 
  mutate(Rx_h2=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^1.3.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_glycoside=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.1.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_diuretic=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.2.1|^2.2.3|^2.2.4"),]$reference_key,1,0)) %>% 
  mutate(Rx_AA=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.3.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_beta=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.4"),]$reference_key,1,0)) %>% 
  mutate(Rx_vasodilator=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.1"),]$reference_key,1,0)) %>% 
  mutate(Rx_aHT=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_alpha=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.4"),]$reference_key,1,0)) %>% 
  mutate(Rx_RAS=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.5.5.[1|2|3]"),]$reference_key,1,0)) %>% 
  mutate(Rx_CCB=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.6.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_anticoagulant=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.8.2"),]$reference_key,1,0)) %>% 
  mutate(Rx_antiplatelet=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.9"),]$reference_key,1,0)) %>% 
  mutate(Rx_LR=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^2.12"),]$reference_key,1,0)) %>% 
  mutate(Rx_NRT=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^4.10"),]$reference_key,1,0)) %>% 
  mutate(Rx_dementia=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^4.11"),]$reference_key,1,0)) %>% 
  mutate(Rx_anticdiabetic=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^6.1.1.[1|2]|^6.1.2.[1|2|3]"),]$reference_key,1,0)) %>% 
  mutate(Rx_corticosteroids=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^6.3.[1|2]"),]$reference_key,1,0)) %>% 
  mutate(Rx_NSAID=ifelse(reference_key%in%drug_baseline[str_which(drug_baseline$therapeutic_classification,"^10.1.1"),]$reference_key,1,0)) %>% 
  mutate(px_orchidectomy=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^62.4"),]$reference_key,1,0)) %>% 
  mutate(px_prostatectomy=ifelse(reference_key%in%Px_baseline[str_which(Px_baseline$code,"^60\\.[2-6]"),]$reference_key,1,0))

adjunct_psa_nona<-adjunct_psa[!is.na(adjunct_psa$unified_results),] %>% 
  filter(index_date>"2013-1-1")
##nha----
adjunct_NHA_psa<-adjunct_psa_nona %>% 
  mutate(treat=ifelse(therapy=="Docetaxel",0,1)) 
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_NHA_psa,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
love.plot(m.out, threshold = .1, which = "both")
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("PFS_nha.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)",ylab="PSA-progression free survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","NHA","KM Docetaxel","KM NHA"),lty=1:2, lwd = 3, col=c("black","black","steelblue","steelblue"))
title(main="NHA vs Docetaxel")
text(1800,0.5,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
dev.off()
#HR
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#1.018 (95%CI 0.906 to 1.139)

matched_data$response<-ifelse((matched_data$response=="1"),1,0)
matched_data$response<-as.factor(matched_data$response)
matched_data$treat<-as.factor(matched_data$treat)
response_table<-table(matched_data$treat,matched_data$response)
chi_test<-prop.test(response_table)
response_test<-prop.test(sum(matched_data$treat==1&matched_data$response==1),sum(matched_data$treat==1))
a<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
response_test<-binom.test(sum(matched_data$treat==0&matched_data$response==1),sum(matched_data$treat==0))
b<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
results<-data.frame(group="NHA",treat=a,control=b,p=round(chi_test$p.value,2))


fit<-survfit(Surv(survtime, death)~1,data=matched_data)
fit
max(matched_data$survtime) #3342
# Summarize the fit
summary_fit <- summary(fit, times = max(matched_data$survtime))  # Use the maximum time in your data as a reference
# The survival probability at the end of the study
survival_at_max_time <- summary_fit$surv
surv_ci_lower <- summary_fit$lower[summary_fit$time == max(matched_data$survtime)]
surv_ci_upper <- summary_fit$upper[summary_fit$time == max(matched_data$survtime)]
print(paste("Survival probability at the end of the study:", round(survival_at_max_time,2)))
#Survival probability at the end of the study: 0.09
print(paste("95% CI for the survival probability at the end of the study: [", 
            round(surv_ci_lower,2), ",", round(surv_ci_upper,2), "]", sep = ""))
#95% CI for the survival probability at the end of the study: [0.02,0.35]

# Summarize the fit
summary_fit <- summary(fit, times = 365.25*5)
# The survival probability at the end of the study
survival_at_max_time <- summary_fit$surv
surv_ci_lower <- summary_fit$lower[summary_fit$time == 365.25*5]
surv_ci_upper <- summary_fit$upper[summary_fit$time == 365.25*5]
print(paste("Survival probability at 5 years:", round(survival_at_max_time,2)))
#Survival probability at the end of the study: 0.2
print(paste("95% CI for the survival probability at 5 years: [", 
            round(surv_ci_lower,2), ",", round(surv_ci_upper,2), "]", sep = ""))
#95% CI for the survival probability at the end of the study: [0.16,0.25]


##abiraterone----
adjunct_AT_psa<-rbind(adjunct_psa_nona[adjunct_psa_nona$therapy=="Abiraterone",],
                  adjunct_psa_nona[adjunct_psa_nona$therapy=="Docetaxel",])
adjunct_AT_psa$treat<-ifelse(adjunct_AT_psa$therapy=="Abiraterone",1,0)
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_AT_psa,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
love.plot(m.out, threshold = .1, which = "both")
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("PFS_abi.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)",ylab="PSA-progression free survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","abiraterone","KM Docetaxel","KM abiraterone"),lty=1:2,lwd = 3, col=c("black","black","steelblue","steelblue"))
title(main="Abiraterone vs Docetaxel")
text(1800,0.5,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
dev.off()
#HR
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#0.955 (95%CI 0.836 to 1.085)

matched_data$response<-ifelse((matched_data$response=="1"),1,0)
matched_data$response<-as.factor(matched_data$response)
matched_data$treat<-as.factor(matched_data$treat)
response_table<-table(matched_data$treat,matched_data$response)
chi_test<-prop.test(response_table)
response_test<-prop.test(sum(matched_data$treat==1&matched_data$response==1),sum(matched_data$treat==1))
a<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
response_test<-binom.test(sum(matched_data$treat==0&matched_data$response==1),sum(matched_data$treat==0))
b<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
results<-rbind(results,
               data.frame(group="abiraterone",treat=a,control=b,p=round(chi_test$p.value,2)))


##enzalutamide----
adjunct_ET_psa<-rbind(adjunct_psa_nona[adjunct_psa_nona$therapy=="Enzalutamide",],
                      adjunct_psa_nona[adjunct_psa_nona$therapy=="Docetaxel",])
adjunct_ET_psa$treat<-ifelse(adjunct_ET_psa$therapy=="Enzalutamide",1,0)
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_ET_psa,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("PFS_enz.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)",ylab="PSA-progression free survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","enzalutamide","KM Docetaxel","KM enzalutamide"),lty=1:2, lwd = 3,col=c("black","black","steelblue","steelblue"))
title(main="Enzalutamide vs Docetaxel")
text(1200,0.55,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
dev.off()
#HR
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#0.878 (95%CI 0.734 to 1.038)"

matched_data$response<-ifelse((matched_data$response=="1"),1,0)
matched_data$response<-as.factor(matched_data$response)
matched_data$treat<-as.factor(matched_data$treat)
response_table<-table(matched_data$treat,matched_data$response)
chi_test<-prop.test(response_table)
response_test<-prop.test(sum(matched_data$treat==1&matched_data$response==1),sum(matched_data$treat==1))
a<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
response_test<-binom.test(sum(matched_data$treat==0&matched_data$response==1),sum(matched_data$treat==0))
b<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
results<-rbind(results,
               data.frame(group="enzalutamide",treat=a,control=b,p=round(chi_test$p.value,2)))

##old-----
adjunct_75_psa<-adjunct_NHA_psa[adjunct_NHA_psa$age_at_diagnosis>=75,]
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_75_psa,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
love.plot(m.out, threshold = .1, which = "both")
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("PFS_old.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)",ylab="PSA-progression free survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","NHA","KM Docetaxel","KM NHA"),lty=1:2,lwd = 3, col=c("black","black","steelblue","steelblue"))
title(main="NHA vs Docetaxel (age>=75)")
text(1050,0.55,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
dev.off()
#HR
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#0.694 (95%CI 0.427 to 1.054)

matched_data$response<-ifelse((matched_data$response=="1"),1,0)
matched_data$response<-as.factor(matched_data$response)
matched_data$treat<-as.factor(matched_data$treat)
response_table<-table(matched_data$treat,matched_data$response)
chi_test<-prop.test(response_table)
response_test<-prop.test(sum(matched_data$treat==1&matched_data$response==1),sum(matched_data$treat==1))
a<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
response_test<-binom.test(sum(matched_data$treat==0&matched_data$response==1),sum(matched_data$treat==0))
b<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
results<-rbind(results,
               data.frame(group="old",treat=a,control=b,p=round(chi_test$p.value,2)))

##young----
adjunct_74_psa<-adjunct_NHA_psa[adjunct_NHA_psa$age_at_diagnosis<75,]
m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
               data = adjunct_74_psa,method = "nearest",caliper = 0.1)
matched_data <- match.data(m.out)
love.plot(m.out, threshold = .1, which = "both")
model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
jpeg("PFS_young.jpeg", width = 1000, height = 900,quality=1000,pointsize = 35)
plot(model, newdata=data.frame(treat=0), xlab="Time since drug dispensing (days)",ylab="PSA-progression free survival",lwd=2)
lines(model, newdata=data.frame(treat=1), lty=2,lwd=2)
lines(survfit(Surv(matched_data$survtime, matched_data$death) ~ matched_data$treat), col="steelblue", lty=1:2,lwd=4)
legend("topright", c("Docetaxel","NHA","KM Docetaxel","KM NHA"),lty=1:2, lwd=3, col=c("black","black","steelblue","steelblue"))
title(main="NHA vs Docetaxel (age<75)")
text(1900,0.5,labels=sprintf("HR %.2f (95%%CI %.2f to %.2f )",
                             exp(coef(model)["treat"]),
                             exp(confint(model)["treat", 1]),
                             exp(confint(model)["treat", 2])))
dev.off()
#HR
sprintf("%.3f (95%%CI %.3f to %.3f)",
        exp(coef(model)["treat"]),
        exp(confint(model)["treat", 1]),
        exp(confint(model)["treat", 2]))
#0.903 (95%CI 0.795 to 1.019)


#Sensitive analysis------
##Remove survtime<30----
data<-list(main=adjunct_NHA_psa,abi=adjunct_AT_psa,enz=adjunct_ET_psa,old=adjunct_75_psa,young=adjunct_74_psa)
shortremove_results<-data.frame(group=NA,rm_30=NA)
for (i in names(data)) {
  adjunct_NHA_30<-data[[i]]%>% 
    filter(survtime>30)
  m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
                 data = adjunct_NHA_30,method = "nearest",caliper = 0.1)
  matched_data <- match.data(m.out)
  model <- stpm2(Surv(survtime,death) ~ treat, data = matched_data)
  nha_30<-sprintf("%.3f (95%%CI %.3f to %.3f)",
                  exp(coef(model)["treat"]),
                  exp(confint(model)["treat", 1]),
                  exp(confint(model)["treat", 2]))
  shortremove_results<-rbind(shortremove_results,
                             data.frame(group=i,
                                        rm_30=sprintf("%.3f (95%%CI %.3f to %.3f)",
                                                      exp(coef(model)["treat"]),
                                                      exp(confint(model)["treat", 1]),
                                                      exp(confint(model)["treat", 2]))))
}
shortremove_results<-shortremove_results[-1,]
##IPTW------
library(ipw)
iptw_results<-data.frame(group=NA,IPTW=NA)
for (i in names(data)) {
  adjunct_name<-data[[i]]
  colnames(adjunct_name)[18:71]<-paste0("cov_", 1:54)
  w1<-ipwpoint(
    exposure = treat,
    family = "binomial",
    link = "logit",
    numerator = ~1,
    denominator = ~cov_1+cov_3+cov_4+cov_5+cov_6+cov_7+cov_8+cov_9+cov_10+cov_11+cov_12+cov_13+cov_14+cov_15+cov_16+cov_17+cov_18+cov_19+cov_20+cov_21+cov_22+cov_23+cov_24+cov_25+cov_26+cov_27+cov_28+cov_29+cov_30+cov_31+cov_32+cov_33+cov_34+cov_35+cov_36+cov_37+cov_38+cov_39+cov_40+cov_41+cov_42+cov_43+cov_44+cov_45+cov_46+cov_47+cov_48+cov_49+cov_50+cov_51+cov_52+cov_53+cov_54,
    data = adjunct_name
  )
  adjunct_name$w1<-w1$ipw.weights
  survival_model<-coxph(Surv(survtime,death)~treat,data=adjunct_name,weights=w1)
  iptw_results<-rbind(iptw_results,
                      data.frame(group=i,
                                 IPTW=sprintf("%.3f (95%%CI %.3f to %.3f)",
                                              exp(coef(survival_model)["treat"]),
                                              exp(confint(survival_model)["treat", 1]),
                                              exp(confint(survival_model)["treat", 2]))))
}
iptw_results<-iptw_results[-1,]


##IPCW------
ipcw_results<-data.frame(group=NA,IPCW=NA)
for (i in names(data)) {
  m.out<-matchit(treat~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_CABG+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_obesity+dx_AF+dx_arrhythmia+dx_cardiomypathy+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_dementia+dx_renal+dx_ulcer+dx_AIDS+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_cv_chemo+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
                 data = data[[i]],method = "nearest",caliper = 0.1)
  matched_data <- match.data(m.out)
  adjunct_IPCW<-adjunct_rk %>% 
    mutate(censor=ifelse(date_of_registered_death==pmin(switch.1,switch.2,date_of_registered_death,endtime,na.rm = T),0,1)) %>% 
    mutate(censor=ifelse(is.na(date_of_registered_death),1,censor)) %>% 
    right_join(matched_data[,c(1,12:72)],by="reference_key")
  c0<-coxph(Surv(survtime,censor)~1,data = adjunct_IPCW)
  c0fit<-summary(survfit(c0),times = adjunct_IPCW$survtime)
  adjunct_IPCW$k0ti<-c0fit$surv
  cz<-coxph(Surv(survtime,censor)~unified_results+age_at_dispensing+dx_IHD+dx_CAD+dx_PVD+px_PCI+px_coronary+dx_HF+dx_liver+dx_HT+dx_DM+dx_lipid+dx_AF+dx_arrhythmia+dx_COPD+dx_asthma+dx_alcohol+dx_smoke+dx_depression+dx_hemiplagia+dx_rheumatic+dx_renal+dx_ulcer+Rx_GnRHa+Rx_GnRHat+Rx_ARB+Rx_chemo+Rx_radiation+Rx_h2+Rx_glycoside+Rx_diuretic+Rx_AA+Rx_beta+Rx_vasodilator+Rx_aHT+Rx_alpha+Rx_RAS+Rx_CCB+Rx_anticoagulant+Rx_antiplatelet+Rx_LR+Rx_NRT+Rx_dementia+Rx_anticdiabetic+Rx_corticosteroids+Rx_NSAID+px_orchidectomy+px_prostatectomy,
            data=adjunct_IPCW)
  process_row <- function(i) {
    datai <-adjunct_IPCW[i, ]
    czfit <- summary(survfit(cz, newdata = datai), times = datai$survtime)
    return(czfit$surv)
  }
  results <- mclapply(1:nrow(adjunct_IPCW), process_row, mc.cores = 10)
  adjunct_IPCW$kzti <- unlist(results)
  adjunct_IPCW<-adjunct_IPCW %>% 
    mutate(ipcw=1/kzti,
           ipcw_s=k0ti/kzti) %>% 
    mutate(ipcw_s_upper=quantile(ipcw_s, 0.95),
           ipcw_s_lower=quantile(ipcw_s, 0.05)) %>%
    mutate(ipcw_s=ifelse(ipcw_s>=ipcw_s_upper,ipcw_s_upper,ipcw_s),
           ipcw_s=ifelse(ipcw_s<=ipcw_s_lower,ipcw_s_lower,ipcw_s))
  survival_model<-coxph(Surv(survtime,death)~treat,data=adjunct_IPCW,weights=ipcw_s)
  ipcw_results<-rbind(ipcw_results,
                      data.frame(group=i,
                                 IPCW=sprintf("%.3f (95%%CI %.3f to %.3f)",
                                              exp(coef(survival_model)["treat"]),
                                              exp(confint(survival_model)["treat", 1]),
                                              exp(confint(survival_model)["treat", 2]))))
  
}
ipcw_results<-ipcw_results[-1,]
sensitive_results<-shortremove_results %>% 
  left_join(iptw_results,by="group") %>% 
  left_join(ipcw_results,by="group")
export(sensitive_results,"sensitive_results_psa.xlsx")

##PSA response table----
matched_data$response<-ifelse((matched_data$response=="1"),1,0)
matched_data$response<-as.factor(matched_data$response)
matched_data$treat<-as.factor(matched_data$treat)
response_table<-table(matched_data$treat,matched_data$response)
chi_test<-prop.test(response_table)
response_test<-prop.test(sum(matched_data$treat==1&matched_data$response==1),sum(matched_data$treat==1))
a<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
response_test<-binom.test(sum(matched_data$treat==0&matched_data$response==1),sum(matched_data$treat==0))
b<-sprintf("%.2f (95%%CI %.2f to %.2f)",
           response_test$estimate,
           response_test$conf.int[1],
           response_test$conf.int[2])
results<-rbind(results,
               data.frame(group="young",treat=a,control=b,p=round(chi_test$p.value,2)))


export(results,"PSA_response_table.xlsx")

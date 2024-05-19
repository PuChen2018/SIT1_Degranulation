library(ggplot2)
library(tidyverse)
# read the file 
rm(list = ls())
files<-list.files("../../../Analysis_killing_assay_incucyte_20240219/",".E[0-9]T[0-9].txt",full.names = FALSE)
list_of_files<-list()
for (file in files){
  if(!file == "20230705_killing_withIL2_4pics_E2T1.txt"){ 
    tbl0<-read.csv2(str_c("../../../Analysis_killing_assay_incucyte_20240219/",file),sep ="\t",skip = 7)
  }else{tbl0<-read.table(str_c("../../../Analysis_killing_assay_incucyte_20240219/",file),header = TRUE,sep="\t")}
  list_of_files[[file]]<-tbl0
}



# combined all the tables into one 
tbl1<-tibble(dir=files)%>%
  mutate(IL2_starvation=ifelse(str_detect(dir,"withIL2"),"withIL2","withoutIL2"))%>%
  mutate(batch=case_when(
    str_detect(dir,"20221020")~1,
    str_detect(dir,"2022111")~2,
    str_detect(dir,"202303")~3,
    str_detect(dir,"2023070")~4
  ) )%>%
  mutate(ratio=case_when(
    str_detect(dir,"E1T1")~"E/T=1:1",
    str_detect(dir,"E2T1")~"E/T=2:1",
    TRUE~"NA"
  ))%>%
  group_by(batch,IL2_starvation,ratio,dir) %>%
  reframe(list_of_files[[dir]])


# reannotate all the datapoint
# all the co-culture condition will be put into the ratio 
tbl2<-tbl1%>%
  pivot_longer(matches('^[A-Z][0-9]+$'))%>%
  rename(GFP_count = value)%>%
  separate_wider_regex(name,c(rows='[A-Z]',column ='\\d+'))%>%
  mutate(E_T_ratio = case_when(
    batch==1&column %in% c(2,3,4,5)~"E/T=1:1",
    batch==1&column%in%c(6,7,8,9)~"E/T=2:1",
    batch==1&column%in%c(10,11)~"target_cells_only",
    batch!=1&rows=="G"~"target_cells_only",
    batch!=1&column==11~"target_cells_only",
    TRUE~ratio
  ))%>% 
  mutate(Aphidicolin_treatment = case_when(
    column %in% c(11) ~ "no_aphidicolin",
    TRUE~"aphidicolin_treated"
  ))%>% ## coculture condition
  mutate(co_culture = ifelse(E_T_ratio=="target_cells_only","target_cells_only","co_culture"))%>%
  ## annotate the individua, batch 1, 2, 3 both the with IL2 and without IL2 were the same but the 
  ## batch 4 has a bit shift 
  mutate(individual=case_when(
    batch!=1& co_culture=="target_cells_only" & Aphidicolin_treatment == "aphidicolin_treated"~"target_only_with_aphi",
    batch!=1& co_culture== "target_cells_only"& Aphidicolin_treatment == "no_aphidicolin"~"target_only_without_aphi",
    batch==1& column %in%c(2,6)~"PID1749",
    batch==1& column %in%c(3,7)~"Patient",
    batch==1& column %in%c(4,8)~"PID1604",
    batch==1& column %in%c(5,9)~ "PID1691",# PID1691#"IGG039"
    batch==2& rows =="B"~"PID1605",
    batch==2& rows =="C"~"PID1749",
    batch==2& rows =="D"~"Patient",
    batch==2& rows =="E"~"PID1604",# =PIG323
    batch==2& rows =="F"~"PID1691",
    batch==3& rows =="B"~"PID1751",
    batch==3& rows =="C"~"Patient",
    batch==3& rows =="D"~"PID1802",
    batch==3& rows =="E"~"PID1605",
    batch==3& rows =="F"~"PIG351",
    batch==4& rows =="B" & IL2_starvation == "withoutIL2"~"PID1605",
    batch==4& rows =="C" & IL2_starvation == "withoutIL2"~"Patient",
    batch==4& rows =="D" & IL2_starvation == "withoutIL2"~"PID1666",
    batch==4& rows =="E" & IL2_starvation == "withoutIL2"~"PID1802",
    batch==4& rows =="F" & IL2_starvation == "withoutIL2"~"PID1820",
    batch==4& rows =="B" & IL2_starvation == "withIL2"~"PID1605",
    batch==4& rows =="C" & IL2_starvation == "withIL2"~"PID1666",
    batch==4& rows =="D" & IL2_starvation == "withIL2"~"Patient",
    batch==4& rows =="E" & IL2_starvation == "withIL2"~"PID1820",
    batch==4& rows =="F" & IL2_starvation == "withIL2"~"PID1802",
    batch ==1& column ==10 ~"target_only_with_aphi",
    batch ==1 & column ==11 ~ "target_only_without_aphi"
  ))%>%
  ## annotate coating 
  mutate(`anti_CD3_µg/ml` = case_when(
    batch==1& rows %in% c("B","C")~ 0,
    batch==1& rows %in% c("D","E")~ 1,
    batch==1& rows %in% c("F","G")~ 5,
    batch!=1& column %in%c(2,3,4)~ 0,
    batch!=1& column %in%c(5,6,7)~ 1,
    batch!=1& column %in%c(8,9,10)~ 5,
    batch!=1& column==11&rows %in%c("B","C")~0,
    batch!=1& column==11&rows %in%c("D","E")~1,
    batch!=1& column==11&rows %in%c("F","G")~5,
    TRUE ~ NA
  )) %>%
  mutate(genotype = ifelse(individual == "Patient","Patient","Controls"))

# clean up the data that has glithes
tbl3<-tbl2%>%filter(GFP_count> 50)

# method two 
Colors2<-c("PID1749"="#DDCC77","Patient"="blue",
           "PID1604"="#CC6677","PID1691"="#117733",
           "target_only_with_aphi"="#882255","target_only_without_aphi"="#44AA99",
           "PID1605"="green","PID1751"="#88CCEE",
           "PID1802"="#AA4499","PIG351"="#661100",
           "PID1666"="#888888","PID1820"="#6699cc")

# normlized the data 
tbl4<-tbl3%>%
  filter(batch == 1& Elapsed >= 0 &Elapsed < 16|
           batch == 2& Elapsed > 0 & Elapsed < 16|
           batch == 3& Elapsed > 0 & Elapsed <16|
           batch == 4& Elapsed > 0 & Elapsed < 16) 

# make plot from the normalized table
for (bth in unique(tbl4$batch)){
  tbl4%>% filter(batch==bth) -> batch_tbl
  batch_tbl%>%group_by(batch,IL2_starvation,E_T_ratio,individual,Elapsed,`anti_CD3_µg/ml`,Aphidicolin_treatment)%>%summarize(GFP_count_mean=mean(GFP_count)) -> means_batch_tbl
  means_batch_tbl%>%group_by(batch, IL2_starvation, E_T_ratio,individual,`anti_CD3_µg/ml`)%>%arrange(Elapsed)%>%
    mutate(normalized_GFP_count_mean=GFP_count_mean/first(GFP_count_mean)) -> normalized_means_batch_tbl
  print( ggplot(NULL)+
           geom_line(data=normalized_means_batch_tbl,aes(x=Elapsed,y=normalized_GFP_count_mean,group = interaction (individual, `anti_CD3_µg/ml`),col= individual))+
           facet_grid(E_T_ratio+Aphidicolin_treatment~IL2_starvation+`anti_CD3_µg/ml`)+
           scale_color_manual(values=Colors2)+
           ggtitle(bth))
}


# plot only for batch that worked 

tbl4%>% 
  filter(batch==2) %>%
  filter(E_T_ratio != "target_cells_only") -> batch_tbl

batch_tbl%>%
  group_by(batch,IL2_starvation,E_T_ratio,individual,Elapsed,`anti_CD3_µg/ml`,Aphidicolin_treatment,genotype)%>%
  summarize(GFP_count_mean=mean(GFP_count)) -> means_batch_tbl

means_batch_tbl%>%
  group_by(batch, IL2_starvation, E_T_ratio,individual,`anti_CD3_µg/ml`,genotype)%>%
  arrange(Elapsed)%>%
  mutate(normalized_GFP_count_mean=GFP_count_mean/first(GFP_count_mean)) -> normalized_means_batch_tbl

p_killing <-ggplot(NULL)+
         geom_line(data=normalized_means_batch_tbl,aes(x=Elapsed,y=normalized_GFP_count_mean,group = interaction (individual, `anti_CD3_µg/ml`),col= genotype))+
         facet_grid(E_T_ratio ~IL2_starvation+`anti_CD3_µg/ml`)+
         scale_color_manual(values=c("Patient"="blue","Controls"= "gray")) +
         theme_classic(base_size = 10)+
         theme(
           panel.background = element_rect(fill = "white", color = "lightgrey"),  # Set background color
           panel.grid.major.y = element_line(size = 0.5, color = "lightgrey"),  # Adjust grid line appearance
           #panel.grid.minor = element_line(size = 0.25, color = "lightgrey"),
           strip.background = element_blank(),
           strip.text.x = element_text(size = rel(1.2)),
           strip.text.y = element_text(size = rel(1.2)),
           text = element_text(family="Times"),
           axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5,size =rel(1)),
           legend.title = element_text(size = rel(1),face = "bold")) +
         guides(fill = guide_legend(
           title = "Samples",
           keywidth = unit(1, "lines"), 
           keyheight = unit(1, "lines") ))+
         labs(y = "normalized GFP count to 0h",
              x = "Elapsed time")
ggsave("killillng_assay.png", plot = p_killing, width = 18, height = 10, units = "cm", dpi = 300)

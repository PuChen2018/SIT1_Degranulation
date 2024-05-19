
library(ggplot2)
library(tidyverse)

files<-list.files(path = "../../Fc_blocker_absent/",recursive = TRUE,full.names = FALSE, ".txt" )
list_of_files<-list()

# read all the files and put them into a list, since some of the columns didnt come out as numbers so there is a step to convert then into numerical data
for (file in files){
  tbl0<-read_tsv(str_c("../../Fc_blocker_absent/",file))
  tbl0_col<-colnames(tbl0)
  newcol<-str_replace(tbl0_col,"/FVD(660)*, (FSC-H subset)*/CD3(-BV605)*, (FSC-H subset)*/Q3: CD8(-PE-cy7)*\\+ , CD4(-BV421)*-","/FVD660-/CD3+/Q3:CD8+CD4-")
  newcol2<-str_replace(newcol, "CD107[A-Za-z0-9_-](FITC)*, CD8(-PE-cy7)* subset ","CD107+CD8+")
  newcol3<-str_replace(newcol2,"Mean \\(CD107-FITC\\)","Mean (CD107a)")
  colnames(tbl0)<-newcol3
  tbl0[2:4]<-lapply(tbl0[2:4],as.numeric)
  print(tbl0_col)
  print(newcol3)
  list_of_files[[file]]<-tbl0
}


# make the master table 
## Table of files
tibble(dir = names(list_of_files) ) %>%
  mutate(IL2_starvation = ifelse(str_detect(dir,"withIL2"),"withIL2","withoutIL2")) %>%
  mutate(batch= ifelse(str_detect(dir,"202210"),1,2))%>%
  mutate(file_name= dir)%>%
  rowwise(IL2_starvation, batch,file_name)%>%
  reframe(list_of_files[[dir]]) -> tbl1


## try to rename and make proper annotation for each treatment 
tbl2<-tbl1 %>%
  filter(str_detect(`Sample:`,"fcs")) %>%
  mutate(row=str_sub(`Sample:`,6,6))%>%
  mutate(column=str_sub(`Sample:`,7,8))%>%
  filter(!(column %in%c("09","10","11","12")&row %in% c("G","H")))%>% # these are the FMO can be deleted 
  mutate(ratio=ifelse(row %in% c("A", "B", "C"),"E/T=1:2","E/T=2:1"))%>%
  mutate(ratio=ifelse(row %in%c("G","H"),"control_rows",ratio))%>%
  mutate(genotype = case_when(
    batch==1 & column %in%c("01","05","09") ~ "patient",
    batch==2 & column %in%c("02","06","10") ~ "patient",
    TRUE ~ "controls"
  ))%>%
  pivot_longer( contains("|") ) %>%
  separate_wider_delim(name, "|",names=c("gating","stat"))%>%
  mutate(gating=str_replace(gating,"Cells/Single Cells/FVD660-/CD3\\+/",""))%>%
  mutate(anti_CD3 =ifelse(column %in% c("01","02","03","04"),"non-coated","1ug/ml"))%>%
  mutate(anti_CD3=ifelse(column %in% c("09","10","11","12"),"5ug/ml",anti_CD3))%>%
  mutate(anti_CD3=ifelse(column %in% c("01","02","03","04")& row %in% c("G","H"), "basal",anti_CD3))%>%
  mutate(anti_CD3=ifelse(column %in% c("05","06","07","08")& row %in% c("G","H"), "PMA+Ino",anti_CD3))%>%
  mutate(anti_CD3 = factor(anti_CD3, levels = c("non-coated","1ug/ml","5ug/ml","basal","PMA+Ino")))


## make the annotation for the individuals 
tbl3_1<-tbl2%>%
  mutate(individual= ifelse(batch==1  & column %in% c("01","05","09"),"Patient","NA"))%>%
  mutate(individual= ifelse(batch==1  & column %in% c("02","06","10"),"PID1749",individual))%>%
  mutate(individual= ifelse(batch==1  & column %in% c("03","07","11"),"PID1604",individual))%>%
  mutate(individual = ifelse(batch==1 & column %in% c("04","08","12"),"PID1691", individual))%>% #IGG039
  mutate(individual = ifelse (batch== 2 & column %in% c("01","05","09"),"PID1749", individual))%>%
  mutate(individual = ifelse (batch== 2 & column %in% c("02","06","10"),"Patient", individual))%>%
  mutate(individual = ifelse (batch== 2 & column %in% c("03","07","11"),"PID1691", individual))%>%
  mutate(individual = ifelse (batch== 2 & column %in% c("04","08","12"),"PID1604", individual))  # PIG323



tbl3<-tbl3_1
tbl3 %>%
  filter(gating =="Q3:CD8+CD4-/CD107+CD8+" & stat == " Freq. of Parent") -> tbl3_CD107a_percen

tbl3%>%
  select(batch,individual,IL2_starvation)%>% 
  distinct()%>% group_by(IL2_starvation,batch)%>% 
  mutate(x.pos= row_number())-> 
  position_integer

Colors2<-c("PID1749"="#DDCC77","Patient"="blue",
           "PID1604"="#CC6677","PID1691"="#117733",
           "target_only_with_aphi"="#882255","target_only_without_aphi"="#44AA99",
           "PID1605"="green","PID1751"="#88CCEE",
           "PID1802"="#AA4499","PIG351"="#661100",
           "PID1666"="#888888","PID1820"="#6699cc")

for (b_no in unique(tbl3_CD107a_percen$batch)){
  tbl3_CD107a_percen%>%
    filter(batch== b_no) ->tbl3_CD107a_percen_batch
  tbl3_CD107a_percen_batch%>% left_join(position_integer) %>% mutate(x=as.integer(factor(anti_CD3))+x.pos/7) -> tbl3_CD107a_percen_batch
  tbl3_CD107a_percen_batch%>% group_by(batch,IL2_starvation,ratio,individual,genotype,x) %>% summarise(value=mean(value))-> means_CD107a_percen_batch
  print(ggplot(NULL)+
          geom_point(aes(x=x,y=value,col=individual,shape = genotype),data=tbl3_CD107a_percen_batch)+
          geom_errorbar(aes(x=x,ymin=value,ymax=value, col = individual), data= means_CD107a_percen_batch) +
          facet_grid(ratio+stat ~ batch + IL2_starvation, scales = "free_y") +
          scale_color_manual(values = Colors2)+
          scale_x_continuous(breaks = (1:5) + 3/7, labels = levels(factor(tbl3_CD107a_percen_batch$anti_CD3))) +
          theme(
            panel.background = element_rect(fill = "white", color = "lightgrey"),  # Set background color
            panel.grid.major = element_line(size = 0.5, color = "lightgrey"),  # Adjust grid line appearance
            panel.grid.minor = element_line(size = 0.25, color = "lightgrey"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          ggtitle(str_c(tbl3_CD107a_percen_batch$gating,"_",tbl3_CD107a_percen_batch$file_name))
        
  )
}  

########################################################################################
#           change to the batch number you want to plot                                #
#                                                                                      #
########################################################################################
# select a batch that you want to make the figure
# put a column that you can indicate the position of the dots 
# if you want to choose a batch, choose here 
tbl3_CD107a_percen%>%
  filter(batch== 2) %>%
  left_join(position_integer) %>%
  mutate(x=as.integer(factor(anti_CD3))+x.pos/7) -> tbl3_CD107a_percen_batch

tbl3_CD107a_percen_batch%>% 
  group_by(batch,IL2_starvation,ratio,individual,genotype,x, anti_CD3,stat) %>% 
  summarise(value_mean=mean(value))-> means_CD107a_percen_batch

# make a plot with more individual details but not necessary need to be shown in the final manuscript 
ggplot(NULL)+
  geom_point(aes(x=x,y=value,col=individual,shape = genotype),data=tbl3_CD107a_percen_batch)+
  geom_errorbar(aes(x=x,ymin=value_mean,ymax=value_mean, col = individual), data= means_CD107a_percen_batch) +
  facet_grid(ratio+stat ~ batch + IL2_starvation, scales = "free_y") +
  scale_color_manual(values = Colors2)+
  scale_x_continuous(breaks = (1:5) + 3/7, labels = levels(factor(tbl3_CD107a_percen_batch$anti_CD3))) +
  theme_classic(base_size = 10)+ 
  theme(
    panel.background = element_rect(fill = "white", color = "lightgrey"),  # Set background color
    panel.grid.major = element_line(size = 0.5, color = "lightgrey"),  # Adjust grid line appearance
    panel.grid.minor = element_line(size = 0.25, color = "lightgrey"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1,size =rel(1))) +
  ggtitle(tbl3_CD107a_percen_batch$gating)

# figure for manuscript testing, dot with mean as bars 
ggplot(NULL)+
  geom_point(aes(x=x,y=value,col= genotype,shape = genotype),data= tbl3_CD107a_percen_batch)+
  geom_errorbar(aes(x=x,ymin=value_mean,ymax=value_mean, col = genotype ), data= means_CD107a_percen_batch) +
  facet_grid(ratio+stat ~ batch + IL2_starvation, scales = "free_y") +
  scale_color_manual(values = c("patient"="blue","controls"="gray"))+
  scale_x_continuous(breaks = (1:5) + 3/7, labels = levels(factor(tbl3_CD107a_percen_batch$anti_CD3))) +
  theme_classic(base_size = 10) + 
  theme(
    panel.background = element_rect(fill = "white", color = "lightgrey"),  # Set background color
    panel.grid.major = element_line(size = 0.5, color = "lightgrey"),  # Adjust grid line appearance
    panel.grid.minor = element_line(size = 0.25, color = "lightgrey"),
    strip.background = element_blank(),
    strip.text.x = element_text(size =rel(1.2)),
    axis.text.x = element_text(angle = 50, vjust = 0.5, hjust = 0.5,size =rel(1.5)),
    
    text = element_text(family="Times")) +
  ggtitle(tbl3_CD107a_percen_batch$gating)


# figure for manuscript testing if I use bar plot 
p_degranulation<-ggplot(data = means_CD107a_percen_batch, aes(x = anti_CD3, y = value_mean,fill = genotype,group = individual)) +
  geom_bar(stat = "identity",position = position_dodge(width = 0.8,preserve = "single"),color = "black",size =0.3,width = 0.8) +
  facet_grid( IL2_starvation ~ratio, scales = "free",space="free")+
  scale_fill_manual(values = c("patient"="blue","controls"="grey")) +
  theme_classic(base_size = 12)+ 
  theme(
    panel.background = element_rect(fill = "white", color = "lightgrey"),  # Set background color
    panel.grid.major.y = element_line(size = 0.5, color = "lightgrey"),  # Adjust grid line appearance
    #panel.grid.minor = element_line(size = 0.25, color = "lightgrey"),
    strip.background = element_blank(),
    #strip.text = element_blank()
    strip.text.x = element_text(size =rel(1)),
    strip.text.y = element_text(size = rel(1.2)),
    text = element_text(family="Times"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5,size =rel(0.8)),
    legend.title = element_text(size = rel(1),face = "bold")) +
  guides(fill = guide_legend(
    title = "Samples",
    keywidth = unit(1, "lines"), 
    keyheight = unit(1, "lines") ))+
  labs(x = "anti-CD3",
       y = "CD107+CD8+% in CD8+")
ggsave("degranulation_single_batch_2.png", plot = p_degranulation, width = 18, height = 10, units = "cm", dpi = 300)
# for batch 1, change the height to 5 since it only has one condition
#template for ggplot to tweak everything for labels 

ggplot(NULL)+
  geom_point(aes(x=x,y=value,col= genotype,shape = genotype),data=tbl3_CD107a_percen_batch)+
  geom_errorbar(aes(x=x,ymin=value,ymax=value, col = genotype ), data= means_CD107a_percen_batch) +
  facet_grid(ratio+stat ~ batch + IL2_starvation, scales = "free_y") +
  scale_color_manual(values = c("patient"="blue","control"="gray"))+
  scale_x_continuous(breaks = (1:5) + 3/7, labels = levels(factor(tbl3_CD107a_percen_batch$anti_CD3))) +
  theme_classic(base_size = 10)+ 
  theme(
    strip.background = element_blank(),  # if you dont want any background for strip label
    panel.background = element_rect(fill = "white", color = "lightgrey"),  # Set background color
    panel.grid.major = element_line(size = 0.5, color = "lightgrey"),  # Adjust grid line appearance
    panel.grid.minor = element_line(size = 0.25, color = "lightgrey"), # Adjust grid line appearance
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size =rel(1)), # the vjust and hjust only has the option as 0, 0.5, 1
    axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 1,size =rel(1)),
    axis.title.y = element_blank(),
    axis.title.x = element_text(hjust =0.5,face = "bold",margin = margin(0,10,0,0),size = rel(1.2),color = 'black'),
    legend.text = element_text(size = rel(0.6)), # the legend text size 
    legend.title = element_text(size = rel(0.9), face = "bold"),
    text = element_text(family="Times")) + # the legend title size 
  guides(fill = guide_legend(
    title = "FDR q value",
    keywidth = unit(0.6, "lines"), 
    keyheight = unit(0.6, "lines") )) + # specs to tweak the legend 
  labs(x= "something",
       y = "something",
       genotype = "samples"
       ) + 
  ggtitle(tbl3_CD107a_percen_batch$gating) 
  
  
  
  # some other ideas of setting the labs in ggplot 
  
  
  labs(
    color = "Genotype",  # Sets the legend title for color
    shape = "Genotype"   # Sets the legend title for shape
  )
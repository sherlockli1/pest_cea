library(ggplot2)
library(rgl)
library(ggrepel)

data<-read.csv("makeplot.csv")
colnames(data)

p1<-ggplot(data, aes(x = Added.PD.Cases.R, y = reorder(Pesticide.Name, Added.PD.Cases.R))) +
  geom_bar(stat = "identity") +
  geom_linerange(aes(xmin = lcl, xmax = ucl))+ 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 14))+
  labs(y = "Pesticides in Residence", x = "Added PD Cases in 20 years")
p1
p3<-ggplot(data, aes(x = Added.PD.Cases.C, y = reorder(Pesticide.Name, Added.PD.Cases.C))) +
  geom_bar(stat = "identity") +
  geom_linerange(aes(xmin = lcl.1, xmax = ucl.1))+ 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 14))+
  labs(y = "Pesticides in Workplace", x = "Added PD Cases in 20 years")
p3

png(filename="Manuscript/figure5_1.png",
    units="in", 
    width=16, 
    height=6, 
    res=300)
print(p1)
dev.off()

png(filename="Manuscript/figure5_3.png",
    units="in", 
    width=16, 
    height=6, 
    res=300)
print(p3)
dev.off()





exposure_res<-read.csv("Manuscript/df_cea_residential_new_psa_mean.csv")
exposure_occ<-read.csv("Manuscript/df_cea_occupational_new_psa_mean.csv")
exposure_res<-exposure_res[-which(exposure_res$v_names_str=="no pest"),]
exposure_occ<-exposure_occ[-which(exposure_occ$'v_names_str'=="no pest"),]


onset_res<-read.csv("Manuscript/df_cea_residential_new_psa_mean_onset.csv")
onset_occ<-read.csv("Manuscript/df_cea_occupational_new_psa_mean_onset.csv")
onset_res<-onset_res[-which(onset_res$'v_names_str'=="no pest"),]
onset_occ<-onset_occ[-which(onset_occ$'v_names_str'=="no pest"),]



p1 <- ggplot(exposure_res, aes(x=Inc_Effect, y=Inc_Cost)) +
  geom_point()+
  scale_size_continuous(range = c(1, 5))+
  labs(title="Incremental Cost Effectiveness Ratio (ICER) \n Compared to No Pesticide Exposure", 
       subtitle = "Residential Location",
       x="Incremental Quality-Adjusted Life-year (QALY)", y="Incremental Cost")+
  theme(panel.background = element_rect(fill = "white"),
        plot.title = element_text(color="black", size=10, face="bold",hjust = 0.5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_vline(xintercept = 0, linetype="dotted", color = "black")+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_text_repel(aes(label = Name), size = 4, color = "black")+ 
  scale_y_continuous(labels = scales::label_number_si())
p1

p2 <- ggplot(exposure_occ, aes(x=Inc_Effect, y=Inc_Cost)) +
  geom_point()+
  scale_size_continuous(range = c(1, 5))+
  labs(title="Incremental Cost Effectiveness Ratio (ICER) \n Compared to No Pesticide Exposure", 
       subtitle = "Workplace",
       x="Incremental Quality-Adjusted Life-year (QALY)", y="Incremental Cost")+
  theme(panel.background = element_rect(fill = "white"),
        plot.title = element_text(color="black", size=10, face="bold",hjust = 0.5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_vline(xintercept = 0, linetype="dotted", color = "black")+
  geom_hline(yintercept=0, linetype="dotted", color = "black")+
  geom_text_repel(aes(label = Name), size = 4, color = "black")+ 
  scale_y_continuous(labels = scales::label_number_si())
p2

png(filename="Manuscript/figure3_res.png",
    units="in", 
    width=10, 
    height=10, 
    res=300)
print(p1)
dev.off()

png(filename="Manuscript/figure3_occ.png",
    units="in", 
    width=10, 
    height=10, 
    res=300)
print(p2)
dev.off()


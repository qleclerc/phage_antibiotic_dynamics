
library(openxlsx)
library(ggplot2)
library(scales)
library(reshape2)
library(RColorBrewer)
library(dplyr)

clean_data = function(data, vol = 20){
  
  for(i in 2:(ncol(data)-1)){
    
    data[,i] = data[,i]*10^(data[,i+1])*(1000/vol)
    
  }
  
  data = data[,c(1,seq(2, ncol(data)-1, 2))]
  
  colnames(data) = gsub("_val", "", colnames(data))
  
  data = melt(data, id.vars = "Concentration")
  colnames(data) = c("Concentration", "time", "cfu")
  
  data$Concentration = data$Concentration
  data$time = as.numeric(as.character(data$time))
  
  data
  
}

palette = c("#00003f", rev(brewer.pal(n = 9, name = "RdBu"))[-5])

all_data = data.frame()

for(bac in c("201kt7", "327", "drp")){
  for(abx in c("ery", "tet")){
    
    filename = grep(abx,list.files(here::here("Data", bac)), value=T)
    
    data = c()
    
    for(i in filename){
      
      data_i = read.xlsx(here::here("Data", bac, i))
      data_i = clean_data(data_i)
      data = rbind(data, data_i)
      
    }
    
    data = data %>%
      group_by(Concentration, time) %>%
      summarise(se = sd(cfu),
                cfu = mean(cfu))
    
    data$abx = abx
    data$bac = bac
    
    all_data = rbind(all_data, data)
  }
}

bac_labs = c("NE201KT7", "NE327", "DRPET1")
names(bac_labs) = c("201kt7", "327", "drp")

abx_labs = c("Erythromycin", "Tetracycline")
names(abx_labs) = c("ery", "tet")

ggplot() +
  geom_line(data = all_data, aes(time, cfu, colour = as.factor(Concentration)), size=1) +
  geom_point(data = all_data, aes(time, cfu, colour = as.factor(Concentration)), size = 3) +
  geom_errorbar(data = all_data, aes(time, cfu, colour = as.factor(Concentration),
                                     ymin = pmax(0,cfu-se), ymax=cfu+se), width = 0.2, size = 1) +
  facet_grid(abx~bac, labeller = labeller(bac = bac_labs, abx = abx_labs)) +
  theme_bw() +
  labs(x = "Time (hours)", y = "cfu per mL", colour = "Antibiotic\nconcentration\n(mg/L):", linetype = "Source:") +
  scale_color_manual(values=palette)+
  coord_cartesian(ylim=c(1e1,1e10))+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))

ggsave(here::here("Figures","suppfig6.png"), dpi = 600)


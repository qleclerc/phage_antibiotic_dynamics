
# bac = "201kt7"
# abx = "tet"

filename = grep(abx,list.files(here::here("Data", bac)), value=T)

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


if(length(filename) == 1){
  
  data = read.xlsx(here::here("Data", bac, filename))
  data = clean_data(data)
  
  ggplot(data, aes(time, cfu, colour = as.factor(Concentration), group = as.factor(Concentration))) +
    geom_line(size=1) +
    geom_point(size = 3) +
    theme_bw() +
    labs(x = "Time (hours)", y = "cfu per mL", colour = "Antibiotic concentration:",
         title = paste0(bac, " ", abx)) +
    scale_color_manual(values=palette)+
    #scale_color_manual(values=palette(9))+
    coord_cartesian(ylim=c(1e1,1e10))+
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(breaks = seq(0,24,4)) +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.text = element_text(size=12),
          strip.text.x = element_text(size=12),
          legend.title = element_text(size=12))
  
} else {
  
  data = c()
  
  for(i in filename){
    
    data_i = read.xlsx(here::here("Data", bac, i))
    data_i = clean_data(data_i)
    data = rbind(data, data_i)
    
  }
  
  data = data %>%
    group_by(Concentration, time) %>%
    summarise(se = sd(cfu),#/sqrt(length(filename)),
              cfu = mean(cfu))
  
  ggplot(data, aes(time, cfu, colour = as.factor(Concentration), group = as.factor(Concentration))) +
    geom_line(size=1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = pmax(0,cfu-se), ymax=cfu+se)) +
    theme_bw() +
    labs(x = "Time (hours)", y = "cfu per mL", colour = "Antibiotic concentration (mg/L):",
         title = paste0(bac, " ", abx)) +
    scale_color_manual(values=palette)+
    #scale_color_manual(values=palette(9))+
    coord_cartesian(ylim=c(1e1,1e10))+
    scale_y_continuous(trans=log10_trans(),
                       breaks=trans_breaks("log10", function(x) 10^x),
                       labels=trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(breaks = seq(0,24,4)) +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.text = element_text(size=12),
          strip.text.x = element_text(size=12),
          legend.title = element_text(size=12))
  
}

write.csv(data, here::here("Data", paste0(bac, "_", abx, "_summary.csv")), row.names =F)
# ggsave(here::here("Figures", paste0(bac, "_", abx,".png")))


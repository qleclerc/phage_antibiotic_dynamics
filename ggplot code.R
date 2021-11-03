
library(ggplot2)
library(scales)
library(RColorBrewer)

accuracy_data = read.csv("accuracy_data.csv")
all_data = read.csv("all_data.csv")

#palette is built using RColorBrewer, but slightly tweaking the order to remove a very light colour in the middle
palette = c("#00003f", rev(brewer.pal(n = 9, name = "RdBu"))[-5])

bac_labs = c("NE201KT7", "NE327", "DRP")
names(bac_labs) = c("201kt7", "327", "drp")

abx_labs = c("Erythromycin", "Tetracycline")
names(abx_labs) = c("ery", "tet")

ggplot(all_data) +
  geom_line(aes(time, cfu, colour = as.factor(Concentration), linetype = source), size = 1) +
  geom_point(aes(time, cfu, colour = as.factor(Concentration)), size = 3) +
  geom_errorbar(aes(time, cfu, colour = as.factor(Concentration),
                    ymin = pmax(0,cfu-se), ymax=cfu+se), width = 0.2, size = 1) +
  facet_grid(abx~bac, labeller = labeller(bac = bac_labs, abx = abx_labs)) +
  theme_bw() +
  labs(x = "Time (hours)", y = "cfu per mL", colour = "Antibiotic\nconcentration\n(mg/L):", linetype = "Source:") +
  scale_color_manual(values = palette)+
  coord_cartesian(ylim = c(1e1,1e10))+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,2)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  geom_text(data = accuracy_data, mapping = aes(x = 1, y = 1e9, label = percent_in))

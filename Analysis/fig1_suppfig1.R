
library(dplyr)

all_params = data.frame()
all_data = data.frame()
all_effects = data.frame()

for(bac in c("327", "201kt7", "drp")){
  for(abx in c("ery", "tet")){
    
    cat("Working on", bac, abx, "\n")
    
    source(here::here("Analysis", "clean_plot_data.R"))
    source(here::here("Analysis", "estimate_parameters.R"))
    
    p = data.frame(bacteria = bac,
                   antibiotic = abx,
                   mu = growth_par[1],
                   Bmax = growth_par[2],
                   kmax = abx_par[1],
                   EC50 = abx_par[2],
                   pow = abx_par[3])
    
    all_params = rbind(all_params, p)
    
    model_results$bac = bac
    model_results$abx = abx
    all_data = rbind(all_data, model_results)
    
    summary_effect$bac = bac
    summary_effect$abx = abx
    all_effects = rbind(all_effects, summary_effect)
    
  }
}

write.csv(all_params, here::here("Parameters", "abx_params.csv"), row.names = F)

accuracy_data = all_data %>%
  dcast(., Concentration+time+bac+abx+se~source, value.var = "cfu") %>%
  group_by(Concentration, time, bac, abx) %>%
  summarise(se = sum(se, na.rm = T),
            Data = sum(Data, na.rm = T),
            model = sum(model, na.rm = T)) %>%
  mutate(in_interval = (model >= (Data-se) & model <= (Data+se))) %>%
  ungroup %>%
  group_by(bac, abx) %>%
  summarise(percent_in = paste0(round(sum(in_interval)/n()*100, 2), "%"))


bac_labs = c("NE201KT7", "NE327", "DRP")
names(bac_labs) = c("201kt7", "327", "drp")

abx_labs = c("Erythromycin", "Tetracycline")
names(abx_labs) = c("ery", "tet")

ggplot() +
  geom_line(data = all_data, aes(time, cfu, colour = as.factor(Concentration), linetype = source), size=1) +
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
  scale_x_continuous(breaks = seq(0,24,2)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  geom_text(data = accuracy_data, mapping = aes(x = 1, y = 1e9, label = percent_in))

ggsave(here::here("Figures","fig1.png"), dpi = 600)


ggplot(all_effects) +
  geom_line(aes(Concentration, effect, colour = "Data"), size = 1) +
  geom_point(aes(Concentration, effect, colour = "Data"), size = 3) +
  geom_line(aes(Concentration, fitted_effect, colour = "Fitted"), size = 1) +
  geom_point(aes(Concentration, fitted_effect, colour = "Fitted"), size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Concentration (mg/L)", y = "Kill rate per hour (relative to growth rate)", colour = "Source:") +
  facet_grid(abx~bac, labeller = labeller(bac = bac_labs, abx = abx_labs)) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

ggsave(here::here("Figures","suppfig1.png"), dpi = 600)

all_data = all_data %>%
  filter(source == "Data")

ggplot() +
  geom_line(data = all_data, aes(time, cfu, colour = as.factor(Concentration)), size=1) +
  geom_point(data = all_data, aes(time, cfu, colour = as.factor(Concentration)), size = 3) +
  geom_errorbar(data = all_data, aes(time, cfu, colour = as.factor(Concentration),
                                     ymin = pmax(0,cfu-se), ymax=cfu+se), width = 0.2, size = 1) +
  facet_grid(abx~bac, labeller = labeller(bac = bac_labs, abx = abx_labs)) +
  theme_bw() +
  labs(x = "Time (hours)", y = "cfu per mL", colour = "Antibiotic\nconcentration\n(mg/L):") +
  scale_color_manual(values=palette)+
  coord_cartesian(ylim=c(1e1,1e10))+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,2)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))

ggsave(here::here("Figures","timekill_curves.png"), dpi = 600)

ggplot(all_effects) +
  geom_line(aes(Concentration, effect, colour = "Data"), size = 1) +
  geom_point(aes(Concentration, effect, colour = "Data"), size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Concentration (mg/L)", y = "Kill rate per hour (relative to growth rate)") +
  facet_grid(abx~bac, labeller = labeller(bac = bac_labs, abx = abx_labs)) +
  theme_bw() +
  guides(colour = F) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

ggsave(here::here("Figures","abx_effect.png"), dpi = 600)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)
library(ggridges)
library(tidyr)
rm(list = setdiff(ls(), "project_root"))
setwd(file.path(project_root, "paper_figures/data/K562_bulk_WGS/"))
output_dir = file.path(project_root, "paper_figures/output/plot_mutations_bulk_WGA/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bulk_data = read.csv("bulk_mutations.csv")

head(bulk_data)

bulk_data = 
  bulk_data %>%
  mutate(Time = substring(Time, 2) %>% as.numeric()) %>% 
  dplyr::select(-Sample) %>%
  pivot_longer(    cols = contains("norm"),
                   names_to = "mutation",
                   values_to = "PMR_value")


ggplot(bulk_data) +
  geom_point(aes(x = Time,
                 y = PMR_value,
                 color = Template)) +
  facet_wrap(~mutation,scales = "free_y")

aavs_average_pmr =
  bulk_data %>% 
  group_by(mutation) %>%
  filter(Template == "aavs") %>%
  summarise(average_PMR = mean(PMR_value))
  
bulk_data %>% 
  left_join(aavs_average_pmr) %>% 
  filter(Template != "aavs") %>%
  filter(mutation == "mut_norm") %>%
  ggplot()+
  geom_point(aes(x = Time,
                 y = PMR_value)) +
  geom_smooth(method = "lm",
              aes(x = Time, 
                  y = PMR_value),se = F) +
    geom_hline(aes(yintercept  = average_PMR),
               linetype = "dotdash") +
  theme_classic() +
  theme(legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)) +
  scale_x_continuous(limits = c(0,15),
                     breaks = c(0,3,6,9,12,15),
                     labels =c("0","3","6","9","12","15")) +
  xlab("Days post editting") +
  ylab("Mutations PMR") 
ggsave(filename = paste(output_dir,"mutations_WGS_all.png",sep = ""),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 2.75)


bulk_data %>%
  left_join(aavs_average_pmr) %>% 
  filter(Template != "aavs") %>% 
  filter(mutation == "TCT.TTT_mut_PMR") %>%
  ggplot()+
  geom_point(aes(x = Time,
                 y = PMR_value)) +
  geom_smooth(method = "lm",
              aes(x = Time, 
                  y = PMR_value),se = F) +
  geom_hline(aes(yintercept  = average_PMR),
             linetype = "dotdash") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_x_continuous(limits = c(0,15),
                     breaks = c(0,3,6,9,12,15),
                     labels =c("0","3","6","9","12","15")) +
  xlab("Days post editting") +
  ylab("Mutations PMR") 
ggsave(filename = paste(output_dir,"mutations_WGS_TCT_TTT.png",sep = ""),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 2.75)

a = bulk_data %>% 
  left_join(aavs_average_pmr) %>% 
  filter(Template != "aavs") %>%
  filter(mutation == "TCT.TTT_norm") 
 
cor(a$mut,a$Time,method = "pearson")

bulk_data %>%
  left_join(aavs_average_pmr) %>% 
  filter(Template != "aavs") %>%
  filter(mutation == "mut_PMR") %>%
  ggplot()+
  geom_point(aes(x = Time,
                 y = PMR_value),
             color = "white") +
  geom_smooth(method = "lm",
              aes(x = Time, 
                  y = PMR_value),se = F) +
  geom_hline(aes(yintercept  = average_PMR),
             linetype = "dotdash",
             color = "#ff1a5e") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_x_continuous(limits = c(0,15),
                     breaks = c(0,3,6,9,12,15),
                     labels =c("0","3","6","9","12","15")) +
  xlab("Days post editting") +
  ylab("Mutations PMR") + 
  theme(
    text = element_text(color = "white"),                # All text elements white
    axis.text = element_text(color = "white", size = 8),          # Axis text white
    axis.title = element_text(color = "white", size = 10),         # Axis titles white
    axis.ticks = element_line(color = "white"),         # Axis ticks white
    axis.line = element_line(color = "white"),          # Axis lines white
    plot.background = element_blank(), #Black plot background
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.ontop = F,
    strip.background = element_blank(),
    strip.text = element_blank(),# Optional: Remove gridlines
    plot.margin = unit(c(1, 1, 1, 1), "cm"))
ggsave(file.path(output_dir, "mutations_WGS_all_legacy.png"),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 2.75)       
                



bulk_data %>%
  left_join(aavs_average_pmr) %>% 
  filter(Template != "aavs") %>% 
  filter(mutation == "TCT.TTT_mut_PMR") %>%
  ggplot()+
  geom_point(aes(x = Time,
                 y = PMR_value),
             color = "white") +
  geom_smooth(method = "lm",
              aes(x = Time, 
                  y = PMR_value),se = F) +
  geom_hline(aes(yintercept  = average_PMR),
             linetype = "dotdash",
             color = "#ff1a5e") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_x_continuous(limits = c(0,15),
                     breaks = c(0,3,6,9,12,15),
                     labels =c("0","3","6","9","12","15")) +
  xlab("Days post editting") +
  ylab("Mutations PMR") + 
  theme(
    text = element_text(color = "white"),                # All text elements white
    axis.text = element_text(color = "white", size = 8),          # Axis text white
    axis.title = element_text(color = "white", size = 10),         # Axis titles white
    axis.ticks = element_line(color = "white"),         # Axis ticks white
    axis.line = element_line(color = "white"),          # Axis lines white
    plot.background = element_blank(), #Black plot background
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.ontop = F,
    strip.background = element_blank(),
    strip.text = element_blank(),# Optional: Remove gridlines
    plot.margin = unit(c(1, 1, 1, 1), "cm"))
ggsave(file.path(output_dir, "mutations_WGS_TCT_TTT_legacy.png"),
       dpi= 600,bg = "transparent",
       height = 2.5, width = 2.75)

  




  
  
  
  
  
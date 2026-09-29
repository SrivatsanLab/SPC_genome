
library(dplyr)
library(ggplot2)
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(file.path(project_root, "paper_figures/data/DNA_protein_retention/Raw_protein_curves"))


getwd()



read.csv('nucleosome_signal.csv')-> nucleosome
read.csv('mono_nucleosome_signal.csv')->mono_nucleosome
read.csv("D4_signal.csv")-> D4
read.csv('Tetrahedron_signal.csv')-> tetrahedron
read.csv('Ribosome_signal.csv')-> ribosome
read.csv('Icosahedron_signal.csv')-> icosahedron
read.csv('D3_signal.csv')->D3 #in downloads not data
read.csv('Thyroglobulin_signal.csv')-> Thyroglobulin
read.csv('RNase_signal.csv')-> RNase
read.csv('Ovalbumin_signal.csv')-> Ovalbumin
read.csv(('Igg_lc_signal.csv'))->Igg_lc
read.csv('Igg_hc_signal.csv')-> Igg_hc
collapsed_mono <- mono_nucleosome %>%
  mutate(collapsed = ceiling(row_number() / 2)) %>%
  group_by(collapsed) %>%
  summarize_all(sum)
collapsed_mono[, -2] -> mono_nucleosome
time<- c(0,5,10,30,60)
as.data.frame(time)->time_df

#changing all the colnmaes so that a combined table is easier to make
colnames(Thyroglobulin)<- c('X','Thyroglobulin')
colnames(RNase)<- c('X','RNase')
colnames(Ovalbumin)<- c('X',"Ovalbumin")
colnames(D3)<- c('X','D3')
colnames(tetrahedron)<-  c("X","tetrahedron")
colnames(ribosome)<- c('X', 'ribosome')
colnames(icosahedron)<- c('X','icosahedron')
colnames(D4)<- c('X','D4')
colnames(Igg_hc)<- c('X','Igg_hc')
colnames(Igg_lc)<- c('X','Igg_lc')
colnames(mono_nucleosome)<- c('X','mono_nucleosome')
colnames(nucleosome) <- c('X','nucleosome')

#creating on combined df of all protein quants 
left_join(icosahedron,tetrahedron, by='X')->i_C_table
head(i_C_table)
left_join(i_C_table,ribosome, by='X')->irc
left_join(irc,D4,by='X')->irc_4
left_join(irc_4,D3,by='X')->irc_4_3

left_join_on_X <- function(df1, df2) {
  left_join(df1, df2, by = "X")
}

left_join_on_X(irc_4_3,Ovalbumin)->results
left_join_on_X(results,RNase)-> results
left_join_on_X(results,mono_nucleosome)-> results
left_join_on_X(results,nucleosome)-> results
left_join_on_X(results,Igg_hc)-> results
left_join_on_X(results,Igg_lc)-> results
left_join_on_X(results,Thyroglobulin)-> results
results-> protein_diffusion
head(protein_diffusion)
bind_cols(protein_diffusion,time_df)-> protein_diffusion_rate
results-> protein_diffusion
#calculating diffusion rates for all columns except for time column and X column 
protein_diffusion_excluded <- protein_diffusion[, !names(protein_diffusion) %in% "time"]
head(protein_diffusion_excluded)
ncol(protein_diffusion)
rate_change <- lapply(protein_diffusion[2:13], function(column) {
  diff(column) / diff(time)
})

head(rate_change)
 print(rate_change)
 as.data.frame(rate_change)-> rate_change_df
 print(rate_change_df)
 
 # Align time and rate of change
 aligned_time <- time_df[-1,]
as.data.frame(aligned_time)-> aligned_time_df
 # Drop the first time value to match lengths
 bind_cols(rate_change_df,aligned_time_df)->rate_change
 head(rate_change)
 
rate_change_df-> rate_change_protein
 
 # Scatter plot using ggplot2
 ggplot(rate_change, aes(x =aligned_time, y = log(max(rate_change)))) +
   geom_point() +
   labs(
     x = "Time",
     y = "rate of protein diffusion"
   ) +
   theme_minimal()
# Convert the list to a dataframe for easier manipulation if needed
rate_change_df <- as.data.frame(rate_change)
head(rate_change_protein)

print(protein_diffusion)
dir.create(file.path(project_root, "paper_figures/output/protein_signal"), recursive = TRUE, showWarnings = FALSE)
write.csv(protein_diffusion, file = file.path(project_root, "paper_figures/output/protein_signal/protein_diffusion.csv"))
write.csv(rate_change_protein, file = file.path(project_root, "paper_figures/output/protein_signal/rate_change_protein.csv"))


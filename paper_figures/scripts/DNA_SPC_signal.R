

library(dplyr)
library(ggplot2)
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(file.path(project_root, "paper_figures/data/DNA_protein_retention/Raw_DNA_curves"))


getwd()



read.csv('20bp_signal.csv')-> twenty
read.csv('35bp_signal.csv')-> thirty_five
read.csv('50bp_signal.csv')-> fifty 
read.csv('75bp_signal.csv')-> seventy_five
read.csv('100bp_signal.csv')-> one_hundred
read.csv('150bp_signal.csv')-> one_fifty
read.csv('200bp_signal.csv')-> two_hundred
read.csv('300bp_signal.csv')-> three_hundred
read.csv('400bp_signal.csv')-> four_hundred
read.csv('500bp_signal.csv')-> five_hundred
read.csv('650bp_signal.csv')-> six_fifty
read.csv('850bp_signal.csv')->eight_fifty
read.csv('1000bp_signal.csv')-> one_thousand
read.csv('1500bp_signal.csv')-> fifteen_hundred
read.csv('15kbp_signal.csv')-> fifteen_thousand
time<- c(0,5,10,15,20,25)

head(fifteen_thousand)
as.data.frame(time)->time_df

#changing all the colnmaes so that a combined table is easier to make
left_join_on_X <- function(df1, df2) {
  left_join(df1, df2, by = "X")
}
left_join_on_X(twenty, thirty_five)->combined
left_join_on_X(combined,fifty)-> combined
left_join_on_X(combined,seventy_five)-> combined
left_join_on_X(combined,one_hundred)-> combined
left_join_on_X(combined,one_fifty)-> combined
left_join_on_X(combined,two_hundred)-> combined
left_join_on_X(combined,three_hundred)-> combined
left_join_on_X(combined,four_hundred)->combined
left_join_on_X(combined, five_hundred)-> combined
left_join_on_X(combined,six_fifty)-> combined
left_join_on_X(combined,eight_fifty)-> combined
left_join_on_X(combined, one_thousand)-> combined
left_join_on_X(combined,fifteen_hundred)-> combined 
left_join_on_X(combined, fifteen_thousand)-> combined_bp

colnames(combined_bp)<- c('X','twenty','thirty_five','fifty','seventy_five','one_hundred',
                          'one_fifty','two_hundred','three_hundred','four_hundred',' five_hundred',
                          'six_fifty','eight_fifty','one_thousand','fifteen_hundred','fifteen_thousand')
head(combined_bp)
combined_bp-> PCR_ladder_diffusion
ncol(combined_bp)

rate_change_bp <- lapply(combined_bp[2:16], function(column) {
  diff(column) / diff(time)
})

head(rate_change_bp)
 
 as.data.frame(rate_change_bp)-> rate_change_bp
 print(rate_change_bp)
 
 # Align time and rate of change
 aligned_time <- time_df[-1,]
as.data.frame(aligned_time)-> aligned_time_df
 # Drop the first time value to match lengths
 bind_cols(rate_change_bp,aligned_time_df)->rate_change_bp
 head(rate_change_bp)
 

 
 # Scatter plot using ggplot2
 ggplot(rate_change_bp, aes(x = aligned_time, y = log(max(rate_change_bp)))) +
   geom_point() +
   labs(
     x = "Time",
     y = "rate of protein diffusion"
   ) +
   theme_minimal()
# Convert the list to a dataframe for easier manipulation if needed
rate_change_df <- as.data.frame(rate_change_bp)
head(rate_change_df)


dir.create(file.path(project_root, "paper_figures/output/dna_signal"), recursive = TRUE, showWarnings = FALSE)
write.csv(PCR_ladder_diffusion, file = file.path(project_root, "paper_figures/output/dna_signal/PCR_ladder_bp.csv"))
write.csv(rate_change_bp, file = file.path(project_root, "paper_figures/output/dna_signal/rate_change_bp.csv"))


#####################################################
# Name: exp6_analysis.R
# Description: Analyze experiment 6 results which
#              focus on showing the time taken for 
#              classifying reads using minimizer
#              indexes
# Date: June 28th, 2022
#####################################################

library(ggplot2)
library(ggpubr)
library(data.table)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################
exp2_data_paths <- c("/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/trial_3/exp2_index_stats.csv")

exp6_data_paths <- c("/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_6/trial_2/exp6_full_index_results.csv",
                     "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_6/trial_2/exp6_total_results.csv")

working_dir <- "/Users/omarahmed/downloads/current_research/spumoni_exps/exp_6/trial_2/plots/"

########################################################################
# Helper methods for generating plots
########################################################################

make_speedup_plot <- function(input_df) {
  plot <- ggplot(input_df, aes(x=w, y=totaltime)) + 
          geom_line(aes(group=type, linetype=type, color=type)) +
          geom_point(aes(group=type, color=type), size=2) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.title=element_text(size=12),
                legend.position= "none", #c(0.75,0.15),
                legend.background = element_rect(linetype="solid", color="black"),
                axis.text=element_text(size=12, color="black")) +
          scale_x_continuous(breaks=seq(0, 32, 2)) +
          scale_y_continuous(breaks=seq(0, 5.0, 0.5)) +
          scale_linetype_discrete(name="Alphabet", labels=c("DNA", "Minimizer")) +
          scale_color_discrete(name="Alphabet", labels=c("DNA", "Minimizer")) +
          #scale_shape_discrete(name="Dataset", labels=c("Ecoli_250", "Ecoli_500", "Salmonella_250", "Salmonella_500")) +
          labs(x="Window size (w)",
               y="Speed-up") 
        return(plot)
}

make_pml_index_plot <- function(input_df) {
  plot <- ggplot(subset(input_df, input_df$name == "ecoli_500"), aes(x=w, y=pml_size)) + 
          geom_line(aes(group=interaction(name,type), color=type, linetype=type))+
          geom_point(aes(group=interaction(name,type), color=type), size=2) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = c(0.25,0.15), 
                legend.background = element_rect(linetype="solid", color="black"),
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.title=element_text(size=12),
                axis.text=element_text(size=12, color="black")) +
          scale_x_continuous(breaks=seq(0, 32, 2)) +
          scale_y_continuous(breaks=seq(0, 1.5, 0.2), limits=c(0,NA)) +
          scale_color_discrete(name="Alphabet", labels=c("DNA", "Minimizer")) +
          scale_linetype_discrete(name="Alphabet", labels=c("DNA", "Minimizer")) +
          #scale_shape_discrete(name="Dataset", labels=c("Ecoli_250", "Ecoli_500", "Salmonella_250", "Salmonella_500")) +
          labs(x="Window size (w)",
               y="Relative Minimizer Index Size") 
        return(plot)
}

########################################################################
# Main method for the code
########################################################################

### Load in data for experiment 6 figure
exp6_full_index_data <- read.csv(exp6_data_paths[1], header=FALSE)
exp6_df <- read.csv(exp6_data_paths[2], header=TRUE)

# Determine total time needed, and use to divide to get speed-up
total_time_with_full_index <- as.double(exp6_full_index_data[,2])
exp6_df$totaltime <- total_time_with_full_index/exp6_df$totaltime

exp6_speedup_plot <- make_speedup_plot(exp6_df)
print(exp6_speedup_plot)
      
### Work on generating the index size plot from experiment 2
pos <- 1
datalist <- list()

for (input_file in exp2_data_paths) {
  datalist[[pos]] <- read.csv(input_file, header=TRUE)
  pos <- pos + 1
}
exp2_full_df <- do.call(rbind, datalist)

# Load in the statistics for full-sized indexes
full_stat_files <- c("/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/trial_3/exp2_full_index_stats.csv")

# Divide all the values by the full-index version
for (stat_file in full_stat_files) {
  curr_stats <- read.csv(stat_file, header=FALSE)
  
  # Read stats for full index
  dataset_name <- curr_stats[,1]
  n <- as.numeric(curr_stats[,2])
  r <- as.numeric(curr_stats[,3])
  n_over_r <- as.double(curr_stats[,4])
  ms_size <- as.numeric(curr_stats[,5])
  slp_size <- as.numeric(curr_stats[,6])
  pml_size <- as.numeric(curr_stats[,7])
  
  # Normalize all the minimizer index stats
  exp2_full_df$n[exp2_full_df$name == dataset_name] <- exp2_full_df$n[exp2_full_df$name == dataset_name]/n
  exp2_full_df$r[exp2_full_df$name == dataset_name] <- exp2_full_df$r[exp2_full_df$name == dataset_name]/r
  exp2_full_df$n_over_r[exp2_full_df$name == dataset_name] <- exp2_full_df$n_over_r[exp2_full_df$name == dataset_name]/n_over_r
  exp2_full_df$ms_size[exp2_full_df$name == dataset_name] <- exp2_full_df$ms_size[exp2_full_df$name == dataset_name]/ms_size
  exp2_full_df$slp_size[exp2_full_df$name == dataset_name] <- exp2_full_df$slp_size[exp2_full_df$name == dataset_name]/slp_size
  exp2_full_df$pml_size[exp2_full_df$name == dataset_name] <- exp2_full_df$pml_size[exp2_full_df$name == dataset_name]/pml_size
}

# Builds the plot showing the reduction in index size
exp2_index_plot <- make_pml_index_plot(exp2_full_df)
print(exp2_index_plot)

### Combine the two plots from above
final_plot <- ggarrange(exp2_index_plot, exp6_speedup_plot, ncol = 2, nrow = 1, labels=c("a","b") )
print(final_plot)

output_name <- paste(working_dir, "exp6_combined_plot.pdf", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="pdf", width=8, height=4)

output_name <- paste(working_dir, "exp6_combined_plot.jpeg", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="jpeg", width=8, height=4)

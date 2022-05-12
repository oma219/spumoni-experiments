#####################################################
# Name: expw_analysis.R
# Description: Generates various plots that show
#              the impact of minimizers on input
#              references when using SPUMONI
# Date: May 3rd, 2022
#####################################################

library(ggplot2)
library(data.table)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################
data_paths <- c("/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2a/exp2_index_stats.csv",
                "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2b/exp2_index_stats.csv",
                "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2c/exp2_index_stats.csv",
                "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2d/exp2_index_stats.csv",
                "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2e/exp2_index_stats.csv",
                "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2f/exp2_index_stats.csv")
working_dir <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/"


########################################################################
# Helper methods for generating plots
########################################################################

make_n_plot <- function(input_df) {
  plot <- ggplot(input_df, aes(x=w, y=n)) + 
          geom_line(aes(color=type, group=interaction(name, type)))+
          geom_point(aes(shape=name, color=type, group=interaction(name, type))) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=12, face="bold"),
                axis.title.x=element_text(size =12),
                axis.title.y=element_text(size=12),
                legend.position = "bottom", 
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.title=element_text(size=12),
                axis.text=element_text(size=12, color="black")) +
          scale_x_continuous(breaks=seq(0, 32, 2)) +
          scale_y_continuous(breaks=seq(0, 1.5, 0.2)) +
          scale_color_discrete(name="Minimizer Type", labels=c("DNA", "Promoted")) +
          scale_shape_discrete(name="Dataset", labels=c("Ecoli_10", "Ecoli_25", "Ecoli_50",
                                                        "Salmonella_10", "Salmonella_25", "Salmonella_50")) +
          labs(x="Window Size (w)",
               y="n - normalized",
               title="Size of Reference File After Different Types of Minimizer Digestion") 
  return(plot)
}

make_nr_plot <- function(input_df) {
  plot <- ggplot(input_df, aes(x=w, y=n_over_r)) + 
          geom_line(aes(group=interaction(name,type), color=type))+
          geom_point(aes(shape=name, color=type, group=interaction(name,type))) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=12, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.title=element_text(size=12),
                axis.text=element_text(size=12, color="black")) +
          scale_x_continuous(breaks=seq(0, 32, 2)) +
          scale_y_continuous(breaks=seq(0, 1.5, 0.2)) +
          scale_color_discrete(name="Minimizer Type", labels=c("DNA", "Promoted")) +
          scale_shape_discrete(name="Dataset", labels=c("Ecoli_10", "Ecoli_25", "Ecoli_50",
                                                        "Salmonella_10", "Salmonella_25", "Salmonella_50")) +
          labs(x="Window Size (w)",
               y="n/r - normalized",
               title="Average Run Length (n/r) for Different Datasets After Minimizer Digestion") 
  return(plot)
}

make_ms_index_plot <- function(input_df) {
  plot <- ggplot(input_df, aes(x=w, y=ms_size)) + 
          geom_line(aes(group=interaction(name,type), color=type))+
          geom_point(aes(shape=name, color=type, group=interaction(name,type))) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=12, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.title=element_text(size=12),
                axis.text=element_text(size=12, color="black")) +
          scale_x_continuous(breaks=seq(0, 32, 2)) +
          scale_y_continuous(breaks=seq(0, 1.5, 0.2)) +
          scale_color_discrete(name="Minimizer Type", labels=c("DNA", "Promoted")) +
          scale_shape_discrete(name="Dataset", labels=c("Ecoli_10", "Ecoli_25", "Ecoli_50",
                                                        "Salmonella_10", "Salmonella_25", "Salmonella_50")) +
          labs(x="Window Size (w)",
               y="MS Index Size - normalized",
               title="Normalized MS Index Size for Different Datasets After Minimizer Digestion") 
  return(plot)
}

make_pml_index_plot <- function(input_df) {
  plot <- ggplot(input_df, aes(x=w, y=pml_size)) + 
          geom_line(aes(group=interaction(name,type), color=type))+
          geom_point(aes(shape=name, color=type, group=interaction(name,type))) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=12, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.title=element_text(size=12),
                axis.text=element_text(size=12, color="black")) +
          scale_x_continuous(breaks=seq(0, 32, 2)) +
          scale_y_continuous(breaks=seq(0, 1.5, 0.2)) +
          scale_color_discrete(name="Minimizer Type", labels=c("DNA", "Promoted")) +
          scale_shape_discrete(name="Dataset", labels=c("Ecoli_10", "Ecoli_25", "Ecoli_50",
                                                        "Salmonella_10", "Salmonella_25", "Salmonella_50")) +
          labs(x="Window Size (w)",
               y="PML Index Size - normalized",
               title="Normalized PML Index Size for Different Datasets After Minimizer Digestion") 
  return(plot)
}


#########################################################################
# Start of the "main" method of code ...
#########################################################################

# Load in the data at all the different window sizes
pos <- 1
datalist <- list()

for (input_file in data_paths) {
  datalist[[pos]] <- read.csv(input_file, header=TRUE)
  pos <- pos + 1
}
full_df <- do.call(rbind, datalist)

# Load in the statistics for full-sized indexes
full_stat_files <- c("/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2a/exp2_full_index_stats.csv",
                     "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2b/exp2_full_index_stats.csv",
                     "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2c/exp2_full_index_stats.csv",
                     "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2d/exp2_full_index_stats.csv",
                     "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2e/exp2_full_index_stats.csv",
                     "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_2/exp_2f/exp2_full_index_stats.csv")

for (stat_file in full_stat_files) {
  curr_stats <- read.csv(stat_file, header=FALSE)
  
  # Read stats for full index
  dataset_name <- curr_stats[,1]
  n <- as.integer(curr_stats[,2])
  r <- as.integer(curr_stats[,3])
  n_over_r <- as.double(curr_stats[,4])
  ms_size <- as.integer(curr_stats[,5])
  slp_size <- as.integer(curr_stats[,6])
  pml_size <- as.integer(curr_stats[,7])
  
  # Normalize all the minimizer index stats
  full_df$n[full_df$name == dataset_name] <- full_df$n[full_df$name == dataset_name]/n
  full_df$r[full_df$name == dataset_name] <- full_df$r[full_df$name == dataset_name]/r
  full_df$n_over_r[full_df$name == dataset_name] <- full_df$n_over_r[full_df$name == dataset_name]/n_over_r
  full_df$ms_size[full_df$name == dataset_name] <- full_df$ms_size[full_df$name == dataset_name]/ms_size
  full_df$slp_size[full_df$name == dataset_name] <- full_df$slp_size[full_df$name == dataset_name]/slp_size
  full_df$pml_size[full_df$name == dataset_name] <- full_df$pml_size[full_df$name == dataset_name]/pml_size
}

# Generate the plots, and save them

# Plot #1
n_plot <- make_n_plot(full_df)
n_plot
output_name <- paste(working_dir, "exp2_normalized_n.png", sep="")
ggsave(output_name, plot=n_plot, dpi=1200, device="jpeg", width=8, height=6)

# Plot #2
nr_plot <- make_nr_plot(full_df)
nr_plot
output_name <- paste(working_dir, "exp2_normalized_n_over_r.png", sep="")
ggsave(output_name, plot=nr_plot, dpi=1200, device="jpeg", width=8, height=6)

# Plot #3
ms_plot <- make_ms_index_plot(full_df)
ms_plot
output_name <- paste(working_dir, "exp2_normalized_ms_size.png", sep="")
ggsave(output_name, plot=ms_plot, dpi=1200, device="jpeg", width=8, height=6)

# Plot #4
pml_plot <- make_pml_index_plot(full_df)
pml_plot
output_name <- paste(working_dir, "exp2_normalized_pml_size.png", sep="")
ggsave(output_name, plot=pml_plot, dpi=1200, device="jpeg", width=8, height=6)

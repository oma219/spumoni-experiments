#####################################################
# Name: exp3_ks_stat_analysis.R
# Description: Generates a dot plot that analyzes
#              the ks-stats that we see from different
#              classes of reads.
# Date: May 15th, 2022
#####################################################

library(ggplot2)
library(data.table)
library(ggpubr)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################

pos_class_reads <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_3c_ksstat_analysis/dna_index/long_positive_reads.fa.ks_stats"
null_class_reads <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_3c_ksstat_analysis/dna_index/long_null_reads.fa.ks_stats"
working_dir <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_3c_ksstat_analysis/plots/dna/"

read_length <- "long"
index_type <- "dna"

########################################################################
# Helper methods for generating plots
########################################################################

make_ks_stat_plot <- function(total_df) {
    # find the threshold
    threshold <- total_df[1,3]

    # create the plot
    plot <- ggplot(total_df, aes(x=ks_stat, y=pos)) + 
      geom_point(aes(color=class, group=class), position=position_jitter(h=0.01,w=0.01), alpha=0.4) +
      theme_bw() +
      theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
            axis.title.x=element_text(size =14),
            axis.title.y=element_text(size=14),
            legend.position = "bottom", 
            legend.text=element_text(size=12),
            legend.box="vertical",
            legend.title=element_text(size=12),
            axis.text=element_text(size=12, color="black")) +
      scale_x_continuous(breaks=seq(0, 1, 0.1)) +
      scale_y_continuous(breaks=seq(0.0, 1.0, 0.1)) +
      scale_color_discrete(name="Read Class", labels=c("Null", "Positive")) +
      geom_vline(xintercept=threshold, linetype="dashed", color = "black", size=1.5) +
      labs(x="KS-statistic",
           y="Position in the Read",
           title="KS-statistics across positive and null long reads using DNA minimizers") 
    return(plot)
}

make_ks_stat_histo_plot <- function(total_df) {
  # find the threshold
  threshold <- total_df[1,3]
  
  # create the plot
  plot <- ggplot(total_df, aes(x=ks_stat)) + 
    geom_histogram(aes(color=class, fill=class, group=class, y=stat(density)), binwidth=0.01, alpha=0.5) +
    geom_density(aes(color=class), size = 1) +
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.position = "bottom", 
          legend.text=element_text(size=12),
          legend.box="vertical",
          legend.title=element_text(size=12),
          axis.text=element_text(size=12, color="black")) +
    geom_vline(xintercept=threshold, linetype="dashed", color = "black", size=1.5) +
    scale_x_continuous(breaks=seq(0, 1, 0.1)) +
    scale_color_discrete(name="Read Class", labels=c("Null", "Positive")) +
    scale_fill_discrete(name="Read Class", labels=c("Null", "Positive")) +
    #geom_vline(xintercept=0.10, linetype="dashed", color = "black", size=1.5) +
    labs(x="KS-statistic",
         y="Density",
         title="Density of KS-statistics for long reads using DNA minimizers") 
  return(plot)
}


#########################################################################
# Start of the "main" method of code ...
#########################################################################

# read in the data for positive/null reads
pos_reads_df <- read.csv(pos_class_reads, header=FALSE, col.names=c("pos", "ks_stat", "threshold"))
null_reads_df <- read.csv(null_class_reads, header=FALSE, col.names=c("pos", "ks_stat", "threshold"))

# sub-sample in the case for long reads, too many data points
pos_reads_df <- pos_reads_df[sample(nrow(pos_reads_df), 100000), ]
null_reads_df <- null_reads_df[sample(nrow(null_reads_df), 100000), ]

# add in attribute identifying the class
pos_reads_df["class"] <- "pos"
null_reads_df["class"] <- "null"

total_df <- rbind(pos_reads_df, null_reads_df)

# Plot 1: Generate plot that looks at KS-statistics for positive/null reads
ks_stat_plot <- make_ks_stat_plot(total_df)
ks_stat_plot

output_name <- paste(working_dir, "exp3_", index_type, "_", read_length, "_ks_stat_points.pdf", sep="")
ggsave(output_name, plot=ks_stat_plot, dpi=800, device="pdf", width=8, height=6)

output_name <- paste(working_dir, "exp3_", index_type, "_", read_length, "_ks_stat_points.jpeg", sep="")
ggsave(output_name, plot=ks_stat_plot, dpi=800, device="jpeg", width=8, height=6)

# Plot 2: Generate histogram to look at the percentage in each bin
ks_stat_histo_plot <- make_ks_stat_histo_plot(total_df)
ks_stat_histo_plot

output_name <- paste(working_dir, "exp3_", index_type, "_", read_length, "_ks_stat_density.pdf", sep="")
ggsave(output_name, plot=ks_stat_histo_plot, dpi=800, device="pdf", width=8, height=6)

output_name <- paste(working_dir, "exp3_", index_type, "_", read_length, "_ks_stat_density.jpeg", sep="")
ggsave(output_name, plot=ks_stat_histo_plot, dpi=800, device="jpeg", width=8, height=6)


# Generate a combined plot
combined_plots <- ggarrange(ks_stat_plot, ks_stat_histo_plot, 
                  labels = c("A", "B"),
                  ncol = 2, nrow = 1)

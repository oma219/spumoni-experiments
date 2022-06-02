#####################################################
# Name: exp3_analysis.R
# Description: Generates various plots that show
#              the shows the binary classification
#              accuracy with different minimizer schemes.
# Date: May 15th, 2022
#####################################################

library(ggplot2)
library(data.table)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################

all_indexes_data_path <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_3d/data/exp3_total_results.csv"
full_index_data_path <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_3d/data/exp3_full_index_results.csv"
working_dir <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_3d/"

########################################################################
# Helper methods for generating plots
########################################################################

make_promotion_acc_plot <- function(input_df, full_index_df) {
    # subset the data as needed
    newdata_df <- subset(input_df, indextype == "promotion")
    
    # extract some statistics on the full FASTA index
    temp_df <- subset(full_index_df, readlength == "long")
    long_read_acc <- temp_df[,"accuracy"]
    temp_df <- subset(full_index_df, readlength == "short")
    short_read_acc <- temp_df[,"accuracy"]

    # create the plot
    plot <- ggplot(newdata_df, aes(x=w, y=accuracy)) + 
      geom_line(aes(color=readlength, group=readlength))+
      geom_point(aes(color=readlength)) +
      theme_bw() +
      theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
            axis.title.x=element_text(size =14),
            axis.title.y=element_text(size=14),
            legend.position = "bottom", 
            legend.text=element_text(size=12),
            legend.box="vertical",
            legend.title=element_text(size=12),
            axis.text=element_text(size=12, color="black")) +
      scale_x_continuous(breaks=seq(6, 32, 2)) +
      scale_y_continuous(breaks=seq(0.60, 1.0, 0.05)) +
      scale_color_discrete(name="Read Length", labels=c("Long", "Short")) +
      geom_hline(yintercept=long_read_acc, linetype="dashed", color="#F8766D", size=0.75) +
      geom_hline(yintercept=short_read_acc, linetype="dashed", color="#00BFC4", size=0.75) +
      labs(x="Window size (w)",
           y="Accuracy",
           title="Binary read classification using promoted minimizer schemes with SPUMONI") 
    return(plot)
}

make_dna_acc_plot <- function(input_df, full_index_df) {
    # subset the data as needed
    newdata_df <- subset(input_df, indextype == "dna")
    
    # extract some statistics on the full FASTA index
    temp_df <- subset(full_index_df, readlength == "long")
    long_read_acc <- temp_df[,"accuracy"]
    temp_df <- subset(full_index_df, readlength == "short")
    short_read_acc <- temp_df[,"accuracy"]
    
    # create the plot
    plot <- ggplot(newdata_df, aes(x=w, y=accuracy)) + 
      geom_line(aes(color=readlength, group=readlength))+
      geom_point(aes(color=readlength)) +
      theme_bw() +
      theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
            axis.title.x=element_text(size =14),
            axis.title.y=element_text(size=14),
            legend.position = "bottom", 
            legend.text=element_text(size=12),
            legend.box="vertical",
            legend.title=element_text(size=12),
            axis.text=element_text(size=12, color="black")) +
      scale_x_continuous(breaks=seq(6, 32, 2)) +
      scale_y_continuous(breaks=seq(0.70, 1.0, 0.05)) +
      scale_color_discrete(name="Read Length", labels=c("Long", "Short")) +
      geom_hline(yintercept=long_read_acc, linetype="dashed", color="#F8766D", size=0.75) +
      geom_hline(yintercept=short_read_acc, linetype="dashed", color="#00BFC4", size=0.75) +
      labs(x="Window size (w)",
           y="Accuracy",
           title="Binary read classification using DNA minimizer schemes with SPUMONI") 
    return(plot)
}

make_accuracy_index_size_plot <- function(total_dataset_df) {
  # create the plot
  plot <- ggplot(total_dataset_df, aes(x=pmlindexsize, y=accuracy)) + 
          geom_line(aes(color=readlength, group=interaction(indextype, readlength)))+
          geom_point(aes(color=readlength, shape=indextype)) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.title=element_text(size=12),
                axis.text=element_text(size=12, color="black")) +
          scale_x_continuous(breaks=seq(8000000, 80000000, 8000000)) +
          scale_y_continuous(breaks=seq(0.60, 1.0, 0.05)) +
          scale_color_discrete(name="Read Length", labels=c("Long", "Short")) +
          scale_shape_discrete(name="Minimizer Type", labels=c("DNA", "Promotion")) +
          #geom_hline(yintercept=long_read_acc, linetype="dashed", color="#F8766D", size=0.75) +
          #geom_hline(yintercept=short_read_acc, linetype="dashed", color="#00BFC4", size=0.75) +
          labs(x="Index Size (bytes)",
               y="Accuracy",
               title="Binary read classification for different minimizer schemes and index sizes") 
  return(plot)
}

#########################################################################
# Start of the "main" method of code ...
#########################################################################

total_dataset_df <- read.csv(all_indexes_data_path, header=TRUE)
full_index_df <- read.csv(full_index_data_path, header=TRUE)

# Generate the plots, and save them

# Plot #1
promotion_acc_plot <- make_promotion_acc_plot(total_dataset_df, full_index_df)
promotion_acc_plot

output_name <- paste(working_dir, "exp3_promotion_min_accuracy.pdf", sep="")
ggsave(output_name, plot=promotion_acc_plot, dpi=1200, device="pdf", width=8, height=6)

output_name <- paste(working_dir, "exp3_promotion_min_accuracy.jpeg", sep="")
ggsave(output_name, plot=promotion_acc_plot, dpi=800, device="jpeg", width=8, height=6)

# Plot #2
dna_acc_plot <- make_dna_acc_plot(total_dataset_df, full_index_df)
dna_acc_plot

output_name <- paste(working_dir, "exp3_dna_min_accuracy.pdf", sep="")
ggsave(output_name, plot=dna_acc_plot, dpi=1200, device="pdf", width=8, height=6)

output_name <- paste(working_dir, "exp3_dna_min_accuracy.jpeg", sep="")
ggsave(output_name, plot=dna_acc_plot, dpi=800, device="jpeg", width=8, height=6)

# Plot #3
acc_index_plot <- make_accuracy_index_size_plot(total_dataset_df)
acc_index_plot

output_name <- paste(working_dir, "exp3_accuracy_index_size.pdf", sep="")
ggsave(output_name, plot=acc_index_plot, dpi=1200, device="pdf", width=8, height=6)

output_name <- paste(working_dir, "exp3_accuracy_index_size.jpeg", sep="")
ggsave(output_name, plot=acc_index_plot, dpi=800, device="jpeg", width=8, height=6)

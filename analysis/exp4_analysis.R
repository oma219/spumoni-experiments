#####################################################
# Name: exp4_analysis.R
# Description: Generates various plots that show
#              the shows the binary classification
#              accuracy with bin sizes for KS-test.
# Date: June 4th, 2022
#####################################################

library(ggplot2)
library(data.table)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################

all_indexes_data_path <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_4/data/exp4_total_results.csv"
working_dir <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_4/"

########################################################################
# Helper methods for generating plots
########################################################################

make_promotion_short_long<- function(input_df) {
    # subset the data as needed
    newdata_df <- subset(input_df, indextype == "promotion")
    newdata_df$binsize <- factor(newdata_df$binsize)

    # create the plot
    plot <- ggplot(newdata_df, aes(x=pmlindexsize, y=accuracy)) + 
      geom_line(aes(color=binsize, group=interaction(binsize, readlength)))+
      geom_point(aes(color=binsize, shape=readlength, group=interaction(binsize, readlength))) +
      theme_bw() +
      theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
            axis.title.x=element_text(size =14),
            axis.title.y=element_text(size=14),
            legend.position = "bottom", 
            legend.text=element_text(size=12),
            legend.box="vertical",
            legend.title=element_text(size=12),
            axis.text=element_text(size=12, color="black")) +
      scale_y_continuous(breaks=seq(0.80, 1.0, 0.05)) +
      scale_x_continuous(breaks=seq(15000000, 45000000, 10000000)) +
      scale_color_discrete(name="Region Size") +
      scale_shape_discrete(name="Read Length", labels=c("Long", "Short")) +
      labs(x="Index Size (bytes)",
           y="Accuracy",
           title="Binary read classification using promoted minimizers with different region sizes") 
    return(plot)
}

make_dna_short_long<- function(input_df) {
  # subset the data as needed
  newdata_df <- subset(input_df, indextype == "dna")
  newdata_df$binsize <- factor(newdata_df$binsize)
  
  # create the plot
  plot <- ggplot(newdata_df, aes(x=pmlindexsize, y=accuracy)) + 
          geom_line(aes(color=binsize, group=interaction(binsize, readlength)))+
          geom_point(aes(color=binsize, shape=readlength, group=interaction(binsize, readlength))) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.title=element_text(size=12),
                axis.text=element_text(size=12, color="black")) +
          scale_y_continuous(breaks=seq(0.80, 1.0, 0.05)) +
          scale_x_continuous(breaks=seq(25000000, 75000000, 10000000)) +
          scale_color_discrete(name="Region Size") +
          scale_shape_discrete(name="Read Length", labels=c("Long", "Short")) +
          labs(x="Index Size (bytes)",
               y="Accuracy",
               title="Binary read classification using DNA minimizers with different region sizes") 
  return(plot)
}

make_promotion_short_dna_short <- function(input_df) {
  # subset the data as needed
  newdata_df <- subset(input_df, readlength == "short")
  newdata_df$binsize <- factor(newdata_df$binsize)
  
  # create the plot
  plot <- ggplot(newdata_df, aes(x=pmlindexsize, y=accuracy)) + 
    geom_line(aes(color=binsize, group=interaction(binsize, indextype)))+
    geom_point(aes(color=binsize, shape=indextype, group=interaction(binsize, indextype))) +
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5, size=12, face="bold"),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.position = "bottom", 
          legend.text=element_text(size=12),
          legend.box="vertical",
          legend.title=element_text(size=12),
          axis.text=element_text(size=12, color="black")) +
    scale_y_continuous(breaks=seq(0.80, 1.0, 0.05)) +
    scale_x_continuous(breaks=seq(15000000, 75000000, 10000000)) +
    scale_color_discrete(name="Region Size") +
    scale_shape_discrete(name="Minimizer Type", labels=c("DNA", "Promoted")) +
    labs(x="Index Size (bytes)",
         y="Accuracy",
         title="Binary classification of short reads using different minimizer schemes with different region sizes") 
  return(plot)
}

make_promotion_long_dna_long <- function(input_df) {
  # subset the data as needed
  newdata_df <- subset(input_df, readlength == "long")
  newdata_df$binsize <- factor(newdata_df$binsize)
  
  # create the plot
  plot <- ggplot(newdata_df, aes(x=pmlindexsize, y=accuracy)) + 
    geom_line(aes(color=binsize, group=interaction(binsize, indextype)))+
    geom_point(aes(color=binsize, shape=indextype, group=interaction(binsize, indextype))) +
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5, size=12, face="bold"),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.position = "bottom", 
          legend.text=element_text(size=12),
          legend.box="vertical",
          legend.title=element_text(size=12),
          axis.text=element_text(size=12, color="black")) +
    scale_y_continuous(breaks=seq(0.80, 1.0, 0.05)) +
    scale_x_continuous(breaks=seq(15000000, 75000000, 10000000)) +
    scale_color_discrete(name="Region Size") +
    scale_shape_discrete(name="Minimizer Type", labels=c("DNA", "Promoted")) +
    labs(x="Index Size (bytes)",
         y="Accuracy",
         title="Binary classification of long reads using different minimizer schemes with different region sizes") 
  return(plot)
}

#########################################################################
# Start of the "main" method of code ...
#########################################################################

### Load the data into a dataframe
total_dataset_df <- read.csv(all_indexes_data_path, header=TRUE)


### Generate the plots, and save them

# Plot #1
plot_one <- make_promotion_short_long(total_dataset_df)
plot_one

output_name <- paste(working_dir, "exp4_promotion_long_short.pdf", sep="")
ggsave(output_name, plot=plot_one, dpi=1200, device="pdf", width=8.5, height=6)

output_name <- paste(working_dir, "exp4_promotion_long_short.jpeg", sep="")
ggsave(output_name, plot=plot_one, dpi=800, device="jpeg", width=8.5, height=6)

# Plot #2
plot_two <- make_dna_short_long(total_dataset_df)
plot_two

output_name <- paste(working_dir, "exp4_dna_long_short.pdf", sep="")
ggsave(output_name, plot=plot_two, dpi=1200, device="pdf", width=8.5, height=6)

output_name <- paste(working_dir, "exp4_dna_long_short.jpeg", sep="")
ggsave(output_name, plot=plot_two, dpi=800, device="jpeg", width=8.5, height=6)

# Plot #3
plot_three <- make_promotion_short_dna_short(total_dataset_df)
plot_three

output_name <- paste(working_dir, "exp4_promotion_short_dna_short.pdf", sep="")
ggsave(output_name, plot=plot_three, dpi=1200, device="pdf", width=8.5, height=6)

output_name <- paste(working_dir, "exp4_promotion_short_dna_short.jpeg", sep="")
ggsave(output_name, plot=plot_three, dpi=800, device="jpeg", width=8.5, height=6)

# Plot #4
plot_four <- make_promotion_long_dna_long(total_dataset_df)
plot_four

output_name <- paste(working_dir, "exp4_promotion_long_dna_long.pdf", sep="")
ggsave(output_name, plot=plot_four, dpi=1200, device="pdf", width=8.5, height=6)

output_name <- paste(working_dir, "exp4_promotion_long_dna_long.jpeg", sep="")
ggsave(output_name, plot=plot_four, dpi=800, device="jpeg", width=8.5, height=6)


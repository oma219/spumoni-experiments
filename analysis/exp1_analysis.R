#####################################################
# Name: exp1_analysis.R
# Description: Generates grouped bar charts for 
#              experiment 1 regarding the document
#              array.
# Date: February 1, 2022
#####################################################

library(ggplot2)
library(data.table)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################
data_paths <- c("/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_1/trial_4/data/illumina_ms_doc_analysis.csv",
                "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_1/trial_4/data/ont_ms_doc_analysis.csv")
output_dir <- "/Users/omarahmed/downloads/current_research/spumoni_exps/exp_1/trial_4/plots/"

dataset_names <- c("E. coli", "S. enterica", "L. monocytogenes", "P. aeruginosa", "B. subtilis", "L. fermentum", "E. faecalis", "S. aureus")
x_labels <- c("Random Reads", "E. coli", "S. enterica", "L. monocytogenes", "P. aeruginosa", "B. subtilis", "L. fermentum", "E. faecalis", "S. aureus")

columns_to_extract <- c("class_1_percent", "class_2_percent", "class_3_percent",
                        "class_4_percent", "class_5_percent", "class_6_percent", "class_7_percent",
                        "class_8_percent")


########################################################################
# Methods for generating grouped bar charts
########################################################################

make_group_bar_plot <- function(input_df, read_type) {
  # Convert to data.table, melt, and then revert to data.frame
  setDT(input_df)
  group_df_melt <- melt(input_df, 
                        measure.vars = columns_to_extract,
                        variable.name = "class", 
                        value.name = "percent")
  setDF(group_df_melt)
  
  cbbPalette <- c( "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73")
  
  # Create the grouped bar-chart
  #read_type <- "Illumina"
  graph_title <- paste("Multi-Class Classification of", read_type, "Reads Using the Document Array")
  plot <- ggplot(group_df_melt, aes(fill=class, x=read_set, y=percent)) + 
          geom_bar(position="dodge", stat="identity")+
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.text.x=element_text(face="italic", size=9),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=12, face="italic"),
                legend.box="vertical",
                legend.title=element_text(size=12),
                axis.text.y=element_text(size=12, color="black")) +
          scale_y_continuous(breaks=seq(0, 1.0, 0.1)) +
          scale_x_discrete(labels=x_labels) +
          labs(x="Simulated Read Set",
               y="Average Ratio of Document Labels Across Read",
               title="") +
          #scale_fill_discrete(name="Species", labels=dataset_names) +
          scale_fill_manual(name="Species", labels=dataset_names, values=cbbPalette)
    return (plot)
}

#########################################################################
# Start of the "main" method of code ...
#########################################################################

# Generates two plots corresponding to short and long reads
pos <- 1
read_types <- c("Illumina", "ONT")

for (input_file in data_paths) {
  input_df <- read.csv(input_file, header=TRUE)
  curr_plot <- make_group_bar_plot(input_df, read_types[pos])
  
  # Save a vector graphic & regular graphic
  output_name <- paste(output_dir, "exp1_plot_", read_types[pos], "_docarray_analysis.pdf", sep="")
  ggsave(output_name, plot=curr_plot, dpi=1200, device="pdf", width=10, height=6)

  output_name <- paste(output_dir, "exp1_plot_", read_types[pos], "_docarray_analysis.jpeg", sep="")
  ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=10, height=6)
  pos <- pos + 1
}


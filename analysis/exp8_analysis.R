#####################################################
# Name: exp8_analysis.R
# Description: Analyze experiment 8 involving looking
#              MS from contaminated assemblies.
#
# Date: July 24th, 2022
#####################################################

library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)

########################################################################
# IMPORTANT: Data paths
########################################################################
data_working_dir <- "/Users/omarahmed/downloads/current_research/spumoni_exps/exp_8/trial_1/data/assembly_1/"
output_dir <- "/Users/omarahmed/downloads/current_research/spumoni_exps/exp_8/trial_1/plots/assembly_1/"

########################################################################
# Helper methods for generating plots
########################################################################

make_length_plot <- function(contig_name, input_df, line_color) {
  x_label <- paste("Position in ", contig_name)
  plot <- ggplot(data=input_df) +
          geom_line(aes(x=pos, y=length), color=line_color, size=0.5) +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =12),
                axis.title.y=element_text(size=12),
                legend.text=element_text(size=12),
                legend.box="horizontal",
                legend.title=element_text(size=12),
                legend.position="bottom",
                axis.text=element_text(size=10, color="black")) +
          #scale_y_continuous(breaks=seq(0, 7000, 1000)) +
          #scale_x_continuous(breaks=seq(0, 40000, 5000)) +
          labs(y="Matching Statistic Length", x=x_label)
  return(plot)         
}

make_single_boxplot <- function(input_df) {
  plot <- ggplot(subset(input_df, name=="All other contigs"), aes(x=name, y=length, fill=status)) +
          geom_boxplot(outlier.shape=4) +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =10),
                axis.title.y=element_text(size=10),
                legend.text=element_text(size=12),
                legend.box="horizontal",
                legend.title=element_text(size=12),
                legend.position="none",
                axis.text=element_text(size=8, color="black"),
                plot.margin = unit(c(0,0,0,0), "cm")) +
        scale_fill_discrete(name="", labels=c("Normal", "Suspicious")) +
        scale_y_continuous(breaks=seq(0, 30, 5)) +
        labs(y="", x="")
  return (plot)
}

make_comparison_plot <- function(input_df) {
  plot <- ggplot(input_df, aes(x=name, y=length, fill=status)) + 
          geom_boxplot(outlier.shape=4) +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =12),
                axis.title.y=element_text(size=12),
                legend.text=element_text(size=12),
                axis.text.x = element_text(angle = 0, vjust=0.5, size=8),
                legend.box="horizontal",
                legend.title=element_text(size=12),
                legend.position=c(0.8, 0.9),
                axis.text.y =element_text(size=10, color="black")) +
        #scale_x_continuous(breaks=seq(0, 32, 2)) +
        #scale_y_continuous(breaks=seq(0, 600, 100)) +
        #scale_linetype_discrete(name="Minimizer Type", labels=c("DNA", "Promoted")) +
        scale_fill_discrete(name="", labels=c("Normal", "Suspicious")) +
        #scale_shape_discrete(name="Dataset", labels=c("Ecoli_250", "Ecoli_500", "Salmonella_250", "Salmonella_500")) +
        labs(y="Matching Statistic Length", x="Contigs") 
  return(plot)
}

########################################################################
# Main method for the code
########################################################################

# Load in sub-sampled data for regular contigs
regular_path <- paste(data_working_dir, "regular_contigs.csv", sep="")
reg_df <- read.csv(regular_path, header=FALSE)

colnames(reg_df) <- c("name", "length")
reg_df["pos"] <- 0
reg_df["status"] <- "normal"
reg_df["name"] <- "All other contigs"

# Load in all values for contigs with most suspicious lengths
top_contig_files <- list.files(path=data_working_dir, pattern="top_contig_[[:digit:]]+_lengths\\.csv$")
df_list <- list()

for (path in top_contig_files) {
  curr_path <- paste(data_working_dir, path, sep="")
  curr_df <- read.csv(curr_path, header=FALSE)
  
  colnames(curr_df) <-  c("name", "pos", "length")
  df_list[[path]] <- curr_df
}
all_top_contigs_df <- do.call("rbind", df_list)
all_top_contigs_df["status"] <- "suspicious"

# Combine reguar and top contigs into one dataframe
total_df <- do.call("rbind", list(reg_df, all_top_contigs_df))

# Generate box plots...
plot_1 <- make_comparison_plot(total_df)
plot_1

plot_2 <- make_single_boxplot(total_df)
plot_2

# Combine the two plots into one ...
combined_plot <- ggdraw() +
                 draw_plot(plot_1) +
                 #draw_plot(plot_2, x = 0.7, y = .65, width = .15, height = .3)
                 draw_plot(plot_2, x=0.12, y=0.5, width=0.2, height=0.4)
combined_plot


# Create a plot showing MS across a suspicious contig ...
curr_path <- paste(data_working_dir, "top_contig_0_lengths.csv", sep="")
suspicious_df <- read.csv(curr_path, header=FALSE)

colnames(suspicious_df) <- c("name", "pos", "length")
contig_name <- suspicious_df[,1][1]

plot_3 <- make_length_plot(contig_name, suspicious_df, "#00BFC4")
plot_3

# Create a plot showing MS across a normal contig ...
curr_path <- paste(data_working_dir, "regular_contig_0_lengths.csv", sep="")
normal_df <- read.csv(curr_path, header=FALSE)

colnames(normal_df) <- c("name", "pos", "length")
contig_name <- normal_df[,1][1]

plot_4 <- make_length_plot(contig_name, normal_df, "#F8766D")
plot_4

# Combine the plots
total_plot <- ggarrange(combined_plot, 
                        ggarrange(plot_3, plot_4, nrow=2, labels=c("b", "c")), 
                        ncol=2, labels=c("a"))
total_plot

# Saving the plot: a vector and non-vector graphic
output_name <- paste(output_dir, "exp8_plot_assembly_1.jpeg", sep="")
ggsave(output_name, plot=total_plot, dpi=800, device="jpeg", width=12, height=6)

output_name <- paste(output_dir, "exp8_plot_assembly_1.pdf", sep="")
ggsave(output_name, plot=total_plot, dpi=800, device="pdf", width=12, height=6)


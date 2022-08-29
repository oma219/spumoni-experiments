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
library(pafr)
library(stringr)

########################################################################
# IMPORTANT: Data paths
########################################################################
data_working_dir <- "/Users/omarahmed/downloads/current_research/spumoni_exps/exp_8/trial_2/data/assembly_1/"
output_dir <- "/Users/omarahmed/downloads/current_research/spumoni_exps/exp_8/trial_2/plots/assembly_1/"

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
          labs(y="Pseudomatching Length", x=x_label)
  return(plot)         
}

make_single_boxplot <- function(input_df) {
  box_plot_stats <- boxplot(subset(input_df, name=="All other contigs")$length, plot=FALSE)$stats
  lower <- box_plot_stats[1] 
  upper <- box_plot_stats[5]
    
  plot <- ggplot(subset(input_df, name=="All other contigs"), aes(x=name, y=length, fill=status)) +
          #geom_boxplot(outlier.shape=NA) +
          geom_boxplot() +
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
        scale_y_continuous(limits=c(lower, upper)) +
        labs(y="", x="")
  return (plot)
}

make_comparison_plot <- function(input_df, lower, upper) {
  plot <- ggplot(input_df, aes(x=name, y=length, fill=status)) + 
          geom_boxplot(outlier.shape=NA) +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=10, face="bold"),
                axis.title.x=element_text(size =12),
                axis.title.y=element_text(size=12),
                legend.text=element_text(size=10),
                axis.text.x = element_text(angle = 0, vjust=0.5, size=6),
                axis.text.y =element_text(size=8, color="black"),
                legend.box="horizontal",
                legend.title=element_text(size=10),
                legend.position=c(0.8, 0.9)) +
        scale_y_continuous(limits=c(lower, upper * 1.1)) +
        scale_fill_discrete(name="", labels=c("Normal", "Suspicious")) +
        labs(y="Pseudomatching Length", x="Contigs") 
  return(plot)
}


process_paf_df <- function(input_df) {

  # create a new dataframe for dot plot
  columns <- c("qname", "x", "y") 
  out_df <- data.frame(qname=character(sum(input_df[,"alen"])),
                       x=numeric(sum(input_df[,"alen"])),
                       y=numeric(sum(input_df[,"alen"])))
  colnames(out_df) = columns

  # Go through each alignment and generate dots for each base
  curr_row <- 1
  for (i in 1:nrow(input_df)) {
    curr_query_seq <- input_df[i,"qname"]
    qstart <- input_df[i,"qstart"]; qend <- input_df[i,"qend"]
    tstart <- input_df[i,"tstart"]; tend <- input_df[i,"tend"]
    cigar <- input_df[i,"cg"]
    curr_x <- qstart; curr_y <- tstart
    
    # Go through each element in CIGAR
    for (x in str_extract_all(c(cigar), "[0-9]+[A-Z]")[[1]]) {
        length <- as.integer(substr(x, 1, nchar(x)-1))
        op <- substr(x, nchar(x), nchar(x))
        
        # Go through each base and add to the overall dataframe
        for (j in 1:length) {
            if (op == "M") {
                out_df[curr_row,] <- c(curr_query_seq, curr_x, curr_y)
                curr_row <- curr_row + 1
                curr_x <- curr_x + 1; curr_y <- curr_y + 1;
            } else if (op == "I") {
                curr_y <- curr_y + length; break;
            } else if (op == "D") {
                curr_x <- curr_x + length; break;
            } else {
              warning("Unrecognized character in the CIGAR string parsing")
            }
        }
    }
    print(curr_query_seq)
    print(i)
  }
  return(out_df)
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

# Find the upper/lower point for
min_value <- .Machine$integer.max
max_value <- 0
for (contig_name in unique(total_df$name)) {
  temp_df <- subset(total_df, name == contig_name)
  box_plot_stats <- boxplot(temp_df$length, plot=FALSE)$stats
  
  max_value <- max(box_plot_stats[5], max_value)
  min_value <- min(box_plot_stats[1], min_value)
}

# Generate box plots...
plot_1 <- make_comparison_plot(total_df, min_value, max_value)
plot_1

plot_2 <- make_single_boxplot(total_df)
plot_2

# Combine the two plots into one ...
combined_plot <- ggdraw() +
                 draw_plot(plot_1) +
                 draw_plot(plot_2, x = 0.7, y = .65, width = .15, height = .3)
                 #draw_plot(plot_2, x=0.12, y=0.5, width=0.2, height=0.4)
combined_plot


####################################################
# EDIT: Commented out these sub-plots
####################################################

# # Create a plot showing MS across a suspicious contig ...
# curr_path <- paste(data_working_dir, "top_contig_0_lengths.csv", sep="")
# suspicious_df <- read.csv(curr_path, header=FALSE)
# 
# colnames(suspicious_df) <- c("name", "pos", "length")
# contig_name <- suspicious_df[,1][1]
# 
# plot_3 <- make_length_plot(contig_name, suspicious_df, "#00BFC4")
# plot_3
# 
# # Create a plot showing MS across a normal contig ...
# curr_path <- paste(data_working_dir, "regular_contig_0_lengths.csv", sep="")
# normal_df <- read.csv(curr_path, header=FALSE)
# 
# colnames(normal_df) <- c("name", "pos", "length")
# contig_name <- normal_df[,1][1]
# 
# plot_4 <- make_length_plot(contig_name, normal_df, "#F8766D")
# plot_4

####################################################
# EDIT: Start of new sub-plot ...
####################################################

# load in paf file
input_paf_file <- "/Users/omarahmed/Downloads/test.paf"
paf_df <- read_paf(input_paf_file)

# remove un-needed alignments
paf_df_after_removal <- subset(paf_df, qname!="ctg7180000000054" |  tname!="NZ_CP086334.1")
paf_df_after_removal <- subset(paf_df_after_removal, qname!="ctg7180000000587" |  tname!="NZ_CP067426.1")
paf_df_after_removal <- subset(paf_df_after_removal, qname!="ctg7180000000651" |  tname!="NZ_CP067426.1")
paf_df_after_removal <- subset(paf_df_after_removal, qname!="ctg7180000000530" |  tname!="NZ_CP067426.1")

# generate data-frame for dot-plot
out_df <- process_paf_df(paf_df_after_removal)

# remove any rows that are empty in dataframe
sus_contigs <- c("ctg7180000000054", "ctg7180000000587", "ctg7180000000651", "ctg7180000000530")
out_df_for_plot <- subset(out_df, qname %in% sus_contigs)

# subset the dataframe for quick plotting
out_df_subset <- out_df_for_plot[seq(1, nrow(out_df_for_plot), 50),]

# make the dot plot
dot_plot <- ggplot(data=out_df_subset, aes(x=as.integer(x), y=as.integer(y))) + 
            geom_point(size=0.1) + 
            labs(x="Contig coordinates", y="Database coordinates") +
            theme_bw() +          
            theme(plot.title=element_text(hjust = 0.5, size=10, face="bold"),
                  axis.title.x=element_text(size =12),
                  axis.title.y=element_text(size=12),
                  legend.text=element_text(size=10),
                  axis.text.x = element_text(angle=0, size=6),
                  axis.text.y =element_text(size=6, color="black"),
                  legend.box="horizontal",
                  legend.title=element_text(size=10)) +
            facet_wrap(vars(qname), scales="free") 
dot_plot


####################################################
# Combine all sub-plots and save them
####################################################

# EDIT: old total plot (save just in case)
# total_plot <- ggarrange(combined_plot, 
#                         ggarrange(plot_3, plot_4, ncol=2, labels=c("b", "c")), 
#                         nrow=2, labels=c("a"))

total_plot <- ggarrange(plot_1, dot_plot, nrow=1, ncol=2, labels=c("a", "b"))
total_plot

# Saving the plot: a vector and non-vector graphic
output_name <- paste(output_dir, "exp8_plot_assembly_1.jpeg", sep="")
ggsave(output_name, plot=total_plot, dpi=800, device="jpeg", width=10, height=5)

output_name <- paste(output_dir, "exp8_plot_assembly_1.pdf", sep="")
ggsave(output_name, plot=total_plot, dpi=800, device="pdf", width=10, height=5)

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
library(ggbeeswarm)

########################################################################
# IMPORTANT: Data paths
########################################################################
data_working_dir <- "/Users/omarahmed/downloads/current_research/spumoni_exps/exp_8/trial_4/data/"
output_dir <- "/Users/omarahmed/downloads/current_research/spumoni_exps/exp_8/trial_4/plots/"

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
        #scale_y_continuous(limits=c(lower, upper)) +
        labs(y="", x="")
  return (plot)
}

make_comparison_plot <- function(input_df, lower, upper) {
  plot <- ggplot(input_df, aes(x=name, y=length, fill=status)) + 
          geom_boxplot(outlier.shape=NA) +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=10, face="bold"),
                axis.title.x=element_text(size=12, vjust=-0.5),
                axis.title.y=element_text(size=12),
                legend.text=element_text(size=10),
                axis.text.x = element_text(angle = 0, vjust=0.5, size=9),
                axis.text.y =element_text(size=8, color="black"),
                legend.box="horizontal",
                legend.title=element_text(size=10),
                legend.position=c(0.8, 0.9)) +
        scale_y_continuous(limits=c(lower, upper * 1.1)) +
        scale_fill_discrete(name="", labels=c("Normal", "Suspicious")) +
        scale_x_discrete(labels=c("All other contigs"="All other contigs\n(n=62482)", 
                                  "ctg7180000000054"="ctg7180000000054\n(n=4658)", 
                                  "ctg7180000000530"="ctg7180000000530\n(n=4671)", 
                                  "ctg7180000000587"="ctg7180000000587\n(n=4723)", 
                                  "ctg7180000000651"="ctg7180000000651\n(n=5476)")) +
        labs(y="Pseudomatching Length", x="Contigs") 
  return(plot)
}

make_comparison_plot_option2 <- function(input_df, lower, upper) {
  plot <- ggplot(input_df, aes(x=name, y=length, fill=status)) + 
          geom_boxplot(outlier.shape=NA) +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=10, face="bold"),
                axis.title.x=element_text(size=12, vjust=-0.5),
                axis.title.y=element_text(size=12),
                legend.text=element_text(size=10),
                axis.text.x = element_text(angle = 0, vjust=0.5, size=9),
                axis.text.y =element_text(size=8, color="black"),
                legend.box="horizontal",
                legend.title=element_text(size=10),
                legend.position=c(0.8, 0.9)) +
          scale_y_continuous(limits=c(lower, upper * 1.1)) +
          scale_fill_discrete(name="", labels=c("Normal", "Suspicious")) +
          scale_x_discrete(labels=c("ctg7180000000220"="ctg7180000000220\n(n=5673)",
                                    "ctg7180000000551"="ctg7180000000551\n(n=31165)", 
                                    "ctg7180000000580"="ctg7180000000580\n(n=189456", 
                                    "ctg7180000000708"="ctg7180000000708\n(n=11120)", 
                                     "ctg7180000000054"="ctg7180000000054\n(n=4658)", 
                                     "ctg7180000000530"="ctg7180000000530\n(n=4671)", 
                                     "ctg7180000000587"="ctg7180000000587\n(n=4723)", 
                                     "ctg7180000000651"="ctg7180000000651\n(n=5476)")) +
          labs(y="Pseudomatching Length", x="Contigs") 
  return(plot)
}



make_comparison_plot_option4 <- function(input_df, lower, upper) {
  plot <- ggplot(input_df, aes(x=name, y=length, fill=status)) + 
          #geom_boxplot(outlier.shape=NA) +
          geom_violin() +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=10, face="bold"),
                axis.title.x=element_text(size=12, vjust=-0.5),
                axis.title.y=element_text(size=12),
                legend.text=element_text(size=10),
                axis.text.x = element_text(angle = 0, vjust=0.5, size=9),
                axis.text.y =element_text(size=8, color="black"),
                legend.box="horizontal",
                legend.title=element_text(size=10),
                legend.position=c(0.8, 0.9)) +
          scale_y_continuous(limits=c(lower, upper * 1.1)) +
          scale_fill_discrete(name="", labels=c("Normal", "Suspicious")) +
          scale_x_discrete(labels=c("All other contigs"="All other contigs\n(n=62482)", 
                                    "ctg7180000000054"="ctg7180000000054\n(n=4658)", 
                                    "ctg7180000000530"="ctg7180000000530\n(n=4671)", 
                                    "ctg7180000000587"="ctg7180000000587\n(n=4723)", 
                                    "ctg7180000000651"="ctg7180000000651\n(n=5476)")) +
          labs(y="Pseudomatching Length", x="Contigs") 
  return(plot)
}



process_paf_df <- function(input_df) {

  # create a new dataframe for dot plot
  columns <- c("qname", "tname","x", "y") 
  out_df <- data.frame(qname=character(sum(input_df[,"alen"])),
                       tname=character(sum(input_df[,"alen"])),
                       x=numeric(sum(input_df[,"alen"])),
                       y=numeric(sum(input_df[,"alen"])))
  colnames(out_df) = columns

  # Go through each alignment and generate dots for each base
  curr_row <- 1
  for (i in 1:nrow(input_df)) {
    curr_query_seq <- input_df[i,"qname"]
    curr_target_seq <- input_df[i,"tname"]
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
                out_df[curr_row,] <- c(curr_query_seq, curr_target_seq, curr_x, curr_y)
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
# Main method for the code - Option #1
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

plot_1_alt <- make_comparison_plot_option4(total_df, min_value, max_value)
plot_1_alt

plot_2 <- make_single_boxplot(total_df)
plot_2

# Combine the two plots into one ...
combined_plot <- ggdraw() +
                 draw_plot(plot_1) +
                 draw_plot(plot_2, x = 0.07, y = .65, width = .15, height = .3)
                 #draw_plot(plot_2, x=0.12, y=0.5, width=0.2, height=0.4)
combined_plot


########################################################################
# Main method for the code - Option #2
########################################################################

# Load in all values for contigs with most suspicious lengths
top_contig_files <- list.files(path=data_working_dir, pattern="top_contig_[[:digit:]]+_lengths\\.csv$")
regular_contig_files <- list.files(path=data_working_dir, pattern="regular_contig_[[:digit:]]+_lengths\\.csv$")

# Load all the 'suspicious' contig files
df_list <- list()
for (path in top_contig_files) {
  curr_path <- paste(data_working_dir, path, sep="")
  curr_df <- read.csv(curr_path, header=FALSE)
  
  colnames(curr_df) <-  c("name", "pos", "length")
  df_list[[path]] <- curr_df
}
all_top_contigs_df <- do.call("rbind", df_list)
all_top_contigs_df["status"] <- "suspicious"

# Load all the 'normal' contig files
df_list <- list()
for (path in regular_contig_files) {
  curr_path <- paste(data_working_dir, path, sep="")
  curr_df <- read.csv(curr_path, header=FALSE)
  
  colnames(curr_df) <-  c("name", "pos", "length")
  df_list[[path]] <- curr_df
}
all_normal_contigs_df <- do.call("rbind", df_list)
all_normal_contigs_df["status"] <- "normal"

# Combine the data-frames
total_df <- do.call("rbind", list(all_normal_contigs_df, all_top_contigs_df))

# Find the upper/lower point for
min_value <- .Machine$integer.max
max_value <- 0
for (contig_name in unique(total_df$name)) {
  temp_df <- subset(total_df, name == contig_name)
  box_plot_stats <- boxplot(temp_df$length, plot=FALSE)$stats
  
  max_value <- max(box_plot_stats[5], max_value)
  min_value <- min(box_plot_stats[1], min_value)
}

# Make boxplot with the 4 suspicous contigs and 4 normal contigs
total_df$name <- factor(total_df$name, levels=c("ctg7180000000220",
                                                "ctg7180000000551", 
                                                "ctg7180000000580", 
                                                "ctg7180000000708", 
                                                "ctg7180000000054", 
                                                "ctg7180000000530", 
                                                "ctg7180000000587", 
                                                "ctg7180000000651"))
plot1 <- make_comparison_plot_option2(total_df, min_value, max_value)
plot1



########################################################################
# Main method for the code - Option #3
########################################################################

# Load the csv file with percentile numbers
percent_df <- read.csv(paste(data_working_dir, "assembly_percentile_report.csv", sep=""), header=TRUE)
colnames(percent_df) <- c("name", "contiglength", "quartileone", "median", "quartilethree", "ninetypercent", "mean")
percent_df["contignum"] <- seq.int(nrow(percent_df))

percent_df["status"] <- "normal"
percent_df[331,]["status"] <- "suspicious"
percent_df[332,]["status"] <- "suspicious"
percent_df[333,]["status"] <- "suspicious"
percent_df[334,]["status"] <- "suspicious"


cbbPalette <- c( "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73")

# Make the plot
plot1_option3 <-  ggplot(percent_df, aes(contiglength, quartileone)) +
                  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=2, fill = "green", alpha = 0.15) +
                  annotate("rect", xmin=-Inf, xmax=Inf, ymin=2, ymax=Inf, fill = "#56B4E9", alpha = 0.15) +
                  annotate("text", x=2.5e7, y=1, label="bold('Normal Contigs')", color="darkgreen", parse = TRUE) +
                  annotate("text", x=2.5e7, y=7, label="bold('Suspicious Contigs')", color="#56B4E9", parse = TRUE) +
                  annotate("text", x=4.15e7, y=2.90, label="PML threshold for contamination", color="black") +
                  annotate("text", x=4.8e6, y=5, label="ctg7180000000530", color="black", size=3.0) +
                  annotate("text", x=4.8e6, y=10, label="ctg7180000000587", color="black", size=3.0) +
                  annotate("text", x=4.8e6, y=11, label="ctg7180000000651", color="black", size=3.0) +
                  annotate("text", x=4.5e6, y=12, label="ctg7180000000054", color="black", size=3.0) +
                  geom_point(aes(fill=status), shape = 21, colour = "black", size = 1.5, stroke = 0.5, alpha=0.5) +
                  theme_classic() +
                  theme(legend.position="none",
                        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
                        aspect.ratio = 1/4,
                        axis.title.x=element_text(size=10),
                        axis.title.y=element_text(size=10)) +
                  geom_hline(yintercept=2, linetype="dashed", color = "black", size=1) +
                  scale_y_continuous(limits=c(0, 12), breaks=seq(0,12,2)) +
                  labs(x="Contig Length (After Digestion)",
                       y="25th percentile of\n PML distribution") +
                  scale_color_discrete(name="Status", labels=c("normal"="Normal", "suspicious"="Suspicious")) +
                  scale_fill_manual(values=c("darkgreen", "#56B4E9"))
                  #scale_fill_manual(values=cbbPalette)
                  #scale_fill_manual(values=("darkgreen", #56B4E9"))
plot1_option3




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
# Start of new sub-plot with dot-plots ...
####################################################

# load in paf file
input_paf_file <- "/Users/omarahmed/Downloads/current_research/spumoni_exps/exp_8/trial_3/test.paf"
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

# change facet titles to include length
contig_names <- list('ctg7180000000054'="ctg7180000000054 (27.1 kb)",
                     'ctg7180000000587'="ctg7180000000587 (22.9 kb)",
                     'ctg7180000000651'="ctg7180000000651 (23.3 kb)",
                     'ctg7180000000530'="ctg7180000000530 (22.8 kb)")

facet_labeller <- function(variable, value){
  return(contig_names[value])
}
contig_names['ctg7180000000054']

# make the dot plot
dot_plot <- ggplot(data=out_df_subset, aes(x=as.integer(x), y=as.integer(y), color=tname)) + 
            geom_point(size=0.1) + 
            labs(x="Contig coordinates", y="Database coordinates") +
            theme_bw() +          
            theme(plot.title=element_text(hjust = 0.5, size=10, face="bold"),
                  axis.title.x=element_text(size =10),
                  axis.title.y=element_text(size=10),
                  axis.text.x = element_text(angle=0, size=8),
                  axis.text.y =element_text(size=8, color="black"),
                  strip.text.x=element_text(size=10),
                  plot.margin = unit(c(0.25, 0.5, 0, 0.5), "cm"),
                  #aspect.ratio = 0.50,
                  legend.box="horizontal",
                  legend.position="bottom",
                  legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
                  legend.text=element_text(size=10),
                  legend.title=element_text(size=10)) +
            guides(colour = guide_legend(override.aes = list(size=2))) +
            facet_wrap(vars(qname), scales="free", labeller=facet_labeller) +
            #scale_color_discrete(name="Accessions", labels=c("NZ_CP067426.1", "NZ_CP086334.1", "NZ_CP091038.1"), values=cbbPalette)
            scale_color_manual(name="Accessions", labels=c("NZ_CP067426.1", "NZ_CP086334.1", "NZ_CP091038.1"), values=cbbPalette)
dot_plot

####################################################
# Combine all sub-plots and save them
####################################################

# EDIT: old total plot (save just in case)
# total_plot <- ggarrange(combined_plot, 
#                         ggarrange(plot_3, plot_4, ncol=2, labels=c("b", "c")), 
#                         nrow=2, labels=c("a"))

total_plot <- ggarrange(plot1_option3, dot_plot, nrow=2, ncol=1, labels=c("a", "b"), heights=c(1, 2))
total_plot

# Saving the plot: a vector and non-vector graphic
output_name <- paste(output_dir, "exp8_plot_assembly_1.jpeg", sep="")
ggsave(output_name, plot=total_plot, dpi=800, device="jpeg", width=7.19, height=6.71)

output_name <- paste(output_dir, "exp8_plot_assembly_1.pdf", sep="")
ggsave(output_name, plot=total_plot, dpi=800, device="pdf", width=7.19, height=6.71)

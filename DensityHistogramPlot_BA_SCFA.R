library(dplyr)

RA = read.csv("~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/Nov2022/Data_tables/MP3_species_ra.csv", header = T)

col_names <- RA[, 1]
t_RA = as.data.frame(t(RA[,-1]))
colnames(t_RA) <- col_names


t_RA <- t_RA %>%
  tibble::rownames_to_column(var = "Study.ID")
t_RA$Study.ID <- gsub("\\.", "-", t_RA$Study.ID)
t_RA$Study.ID <- gsub("IBS-040_02262018", "IBS-040", t_RA$Study.ID)
t_RA$Study.ID <- gsub("IBS-038", "IBS-0038", t_RA$Study.ID)
t_RA$Study.ID <- gsub("IBS-057", "IBS - 057", t_RA$Study.ID)

Demo = read.csv("~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/All/Clinical_Data_2.4.2022_baseline_demographic.csv", header = T)
BA_SCFA = read.csv("~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/All/Clinical_Data_2.4.2022_BA_SCFA.csv", header = T)
Demo_SCFA=merge.data.frame(BA_SCFA, Demo,by = "Study.ID", all.x = TRUE)
#Demo_SCFA_seq = merge.data.frame(t_RA, Demo_SCFA,by = "Study.ID", all.x = TRUE)[-(1:2),-(2:464)]
Demo_SCFA_seq = merge.data.frame(t_RA, Demo_SCFA,by = "Study.ID", all.x = TRUE)[-(1:2),]

# list of high abundance
high_abundance = read.csv("~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/Nov2022/Data_tables/MP3_species_ra_mean_sd0.1.csv")
pattern <- "^(.*?)_(mean|sd)$"

# Use sub to extract the desired part of the column names
high_abundance_bac_list <- unique(sub(pattern, "\\1", colnames(high_abundance)))
Abund=Demo_SCFA_seq[,high_abundance_bac_list[-1]]
Abund_Demo_SCFA = cbind(Demo_SCFA_seq[,1], Demo_SCFA_seq[,high_abundance_bac_list[-1]],Demo_SCFA_seq[,465:493])
colnames(Abund_Demo_SCFA)[1] = "Study.ID"

#Display the BA and SCFA density in three cohort groups. (summary ppt page 2)
# Create histograms with density curves using ggplot 
library(ggplot2)
library(gridExtra)


p_list <- lapply(109:113, function(i) {
  ggplot(Abund_Demo_SCFA, aes(x = Abund_Demo_SCFA[,i], fill = Type.of.IBS)) +
    geom_histogram(aes(y = ..density..), position = "dodge", binwidth = 1, alpha = 0.5) +
    geom_density(alpha = 0.5, aes(color = Type.of.IBS), size = 1) +
    scale_fill_manual(values = c("HV" = "green", "IBS - C" = "blue", "IBS - D" = "red")) +
    scale_color_manual(values = c("HV" = "green", "IBS - C" = "blue", "IBS - D" = "red")) +
    theme(
      legend.text = element_text(size = 8),  # Adjust legend text size
      plot.title = element_text(size = 14),   # Adjust title text size
      legend.title = element_text(size = 8), # Adjust legend title size
      axis.title.x = element_text(size = 12), # Adjust x-axis title size
      axis.title.y = element_text(size = 12), # Adjust y-axis title size
      axis.text.x = element_text(size = 10),  # Adjust x-axis text size
      axis.text.y = element_text(size = 10),   # Adjust y-axis text size
      legend.key.size = unit(0.5, "lines")    # Adjust legend key size
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    #labs(title = "Histograms with Density Curves by Cohort", x = colnames(Abund_Demo_SCFA)[i], y = "Density")
  labs(x = colnames(Abund_Demo_SCFA)[i], y = "Density")
  
  })

# Arrange the plots in a 3x2 layout
grid.arrange(grobs = p_list, ncol = 2)

# Save the plot as a PDF with increased width
pdf(file = "~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/April2024/Density_SCFA_BA_inIBSType_all_test2.pdf", width = 14, height = 8)
grid.arrange(grobs = p_list, ncol = 2)
dev.off()

#Bile Acid and Total Bile Acid (log sclae)
p1 = ggplot(Abund_Demo_SCFA, aes(x = Abund_Demo_SCFA[,107], fill = Type.of.IBS)) +
  geom_histogram(aes(y = ..density..), position = "dodge", binwidth = 1, alpha = 0.5) +
  geom_density(alpha = 0.5, aes(color = Type.of.IBS), size = 1) +
  scale_fill_manual(values = c("HV" = "green", "IBS - C" = "blue", "IBS - D" = "red")) +
  scale_color_manual(values = c("HV" = "green", "IBS - C" = "blue", "IBS - D" = "red")) +
  theme(
    legend.text = element_text(size = 8),  # Adjust legend text size
    plot.title = element_text(size = 14),   # Adjust title text size
    legend.title = element_text(size = 8), # Adjust legend title size
    axis.title.x = element_text(size = 12), # Adjust x-axis title size
    axis.title.y = element_text(size = 12), # Adjust y-axis title size
    axis.text.x = element_text(size = 10),  # Adjust x-axis text size
    axis.text.y = element_text(size = 10),   # Adjust y-axis text size
    legend.key.size = unit(0.5, "lines")    # Adjust legend key size
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Histograms with Density Curves by Cohort", x = colnames(Abund_Demo_SCFA)[107], y = "Density")

p2 = ggplot(Abund_Demo_SCFA, aes(x = log(Abund_Demo_SCFA[,108]), fill = Type.of.IBS)) +
  geom_histogram(aes(y = ..density..), position = "dodge", binwidth = 1, alpha = 0.5) +
  geom_density(alpha = 0.5, aes(color = Type.of.IBS), size = 1) +
  scale_fill_manual(values = c("HV" = "green", "IBS - C" = "blue", "IBS - D" = "red")) +
  scale_color_manual(values = c("HV" = "green", "IBS - C" = "blue", "IBS - D" = "red")) +
  theme(
    legend.text = element_text(size = 8),  # Adjust legend text size
    plot.title = element_text(size = 14),   # Adjust title text size
    legend.title = element_text(size = 8), # Adjust legend title size
    axis.title.x = element_text(size = 12), # Adjust x-axis title size
    axis.title.y = element_text(size = 12), # Adjust y-axis title size
    axis.text.x = element_text(size = 10),  # Adjust x-axis text size
    axis.text.y = element_text(size = 10),   # Adjust y-axis text size
    legend.key.size = unit(0.5, "lines")    # Adjust legend key size
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Histograms with Density Curves by Cohort", x = paste0(colnames(Abund_Demo_SCFA)[108],"   (log scale)"), y = "Density")

# Arrange the plots side by side
grid.arrange(p1, p2, ncol = 2)

# Save the arranged plots as a PDF with increased width
pdf(file = "~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/April2024/Density_BA_inIBSType2.pdf", width = 14, height = 4)
grid.arrange(p1, p2, ncol = 2)
dev.off()



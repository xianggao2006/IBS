
library(dplyr)
RA = read.csv("MP3_species_ra.csv", header = T)
col_names <- RA[, 1]
t_RA = as.data.frame(t(RA[,-1]))
colnames(t_RA) <- col_names


t_RA <- t_RA %>%
  tibble::rownames_to_column(var = "Study.ID")
t_RA$Study.ID <- gsub("\\.", "-", t_RA$Study.ID)
t_RA$Study.ID <- gsub("IBS-040_02262018", "IBS-040", t_RA$Study.ID)
t_RA$Study.ID <- gsub("IBS-038", "IBS-0038", t_RA$Study.ID)
t_RA$Study.ID <- gsub("IBS-057", "IBS - 057", t_RA$Study.ID)

Demo = read.csv("Clinical_Data_2.4.2022_baseline_demographic.csv", header = T)
BA_SCFA = read.csv("Clinical_Data_2.4.2022_BA_SCFA.csv", header = T)
Demo_SCFA=merge.data.frame(BA_SCFA, Demo,by = "Study.ID", all.x = TRUE)
#Demo_SCFA_seq = merge.data.frame(t_RA, Demo_SCFA,by = "Study.ID", all.x = TRUE)[-(1:2),-(2:464)]
Demo_SCFA_seq = merge.data.frame(t_RA, Demo_SCFA,by = "Study.ID", all.x = TRUE)[-(1:2),]

TransitTime = read.csv("Clinical_Data_2.4.2022_ColonicTransit.csv", header = T)
TransitTime_cleaned <- TransitTime[!is.na(TransitTime$Transit.Time.in.Days) & TransitTime$Transit.Time.in.Days != 0, ]
Demo_SCFA_seq_Transit = merge(Demo_SCFA_seq, TransitTime_cleaned, by = "Study.ID", all.x = FALSE, all.y = FALSE)

high_abundance = read.csv("MP3_species_ra_mean_sd0.1.csv")
pattern <- "^(.*?)_(mean|sd)$"

high_abundance_bac_list <- unique(sub(pattern, "\\1", colnames(high_abundance)))
Abund = Demo_SCFA_seq_Transit[,high_abundance_bac_list[-1]]
Abund_Demo_SCFA_Transit= cbind(Demo_SCFA_seq_Transit[,1], Demo_SCFA_seq_Transit[,high_abundance_bac_list[-1]],Demo_SCFA_seq_Transit[,465:501])
colnames(Abund_Demo_SCFA_Transit)[1] = "Study.ID"

###CCA (summary ppt page 3 and 4)
library(vegan)


######final results ########
# Running partial CCA
# BacteriaAbundanceMatrix ~ SCFA variables conditioned on "transit time"
result_pCCA <- cca(Abund ~ Acetic.Acid+Propionic.Acid+ Butyric.Acid +Total.SCFA + Acetic.to.Butyrate.ratio + Condition(Transit.Time.in.Days), data = Abund_Demo_SCFA_Transit)
permutest(result_pCCA, permutations = how(nperm = 999), model = "reduced")


cor_data <- Abund_Demo_SCFA_Transit[, c("Acetic.Acid", "Propionic.Acid", "Butyric.Acid", "Total.SCFA", "Acetic.to.Butyrate.ratio")]
cor(cor_data)

RsquareAdj(result_pCCA)

scfa_scores <- scores(result_pCCA, display = "bp")  # bp: biplot scores
species_scores <- scores(result_pCCA, display = "species")

sorted_bac_list <- list()

# Iterate over each row in the scfa_scores DataFrame
for (scfa_name in rownames(scfa_scores)) {
  env_arrow <- scfa_scores[scfa_name, ]
  unit_env_arrow <- env_arrow / sqrt(sum(env_arrow^2))
  species_scores_matrix <- as.matrix(species_scores)
  unit_env_vector <- as.vector(matrix(unit_env_arrow, ncol = 1))
  row_magnitudes <- sqrt(rowSums(species_scores_matrix^2))
  cosine_of_angle <- (species_scores_matrix %*% unit_env_vector) / row_magnitudes
  projections <- as.matrix(species_scores) %*% matrix(unit_env_arrow, ncol = 1)
  angles_in_degrees <- acos(cosine_of_angle) * (180 / pi)

  results <- as.data.frame(species_scores)
  results$Angle <- angles_in_degrees
  
  # Filter species based on angles
  species_filtered <- results[results$Angle < 30 | results$Angle > 330 | (results$Angle > 150 & results$Angle < 210), ]
  
  species_filtered$Projections <- projections[rownames(species_filtered), , drop = FALSE]
  ranked_species_filtered <- species_filtered[order(-abs(species_filtered$Projections)), ]
  #sorted_bac_with_scfa <- ranked_species_filtered[abs(ranked_species_filtered$Projections) >= 0.5, ]
  sorted_bac_with_scfa <- ranked_species_filtered
  sorted_bac_with_scfa$SCFA <- scfa_name
  sorted_bac_with_scfa$Species <- rownames(sorted_bac_with_scfa)
  sorted_bac_list[[scfa_name]] <- sorted_bac_with_scfa
}

sorted_bac_with_scfa <- do.call(rbind, sorted_bac_list)
sorted_bac_with_scfa_All <- sorted_bac_with_scfa[, c( "SCFA", "Species", "Angle", "Projections")]

write.csv(sorted_bac_with_scfa_All, "Bac_SCFA_BA_CCA_CosineProjection/Sorted_Bac_with_SCFA_pCCA_All.csv", row.names = FALSE)
unique_species <- unique(sorted_bac_with_scfa_All$Species)


print(unique_species)
write.csv(unique_species, "TopRankBac_SCFA_pCCA_All_nofilter.csv", row.names = FALSE)

plot(result_pCCA, scaling = 2)

species_scores <- scores(result_pCCA, display = "species")
label_indices <- rownames(species_scores) %in% unique_species

text(species_scores[label_indices, ], labels = rownames(species_scores)[label_indices], cex = 0.4, pos = 3)
dev.copy(png, file = "Bacteria_AllSCFA_ConditionalOnTT_partialCCA_topRankSpeceisNameOnly_Allcohorts_TriPlot.png")
dev.off() 



#### SCFA CCA to each cohorts (HV, IBS-C, IBS-D) (ppt page 6 -7)
HV = Abund_Demo_SCFA_Transit[Abund_Demo_SCFA_Transit$Type.of.IBS == "HV",]
result_pCCA <- cca(HV[,2:106] ~ Acetic.Acid+Propionic.Acid+ Butyric.Acid +Total.SCFA + Acetic.to.Butyrate.ratio + Condition(Transit.Time.in.Days), data = HV)
permutest(result_pCCA, permutations = how(nperm = 999), model = "reduced")

summary(result_pCCA)
RsquareAdj(result_pCCA)$r.squared
scores(result_pCCA)

scfa_scores <- scores(result_pCCA, display = "bp")  # bp: biplot scores
species_scores <- scores(result_pCCA, display = "species")

sorted_bac_list <- list()

for (scfa_name in rownames(scfa_scores)) {
  env_arrow <- scfa_scores[scfa_name, ]
  unit_env_arrow <- env_arrow / sqrt(sum(env_arrow^2))
  species_scores_matrix <- as.matrix(species_scores)
  unit_env_vector <- as.vector(matrix(unit_env_arrow, ncol = 1))
  row_magnitudes <- sqrt(rowSums(species_scores_matrix^2))
  cosine_of_angle <- (species_scores_matrix %*% unit_env_vector) / row_magnitudes
  projections <- as.matrix(species_scores) %*% matrix(unit_env_arrow, ncol = 1)
  angles_in_degrees <- acos(cosine_of_angle) * (180 / pi)

  results <- as.data.frame(species_scores)
  results$Angle <- angles_in_degrees
  species_filtered <- results[results$Angle < 30 | results$Angle > 330 | (results$Angle > 150 & results$Angle < 210), ]
  species_filtered$Projections <- projections[rownames(species_filtered), , drop = FALSE]
  
  ranked_species_filtered <- species_filtered[order(-abs(species_filtered$Projections)), ]
  #sorted_bac_with_scfa <- ranked_species_filtered[abs(ranked_species_filtered$Projections) >= 0.5, ]
  sorted_bac_with_scfa <- ranked_species_filtered

  sorted_bac_with_scfa$SCFA <- scfa_name
  sorted_bac_with_scfa$Species <- rownames(sorted_bac_with_scfa)
  sorted_bac_list[[scfa_name]] <- sorted_bac_with_scfa
}

sorted_bac_with_scfa <- do.call(rbind, sorted_bac_list)
sorted_bac_with_scfa_HV <- sorted_bac_with_scfa[, c("SCFA","Species", "Angle", "Projections")]


write.csv(sorted_bac_with_scfa_HV, "Bac_SCFA_BA_CCA_CosineProjection/Sorted_Bac_with_SCFA_pCCA_HV.csv", row.names = FALSE)
unique_species <- unique(sorted_bac_with_scfa_HV$Species)

print(unique_species)
write.csv(unique_species, "Bac_SCFA_BA_CCA_CosineProjection/TopRankBac_SCFA_pCCA_HV.csv", row.names = FALSE)


plot(result_pCCA, display = c("species", "bp"), scaling = 2)
species_scores <- scores(result_pCCA, display = "species")
label_indices <- rownames(species_scores) %in% unique_species

text(species_scores[label_indices, ], labels = rownames(species_scores)[label_indices], cex = 0.4, pos = 3)
dev.copy(png, file = "Bac_SCFA_BA_CCA_CosineProjection/Bacteria_AllSCFA_ConditionalOnTT_partialCCA_topRankSpeceisNameOnly_HV.png")
dev.off() 


######## IBS-C
IBSC = Abund_Demo_SCFA_Transit[Abund_Demo_SCFA_Transit$Type.of.IBS == "IBS - C",]

# Running partial CCA
result_pCCA <- cca(IBSC[,2:106] ~ Acetic.Acid+Propionic.Acid+ Butyric.Acid +Total.SCFA + Acetic.to.Butyrate.ratio + Condition(Transit.Time.in.Days), data = IBSC)
permutest(result_pCCA, permutations = how(nperm = 999), model = "reduced")

summary(result_pCCA)
RsquareAdj(result_pCCA)
scores(result_pCCA)

scfa_scores <- scores(result_pCCA, display = "bp")  # bp: biplot scores
species_scores <- scores(result_pCCA, display = "species")
sorted_bac_list <- list()

for (scfa_name in rownames(scfa_scores)) {
  env_arrow <- scfa_scores[scfa_name, ]
  unit_env_arrow <- env_arrow / sqrt(sum(env_arrow^2))
  
  species_scores_matrix <- as.matrix(species_scores)
  unit_env_vector <- as.vector(matrix(unit_env_arrow, ncol = 1))
  row_magnitudes <- sqrt(rowSums(species_scores_matrix^2))
  cosine_of_angle <- (species_scores_matrix %*% unit_env_vector) / row_magnitudes
  projections <- as.matrix(species_scores) %*% matrix(unit_env_arrow, ncol = 1)
  angles_in_degrees <- acos(cosine_of_angle) * (180 / pi)
  
  results <- as.data.frame(species_scores)
  results$Angle <- angles_in_degrees
  species_filtered <- results[results$Angle < 30 | results$Angle > 330 | (results$Angle > 150 & results$Angle < 210), ]
  species_filtered$Projections <- projections[rownames(species_filtered), , drop = FALSE]
  
  ranked_species_filtered <- species_filtered[order(-abs(species_filtered$Projections)), ]
  
  # Filter species with absolute projections greater than or equal to 0.5
  #sorted_bac_with_scfa <- ranked_species_filtered[abs(ranked_species_filtered$Projections) >= 0.5, ]
  sorted_bac_with_scfa <- ranked_species_filtered
  
  sorted_bac_with_scfa$SCFA <- scfa_name
  sorted_bac_with_scfa$Species <- rownames(sorted_bac_with_scfa)
  sorted_bac_list[[scfa_name]] <- sorted_bac_with_scfa
}


sorted_bac_with_scfa <- do.call(rbind, sorted_bac_list)
sorted_bac_with_scfa_IBSC <- sorted_bac_with_scfa[, c("SCFA", "Species", "Angle", "Projections")]


#write.csv(sorted_bac_with_scfa_IBSC, "Bac_SCFA_BA_CCA_CosineProjection/Sorted_Bac_with_SCFA_pCCA_IBSC.csv", row.names = FALSE)
write.csv(sorted_bac_with_scfa_IBSC, "Bac_SCFA_BA_CCA_CosineProjection/Sorted_Bac_with_SCFA_pCCA_IBSC_nofilter.csv", row.names = FALSE)

unique_species <- unique(sorted_bac_with_scfa_IBSC$Species)

print(unique_species)

#write.csv(unique_species, "Bac_SCFA_BA_CCA_CosineProjection/TopRankBac_SCFA_pCCA_IBSC.csv", row.names = FALSE)
write.csv(unique_species, "Bac_SCFA_BA_CCA_CosineProjection/TopRankBac_SCFA_pCCA_IBSC_nofilter.csv", row.names = FALSE)

plot(result_pCCA, display = c("species", "bp"), scaling = 2)

species_scores <- scores(result_pCCA, display = "species")
label_indices <- rownames(species_scores) %in% unique_species

text(species_scores[label_indices, ], labels = rownames(species_scores)[label_indices], cex = 0.4, pos = 3)
dev.copy(png, file = "~/Bac_SCFA_BA_CCA_CosineProjection/Bacteria_AllSCFA_ConditionalOnTT_partialCCA_topRankSpeceisNameOnly_IBSC.png")
dev.off() 



##### IBS_D
IBSD = Abund_Demo_SCFA_Transit[Abund_Demo_SCFA_Transit$Type.of.IBS == "IBS - D",]
result_pCCA <- cca(IBSD[,2:106] ~ Acetic.Acid+Propionic.Acid+ Butyric.Acid +Total.SCFA + Acetic.to.Butyrate.ratio + Condition(Transit.Time.in.Days), data = IBSD)
permutest(result_pCCA, permutations = how(nperm = 999), model = "reduced")


summary(result_pCCA)
RsquareAdj(result_pCCA)
scores(result_pCCA)

scfa_scores <- scores(result_pCCA, display = "bp")  # bp: biplot scores
species_scores <- scores(result_pCCA, display = "species")

sorted_bac_list <- list()

# Iterate over each row in the scfa_scores DataFrame
for (scfa_name in rownames(scfa_scores)) {
  env_arrow <- scfa_scores[scfa_name, ]
  unit_env_arrow <- env_arrow / sqrt(sum(env_arrow^2))
  species_scores_matrix <- as.matrix(species_scores)
  unit_env_vector <- as.vector(matrix(unit_env_arrow, ncol = 1))
  row_magnitudes <- sqrt(rowSums(species_scores_matrix^2))
  cosine_of_angle <- (species_scores_matrix %*% unit_env_vector) / row_magnitudes
  projections <- as.matrix(species_scores) %*% matrix(unit_env_arrow, ncol = 1)
  angles_in_degrees <- acos(cosine_of_angle) * (180 / pi)
 
  results <- as.data.frame(species_scores)
  results$Angle <- angles_in_degrees
  
 
  species_filtered <- results[results$Angle < 30 | results$Angle > 330 | (results$Angle > 150 & results$Angle < 210), ]
  species_filtered$Projections <- projections[rownames(species_filtered), , drop = FALSE]
  ranked_species_filtered <- species_filtered[order(-abs(species_filtered$Projections)), ]
  # Filter species with absolute projections greater than or equal to 0.5
  #sorted_bac_with_scfa <- ranked_species_filtered[abs(ranked_species_filtered$Projections) >= 0.5, ]
  sorted_bac_with_scfa <- ranked_species_filtered

  sorted_bac_with_scfa$SCFA <- scfa_name
  sorted_bac_with_scfa$Species <- rownames(sorted_bac_with_scfa)
  sorted_bac_list[[scfa_name]] <- sorted_bac_with_scfa
}

sorted_bac_with_scfa <- do.call(rbind, sorted_bac_list)
sorted_bac_with_scfa_IBSD <- sorted_bac_with_scfa[, c("SCFA", "Species", "Angle", "Projections")]
#write.csv(sorted_bac_with_scfa_IBSD, "Bac_SCFA_BA_CCA_CosineProjection/Sorted_Bac_with_SCFA_pCCA_IBSD.csv", row.names = FALSE)
write.csv(sorted_bac_with_scfa_IBSD, "Bac_SCFA_BA_CCA_CosineProjection/Sorted_Bac_with_SCFA_pCCA_IBSD.csv_nofilter", row.names = FALSE)

unique_species <- unique(sorted_bac_with_scfa_IBSD$Species)

# Print the unique species
print(unique_species)

#write.csv(unique_species, "~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/April2024/Bac_SCFA_BA_CCA_CosineProjection/TopRankBac_SCFA_pCCA_IBSD.csv", row.names = FALSE)
write.csv(unique_species, "~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/April2024/Bac_SCFA_BA_CCA_CosineProjection/TopRankBac_SCFA_pCCA_IBSD_nofilter.csv", row.names = FALSE)


plot(result_pCCA, display = c("species", "bp"), scaling = 2)
species_scores <- scores(result_pCCA, display = "species")
label_indices <- rownames(species_scores) %in% unique_species

text(species_scores[label_indices, ], labels = rownames(species_scores)[label_indices], cex = 0.4, pos = 3)
dev.copy(png, file = "Bac_SCFA_BA_CCA_CosineProjection/Bacteria_AllSCFA_ConditionalOnTT_partialCCA_topRankSpeceisNameOnly_IBSD.png")
dev.off() 



######## Generate file for each SCFA in order to do heatmap
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)

dataframes <- list(
 #All = sorted_bac_with_scfa_All,
  HV = sorted_bac_with_scfa_HV,
  IBSC = sorted_bac_with_scfa_IBSC,
  IBSD = sorted_bac_with_scfa_IBSD
)

scfas <- c("Acetic.Acid", "Propionic.Acid", "Butyric.Acid", "Acetic.to.Butyrate.ratio")

filter_and_combine_by_scfa <- function(scfa, dfs, cohort_names) {
  combined_list <- list()
  
  for (i in seq_along(dfs)) {
    df <- dfs[[i]]
    cohort <- cohort_names[i]
    
    # Filter rows starting with the specific SCFA
    filtered_df <- df[grepl(paste0("^", scfa), rownames(df)), ]
    
    # Add Cohort column
    filtered_df$Cohort <- cohort
    
    combined_list[[cohort]] <- filtered_df
  }
  
  combined_df <- do.call(rbind, combined_list)
  rownames(combined_df) <- NULL
  combined_df <- combined_df %>% select(-Angle)
  return(combined_df)
}


for (scfa in scfas) {
  combined_df <- filter_and_combine_by_scfa(scfa, dataframes, names(dataframes))
  #write.csv(combined_df, paste0("Bac_SCFA_BA_CCA_CosineProjection/",scfa, "_combined_3cohorts.csv"), row.names = TRUE)
  write.csv(combined_df, paste0("Bac_SCFA_BA_CCA_CosineProjection/",scfa, "_nofilter_combined_3cohorts.csv"), row.names = TRUE)
  
  print(paste("File for", scfa, "has been saved as", paste0(scfa, "_combined_3cohorts.csv")))
  
  data = combined_df 
  heatmap_data <- dcast(data, Species ~ Cohort, value.var = "Projections")
  heatmap_data[is.na(heatmap_data)] <- 0
  
  # Extract the CCA values and convert to a matrix for the heatmap
  cca_matrix <- as.matrix(heatmap_data[, -1])  # Exclude the taxa column for the matrix
  heatmap_data[-1] <- sapply(heatmap_data[-1], as.numeric)
  heatmap_data[is.na(heatmap_data)] <- 0 # Set row names of the matrix to the taxa
  
  long_heatmap_data <- melt(heatmap_data, id.vars = 'Species', variable.name = 'Cohort', value.name = 'Projections')
  p = ggplot(long_heatmap_data, aes(x = Cohort, y = Species, fill = Projections)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Adjust x-axis text size here
      axis.text.y = element_text(size = 7)  # Adjust y-axis text size here
    ) +
  labs(title = paste("Heatmap for", scfa))  # Add title here
  p
  ggsave(filename = paste0("~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/April2024/Bac_SCFA_BA_CCA_CosineProjection/", scfa, "_nofilter_3cohorts_heatmap.pdf"), plot = p, width = 10, height = 8)
  #ggsave(filename = paste0("~/OneDrive - Loyola University Chicago/Research/Shin_R03_06032021/April2024/Bac_SCFA_BA_CCA_CosineProjection/", scfa, "_heatmap_3cohorts.pdf"), plot = p, width = 10, height = 8)
  }



########### end of script ##########


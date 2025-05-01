library(dplyr)

# Use relative abundance to calculate the mean and sd of each cohort
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

Demo_seq = merge.data.frame(t_RA, Demo,by = "Study.ID", all = FALSE)
dim(Demo_seq)
###[1]  58 485
Demo_seq[Demo_seq$Control.patients.for.IBS.study == "Yes", ]$Study.ID  #17 control
Demo_seq[Demo_seq$Type.of.IBS == "IBS - C", ]$Study.ID  #15 IBS-C
Demo_seq[Demo_seq$Type.of.IBS == "IBS - D", ]$Study.ID #26 IBS-D

# Filter the dataframe to keep only rows with "Type of IBS" equal to "IBS - C"
Demo_seq_IBSC <- Demo_seq %>%
   filter(`Type.of.IBS` == "IBS - C")
Demo_seq_IBSD <- Demo_seq %>%
   filter(`Type.of.IBS` == "IBS - D")
Demo_seq_HV <- Demo_seq %>%
   filter(`Control.patients.for.IBS.study` == "Yes")

# Calculate the average and standard deviation for each bacterial column
HV_numeric_columns <- Demo_seq_HV[,1:464] %>%
   select(-c("Study.ID")) %>%
   lapply(as.numeric) %>%
   as.data.frame()

# Calculate the mean and standard deviation for each numeric column
HV_result <- HV_numeric_columns %>%
   summarise(across(everything(), list(mean = mean, sd = sd), na.rm = TRUE))


# Calculate the average and standard deviation for each bacterial column
IBSC_numeric_columns <- Demo_seq_IBSC[,1:464] %>%
   select(-c("Study.ID")) %>%
   lapply(as.numeric) %>%
   as.data.frame()

# Calculate the mean and standard deviation for each numeric column
IBSC_result <- IBSC_numeric_columns %>%
   summarise(across(everything(), list(mean = mean, sd = sd), na.rm = TRUE))

# Calculate the average and standard deviation for each bacterial column
IBSD_numeric_columns <- Demo_seq_IBSD[,1:464] %>%
   select(-c("Study.ID")) %>%
   lapply(as.numeric) %>%
   as.data.frame()

# Calculate the mean and standard deviation for each numeric column
IBSD_result <- IBSD_numeric_columns %>%
   summarise(across(everything(), list(mean = mean, sd = sd), na.rm = TRUE))



#### USE ALR transformed abundance for GLM analyis 

ALR = read.csv("MP3_species_count.ALR.csv", header = T)

col_names <- ALR[, 1]
t_ALR = as.data.frame(t(ALR[,-1]))
colnames(t_ALR) <- col_names

t_ALR <- t_ALR %>%
  tibble::rownames_to_column(var = "Study.ID")
t_ALR$Study.ID <- gsub("\\.", "-", t_ALR$Study.ID)
t_ALR$Study.ID <- gsub("IBS-040_02262018", "IBS-040", t_ALR$Study.ID)
t_ALR$Study.ID <- gsub("IBS-038", "IBS-0038", t_ALR$Study.ID)
t_ALR$Study.ID <- gsub("IBS-057", "IBS - 057", t_ALR$Study.ID)


Demo = read.csv("Clinical_Data_2.4.2022_baseline_demographic.csv", header = T)
Diet = read.csv("Clinical_Data_2.4.2022_baseline_diet.csv", header = T)
Data = merge.data.frame(Diet, Demo,by = "Study.ID", all.x = TRUE)

Demo_seq = merge.data.frame(t_ALR, Demo,by = "Study.ID", all = FALSE)

# list of high abundance
high_abundance = read.csv("MP3_species_ra_mean_sd0.1.csv")
pattern <- "^(.*?)_(mean|sd)$"

# Use sub to extract the desired part of the column names
high_abundance_bac_list <- sub(pattern, "\\1", colnames(high_abundance))

Demo_seq$Type.of.IBS = as.factor(Demo_seq$Type.of.IBS)
Demo_seq <- Demo_seq[Demo_seq$BMI..kg.m2. != "#DIV/0!", ]
Demo_seq$BMI..kg.m2. = as.numeric(Demo_seq$BMI..kg.m2.)


block = list()
P_value = list()

for (i in 2:464) {
  if( colnames(Demo_seq)[i] %in% high_abundance_bac_list) {
    model = glm(Demo_seq[,i] ~ Demo_seq$Type.of.IBS + Demo_seq$Age + Demo_seq$Gender+ Demo_seq$BMI..kg.m2., family = gaussian)
    newblock = cbind(coef(summary(model)),colnames(Demo_seq)[i])
    model2 = glm(Demo_seq[,i] ~ relevel(Demo_seq$Type.of.IBS, "IBS - C") + Demo_seq$Age + Demo_seq$Gender+ Demo_seq$BMI..kg.m2., family = gaussian)
    newblock2 = cbind(coef(summary(model2)),colnames(Demo_seq)[i])
    block<- rbind(newblock, newblock2)
    rows_to_keep <- c(2:3, 9)
    subset_block <- block[rows_to_keep, ]
    P_value = rbind(P_value, subset_block)
  }
}
Level = rownames(P_value)
P_value = cbind(Level, P_value)
P_value = as.data.frame(P_value)
colnames(P_value) = c("Level", "Estimate","Std", "t_value", "Pvalue", "Taxa")
rownames(P_value) <- NULL
P_value$Taxa = as.character(P_value$Taxa)
P_value$Level = as.character(P_value$Level)
P_value$Estimate = as.numeric(P_value$Estimate)
P_value$Std = as.numeric(P_value$Std)
P_value$t_value = as.numeric(P_value$t_value)
P_value$Pvalue = as.numeric(P_value$Pvalue)


# add the RA mean and sd 
t_high_abundance = as.data.frame(t(high_abundance[,-1]))
col_names = high_abundance[,1]
colnames(t_high_abundance) = col_names
odd_rows <- t_high_abundance[seq(1, nrow(t_high_abundance), by = 2), ]
even_rows <-t_high_abundance[seq(2, nrow(t_high_abundance), by = 2), ]

# Concatenate the even rows to the right of the odd rows
combined_df <- cbind(odd_rows, even_rows)
Taxa = sub(pattern, "\\1", rownames(combined_df))
combined_df2  = cbind(Taxa,combined_df)
final = merge.data.frame(P_value,combined_df2, by = "Taxa")
colnames(final) = c("Taxa", "Level", "Estimate", "Std", "t_value", "Pvalue", "HV_mean","IBS-C_mean", "IBS-D_mean", "HV_sd", "IBS-C_sd", "IBS-D_sd")
write.csv(final, "MetagenomicsALRTaxaGLM_IBS_age_gender_IBM_PValue_highAbundance0.1.csv", row.names = F)


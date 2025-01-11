#===============================================================================
# The SphingoHIIT study: linear mixed model
# Author: Nadia Weber, Seraina Fische, Justin Carrard, & Denis Infanger
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and import library
#-------------------------------------------------------------------------------

library(readxl)
library(car)
library(broom)
library(broom.mixed)
library(dplyr)
library(tidyverse)
library(ggdendro)
library(egg)
library(ggplot2)
library(reshape2)
library(lme4)
library(texreg)
library(nlme)
library(tidyr)
library(magrittr)
library(pander)
library(emmeans)
library(performance)
library(stringr)
library(writexl)
library(Matrix)
# Install the mice package if not already installed
if (!requireNamespace("mice", quietly = TRUE)) install.packages("mice")
library(mice)

#install.packages("multcomp")
library(multcomp)

# Install Cairo if not installed
if (!requireNamespace("Cairo", quietly = TRUE)) install.packages("Cairo")

rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

# Paths
setwd("") # Main path

data_path <- "./data" # Path for data
graphics_path <- "./output/graphics" # Path for graphics
text_path <- "./output/text" # Path for text-output

#-------------------------------------------------------------------------------
# Import data
#-------------------------------------------------------------------------------

# Import data
dat.all <- (read_excel(paste0(data_path, "/", "")))

# View(dat)
head(dat.all)
str(dat.all)

#-------------------------------------------------------------------------------
#Define time_points and lipid_values
#-------------------------------------------------------------------------------

# Defining dat_1
dat_1 <- dat.all %>% pivot_longer(cols=c(39 : 414), 
                                  names_to = 'time_point',
                                  values_to = 'lipid_values')

# Split the column "time_point" into one column "time point" and one column "lipid type"
dat_2 <- separate_wider_delim(cols = time_point, names = c("lipid_type", "time_point"), delim = "_", data = dat_1)
 
dat_3 <- spread(data = dat_2, key = "lipid_type", value = "lipid_values")

dat_3 <- dat_3 %>% 
  mutate(time_point = factor(
    time_point
    , levels = c("1", "2", "3", "4", "5", "6", "7", "8")
    , labels = c("1", "2", "3", "4", "5", "6", "7", "8")
  ))

#-------------------------------------------------------------------------------
# Create variables for interaction
#-------------------------------------------------------------------------------

dat_3 <- dat_3 %>%
  mutate(
    baseline = case_when(
      time_point %in% c("1", "2", "3")  ~ 0,
      TRUE ~ 1
    )
    , time_point_num = case_when(
      time_point %in% c("1", "2", "3")  ~ 0,
      time_point %in% "4" ~ 4,
      time_point %in% "5" ~ 5,
      time_point %in% "6" ~ 6,
      time_point %in% "7" ~ 7,
      time_point %in% "8" ~ 8
    )
    , group_num = case_when(
      Group %in% "HIIT" ~ 1,
      TRUE ~ 0
    )
    , int = factor(group_num*time_point_num*baseline) # 1 = HIIT, timepoints after baseline
  )

dat_3$Sex <- factor(dat_3$Sex)

# Create a new variable for intervention phase
dat_3 <- dat_3 %>%
  mutate(intervention_phase = ifelse(time_point %in% c("1", "2", "3"), "Pre-intervention", "Post-intervention"))

#-------------------------------------------------------------------------------
# Compute the sphingolipid ratios contained in the CERT1 score
#-------------------------------------------------------------------------------

# Add ratio columns
dat_3 <- dat_3 %>%
  mutate(
    Ratio_SL7_SL15 = SL7 / SL15,    # Ratio of SL7 to SL15
    Ratio_SL10_SL15 = SL10 / SL15,  # Ratio of SL10 to SL15
    Ratio_SL16_SL15 = SL16 / SL15   # Ratio of SL16 to SL15
  )

#-------------------------------------------------------------------------------
# List of column names that need to be converted to numeric
#-------------------------------------------------------------------------------

columns_to_convert <- c(
  "SL1", "SL2", "SL3", "SL4", "SL5", "SL6", "SL7", "SL8", "SL9", 
  "SL10", "SL11", "SL12", "SL13", "SL14", "SL15", "SL16", "SL17", 
  "SL18", "SL19", "SL20", "SL21", "SL22", "SL23", "SL24", "SL25", 
  "SL26", "SL27", "SL28", "SL29", "SL30", "SL31", "SL32", "SL33", 
  "SL34", "SL35", "SL36", "SL37", "SL38", "SL39", "SL40", "SL41", 
  "SL42", "SL43", "SL44", "SL45", "SL46", "SL47",
  "Ratio_SL7_SL15",
  "Ratio_SL10_SL15",
  "Ratio_SL16_SL15")

# Convert each specified column to numeric
dat_3[columns_to_convert] <- lapply(dat_3[columns_to_convert], function(x) as.numeric(as.character(x)))

#-------------------------------------------------------------------------------
# Impute missing Physical Activity data (should be done once only, subsequently
# go directly to row 291)
#-------------------------------------------------------------------------------

# Check for constant variables

constant_variables <- dat_3 %>% 
  group_by(ID) %>% 
  summarise(across(everything(), n_distinct)) %>% 
  select_if(function(.) all(. == 1)) %>% 
  ungroup()

changing_variables <- dat_3 %>% 
  group_by(ID) %>% 
  dplyr::summarise(across(everything(), n_distinct)) %>% 
  select_if(function(.) any(. > 1)) %>% 
  ungroup()

id_cols <- c("ID", names(constant_variables))

values_cols <- names(changing_variables)
values_cols <- values_cols[-which(values_cols %in% "ID")]

# Convert from long to wide format for imputation

dat_wide <- dat_3 %>%
  pivot_wider(
    id_cols = all_of(id_cols)
    , names_from = "time_point"
    , values_from = all_of(values_cols)
  )

pred <- quickpred(
  dat_wide
  , include = c(
    "TPA_min_day"
    , "VPA_min_day"
    , "VPA_min_day_one_min_bouts"
  )
  , exclude = c(
    "ID"
    , "Date_of_Birth"
    , "group_num"
    , paste0("baseline_", 1:8, collapse = ",")
    , paste0("time_point_num_", 1:8, collapse = ",")
    , paste0("time_point_", 1:8, collapse = ",")
  )
  , mincor = 0.6
  , method = "pearson"
)

table(rowSums(pred)) # Table of number of predictors
mean(rowSums(pred)) # Mean number of predictors

#-------------------------------------------------------------------------------
# Initialize imputation
#-------------------------------------------------------------------------------

imp_wide0 <- mice(dat_wide, maxit = 0, predictorMatrix = pred)

#-------------------------------------------------------------------------------
# See events/errors
#-------------------------------------------------------------------------------

imp_wide0$loggedEvents

#-------------------------------------------------------------------------------
# Inspect predictor matrix
#-------------------------------------------------------------------------------

# Only impute "TPA_min_day", "VPA_min_day", "VPA_min_day_one_min_bouts"

meth <- imp_wide0$method

meth[!names(meth) %in% c("TPA_min_day", "VPA_min_day", "VPA_min_day_one_min_bouts")] <- ""

meth[c("TPA_min_day", "VPA_min_day", "VPA_min_day_one_min_bouts")]

table(meth)

#-------------------------------------------------------------------------------
# Impute data
#-------------------------------------------------------------------------------

dat_imp_wide <- futuremice(dat_wide, m = 50, predictorMatrix = pred, method = meth, maxit = 10, parallelseed = 142857, n.core = parallelly::availableCores() - 1)

# Diagnostics

plot(dat_imp_wide, y = TPA_min_day + VPA_min_day + VPA_min_day_one_min_bouts ~ .it | .ms)

mice::stripplot(dat_imp_wide, TPA_min_day~.imp)
mice::stripplot(dat_imp_wide, VPA_min_day_one_min_bouts~.imp)
mice::stripplot(dat_imp_wide, VPA_min_day~.imp)

imp_comp <- mice::complete(dat_imp_wide, "long", include = TRUE)

#-------------------------------------------------------------------------------
# Convert imputed datasets from wide to long for modelling
#-------------------------------------------------------------------------------

var_names_long <- values_cols

varying_list <- vector("list", length = length(var_names_long))
names(varying_list) <- var_names_long

for (i in seq_along(var_names_long)) {
  varying_list[[i]] <- paste0(var_names_long[i], "_", 1:8)
}

imp_comp_long <- plyr::ddply(.data = imp_comp, .variables = ".imp", .fun = function(df) {
  
  # df <- subset(imp_comp, .imp == 1)
  
  df_long <- reshape(
    df
    , direction = "long"
    , idvar = c("ID")
    , varying = varying_list
    , v.names = var_names_long
  )
  
  # df_long$time_cat <- factor(df_long$time, levels = 1:4)
  
  df_long
  
})

imp_comp_long <- imp_comp_long %>% 
  arrange(.imp, .id, time_point)

names(imp_comp_long)

with(imp_comp_long, table(ID, .id))

# Convert to mids object

dat_imp_long <- as.mids(imp_comp_long, .imp = ".imp", .id = ".pid")

# Save imputed dataset for later retrieval in order to impute missing data only once

saveRDS(dat_imp_long, paste0(data_path, "/", "dat_imputed_long_2024-12-03.rds"))

# Load saved imputed dataset

dat_imp <- readRDS(paste0(data_path, "/", "dat_imputed_long_2024-12-03.rds"))

#-------------------------------------------------------------------------------
#Descriptive statistics
#-------------------------------------------------------------------------------

# Replace sphingolipid IDs with corresponding sphingolipid names
sphingo_names <- c(
  "SL1" = "Cer(d18:0/22:0) or Cer(m18:1/22:0)",
  "SL2" = "Cer(m18:1/24:1)",
  "SL3" = "1-deoxysphinganine(m18:0)",
  "SL4" = "1-deoxymethylsphinganine(m17:0)",
  "SL5" = "Cer(d18:1/12:0)",
  "SL6" = "Cer(d18:1/14:0)",
  "SL7" = "Cer(d18:1/16:0)",
  "SL8" = "Cer(d18:2/16:0)",
  "SL9" = "Cer(d18:1/17:0)",
  "SL10" = "Cer(d18:1/18:0)",
  "SL11" = "Cer(d18:1/18:1)",
  "SL12" = "Cer(d18:1/20:0)",
  "SL13" = "Cer(d18:1/22:0)",
  "SL14" = "Cer(d17:1/24:0)",
  "SL15" = "Cer(d18:1/24:0)",
  "SL16" = "Cer(d18:1/24:1)",
  "SL17" = "Cer1P(d18:1/16:0)",
  "SL18" = "Cer1P(d18:1/24:0)",
  "SL19" = "Cer(d18:0/16:0)",
  "SL20" = "Cer(d18:0/18:0)",
  "SL21" = "Cer(d18:0/20:0)",
  "SL22" = "Cer(d18:0/24:0)",
  "SL23" = "Cer(d18:0/24:1)",
  "SL24" = "HexCer(d18:1/16:0)",
  "SL25" = "HexCer(d18:1/18:0)",
  "SL26" = "HexCer(d18:1/18:1)",
  "SL27" = "HexCer(d18:1/24:1)",
  "SL28" = "HexSph(d18:1)",
  "SL29" = "Hex2Cer(d18:1/16:0)",
  "SL30" = "Hex2Cer(d18:1/24:0)",
  "SL31" = "Hex2Cer(d18:1/24:1)",
  "SL32" = "N-nervonoyl-1-deoxysphinganine",
  "SL33" = "N-nervonoyl-1-deoxysphingosine",
  "SL34" = "N-nervonoyl-1-desoxymethylsphinganine",
  "SL35" = "N-palmitoyl-1-deoxysphinganine",
  "SL36" = "SM(d18:1/12:0)",
  "SL37" = "SM(d18:1/16:0)",
  "SL38" = "SM(d18:1/17:0)",
  "SL39" = "SM(d18:1/18:0)",
  "SL40" = "SM(d18:1/18:1)",
  "SL41" = "SM(d18:1/22:0)",
  "SL42" = "SM(d18:1/24:0)",
  "SL43" = "SM(d18:1/24:1)",
  "SL44" = "Spa(d18:0)",
  "SL45" = "Sph(d18:1)",
  "SL46" = "Sph1P(d17:1)",
  "SL47" = "Sph1P(d18:1)",
  "Ratio_SL7_SL15" = "Cer(d18:1/16:0) to Cer(d18:1/24:0)",
  "Ratio_SL10_SL15" = "Cer(d18:1/18:0) to Cer(d18:1/24:0)",
  "Ratio_SL16_SL15" = "Cer(d18:1/24:1) to Cer(d18:1/24:0)"
)

# Define the new labels for the time points
time_point_labels <- c("B1", "B2", "B3", "2min", "15min", "30min", "60min", "24h")

# Define the new labels for the facet groups
facet_labels <- c("control" = "Control group", "HIIT" = "HIIT group")

variables_to_plot <- c(
  "SL1", "SL2", "SL3", "SL4", "SL5", "SL6", "SL7", "SL8", 
  "SL9", "SL10", "SL11", "SL12", "SL13", "SL14", "SL15",
  "SL16", "SL17", "SL18", "SL19", "SL20", "SL21", "SL22",
  "SL23", "SL24", "SL25", "SL26", "SL27", "SL28", "SL29",
  "SL30", "SL31", "SL32", "SL33", "SL34", "SL35", "SL36",
  "SL37", "SL38", "SL39", "SL40", "SL41", "SL42", "SL43",
  "SL44", "SL45", "SL46", "SL47",
  "Ratio_SL7_SL15",
  "Ratio_SL10_SL15",
  "Ratio_SL16_SL15")

# Loop over each variable and create a plot
for(i in seq_along(variables_to_plot)) {
    variable <- variables_to_plot[i]

    # Calculate max y-value for the annotations
    max_y_value <- max(dat_3[[variable]], na.rm = TRUE)
    
  plot <- ggplot(dat_3, aes(x = time_point, y = !!sym(variable), group = ID)) +
    geom_point(aes(colour = Group), position = position_dodge(width = 0.2)) +
    geom_line(aes(colour = Group), position = position_dodge(width = 0.2)) +
    geom_vline(xintercept = 3.5, linewidth = 0.5) +
    stat_summary(aes(group = NULL), fun = "mean", geom = "point", pch = 15, size = 4) +
    stat_summary(aes(group = Group), fun = "mean", geom = "line", linewidth = 1, linetype = 2) +
    scale_y_continuous(
      #trans = "log2"
      ) +
    scale_x_discrete(labels = time_point_labels) +
    #ylab(paste0("Molar fraction [unitless]")) +
    #ylab(paste0("Concentration [nM]")) +
    ylab(paste0("Ratio")) +
    xlab(paste0("Time points")) +
    ggtitle(paste0(sphingo_names[i])) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Group, labeller = as_labeller(facet_labels)) +
    theme(legend.position = "none",
          strip.background = element_rect(colour = "black", fill="white"),
          strip.text = element_text(face = "bold")) +
    annotate("text", x = 2, y = max_y_value * 1.2, label = "Pre-intervention", vjust = 2, size = 4, color = "black") +
    annotate("text", x = 6, y = max_y_value * 1.2, label = "Post-intervention", vjust = 2, size = 4, color = "black")

  # Print or save the plot
  print(plot)
  ggsave(paste("plot_", variable, ".png", sep=""), plot)
}

#-------------------------------------------------------------------------------
# Linear mixed models with a random effect for each subject
#-------------------------------------------------------------------------------

dat_3_data_frame <- data.frame(dat_3)

outcome_variables <- c(
  "SL1", "SL2", "SL3", "SL4", "SL5", "SL6", "SL7", "SL8", 
  "SL9", "SL10", "SL11", "SL12", "SL13", "SL14", "SL15",
  "SL16", "SL17", "SL18", "SL19", "SL20", "SL21", "SL22",
  "SL23", "SL24", "SL25", "SL26", "SL27", "SL28", "SL29",
  "SL30", "SL31", "SL32", "SL33", "SL34", "SL35", "SL36",
  "SL37", "SL38", "SL39", "SL40", "SL41", "SL42", "SL43",
  "SL44", "SL45", "SL46", "SL47",
  "Ratio_SL7_SL15",
  "Ratio_SL10_SL15",
  "Ratio_SL16_SL15")

# Create an empty data frame to store results
results_df <- data.frame()

for (i in outcome_variables) {
  
  cat("Calculating models for lipid:", i, "\n")
  
  form <- paste0(i, "~ time_point + int + Sex + Percent_Body_Fat + VO2peak_ml_min_kg + TPA_min_day + (1|ID)")
  
  model_fit <- with(dat_imp, lmer(format(form), REML = FALSE))
  
  model_fit.emm <- emmeans(model_fit, "int")
  
  contrasts_summary <- as.data.frame(summary(emmeans::contrast(model_fit.emm, "trt.vs.ctrl"), by = NULL, adjust = "mvt", infer = c(TRUE, TRUE)))
  
  #Extract relevant information from contrast_summary
  for (j in 1:nrow(contrasts_summary)) {
    temp_result <- data.frame(
      Variable = i,
      Contrast = contrasts_summary[j, "contrast"],  # Adjust column names as per your contrast_summary structure
      Estimate = contrasts_summary[j, "estimate"],
      StdError = contrasts_summary[j, "SE"],
      LowerCI = contrasts_summary[j, "lower.CL"],
      UpperCI = contrasts_summary[j, "upper.CL"],
      DF = contrasts_summary[j, "df"],
      TValue = contrasts_summary[j, "t.ratio"],
      PValue = contrasts_summary[j, "p.value"]
      )

    # Append this iteration's results to the results data frame
    results_df <- rbind(results_df, temp_result)
  }
}

#------------------------------------------------------------------------------------------
# Adjust p-values for multiple testing
#------------------------------------------------------------------------------------------

#Adjust p-values
hist(results_df$PValue, main = "Distribution of Raw P-Values", xlab = "P-Value", breaks = 20)
results_df <- results_df[order(results_df$PValue),]
results_df$BH <- p.adjust(results_df$PValue, method = "BH")

#Categorise BH p-values
results_df$BH_cat <- NA

results_df$BH_cat[results_df$BH > 0.05] <- "> 0.05"
results_df$BH_cat[results_df$BH <= 0.05 & results_df$BH > 0.01] <- "≤ 0.05"
results_df$BH_cat[results_df$BH <= 0.01 & results_df$BH > 0.001] <- "≤ 0.01"
results_df$BH_cat[results_df$BH <= 0.001 & results_df$BH > 0.0001] <- "≤ 0.001"
results_df$BH_cat[results_df$BH <= 0.0001] <- "≤ 0.0001"

results_df$BH_cat <- factor(results_df$BH_cat, levels = c("> 0.05", "≤ 0.05", "≤ 0.01", "≤ 0.001", "≤ 0.0001"))

#Order by decreasing estimates
results_df <- results_df[order(-results_df$Estimate),]


# Replace sphingolipid IDs with corresponding sphingolipid names
sphingolipid_names <- c(
  "SL1" = "Cer(d18:0/22:0) or Cer(m18:1/22:0)",
  "SL2" = "Cer(m18:1/24:1)",
  "SL3" = "1-deoxysphinganine(m18:0)",
  "SL4" = "1-deoxymethylsphinganine(m17:0)",
  "SL5" = "Cer(d18:1/12:0)",
  "SL6" = "Cer(d18:1/14:0)",
  "SL7" = "Cer(d18:1/16:0)",
  "SL8" = "Cer(d18:2/16:0)",
  "SL9" = "Cer(d18:1/17:0)",
  "SL10" = "Cer(d18:1/18:0)",
  "SL11" = "Cer(d18:1/18:1)",
  "SL12" = "Cer(d18:1/20:0)",
  "SL13" = "Cer(d18:1/22:0)",
  "SL14" = "Cer(d17:1/24:0)",
  "SL15" = "Cer(d18:1/24:0)",
  "SL16" = "Cer(d18:1/24:1)",
  "SL17" = "Cer1P(d18:1/16:0)",
  "SL18" = "Cer1P(d18:1/24:0)",
  "SL19" = "Cer(d18:0/16:0)",
  "SL20" = "Cer(d18:0/18:0)",
  "SL21" = "Cer(d18:0/20:0)",
  "SL22" = "Cer(d18:0/24:0)",
  "SL23" = "Cer(d18:0/24:1)",
  "SL24" = "HexCer(d18:1/16:0)",
  "SL25" = "HexCer(d18:1/18:0)",
  "SL26" = "HexCer(d18:1/18:1)",
  "SL27" = "HexCer(d18:1/24:1)",
  "SL28" = "HexSph(d18:1)",
  "SL29" = "Hex2Cer(d18:1/16:0)",
  "SL30" = "Hex2Cer(d18:1/24:0)",
  "SL31" = "Hex2Cer(d18:1/24:1)",
  "SL32" = "N-nervonoyl-1-deoxysphinganine",
  "SL33" = "N-nervonoyl-1-deoxysphingosine",
  "SL34" = "N-nervonoyl-1-desoxymethylsphinganine",
  "SL35" = "N-palmitoyl-1-deoxysphinganine",
  "SL36" = "SM(d18:1/12:0)",
  "SL37" = "SM(d18:1/16:0)",
  "SL38" = "SM(d18:1/17:0)",
  "SL39" = "SM(d18:1/18:0)",
  "SL40" = "SM(d18:1/18:1)",
  "SL41" = "SM(d18:1/22:0)",
  "SL42" = "SM(d18:1/24:0)",
  "SL43" = "SM(d18:1/24:1)",
  "SL44" = "Spa(d18:0)",
  "SL45" = "Sph(d18:1)",
  "SL46" = "Sph1P(d17:1)",
  "SL47" = "Sph1P(d18:1)",
  "Ratio_SL7_SL15" = "Cer(d18:1/16:0) to Cer(d18:1/24:0)",
  "Ratio_SL10_SL15" = "Cer(d18:1/18:0) to Cer(d18:1/24:0)",
  "Ratio_SL16_SL15" = "Cer(d18:1/24:1) to Cer(d18:1/24:0)"
)

# Replace the IDs in results_df$Variable with the sphingolipid names
results_df <- results_df %>%
  mutate(Variable = recode(Variable, !!!sphingolipid_names))

# Check the results to ensure replacements worked as expected
head(results_df)

# Rename the contrasts analysed with more descriptive labels
contrast_labels <- c(
  "int4 - int0" = "2min post-intervention vs. baseline",
  "int5 - int0" = "15min post-intervention vs. baseline",
  "int6 - int0" = "30min post-intervention vs. baseline",
  "int7 - int0" = "60min post-intervention vs. baseline",
  "int8 - int0" = "24h post-intervention vs. baseline"
)

# Replace the values in results_df$Contrast with the descriptive labels
results_df <- results_df %>%
  mutate(Contrast = recode(Contrast, !!!contrast_labels))

# Check the updated results_df
head(results_df)

#-------------------------------------------------------------------------------
# Creating a forest plot to visualise the results
#-------------------------------------------------------------------------------

# Arrange data for a cleaner display
results_df <- results_df %>%
  arrange(Variable, desc(Estimate))

# Plot
p <- 
  ggplot(results_df, aes(x = Estimate, y = reorder(paste(Variable, Contrast, sep = ": "), Estimate))) +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Effect of a single HIIT session (vs. physical rest) on sphingolipid species at different time points post-intervention",
    x = "Estimate and 95% confidence intervals",
    y = "Sphingolipids at post-intervention time points vs. baseline"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
p

# Save plot with CairoPDF to support special characters
Cairo::CairoPDF(file = paste(graphics_path, "forest_plot.pdf", sep = "/"), width = 20, height = 40, pointsize = 12)  # Adjust pointsize if needed
print(p)
dev.off()

#rename columns
names(results_df)[1] <- "Sphingolipid"
names(results_df)[2] <- "Contrast"
names(results_df)[3] <- "Estimate"
names(results_df)[4] <- "Standard error"
names(results_df)[5] <- "95% CI lower bound for β coefficient"
names(results_df)[6] <- "95% CI higher bound for β coefficient"
names(results_df)[7] <- "Degrees of Freedom"
names(results_df)[8] <- "t-values"
names(results_df)[9] <- "p-values"
names(results_df)[10] <- "BH p-values"
names(results_df)[11] <- "Categorical BH p-value"

#-------------------------------------------------------------------------------
# Export results as xlsx
#-------------------------------------------------------------------------------

write_xlsx(results_df,"")


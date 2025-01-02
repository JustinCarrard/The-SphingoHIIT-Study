#===============================================================================
# The SphingoHIIT study: linear mixed model
# Author: Nadia Weber, Seraina Fische, Justin Carrard, & Denis Infanger
#===============================================================================
#-------------------------------------------------------------------------------
# Change format of the data
#-------------------------------------------------------------------------------

# Load necessary libraries
library(readxl)
library(tidyr)
library(dplyr)
library(writexl)

# Paths
setwd("") # Main path

data_path <- "./data" # Path for data
graphics_path <- "./output/graphics" # Path for graphics
text_path <- "./output/text" # Path for text-output

# import data raw
data <- (read_excel(paste0(data_path, "/", "SphingoHIIT_raw_data.xlsx")))

# Calculate the sum of sphingolipids for each time point
data$Total_Sphingolipids <- rowSums(data[ , 3:ncol(data)])

# Normalize sphingolipid concentrations 
# (only if normalized concentrations should be analysed, otherwise skip this part)
normalized_data <- data
normalized_data[ , 3:(ncol(data)-1)] <- data[ , 3:(ncol(data)-1)] / data$Total_Sphingolipids

# Replace time points from 'a'-'h' to '1'-'8'
normalized_data <- normalized_data %>%
  mutate(Time_point = recode(Time_point, 
                             "a" = "1", "b" = "2", "c" = "3", "d" = "4", 
                             "e" = "5", "f" = "6", "g" = "7", "h" = "8"))

# Reshape normalized data to long format with one row per sphingolipid per time point
normalized_long_data <- normalized_data %>%
  pivot_longer(
    cols = starts_with("SL"),  # Select columns for sphingolipid concentrations
    names_to = "sphingolipid",  # Name for the sphingolipid column
    values_to = "concentration"  # Name for the concentration values
  )

# Remove the Total_Sphingolipids and Time_points columns
normalized_long_data$Total_Sphingolipids <- NULL

# Create unique column names for each sphingolipid-time combination
normalized_long_data <- normalized_long_data %>%
  unite("sphingolipid_time", sphingolipid, Time_point, sep = "_") 

# Pivot data back to wide format to create separate columns for each sphingolipid-time combination
normalized_wide_data <- normalized_long_data %>%
  pivot_wider(
    names_from = sphingolipid_time,
    values_from = concentration
  )

# Custom ordering function to sort in desired order
ordered_columns <- c("ID", 
                     paste0("SL", rep(1:47, each = 8), "_", rep(1:8, times = 47)))

# Reorder columns in wide_data based on ordered_columns
normalized_wide_data <- normalized_wide_data %>%
  select(all_of(ordered_columns))

# View the final data
print(normalized_wide_data)

# Export the final data as .xlsx
write_xlsx(wide_data,"")

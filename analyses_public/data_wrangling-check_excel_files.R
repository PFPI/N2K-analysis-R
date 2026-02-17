# Install the package (only need to do this once)
install.packages("readxl")

# Load the package for your current session
library(readxl)
library(tidyverse)

# Define the path to your Excel file
file_path <- "imports/FullEU_PointsMT1kSqM_N2K_Intersected_Export.csv"

# Read the data from the first sheet into a variable called 'df'
df <- read.csv(file_path)

# View the first few rows of your new data frame
head(df)
str(df)
intersected_df <- df

file_path_2 <- "imports/Vectorized_FullEU_LossYear2001_For_Cleanup_Export.csv"
vectorized_df <- read.csv(file_path_2)
head(vectorized_df)
str(vectorized_df)

## use intersected, it has site info. 
## get the sitecodes we want to pull from because those are correct
sitecodes_file <- "importsisual_verification_points_all_regions.csv" #nolint
sitecodes_df <- read.csv(sitecodes_file)
head(sitecodes_df)
sitecodes_to_filter_by <- unique(sitecodes_df$SITECODE)
str(sitecodes_to_filter_by)

filtered_df <- intersected_df %>% filter(SITECODE %in% sitecodes_to_filter_by)
str(filtered_df)

sliced_df <- filtered_df %>% group_by(SITECODE) %>% slice_sample(n = 20) %>% ungroup()
head(sliced_df)
write.csv(sliced_df, "csvs/new_sampled_sites.csv")
unique(filtered_df$MS)
unique(sitecodes_df$COUNTRY)

## missing Austria, Finland, Ireland, Lithuania, Poland
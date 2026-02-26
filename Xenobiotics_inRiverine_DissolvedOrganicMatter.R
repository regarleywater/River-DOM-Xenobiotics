### Code to assess the contribution of xenobiotic contaminants to detected dissolved organic matter by LC-MS/MS. 
### Note: Code replicates main Hitchhikers script for calling files, blank removal and data set selection for data analysis, which is linked below
### https://github.com/Functional-Metabolomics-Lab/FBMN-STATS/blob/main/R/.ipynb_checkpoints/Stats_Untargeted_Metabolomics-checkpoint.ipynb


# The main components of this code are: (1) load the libraries and (2) call the data files and process them for analysis. 
# Subsequently (3) the non-target data is analyzed to determine (3.1) the percentage of xenobiotic features 
# contribute to detected dissolved organic matter (DOM) peak area and (3.2) how xenobiotic features 
# correspond to population density and alter the compositional character of DOM.


##### 1.  LOAD LIBRARIES #####
# Global settings for plot size in the output cell:
options(repr.plot.width=10, repr.plot.height=10, res=600) # the parameters: width, height & resolution can be changed

library("tidyverse")
library("ggdist")
library("ggplot2")
library("cowplot")
library("dplyr")
library("igraph")
library("grid")
library("gtable")


##### 2.  DATA PROCESSING #####
##### 2.1 CALL FILES and SET UP FILES for ANALYSIS #####

# Set directory
Directory <- "C:/Users/Lana/Documents/Projects/R_projects/DataAnalysis/MassSpec_postFBMN/NeckarSpree2_BOTH/Xenobiotics_Riverine_DissolvedOrganicMatter"
setwd(Directory)

# Create folders for output if not already there
dir.create(file.path(Directory, "02_Output_DataProducts"))
dir.create(file.path(Directory, "03_Output_Figures"))

# Call file with supplemental data that can go in scatter plots 
fileName_xenobiotic <- "C:/Users/Lana/Documents/Projects/R_projects/DataAnalysis/MassSpec_postFBMN/NeckarSpree2_BOTH/Xenobiotics_Riverine_DissolvedOrganicMatter/01_Data/S7a_Features_Groups_Assignment.csv"
df_csv_xenobiotic <- read.csv(fileName_xenobiotic,  head=TRUE, sep=",")
colnames(df_csv_xenobiotic)

# Call file with MN results
fileName_edges <- "C:/Users/Lana/Documents/Projects/R_projects/DataAnalysis/MassSpec_postFBMN/NeckarSpree2_BOTH/Xenobiotics_Riverine_DissolvedOrganicMatter/01_Data/S5_Features_Edges_MolecularNetwork.csv"
df_csv_edges <- read.csv(fileName_edges,  head=TRUE, sep=",")
dim(df_csv_edges)

# Call file with MN results
fileName_SIRIUS_formulae <- "C:/Users/Lana/Documents/Projects/R_projects/DataAnalysis/MassSpec_postFBMN/NeckarSpree2_BOTH/Xenobiotics_Riverine_DissolvedOrganicMatter/01_Data/S6_Features_Annotations_SIRIUS.csv"
df_csv_SIRIUS_formulae <- read.csv(fileName_SIRIUS_formulae,  head=TRUE, sep=",")
dim(df_csv_SIRIUS_formulae)

# Identify names of feature quant and metadata tables in directory
input_str <- "S2_Features_PeakAreas_Raw.csv,S1_Samples_Metadata.tsv,S4_Features_Annotations_GNPS.tsv"

# Read the files indicated by the indecies 
input <- (strsplit(input_str, ",")[[1]])

# Read feature quant table
first_line <- readLines(file.path(Directory, "01_Data", input[1]), n = 1)
if (length(strsplit(first_line, ';')[[1]]) > 1) {
  ft <- read.csv(file.path(Directory, "01_Data", input[1]), header = T, check.names = F, sep = ';') # in case, ';' is the separator
} else {
  ft <- read.csv(file.path(Directory, "01_Data", input[1]), header = T, check.names = F)
}

# Read metadata
md <- read.csv(file.path(Directory, "01_Data", input[2]), header = T, check.names = F, sep = '\t') # mention seperator as "/t"(tab-separated) in case of txt or tsv fi

### Annotations ###
# Load GNPS annotations:
an_gnps <- read.csv(file.path(Directory, "01_Data", input[3]), header = T, check.names = F, sep = '\t') 

# Check table dimensions
dim(ft) 
dim(md)
dim(an_gnps) 

### Reformat files ###
# Function to summarize metadata
InsideLevels <- function(metatable){
  LEVELS <- c() #creating empty vector to store information
  typ <-c()
  COUNT <- c()
  for(i in 1:ncol(metatable)){ # for each metadata column
    temp <- as.data.frame(table(metatable[,i])) #table function gives the category in each column and the count of each category
    x <- temp$Var1 #getting the name of each category in every column
    if(is.double(metatable[,i])==T){ # for numeric columns in metadata table, round the category values
      x=round(as.double(x),2)
    } 
    LEVELS <- rbind(LEVELS,toString(x)) # adding all the category values in a row
    COUNT <- rbind(COUNT,toString(temp$Freq)) # getting the frequency of each level in every column
    typ <- rbind(typ,class(metatable[,i])) # getting the class of each column
  }
  out <- data.frame(INDEX=c(1:ncol(metatable)), #creating an output dataframe with 1st column as INDEX
                    ATTRIBUTES=colnames(metatable), #2nd column ATTRIBUTES will be the column name of metadata table
                    LEVELS, #3rd column LEVELS will give the different categories in each ATTRIBUTE
                    COUNT, #4th column COUNT will give the number of files present with each category
                    'ATTRIBUTE_CLASS'=typ, row.names=NULL) #Final column indicating the Class or datatype of each ATTTRIBUTE
  return(out)
}

# Check out metadata dimensions & content
ncol(md) #number of columnns of metadata
InsideLevels(md[, 2:ncol(md)]) #excluding 1st filename

# Check compatibility of annotations with feature table ids by checking their class
identical(class(ft$`row ID`),class(an_gnps$`#Scan#`))

# Arrange annotations by Scan number
an_final <- an_gnps %>% arrange(`#Scan#`) #arranging by scan ID

# Reformat GNPS annotations
# Function to compare between gnps_compound_name and its library match (analog)
combine_names <- function(compound_name) {
  return(paste(compound_name))
}

# Consolidate multiple annotations for a single '#Scan#' into one combined name
an_final_single <- an_final %>%
  group_by(`#Scan#`) %>%
  summarise(Combined_Name = combine_names(Compound_Name[1])) %>%
  ungroup() %>%
  as.data.frame()

# Merge annotations with feature table
ft_an <- merge(ft, an_final_single, by.x="row ID", by.y="#Scan#", all.x= TRUE) 
head(ft_an, 2)
dim(ft_an) 

# Arrange feature table and metadata in same order
# Create duplicate (working) files
new_ft <- ft #storing the files under different names to preserve the original files
new_md <- md

# Clean the new files
colnames(new_ft) <- gsub(' Peak area','',colnames(new_ft)) #removing Peak area extensions from the column names of ft
new_ft <- new_ft[order(new_ft$`row ID`),,drop=F] #arranging the rows of ft file by  by ascending order of row ID
new_ft <- new_ft[,colSums(is.na(new_ft))<nrow(new_ft)] #removing if any NA columns present in the ft file,
new_md <- new_md[,colSums(is.na(new_md))<nrow(new_md)] #removing if any NA columns present in the md file,
new_md <- new_md[apply(new_md != "", 1, any), ] # Removing rows that are completely filled with empty strings,
new_md <- new_md[, apply(new_md != "", 2, any)] # Removing columns that are completely filled with empty strings

# Remove the (front & tail) spaces, if any present, from the filenames of md,
new_md$filename <- trimws(new_md$filename, which = c("both"))
rownames(new_md) <- new_md$filename
new_md <- new_md[, -which(names(new_md) == "filename")]

# Update row names of feature table
if(exists("ft_an")){identical(ft_an$`row ID`,new_ft$`row ID`)} #should return TRUE if you have annotation file

# Changing the row names of the files into the combined name as "XID_mz_RT":
rownames(new_ft) <- paste(paste0("X",new_ft$`row ID`),
                          round(new_ft$`row m/z`,digits = 3),
                          round(new_ft$`row retention time`,digits = 3),
                          if(exists("ft_an")){ft_an$Combined_Name}, 
                          sep = '_') 

# Remove the trailing underscore at rownames
rownames(new_ft) <- sub("_$", "", rownames(new_ft))

# In the feature table, identify which columns correspond to samples
# Check if columns contain 'mzXML' or 'mzML' extensions
if (any(grepl('.mzML', colnames(new_ft)))) {
  # Picking only the files with column names containing 'mzXML' or 'mzML'
  new_ft <- new_ft[, grepl('.mzML', colnames(new_ft))]
  
  # Message if both .mzXML and .mzML files are present
  if (any(grepl('.mzXML$', colnames(new_ft))) && any(grepl('.mzML$', colnames(new_ft)))) {
    print("Both .mzXML and .mzML file types are present in the data")
  }
} else {
  # Ask the user for the extension if neither 'mzXML' nor 'mzML' is found
  
  your_extension <- readline('If your file extension is not .mzML or .mzXML, enter your extension (e.g., ".txt"): ')
  new_ft <- new_ft[, grepl(your_extension, colnames(new_ft))]
}

# Checking the files again to see if the above changes have been made
dim(new_ft)
dim(new_md)

# Check overlap between the feature table and metadata
new_ft <- new_ft[,order(colnames(new_ft)), drop=F] #ordering the ft by its column names
new_md <- new_md[order(rownames(new_md)),, drop=F] #ordering the md by the 1st column filename

# Check how many files in the metadata are also present in the feature table
table(rownames(new_md) %in% colnames(new_ft))

# Check if the sample names the same
identical(rownames(new_md), colnames(new_ft))

# Chech which file names in the metadata are not in the feature table
setdiff(rownames(new_md),colnames(new_ft))
# print(colnames(new_ft)) # uncomment to check the column names of new_ft

# Checking the dimensions of our new ft and md:
cat("The number of rows and columns in our original ft is:",dim(ft),"\n")
cat("The number of rows and columns in our new ft is:",dim(new_ft),"\n")
cat("The number of rows and columns in our new md is:",dim(new_md))

# Remove rows that sum to 0 (i.e., no remaining samples have peak areas > 0)
new_ft <- new_ft[rowSums(new_ft != 0) > 0,]
dim(new_ft)
cat("The number of rows and non-zero columns in our new ft is:",dim(new_ft))

# Create "Figures" folder if not already there
dir.create(file.path(Directory, "Figures"))

# Transpose and merge feature table and metadata
ft_t <- as.data.frame(t(new_ft)) #transposing the ft
ft_t <- ft_t %>% mutate_all(as.numeric)  #converting all values to numeric
identical(rownames(new_md),rownames(ft_t)) #should return TRUE now

# Check out the results
head(ft_t,3)

# Merging metadata (new_md) and transposed feature table based on the sample names
ft_merged_with_md <- merge(new_md, ft_t, by= 0, all.x=TRUE) #by =0 indicates the rownames of new_md and ft_t
head(ft_merged_with_md, 3)


##### 2.2 BLANK REMOVAL and SUBTRACTION #####

# Get the index levels in the data
InsideLevels(new_md)

# Enter indecies
# Sample_attribute <- as.numeric(readline('Enter the index number of the attribute containing sample and blanks information: '))
sample_attribute <- as.numeric("5")
unique_sampletypes <- unique(new_md[, sample_attribute])

# Display the unique sample types along with their index
cat("Available sample types:/n")
print(unique_sampletypes)

# Provide index number of blank label
blank_ID_str <- readline('Enter the index number(s) of the blanks: ')
#blank_ID <- as.numeric(strsplit(blank_ID_str, ",")[[1]])
blank_ID <- as.numeric("1")

# Provide index number of sample label
sample_ID_str <- readline('Enter the index number(s) of the samples: ')
# sample_ID <- as.numeric(strsplit(sample_ID_str, ",")[[1]])
sample_ID <- as.numeric("2")

# Filtering the rows from metadata with the condition = blank and sample
md_Blank <- new_md[new_md[, sample_attribute] %in% unique_sampletypes[blank_ID],]
md_Samples <- new_md[new_md[, sample_attribute] %in% unique_sampletypes[sample_ID],]

# Getting the corresponding rows from ft_t
Blank <- ft_t[which(rownames(ft_t) %in% (rownames(md_Blank))), , drop=F]
Samples <- ft_t[which(rownames(ft_t) %in% (rownames(md_Samples))), , drop=F]
InsideLevels(new_md)

# Check out results
dim(Blank) 
dim(Samples)

# Set Blank cutoff
# When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
# Cutoff <- as.numeric(readline('Enter Cutoff value between 0.1 & 1:')) # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3
Cutoff <- as.numeric(0.1)

# Getting mean for every feature in blank and Samples in a data frame named 'Avg_ft'
Avg_ft <- data.frame(Avg_blank=colMeans(Blank, na.rm= F)) # set na.rm = F to check if there are NA values. When set as T, NA values are changed to 0
Avg_ft$Avg_samples <- colMeans(Samples, na.rm= F) # adding another column 'Avg_samples' for feature means of samples

# Getting the ratio of blank vs Sample
Avg_ft$Ratio_blank_Sample <- (Avg_ft$Avg_blank+1)/(Avg_ft$Avg_samples+1)

# Creating a bin with 1s when the ratio>Cutoff, else put 0s
Avg_ft$Bg_bin <- ifelse(Avg_ft$Ratio_blank_Sample > Cutoff, 1, 0 )

# Calculating the number of background features and features present
print(paste("Total no.of features:",nrow(Avg_ft)))
print(paste("No.of Background or noise features:",sum(Avg_ft$`Bg_bin` ==1,na.rm = T)))
print(paste("No.of features after excluding noise:",(ncol(Samples) - sum(Avg_ft$`Bg_bin` ==1,na.rm = T))))

# Create df with blank removed features
blk_rem_1 <- merge(as.data.frame(t(Samples)), Avg_ft, by=0) %>%
  filter(Bg_bin == 0) %>% #picking only the features
  select(-c(Avg_blank,Avg_samples,Ratio_blank_Sample,Bg_bin)) %>% #removing the last 4 columns
  column_to_rownames(var="Row.names") 

# Transpose blank removed features df
blk_rem <- as.data.frame(t(blk_rem_1))
  
# Perform blank subtraction for features that remain but have non-zero peak area values in the blanks. 
# The maximum peak area from the blanks is subtracted from all non-blank samples. 
# Add Max_blank column to Avg_ft
Avg_ft$Max_blank <- apply(t(Blank), 1, max, na.rm = TRUE)

# Correctly match feature names using row names from Avg_ft
Blank_max_filtered <- Avg_ft[names(blk_rem), "Max_blank", drop = FALSE]$Max_blank

# Subtract Max blank values from Samples, ensuring no negative values
blk_rem_sub <- sweep(blk_rem, 2, Blank_max_filtered, FUN = "-")
blk_rem_sub[blk_rem_sub < 0] <- 0  # Set negative values to zero

# Determine the limit of detection (LOD)
Cutoff_LOD <- round(min(blk_rem[blk_rem > 0], na.rm = TRUE))
print(paste0("The limit of detection (LOD) is: ", Cutoff_LOD))

# Apply LOD correction
blk_rem_sub[blk_rem_sub > 0 & blk_rem_sub < Cutoff_LOD] <- 0

# View the final processed dataset
print(dim(blk_rem_sub))

# Ensure the final dataframe is formatted correctly
blk_rem <- as.data.frame(blk_rem)
blk_rem_sub <- as.data.frame(blk_rem_sub)

# Review feature tables without blanks
dim(blk_rem)
dim(blk_rem_sub)

# Review metadata without the blanks info
dim(md_Samples)

# Define output file path & write the .csv file
output_path <- file.path(Directory, "02_Output_DataProducts", paste0('S3_Feature_PeakAreas_BlkRemSub_with_cutoff',Cutoff,'.csv'))
write.csv(t(blk_rem_sub), output_path, row.names = TRUE)

##### 2.3 SELECT DATA for ANALYSIS #####

# Select quant table
blk_rem_use <- blk_rem_sub
dim(blk_rem_use)

# Select and format metadata
md_Samples_use <- md_Samples
md_Samples_use$Sample <- rownames(md_Samples_use)
dim(md_Samples_use)

##### 2.4 SUBSET XENOBIOTIC FEATURES USING ANNOTATIONS and MOLECULAR NETWORKS #####

# Create a copy of the original feature group table
features_use <- df_csv_xenobiotic

# Define the group columns
group_cols <- c(paste0("Group", 0:2), "SubCategory")

# Find the rows where cosine is less than 0.7 and the MZerror is greater than 10 - the annotations will not be used from these features
rows_to_blank <- which(!is.na(features_use$Cosine) & features_use$Cosine < 0.7 | !is.na(features_use$Mzerror) & features_use$Mzerror >= 10)

# Check out which rows will be set to blank
features_use[rows_to_blank, c("FeatureID", "Cosine", "Mzerror")]

# Set the group column values to blank if a feature's cosine or Mzerror values do not meet the constraints
features_use[rows_to_blank, group_cols] <- ""

# View result
dim(features_use)
names(features_use)
str(features_use[group_cols])
length(rows_to_blank)

# Filter the features based on if they have an annotation category
features_filtered <- features_use %>%
  filter(Group0 %in% c("a"))

# Extract FeatureIDs present in the feature table after blank removal and subtraction, replicate averaging and scaling
blk_feature_ids <- gsub(".*X([0-9]+_[a-zA-Z]+).*$", "\\1", colnames(blk_rem_use))

# Filter df_csv_edges to keep only rows where BOTH nodes are in blk_rem
df_csv_edges_use <- df_csv_edges %>%
  filter(node1 %in% blk_feature_ids & node2 %in% blk_feature_ids)
dim(df_csv_edges_use)

# Initialize all features df
all_features_df <- data.frame(
  FeatureID = character(),
  Parent = character(),
  Source = character(),
  stringsAsFactors = FALSE
)

# Add filtered features with their subsets
all_features_df <- features_filtered %>%
  select(FeatureID, Group0) %>%
  rename(Source = Group0) %>%
  mutate(Parent = FeatureID)

# For each xenobiotic feature, get first neighbors
for (feature in features_filtered$FeatureID) {

  connected_node2 <- df_csv_edges_use %>%
    filter(node1 == feature) %>%
    select(node2) %>%
    rename(FeatureID = node2) %>%
    mutate(Source = "MN", Parent = feature)

  connected_node1 <- df_csv_edges_use %>%
    filter(node2 == feature) %>%
    select(node1) %>%
    rename(FeatureID = node1) %>%
    mutate(Source = "MN", Parent = feature)

  connected_features <- bind_rows(connected_node1, connected_node2)
  all_features_df <- bind_rows(all_features_df, connected_features)
}

# Collect all IDs found so far
combined_feature_ids <- all_features_df$FeatureID

# Create new df to use
all_features_df_1 <- all_features_df

### Use igraph to add nth MN and cluster info ###
# Ensure only 2 columns are passed to define edges
edge_df <- df_csv_edges_use[, c("node1", "node2")]

# Build the graph cleanly
g <- graph_from_data_frame(edge_df, directed = FALSE)

# Check vertex names
head(V(g)$name)

# Get xenobiotic features from the combined df
xenobiotic_ids <- all_features_df_1 %>%
  filter(Source %in% c("a")) %>%
  pull(FeatureID) %>%
  intersect(V(g)$name)  # make sure they exist in graph

# Get component membership
components_info <- components(g)
node_to_component <- components_info$membership

# Add cluster ID to all_features_df_1
all_features_df_1$ClusterID <- node_to_component[match(all_features_df_1$FeatureID, names(node_to_component))]

# Identify which clusters contain xenobiotic features
xenobiotic_clusters <- unique(all_features_df_1$ClusterID[all_features_df_1$Source == "a"])

# Now find nodes in those xenobiotic clusters
nodes_in_xenobiotic_clusters <- names(node_to_component)[node_to_component %in% xenobiotic_clusters]

# Find nth MN: in xenobiotic clusters, but not xenobiotic or 1st neighbors
nth_mn <- setdiff(nodes_in_xenobiotic_clusters, all_features_df_1$FeatureID[all_features_df_1$Source %in% c("a", "MN")])

# Add these to the dataframe
nth_mn_df <- data.frame(
  FeatureID = nth_mn,
  Source = "nthMN",
  Parent = nth_mn
)

# Merge all features together and merge with group info.
all_features_df_2 <- bind_rows(all_features_df_1, nth_mn_df)

# Collect all features into one df with their group info
combined_df <- full_join(all_features_df_2,
                         features_use[, c("FeatureID", "Feature", "SubCategory", "Group0", "Group1", "Group2", "Mode")],
                         by = "FeatureID")

# Create a lookup table for parent feature group info
parent_groups <- features_use %>%
  select(FeatureID, SubCategory, Group0, Group1, Group2) %>%
  rename(Parent = FeatureID,
         SubCategory_parent = SubCategory,
         Group0_parent = Group0,
         Group1_parent = Group1,
         Group2_parent = Group2)

# Join parent info to combined_df
combined_df_joined <- combined_df %>%
  left_join(parent_groups, by = "Parent")

# Replace Group0–2 values with parent values if Source == "MN"
all_features_df_use <- combined_df_joined %>%
  mutate(
    SubCategory = ifelse(Source == "MN" & !is.na(SubCategory_parent), SubCategory_parent, SubCategory),
    Group0 = ifelse(Source == "MN" & !is.na(Group0_parent), Group0_parent, Group0),
    Group1 = ifelse(Source == "MN" & !is.na(Group1_parent), Group1_parent, Group1),
    Group2 = ifelse(Source == "MN" & !is.na(Group2_parent), Group2_parent, Group2)
  ) %>%
  select(-Group0_parent, -Group1_parent, -Group2_parent)

# Add cluster ID to all_features_df_use
all_features_df_use$ClusterID <- node_to_component[match(all_features_df_use$FeatureID, names(node_to_component))]

# Label all others explicitly as "unknown" if still NA
all_features_df_use$Group0[all_features_df_use$Group0 == ""] <- "unknown"
all_features_df_use$Group1[all_features_df_use$Group1 == ""] <- "unknown"
all_features_df_use$Group2[all_features_df_use$Group2 == ""] <- "unknown"
all_features_df_use$Group1[all_features_df_use$Group0 == "biological"] <- "biological compound"
all_features_df_use$Group2[all_features_df_use$Group0 == "biological"] <- "biological compound"

# Add priority
print(unique(all_features_df_use$Source))
all_features_df_use <- all_features_df_use %>%
  mutate(Priority = case_when(
    Source == "a" ~ 1,
    Source == "MN" ~ 2,
    Source == "nthMN" ~ 3
  )) %>%
  arrange(FeatureID, Priority) %>%
  distinct(FeatureID, .keep_all = TRUE) %>%
  select(-Priority)

# Count the occurrences of each bin using table()
summary <- table(all_features_df_use$Source)
print(summary)

# For Mode = positive
summary_pos <- table(all_features_df_use$Source[all_features_df_use$Mode == "pos"])
print(summary_pos)

# For Mode = negative
summary_neg <- table(all_features_df_use$Source[all_features_df_use$Mode == "neg"])
print(summary_neg)

##### 3.  DATA ANALYSIS #####
##### 3.1 DATA ANALYSIS SET-UP #####

# Prep inputs
group_subset <- c("a")
source_subset <- c("nthMN")

# Create copy of quant data frame
feature_table <- blk_rem_use

#Create copy of metadata data frame
metadata_table <- md_Samples_use
metadata_table$Sample <- rownames(metadata_table)

# Compile data for analysis
feature_table_long <- feature_table %>%
  rownames_to_column(var = "Sample") %>%  # Keep sample names
  pivot_longer(cols = -Sample, names_to = "Feature", values_to = "PeakArea") %>%
  mutate(
    #Mode = ifelse(grepl("_pos", Feature), "pos", "neg"),
    FeatureID = sub(".*X([0-9]+_[a-zA-Z]+).*$", "\\1", Feature)
  ) %>%
  left_join(metadata_table %>% mutate(Sample = rownames(metadata_table)), by = "Sample") %>%
  left_join(all_features_df_use  %>% select(-Feature), by = "FeatureID") %>%
  left_join(df_csv_SIRIUS_formulae, by = "FeatureID")

# PAG sample calculation
pag_summary <- feature_table_long %>%
  group_by(Mode, Sample) %>%
  summarise(
    TotalPeakAreaPAG = sum(PeakArea[SubCategory == "PAG"], na.rm = TRUE),
    TotalPeakAreaAll = sum(PeakArea, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    PAG_Percentage = (TotalPeakAreaPAG / TotalPeakAreaAll) * 100
  )

# Write results
# Define output file path & write the csv file
output_path <- file.path(Directory, "02_Output_DataProducts", paste0("S9_Samples_PAGPercentage.csv"))
write.csv(pag_summary, output_path, row.names = FALSE)

# Identify samples with >=10% PAG in POS mode
samples_with_PAG <- pag_summary %>%
  filter(Mode == "pos", round(PAG_Percentage, 0) >= 10) %>%
  pull(Sample)

# Merge PAG percentages into main table
feature_table_long_pag <- feature_table_long %>%
  left_join(
    pag_summary %>% select(Sample, Mode, PAG_Percentage),
    by = c("Sample", "Mode")
  )

# Label samples based on POS mode only
feature_table_long <- feature_table_long_pag %>%
  mutate(
    ATTRIBUTE_RiverGroup_PAG = ifelse(
      Sample %in% samples_with_PAG,
      paste0(ATTRIBUTE_RiverGroup, "_PAG"),
      ATTRIBUTE_RiverGroup
    )
  )

# Calculate percent contribution per sample
feature_table_with_pct <- feature_table_long %>%
  group_by(Sample) %>%
  mutate(TotalPeakArea = sum(PeakArea, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(PercentOfSample = ifelse(TotalPeakArea > 0, (PeakArea / TotalPeakArea) * 100, 0))

# Summarize by FeatureID and Mode
feature_stats <- feature_table_with_pct %>%
  group_by(FeatureID) %>%
  summarise(
    MeanPeakArea = mean(PeakArea, na.rm = TRUE),
    SDPeakArea = sd(PeakArea, na.rm = TRUE),
    NonZeroCount = sum(PeakArea > 0, na.rm = TRUE),
    MeanPercent = mean(PercentOfSample, na.rm = TRUE),
    SDPercent = sd(PercentOfSample, na.rm = TRUE),
    .groups = "drop"
  )

# Join back to original long table
feature_table_long <- feature_table_long %>%
  left_join(feature_stats, by = c("FeatureID"))

# Make summary table for all features in each mode
feature_summary_table <- feature_table_long %>%
  select(FeatureID, Feature, 
         Mode,Source, Parent, Group0,NPC.pathway, NPC.pathway.Probability,
         NonZeroCount, MeanPeakArea, SDPeakArea, MeanPercent, SDPercent) %>%
  distinct


# Write results
# Define output file path & write the csv file
output_path <- file.path(Directory, "02_Output_DataProducts", paste0("S7b_Features_Info.csv"))
write.csv(feature_summary_table, output_path, row.names = FALSE)

# Merge PAG percentages into metadata table
metadata_table_pag <- metadata_table %>%
  left_join(
    pag_summary %>%
      select(Sample, Mode, PAG_Percentage),
    by = "Sample"
  )

# Create new ATTRIBUTE_RiverGroup_PAG column
metadata_table <- metadata_table_pag %>%
  mutate(
    ATTRIBUTE_RiverGroup_PAG = ifelse(
      Sample %in% samples_with_PAG,
      paste0(ATTRIBUTE_RiverGroup, "_PAG"),
      ATTRIBUTE_RiverGroup
    )
  )
             
              
##### 3.2 XENOBIOTIC PERCENTAGES of DOM (SI & FIGURE 1) #####                                              

# Compute percentage of xenobiotic features per sample and mode
xeno_summary <- feature_table_long %>%
  filter(!is.na(Mode)) %>%
  group_by(Mode, Sample, Source, Group0) %>%
  summarise(TotalFeaturesXeno = n(), 
            TotalPeakAreaXeno = sum(PeakArea, na.rm = TRUE),  .groups = "drop") %>%
  group_by(Mode, Sample) %>%
  mutate(Percentage = (TotalPeakAreaXeno / sum(TotalPeakAreaXeno, na.rm = TRUE)) * 100) %>%
  ungroup()

# Sum xenobiotic contribution to samples for each feature confidence category
xeno_summary_use <- xeno_summary %>%
  group_by(Sample, Mode) %>%
  summarise(
    Category_Anno = sum(Percentage[Source == "a" ], na.rm = TRUE),
    Category_Anno1 = sum(Percentage[Source %in% c("a", "MN")], na.rm = TRUE),
    Category_Anno1n = sum(Percentage[Source %in% c("a", "MN", "nthMN")], na.rm = TRUE),
    .groups = "drop"
  )

# Average across ionization modes
xeno_summary_both <- xeno_summary_use %>%
  group_by(Sample) %>%
  summarise(
    Mode = "both",
    Category_Anno = mean(Category_Anno, na.rm = TRUE),
    Category_Anno1 = mean(Category_Anno1, na.rm = TRUE),
    Category_Anno1n = mean(Category_Anno1n, na.rm = TRUE),
    .groups = "drop"
  )

# Combine the averaged ("both") data with the original per-mode data
xeno_summary_final <- bind_rows(xeno_summary_use, xeno_summary_both)

# Summarize the xenobiotic percentage across all samples, by annotation category
avg_xeno_contrib <- xeno_summary_final %>%
  group_by(Mode) %>%
  summarise(
    avg_Category_Anno   = mean(Category_Anno, na.rm = TRUE),
    avg_Category_Anno1  = mean(Category_Anno1, na.rm = TRUE),
    avg_Category_Anno1n = mean(Category_Anno1n, na.rm = TRUE),
    .groups = "drop"
  )

# Reshape to long format for plotting
xeno_long_1 <- xeno_summary_final %>%
  pivot_longer(cols = starts_with("Category"),
               names_to = "Category",
               values_to = "Percentage")

# Clean category names
xeno_long_1$Category <- recode(xeno_long_1$Category,
                             "Category_Anno" = "Anno.",
                             "Category_Anno1" = "Anno. + 1st",
                             "Category_Anno1n" = "Anno. + 1st + nth")

# Convert to factor with ordered levels
xeno_long_1$Category <- factor(xeno_long_1$Category, levels = c(
  "Anno.", "Anno. + 1st", "Anno. + 1st + nth"))

# Join to specific columns from the metadata
xeno_long_2 <- xeno_long_1 %>%
  left_join(
    metadata_table %>%
      select(Sample, ATTRIBUTE_RiverGroup) %>%
      distinct(Sample, .keep_all = TRUE),  # keep one row per sample
    by = "Sample"
  )

# Create new ATTRIBUTE_RiverGroup_PAG column
xeno_long <- xeno_long_2 %>%
  mutate(
    ATTRIBUTE_RiverGroup_PAG = ifelse(
      Sample %in% samples_with_PAG,
      paste0(ATTRIBUTE_RiverGroup, "_PAG"),
      ATTRIBUTE_RiverGroup
    )
  )

# Create one plot per ionization mode
# Positive mode
pos <- ggplot(xeno_long[xeno_long$Category %in% c("Anno.", "Anno. + 1st", "Anno. + 1st + nth") & xeno_long$Mode == "pos",], aes(x = Category, y = Percentage, fill = Category)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0,
               justification = -0.4, point_colour = NA, alpha = 0.4) +  # violin half 
  geom_jitter(aes(color = ATTRIBUTE_RiverGroup_PAG), width = 0.15, size = 1, shape = 16, alpha = 0.7) +  # points by River
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_nudge(x = 0.22), alpha = 1) +
  scale_fill_manual(values = c("darkgray", "darkgray", "darkgray")) +
  scale_color_manual(values = c(
    "Neckar" ="maroon3",
    "Neckar_PAG" = "goldenrod",
    "Spree" = "darkblue"
  )) +
  labs(x = "Xenobiotic feature confidence category", 
       y = "Samples' Xenobiotic peak area percentage (%)",
       title = "Positive Mode") +
  theme_bw(base_size = 7) +
  theme(legend.position = "right",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7, face = "bold"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, face = "bold"),
        #axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = "none")  # Remove fill legend
pos

# Negative mode
neg <- ggplot(xeno_long[xeno_long$Category %in% c("Anno.", "Anno. + 1st", "Anno. + 1st + nth") & xeno_long$Mode == "neg",], aes(x = Category, y = Percentage, fill = Category)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0,
               justification = -0.4, point_colour = NA, alpha = 0.4) +  # violin half 
  geom_jitter(aes(color = ATTRIBUTE_RiverGroup_PAG), width = 0.15, size = 1, shape = 16, alpha = 0.7) +  # points by River
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_nudge(x = 0.22), alpha = 1) +
  scale_fill_manual(values = c("darkgrey", "darkgray", "darkgray")) +
  scale_color_manual(values = c(
    "Neckar" ="maroon3",
    "Neckar_PAG" = "goldenrod",
    "Spree" = "darkblue"
  )) +
  labs(x = "Xenobiotic feature confidence category", 
       y = "Samples' Xenobiotic peak area percentage (%)",
       title = "Negative Mode") +
  theme_bw(base_size = 7) +
  theme(legend.position = "right",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7, face = "bold"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, face = "bold"),
        #axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = "none")  # Remove fill legend
neg

# Positive and Negative modes - averaged
both <- ggplot(xeno_long[xeno_long$Category %in% c("Anno.", "Anno. + 1st", "Anno. + 1st + nth") & xeno_long$Mode == "both",], aes(x = Category, y = Percentage, fill = Category)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0,
               justification = -0.4, point_colour = NA, alpha = 0.4) +  # violin half 
  geom_jitter(aes(color = ATTRIBUTE_RiverGroup_PAG), width = 0.15, size = 1, shape = 16, alpha = 0.7) +  # points by River
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_nudge(x = 0.22), alpha = 1) +
  scale_fill_manual(values = c("darkgray", "darkgray", "darkgray")) +
  scale_color_manual(values = c(
    "Neckar" ="maroon3",
    "Neckar_PAG" = "goldenrod",
    "Spree" = "darkblue"
  )) +
  labs(x = "Xenobiotic feature confidence category", 
       y = "Samples' Xenobiotic peak area percentage (%)",
       title = "Positive & Negative Modes - Averaged",
       color = "Samples") +
  theme_bw(base_size = 7) +
  theme(legend.position = "right",
        #axis.text.x = element_blank(),
        axis.text = element_text(size = 7, face = "bold"),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = 7, face = "bold"),
        #axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = "none")  # Remove fill legend
both

# Set pdf dimensions
h <- 3
w <- 4

# Save plot
output_path <- file.path(Directory, "03_Output_Figures", paste0("Fig1B_RainPlots_XenobioticFeatureConfidence_pos.pdf"))
pdf(file = output_path, height = h, width = w)
print(pos)
dev.off()

# Save plot
output_path <- file.path(Directory, "03_Output_Figures", paste0("Fig1B_RainPlots_XenobioticFeatureConfidence_neg.pdf"))
pdf(file = output_path, height = h, width = w)
print(neg)
dev.off()

# Save plot
output_path <- file.path(Directory, "03_Output_Figures", paste0("Fig1B_RainPlots_XenobioticFeatureConfidence_both.pdf"))
pdf(file = output_path, height = h, width = w)
print(both)
dev.off()

# Write results
# Define output file path & write the csv file
output_path <- file.path(Directory, "02_Output_DataProducts", paste0("S8_Samples_XenobioticPercentage.csv"))
write.csv(xeno_long, output_path, row.names = FALSE)


##### 3.3 NPC PERCENTAGES of DOM (SI & FIGURE 2) #####

# Calculate the sample peak area values
xenobiotic_feature_plot_data <- feature_table_long %>%
  filter(!is.na(Mode))  %>%
  group_by(Mode, Sample, NPC.pathway) %>%
  summarise(
    TotalPeakArea = sum(PeakArea, na.rm = TRUE),
    # Total peak area of Xenobiotic contaminants in each feature category
    XenobioticPeakArea = sum(PeakArea[Group0 %in% group_subset & !Source %in% source_subset], na.rm = TRUE),
    Sample_Type = unique(ATTRIBUTE_Type_2),
    RiverGroup = unique(ATTRIBUTE_RiverGroup_PAG),
    .groups = "drop"
  ) %>%
  # Calculate the percentage of xenobiotic contaminants in the sample
  group_by(Sample, Mode) %>%
  mutate(
    XenobioticPercentSample = (sum(XenobioticPeakArea) / sum(TotalPeakArea)) * 100,
    # Calculate the percentage of feature group contribution to the sample
    FeatureGroupPercentSample = (TotalPeakArea / sum(TotalPeakArea)) * 100,
    # Calculate the percentage of xenobiotic contribution to each feature group
    XenobioticPercentGroup = (XenobioticPeakArea / TotalPeakArea) * 100
  ) %>%
  ungroup()
  
                                            
# Sample-level totals
xenobiotic_percent_sample <- feature_table_long %>%
  group_by(Sample, Mode) %>%
  summarise(
    TotalPeakArea = sum(PeakArea, na.rm = TRUE),
    XenobioticPeakArea = sum(PeakArea[Group0 %in% group_subset & !Source %in% source_subset], na.rm = TRUE),
    XenobioticPercentSample = (XenobioticPeakArea / TotalPeakArea) * 100,
    Sample_Type = unique(ATTRIBUTE_Type_2),
    RiverGroup = unique(ATTRIBUTE_RiverGroup_PAG),
    .groups = "drop"
  )

# Group-level xenobiotic peak areas per sample
xenobiotic_contrib_by_group <- feature_table_long %>%
  filter(Group0 %in% group_subset & !Source %in% source_subset) %>%  # keep only xenobiotic features
  group_by(Sample, Mode, NPC.pathway) %>%
  summarise(
    GroupXenobioticPeakArea = sum(PeakArea, na.rm = TRUE),
    .groups = "drop"
  )

# Join to get percentage of sample total
xenobiotic_feature_plot_data_XenoOnly <- xenobiotic_contrib_by_group %>%
  left_join(xenobiotic_percent_sample, by = c("Sample", "Mode")) %>%
  mutate(
    GroupXenobioticPercentOfSample = (GroupXenobioticPeakArea / TotalPeakArea) * 100
  )

xenobiotic_feature_plot_data <- left_join(xenobiotic_feature_plot_data, 
                                      xenobiotic_feature_plot_data_XenoOnly[,c("Sample", "Mode", "NPC.pathway", "GroupXenobioticPercentOfSample")], 
                                      by = c("Sample", "Mode", "NPC.pathway"))

# Write results
# Define output file path & write the csv file
output_path <- file.path(Directory, "02_Output_DataProducts", paste0("S10_Samples_NPCPercentage.csv"))
write.csv(xenobiotic_feature_plot_data, output_path, row.names = FALSE)

# Initialize an empty list to store the plots
NPC_plots <- list()

# Set colors for catchments (Neckar, Neckar_PAG, Spree)
catchments <- c("maroon3", "goldenrod", 'darkblue')

# Plot the percentages for each NPC pathway
for (m in unique(xenobiotic_feature_plot_data$Mode)){
  print(paste(m))
  data_1 <-  xenobiotic_feature_plot_data[xenobiotic_feature_plot_data$Mode == m,]
  
  for (npc in unique(data_1$NPC.pathway)) {
    npc_lab <- npc
    if (is.na(npc_lab) || npc_lab == "") npc_lab <- "Unassigned"
    print(paste(npc_lab))
    data <- data_1[data_1$NPC.pathway == npc,]
    box <- ggplot(data, aes(x = XenobioticPercentSample, y = FeatureGroupPercentSample, shape = RiverGroup, color = RiverGroup)) +
       
      # Hollow points: FeatureGroupPercentSample
      geom_point(data = data,
                 aes(x = XenobioticPercentSample,
                     y = FeatureGroupPercentSample,
                     color = RiverGroup,
                     shape = Sample_Type),
                 fill = NA, size = 3, stroke = 1) +

      # Filled points: GroupXenobioticPercentOfSample
      geom_point(data = data,
                 aes(x = XenobioticPercentSample,
                     y = GroupXenobioticPercentOfSample,
                     fill = RiverGroup,
                     shape = Sample_Type),
                 size = 3, stroke = 0.5, color = "gray75")  +
      scale_color_manual(values = catchments) +  # Update to match RiverGroup
      scale_fill_manual(values = catchments) +
      scale_shape_manual(values = c(21, 21, 21, 21, 21)) +
      labs(title = npc_lab, 
           x = "Xeno Perc", 
           y = "Percentage", 
           fill = "River") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 14),
        panel.grid.minor = element_blank()
      )
    
    NPC_plots[[paste(m, npc,"_perc", sep = "")]] <- box
  }}

# Plot scatter plots - negative
Boxes_long_perc_n <- plot_grid(plotlist = NPC_plots[c(2:8,1)], ncol = 1, nrow = 8, rel_widths = c(1), rel_heights = c(1), label_size = 12, label_y = 0.9,
                               align = "v")
# Add labels
Boxes_long_perc_n_final <-
  ggdraw() +
  draw_label("Negative Mode", fontface = "bold", size = 16, x = 0.5, y = 0.98, hjust = 0.5) +
  draw_label("Samples' Xenobiotic Percentage (%)", x = 0.5, y = 0.02, hjust = 0.5, size = 14) +
  draw_label("Samples' NPC Pathway Percentage (%)", x = 0.04, y = 0.5, angle = 90, hjust = 0.5, size = 14) +
  draw_plot(Boxes_long_perc_n, x = 0.07, y = 0.05, width = 0.9, height = 0.9)

# Plot scatter plots - positive
Boxes_long_perc_p <- plot_grid(plotlist = NPC_plots[c(10:16,9)], ncol = 1, nrow = 8,rel_widths = c(1), rel_heights = c(1), label_size = 12, label_y = 0.9,
                               align = "v")

# Add labels
Boxes_long_perc_p_final <-
  ggdraw() +
  draw_label("Positive Mode", fontface = "bold", size = 16, x = 0.5, y = 0.98, hjust = 0.5) +
  draw_label("Samples' Xenobiotic Percentage (%)", x = 0.5, y = 0.02, hjust = 0.5, size = 14) +
  draw_label("Samples' NPC Pathway Percentage (%)", x = 0.04, y = 0.5, angle = 90, hjust = 0.5, size = 14) +
  draw_plot(Boxes_long_perc_n, x = 0.07, y = 0.05, width = 0.9, height = 0.9)

# Set pdf dimensions
h <- 16
w <- 5.4

# Plot the entire plot
pdf(file = file.path(Directory, "03_Output_Figures", paste("FigS4A_ScatterPlots_NPC_and_XenobioticFeatures_Anno1_pos.pdf")), height = h, width = w)
plot_grid(Boxes_long_perc_p_final, ncol = 1, nrow = 1, rel_widths = c(1,1), rel_heights = c(1,5), label_size = 11, label_y = 0.9)
dev.off()         

# Plot the entire plot
pdf(file = file.path(Directory, "03_Output_Figures", paste("FigS4B_ScatterPlots_NPC_and_XenobioticFeatures_Anno1_neg.pdf")), height = h, width = w)
plot_grid(Boxes_long_perc_n_final, ncol = 1, nrow = 1, rel_widths = c(1,1), rel_heights = c(1,4), label_size = 11, label_y = 0.9)
dev.off()     


### Plot FIGURE 2 ###

# Build plotting datasets
# Data - Pop density per sample+mode
pop_density_df <- feature_table_long %>%
  filter(!is.na(Mode)) %>%
  group_by(Sample, Mode) %>%
  summarise(
    PopDensity = dplyr::first(na.omit(ATTRIBUTE_PopulationDensity)),
    .groups = "drop"
  ) %>%
  left_join(
    xenobiotic_percent_sample %>% select(Sample, Mode, XenobioticPercentSample, Sample_Type, RiverGroup),
    by = c("Sample", "Mode")
  ) %>%
  mutate(PopDensity = as.numeric(PopDensity)) %>%
  mutate(PointType = "Samples")

# Exclude the two high pop-density points (Wuhle) from the trendline
pop_density_smooth_df <- pop_density_df %>%
  filter(!is.na(PopDensity), PopDensity <= 700,
         !is.na(XenobioticPercentSample))

# Helper: robust pathway matching
is_alk  <- function(x) grepl("alkaloid", tolower(x))
is_terp <- function(x) grepl("terpen",  tolower(x))

# Data - Alkaloids (ALL features vs XENO-only)
alk_all <- xenobiotic_feature_plot_data %>%
  filter(!is.na(NPC.pathway), is_alk(NPC.pathway)) %>%
  group_by(Sample, Mode) %>%
  summarise(
    XenobioticPercentSample = dplyr::first(XenobioticPercentSample),
    Y = sum(FeatureGroupPercentSample, na.rm = TRUE),
    RiverGroup = dplyr::first(RiverGroup),
    .groups = "drop"
  ) %>%
  mutate(PointType = "All features\nin Sample")

alk_xeno <- xenobiotic_feature_plot_data %>%
  filter(!is.na(NPC.pathway), is_alk(NPC.pathway)) %>%
  group_by(Sample, Mode) %>%
  summarise(
    XenobioticPercentSample = dplyr::first(XenobioticPercentSample),
    Y = sum(GroupXenobioticPercentOfSample, na.rm = TRUE),
    RiverGroup = dplyr::first(RiverGroup),
    .groups = "drop"
  ) %>%
  mutate(PointType = "Xenobiotic features\nin Sample")

alk_df_long <- bind_rows(alk_all, alk_xeno)

# Data - Terpenoids (ALL features vs XENO-only)
terp_all <- xenobiotic_feature_plot_data %>%
  filter(!is.na(NPC.pathway), is_terp(NPC.pathway)) %>%
  group_by(Sample, Mode) %>%
  summarise(
    XenobioticPercentSample = dplyr::first(XenobioticPercentSample),
    Y = sum(FeatureGroupPercentSample, na.rm = TRUE),
    RiverGroup = dplyr::first(RiverGroup),
    .groups = "drop"
  ) %>%
  mutate(PointType = "All features\nin Sample")

terp_xeno <- xenobiotic_feature_plot_data %>%
  filter(!is.na(NPC.pathway), is_terp(NPC.pathway)) %>%
  group_by(Sample, Mode) %>%
  summarise(
    XenobioticPercentSample = dplyr::first(XenobioticPercentSample),
    Y = sum(GroupXenobioticPercentOfSample, na.rm = TRUE),
    RiverGroup = dplyr::first(RiverGroup),
    .groups = "drop"
  ) %>%
  mutate(PointType = "Xenobiotic features\nin Sample")

terp_df_long <- bind_rows(terp_all, terp_xeno)

# Colors (river groups)
catchments <- c("maroon3", "goldenrod", "darkblue")
names(catchments) <- c("Neckar", "Neckar_PAG", "Spree") 

# Shapes (point types)
shape_map <- c(
  "Samples" = 17,
  "All features\nin Sample" = 1,
  "Xenobiotic features\nin Sample" = 16
)


# Create plot components
theme_base <- function(show_x = FALSE) {
  theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = if (show_x) element_text(size = 8) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.text.y  = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 8)
    )
}

# Pop density: points colored by RiverGroup, SINGLE pooled trendline (using filtered data)
make_pop_plot_pooled <- function(df_points, df_smooth, mode_value, show_x = FALSE) {
  pts <- df_points %>%
    filter(Mode == mode_value) %>%
    filter(!is.na(PopDensity), !is.na(XenobioticPercentSample))
  
  sm  <- df_smooth %>%
    filter(Mode == mode_value)
  
  ggplot(pts, aes(x = XenobioticPercentSample, y = PopDensity)) +
    geom_point(aes(color = RiverGroup, shape = PointType), size = 3) +
    geom_smooth(
      inherit.aes = FALSE,
      data = sm,
      aes(x = XenobioticPercentSample, y = PopDensity),
      method = "lm", se = FALSE, linewidth = 0.6, color = "black"
    ) +
    scale_color_manual(values = catchments, name = "River group") +
    scale_shape_manual(values = shape_map, name = "Point type") +
    coord_cartesian(ylim = c(0, 700)) +
    ylab("Upstream catchment\npopulation density") +
    theme_base(show_x = show_x)
}

# Alk/Terp: points colored by RiverGroup, SINGLE pooled trendline per PointType
make_pct_plot_pooled <- function(df_long, ylab_txt, mode_value, show_x = FALSE) {
  df2 <- df_long %>%
    filter(Mode == mode_value) %>%
    filter(!is.na(Y), !is.na(XenobioticPercentSample))
  
  ggplot(df2, aes(x = XenobioticPercentSample, y = Y)) +
    geom_point(aes(color = RiverGroup, shape = PointType), size = 3) +
    
    # pooled trend for ALL features (solid)
    geom_smooth(
      data = df2 %>% filter(PointType == "All features\nin Sample"),
      inherit.aes = FALSE,
      aes(x = XenobioticPercentSample, y = Y),
      method = "lm", se = FALSE, linewidth = 0.6, color = "black"
    ) +
    
    # pooled trend for XENO-only (dashed)
    geom_smooth(
      data = df2 %>% filter(PointType == "Xenobiotic features\nin Sample"),
      inherit.aes = FALSE,
      aes(x = XenobioticPercentSample, y = Y),
      method = "lm", se = FALSE, linewidth = 0.6, linetype = "dashed", color = "black"
    ) +
    
    scale_color_manual(values = catchments, name = "River group") +
    scale_shape_manual(values = shape_map, name = "Point type") +
    ylab(ylab_txt) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_base(show_x = show_x)
}

# Build the 2-col x 3-row panel plot (Pos left, Neg right)
modes <- unique(xenobiotic_feature_plot_data$Mode)
mode_pos <- modes[grepl("pos", tolower(modes))][1]
mode_neg <- modes[grepl("neg", tolower(modes))][1]

# Top row: Alkaloids
p_alk_pos  <- make_pct_plot_pooled(alk_df_long,  "N-containing [Alkaloid]\nfeatures peak area %",  mode_pos, show_x = FALSE)
p_alk_neg  <- make_pct_plot_pooled(alk_df_long,  "N-containing [Alkaloid]\nfeatures peak area %",  mode_neg, show_x = FALSE) +
  theme(axis.title.y = element_blank())

# Middle row: Terpenoids
p_terp_pos <- make_pct_plot_pooled(terp_df_long, "Terpenoid\nfeatures peak area %", mode_pos, show_x = TRUE)
p_terp_neg <- make_pct_plot_pooled(terp_df_long, "Terpenoid\nfeatures peak area %", mode_neg, show_x = TRUE) +
  theme(axis.title.y = element_blank())

# Bottom row: Pop density
p_pop_pos  <- make_pop_plot_pooled(pop_density_df, pop_density_smooth_df, mode_pos, show_x = FALSE)
p_pop_neg  <- make_pop_plot_pooled(pop_density_df, pop_density_smooth_df, mode_neg, show_x = FALSE) +
  theme(axis.title.y = element_blank())

# Extract the legends to include in the plot
get_legend_grob <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(idx) == 0) idx <- grep("guide", sapply(g$grobs, function(x) x$name))
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

# SHAPE legend
legend_shape_plot <- ggplot(
  data.frame(PointType = factor(names(shape_map), levels = names(shape_map)),
             x = seq_along(names(shape_map)), y = 1),
  aes(x, y, shape = PointType)) +
  geom_point(size = 3, color = "black") +
  scale_shape_manual(values = shape_map, name = "Point type", drop = FALSE) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 8))

legend_shape <- get_legend_grob(legend_shape_plot)

# COLOR legend
legend_color_plot <- ggplot(
  data.frame(RiverGroup = factor(names(catchments), levels = names(catchments)),
             x = seq_along(names(catchments)), y = 1),
  aes(x, y, color = RiverGroup)) +
  geom_point(size = 3) +
  scale_color_manual(values = catchments, name = "River group", drop = FALSE) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 8))

legend_color <- get_legend_grob(legend_color_plot)

legend <- plot_grid(
  ggdraw() + draw_grob(legend_shape),
  ggdraw() + draw_grob(legend_color),
  ncol = 1,
  rel_heights = c(1, 1)
)

# Assemble the plot with figure letters, column titles, and axes labels
strip_legend <- function(p) p + theme(legend.position = "none")

panel <- plot_grid(
  strip_legend(p_pop_pos),  strip_legend(p_pop_neg),
  strip_legend(p_alk_pos),  strip_legend(p_alk_neg),
  strip_legend(p_terp_pos), strip_legend(p_terp_neg),
  ncol = 2,
  align = "hv",
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 10,
  label_fontface = "bold",
  label_x = 0.07
  
)

panel_labeled <- ggdraw() +
  draw_plot(panel, x = 0.02, y = 0.04, width = 0.96, height = 0.93) +
  draw_label("Positive Mode", fontface = "bold", size = 10, x = 0.29, y = 0.996, hjust = 0.5, vjust = 1) +
  draw_label("Negative Mode", fontface = "bold", size = 10, x = 0.78, y = 0.996, hjust = 0.5, vjust = 1) +
  # duplicated x-axis titles (one per column)
  draw_label("Xenobiotic features peak area %", size = 8, x = 0.29, y = 0.03, hjust = 0.5) +
  draw_label("Xenobiotic features peak area %", size = 8, x = 0.78, y = 0.03, hjust = 0.5)
 

panel_final <- plot_grid(
  panel_labeled,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.09)
)

# Save the plot
h <- 6
w <- 6.5
pdf(
  file = file.path(Directory, "03_Output_Figures",
                   "Fig2_ScatterPlots_PopDen_NPC_and_XenobioticFeatures_Anno1.pdf"),
  height = h, width = w
)
print(panel_final)
dev.off()


### Correlations & Linear Regressions ###

# Functions to handle NAs and unexpected data

safe_cor <- function(x, y, method) {
  tryCatch({
    cor.test(x, y, method = method)
  }, error = function(e) {
    list(estimate = NA, p.value = NA)
  })
}

safe_lm <- function(formula, data) {
  tryCatch({
    model <- lm(formula, data = data)
    s <- summary(model)
    list(
      slope = coef(model)[2],
      intercept = coef(model)[1],
      r2 = s$r.squared,
      p = coef(s)[2, 4]
    )
  }, error = function(e) {
    list(slope = NA, intercept = NA, r2 = NA, p = NA)
  })
}

safe_shapiro <- function(x) {
  tryCatch({
    shapiro.test(x)$p.value
  }, error = function(e) NA)
}

# Run NPC-pathway correlations & linear model
river_groups <- list(
  All = xenobiotic_feature_plot_data,
  Non_PAG = xenobiotic_feature_plot_data[xenobiotic_feature_plot_data$RiverGroup %in% c("Spree", "Neckar"),],
  Spree = xenobiotic_feature_plot_data[xenobiotic_feature_plot_data$RiverGroup == "Spree",],
  Neckar = xenobiotic_feature_plot_data[xenobiotic_feature_plot_data$RiverGroup == "Neckar",],
  Neckar_PAG = xenobiotic_feature_plot_data[xenobiotic_feature_plot_data$RiverGroup == "Neckar_PAG",]
)

npc_results <- list()

for (river_name in names(river_groups)) {
  data_group <- river_groups[[river_name]]
  
  for (npc in unique(data_group$NPC.pathway)) {
    data_npc <- data_group[data_group$NPC.pathway == npc, ]
    
    for (m in unique(data_npc$Mode)) {
      data <- data_npc[data_npc$Mode == m, ]
      
      if (nrow(data) < 3) {
        measures <- c("TotalPeakAreaOfGroup", "FeatureGroupPercentOfSample", "GroupXenobioticPercentOfSample")
        for (meas in measures) {
          npc_results[[paste("NPC", river_name, npc, m, meas, sep = "_")]] <- data.frame(
            Analysis = "NPC",
            River = river_name,
            NPC_pathway = npc,
            Mode = m,
            Measure = meas,
            Shaprio_Xenobiotic_p = NA,
            Shaprio_p = NA,
            Pearson_r = NA, Pearson_p = NA,
            Spearman_rho = NA, Spearman_p = NA,
            Kendall_tau = NA, Kendall_p = NA,
            LM_Slope = NA, LM_Intercept = NA, LM_R2 = NA, LM_p = NA
          )
        }
        next
      }
      
      shapiro_xenobiotic_p <- safe_shapiro(data$XenobioticPercentSample)
      shapiro_pa_p         <- safe_shapiro(data$TotalPeakArea)
      shapiro_perc_p       <- safe_shapiro(data$FeatureGroupPercentSample)
      shapiro_xeno_perc_p  <- safe_shapiro(data$GroupXenobioticPercentOfSample)
      
      # TotalPeakArea (called TotalPeakAreaOfGroup in output)
      pearson_pa  <- safe_cor(data$XenobioticPercentSample, data$TotalPeakArea, "pearson")
      spearman_pa <- safe_cor(data$XenobioticPercentSample, data$TotalPeakArea, "spearman")
      kendall_pa  <- safe_cor(data$XenobioticPercentSample, data$TotalPeakArea, "kendall")
      lm_pa       <- safe_lm(TotalPeakArea ~ XenobioticPercentSample, data)
      
      npc_results[[paste("NPC", river_name, npc, m, "TotalPeakAreaOfGroup", sep = "_")]] <- data.frame(
        Analysis = "NPC",
        River = river_name, NPC_pathway = npc, Mode = m, Measure = "TotalPeakAreaOfGroup",
        Shaprio_Xenobiotic_p = shapiro_xenobiotic_p, Shaprio_p = shapiro_pa_p,
        Pearson_r = pearson_pa$estimate, Pearson_p = pearson_pa$p.value,
        Spearman_rho = spearman_pa$estimate, Spearman_p = spearman_pa$p.value,
        Kendall_tau = kendall_pa$estimate, Kendall_p = kendall_pa$p.value,
        LM_Slope = lm_pa$slope, LM_Intercept = lm_pa$intercept, LM_R2 = lm_pa$r2, LM_p = lm_pa$p
      )
      
      # FeatureGroupPercentSample
      pearson_perc  <- safe_cor(data$XenobioticPercentSample, data$FeatureGroupPercentSample, "pearson")
      spearman_perc <- safe_cor(data$XenobioticPercentSample, data$FeatureGroupPercentSample, "spearman")
      kendall_perc  <- safe_cor(data$XenobioticPercentSample, data$FeatureGroupPercentSample, "kendall")
      lm_perc       <- safe_lm(FeatureGroupPercentSample ~ XenobioticPercentSample, data)
      
      npc_results[[paste("NPC", river_name, npc, m, "FeatureGroupPercentOfSample", sep = "_")]] <- data.frame(
        Analysis = "NPC",
        River = river_name, NPC_pathway = npc, Mode = m, Measure = "FeatureGroupPercentOfSample",
        Shaprio_Xenobiotic_p = shapiro_xenobiotic_p, Shaprio_p = shapiro_perc_p,
        Pearson_r = pearson_perc$estimate, Pearson_p = pearson_perc$p.value,
        Spearman_rho = spearman_perc$estimate, Spearman_p = spearman_perc$p.value,
        Kendall_tau = kendall_perc$estimate, Kendall_p = kendall_perc$p.value,
        LM_Slope = lm_perc$slope, LM_Intercept = lm_perc$intercept, LM_R2 = lm_perc$r2, LM_p = lm_perc$p
      )
      
      # GroupXenobioticPercentOfSample
      pearson_xeno_perc  <- safe_cor(data$XenobioticPercentSample, data$GroupXenobioticPercentOfSample, "pearson")
      spearman_xeno_perc <- safe_cor(data$XenobioticPercentSample, data$GroupXenobioticPercentOfSample, "spearman")
      kendall_xeno_perc  <- safe_cor(data$XenobioticPercentSample, data$GroupXenobioticPercentOfSample, "kendall")
      lm_xeno_perc       <- safe_lm(GroupXenobioticPercentOfSample ~ XenobioticPercentSample, data)
      
      npc_results[[paste("NPC", river_name, npc, m, "GroupXenobioticPercentOfSample", sep = "_")]] <- data.frame(
        Analysis = "NPC",
        River = river_name, NPC_pathway = npc, Mode = m, Measure = "GroupXenobioticPercentOfSample",
        Shaprio_Xenobiotic_p = shapiro_xenobiotic_p, Shaprio_p = shapiro_xeno_perc_p,
        Pearson_r = pearson_xeno_perc$estimate, Pearson_p = pearson_xeno_perc$p.value,
        Spearman_rho = spearman_xeno_perc$estimate, Spearman_p = spearman_xeno_perc$p.value,
        Kendall_tau = kendall_xeno_perc$estimate, Kendall_p = kendall_xeno_perc$p.value,
        LM_Slope = lm_xeno_perc$slope, LM_Intercept = lm_xeno_perc$intercept, LM_R2 = lm_xeno_perc$r2, LM_p = lm_xeno_perc$p
      )
    }
  }
}

correlation_df_npc <- do.call(rbind, npc_results)

# Run Population Density correlations  & linear model
pop_cutoff <- 700 # Use cutoff to exclude the Wuhle (pop. density > 4.8k)

pop_density_df <- feature_table_long %>%
  filter(!is.na(Mode)) %>%
  group_by(Sample, Mode) %>%
  summarise(
    PopDensity = dplyr::first(na.omit(ATTRIBUTE_PopulationDensity)),
    .groups = "drop"
  ) %>%
  left_join(
    xenobiotic_percent_sample %>% select(Sample, Mode, XenobioticPercentSample, RiverGroup),
    by = c("Sample", "Mode")
  ) %>%
  mutate(PopDensity = as.numeric(PopDensity))

pop_df_use <- pop_density_df %>%
  filter(!is.na(PopDensity), !is.na(XenobioticPercentSample)) %>%
  filter(PopDensity <= pop_cutoff)

pop_groups <- list(
  All = pop_df_use,
  Non_PAG = pop_df_use[pop_df_use$RiverGroup %in% c("Spree", "Neckar"), ],
  Spree = pop_df_use[pop_df_use$RiverGroup == "Spree", ],
  Neckar = pop_df_use[pop_df_use$RiverGroup == "Neckar", ],
  Neckar_PAG = pop_df_use[pop_df_use$RiverGroup == "Neckar_PAG", ]
)

pop_results <- list()

for (river_name in names(pop_groups)) {
  data_group <- pop_groups[[river_name]]
  
  for (m in unique(data_group$Mode)) {
    data <- data_group[data_group$Mode == m, ]
    
    if (nrow(data) < 3) {
      pop_results[[paste("POP", river_name, m, sep = "_")]] <- data.frame(
        Analysis = "POP",
        River = river_name,
        NPC_pathway = "PopDensity",
        Mode = m,
        Measure = "PopDensity",
        Shaprio_Xenobiotic_p = NA,
        Shaprio_p = NA,
        Pearson_r = NA, Pearson_p = NA,
        Spearman_rho = NA, Spearman_p = NA,
        Kendall_tau = NA, Kendall_p = NA,
        LM_Slope = NA, LM_Intercept = NA, LM_R2 = NA, LM_p = NA
      )
      next
    }
    
    shapiro_xenobiotic_p <- safe_shapiro(data$XenobioticPercentSample)
    shapiro_pop_p        <- safe_shapiro(data$PopDensity)
    
    pearson_pop  <- safe_cor(data$XenobioticPercentSample, data$PopDensity, "pearson")
    spearman_pop <- safe_cor(data$XenobioticPercentSample, data$PopDensity, "spearman")
    kendall_pop  <- safe_cor(data$XenobioticPercentSample, data$PopDensity, "kendall")
    
    lm_pop <- safe_lm(PopDensity ~ XenobioticPercentSample, data)
    
    pop_results[[paste("POP", river_name, m, "PopDensity", sep = "_")]] <- data.frame(
      Analysis = "POP",
      River = river_name,
      NPC_pathway = "PopDensity",
      Mode = m,
      Measure = "PopDensity",
      Shaprio_Xenobiotic_p = shapiro_xenobiotic_p,
      Shaprio_p = shapiro_pop_p,
      Pearson_r = pearson_pop$estimate, Pearson_p = pearson_pop$p.value,
      Spearman_rho = spearman_pop$estimate, Spearman_p = spearman_pop$p.value,
      Kendall_tau = kendall_pop$estimate, Kendall_p = kendall_pop$p.value,
      LM_Slope = lm_pop$slope, LM_Intercept = lm_pop$intercept, LM_R2 = lm_pop$r2, LM_p = lm_pop$p
    )
  }
}

correlation_df_pop <- do.call(rbind, pop_results)

# Combine the NPC and Pop density correlation results into one csv file
correlation_df_all <- dplyr::bind_rows(correlation_df_npc, correlation_df_pop)

# Write results
# Define output file path & write the csv file
output_path_all <- file.path(Directory, "02_Output_DataProducts", "S11_Samples_Correlations_LinearReg_Anno1Xenobiotics.csv")
write.csv(correlation_df_all, output_path_all, row.names = FALSE)




                                  
                    
                     
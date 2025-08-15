library(readr)
library(dplyr)
library(stats)

# import the txt file and convert to data frame
data <- read.delim("Dros_H3K9ac_bulkChIC_Analysis/push_results/integration_counts.txt", header = FALSE)

col_names <- (c(
  "sample", "total_reads", "chr", 
  "start_pos", "end_pos", "site",
  "overlap_reads", "bases", "coverage", 
  "overlap_fraction"
))

# assign column names
colnames(data) <- col_names

# Define the controls vs the experimental groups
groupings <- c("cntr", "cntr", "exp", "exp")

# Add the groupings to the data frame
df <- data %>%
  mutate(group = rep(groupings, each = nrow(data) / length(groupings)))

# Create new df reshaped 
df_2 <- df %>% 
  select(c("total_reads", "overlap_reads", "group")) %>% 
  group_by(group) %>%
  summarise(
    total_reads = sum(total_reads),
    overlap_reads = sum(overlap_reads)
  ) %>% 
  mutate(
    nnoverlap_reads = total_reads - overlap_reads
  ) %>% 
  select(-total_reads) 

# Covnert to 2x2 matrix
contingency_table <- as.matrix(df_2[, c("overlap_reads", "nnoverlap_reads")])
rownames(contingency_table) <- df_2$group

print(contingency_table)


# Perform fishers exact test
fisher.test(contingency_table)
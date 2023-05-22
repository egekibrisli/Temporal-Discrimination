# Read the data
all_data <- read.csv("C:/Users/PC/Dropbox/eloras_rats/data/deneme/all_data_step6exp_filtered.csv")

# Create four subsets based on the given conditions
NARO_LM <- subset(all_data, cond == 4 & group_width == "NARO" & group_dir == "LM")
NARO_ML <- subset(all_data, cond == 4 & group_width == "NARO" & group_dir == "ML")
WIDE_LM <- subset(all_data, cond == 4 & group_width == "WIDE" & group_dir == "LM")
WIDE_ML <- subset(all_data, cond == 4 & group_width == "WIDE" & group_dir == "ML")

# Count the number of 1s in SE and ME columns for each group
se_counts <- c(sum(NARO_LM$SE == 1, na.rm = TRUE), sum(NARO_ML$SE == 1, na.rm = TRUE),
               sum(WIDE_LM$SE == 1, na.rm = TRUE), sum(WIDE_ML$SE == 1, na.rm = TRUE))
me_counts <- c(sum(NARO_LM$ME == 1, na.rm = TRUE), sum(NARO_ML$ME == 1, na.rm = TRUE),
               sum(WIDE_LM$ME == 1, na.rm = TRUE), sum(WIDE_ML$ME == 1, na.rm = TRUE))

# Combine SE and ME counts
counts <- rbind(se_counts, me_counts)

# Create a bar plot for SE and ME counts
barplot(counts, beside = TRUE, legend.text = c("SE", "ME"),
        names.arg = c("NARO LM", "NARO ML", "WIDE LM", "WIDE ML"),
        xlab = "Group", ylab = "Count", main = "SRNH and MRNH Counts for Normal Durations")
#install.packages("ggthemes")
# Read data from the CSV file
data <- read.csv("/Users/PC/Dropbox/eloras_rats/data/all_data_step6exp_filtered.csv")

# Filter the data for the four groups and SE/ME values smaller than or equal to 15
naro_lm <- subset(data, group_width == "NARO" & group_dir == "LM" & cond == 4)
naro_ml <- subset(data, group_width == "NARO" & group_dir == "ML" & cond == 4)
wide_lm <- subset(data, group_width == "WIDE" & group_dir == "LM" & cond == 4)
wide_ml <- subset(data, group_width == "WIDE" & group_dir == "ML" & cond == 4)

# Calculate the mean and standard deviation of SE and ME for each group
mean_SE <- sapply(list(naro_lm, naro_ml, wide_lm, wide_ml), function(group) mean(group$SE, na.rm = TRUE))
mean_ME <- sapply(list(naro_lm, naro_ml, wide_lm, wide_ml), function(group) mean(group$ME, na.rm = TRUE))
sd_SE <- sapply(list(naro_lm, naro_ml, wide_lm, wide_ml), function(group) sd(group$SE, na.rm = TRUE))
sd_ME <- sapply(list(naro_lm, naro_ml, wide_lm, wide_ml), function(group) sd(group$ME, na.rm = TRUE))

# Perform statistical tests for significance between SE and ME values in each group
p_values <- sapply(list(naro_lm, naro_ml, wide_lm, wide_ml), function(group) {
  if (length(group$SE) == length(group$ME)) {
    t_test <- t.test(group$SE, group$ME, paired = TRUE)
    return(t_test$p.value)
  } else {
    wilcox_test <- wilcox.test(group$SE, group$ME, paired = TRUE)
    return(wilcox_test$p.value)
  }
})

# Create a data frame for mean values, standard deviations, and p-values
df <- data.frame(
  Groups = c("NARO LM", "NARO ML", "WIDE LM", "WIDE ML"),
  SE_mean = mean_SE,
  ME_mean = mean_ME,
  SE_sd = sd_SE,
  ME_sd = sd_ME,
  p_value = p_values
)

# Create a function to add significance asterisks based on p-values
add_significance <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Add significance levels to the data frame
df$significance <- sapply(df$p_value, add_significance)

# Create a side-by-side bar graph for mean SE and ME values by group with error bars and significance levels
library(ggplot2)
library(ggthemes)

df_long <- tidyr::gather(df, variable, value, -Groups, -SE_sd, -ME_sd, -p_value, -significance)

# Define custom colors
colors <- c("#377EB8", "#E41A1C")

# Define custom theme
custom_theme <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Plot the graph
ggplot(df_long, aes(x = Groups, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(ymin = value - ifelse(variable == "SE_mean", SE_sd, ME_sd),
                    ymax = value + ifelse(variable == "SE_mean", SE_sd, ME_sd)),
                position = position_dodge(width = 0.9), width = 0.2) +
  labs(title = "Mean Short and Medium Duration Expectancy Latencies by Groups For Normal Durations", y = "Latencies (s)", fill = "Variable") +
  scale_fill_manual(values = colors) +
  custom_theme +
  geom_text(data = subset(df_long, variable == "ME_mean"),
            aes(label = significance),
            position = position_dodge(width = 0.9),
            vjust = -1.5)
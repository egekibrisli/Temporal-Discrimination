#install.packages("ggsignif")
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggsignif)

# Read the data from csv file
df <- read.csv("/Users/PC/Dropbox/eloras_rats/data/first_half_step7_columned.csv", check.names = FALSE) %>%
  setNames(make.names(names(.)))

# Create a new column for accuracy based on the value of correct_entry column
df <- df %>%
  mutate(accuracy = ifelse(X.correct_entry. == 1, "Correct", "Wrong"))

# Create a new column for response type
df <- df %>%
  mutate(response_type = paste0(X.group_width., " ", X.group_dir.))

# Filter out rows with missing values
df <- df %>%
  filter(!is.na(X.response_lat.), !is.na(X.correct_entry.), !is.na(first_check_lat))

# Calculate mean and standard deviation for each response_type and accuracy combination
summary_data <- df %>%
  group_by(response_type, accuracy) %>%
  summarise(mean_response_lat = mean(X.response_lat.), sd_response_lat = sd(X.response_lat.))

# Perform pairwise t-tests between "Wrong" and "Correct" groups for each subgroup
p_values <- df %>%
  filter(response_type %in% c("NARO LM", "WIDE LM", "NARO ML", "WIDE ML")) %>%
  group_by(response_type) %>%
  summarise(p_value = t.test(X.response_lat. ~ accuracy, data = cur_data())$p.value)

# Plot the response latency for each group
ggplot(df, aes(x = response_type, y = X.response_lat., fill = accuracy)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.9), width = 0.4, size = 0.8) +
  geom_point(aes(y = X.response_lat.), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), shape = 1, size = 0.35, color = "black") +
  labs(x = "Response type", y = "Response Latency (s)", fill = "Accuracy") +
  ggtitle("Response Time Latencies by Subgroup For First Half Lagged Durations") +
  theme_classic() +
  theme(legend.position = "top") +
  geom_signif(
    data = p_values,
    aes(xmin = response_type, xmax = response_type, y_position = 1.5, annotations = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "")))),
    manual = TRUE,
    inherit.aes = FALSE,
    textsize = 4,
    color = "black",
    vjust = -1
  ) +
  coord_cartesian(ylim = c(0, 2))
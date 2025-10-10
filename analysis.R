# BB2 Data Analysis Script
# Author: Research Team
# Date: 2024
# Description: Analysis of experimental data for BB2 study

# Load required libraries
library(tidyverse)
library(ggplot2)

# Read the data files
data1 <- read.csv("data1.csv")
data2 <- read.csv("data2.csv")

# Summary statistics for data1
cat("Summary Statistics for Data1:\n")
print(summary(data1))

# Group analysis by condition
cat("\nMean scores by condition:\n")
condition_summary <- data1 %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_score = mean(score),
    sd_score = sd(score),
    se_score = sd(score) / sqrt(n())
  )
print(condition_summary)

# Visualization: Scores by condition
ggplot(data1, aes(x = condition, y = score, fill = condition)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Score Distribution by Condition",
    x = "Condition",
    y = "Score"
  ) +
  theme(legend.position = "none")

ggsave("scores_by_condition.png", width = 8, height = 6)

# Analysis of response times from data2
cat("\nResponse Time Analysis:\n")
rt_summary <- data2 %>%
  group_by(participant) %>%
  summarise(
    mean_rt = mean(response_time),
    accuracy_rate = mean(accuracy)
  )
print(rt_summary)

# Statistical test: t-test comparing conditions
cat("\nT-test Results:\n")
control_scores <- data1 %>% filter(condition == "control") %>% pull(score)
experimental_scores <- data1 %>% filter(condition == "experimental") %>% pull(score)

t_test_result <- t.test(control_scores, experimental_scores)
print(t_test_result)

# Effect size (Cohen's d)
cohens_d <- (mean(experimental_scores) - mean(control_scores)) / 
  sqrt((sd(experimental_scores)^2 + sd(control_scores)^2) / 2)
cat("\nCohen's d effect size:", round(cohens_d, 3), "\n")

cat("\nAnalysis complete!\n")

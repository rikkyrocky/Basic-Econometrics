library(ggplot2)


data <- read.csv("births_tidy_2019 (1).csv")


data <- subset(data, dbwt != 9999 & cig_0 != 99)


model <- lm(dbwt ~ cig_0, data = data)
summary(model)


ggplot(data, aes(x = cig_0, y = dbwt)) +
  stat_bin2d(bins = 50) +   # heatmap of counts
  scale_fill_gradient(low = "lightblue", high = "darkred", trans = "log") +
  geom_smooth(method = "lm", color = "black", se = FALSE, size = 1.2) +  # regression line
  labs(
    title = "Birthweight vs Cigarettes Smoked Before Pregnancy",
    x = "Cigarettes Smoked Before Pregnancy",
    y = "Birthweight (grams)",
    fill = "Count"
  ) +
  theme_minimal()


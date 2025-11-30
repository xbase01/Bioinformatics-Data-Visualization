# =============================================================================
# Data: Life Satisfaction, Life Expectancy, and Economic Indicators
# =============================================================================

# Set working directory
setwd("C:/Users/judo4/Documents/Data_Visualization")

# Load required libraries
library(ggplot2)
library(countrycode)
library(ggrepel)      # For better text label positioning
library(scales)       # For better axis formatting
library(viridis)      # For colorblind-friendly palettes
library(patchwork)    # For combining plots
library(extrafont)    # For publication fonts (optional)

# Read datasets
life_sat_vs_exp <- read.csv("life-satisfaction-vs-life-expectancy.csv")
gdp_vs_happiness <- read.csv("gdp-vs-happiness.csv")
countries_meta <- read.table("countries_metadata.txt", header = TRUE, sep = "\t")

# Clean column names
colnames(life_sat_vs_exp)[1] <- "Country"
colnames(life_sat_vs_exp)[2] <- "Country_Code"
colnames(life_sat_vs_exp)[4] <- "Life_expectancy_at_birth"
colnames(life_sat_vs_exp)[5] <- "Life_satisfaction"
colnames(life_sat_vs_exp)[6] <- "Population"

# =============================================================================
# DATA PREPARATION
# =============================================================================

# Prepare 2021 data
life_sat_vs_exp2021 <- na.omit(subset(life_sat_vs_exp, Year == 2021))
life_sat_vs_exp2021$Continent <- countrycode(
  sourcevar = life_sat_vs_exp2021$Country, 
  origin = "country.name", 
  destination = "continent"
)
life_sat_vs_exp2021$Continent[life_sat_vs_exp2021$Country == "Kosovo"] <- "Europe"

# Prepare 2017-2021 data
life_sat_vs_exp2017_2021 <- na.omit(subset(life_sat_vs_exp, Year %in% 2017:2021))
life_sat_vs_exp2017_2021$Continent <- countrycode(
  sourcevar = life_sat_vs_exp2017_2021$Country,
  origin = "country.name",
  destination = "continent"
)
life_sat_vs_exp2017_2021$Continent[life_sat_vs_exp2017_2021$Country == "Kosovo"] <- "Europe"

# =============================================================================
# FIGURE 1: Life Satisfaction vs Life Expectancy (2021)
# Publication-ready bubble plot with selective labeling
# =============================================================================

# Identify outliers and interesting cases for labeling
life_sat_vs_exp2021 <- life_sat_vs_exp2021[order(-life_sat_vs_exp2021$Population), ]
life_sat_vs_exp2021$label_flag <- FALSE

# Label top 15 most populous countries and extreme cases
life_sat_vs_exp2021$label_flag[1:15] <- TRUE
life_sat_vs_exp2021$label_flag[
  life_sat_vs_exp2021$Life_satisfaction > 7.5 | 
    life_sat_vs_exp2021$Life_satisfaction < 4 |
    life_sat_vs_exp2021$Life_expectancy_at_birth < 60
] <- TRUE

# Create labeled subset
labels_df <- subset(life_sat_vs_exp2021, label_flag == TRUE)

# Generate Figure 1
fig1 <- ggplot(life_sat_vs_exp2021, 
               aes(x = Life_satisfaction, 
                   y = Life_expectancy_at_birth, 
                   color = Continent, 
                   size = Population)) +
  geom_point(alpha = 0.6, stroke = 0.5) +
  geom_text_repel(
    data = labels_df,
    aes(label = Country),
    size = 2.5,
    max.overlaps = 20,
    segment.size = 0.2,
    segment.color = "grey50",
    box.padding = 0.5,
    point.padding = 0.3,
    force = 2,
    seed = 42
  ) +
  scale_color_viridis_d(option = "turbo", end = 0.9) +
  scale_size_continuous(
    range = c(1, 20),
    breaks = c(1e7, 1e8, 5e8, 1e9),
    labels = c("10M", "100M", "500M", "1B"),
    name = "Population"
  ) +
  scale_x_continuous(limits = c(2.5, 8), breaks = seq(3, 8, 1)) +
  scale_y_continuous(limits = c(50, 90), breaks = seq(50, 90, 10)) +
  labs(
    x = "Life Satisfaction (Cantril Ladder score, 0-10)",
    y = "Life Expectancy at Birth (years)",
    color = "Continent"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    legend.position = "right",
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "lines"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4), order = 1),
    size = guide_legend(order = 2)
  )

# Save Figure 1
ggsave("Figure1_LifeSat_vs_LifeExp_2021.png", 
       plot = fig1, 
       width = 8, 
       height = 6, 
       dpi = 600, 
       bg = "white")

ggsave("Figure1_LifeSat_vs_LifeExp_2021.pdf", 
       plot = fig1, 
       width = 8, 
       height = 6, 
       device = "pdf")

print(fig1)

# =============================================================================
# FIGURE 2: Life Expectancy Distribution by Continent (2017-2021)
# Publication-ready violin plots with statistical annotations
# =============================================================================

# Calculate summary statistics for annotations
continent_stats <- aggregate(
  Life_expectancy_at_birth ~ Continent + Year,
  data = life_sat_vs_exp2017_2021,
  FUN = function(x) c(mean = mean(x), median = median(x), n = length(x))
)

fig2 <- ggplot(life_sat_vs_exp2017_2021, 
               aes(x = as.factor(Year), 
                   y = Life_expectancy_at_birth, 
                   fill = as.factor(Year))) +
  geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, 
               alpha = 0.8, 
               outlier.shape = 21, 
               outlier.size = 1.5,
               outlier.stroke = 0.3,
               show.legend = FALSE) +
  stat_summary(fun = mean, 
               geom = "point", 
               shape = 23, 
               size = 2.5, 
               fill = "white",
               color = "black",
               stroke = 0.5) +
  facet_wrap(~Continent, 
             ncol = 3, 
             scales = "free_y",
             labeller = labeller(Continent = function(x) paste0(x))) +
  scale_fill_viridis_d(option = "mako", 
                       begin = 0.2, 
                       end = 0.8,
                       name = "Year") +
  scale_y_continuous(breaks = seq(50, 90, 10)) +
  labs(
    x = "Year",
    y = "Life Expectancy at Birth (years)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.6),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save Figure 2
ggsave("Figure2_LifeExp_Distribution_2017-2021.png", 
       plot = fig2, 
       width = 10, 
       height = 7, 
       dpi = 600, 
       bg = "white")

ggsave("Figure2_LifeExp_Distribution_2017-2021.pdf", 
       plot = fig2, 
       width = 10, 
       height = 7, 
       device = "pdf")

print(fig2)

# =============================================================================
# ADDITIONAL ANALYSIS: Correlation statistics
# =============================================================================

# Calculate Pearson correlation
cor_test <- cor.test(life_sat_vs_exp2021$Life_satisfaction, 
                     life_sat_vs_exp2021$Life_expectancy_at_birth,
                     method = "pearson")

cat("\n=============================================================================\n")
cat("STATISTICAL SUMMARY FOR MANUSCRIPT\n")
cat("=============================================================================\n\n")
cat(sprintf("Pearson correlation coefficient (r): %.3f\n", cor_test$estimate))
cat(sprintf("95%% CI: [%.3f, %.3f]\n", cor_test$conf.int[1], cor_test$conf.int[2]))
cat(sprintf("p-value: %.2e\n", cor_test$p.value))
cat(sprintf("Sample size (n): %d countries\n\n", nrow(life_sat_vs_exp2021)))

# Summary by continent
cat("Life Expectancy by Continent (2021):\n")
print(aggregate(Life_expectancy_at_birth ~ Continent, 
                data = life_sat_vs_exp2021, 
                FUN = function(x) c(Mean = mean(x), SD = sd(x), N = length(x))))

cat("\n=============================================================================\n")
cat("Figures saved successfully in working directory\n")
cat("=============================================================================\n")
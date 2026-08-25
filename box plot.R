library(ggpubr)

# Find the original signal files
signal_paths <- list.files(
  path = "results",
  pattern = "^kinome-signal_.*-STK\\.csv$",
  full.names = TRUE
)

# Combine files and remove observations repeated across comparisons
all_signal_data <- signal_paths |>
  map_dfr(
    ~ read_csv(.x, show_col_types = FALSE)
  ) |>
  distinct(
    Group,
    SampleName,
    Peptide,
    .keep_all = TRUE
  )

# Calculate one overall mean slope per biological sample
plot_data <- all_signal_data |>
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    .by = c(Group, SampleName)
  ) |>
  mutate(
    Group = factor(
      Group,
      levels = c(
        "ShSr-C",
        "ShHSFi-C",
        "ShSr-NPTX",
        "ShHSFi-NPTX"
      )
    )
  )
comparisons_to_show <- list(
  c("ShHSFi-C",    "ShSr-C"),
  c("ShHSFi-NPTX", "ShHSFi-C"),
  c("ShHSFi-NPTX", "ShSr-C"),
  c("ShHSFi-NPTX", "ShSr-NPTX"),
  c("ShSr-NPTX",   "ShSr-C")
)

group_violin_plot <- ggplot(
  plot_data,
  aes(
    x = Group,
    y = mean_slope,
    fill = Group
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.55,
    color = "grey30"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.07,
    size = 3,
    alpha = 0.9
  ) +
  stat_compare_means(
    comparisons = comparisons_to_show,
    method = "t.test",
    label = "p.format",
    hide.ns = FALSE,
    tip.length = 0.01,
    step.increase = 0.08,
    size = 3
  ) +
  scale_fill_brewer(
    palette = "Set2",
    guide = "none"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.60))
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Overall Phosphorylation Signal by Group",
    x = NULL,
    y = "Mean peptide slope per biological sample",
    caption = paste(
      "Points represent biological samples.",
      "Brackets show unpaired Welch’s t-test p-values."
    )
  )

group_violin_plot

ggsave(
  "results/overall_mean_group_violin_with_brackets.png",
  group_violin_plot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "results/overall_mean_group_violin_with_brackets.svg",
  group_violin_plot,
  width = 9,
  height = 7
)

#Box plot 

group_boxplot <- ggplot(
  plot_data,
  aes(
    x = Group,
    y = mean_slope,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    alpha = 0.65,
    color = "grey30"
  ) +
  geom_jitter(
    width = 0.07,
    size = 3,
    alpha = 0.9
  ) +
  stat_compare_means(
    comparisons = comparisons_to_show,
    method = "t.test",
    label = "p.format",
    hide.ns = FALSE,
    tip.length = 0.01,
    step.increase = 0.08,
    size = 3
  ) +
  scale_fill_brewer(
    palette = "Set2",
    guide = "none"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.60))
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Overall Phosphorylation Signal by Group in STK",
    x = NULL,
    y = "Mean peptide slope per biological sample",
    caption = paste(
      "Points represent biological samples.",
      "Boxes show the median and interquartile range.",
      "Brackets show unpaired Welch’s t-test p-values."
    )
  )

group_boxplot

ggsave(
  "results/overall_mean_group_boxplot_with_brackets STK.png",
  group_boxplot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "results/overall_mean_group_boxplot_with_brackets STK.svg",
  group_boxplot,
  width = 9,
  height = 7
)

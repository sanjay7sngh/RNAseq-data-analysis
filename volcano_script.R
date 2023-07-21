##making volcano plot with the help of tidyverse package
test_res <- test_res %>%
  rownames_to_column(var = "id") %>%
  mutate(threshold = pval < 0.05)
##upper part not functional yet

ggplot(res_fenb) +
  geom_point(aes(x = log2FoldChange, sy = -log10(pval),
                 color = pval<0.05)) +
  xlab("log2 fold change") +
  ylab("-log10 p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

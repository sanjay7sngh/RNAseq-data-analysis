### boxplot woith selection with grepl
cotton_pheLong <- read_csv("cotton_pheLong.csv")
## subset
x <- cotton_pheLong[grepl("_pau", cotton_pheLong$location),]

## general boxplot
ggplot(x, aes(location, value)) + geom_boxplot(fill = "#6AD5E7", notch = TRUE) + facet_wrap(~trait, scales = "free", nrow = 1) + theme_cowplot()

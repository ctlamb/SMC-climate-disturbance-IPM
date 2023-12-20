library(hrbrthemes)
library(ggplot2)

theme.custom <- list(  theme_ipsum() +
                         theme(
                           plot.title = element_text(size = 18),
                           axis.title.y = element_text(size = 15),
                           axis.title.x = element_text(size = 15),
                           axis.text.x = element_text(size = 12),
                           axis.text.y = element_text(size = 12),
                           strip.text.x = element_text(size = 14),
                           strip.text.y = element_text(size = 14),
                           legend.text = element_text(size = 12),
                           legend.title = element_text(size = 14),
                           legend.position = "bottom"
                         )
)
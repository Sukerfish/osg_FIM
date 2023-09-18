### overall project ggplot theme ###

osg_theme <- theme_bw() +   
  theme(text = element_text(family="serif"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "gray93"),
        title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

## save it as RDS
osg_theme %>% saveRDS('osg_theme.rds')
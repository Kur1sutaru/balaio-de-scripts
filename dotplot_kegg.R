library(ggplot2)
mydata<-KEGG_Enrichment_Results_RS.OB
ggplot(mydata,
       aes_string(
         x = 'Count',
         y = "Description",
         color = 'pvalue'
       )
) +
  geom_point(size=5) +
  scale_color_continuous(
    low = "red",
    high = "blue",
    name = 'pvalue',
    guide = guide_colorbar(reverse = TRUE)
  ) +
  ylab(NULL) + theme(axis.text.y = element_text(
    size = 12,
    angle = 0,
    hjust = 1,
    vjust = 0,
    face = "plain"
  ))  

 ggplot( data, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= abs(NES) > 2.5) ) +
  coord_flip() +
  labs( x="Normalized Enrichment Score",
       title="Hallmark & KEGG pathways NES from GSEA") + 
	   theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())  +
		 geom_text(aes(label = pathway), hjust = ifelse(data$NES > 0, 1, 0), color = "black", y=0 ) +
		 theme(panel.grid.major = element_blank(), 
		       panel.grid.minor = element_blank(),
               panel.background = element_blank())

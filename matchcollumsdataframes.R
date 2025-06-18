# Match two columns of different dataframes and save in a new dataframe

setwd("D:/baylor/HumanAnteriorSegment/downstream_analysis/degs/maintables/combined_degs_pathways")

# Upper case all the info
ciliarykegg = toupper(ciliarytop100geneskegg)


# To macth the columns
genelistWangKdegs <- merge(WangKdegs, genelist,
                           by.x = "Feature.Name", by.y = "ï..Gene" )

write.csv(genelistWangKdegs, "genelistWangKdegs.csv")


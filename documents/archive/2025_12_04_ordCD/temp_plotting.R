nodes <- names(bn_structure$nodes)
adj_matrix <- as.matrix(table(
  factor(bn_structure$arcs[,1], levels = nodes), 
  factor(bn_structure$arcs[,2], levels = nodes)
  )
)

net <- network::network(bn_structure$arcs, directed = TRUE)
cols <- RColorBrewer::brewer.pal(9, "Set1")
names(cols) <- unique(c("material", "caregiving", "childcare", "SES", "material", 
                        "material", "transport", "internet", "psychosocial", "psychosocial", 
                        "psychosocial", "legal", "SES", "outcome", "outcome"))

ggnet2(net, color = c("material", "caregiving", "childcare", "SES", "material", 
                      "material", "transport", "internet", "psychosocial", "psychosocial", 
                      "psychosocial", "legal", "SES", "outcome", "outcome"), 
       label = TRUE, 
       palette = cols)


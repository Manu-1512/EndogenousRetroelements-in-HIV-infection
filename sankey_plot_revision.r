########### Sankey diagram
setwd("E:/17July24")
set.seed(0)
library(dplyr)
library(networkD3)

data <- read.delim(file = "Metadata_Elite_comparison_genes_and_TE_loci.tsv", header = TRUE, sep = "\t")
df1 <- data %>% select(genes, cluster)
df2 <- data %>% select(genes.1, cluster.1)
colnames(df2) <- c("genes", "cluster")

levels = 0:3

# Set the cluster column as a factor with the specified levels
df1$cluster <- factor(df1$cluster, levels = levels)
df2$cluster <- factor(df2$cluster, levels = levels)

# Find shared genes between datasets for each cluster pair
shared_genes <- df1 %>%
  inner_join(df2, by = "genes") %>%
  group_by(cluster.x, cluster.y) %>%
  summarise(shared_genes = n()) %>%
  ungroup() %>%
  rename(df1_cluster = cluster.x, df2_cluster = cluster.y)

write.table(shared_genes, "shared_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

links <- shared_genes %>%
  mutate(source = paste0("df1_", df1_cluster),
         target = paste0("df2_", df2_cluster)) %>%
  select(source, target, shared_genes) %>%
  rename(value = shared_genes)

# Create nodes for the Sankey diagram
nodes <- data.frame(name = unique(c(links$source, links$target)))

# Create a mapping of node names to node IDs
nodes$id <- 0:(nrow(nodes) - 1)
nodes_map <- setNames(nodes$id, nodes$name)

# Update links with node IDs
links$source <- nodes_map[links$source]
links$target <- nodes_map[links$target]

# Add a color column to links based on source and target comparison
links$color <- ifelse(links$source == links$target, "red", "steelblue")

# Create the Sankey diagram
sankey <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  units = "Genes",
  fontSize = 12,
  nodeWidth = 30
)

# Save the Sankey diagram as an HTML file
saveNetwork(sankey, "SankeyDiagram_auto_order.html")

webshot::webshot("SankeyDiagram_auto_order.html","SankeyDiagram_auto_order.html.png", vwidth = 1000, vheight = 900, zoom = 15)


# Create the Sankey diagram
sankey2 <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  units = "Genes",
  fontSize = 12,
  nodeWidth = 30,
  iterations = 0 # To change the order
)

# Save the Sankey diagram as an HTML file
saveNetwork(sankey2, "SankeyDiagram_custom_order.html")

webshot::webshot("SankeyDiagram_custom_order.html","SankeyDiagram_custom_order.html.png", vwidth = 1000, vheight = 900, zoom = 15)



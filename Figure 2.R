### This code produces Figure 2, and associated analyses, in: 
### Hooks AP, Burgess SC. Variation in polyandry, reproductive output, 
# and within-brood genetic diversity in a marine snail population 
# across seasons and years. Marine Ecology Progress Series
# Code written by Alex Hooks, hooksap@gmail.com
# with input from Scott Burgess, sburgess@bio.fsu.edu

# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"

# Load required libraries
library('tidyverse')
library('tidygraph') 
library('ggraph')


# Import data
BestCluster_2018 <- read.csv("BestCluster_2018.csv",header=T) 
BestCluster_2020 <- read.csv("BestCluster_2020.csv",header=T) 

# Prepare 2018 data
network_dat_2018 <- BestCluster_2018 %>% 
  filter(Probability > 0) %>% 
  group_by(MotherID) %>% 
  count(FatherID=FatherID) %>% 
  ungroup()

nodes_2018 <- c(network_dat_2018$MotherID, network_dat_2018$FatherID) %>% 
  unique() %>%
  tibble(label = .) %>%
  rowid_to_column("id")
nodes_2018$ParentCat <- ifelse(grepl("parent",nodes_2018$label),"parent",nodes_2018$label)
nodes_2018$ParentColor <- ifelse(nodes_2018$ParentCat=="parent",2,10)

edges_2018 <- network_dat_2018 %>%
  left_join(nodes_2018, by = c("MotherID"="label")) %>%
  rename(from = "id") %>%
  left_join(nodes_2018, by = c("FatherID"="label")) %>%
  rename("to" = "id") %>%
  select(from, to, n)

graph_tidy_2018 <- tbl_graph(nodes = nodes_2018, edges = edges_2018, directed = F)


# Make plot for 2018
quartz(width=11,height=5)
graph_tidy_2018 %>%
  mutate(Centrality = centrality_authority()) %>%
  ggraph(layout = "auto") + 
  geom_node_point(aes(size=2, colour = ParentColor), show.legend = F) +
  geom_edge_link(aes(width = n), alpha = 0.6, show.legend = F) + 
  scale_edge_width(range = c(0.2, 3)) +
  geom_node_text(aes(label = label), max.overlaps=109, repel = TRUE)




# Prepare 2020 data
network_dat_2020 <- BestCluster_2020 %>% 
  filter(Probability > 0) %>% 
  group_by(MotherID) %>% 
  count(FatherID=FatherID) %>% 
  ungroup()

nodes_2020 <- c(network_dat_2020$MotherID, network_dat_2020$FatherID) %>% 
  unique() %>%
  tibble(label = .) %>%
  rowid_to_column("id")
nodes_2020$ParentCat <- ifelse(grepl("parent",nodes_2020$label),"parent",nodes_2020$label)
nodes_2020$ParentColor <- ifelse(nodes_2020$ParentCat=="parent",2,10)

edges_2020 <- network_dat_2020 %>%
  left_join(nodes_2020, by = c("MotherID"="label")) %>%
  rename(from = "id") %>%
  left_join(nodes_2020, by = c("FatherID"="label")) %>%
  rename("to" = "id") %>%
  select(from, to, n)

graph_tidy_2020 <- tbl_graph(nodes = nodes_2020, edges = edges_2020, directed = F)


# Make plot for 2020
quartz(width=11,height=5)
graph_tidy_2020 %>%
  mutate(Centrality = centrality_authority()) %>%
  ggraph(layout = "auto") + 
  geom_node_point(aes(size=2, colour = ParentColor), show.legend = F) +
  geom_edge_link(aes(width = n), alpha = 0.6, show.legend = F) + 
  scale_edge_width(range = c(0.2, 3)) +
  geom_node_text(aes(label = label), max.overlaps=109, repel = TRUE)

library(tidyverse)
# library(msa)

# fasta <- readAAStringSet("data/Proteom/Anaerostipes_caccae.faa", format = "fasta")
# fasta_align <- msa(fasta, verbose = TRUE)

cofactors <- c("CO2", "H(+)", "AMP", "ATP", "ADP", "GDP", "NAD(+)", "NADH", "NADP(+)", "NADPH", "FAD", "FADH2", "UTP", "CTP", "heme b", "CoA", "FMN", "H2O", "NH4(+)", "phosphate", "hydrogen sulfide")

flux_node_info_biomass_final <- tibble()
flux_graph_info_final <- tibble()
metabolites_filtered_final <- tibble()
cofac_filtered_final <- tibble()
reac_filtered_final <- tibble()

for (file in list.files("output/", pattern = "flux_graph_info", recursive = TRUE)) {
    organism <- dirname(file)
    flux_graph_info <- read_delim(file.path("output/", file), delim = ",")
    flux_graph_info$group <- organism
    flux_graph_info_final <- bind_rows(flux_graph_info_final, flux_graph_info)
}

read_node_csv <- function(file) {
    read_delim(file.path("output/", file), delim = ",")
}

for (file in list.files("output/", pattern = "flux_node_info", recursive = TRUE)) {
    if (grepl("csv", file)) {
        organism <- dirname(file)
        flux_node_info <- try(read_node_csv(file))
        if (inherits(flux_node_info, "try-error")) {
            print(str_glue("read csv: {organism} failed!"))
        } else {
            flux_node_info$group <- organism

            flux_node_info_cofac <- flux_node_info %>%
                filter(Name %in% cofactors) %>%
                filter(Reaction == FALSE) %>%
                group_by(group) %>%
                top_n(100, wt = Flux)

            cofac_filtered <- flux_node_info_cofac %>%
                count(Name)
            cofac_filtered_final <- bind_rows(cofac_filtered_final, cofac_filtered)

            flux_node_info_reac <- flux_node_info %>%
                filter(Reaction == TRUE) %>%
                group_by(group) %>%
                top_n(100, wt = Flux)

            reac_filtered <- flux_node_info_reac %>%
                count(Name)
            reac_filtered_final <- bind_rows(reac_filtered_final, reac_filtered)

            flux_node_info_metabolites <- flux_node_info %>%
                filter(!Name %in% cofactors) %>%
                filter(Reaction == FALSE) %>%
                group_by(group) %>%
                top_n(100, wt = Flux)

            metabolites_filtered <- flux_node_info_metabolites %>%
                count(Name)
            metabolites_filtered_final <- bind_rows(metabolites_filtered_final, metabolites_filtered)

            flux_node_info_biomass <- flux_node_info %>%
                filter(Name == "biomass") %>%
                group_by(Protein, group) %>%
                summarise(mean_flux = mean(Flux))
            flux_node_info_biomass_final <- bind_rows(flux_node_info_biomass_final, flux_node_info_biomass)
        }
    }
}
# flux_graph_info_final %>%
#   select(Number_of_nodes_bio, Number_of_edges_bio, Number_of_nodes, Number_of_edges, group) %>%
#   group_by(group) %>%
#   mutate(mean_nodes = mean(Number_of_nodes), mean_edges = mean(Number_of_edges)) %>%
#   select(Number_of_nodes_bio, Number_of_edges_bio, mean_nodes, mean_edges, group) %>%
#   distinct() %>%
#   pivot_longer(cols = c("mean_nodes", "Number_of_nodes_bio"), values_to = "n_nodes", names_to = "node_group") %>%
#   pivot_longer(cols = c("mean_edges", "Number_of_edges_bio"), values_to = "n_edges", names_to = "edge_group" ) %>%
#   ggplot(aes(x))

flux_graph_info_final$medium <- ifelse(
    grepl("adam", flux_graph_info_final$group, ignore.case = TRUE),
    "adam",
    ifelse(
        grepl("cimIV", flux_graph_info_final$group, ignore.case = TRUE),
        "cimIV",
        "none"
    )
)

flux_graph_info_final$group <- str_remove(flux_graph_info_final$group, "_adam")
flux_graph_info_final$group <- str_remove(flux_graph_info_final$group, "_cimIV")

p <- flux_graph_info_final %>%
    select(Number_of_nodes_bio, Number_of_nodes, group, medium) %>%
    group_by(medium, group) %>%
    mutate(mean_nodes = mean(Number_of_nodes)) %>%
    select(Number_of_nodes_bio, mean_nodes, group, medium) %>%
    distinct() %>%
    pivot_longer(cols = c("mean_nodes", "Number_of_nodes_bio"), values_to = "n_nodes", names_to = "node_group") %>%
    ggplot(aes(x = group, y = n_nodes, fill = node_group)) +
    facet_wrap(~medium, ncol = 2) +
    theme_bw() +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_viridis_d(alpha = 0.6, name = "Treatment", labels = c("Flux analysis filtered", "D-glucose-AS-pw")) +
    # scale_fill_discrete() +
    ggtitle("Reconstructed AS biosythesis pathways") +
    ylab("Number nodes in graph") +
    xlab("Organism") +
    theme(
        axis.text.x = element_text(size = 9, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 9, angle = 0),
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 9, hjust = 1)
    )
ggsave("output/flux_graph.pdf", p, dpi = 300)

p <- reac_filtered_final %>%
    ggplot(aes(y = n, x = reorder(Name, -n), fill = group)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d(alpha = 0.6) +
    ggtitle("Top 100 reactions for each group of the first 500 proteins in their corresponding proteome") +
    ylab("Number of times reaction was in top 100") +
    xlab("Reactions") +
    theme(
        axis.text.x = element_text(size = 4, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 9, angle = 0),
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 9, hjust = 1)
    )
ggsave("output/reac_filtered_final.pdf", p, dpi = 300)

p <- cofac_filtered_final %>%
    ggplot(aes(y = n, x = reorder(Name, -n), fill = group)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d(alpha = 0.6) +
    ggtitle("Top 100 cofactors for each group of the first 500 proteins in their corresponding proteome") +
    ylab("Number of times cofactors was in top 100") +
    xlab("Cofactors") +
    theme(
        axis.text.x = element_text(size = 9, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 9, angle = 0),
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 9, hjust = 1)
    )
ggsave("output/cofac_filtered_final.pdf", p, dpi = 300)

p <- metabolites_filtered_final %>%
    ggplot(aes(y = n, x = reorder(Name, -n), fill = group)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d(alpha = 0.6) +
    ggtitle("Top 100 metabolites for each group of the first 500 proteins in their corresponding proteome") +
    ylab("Number of times metabolite was in top 100") +
    xlab("Metabolte without cofactors") +
    theme(
        axis.text.x = element_text(size = 9, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 9, angle = 0),
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 9, hjust = 1)
    )
ggsave("output/metabolites_filtered_final.pdf", p, dpi = 300)

p <- flux_node_info_biomass_final %>%
    ggplot(aes(y = mean_flux, x = group, fill = group)) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_viridis_d(alpha = 0.6) +
    theme(
        legend.position = "none",
        plot.title = element_text(size = 11)
    ) +
    ggtitle("Mean computed sythesis flux of first 500 preoteins from corresponding proteom") +
    ylab("Mean flux of protein sythesis") +
    xlab("Organisms by medium") +
    theme(
        axis.text.x = element_text(size = 9, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 9, angle = 0),
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 9, hjust = 1)
    )
ggsave("output/flux_node_info_biomass_final.pdf", p, dpi = 300)

#
# flux_graph_info_cimIV %>%
#   select(Number_of_nodes) %>%
#   ggplot(aes(x = Number_of_nodes)) +
#     geom_histogram(binwidth = 2, fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
#     ggtitle("Bin size = 3") +
#     theme_bw() +
#     theme(plot.title = element_text(size = 15))
#
# df_mean_metabolites_min <- flux_node_info %>%
#   filter(Protein %in% min_prots) %>%
#   filter(Reaction == FALSE) %>%
#   group_by(Protein, Name) %>%
#   summarise(mean_flux = mean(Flux)) %>%
#   arrange(desc(mean_flux)) %>%
#   slice(1:50)
#
# df_mean_reactions_min <- flux_node_info %>%
#   filter(Protein %in% min_prots) %>%
#   filter(Reaction == TRUE) %>%
#   group_by(Protein, Name) %>%
#   summarise(mean_flux = mean(Flux)) %>%
#   arrange(desc(mean_flux)) %>%
#   slice(1:50)
#
# df_mean_metabolites_max <- flux_node_info %>%
#   filter(Protein %in% max_prots) %>%
#   filter(Reaction == FALSE) %>%
#   group_by(Protein, Name) %>%
#   summarise(mean_flux = mean(Flux)) %>%
#   arrange(desc(mean_flux)) %>%
#   slice(1:50)
#
# df_mean_reactions_max <- flux_node_info %>%
#   filter(Protein %in% max_prots) %>%
#   filter(Reaction == TRUE) %>%
#   group_by(Protein, Name) %>%
#   summarise(mean_flux = mean(Flux)) %>%
#   arrange(desc(mean_flux)) %>%
#   slice(1:50)
#
# df_mean <- df_mean_metabolites_min
#
# # perform hierarchical clustering on the mean flux values
# d <- dist(df_mean$mean_flux)
# hc <- hclust(d)
#
# # reorder rows and columns of the data frame based on the clustering
# df_mean <- df_mean %>%
#   mutate(row_order = hc$order[match(Name, names(hc$order))],
#          col_order = hc$order[match(Protein, names(hc$order))]) %>%
#   arrange(row_order, col_order) %>%
#   select(-row_order, -col_order)
#
# # create heatmap plot
# ggplot(df_mean, aes(x = Protein, y = Name, fill = mean_flux)) +
#   labs(x = "Protein", y = "Name", fill = "Mean Flux") +
#   geom_tile(colour= "black", linewidth=0.0) +
#   scale_fill_continuous(na.value = "black") +
#   theme_bw()+
#   theme(
#     panel.background = element_rect(fill = "black"),
#     panel.grid = element_line(color = "black"))

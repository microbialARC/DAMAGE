library(ggplot2)
library(dplyr)

# Import from snakemake
cluster_alleles_summary <- readRDS(snakemake@input[["cluster_alleles_summary"]])
# Function to generate plot
get_alleles_plot <- function(alleles_plot_df){
  
  plot <-  ggplot(alleles_plot_df,
                  aes(x=cluster,
                      y=num_alleles)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = function(x) c(0, max(x) * 1.05)) +
    geom_bar(fill = "#86C2E8",
             stat = "identity") + 
    geom_text(aes(label = num_alleles),
              color = "#305583",
              vjust = -0.5,
              size = 5,
              fontface = "bold") + 
    labs(x = "Cluster",
         y = "Alleles") +
    theme(axis.text.x = element_text(face = "bold",
                                     angle = 0,
                                     size = 15),
          axis.text.y = element_text(face = "bold",
                                     size = 25),
          axis.title.x = element_text(face = "bold", size = 30),
          axis.title.y = element_text(face = "bold", size = 30),
          plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          legend.key.size = unit(1,"cm"),
          axis.line.x = element_line(linewidth = 1.5),
          axis.line.y = element_line(linewidth = 1.5),
          axis.ticks.length=unit(0.3,"cm"))
  return(plot)
}

# Ortholog alleles 
ortholog_alleles_plot_df <- data.frame(cluster = sapply(cluster_alleles_summary[["ortholog"]], function(cluster) cluster$Cluster),
                                       num_alleles = sapply(cluster_alleles_summary[["ortholog"]], function(cluster) length(cluster$Allele_ID)))

ortholog_alleles_plot_df$cluster <- factor(ortholog_alleles_plot_df$cluster,
                                           levels = ortholog_alleles_plot_df$cluster[order(as.integer(ortholog_alleles_plot_df$cluster))])

ortholog_alleles_plot <- get_alleles_plot(ortholog_alleles_plot_df)

pdf(file=snakemake@output[["ortholog_alleles_visual"]],
    width=22.5,
    height=8)
print(ortholog_alleles_plot)
dev.off()


# Exclusive alleles
#exclusive_alleles_plot_df <- data.frame(cluster = sapply(cluster_alleles_summary[["exclusive"]], function(cluster) cluster$Cluster),
#                                       num_alleles = sapply(cluster_alleles_summary[["exclusive"]], function(cluster) length(na.omit(cluster$Allele_ID))))

#exclusive_alleles_plot_df$cluster <- factor(exclusive_alleles_plot_df$cluster,
#                                           levels = exclusive_alleles_plot_df$cluster[order(as.integer(exclusive_alleles_plot_df$cluster))])

#exclusive_alleles_plot <- get_alleles_plot(exclusive_alleles_plot_df)

#pdf(file=snakemake@output[["exclusive_alleles_visual"]],
#    width=22.5,
#    height=8)
#print(exclusive_alleles_plot)
#dev.off()


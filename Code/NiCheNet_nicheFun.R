
nicheFun = function(background_expressed_genes,
                    expressed_genes_sender,
                    expressed_genes_receiver,
                    geneset_oi,
                    ligand_target_matrix,
                    lr_network, 
                    weighted_networks,
                    weighted_networks_lr,
                    dotplot_object
){
  ## potatial ligand
  ligands = unique(lr_network$from)
  receptors = unique(lr_network$to)
  
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  ## ligand activity
  ligand_activities = predict_ligand_activities(
    geneset = geneset_oi, 
    background_expressed_genes = background_expressed_genes, 
    ligand_target_matrix = ligand_target_matrix, 
    potential_ligands = potential_ligands)
  
  ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
  
  ## ligand-target heatmap
  best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi, 
    ligand_target_matrix = ligand_target_matrix, 
    n = 200) %>% bind_rows() %>% drop_na()
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                   ligand_target_matrix = ligand_target_matrix, 
                                                                   cutoff = 0.33)
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", 
                                                                      color = "purple",legend_position = "top", x_axis_position = "top",
                                                                      legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + 
    scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  
  ## top ligands
  ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  
  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", 
                                                                legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + 
    theme(legend.text = element_text(size = 9))
  
  ## ligand-receptor heatmap
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  vis_ligand_receptor_network = as.matrix(data.frame(vis_ligand_receptor_network)[,intersect(rownames(vis_ligand_pearson),colnames(vis_ligand_receptor_network))])

  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

  ## dotplot of ligands
  rotated_dotplot = DotPlot(dotplot_object, 
                          features = gsub('\\.','-',rownames(vis_ligand_pearson)))+RotatedAxis()+coord_flip()

  ## combined plot
  row1 = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", 
                          axis.ticks = element_blank(), 
                          axis.title.x = element_text(size = 12), 
                          axis.text.y = element_text(face = "italic", size = 9),
                          axis.text.x = element_text(size = 9,  angle = 45,hjust = 0)) + 
    ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_receptor_network+ theme(legend.position = "none"),
  #align = "hv",
  ncol = 3,
  rel_widths = c(ncol(vis_ligand_pearson)+8, 12,  nrow(vis_ligand_receptor_network)-5))

  legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  cowplot::plot_grid(ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
                     ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_receptor_network)),
                     ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
                     ncol = 1),
  nrow=1,
  align = "hv")

  row2 = cowplot::plot_grid(
  p_ligand_target_network + theme(legend.position = "none") + ylab(""),
  legends,
  nrow=1,
  rel_widths = c(ncol(vis_ligand_target),12))
  
  combined_plot = cowplot::plot_grid(row1, row2,rel_heights = c(14,12), nrow = 2, align = "hv")

  return(list(row1 = row1,row2 = row2, combined = combined_plot))
  
}




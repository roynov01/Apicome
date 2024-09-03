matnorm = function(df) {
  # normilizes a dataframe, each value is devidded by the sum of its column
  ind = unlist(lapply(df, is.numeric), use.names = FALSE)  
  numeric_df = df[,ind]
  sum_col = unname(colSums(numeric_df,na.rm = T))
  normilized_df =  as.data.frame(t(t(numeric_df) / sum_col))
  final_df = data.frame(df[,!ind],normilized_df)
  colnames(final_df) = c(colnames(df)[!ind],colnames(df)[ind])
  return (final_df)
}

merge_quick = function(df_list,by=NULL,all=F,suffixes=c("_x","_y")){
  # same as "merge", but for a list of dataframes. if by=NULL, will use rownames
  if (!is.list(df_list) || length(df_list) < 2) {stop("Please provide a list containing at least two data frames.")}
  if (is.null(by)){
    by="temp42"
    df_list = lapply(df_list, function(df) {
      df$temp42 = rownames(df)
      return(df)})
  }
  merged_df = Reduce(function(x,y) {merge(x,y,by=by,all=all,suffixes=suffixes)},df_list)
  if (by=="temp42"){
    rownames(merged_df) = merged_df$temp42
    merged_df$temp42 = NULL
  }
  return(merged_df)
}

remove_duplicated_rows = function(df, column_containing_duplicates) {
  library(data.table)
  dt = data.table(df)
  dt = unique(dt, by = column_containing_duplicates)
  return (as.data.frame(dt))
}

averege_duplicated_rows = function(df, column_containing_duplicates) {
  return (aggregate( . ~ eval(parse(text=column_containing_duplicates)), df, mean, na.rm=T, na.action=NULL))
}

max_duplicated_rows = function(df, column_name) {
  df$medians = apply(df[sapply(df, is.numeric)],1,median)
  df = df[order(df$medians, decreasing = TRUE), ]
  res = df[!duplicated(df[,column_name]), ]
  res$medians = NULL
  return (res)
}

pca = function(df,k_means=NULL,first_pc=1,title="PCA",number_of_genes=20) {
  #' creates a PCA of a dataframe. rownames are genes, colnames are samples.
  #' k_means = how many clusters to make(1<k_means<ncol(df). Will compute silhouette score. if NULL - will not cluster.
  #' first_pc = will show the PC, and +1 from it
  #' number_of_genes = for variable genes plots
  #' returns: plots and data (PCA, elbow), and silhouette.
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  
  if (first_pc >= ncol(df)){stop(paste0("no PC ",first_pc+1," possible in a data of ",ncol(df)," samples"))}
  
  genes_sum_expression = rowSums(df)
  df_filtered = df[genes_sum_expression>0,]
  pca_result = prcomp(t(df_filtered),center=T,scale.=T)
  var_genes=data.frame(pca_result$rotation)
  pca_data = data.frame(pca_result$x)
  
  var_explained = pca_result$sdev^2
  var_explained = var_explained / sum(var_explained)
  elbow_data = data.frame(PC = 1:length(var_explained), Variance = var_explained)
  
  if (!is.null(k_means)){
    if (k_means >= nrow(pca_data)){stop("k_means needs to be lower than number of samples")}
    library(cluster)
    kmeans_result = kmeans(pca_result$x[, first_pc:first_pc+1], centers = k_means)
    dissimilarity_matrix = dist(pca_data)
    pca_data$cluster = as.factor(kmeans_result$cluster)
    silhouette_scores = silhouette(kmeans_result$cluster, dissimilarity_matrix)
    silhouette_score = summary(silhouette_scores)$avg.width
  } else {
    pca_data$cluster = "1"
  }
  pca_data$sample = rownames(pca_data)
  x = paste0("PC",first_pc)
  y = paste0("PC",first_pc+1)
  xlab = paste(paste0("PC ",first_pc," (",signif(var_explained[first_pc]*100,3),"%)"))
  ylab = paste(paste0("PC ",first_pc+1," (",signif(var_explained[first_pc+1]*100,3),"%)"))
  pca_plot = ggplot(pca_data, aes(x=!!sym(x),y=!!sym(y))) +
    geom_point(aes(color=cluster), size=3) +
    geom_text_repel(aes(label=sample)) +
    theme_bw() +
    theme(panel.grid=element_blank()) + 
    scale_color_brewer(palette = "Set1") +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title)
  elbow_plot = ggplot(elbow_data, aes(x=PC, y=Variance)) +
    geom_point() +
    theme_bw() +
    geom_line() +
    xlab("Principal Component") +
    ylab("Proportion of Variance Explained") +
    ggtitle("Elbow Plot")
  
  var_genes$gene = rownames(var_genes)
  df1 = var_genes[order(var_genes[,x]),c("gene",x)]
  df1 = rbind(head(df1, number_of_genes), tail(df1, number_of_genes))
  df2 = var_genes[order(var_genes[,y]),c("gene",y)]
  df2 = rbind(head(df2, number_of_genes), tail(df2, number_of_genes))
  
  genes1 = ggplot(df1, aes(x = reorder(gene,!!sym(x)),!!sym(x))) +
    geom_bar(stat = "identity",fill="blue") +
    geom_hline(yintercept=0, color="black") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid=element_blank()) +
    xlab("")
  genes2 = ggplot(df2, aes(x = reorder(gene,!!sym(y)),!!sym(y))) +
    geom_bar(stat = "identity",fill="blue") +
    geom_hline(yintercept=0, color="black") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid=element_blank()) +
    xlab("")
  plot_genes = plot_grid(genes1,genes2,ncol=1)
  
  if (!is.null(k_means)){
    lab = paste0("Silhouette mean score: ",as.character(round(silhouette_score,2)))
    pca_plot = pca_plot + labs(subtitle=lab)
    return(list(plot_pca=pca_plot, plot_elbow=elbow_plot, pca_df=pca_data, 
                variance=elbow_data, silhouette=silhouette_scores, 
                var_genes=var_genes[,!grepl("gene",colnames(var_genes))],
                plot_genes=plot_genes))
  }
  return(list(plot_pca=pca_plot, plot_elbow=elbow_plot, pca_df=pca_data, variance=elbow_data, 
              var_genes=var_genes[,!grepl("gene",colnames(var_genes))],plot_genes=plot_genes))
}

find_gene = function(data,part_name,gene=F) {
  # finds genes based on regex of the gene name or the full protein name
  # requires metadata_features to be loaded
  part_name = toupper(part_name)
  if (gene==F){
    names = metadata_features$Protein_names[grepl(part_name, toupper(metadata_features$Protein_names))]
    genes = metadata_features$gene[grepl(part_name, toupper(metadata_features$Protein_names))]
  } else {
    names = metadata_features$Protein_names[grepl(part_name, toupper(metadata_features$gene))]
    genes = metadata_features$gene[grepl(part_name, toupper(metadata_features$gene))]
  }
  return (intersect(rownames(data),genes))
}

kegg_get_genes = function (pathway="hsa04975"){
  # get gene names of a pathway
  library(KEGGREST)
  genes = keggGet(pathway)[[1]]$GENE
  genes = data.frame(matrix(genes,ncol=2, byrow=TRUE))
  colnames(genes) = c("id","name")
  genes = separate(genes, name, c('gene', 'name'),sep=";")
  return (genes)
}

kegg_get_pathways = function(partial_name=NA,organism="human",type="pathway"){
  # get pathway ids of the organism
  library(KEGGREST)
  # to get other types: listDatabases()
  if (organism=="human"){organism="hsa"}
  if (organism=="mouse"){organism="mmu"}
  pathways  = keggList("pathway", organism)
  pathways = data.frame(name=unname(pathways),id=names(pathways))
  if (is.na(partial_name)){
    return(pathways)
  } else {
    partial_name = toupper(partial_name)
    upper_pathways = toupper(pathways$name)
    return(pathways[grepl(partial_name,upper_pathways),])
  }
}

kegg_draw_pathway = function(gene_names,values,pathway.id=NULL,
                             colors=c("red","yellow","green"),output_dir=".",
                             filename="pathview",limit=NULL){
  # draw the pathway, and color by a value.
  # data: data.matrix with rownames = entrez_gene_id, and a column with values
  # pathway id example: "04975"
  library(pathview)
  organism="hsa"
  if (is.null(pathway.id)){stop("\n# # # # # #\nTo find a pathway.id, (for example 04975), use kegg_get_pathways('protein') or other regex\n# # # # # #")}
  input = data.frame(gene=gene_names,value=values)
  pathway_data = merge(input[,c("gene","value")],translation_human[,c("external_gene_name","entrezgene_id")],by.y="external_gene_name",by.x="gene")
  pathway_data = pathway_data[!is.na(pathway_data$entrezgene_id),]
  pathway_data$entrezgene_id = as.character(pathway_data$entrezgene_id)
  rownames(pathway_data) = pathway_data$entrezgene_id
  if (is.null(limit)){
    current_genes = kegg_get_genes(paste0(organism,pathway.id))$gene
    current_genes = pathway_data[pathway_data$gene %in% current_genes,]
    limit=list(gene=c(min(current_genes$value),max(current_genes$value)),cpd=1)
  } else {limit=list(gene=limit,cpd=1)}
  pathway_data = data.matrix(pathway_data[,c("value","gene")])
  pathview(gene.data=pathway_data[,1],pathway.id=pathway.id,species=organism,kegg.native=T,
           low=list(gene=colors[1],cpd="black"),mid=list(gene=colors[2],cpd="black"),high =list(gene=colors[3],cpd="black"),
           out.suffix=filename,kegg.dir=output_dir,limit=limit)
}


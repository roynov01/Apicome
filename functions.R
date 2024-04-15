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


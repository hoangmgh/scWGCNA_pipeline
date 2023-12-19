preprocess_input=function (meta_file, count_file, cluster_label = "cluster", sample_label = "sample", 
    vars_to_keep = c("sample", "immune_infiltration_level", "Condition"), 
    normalize = FALSE, cluster_to_use = c("c5", "c1"), normalize_method = "edgeR", 
    remove_genes = Ig_rm, vst = TRUE,remove_missing=NA
                          ) 
{
    "
    input:
        meta_file: a metadata file for all population of cells, with
                    cluster_label: the columns correspond to cell-type (cluster)
                    sample_label: the column of sampleID
                    vars_to_keep: are the variables to be retained for further analyses
        normalize:whether to normalize the count data
        normalize_method: using existing package for normalization, such as edgeR or deSeq2
        vst: whether to use variance stabilizing transformation alongside with deseq2 normalization
        cluster_to_use: indicating the cell population to be used, default to all existing populations 
        remove_low_count: remove genes with a total sum of counts surpassing this value, or NA if use all available genes
        removes_genes: remove genes from a given blacklist, or NA to keep the existing set of genes
    output:
        a list containing raw_data, normalized data, and metadata for all clusters specified.  
    "
    print("1. read in metadata:")
    ################################################
    dat_meta <<- data.frame(as.matrix(fread(meta_file), rownames = 1))
    clusters = unique(dat_meta[[cluster_label]])
    if (length(cluster_to_use) == 1 & cluster_to_use[1] == "all") {
        print("use all clusters")
        cluster_to_use = clusters
    }
    else {
        if (any(cluster_to_use %in% dat_meta[[cluster_label]])) {
            print("one of cluster is not in datameta")
        }
        
    }
    print(paste0("using a total of:", length(cluster_to_use), 
        " unique clusters"))
    ################################################
    dat_meta_splitted = lapply(clusters, function(i) {
        dat_meta[which(dat_meta[[cluster_label]] == i), ]
    })
    names(dat_meta_splitted) = clusters
    print("2. sanity check  metadata:")
    if (sum(!table(dat_meta[[cluster_label]], dat_meta[[sample_label]]) <= 
        1)) {
        print("more than 2 identical samples  exist in one cluster")
        break
    }
    if ( !sample_label %in% colnames(dat_meta_splitted)) {
        print("sample label must exist in metadata")
    }
    if ( !all(vars_to_keep %in% colnames(dat_meta_splitted))) {
        print("all variables used  must exist in metadata")
    }
    if (!sample_label %in% vars_to_keep ) {
        print("sample label and cluster label must be in vars_to_keep")
        break
    }
    
    dat =as.matrix(fread(count_file), rownames = 1)
   
    if (length(remove_genes)>1 & !is.na(remove_genes[1])){
        print("3. remove gene in the blacklist:")
        dat = dat[which(!rownames(dat) %in% remove_genes), ]       
    }
      print(paste0("remaining genes", nrow(dat)))
    dat = lapply(clusters, function(i) t(dat[, which(dat_meta[[cluster_label]] == 
        i)]))
                 
    names(dat) = clusters
    if (!is.na(remove_missing)){
        print("removing genes with low counts overall:")

        keep = lapply(clusters, function(i) colnames(dat[[i]])[which(colSums(dat[[i]],na.rm=T) >=remove_missing)])
        print(paste0("remaining genes: ", length(keep[[1]])))
    }
    else{
    print(clusters)
    keep = lapply(clusters, function(i) colnames(dat[[i]]))  
    }
 
    names(keep) = clusters
    print("replace NA with 0 for calculating CPM:")
    na_mask = lapply(clusters, function(i) apply(dat[[i]],2,is.na))
    names(na_mask)=clusters
    for (i in clusters){
           dat[[i]][na_mask[[i]]]=0
    }
 
    for (i in seq_along(dat_meta_splitted)) rownames(dat_meta_splitted[[i]]) <- dat_meta_splitted[[i]][[sample_label]]
    for (i in seq_along(dat_meta_splitted)) rownames(dat[[i]]) <- rownames(dat_meta_splitted[[i]])
    if (normalize) {
        if (normalize_method == "deseq2") {
            print("using deseq2 normalization")
            dds = lapply(cluster_to_use, function(i) DESeqDataSetFromMatrix(countData = t(dat[[i]]), 
                colData = data.frame(dat_meta_splitted[[i]]), 
                design = ~1))
            if (vst == TRUE) {
                print("also using VST")
                vsd <- lapply(dds, function(i) vst(i, blind = FALSE))
                mat <- lapply(vsd, function(i) t(assay(i)))
                names(mat) = cluster_to_use

                for (i in clusters){
                    mat[[i]][na_mask[[i]]]=NA
                    dat[[i]][na_mask[[i]]]=NA
                }

                multiExpr = lapply(cluster_to_use, function(i) list(data = mat[[i]][, 
                  keep[[i]]]))
                                   
                names(multiExpr) = cluster_to_use
            }
            else {
                names(dds) = cluster_to_use
                dds = lapply(dds, estimateSizeFactors)
                multiExpr = lapply(cluster_to_use, function(i) list(data = t(log10(counts(dds[[i]], 
                  normalized = T) + 1))))
                for (i in clusters){
                    multiExpr[[i]]$data[na_mask[[i]]]=NA
                    dat[[i]][na_mask[[i]]]=NA
                }

                multiExpr = lapply(cluster_to_use, function(i) {
                    selected=colnames(multiExpr[[i]]$data) %in% keep[[i]]
                return(list(data=multiExpr[[i]]$data[,selected]))
                })
                names(multiExpr) = cluster_to_use
            }
                                   
        }
                                   
        else if (normalize_method == "edgeR") {
            print("using edgeR normalization:")
            multiExpr = lapply(cluster_to_use, function(i) {
            list(data = t(cpm(t(dat[[i]]), log = T)))})
            names(multiExpr) = cluster_to_use
            for (i in clusters){
               multiExpr[[i]]$data[na_mask[[i]]]=NA
               dat[[i]][na_mask[[i]]]=NA
            }

            multiExpr = lapply(cluster_to_use, function(i) {
                selected=colnames(multiExpr[[i]]$data) %in% keep[[i]]
                return(list(data=multiExpr[[i]]$data[,selected]))
            }
            )
            names(multiExpr) = cluster_to_use

        }        
    }
    else {
        multiExpr = lapply(cluster_to_use, function(i) list(data = dat[[i]][, 
         keep[[i]]]))
        names(multiExpr) = cluster_to_use
    }
    return(list(raw_data=dat, normalized_data=multiExpr, metadata=dat_meta_splitted))
}


find_power=function(multiExpr,which_cluster="all",meta,blockSize=30000,verbose=2,corFnc=WGCNA::cor,networkType="signed") {
    "
    parameter:
        multiExpr: a list of normalized datamatrix created from calling preprocess_input()
        which_cluster: which cluster to find power for. default to all
        blockSize: number of genes to be processed in parallel
        corFnc: correlation function to be used, default to WGCNA s cor 
    output:
        powerTables:a list of tables 
    "
    
    if (length(which_cluster)==1 &which_cluster=="all"){
        nSet=length(multiExpr)
        which_cluster=names(multiExpr)
    }
    else {
    nSet=length(which_cluster)
    }
    
    powers <<- c(seq(2,10,by=1), seq(12,20, by=2));
    # Initialize a list to hold the results of scale-free analysis
    powerTables <<- vector(mode = "list", length = nSet);
    # Call the network topology analysis function for each set in turn
    for (set in seq_along(which_cluster)){
      powerTables[[set]] <<-list(data = pickSoftThreshold(multiExpr[[which_cluster[set]]]$data, powerVector=powers,
                                                         verbose = 2,corFnc=corFnc,blockSize=30000,networkType=networkType)[[2]])
    }
    collectGarbage();
    # Plot the results:
    # Will plot these columns of the returned scale free analysis tables
    names(powerTables)=which_cluster
    return(powerTables)
}

construct_network=function(multiExpr,which_cluster="all",nThreads=30,power_list,minModuleSize = 50, reassignThreshold = 0, 
                           mergeCutHeight = .15,numericLabels = TRUE, pamRespectsDendro = TRUE, saveTOMs = FALSE,
                      verbose = 3,maxBlockSize=30000,deepSplit=4,  checkMissingData = TRUE,
                      corFnc=WGCNA::cor,save_output=NA) {
if (length(which_cluster)==1 &which_cluster=="all") cluster_to_use= names(multiExpr)
else cluster_to_use=which_cluster
    
net=lapply(cluster_to_use,function(i) {
     blockwiseModules(step1[[2]][[i]]$data,nThreads=30, power =16,
                      minModuleSize = 50, reassignThreshold = 0, mergeCutHeight = .15, 
                      numericLabels = TRUE, pamRespectsDendro = TRUE, saveTOMs = FALSE,
                      verbose = 3,maxBlockSize=30000,deepSplit=4,  checkMissingData = TRUE,
                      corFnc=cor)})
names(net)=names(cluster_to_use)
if (!is.na(save_output)){
saveRDS(net,save_output)
}
return(net)
}


analyze_network=function (net, multiExpr, meta, use = "pairwise.complete.obs", 
    scale = T, summary = "pc") 
{
    "
    parameter:
        net: the list of network constructed from the previous step (construct_network)
        multiExpr:the list of normalized matrices from step1 (preprocess_input)
        meta: the metadata obtained from step1
        use: this will be feed into the correlation function. for example use='p' means use the none NA entries
        scale: whether to scale the data before correlation calculation. this is advised for using the mean trend strategy (summary ='mean')
        summary: whether to take the average value of the genes within each module ('mean'), or use WGCNA  module eigengene ('pc') 


    output: 
        cor_df: a dataframe storing  correlation coefficients between each module and variables of interest 
        p_df: a dataframe storing  correlation p-value between each module and variables of interest 
        group_gene: a named lists where each element is a list of genes per module
        geneModuleMembership: the membership of each gene to all modules
        mergedColors: the assignment of each gene to their module (color as the module nmae)

    "
    mergedColors = lapply(net, function(i) labels2colors(i$colors))
    names(mergedColors) =  names(net)
    geneModuleMembership = lapply(names(net), function(i) as.data.frame(WGCNA::cor(multiExpr[[i]]$data, 
        net[[i]]$MEs, use = use,method="spearman")))
    names(geneModuleMembership) = names(net)
    group_gene = lapply(names(net), function(i) split(colnames(multiExpr[[i]]$data), 
        mergedColors[[i]]))
    MMPvalue = lapply(names(net), function(i) as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership[[i]]), 
        ncol(net[[i]]))))
    print("test")
    names(group_gene) = names(net)
    names(MMPvalue) = names(net)
    names(mergedColors) = names(net)
    groups = names(net)
    cor_result = lapply(names(net), function(i) {
        nGenes = ncol(multiExpr[[i]]$data)
        nSamples = nrow(multiExpr[[i]]$data)
        if (scale) {
            tmp_data = scale(multiExpr[[i]]$data)
        }
        else {
            tmp_data = (multiExpr[[i]]$data[rownames(meta[[i]]), 
                ])
        }
        print(dim(tmp_data))
        if (summary == "mean") 
            MEs = orderMEs(moduleMeangenes(tmp_data, mergedColors[[i]]))
        else if (summary =="pc") MEs = orderMEs(moduleEigengenes(tmp_data, mergedColors[[i]])$eigengenes)
        moduleTraitCor = WGCNA::cor(MEs, meta[[i]], use = use, 
            method = "pearson")
        moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 
            nSamples)
        rownames(moduleTraitCor) = paste0(i, "_", rownames(moduleTraitCor))
        rownames(moduleTraitPvalue) = paste0(i, "_", rownames(moduleTraitPvalue))
        return(list(Cor = moduleTraitCor, p = moduleTraitPvalue))
    })
    cor_df = data.frame(do.call(plyr::rbind.fill, lapply(cor_result, 
        function(i) data.frame(i$Cor, rownames = rownames(i$Cor)))))
    rownames(cor_df) = unlist(sapply(cor_result, function(i) rownames(i$Cor)))
    p_df = data.frame(do.call(plyr::rbind.fill, lapply(cor_result, 
        function(i) data.frame(i$p, rownames = rownames(i$p)))))
    rownames(p_df) = unlist(sapply(cor_result, function(i) rownames(i$p)))
    return(list(cor_df = cor_df, p_df = p_df,group_gene=group_gene,geneModuleMembership=geneModuleMembership,
               mergedColors=mergedColors))
}





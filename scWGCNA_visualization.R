generate_trend_plots=function (fu_df, multiExpr, group_gene, mergedColors, response_var, 
    condition_label = "Condition", meta, sample_label = "sample", fit_curve="gam",
    xlabel = "% of immune cells found in tissue", ylabel = "scaled log10 CPM", 
    save_dir, w = 15, h = 15,loading_thres=.7,condition_palette=c("#c1c3c5", "#51779c", "#dd7784"),summary="mean") 
{
    "
    This function will generate trend plots with plot_trend() in bulk for a given set of clusters and modules and save this in save_dir
    it will also create a excel sheet for each module, storing gene name, loading, and the trait-important correlation 

    input:
        fu_df: a followup_df from analyze_network() or a dataframe storing clusters and modules name in $cluster and $module
        multiExpr: list of normalized expression for clusters generated from step 1 in preprocess_input(), or you can create a list of the form:
                    multiExpr[[cluster]]$data = normalized_matrix
        group_gene: a grouping of the genes for each cluster, group_gene[[cluster]] will contain lists of modules, 
                        and group_gene[[cluster]][[module]] is a set of genes in the module
                        one can also use the group_gene slot from analyze_network() as well
        mergedColors: a list of color assignment for all genes for each cell cluster. length(mergedColors[[cluster]]) == ncol(multiExpr[[cluster]])
        response_var: the response variable to be plot on the y-axis. this must be in the metadata df
        condition_label: label for the Condition
        meta: metadata df taken from step1, preprocess_input()
        sample_label: label for the sampleId in metadata
        loading_thres: only plot genes with a modulemembership higher than this value
        save_dir: a directory for the plot and excel sheet
        fit_curve: whether to fit a gam, lm, or None to the genes    
    output:
        a list of datafrane containing the loading and gene-trait correlation for each module  given
    "
    library(ggplot2)
    library(gridExtra)
    pp = lapply(1:nrow(fu_df), function(line) {
        cluster = fu_df[line, ]$cluster
        module = fu_df[line, ]$module
        p = plot_trend(step1[[2]], cluster, module, group_gene, 
            mergedColors, meta, response_var, condition_label, 
            sample_label, xlabel = xlabel, ylabel = ylabel, title = paste0(cluster, 
                ":", module), loading_thres = loading_thres, 
            fit_curve = fit_curve,condition_palette = condition_palette,summary=summary)
    })
    ggsave(filename = paste0(save_dir, "/trend_plots.pdf"), plot = marrangeGrob(lapply(pp, 
        "[[", "plot"), nrow = 1, ncol = 1), width = w, height = h)
    test = lapply(pp, function(res) data.frame(list(loading = t(res$loading), 
        importance = t(res$importance))))
    names(test) = paste0(fu_df$cluster, "_", fu_df$module)
    names(pp) = names(test)
    library("openxlsx")
    wb <- createWorkbook()
    for (i in names(test)) {
        test[[i]] = test[[i]][pp[[i]]$genes, ]
        test[[i]]$Gene = rownames(test[[i]])
        addWorksheet(wb, i)
        writeData(wb, i, test[[i]], startRow = 1, startCol = 1)
    }
    saveWorkbook(wb, file = paste0(save_dir, "/module_info.xlsx"), 
        overwrite = TRUE)
    return(test)
}   



plot_trend=function(multiExpr,cluster,module,group_gene,mergedColors,meta,response_var,
                    condition_var,sample_var,loading_thres=.7,
                   condition_palette=c("#c1c3c5","#51779c","#dd7784"),
                   title="c5's MEblue genes plotted in c5",
                   xlabel="",ylabel="",fit_curve="lm",summary="pc"){
    "
    input:
        multiExpr: list of normalized expression for clusters generated from step 1 in preprocess_input(), or you can create a list of the form:
                    multiExpr[[cluster]]$data = normalized_matrix
        cluster: the cell population to be used with group_gene
        module: the genes from the module of the given cluster to plot 
        group_gene: a grouping of the genes for each cluster, group_gene[[cluster]] will contain lists of modules, 
                        and group_gene[[cluster]][[module]] is a set of genes in the module
                        one can also use the group_gene slot from analyze_network() as well
        mergedColors: a list of color assignment for all genes for each cell cluster. length(mergedColors[[cluster]]) == ncol(multiExpr[[cluster]])
        response_var: the response variable to be plot on the y-axis. this must be in the metadata df
        condition_label: label for the Condition
        meta: metadata df taken from step1, preprocess_input()
        sample_label: label for the sampleId in metadata
        loading_thres: only plot genes with a modulemembership higher than this value
        save_dir: saving directory for the plot
        fit_curve: whether to fit a gam, lm, or None to the genes    
        summary: whether to use WGCNA's eigengene or module's mean of scaled genes to summarize the trend 
    output:
        a list of datafrane containing the loading and gene-trait correlation for each module given
    "        

    genes=group_gene[[cluster]][[module]]
    if (summary=="pc")   MEs=  moduleEigengenes(multiExpr[[cluster]]$data, mergedColors[[cluster]])$eigengenes
    else MEs = orderMEs(moduleMeangenes(scale(multiExpr[[cluster]]$data), mergedColors[[i]]))

    loading=cor(MEs[,paste0("ME",module)],multiExpr[[cluster]]$data,method ="spearman",use="p")
    importance=unlist(cor(as.numeric(meta[[cluster]][[response_var]]),multiExpr[[cluster]]$data,method ="spearman",use="p"))
    selected=colnames(loading)[which(abs(loading)>=loading_thres & colnames(loading) %in% genes)]
    importance_selected=importance[which(abs(loading)>=loading_thres& colnames(loading) %in% genes)]
    names(importance_selected)=selected
    module_selected=as.data.frame(scale(multiExpr[[cluster]]$data[,selected]))
    module_selected$mean=apply(module_selected,1,mean)

    
    module_selected$sample=(as.character(meta[[cluster]][[sample_var]]))
    module_selected$Condition=(as.character(meta[[cluster]][[condition_var]]))
    module_selected[[response_var]]=(as.character(meta[[cluster]][[response_var]]))
    module_selected=melt(module_selected)
    module_selected$variable=as.factor(module_selected$variable)
    module_selected[[response_var]]=(as.numeric(module_selected[[response_var]]))
    module_selected$Condition=as.factor(module_selected$Condition)
    module_selected$sample=as.factor(module_selected$sample)

if (fit_curve=="gam") {
gp=ggplot(data=module_selected,
       aes_string(x=response_var, y="value",group="variable",label="sample",fill="Condition")) +
        scale_fill_manual(values=condition_palette)+
   geom_line(size=0.5, alpha=0.3,col="grey")+
        geom_point(size=0.01,col="black",alpha=0.2)+
geom_smooth(data = module_selected[module_selected$variable=="mean",],method='gam',se=F,col="blue")+
geom_point(data = module_selected[module_selected$variable=="mean",],aes_string(x=response_var,y="value"),size=5,col="black",pch=21, size=5)+
        xlab(xlabel)+ylab(ylabel)+
        ggtitle(title) +theme_bw() +
  theme(plot.title = element_text(hjust = 1))+ 
        geom_text(nudge_y=.2,nudge_x=.2,aes(label=ifelse(variable=="mean",as.character(sample),'')))
 }
else if (fit_curve=="lm"){
    gp=ggplot(data=module_selected,
       aes_string(x=response_var, y="value",group="variable",label="sample",fill="Condition")) +
        scale_fill_manual(values=condition_palette)+
   geom_line(size=0.5, alpha=0.3,col="grey")+
        geom_point(size=0.01,col="black",alpha=0.2)+
geom_smooth(data = module_selected[module_selected$variable=="mean",],method='lm',formula=y~x,se=F)+
geom_point(data = module_selected[module_selected$variable=="mean",],aes_string(x=response_var,y="value"),size=5,col="black",pch=21, size=5)+
        xlab(xlabel)+ylab(ylabel)+
        ggtitle(title) +theme_bw() +
  theme(plot.title = element_text(hjust = 1))+ 
        geom_text(nudge_y=.2,nudge_x=.2,aes(label=ifelse(variable=="mean",as.character(sample),'')))
    
}
else{
     gp=ggplot(data=module_selected,
       aes_string(x=response_var, y="value",group="variable",label="sample",fill="Condition")) +
        scale_fill_manual(values=condition_palette)+
   geom_line(size=0.5, alpha=0.3,col="grey")+
        geom_point(size=0.01,col="black",alpha=0.2)+
geom_point(data = module_selected[module_selected$variable=="mean",],aes_string(x=response_var,y="value"),size=5,col="black",pch=21, size=5)+
        xlab(xlabel)+ylab(ylabel)+
        ggtitle(title) +theme_bw() +
  theme(plot.title = element_text(hjust = 1))+ 
        geom_text(nudge_y=.2,nudge_x=.2,aes(label=ifelse(variable=="mean",as.character(sample),'')))
}
    return(list("plot"=gp,"loading"=loading,"importance"=importance,"genes"=genes,"df"=module_selected)) 
}            


plot_each_gene=function(multiExpr,cluster,module,group_gene,mergedColors,meta,response_var,
                    condition_var,sample_var,loading_thres=.7,
                   condition_palette=c("#c1c3c5","#51779c","#dd7784"),
                   title="c5's MEblue genes plotted in c5",
                   xlabel="",ylabel="",save_fig=NA,w=10,h=10,fit_curve="lm"){
    "
    input:
        multiExpr: list of normalized expression for clusters generated from step 1 in preprocess_input(), or you can create a list of the form:
                    multiExpr[[cluster]]$data = normalized_matrix
        cluster: the cell population to be used with group_gene
        module: the genes from the module of the given cluster to plot 
        group_gene: a grouping of the genes for each cluster, group_gene[[cluster]] will contain lists of modules, 
                        and group_gene[[cluster]][[module]] is a set of genes in the module
                        one can also use the group_gene slot from analyze_network() as well
        mergedColors: a list of color assignment for all genes for each cell cluster. length(mergedColors[[cluster]]) == ncol(multiExpr[[cluster]])
        response_var: the response variable to be plot on the y-axis. this must be in the metadata df
        condition_label: label for the Condition
        meta: metadata df taken from step1, preprocess_input()
        sample_label: label for the sampleId in metadata
        loading_thres: only plot genes with a modulemembership higher than this value
        save_fig: saving directory for the plot or returning the ggplot list, if NA
        fit_curve: whether to fit a gam, lm, or None to the genes (default = 'lm')
        summary: whether to use WGCNA's eigengene or module's mean of scaled genes to summarize the trend 
        w: width of the plot (default = 10)
        h: height of the plot (default = 10)
    
    output:
        a list of datafrane containing the loading and gene-trait correlation for each module given
    "       
                    
    genes=group_gene[[cluster]][[module]]
    print(length(genes))   
    MEs=  moduleEigengenes(scale(multiExpr[[cluster]]$data), mergedColors[[cluster]])$eigengenes
    loading=cor(MEs[,paste0("ME",module)],multiExpr[[cluster]]$data,method ="spearman",use="p")
    importance=unlist(cor(as.numeric(meta[[cluster]][[response_var]]),multiExpr[[cluster]]$data,method="spearman",use="p"))
    selected=colnames(loading)[which(abs(loading)>=loading_thres & colnames(loading) %in% genes)]
    importance_selected=importance[which(abs(loading)>=loading_thres & colnames(loading) %in% genes)]
    names(importance_selected)=selected
    module_selected=as.data.frame(scale(multiExpr[[cluster]]$data[,selected]))
    module_selected=module_selected[,order(importance_selected,decreasing=T)]
    importance_selected=importance_selected[order(importance_selected,decreasing=T)]    

    module_selected$mean_trend=apply(module_selected,1,mean)

    module_selected$sample=(as.character(meta[[cluster]][[sample_var]]))
    module_selected$Condition=(as.character(meta[[cluster]][[condition_var]]))
    module_selected[[response_var]]=(as.character(meta[[cluster]][[response_var]]))

    module_selected=melt(module_selected)
    module_selected[[response_var]]=(as.numeric(module_selected[[response_var]]))
    module_selected$variable=as.factor(module_selected$variable)
    if (fit_curve=="lm"){
    gps=lapply(unique(module_selected$variable),function(i) {
        return(ggplot(data=module_selected[module_selected$variable==i,],
        aes_string(x=response_var, y="value",label="sample")) +
            scale_fill_manual(values=condition_palette)+
            geom_point(size=2,col="black")+geom_smooth(method="lm",formula=y~x)+
            xlab("")+ylab("")+ ggtitle(paste0(i," r=",round(importance_selected[i],2))) +theme_bw() +
    theme(plot.title = element_text(hjust = 1)))})}
    else {
    gps=lapply(unique(module_selected$variable),function(i) {
        return(ggplot(data=module_selected[module_selected$variable==i,],
        aes_string(x=response_var, y="value",label="sample_var")) +
            scale_fill_manual(values=condition_palette)+
            geom_point(size=2,col="black")+geom_smooth(method="gam")+
            xlab("")+ylab("")+ ggtitle(paste0(i," r=",round(importance_selected[i],2))) +theme_bw() +
    theme(plot.title = element_text(hjust = 1)))})  
    }
        if (!is.na(save_fig)) {
    ggsave(filename = save_fig,
    plot = marrangeGrob(gps, nrow=3, ncol=3,bottom=xlabel,left=ylabel), 
    width = w, height = h)
        }

        return(gps)
}
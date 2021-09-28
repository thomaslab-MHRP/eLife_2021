doFindMarkers = function(seuratobject, min.pct = .1, logfc.threshold = .25, ident.1, ident.2, group.by = "", abs_logFC = FALSE, test.use = "wilcox"){

    # min.pct, test.use, and logfc.threshold defaults are the same as the Seurat defaults for this function.

    if(group.by != ""){	
        result.findmarkers = FindMarkers(seuratobject, only.pos = FALSE, min.pct = min.pct, logfc.threshold = logfc.threshold, ident.1 = ident.1, ident.2 = ident.2, group.by = group.by, test.use = test.use)
    } else {
        result.findmarkers = FindMarkers(seuratobject, only.pos = FALSE, min.pct = min.pct, logfc.threshold = logfc.threshold, ident.1 = ident.1, ident.2 = ident.2, test.use = test.use)
    }

    if(abs_logFC){
        result.findmarkers$cluster = ifelse(result.findmarkers$avg_logFC > 0, ident.1, ident.2)
        pct1new = ifelse(result.findmarkers$avg_logFC < 0, result.findmarkers$pct.2, result.findmarkers$pct.1)
        pct2new = ifelse(result.findmarkers$avg_logFC < 0, result.findmarkers$pct.1, result.findmarkers$pct.2)
        result.findmarkers$pct.1 = pct1new
        result.findmarkers$pct.2 = pct2new
        result.findmarkers$orig_logFC = result.findmarkers$avg_logFC
        result.findmarkers$avg_logFC = abs(result.findmarkers$avg_logFC)
        result.findmarkers$gene = rownames(result.findmarkers)
    }

    return(result.findmarkers)
}
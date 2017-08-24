## Combined analysis method function

Combined_method <- function(RNAseqcount, label, alpha){

    ## Filtering

    keep <- rowSums(RNAseqcount[, -1])>=10
    nkeep <- sum(keep)
    RNAseqcountnew <- RNAseqcount[keep,]

    genename <- RNAseqcount[, 1]
    
    ## edgeR Exact Test 3.12.0

    design <- cbind(Grp1=1,Grp2vs1 = label)

    edgeR.dgelist <- DGEList(counts = RNAseqcountnew[, -1], group = factor(label))
    edgeR.dgelist <- calcNormFactors(edgeR.dgelist, method = "TMM")

    count.norm <- round(edgeR.dgelist$counts, 0)

    ## DESeq2  1.10.1

    group.con <- data.frame(cbind(c(rep("control", sum(label == 0)), rep("treatment", sum(label == 1))), 
                                  c(rep("single-read", length(label)))))
    colnames(group.con) <- c("condition", "type")
    rownames(group.con) <- c(paste("control", 1:sum(label == 0)), paste("treatment", 1:sum(label == 1)))
    colnames(count.norm) <- rownames(group.con)
    DESeq2.dds <- DESeqDataSetFromMatrix(countData = count.norm, colData = group.con, design = ~condition)
    DESeq2.test <- DESeq(DESeq2.dds, quiet = TRUE)
    DESeq2.pval <- results(DESeq2.test)$pvalue
    DESeq2.adjp <- p.adjust(DESeq2.pval, "BH")
    DESeq2.count<- cbind(RNAseqcountnew, results(DESeq2.test), DESeq2.adjp)

    DESeq2.sig.gene <- as.numeric(rownames(RNAseqcountnew)[which(DESeq2.adjp <= alpha)])

    ## EBSeq  1.10.0

    sizes <- MedianNorm(count.norm)
    EBSeq.out <- EBTest(Data = count.norm, Conditions = factor(label), sizeFactors = sizes, maxround = 5)
    EBSeq.result <- GetDEResults(EBSeq.out, FDR = alpha)
    EBSeq.status <- EBSeq.result$Status
    sum(EBSeq.status == "DE", na.rm=TRUE)
    EBSeq.adjp <- EBSeq.result$PPMat[, 1]

    EBSeq.sig.gene <- as.numeric(rownames(RNAseqcountnew)[which(EBSeq.status == "DE")])

    ## SAMSeq   2.0
    
    label2<-c(rep(1, sum(label == 0)), rep(2, sum(label == 1)))
    SAMSeq.test <- SAMseq(count.norm, label2, resp.type = "Two class unpaired", 
                          geneid = rownames(count.norm), genenames = rownames(count.norm), nperms = 100, 
                          nresamp = 20, fdr.output = alpha)
    SAMSeq.result <- rbind(SAMSeq.test$siggenes.table$genes.up, SAMSeq.test$siggenes.table$genes.lo)
    SAMSeq.s <- strsplit(SAMSeq.result[, 1], "[^[:digit:]]")
    SAMSeq.solution <- as.numeric(unlist(SAMSeq.s))
    SAMSeq.solution <- unique(SAMSeq.solution[!is.na(SAMSeq.solution)])

    SAMSeq.sig.gene <- SAMSeq.solution

    ## NOISeq    2.14.1
    
    NOISeq.data <- readData(count.norm, factors = group.con)
    NOISeq.result <- noiseqbio(NOISeq.data, norm = "tmm", factor = "condition", lc = 1, filter = 0)
    NOISeq.DE <- degenes(NOISeq.result, q = 1 - alpha, M = NULL)
    NOISeq.s <- strsplit(rownames(NOISeq.DE), "[^[:digit:]]")
    NOISeq.solution <- as.numeric(unlist(NOISeq.s))
    NOISeq.solution <- unique(NOISeq.solution[!is.na(NOISeq.solution)])

    NOISeq.sig.gene <- NOISeq.solution

    common.sig <- Reduce(intersect, list(DESeq2.sig.gene, EBSeq.sig.gene, SAMSeq.sig.gene, NOISeq.sig.gene))

    common.select <- DESeq2.count[which(common.sig %in% DESeq2.sig.gene), ]

    mylist <- list(DESeq2.sig = genename[DESeq2.sig.gene], EBSeq.sig = genename[EBSeq.sig.gene], 
                   SAMSeq.sig = genename[SAMSeq.sig.gene], NOISeq.sig = genename[NOISeq.sig.gene], 
                   Combined.sig = genename[common.sig], Combined.sig.table = common.select)

    return(mylist)

}


This is a R package named Combine which use a combined method for RNA-Seq differentiall analysis to achieve close to nominal FDR control when sample size is small such as three in each group.

Package: Combine
Type: Package
Title: Combine: An ensemble method for RNA-Seq differential analysis 
Version: 0.1.0
Author: Dongmei Li
Maintainer: Dongmei Li <Dongmei_Li@urmc.rochester.edu>
Depends: R (>= 3.3.0), edgeR, DESeq2, EBSeq, samr, NOISeq
Description: This ensemble method use four current methods (DESeq2, EBSeq, SAMSeq, NOISeq) with equal weight on their 
             resulting FDR-adjusted P-values to identify significant genes from RNA-Seq data.
License:  GPL (>= 2)
Encoding: UTF-8
LazyData: true

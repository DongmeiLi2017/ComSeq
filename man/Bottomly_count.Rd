\name{Bottomly_count}
\docType{data}
\alias{Bottomly_count}
\title{Evaluating gene expression in C57BL/6J and DBA/2J mouse striatum using RNA-Seq data}
\description{
The RNA-Seq data includes 10 B6 mouse strains and 11 D2 strains with an average of 22 million short sequencing reads generated per sample. 
The reads were aligned to the mouse reference genome, and counts from 16,183 Ensembl genes for both strains were generated after 
preprocessing. The Bottomly_count data in the ComSeq package include data from first 2000 genes in all 21 mouse samples.
}
\usage{Bottomly_count}
\value{
The Bottomly_count RNA-Seq example data set contains the following information:
 \item{gene}{Name of gene}
 \item{SRX}{Name of mouse strain: First 10 are B6 strains and remaining 11 are D2 strains}
}

\format{A matrix containing 2000 rows and 22 columns with each row denoting a gene and each column denoting a sample except the first column denoting the gene names.}
\source{http://bowtie-bio.sourceforge.net/recount/}
\references{ 
Bottomly D, Walter NAR, Hunter JE, Darakjian P, Kawane S, et al. (2011) Evaluating Gene Expression in C57BL/6J and DBA/2J Mouse Striatum Using RNA-Seq and Microarrays. PLOS ONE 6(3): e17820. https://doi.org/10.1371/journal.pone.0017820}
\keyword{datasets}

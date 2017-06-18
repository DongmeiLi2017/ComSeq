\name{Combined_method}
\alias{Ensemble_method}

\title{Combined_method: a R function for RNA-Seq data differential analysis
}
\description{This ensemble method combined four current methods (DESeq2, EBSeq, SAMSeq, NOISeq) with equal weight on their resulting FDR-adjusted P-values to identify significantly differentially expressed genes from RNA-Seq data.
}
\usage{
Combined(RNAseqcount, label, alpha)
}

\arguments{
  \item{RNASeqcount}{
   The input RNA-Seq count matrix with rows and columns denoting features and samples, respectively. The first column of the matrix denotes names of features.
}
  \item{label}{
   A vector for group notation such as 1s denote treatment group and 0s denote control group
}

  \item{alpha}{
   The signifiance level
}
}

\details{The Combined method combined the analysis results from four current popular RNA-Seq differential analysis methods (DESeq2, EBSeq, SAMSeq, NOISeq) and declare a gene is significantly differentially expressed only when the gene is identified by all four methods.}

\value{
Combined_method produces a named list with the following components:
\item{DESeq2.sig}{Feature names identified by DESeq2 differential analysis method}
\item{EBSeq.sig}{Feature names identified by EBSeq differential analysis method}
\item{SAMSeq.sig}{Feature names identified by SAMSeq differential analysis method}
\item{NOISeq.sig}{Feature names identified by NOISeq differential analysis method}
\item{Combined.sig}{Feature names identified by Combined differential analysis method}
\item{Combined.sig.table}{Significant feature tables from combined method including original feature counts from all samples, test statistics from DESeq2 method, and adjusted P-values from DESeq2 method}
}

\references{
Love MI, Huber W and Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, pp. 550. doi: 10.1186/s13059-014-0550-8.
Leng N and Kendziorski C (2015). EBSeq: An R package for gene and isoform differential expression analysis of RNA-seq data. R package version 1.14.0.
Tarazona S, Garcia-Alcalde F, Dopazo J, Ferrer A and Conesa A (2011). “Differential expression in RNA-seq: a matter of depth.” Genome Research, 21(12), pp. 4436.
}
\author{Dongmei Li}

\seealso{
The DESeq2, EBSeq, SAMSeq, and NOISeq packages.
}

\examples{
library(ComSeq)
system.file("data", package = "ComSeq")
data("Bottomly_count")
summary(Bottomly_count)
label <- c(rep(0, 10), rep(1, 11))
Result <- Combined_method(RNAseqcount = Bottomly_count, label = label, alpha = 0.05)
summary(Result)
Result
}

\keyword{ RNA-Seq }
\keyword{ Differential Analysis }

#' Eval_Compactness.R
#'
#' Produces a heatmap of an autocorrelation matrix of genes involved in a gene signature, which can be used to evaluate signature compactness.
#' @param sig_dir This is the signature directory file location
#' @param genes_dir This is the gene expression file location
#' @keywords compactness
#' @export
#' @examples
#' eval_compactness()

#


eval_compactness <- function(sig_dir, genes_dir){
	library(gplots)

	gene_sig <- read.table(sig_dir,stringsAsFactors=F, quote="",sep="\t",header=F)
	gene_sig <- gene_sig[,1]
	gene_sig <- as.character(gene_sig)
	
	#load the expression matrix
	mRNA_expr_matrix <- read.table(genes_dir, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(mRNA_expr_matrix) <- gsub('[.]','-',colnames(mRNA_expr_matrix))

	autocors <- cor(t(mRNA_expr_matrix[gene_sig,]),method='spearman')
	pdf('sig_eval_compactness.pdf')
	heatmap.2( autocors,
		col = redgreen(100),#colorpanel(100,"red","yellow","green"),
		trace = "none", 
		xlab = "Gene ID",
		na.color="blue",
		labRow=F,
		main = paste("Autocorrelation of gene signature"),
		dendrogram = "none",
		symbreaks = FALSE,
		Rowv = T,Colv=T )
	dev.off()
}
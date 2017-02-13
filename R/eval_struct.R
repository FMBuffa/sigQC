#' Eval_Struct.R
#'
#' Evaluates the structure of the signature when considering expression across patients. Enables visual clustering of signature genes/patients.
#' @param sig_dir This is the signature directory file location
#' @param genes_dir This is the gene expression file location
#' @keywords expression
#' @export
#' @examples
#' eval_expr()

#

eval_struct <- function(sig_dir, genes_dir, logT=FALSE){
	library(gplots)
	#load the gene signature
	gene_sig <- read.table(sig_dir,stringsAsFactors=F, quote="",sep="\t",header=F)
	gene_sig <- gene_sig[,1]
	gene_sig <- as.character(gene_sig)
	
	#load the expression matrix
	mRNA_expr_matrix <- read.table(genes_dir, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(mRNA_expr_matrix) <- gsub('[.]','-',colnames(mRNA_expr_matrix))
	sig_scores <- as.matrix(mRNA_expr_matrix[gene_sig,])
	if (!logT){
		sig_scores <- log2(sig_scores)
	}
	sig_scores[!(is.finite(sig_scores))] <- NA
	#if user wants a barchart made then do this:
	
	pdf('sig_eval_struct.pdf')
	heatmap.2(t(sig_scores),
		col = redgreen(100),#colorpanel(100,"red","yellow","green"),
		trace = "none", 
		xlab = "Gene ID",
		na.color="blue",
		labRow=F,
		main = paste("Gene signature expression"),
		dendrogram = "col",
		symbreaks = F,
		Rowv = T,Colv=T )
	dev.off()
	
}
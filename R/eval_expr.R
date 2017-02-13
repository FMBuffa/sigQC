#' Eval_Expr.R
#'
#' Evaluates the expression of the signature genes across samples
#' @param sig_dir This is the signature directory file location
#' @param genes_dir This is the gene expression file location
#' @keywords expression
#' @export
#' @examples
#' eval_expr()

#

eval_expr <- function(sig_dir, genes_dir, make_chart=FALSE){
	#load the gene signature
	gene_sig <- read.table(sig_dir,stringsAsFactors=F, quote="",sep="\t",header=F)
	gene_sig <- gene_sig[,1]
	gene_sig <- as.character(gene_sig)
	
	#load the expression matrix
	mRNA_expr_matrix <- read.table(genes_dir, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(mRNA_expr_matrix) <- gsub('[.]','-',colnames(mRNA_expr_matrix))

	#calculate the porportion of nonzero expression data in the matrix
	genes_expr <- mRNA_expr_matrix[gene_sig,]
	gene_expr_vals <- 1 - (rowSums(genes_expr==0) / length(colnames(genes_expr)))
	gene_expr_vals <- sort(gene_expr_vals)
	
	#if user wants a barchart made then do this:
	if (make_chart){
		pdf('sig_expr_eval.pdf')
		bar_expr <- barplot(gene_expr_vals,xlab="Signature Gene IDs",ylab="Proportion with non-zero expression",main="Signature gene expression across samples", axisnames=F,axis=F)
		text(bar_expr, par("usr")[3], labels = names(gene_expr_vals), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
		axis(2)
		#or we could do a density plot if there are a lot of genes
		plot(density(gene_expr_vals),main="Signature gene expression across samples",ylab="Density")
		dev.off()
	}

	print(paste0("Min expression occurs for gene ID ",names(gene_expr_vals)[1], ", with ", round(gene_expr_vals[1]*100,2),"% non-zero expression."))
}
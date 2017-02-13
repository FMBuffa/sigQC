#' Eval_Var.R
#'
#' Evaluates the variance of the signature genes across samples
#' @param sig_dir This is the signature directory file location
#' @param genes_dir This is the gene expression file location
#' @keywords expression
#' @export
#' @examples
#' eval_var()

#

eval_var <- function(sig_dir, genes_dir, make_chart=FALSE){
	#library(ggplot2)
	#load the gene signature
	gene_sig <- read.table(sig_dir,stringsAsFactors=F, quote="",sep="\t",header=F)
	gene_sig <- gene_sig[,1]
	gene_sig <- as.character(gene_sig)
	
	#load the expression matrix
	mRNA_expr_matrix <- read.table(genes_dir, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(mRNA_expr_matrix) <- gsub('[.]','-',colnames(mRNA_expr_matrix))

	#calculate the coefficient of variation
	coeff_of_var <- apply(mRNA_expr_matrix,1,function(x) sd(as.numeric(x),na.rm=T) / mean(as.numeric(x),na.rm=T))
	coeff_of_var_gene_sig <- coeff_of_var[gene_sig]
	if (make_chart){
		pdf('sig_expr_var.pdf')
		# par(lwd=2)
		# plot(density(na.omit(coeff_of_var)), ylim = c(0, max(max(density(na.omit(coeff_of_var))$y),max(density(na.omit(coeff_of_var_gene_sig))$y))),
		# 	main='Comparison of coefficients of variation',
		# 	xlab='Coefficient of variation')
		# lines(density(na.omit(coeff_of_var_gene_sig)),col='red')
		boxplot(coeff_of_var,coeff_of_var_gene_sig,log="y",
			names=c('All Genes','Gene Signature'),
			ylab='Coefficient of Variation',
			main='Variance of signature genes vs. all genes')
		dev.off()

	}

}
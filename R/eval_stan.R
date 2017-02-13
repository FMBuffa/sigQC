#' Eval_Stan.R
#'
#' Evaluates the effect of standardisation of data on the gene signature score.
#' @param sig_dir This is the signature directory file location
#' @param genes_dir This is the gene expression file location
#' @keywords expression
#' @export
#' @examples
#' eval_stan()



eval_stan <- function(sig_dir, genes_dir, logT = FALSE, make_chart=FALSE){
	library(GGally)

	gene_sig <- read.table(sig_dir,stringsAsFactors=F, quote="",sep="\t",header=F)
	gene_sig <- gene_sig[,1]
	gene_sig <- as.character(gene_sig)
	
	#load the expression matrix
	mRNA_expr_matrix <- read.table(genes_dir, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(mRNA_expr_matrix) <- gsub('[.]','-',colnames(mRNA_expr_matrix))

	if (logT){
		mRNA_expr_matrix <- log2(mRNA_expr_matrix)
	}
	mRNA_expr_matrix[!(is.finite(as.matrix(mRNA_expr_matrix)))] <- NA

	z_transf_mRNA <- mRNA_expr_matrix[gene_sig,]
	for (i in 1:length(gene_sig)){
			z_transf_mRNA[gene_sig[i],] <- (as.numeric(z_transf_mRNA[gene_sig[i],]) - mean(as.numeric(z_transf_mRNA[gene_sig[i],]),na.rm=T)) / sd(as.numeric(z_transf_mRNA[gene_sig[i],]),na.rm=T)
	}
	z_transf_scores <- apply(z_transf_mRNA[gene_sig,],2,function(x) median(na.omit(x)))


	med_scores <- apply(mRNA_expr_matrix[gene_sig,],2,function(x) median(na.omit(x)))
	#mean_scores <- apply(mRNA_expr_matrix[gene_sig,],2,function(x) mean(na.omit(x)))
	#pca1_scores <- prcomp(na.omit(t(mRNA_expr_matrix[gene_sig,])),retx=T)
	#pca1_scores <- pca1_scores$x[,1]

#	common_score_cols <- intersect(names(med_scores),intersect(names(mean_scores),names(pca1_scores)))
	combined_scores <- as.data.frame(cbind(med_scores,z_transf_scores))
	#combined_scores <- as.data.frame(cbind(med_scores[common_score_cols],mean_scores[common_score_cols],pca1_scores[common_score_cols]))
	colnames(combined_scores) <- c('Median',"Z_Median")
	pdf('sig_eval_stan.pdf')
	plots <- ggpairs(combined_scores,upper=list(continuous = wrap('cor', method = "spearman")),title='Standardisation comparison')
	print(plots)
	dev.off()


}
#' make_all_plots.R
#'
#' Makes the plots in each of the subfunctions.
#' @param sig_dir This is the signature directory file location
#' @param genes_dir This is the gene expression file location
#' @keywords expression
#' @export
#' @examples
#' make_all_plots()



make_all_plots <- function(sig_dir, genes_dir, logT = FALSE){

	gene_sig <- read.table(sig_dir,stringsAsFactors=F, quote="",sep="\t",header=F)
	gene_sig <- gene_sig[,1]
	gene_sig <- as.character(gene_sig)
	
	#load the expression matrix
	mRNA_expr_matrix <- read.table(genes_dir, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(mRNA_expr_matrix) <- gsub('[.]','-',colnames(mRNA_expr_matrix))

	eval_var_loc(gene_sig, mRNA_expr_matrix)
	eval_expr_loc(gene_sig,mRNA_expr_matrix)
	compare_metrics_loc(gene_sig,mRNA_expr_matrix,logT)
	eval_stan_loc(gene_sig,mRNA_expr_matrix,logT)
	eval_struct_loc(gene_sig,mRNA_expr_matrix,logT)
}

eval_stan_loc <- function(gene_sig, mRNA_expr_matrix, logT = FALSE){
	library(GGally)

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
	
	combined_scores <- as.data.frame(cbind(med_scores,z_transf_scores))
	colnames(combined_scores) <- c('Median',"Z_Median")
	pdf('sig_eval_stan.pdf')
	plots <- ggpairs(combined_scores,upper=list(continuous = wrap('cor', method = "spearman")),title='Standardisation comparison')
	print(plots)
	dev.off()
}


eval_struct_loc <- function(gene_sig, mRNA_expr_matrix, logT=FALSE){
	library(gplots)
	
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

eval_var_loc <- function(gene_sig, mRNA_expr_matrix){
	
	#calculate the coefficient of variation
	coeff_of_var <- apply(mRNA_expr_matrix,1,function(x) sd(as.numeric(x),na.rm=T) / mean(as.numeric(x),na.rm=T))
	coeff_of_var_gene_sig <- coeff_of_var[gene_sig]
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

eval_expr_loc <- function(gene_sig, mRNA_expr_matrix){

	#calculate the porportion of nonzero expression data in the matrix
	genes_expr <- mRNA_expr_matrix[gene_sig,]
	gene_expr_vals <- 1 - (rowSums(genes_expr==0) / length(colnames(genes_expr)))
	gene_expr_vals <- sort(gene_expr_vals)
	
	#if user wants a barchart made then do this:
	
		pdf('sig_expr_eval.pdf')
		bar_expr <- barplot(gene_expr_vals,xlab="Signature Gene IDs",ylab="Proportion with non-zero expression",main="Signature gene expression across samples", axisnames=F,axis=F)
		text(bar_expr, par("usr")[3], labels = names(gene_expr_vals), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
		axis(2)
		#or we could do a density plot if there are a lot of genes
		plot(density(gene_expr_vals),main="Signature gene expression across samples",ylab="Density")
		dev.off()


	print(paste0("Min expression occurs for gene ID ",names(gene_expr_vals)[1], ", with ", round(gene_expr_vals[1]*100,2),"% non-zero expression."))
}


compare_metrics_loc <- function(gene_sig, mRNA_expr_matrix, logT = FALSE){
	library(GGally)

	if (logT){
		mRNA_expr_matrix <- log2(mRNA_expr_matrix)
	}

	mRNA_expr_matrix[!(is.finite(as.matrix(mRNA_expr_matrix)))] <- NA

	med_scores <- apply(mRNA_expr_matrix[gene_sig,],2,function(x) median(na.omit(x)))
	mean_scores <- apply(mRNA_expr_matrix[gene_sig,],2,function(x) mean(na.omit(x)))
	pca1_scores <- prcomp(na.omit(t(mRNA_expr_matrix[gene_sig,])),retx=T)
	pca1_scores <- pca1_scores$x[,1]

	common_score_cols <- intersect(names(med_scores),intersect(names(mean_scores),names(pca1_scores)))

	combined_scores <- as.data.frame(cbind(med_scores[common_score_cols],mean_scores[common_score_cols],pca1_scores[common_score_cols]))
	colnames(combined_scores) <- c('Median','Mean','PCA1')
	pdf('sig_compare_metrics.pdf')
	plots <- ggpairs(combined_scores,upper=list(continuous = wrap('cor', method = "spearman")),title='Scoring metric comparison')
	print(plots)
	dev.off()
}

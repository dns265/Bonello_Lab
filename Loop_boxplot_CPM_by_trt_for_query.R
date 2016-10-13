#Required packages
#install.packages("ggplot2")
#library(ggplot2)
#library(reshape2)

#source script
  #source("~/David/Ohio State/Research/Fertigator/2013/Transcriptomics/Loop_boxplot_CPM_by_trt_for_query.R")

#this assumes you have just run the R script 'input_for_boxplots_from_annot_query.R' and have just specified output_folder and generated modified .csv file using Excell
#execute script by typing:
  #boxplot_by_gene(output_folder)

#creating graphing fuction
boxplot_by_gene<- function(results, na.rm=TRUE,... ) {
  #importing CPM data arranged with samples in rows and genes in columns
  CPM_table <- read.csv(file.path(output_folder,"query_mw_expression_transpose.csv"))
  #reshaping data to have all RPKM values in one column and adding a column called "Gene" 
  long_data= melt(CPM_table, id=c("TRT", "Sample"), variable.name = "gene", value.name = "CPM")
  #create list of genes in data to loop over
  gene_list<- unique(long_data$gene)
  
  #create for loop to produce ggplot graphs
  for (i in seq_along(gene_list)) {
    
    #create plot for each gene in df
    plot<- 
      ggplot(subset(long_data, long_data$gene==gene_list[i]),
             aes(x=TRT, y=CPM)) + geom_boxplot() +
              labs(x="Treatment") + 
              geom_jitter(aes(colour=TRT), size=5, show.legend=FALSE)+
              geom_text(aes(label=Sample), size=5, show.legend = FALSE)+
            ggtitle(paste(gene_list[i], sep=''))

#save plots as .pdf
ggsave(plot, file=paste(results,"/",gene_list[i], ".pdf", sep=''), scale=2)

#print plots to screen
#print(plot)

  }
}
# https://stackoverflow.com/questions/34384907/how-can-put-multiple-plots-side-by-side-in-shiny-r
# https://excelquick.com/r-shiny/selectinput-dependent-on-another-input/
# https://shiny.rstudio.com/reference/shiny/latest/updateSelectInput.html
# https://shiny.rstudio.com/gallery/selectize-rendering-methods.html
# https://shiny.rstudio.com/articles/layout-guide.html

shinyClust <- function(ClusterInput = NULL, G.df = NULL){

# load libraries
	library(dplyr)
	library(shinyBS)
	library(shiny)
	library(dtwclust)
ui <- fluidPage(
  ## title
  titlePanel("Clustering analysis of clonotype lineages"),
  
    fluidRow(
      column(3,
        wellPanel(
            bsCollapse(id="clusterInfo", open="Clustering Details",
                   
                   bsCollapsePanel("Clustering Details",
                                   numericInput(inputId = 'clustNum',
												label = 'Number of clusters(k)',
												min = 2,value = 5),
                                   radioButtons(inputId = "clustertype",
                                               label = "Type of clustering",
                                               choices = c("Partitional", "Hierarchical"),
                                               selected = 'Hierarchical'),
									uiOutput('clusterdetails')
											
								)                                  
                   ),
				   hr(),
				   # this one depends on the value of clustNum
				   sliderInput(inputId = 'clusterIndex',label = 'specify a cluster to show plots',min=1,max=7,value=1)    
        ) #wellPanel
      ),
      column(4,plotOutput(outputId = 'plotL')),
	  column(5,plotOutput(outputId ='plotUMAP'))
			
	  )
   )

	   
server <- function(input,output,session) {

		

		
		output$clusterdetails <-renderUI({
			switch(input$clustertype,
					'Hierarchical' = selectInput(inputId = 'lineage',
											label = 'hc.lineage',
											choices = c('complete','average','single','ward.D'),
											selected = 'average') ,
					'Partitional' = list(numericInput(inputId = 'iterMax',
												label = 'iterMax',
												min = 5, max=100, value = 20),
										selectInput(inputId = 'centroid',
											    label = 'centroid',
												choices = c('pam','dba','mean','median'),
												selected = 'pam'),
										numericInput(inputId = 'setSeed',
												label = 'seed',
												value = 123456)))		
		})
		
	
		
        ClusterID <- reactive(input$clusterIndex)
	
		df <- ClusterInput@lineages
		k <- reactive(input$clustNum)
		NUM <- seq(1,length(df))
		## get xlim and ylim
		
		UMAP1_min <- floor(min(sapply(df,function(x) min(x$UMAP_1))))-1
		UMAP1_max <- ceiling(max(sapply(df,function(x) max(x$UMAP_2))))+1
		UMAP2_min <- floor(min(sapply(df,function(x) min(x$UMAP_2))))-1
		UMAP2_max <- ceiling(max(sapply(df,function(x) max(x$UMAP_2))))+1
		NUM <-seq(1,length(df))
		
		plotUMAPInput <-function(){
			clusterRST <- reactive({
				switch(input$clustertype,
					'Hierarchical' = tsclust(df,type='h',k=input$clustNum, hierarchical_control(method=input$lineage)),
					'Partitional' = tsclust(df,type='p',k=input$clustNum,seed=input$setSeed,iter.max = input$iterMax,centroid=input$centroid)
				)
			})
		
			## UMAP 
			CLONES <- names(which(clusterRST()@cluster ==input$clusterIndex))
			INDEX <- which(G.df$cdr3 %in% CLONES)
			G.df.highlight <- G.df[INDEX,]
			ggplot(G.df,aes(UMAP1,UMAP2,color=G))+geom_point(alpha=0.3,shape=1,size=0.3)+
			geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=1.2,shape=17)+theme_bw()
			
		}
		
		output$UMAP <-renderPlot({
			print(plotUMAPInput())
		})
		
		plotLinInput <-function(){
			clusterRST <- reactive({
				switch(input$clustertype,
					'Hierarchical' = tsclust(df,type='h',k=input$clustNum, hierarchical_control(method=input$lineage)),
					'Partitional' = tsclust(df,type='p',k=input$clustNum,seed=input$setSeed,iter.max = input$iterMax,centroid=input$centroid)
				)
			})
			
			INDEX <-which(clusterRST()@cluster== input$clusterIndex)
			INDEX2 <-NUM[-INDEX]
			
			## plot lineage
			plot(df[[1]]$UMAP_1,df[[1]]$UMAP_2,type='l',xlab='UMAP1',col='grey',ylab='UMAP2',xlim=c(UMAP1_min,UMAP1_max),ylim=c(UMAP2_min,UMAP2_max),pch=19,lwd=1.1)
			text(df[[1]]$UMAP_1,df[[1]]$UMAP_2,labels=rownames(df[[i]]),cex=0.4,col='grey')

			for (i in INDEX2){
				lines(df[[i]]$UMAP_1,df[[i]]$UMAP_2,type='l',col="grey",pch=19,cex=0.4,lwd=1.1)
				text(df[[i]]$UMAP_1,df[[i]]$UMAP_2,labels=rownames(df[[i]]),cex=0.4,col='grey')
			}

			for (j in INDEX){
				print(j)
				lines(df[[j]][,1],df[[j]][,2],type='l',col="blue",pch=19,cex=0.4,lwd=1.2)
				text(df[[j]][,1],df[[j]][,2],labels=rownames(df[[j]]),cex=0.4)
			}
		}
		
		
		## plot component lineages
		output$plotL <- renderPlot({
			print(plotLinInput())
		})
	}	
	
	shinyApp(ui = ui, server = server)
}

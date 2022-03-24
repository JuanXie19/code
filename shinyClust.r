shinyClust <- function(ClusterInput = NULL, G.df = NULL){

# load libraries
	library(dplyr)
	library(shinyBS)
	library(shiny)
	library(dtwclust)

ui <- fluidPage(

  # App title ----
  titlePanel("Clustering analysis of clonotype lineages"),

  fluidRow(
	 column(3,
		wellPanel(	

		bsCollapse(id="clusterInfo", open="Clustering Details",

                   bsCollapsePanel("Clustering Details",
									numericInput(inputId = 'clustNum',
												label = 'Desired number of clusters',
												min = 2, value = 2),
									numericInput(inputId = 'setSeed',
												label = 'Random seed',value=123456),
								
                                   radioButtons(inputId = "clustertype",
                                               label = "Type of clustering",
                                               choices = c("Partitional", "Hierarchical"),
                                               selected = 'Hierarchical'),
									uiOutput("clusterdetails")
														
								)                                  
                   ),
				   hr(),
				   sliderInput(inputId ='clusterIndex',label = 'select a cluster to show plots',min=1,max=7,value = 1,step=1),
				   hr(),
				   radioButtons(inputId='radio',
									  label = 'Showing option for lineages in a cluster',
									  choices = list('show all lineages'=1,'show representative'=2),									
									  selected = 1
					)
				)),
		column(4,
				h4('plot 1'),
				plotOutput('plotLin'),
				downloadButton('downloadPlot1',label = 'Download plot1')),
		column(5,
				h4('plot 2'),
				plotOutput('plotUMAP'),
				downloadButton('downloadPlot2',label = 'Download plot2'))

	)
)
		   



	server <- function(input,output,session) {
	
		output$clusterdetails <- renderUI({
			if(is.null(input$clustertype))
				return()
				
			switch(input$clustertype,
				'Hierarchical' = selectInput(inputId = 'lineage',
												label = 'hc.lineage',
												choices = c('complet','average','single','ward.D'),
												selected = 'average'),
				'Partitional' = list(numericInput(inputId = 'iterMax',
												label = 'iterMax',
												min = 5, max=100, value = 20),
										selectInput(inputId = 'centroid',
											    label = 'centroid',
												choices = c('pam','dba','mean','median'),
												selected = 'pam')))			
		})

		df <- ClusterInput@lineages

		itermax <- reactive(input$iterMax)
		centroid <- reactive(input$centroid)
		seed <- reactive(input$setSeed)
		k <- reactive(input$clustNum)

		method <- reactive(input$lineage)
		
		
		observeEvent(input$clustNum,{
			updateSliderInput(session,'clusterIndex',max = input$clustNum)
		})
		
		clusterID <- reactive(input$clusterIndex)
		## get xlim and ylim

		UMAP1_min <- floor(min(sapply(df,function(x) min(x$UMAP_1))))-1
		UMAP1_max <- ceiling(max(sapply(df,function(x) max(x$UMAP_2))))+1
		UMAP2_min <- floor(min(sapply(df,function(x) min(x$UMAP_2))))-1
		UMAP2_max <- ceiling(max(sapply(df,function(x) max(x$UMAP_2))))+1
		NUM <-seq(1,length(df))


		UMAP <- function(){
		## clustering 
			clusterRST <- reactive({
				switch(input$clustertype,
					'Hierarchical' = tsclust(df,type='h',k=k(),seed=seed(),control=hierarchical_control(method=method())),
					'Partitional' = tsclust(df,type='p',k=k(),seed=seed(),iter.max=itermax(),centroid=centroid())
			
				)
			})
			
		CLONES <- names(which(clusterRST()@cluster==clusterID()))
		INDEX <- which(G.df$cdr3 %in% CLONES)
		G.df.highlight <- G.df[INDEX,]
		ggplot(G.df,aes(UMAP1,UMAP2,color=G))+geom_point(alpha=0.3,shape=1,size=0.3)+
				geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=1.2,shape=17)+theme_bw()		
		}
		
		output$plotUMAP <- renderPlot({
			print(UMAP())
		
		})
		
		output$downloadPlot1 <- downloadHandler(
			filename = function(){
			'plot1.pdf'},
			content = function(file){
				pdf(file)
				print(UMAP())
				dev.off()
			}		
		)
		
		
		LIN <- function(){
			clusterRST <- reactive({
				switch(input$clustertype,
					'Hierarchical' = tsclust(df,type='h',k=k(),seed=seed(),control=hierarchical_control(method=method())),
					'Partitional' = tsclust(df,type='p',k=k(),seed=seed(),iter.max=itermax(),centroid=centroid())
			
				)
			})
			
			INDEX <- which(clusterRST()@cluster == clusterID())
			INDEX2 <- NUM[-INDEX]
			CENTROID <- as.data.frame(clusterRST()@centroids[[clusterID()]])
			
			plot(df[[1]]$UMAP_1,df[[1]]$UMAP_2,type='l',xlab='UMAP1',col='grey',ylab='UMAP2',xlim=c(UMAP1_min,UMAP1_max),ylim=c(UMAP2_min,UMAP2_max),pch=19,lwd=1.1)
			text(df[[1]]$UMAP_1,df[[1]]$UMAP_2,labels=rownames(df[[1]]),cex=0.4,col='grey')

			for (i in INDEX2){
				lines(df[[i]]$UMAP_1,df[[i]]$UMAP_2,type='l',col="grey",pch=19,cex=0.4,lwd=1.1)
				text(df[[i]]$UMAP_1,df[[i]]$UMAP_2,labels=rownames(df[[i]]),cex=0.4,col='grey')
			}

			for (j in INDEX){
				print(j)
				lines(df[[j]][,1],df[[j]][,2],type='l',col="blue",pch=19,cex=0.4,lwd=1.2)
				text(df[[j]][,1],df[[j]][,2],labels=rownames(df[[j]]),cex=0.4)
			}
			lines(CENTROID[,1],CENTROID[,2],type='l',col='red',pch=19,cex=0.4,lwd=1.2)
			text(CENTROID[,1],CENTROID[,2],labels=rownames(CENTROID),cex=0.4)

		}
		
		
		output$plotLin <- renderPlot({
			print(LIN())
		
		})
		
		output$downloadPlot2 <- downloadHandler(
			filename = function(){
				paste('cluster',input$clusterIndex,'_plot2.pdf',sep='')
			},
			content = function(file){
				pdf(file)
				print(LIN())
				dev.off()
			
			}
		)
		
	}

	shinyApp(ui = ui, server = server)
}
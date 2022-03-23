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
                                   
                                   radioButtons(inputId = "clustertype",
                                               label = "Type of clustering",
                                               choices = c("Partitional", "Hierarchical"),
                                               selected = 'Hierarchical'),
									numericInput(inputId = 'clustNum',
												label = 'Number of clusters(k)',
												min = 2,value = 5),
							
									## only show this if clustering type is partitional
									conditionalPanel(
										condition = "input.clustertype=='Partitional'",
										numericInput(inputId = 'iterMax',
												label = 'iterMax',
												min = 5, max=100, value = 20),
										selectInput(inputId = 'centroid',
											    label = 'centroid',
												choices = c('pam','dba','mean','median'),
												selected = 'pam'),
										numericInput(inputId = 'setSeed',
												label = 'seed',
												value = 123456)
									),							
									
									## only show this if clustering type is hierarchicha
									conditionalPanel(
										condition = "input.clustertype=='Hierarchical'",							
										selectInput(inputId = 'lineage',
												label = 'hc.lineage',
												choices = c('complet','average','single','ward.D'),
												selected = 'average')
									)									
								)                                  
                   ),
				   hr(),
				   # this one depends on the value of clustNum
				   sliderInput(inputId = 'clusterIndex',label = 'specify a cluster to show plots',min=1,max=7,value=1)    
        ) #wellPanel
      ),
      column(9,
		helpText('TEST'),
		#textOutput(outputId = 'desc'),
		fluidRow(
			splitLayout(cellWidths=c('50%','50%')),
			conditionalPanel(
				condition = "input.clustertype=='Hierarchical'",
				plotOutput('plotL_hc')),
				conditionalPanel(
				condition = "input.clustertype=='Partitional'",
				plotOutput('plotL_p')),
				conditionalPanel(
				condition = "input.clustertype=='Hierarchical'",
				plotOutput('plotUMAP_hc')),
				conditionalPanel(
				condition = "input.clustertype=='Partitional'",
				plotOutput('plotUMAP_p'))			
		
		)
			
	  )
   )
  )

	   

server <- function(input,output,session) {
		
		df <- ClusterInput@lineages
		
		itermax <- reactive(input$itermax)
		centroid <- reactive(input$centroid)
		seed <- reactive(input$setSeed)
		k <- reactive(input$clustNum)
	
		method <- reactive(input$lineage)
		
		if(FALSE){observeEvent(input$clustNum,{
			updateSliderInput(session,'clusterIndex',max = input$clustNum)
	
		})}
		
		
		
        ClusterID <- reactive(input$clusterIndex)
	
		## get xlim and ylim
		
		UMAP1_min <- floor(min(sapply(df,function(x) min(x$UMAP_1))))-1
		UMAP1_max <- ceiling(max(sapply(df,function(x) max(x$UMAP_2))))+1
		UMAP2_min <- floor(min(sapply(df,function(x) min(x$UMAP_2))))-1
		UMAP2_max <- ceiling(max(sapply(df,function(x) max(x$UMAP_2))))+1
		NUM <-seq(1,length(df))
		
	
		## plot component cells on UMAP
		output$plotUMAP_p <- renderPlot({
		
		  part.dtw <-tsclust(df,type='p',k=k(),seed=seed(),iter.max=itermax(),centroid=centroid())
		  CLONES <- names(which(part.dtw@cluster == ClusterID()))
		  INDEX <- which(G.df$cdr3 %in% CLONES)
		  G.df.highlight <-G.df[INDEX,]
		  ggplot(G.df,aes(UMAP1,UMAP2,color=G))+geom_point(alpha=0.3,shape=1,size=0.3)+
			geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=1.2,shape=17)+theme_bw()
		},height='50%',width='60%')
	
		output$plotUMAP_hc <-renderPlot({
			hc.dtw <-tsclust(df,type='h',k=k(),control=hierarchical_control(method=method()))
			CLONES <- names(which(hc.dtw@cluster == ClusterID()))
			INDEX <- which(G.df$cdr3 %in% CLONES)
			G.df.highlight <-G.df[INDEX,]
			ggplot(G.df,aes(UMAP1,UMAP2,color=G))+geom_point(alpha=0.3,shape=1,size=0.3)+
				geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=1.2,shape=17)+theme_bw()
		},height='50%',width='60%')
		
		## plot component lineages
		output$plotL_p <- renderPlot({
		
		  part.dtw <-tsclust(df,type='p',k=k(),seed=seed(),iter.max=itermax(),centroid=centroid())
		  
		  INDEX <-which(part.dtw@cluster== ClusterID())
		  INDEX2 <-NUM[-INDEX]
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
		},height='50%',width='50%')
	
		output$plotL_hc <-renderPlot({
			hc.dtw <-tsclust(df,type='h',k=k(),control=hierarchical_control(method=method()))
			INDEX <-which(hc.dtw@cluster== ClusterID())
			INDEX2 <-NUM[-INDEX]
			plot(df[[1]]$UMAP_1,df[[1]]$UMAP_2,type='l',xlab='UMAP1',col='grey',ylab='UMAP2',xlim=c(UMAP1_min,UMAP1_max),ylim=c(UMAP2_min,UMAP2_max),pch=19,lwd=1.1)
			text(df[[1]]$UMAP_1,df[[1]]$UMAP_2,labels=rownames(df[[i]]),cex=0.4,col='grey')

			for (i in INDEX2){
				lines(df[[i]]$UMAP_1,df[[i]]$UMAP_2,type='l',col="grey",pch=19,cex=0.4,lwd=1.1)
				text(df[[i]]$UMAP_1,df[[i]]$UMAP_2,labels=rownames(df[[i]]),cex=0.4,col='grey')
			}

			for (j in INDEX){
				lines(df[[j]][,1],df[[j]][,2],type='l',col="blue",pch=19,cex=0.4,lwd=1.2)
				text(df[[j]][,1],df[[j]][,2],labels=rownames(df[[j]]),cex=0.4)
			}			
			
		},height='50%',width='50%')
	
		
		
		
	}
	
	shinyApp(ui = ui, server = server)
}
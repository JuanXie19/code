shinyClone <- function(CombinedDataSet =NULL){

# load libraries
	library(Seurat)
	library(dplyr)
	library(DT)
# extract needed information
	UMAP1 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,1]
	UMAP2 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,2]
	cdr3 <- CombinedDataSet@seurat$cdr3
	G <- CombinedDataSet@seurat$clusters
	G.df <- data.frame(UMAP1,UMAP2,G,cdr3)
		
	CLONE <- unique(G.df$cdr3)
	n <-length(CLONE)
	
	## find each clonotype involve how many cells
	
	cdr3.df <- data.frame(rbind(table(G.df$cdr3,G.df$G)))
	cdr3.df$cell.num <- rowSums(cdr3.df)
	cdr3.df$cluster.num <-apply(cdr3.df,1, function(x) sum(x>0)-1) # count how many clusters it spans
	
	ui <-fluidPage(
		titlePanel('Clonoal distribution on UMAP'),
		sidebarLayout(
		 sidebarPanel(
		 textInput(inputId = 'cloneName',label = strong('Clonotype'),value = 'CASSFAGGASYEQYF'),
		 downloadButton('downloadPlot',label='Download Plot')
		 ),
		 mainPanel(
			tabsetPanel(
				tabPanel('table',
					DT::dataTableOutput('clonetable')
					),
				tabPanel('plot',
					textOutput(outputId = 'desc'),
					plotOutput(outputId = 'UMAP'))
			)
			
			
			
		 )
		)
	)
	
	server <- function(input,output){
		cloneInd <- reactive({
			req(input$cloneName)
			validate(need (input$cloneName %in% CLONE,'Error:invalid clonotype. Please double check'))
			#validate(need(is.numeric(input$cloneIndex),'Error: cloneIndex should be a integer'))		
			INDEX <-which(G.df$cdr3 == input$cloneName)
			return(INDEX)
		})
		
		
		output$clonetable <- DT::renderDataTable({cdr3.df})
		
		## description of the plotted clonotype
		output$desc <-renderText({
			G.df.highlight <- G.df[cloneInd(),]
			paste('This clonotype consists of ',length(cloneInd()),'cells','spanning',length(unique(G.df.highlight$G)),'cluster/clusters')
		})
		
		## plot the distribution
		plotUMAP <-function(){
			if(length(cloneInd())<3){
				print('too few cells')
			}
			
			G.df.highlight <- G.df[cloneInd(),]
			
			ggplot(G.df,aes(UMAP1,UMAP2,color=G))+geom_point(alpha=0.4,shape=1,size=0.6)+
				geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=2,shape=17)+
				theme_bw()+labs(title=CLONE[cloneInd()])
		
		}
		
		output$UMAP <-renderPlot({
			print(plotUMAP())
		})
		output$downloadPlot <-downloadHandler(
			filename = function(){
			paste(CLONE[cloneInd()],'_UMAP.pdf',sep='')},
			content = function(file){
			 pdf(file)
			 print(plotUMAP())
			 dev.off()
			}
		)
		
	}
	## create shiny object
	shinyApp(ui = ui, server = server)
}
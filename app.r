library(shiny)
library(shinyBS)
library(dtwclust)
library(dplyr)

ui <- fluidPage(

  # App title ----
  titlePanel("Clustering analysis of clonotype lineages"),

  sidebarLayout(
	wellPanel(
		fileInput(inputId = 'upload',NULL,
					label= 'Upload a list of clonoetype lineages',
					buttonLabel = 'Choose File',multiple=FALSE,
					placeholder = 'No file selected',
					accept = '.csv'),
		hr(),
	
		bsCollapse(id="clusterInfo", open="Clustering Details",
                   
                   bsCollapsePanel("Clustering Details",
                                   
                                   radioButtons(inputId = "clustertype",
                                               label = "Type of clustering",
                                               choices = c("Partitional", "Hierarchical"),
                                               selected = 'Hierarchical'),
									sliderInput(inputId = 'clustNum',
												label = 'Number of clusters',
												min = 1,max = 8,value = 1),
							
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
                   )
				),
	mainPanel(
		tabsetPanel(
			tabPanel('Plot',
				conditionalPanel(
				condition = "input.clustertype=='Hierarchical'",
				plotOutput('plot_hc')),
				conditionalPanel(
				condition = "input.clustertype=='Partitional'",
				plotOutput('plot_p'))
			),
			tabPanel('Summary',
			conditionalPanel(
			condition = "input.clustertype=='Hierarchical'",
			verbatimTextOutput('summaryHC')),
			conditionalPanel(
			condition="input.clustertype=='Partitional'",
			verbatimTextOutput('summaryP'))
			),
			tabPanel('Cluster label',
			conditionalPanel(
			condition = "input.clustertype=='Hierarchical'",
			verbatimTextOutput('labelHC')),
			conditionalPanel(
			condition="input.clustertype=='Partitional'",
			verbatimTextOutput('labelP'))
			)
		
		
		)
		
    )		   
  )
 )


server <- function(input,output) {
	my_dtw <-reactive({
		file <-input$upload
		
		ext <- tools::file_ext(file$datapath)
		req(file)
    
		validate(need(ext == "csv", "Please upload a csv file"))
		df <-read.csv(file$datapath,header=TRUE)
		# split the dataframe into list of dataframe
		temp <- df%>%group_split(names)
		NAMES <- rep(NA,length(temp))
		for (i in 1:length(temp)){
		 NAMES[i] <-unique(temp[[i]]$names)
		}
		names(temp) <-NAMES
		temp <- lapply(temp,function(x) x[,c('UMAP_1','UMAP_2')])
		return(temp)
		})
	
	
		
	

	
	itermax <-reactive(input$itermax)
	centroid <-reactive(input$centroid)
	seed <-reactive(input$setSeed)
	k <-reactive(input$clustNum)
	
	
	method <-reactive(input$lineage)

	
	output$plot_p <- renderPlot({
		
		  part.dtw <-tsclust(my_dtw(),type='p',k=k(),seed=12345,iter.max=itermax(),centroid=centroid())
		  plot(part.dtw)
	})
	
	output$plot_hc <-renderPlot({
		hc.dtw <-tsclust(my_dtw(),type='h',k=k(),control=hierarchical_control(method=method()))
		plot(hc.dtw,type='series')
	})
	
	output$summaryP <-renderPrint({
	
		part.dtw <-tsclust(my_dtw(),type='p',k=k(),seed=12345,iter.max=itermax(),centroid=centroid())
		part.dtw
	})
	
	output$summaryHC <-renderPrint({
		hc.dtw <-tsclust(my_dtw(),type='h',k=k(),control=hierarchical_control(method=method()))
		hc.dtw
	
	})
	
	output$labelP <-renderPrint({
	
		part.dtw <-tsclust(my_dtw(),type='p',k=k(),seed=12345,iter.max=itermax(),centroid=centroid())
		part.dtw@cluster
	})
	
	output$labelHC <-renderPrint({
		hc.dtw <-tsclust(my_dtw(),type='h',k=k(),control=hierarchical_control(method=method()))
		hc.dtw@cluster
	
	})
}

shinyApp(ui = ui, server = server)
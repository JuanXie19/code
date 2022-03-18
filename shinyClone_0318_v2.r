library(shiny)

library(DT)

ui <-fluidPage(
  fluidRow(
    column(6,DT::dataTableOutput('clonetable')),
    column(6,plotOutput('UMAP'))
  )
)

server <-function(input,output,session){
  output$clonetable = DT::renderDataTable(cdr3.df,server=FALSE)
  
  observe({
    req(input$clonetable_rows_selected)
    s = input$clonetable_rows_selected
    output$UMAP <-renderPlot({
      if (length(s)==1){
        CLONE <- rownames(cdr3.df)[s]
        INDEX <- which(G.df$cdr3 == CLONE)
        G.df.highlight <- G.df[INDEX,]
        ggplot(G.df,aes(UMAP1,UMAP2,color=G))+geom_point(alpha=0.4,shape=1,size=0.6)+
          geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=2,shape=17)+
          theme_bw()+labs(title=CLONE)
      }elseif (length(s)>1){		
        par(mfrow = c(2,2))
        for (j in 1:length(s)){
          
          CLONE <- rownames(cdr3.df)[s[j]]
          INDEX <- which(G.df$cdr3 == CLONE)
          G.df.highlight <- G.df[INDEX,]
          print(
            ggplot(G.df,aes(UMAP1,UMAP2,color=G))+geom_point(alpha=0.4,shape=1,size=0.6)+
              geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=2,shape=17)+
              theme_bw()+labs(title=CLONE)
          )
        }		
      }else{
        ggplot(G.df,aes(UMAP1,UMAP2,color=G))+geom_point(alpha=0.4,shape=1,size=0.6)
      }
      
    })
    
    
    
  })
  
  
  
}

## create shiny object
shinyApp(ui = ui, server = server)

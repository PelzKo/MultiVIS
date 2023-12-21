#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(data.table)
library(shinyFiles)
library(shinyalert)

options(shiny.maxRequestSize=200*1024^2, scipen=999)


# Define UI for application that draws a histogram
ui <- fluidPage(
  tabsetPanel(
    tabPanel("Overview",
             shinyDirButton("dir", "Chose directory", "Upload"),
             dataTableOutput("pval_list")
             ),
    tabPanel("Detail",
             sidebarPanel(
               selectInput(inputId = "method", label = "Visualization method", choices = c("TSNE 2D", "TSNE 3D", "UMAP 2D", "UMAP 3D")),
               selectInput(inputId = "enrichment", label = "Color based on", choices = c("type", "biodomain", "score", "trait_assoc", "49", "24", "12", "chromosome")),
               conditionalPanel(condition = "input.enrichment == 'biodomain'",
                                selectInput(inputId = "biodom", label = "Biodomain", choices = c("LOADING"))),
               conditionalPanel(condition = "input.enrichment == 'score'",
                                selectInput(inputId = "scores", label = "Score", choices = c("LOADING"))),
               selectInput(inputId = "dataset", label = "Dataset", choices = c("gm_brain_set","gm_brain_excl_coreg_set")),
               selectInput(inputId = "dims", label = "Dimensions", choices = c(25, 75, 105, 128, 150), selected = 128),
               selectInput(inputId = "nodes", label = "Closest nodes", choices = c(200, 300, 400), selected = 300),
               selectInput(inputId = "samples", label = "Negative samples", choices = c(5, 10, 20), selected = 10),
               selectInput(inputId = "learning", label = "Learning rate", choices = c(0.001, 0.01, 0.1), selected = 0.01),
               selectInput(inputId = "steps", label = "Steps", choices = c(100000, 100000000), selected = 100000000)
             ),
             mainPanel(plotlyOutput("plot")))
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
    wd_root <- "C:/Users/Konstantin/Downloads"
    #wd_root <- "D:"
    shinyDirChoose(input, 'dir', roots = c(wd=wd_root))
    wd <- reactive({
      if ("list" %in% class(input$dir)){
        return(paste0(wd_root, paste(input$dir$path, collapse = "/")))
      }
      return("")
    })
    
    load_data <- reactive({
      path <- paste0(wd(), "/defcut_anno_", input$dataset, "_", input$dims, "_", input$nodes, "_",
                     input$samples, "_", input$learning, "_", prettyNum(input$steps, scientific=FALSE),
                     "_clust_proj.tsv")
      if (!file.exists(path)){
        # If steps are 200,000,000 or 400,000,000 it is expected that they just have a
        # specific combination of parameters
        if (input$steps==200000000 | input$steps==400000000){
          shinyalert("File not found", paste0(path, " was not found on your system. ",
                                              "This might be because steps 200,000,000 ",
                                              "and 400,000,000 can only be used with ",
                                              "dim=150, nodes=200, samples=5, ",
                                              "learning=0.01 or learning=0.001"),
                     type = "error")
        } else {
          shinyalert("File not found", paste0(path, " was not found on your system."),
                     type = "error")
        }
        return(plot_ly())
      }
      data <- data.table::fread(path)
      
      # Clean the data a bit: Split biodomains and transform them and cluster information into factors
      
      data[biodomain=="",biodomain:="No_BIODOM"]
      split <- data[,tstrsplit(biodomain,"\\|")]
      split[, individual := 1:nrow(data)]
      biodom <- melt(split, id.vars = 'individual',na.rm = TRUE)[, variable := NULL]
      biodom <- biodom[order(individual)]
      biodom_one_hot <- dcast(biodom, individual ~ value, fun = length)
      data <- cbind(data, biodom_one_hot)[, individual := NULL]
      data[No_BIODOM==1,unique(biodom$value)] <- NA
      data[is.na(No_BIODOM), No_BIODOM:=1]
      # Columns have to be transformed into characters so they can be overwritten
      data[ ,         # Modifying classes of particular variables
                   (unique(biodom$value)) := lapply(.SD, as.character),
                   .SDcols = unique(biodom$value)]
      for(col in unique(biodom$value)) set(data, i=which(data[[col]]==0), j=col, value="Not part of Biodom")
      for(col in unique(biodom$value)) set(data, i=which(data[[col]]==1), j=col, value="Part of Biodom")
      for(col in unique(biodom$value)) set(data, i=which(is.na(data[[col]])), j=col, value="No Info")
      
      factor_cols <- c("12", "24", "49", unique(biodom$value))
      data[ ,         # Modifying classes of particular variables
            (factor_cols) := lapply(.SD, as.factor),
            .SDcols = factor_cols]
      # Join with additional scores
      ensembl_mapping <- fread(paste0(wd(), "/adatlas_id_to_ensembl.csv"))
      ensembl_mapping[,n.ensembl_id:=str_sub(n.ensembl_id, 3, -3)]
      new_scores <- fread(paste0(wd(), "/20231121_AD_gene_scores_updated.csv"))
      new_scores <- new_scores[match(ensembl_mapping$n.ensembl_id, sid)]
      new_scores$ad_ids <- ensembl_mapping$`ID(n)`
      new_scores <- new_scores[!is.na(sid)]
      # The names are exactly the same and will otherwise get a .y suffix
      new_scores[,name:=NULL]
      
      data <- merge(data, new_scores, by.x="id", by.y="ad_ids", all.x=TRUE)
      
      
      
      # Update the choices of biodom
      updateSelectInput(
        inputId = "biodom",
        choices = unique(biodom$value)
      )
      # Update the choices of scores
      updateSelectInput(
        inputId = "scores",
        choices = colnames(data)[grepl("score", colnames(data), fixed = TRUE) | colnames(data) %in% c("Overall", "GeneticsScore", "OmicsScore")]
      )
      return(data)
    })
    
    output$pval_list <- renderDataTable({
      path <- paste0(wd(), "/all_enrichments_defcut.tsv")
      if(!file.exists(path)){
        return(data.table())
      }
      all_p_vals <- data.table::fread(path)
      all_p_vals <- all_p_vals[parameters!="parameters"]
      numeric_cols <- colnames(all_p_vals)[!colnames(all_p_vals) %in% c('biodomain', 'parameters')]
      all_p_vals[ ,                    # Modifying classes of particular variables
                  (numeric_cols) := lapply(.SD, as.numeric),
                  .SDcols = numeric_cols]
      return(all_p_vals)
    })
    
    
    output$plot <- renderPlotly({
      z_dim <- NULL
      
      if (input$method =="TSNE 2D"){
        x_dim <- "TSNE_1_2D"
        y_dim <- "TSNE_2_2D"
        x_dim_title <- "t-SNE dimension one"
        y_dim_title <- "t-SNE dimension two"
      } else if (input$method =="TSNE 3D"){
        x_dim <- "TSNE_1_3D"
        y_dim <- "TSNE_2_3D"
        z_dim <- "TSNE_3_3D"
        x_dim_title <- "t-SNE dimension one"
        y_dim_title <- "t-SNE dimension two"
        z_dim_title <- "t-SNE dimension three"
      } else if (input$method =="UMAP 2D"){
        x_dim <- "UMAP_1_2D"
        y_dim <- "UMAP_2_2D"
        x_dim_title <- "UMAP dimension one"
        y_dim_title <- "UMAP dimension two"
      } else if (input$method =="UMAP 3D"){
        x_dim <- "UMAP_1_3D"
        y_dim <- "UMAP_2_3D"
        z_dim <- "UMAP_3_3D"
        x_dim_title <- "UMAP dimension one"
        y_dim_title <- "UMAP dimension two"
        z_dim_title <- "UMAP dimension three"
      } else {
        print(paste0("Unknown method ", input$method))
        
      }
      data <- load_data()
      if (length(unique(data[,get(x_dim)])) == 1 || is.na(unique(data[,get(x_dim)]))){
        shinyalert("No values", paste0("There was an error in calculating the values for ",
                                       input$method, ". Please choose another visualization method."),
                   type = "error")
        return(plot_ly())
      }
      
      color_by <- input$enrichment
      if(input$enrichment == "biodomain"){
        color_by <- input$biodom
      }
      if (nrow(unique(data[,..color_by]))>350){
        shinyalert("Too many groups", paste0("You are trying to color ", nrow(unique(data[,..color_by])),
                                             " different groups which is more than the maximum of 350. ",
                                             "Please choose a different property to color by."),
                   type = "error")
        return(plot_ly())
      }
      
      if (is.null(z_dim)){
        fig <- plot_ly(data, x = ~get(x_dim), y = ~get(y_dim), color = ~get(color_by))#, colors = c('#BF382A', '#0C4B8E'))
        fig <- fig %>% add_markers(marker = list(size = 2))
        fig <- fig %>% layout(xaxis = list(title = x_dim_title),
                              yaxis = list(title = y_dim_title),
                              legend=list(title=list(text=paste0("<b> ", color_by, " </b>"))))
      } else {
        fig <- plot_ly(data, x = ~get(x_dim), y = ~get(y_dim), z = ~get(z_dim), color = ~get(color_by), marker = list(size = 10))#, colors = c('#BF382A', '#0C4B8E'))
        fig <- fig %>% add_markers(marker = list(size = 3))
        fig <- fig %>% layout(scene = list(xaxis = list(title = x_dim_title),
                                           yaxis = list(title = y_dim_title),
                                           zaxis = list(title = z_dim_title)),
                              legend=list(title=list(text=paste0("<b> ", color_by, " </b>"))))
      }
      return(fig)
    })
    
    
    #output$biodom_select <- renderUI({
    #  if(input$enrichment == "biodomain"){
    #    selectInput(inputId = "biodom", label = "Biodomain", choices = biodomains())
    #  } else {
    #    NULL
    #  }
    #})
}

# Run the application 
shinyApp(ui = ui, server = server)

############################
######## Libraries #########
############################

library(shiny)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(plotly)
library(ggrepel)
library(kohonen)
library(effectsize)
library(tibble)
library(tidyr)
library(rMQanalysis)
library(shinyBS)
library(DT)
library(mailtoR)
library(viridis)
library(pheatmap)
library(ggpointdensity)
set.seed(666)

############################
####### Variables ##########
############################

load("./InputFiles/01_LInfantum_Data.RData")
load("./InputFiles/02_LMajor_Data.RData")
load("./InputFiles/03_LMexicana_Data.RData")
load("./InputFiles/04_Clustering_Data.RData")
load("./InputFiles/05_Orthology_Data.RData")

############################
####### Functions ##########
############################

createLink <- function(species, val) {
  if (species == "Mouse") {
    val_sub <- gsub(pattern = "\\..*",replacement = "",x = val)
    sprintf('<a href="https://mart.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=%s" target="_blank" class="btn btn-link" role="button">Ensembl</a>',val_sub)
  }else{
    sprintf('<a href="https://tritrypdb.org/tritrypdb/app/record/gene/%s" target="_blank" class="btn btn-link" role="button">TriTrypDB</a>',val)
  }
  
}

############################
######### Script ###########
############################

server <- function(input, output, session) {
  
  #### Data functions ####
  
  # Modify the different data sets in function of
  # user's input
  
  ##### L. infantum ####
  
  ## Main data table
  
  infantum_MainDataTable1 <- reactive({
    my_main_table_data <- df_infantum 
    if(input$infantum_ui_datasetID != paste0(unique(df_infantum$Species),collapse = " & ")) {
      my_main_table_data <- df_infantum %>% filter(Species == input$infantum_ui_datasetID)
    }
    if(length(input$infantum_ui_timepointID) > 0) {
      my_main_table_data <- my_main_table_data %>% filter(TimePoint %in% input$infantum_ui_timepointID)
    }
    my_main_table_data
  })
  
  infantum_MainDataTable2 <- reactive({
    my_main_table_data <- infantum_MainDataTable1()
    my_main_table_data <- my_main_table_data %>% filter(Majority.protein.IDs %in% input$infantum_ui_selectedIDs) 
    my_main_table_data
  })
  
  infantum_MainDataTable3 <- reactive({
    my_main_table_data <- infantum_MainDataTable1()
    my_main_table_data$TimePoint <- factor(x = my_main_table_data$TimePoint, levels = unique(my_main_table_data$TimePoint))
    my_main_table_data <- my_main_table_data %>%
      group_by(Gene.Name,Majority.protein.IDs,TimePoint,Species)%>%
      summarize("mean(LFQ)" = mean(Value, na.rm=TRUE))
    my_main_table_data$Database  <- "Jnk"
    my_main_table_data[my_main_table_data$Species == "M. musculus",]$Database <- createLink("Mouse", my_main_table_data[my_main_table_data$Species == "M. musculus",]$Majority.protein.IDs)
    my_main_table_data[my_main_table_data$Species != "M. musculus",]$Database <- createLink("Leishmania", my_main_table_data[my_main_table_data$Species != "M. musculus",]$Majority.protein.IDs)
    colnames(my_main_table_data)[1] <- "Gene Name"
    colnames(my_main_table_data)[2] <- "Protein ID"
    my_main_table_data
  })
  
  ## Test data table
  
  infantum_TestDataTable1 <- reactive({
    my_test_table_data <- test_infantum 
    if(input$infantum_ui_datasetID != paste0(unique(df_infantum$Species),collapse = " & ")) {
      my_test_table_data <- test_infantum %>% filter(Species == input$infantum_ui_datasetID)
    }
    if(length(input$infantum_ui_timepointID) > 0) {
      my_test_table_data <- my_test_table_data %>% filter(group1 %in% input$infantum_ui_timepointID)
      my_test_table_data <- my_test_table_data %>% filter(group2 %in% input$infantum_ui_timepointID)
    }
    my_test_table_data
  })
  
  infantum_TestDataTable2 <- reactive({
    my_test_table_data <- infantum_TestDataTable1()
    my_test_table_data <- my_test_table_data %>% filter(Majority.protein.IDs %in% input$infantum_ui_selectedIDs) 
    my_test_table_data
  })
  
  ## Exploratory data
  
  infantum_FullDataTable1 <- reactive({
    my_full_table_data <- full_infantum 
    if(input$infantum_ui_datasetID != paste0(unique(df_infantum$Species),collapse = " & ")) {
      my_full_table_data <- full_infantum %>% filter(Species == input$infantum_ui_datasetID)
    }
    if(length(input$infantum_ui_timepointID) > 0) {
      
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
      my_pattern <- paste0(collapse = "|",  
                           c(paste0(collapse = "|",
                                    colnames(full_infantum)[1:11]),
                             paste0(collapse = "|",
                                    infantum_time_df_ui$time_point)))
      
      my_full_table_data <- my_full_table_data[,grep(pattern = my_pattern,x = colnames(my_full_table_data))]
    
    }
    my_full_table_data
  })

  infantum_FullDataTable2 <- reactive({
    full_infantum2 <- infantum_FullDataTable1()
    if (length(input$infantum_ui_selectedIDs) != 0) {
      full_infantum2 <- full_infantum2 %>% filter(Majority.protein.IDs %in% input$infantum_ui_selectedIDs)
      full_infantum2
    }else{
      full_infantum2 <- data.frame()
    }
  })
  
  ##### L. major ####
  
  ## Main data table
  
  major_MainDataTable1 <- reactive({
    my_main_table_data <- df_major 
    if(input$major_ui_datasetID != paste0(unique(df_major$Species),collapse = " & ")) {
      my_main_table_data <- df_major %>% filter(Species == input$major_ui_datasetID)
    }
    if(length(input$major_ui_timepointID) > 0) {
      my_main_table_data <- my_main_table_data %>% filter(TimePoint %in% input$major_ui_timepointID)
    }
    my_main_table_data
  })
  
  major_MainDataTable2 <- reactive({
    my_main_table_data <- major_MainDataTable1()
    my_main_table_data <- my_main_table_data %>% filter(Majority.protein.IDs %in% input$major_ui_selectedIDs) 
    my_main_table_data
  })
  
  major_MainDataTable3 <- reactive({
    my_main_table_data <- major_MainDataTable1()
    my_main_table_data$TimePoint <- factor(x = my_main_table_data$TimePoint, levels = unique(my_main_table_data$TimePoint))
    my_main_table_data <- my_main_table_data %>%
      group_by(Gene.Name,Majority.protein.IDs,TimePoint,Species)%>%
      summarize("mean(LFQ)" = mean(Value, na.rm=TRUE))
    my_main_table_data$Database  <- "Jnk"
    my_main_table_data[my_main_table_data$Species == "M. musculus",]$Database <- createLink("Mouse", my_main_table_data[my_main_table_data$Species == "M. musculus",]$Majority.protein.IDs)
    my_main_table_data[my_main_table_data$Species != "M. musculus",]$Database <- createLink("Leishmania", my_main_table_data[my_main_table_data$Species != "M. musculus",]$Majority.protein.IDs)
    colnames(my_main_table_data)[1] <- "Gene Name"
    colnames(my_main_table_data)[2] <- "Protein ID"
    my_main_table_data
  })
  
  ## Test data table
  
  major_TestDataTable1 <- reactive({
    my_test_table_data <- test_major 
    if(input$major_ui_datasetID != paste0(unique(df_major$Species),collapse = " & ")) {
      my_test_table_data <- test_major %>% filter(Species == input$major_ui_datasetID)
    }
    if(length(input$major_ui_timepointID) > 0) {
      my_test_table_data <- my_test_table_data %>% filter(group1 %in% input$major_ui_timepointID)
      my_test_table_data <- my_test_table_data %>% filter(group2 %in% input$major_ui_timepointID)
    }
    my_test_table_data
  })
  
  major_TestDataTable2 <- reactive({
    my_test_table_data <- major_TestDataTable1()
    my_test_table_data <- my_test_table_data %>% filter(Majority.protein.IDs %in% input$major_ui_selectedIDs) 
    my_test_table_data
  })
  
  ## Exploratory data
  
  major_FullDataTable1 <- reactive({
    my_full_table_data <- full_major 
    if(input$major_ui_datasetID != paste0(unique(df_major$Species),collapse = " & ")) {
      my_full_table_data <- full_major %>% filter(Species == input$major_ui_datasetID)
    }
    if(length(input$major_ui_timepointID) > 0) {
      
      major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
      my_pattern <- paste0(collapse = "|",  
                           c(paste0(collapse = "|",
                                    colnames(full_major)[1:11]),
                             paste0(collapse = "|",
                                    major_time_df_ui$time_point)))
      
      my_full_table_data <- my_full_table_data[,grep(pattern = my_pattern,x = colnames(my_full_table_data))]
      
    }
    my_full_table_data
  })
  
  major_FullDataTable2 <- reactive({
    full_major2 <- major_FullDataTable1()
    if (length(input$major_ui_selectedIDs) != 0) {
      full_major2 <- full_major2 %>% filter(Majority.protein.IDs %in% input$major_ui_selectedIDs)
      full_major2
    }else{
      full_major2 <- data.frame()
    }
  })
  
  ##### L. mexicana ####
  
  ## Main data table
  
  mexicana_MainDataTable1 <- reactive({
    my_main_table_data <- df_mexicana 
    if(input$mexicana_ui_datasetID != paste0(unique(df_mexicana$Species),collapse = " & ")) {
      my_main_table_data <- df_mexicana %>% filter(Species == input$mexicana_ui_datasetID)
    }
    if(length(input$mexicana_ui_timepointID) > 0) {
      my_main_table_data <- my_main_table_data %>% filter(TimePoint %in% input$mexicana_ui_timepointID)
    }
    my_main_table_data
  })
  
  mexicana_MainDataTable2 <- reactive({
    my_main_table_data <- mexicana_MainDataTable1()
    my_main_table_data <- my_main_table_data %>% filter(Majority.protein.IDs %in% input$mexicana_ui_selectedIDs) 
    my_main_table_data
  })
  
  mexicana_MainDataTable3 <- reactive({
    my_main_table_data <- mexicana_MainDataTable1()
    my_main_table_data$TimePoint <- factor(x = my_main_table_data$TimePoint, levels = unique(my_main_table_data$TimePoint))
    my_main_table_data <- my_main_table_data %>%
      group_by(Gene.Name,Majority.protein.IDs,TimePoint,Species)%>%
      summarize("mean(LFQ)" = mean(Value, na.rm=TRUE))
    my_main_table_data$Database  <- "Jnk"
    my_main_table_data[my_main_table_data$Species == "M. musculus",]$Database <- createLink("Mouse", my_main_table_data[my_main_table_data$Species == "M. musculus",]$Majority.protein.IDs)
    my_main_table_data[my_main_table_data$Species != "M. musculus",]$Database <- createLink("Leishmania", my_main_table_data[my_main_table_data$Species != "M. musculus",]$Majority.protein.IDs)
    colnames(my_main_table_data)[1] <- "Gene Name"
    colnames(my_main_table_data)[2] <- "Protein ID"
    my_main_table_data
  })
  
  ## Test data table
  
  mexicana_TestDataTable1 <- reactive({
    my_test_table_data <- test_mexicana 
    if(input$mexicana_ui_datasetID != paste0(unique(df_mexicana$Species),collapse = " & ")) {
      my_test_table_data <- test_mexicana %>% filter(Species == input$mexicana_ui_datasetID)
    }
    if(length(input$mexicana_ui_timepointID) > 0) {
      my_test_table_data <- my_test_table_data %>% filter(group1 %in% input$mexicana_ui_timepointID)
      my_test_table_data <- my_test_table_data %>% filter(group2 %in% input$mexicana_ui_timepointID)
    }
    my_test_table_data
  })
  
  mexicana_TestDataTable2 <- reactive({
    my_test_table_data <- mexicana_TestDataTable1()
    my_test_table_data <- my_test_table_data %>% filter(Majority.protein.IDs %in% input$mexicana_ui_selectedIDs) 
    my_test_table_data
  })
  
  ## Exploratory data
  
  mexicana_FullDataTable1 <- reactive({
    my_full_table_data <- full_mexicana 
    if(input$mexicana_ui_datasetID != paste0(unique(df_mexicana$Species),collapse = " & ")) {
      my_full_table_data <- full_mexicana %>% filter(Species == input$mexicana_ui_datasetID)
    }
    if(length(input$mexicana_ui_timepointID) > 0) {
      
      mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
      my_pattern <- paste0(collapse = "|",  
                           c(paste0(collapse = "|",
                                    colnames(full_mexicana)[1:11]),
                             paste0(collapse = "|",
                                    mexicana_time_df_ui$time_point)))
      
      my_full_table_data <- my_full_table_data[,grep(pattern = my_pattern,x = colnames(my_full_table_data))]
      
    }
    my_full_table_data
  })
  
  mexicana_FullDataTable2 <- reactive({
    full_mexicana2 <- mexicana_FullDataTable1()
    if (length(input$mexicana_ui_selectedIDs) != 0) {
      full_mexicana2 <- full_mexicana2 %>% filter(Majority.protein.IDs %in% input$mexicana_ui_selectedIDs)
      full_mexicana2
    }else{
      full_mexicana2 <- data.frame()
    }
  })
  
  ##### Clustering ####
  
  ## Main data table
  
  clustering_MainDataTable1 <- reactive({
    
    infantum_data_table <- infantum_MainDataTable3()
    
    infantum_data_table$Experiment <- "L. infantum"
    
    major_data_table <- major_MainDataTable3()
    
    major_data_table$Experiment <- "L. major"
    
    mexicana_data_table <- mexicana_MainDataTable3()
    
    mexicana_data_table$Experiment <- "L. mexicana"
    
    clustering_data_table <- rbind(infantum_data_table, major_data_table, mexicana_data_table)
    
    if (input$clustering_ui_datasetID == "Leishmania spp.") {
      
      clustering_data_table <- clustering_data_table %>% filter(`Protein ID` %in% cluster_leishmania_df$`Protein ID`)
      
      clustering_data_table <- left_join(clustering_data_table, cluster_leishmania_df, by= "Protein ID")
      
      clustering_data_table <- arrange(clustering_data_table, Cluster, `Gene Name`)
      
      # clustering_data_table <- clustering_data_table %>% dplyr::select(!c(Experiment))
      
      clustering_data_table <- clustering_data_table[,c(8,1,2,9,7,3:6)]

    }
    
    if (input$clustering_ui_datasetID == "M. musculus") {
      
      clustering_data_table$Filter_ID <- paste0(clustering_data_table$`Protein ID`, "_", clustering_data_table$Experiment)
      
      clustering_data_table <- clustering_data_table %>% filter(`Filter_ID` %in% cluster_mouse_df$`Filter_ID`)
      
      clustering_data_table <- left_join(clustering_data_table, cluster_mouse_df, by= "Filter_ID")
      
      clustering_data_table <- arrange(clustering_data_table, Cluster, `Gene Name`)
      
      clustering_data_table <- clustering_data_table %>% dplyr::select(!c(Filter_ID))
      
      clustering_data_table <- clustering_data_table[,c(8,1,2,7,3:6)]
      
    } 
    
    clustering_data_table
    
  })
  
  clustering_MainDataTable2 <- reactive({
    clustering_data_table <- clustering_MainDataTable1()
    clustering_data_table <- clustering_data_table %>% filter(`Protein ID` %in% input$clustering_ui_selectedIDs) 
    clustering_data_table
  })
  
  #### Output functions ####
  
  #### L. infantum #### 
  
  #### Main table output  
  
  output$infantum_MainDataTable <- renderDataTable({
    infantum_MainDataTable3() %>%
      datatable(options = list(pageLength = 7, autoWidth = TRUE), 
                rownames = F,escape = FALSE) %>%
      formatRound("mean(LFQ)", digits = 2)
  })
  
  infantum_MainDataTable_proxy <- dataTableProxy('infantum_MainDataTable')
  
  #### Exploratory analysis output 
  
  output$infantum_PCA <- renderPlot({
    
    # Load and format the data
    
    if(length(input$infantum_ui_timepointID) > 0) {
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
    }else{
      infantum_time_df_ui <- infantum_time_df
    }
    plot_df_infantum <- infantum_FullDataTable1()
    pca_columns <- grep(pattern = "^imputed.log2.LFQ.intensity.", colnames(plot_df_infantum))
    pca <- prcomp(t(na.omit(plot_df_infantum[,pca_columns])), scale.=TRUE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    scores <- data.frame(pca$x[,c("PC1","PC2")])
    scores$PC2 <- scores$PC2 * -1
    scores$TimePoint <- factor(rep(infantum_time_df_ui$hour, each = 4), levels = infantum_time_df_ui$hour)
    scores$Label <- gsub(".*_", replacement = "", x = gsub(pattern = "^imputed.log2.LFQ.intensity.",replacement = "",rownames(scores)))
    
    if (length(infantum_time_df_ui$hour)>1) {
      colors_infantum <- brewer.pal(n = 9,name = "Greys")[3:(length(infantum_time_df_ui$hour)+2)]
    }else{
      colors_infantum <- brewer.pal(n = 9,name = "Greys")[3:9]
    }
    
    # Input plot
    
    infantum_PCA <- ggplot(data=scores,aes(x=PC1, y=PC2, colour = TimePoint))+
      geom_hline(yintercept = 0,color = "#A3C585", linewidth = 1)+
      geom_vline(xintercept = 0, color = "#A3C585", linewidth = 1)+
      stat_ellipse(geom = "polygon", alpha = 0.15, aes(fill = TimePoint))+
      geom_point(size=3.5)+
      geom_text_repel(mapping = aes(label = Label), 
                      size = 5,
                      fontface = 'bold', 
                      color = 'black')+
      theme_minimal()+
      xlab(paste("PC1: ", round(percentVar[1] * 100,digits = 2), "%"))+
      ylab(paste("PC2: ", round(percentVar[2] * 100, digits = 2), "%"))+
      scale_y_continuous(limits = c(-150,175),breaks = seq(-150, 150, by = 50))+
      # scale_x_continuous(limits = c(-100,100),breaks = seq(-100, 100, by = 25))+
      scale_color_manual(name = "Hours\npost-infection", values = colors_infantum)+
      scale_fill_manual(name = "Hours\npost-infection", values = colors_infantum)+
      guides(colour = guide_legend(nrow = 1))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
      theme(legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 14),
            legend.key.size = unit(1, 'cm'),
            legend.position = "top")
    
    print(infantum_PCA)
    
  })
  
  output$infantum_HeatMap <- renderPlot({
    
    # Load and format the data
    
    if(length(input$infantum_ui_timepointID) > 0) {
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
    }else{
      infantum_time_df_ui <- infantum_time_df
    }
    plot_df_infantum <- infantum_FullDataTable1()
    lfq <- plot_df_infantum[, grep("^imputed.log2.LFQ.intensity.", names(plot_df_infantum))]
    colnames(lfq) <- gsub("^imputed.log2.LFQ.intensity.", "", colnames(lfq))
    pattern <- unique(gsub(pattern = "_.*", replacement = "", colnames(lfq)))
    
    for (i in 1:length(infantum_time_df_ui$hour)) {
      
      colnames(lfq) <- gsub(pattern = pattern[i], replacement = infantum_time_df_ui$hour[i], colnames(lfq))
      
    }
    
    group <- as.factor(gsub("_[1-4]$", "", colnames(lfq)))
    lfq.avg <- t(apply(lfq, 1, function(x) tapply(x, group, mean,na.rm = T)))
    v <- apply(lfq, 1, sd)
    
    corr <- apply(lfq[rev(order(v)), ], 2, function(x) {
      apply(lfq[rev(order(v)), ], 2, function(y) {
        cor(x, y, use = "na.or.complete")
      })
    })
    
    if (length(infantum_time_df_ui$hour)>0) {
      
      Exp <- rep("L. infantum", length(infantum_time_df_ui$hour), each = 4)
      HPI <- rep(as.character(infantum_time_df_ui$hour), each = 4)
      
    }else{
      
      Exp <- rep("L. infantum", length(as.character(hours)), each = 4)
      HPI <- rep(as.character(hours), each = 4)
      
    }
    
    df <- data.frame(HPI,Exp)
    
    rownames(df) <- colnames(corr)
    
    Exp <- "#A3C585"
    names(Exp) <- unique(df$Exp)
    
    HPI <- brewer.pal(n = 9,name = "Greys")[3:(length(infantum_time_df_ui$hour)+2)]
    names(HPI) <- unique(df$HPI)
    
    anno_colors <- list(HPI = HPI,
                        Exp = Exp)
    
    # Input plot
    
    hm_pal <- viridis(n = 15,option = "G",direction = -1)

    # infantum_HeatMap <- gplots::heatmap.2(corr,
    #                                      trace="none",
    #                                      Colv=T,
    #                                      Rowv=T,
    #                                      dendrogram="column",
    #                                      density.info="density",
    #                                      srtCol=45,
    #                                      margins=c(10, 10),
    #                                      key.xlab="Pearson's correlation coefficient",
    #                                      key.title="",
    #                                      keysize=1.5,
    #                                      col=hm_pal)
    
    infantum_heatmap <- pheatmap(corr,
                                 col=hm_pal,
                                 cluster_cols = F,
                                 cluster_rows = F,
                                 show_rownames = F,
                                 show_colnames = F,
                                 # trace = "none",
                                 annotation_col = df,
                                 annotation_row = df,
                                 annotation_colors = anno_colors,
                                 annotation_names_row = T,
                                 annotation_legend = F,
                                 border_color=NA,
                                 fontsize = 14)
    dev.off()
    print(infantum_heatmap)
                
  })
  
  #### Dynamicity analysis output
  
  output$infantum_ScatterPlot <- renderPlotly({
    
  if (length(input$infantum_ui_timepointID) > 1) {
    
    # Settings:
    
    viridis_option <- "G"
    
    exp_color <- "#A3C585"
    
    # Proteins to highlight:
    
    if(length(input$infantum_ui_DynamicControls) == 2) {
      
      highlight <- c("LINF_050008500", "LINF_230005400", "LINF_290017500", 
                     "LINF_040017600",
                     "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126",
                     "ENSMUSG00000004677")
    }
    
    if(length(input$infantum_ui_DynamicControls) == 1) {
      
      if(input$infantum_ui_DynamicControls == "Dynamic proteins") {
        
        highlight <- c("LINF_050008500", "LINF_230005400", "LINF_290017500",
                       "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126")
      }
      
      if(input$infantum_ui_DynamicControls == "Stable proteins") {
        
        highlight <- c("LINF_040017600",
                       "ENSMUSG00000004677")
      }
      
    }
    
    if(length(input$infantum_ui_DynamicControls) == 0) {
      
      highlight <- NULL
    }
    
    UI_ID <- input$infantum_ui_selectedIDs
    
    highlight <- c(highlight, unique(UI_ID))
    
    highlight <- ifelse(test = grepl(pattern = "ENSMUS",x = highlight),
                        yes = gsub(pattern = "\\..*", replacement = "",highlight),
                        no = highlight) 
    
    # gini function
    
    gini <- function(x) {
      x <- x - min(df_means)
      sum(sapply(1:length(x), function(i) sum(abs(x[i] - x)))) / (2 * length(x)^2 * mean(x))
    }
    
    # gini threshold
    
    GINI_THRES <- input$infantum_Dynamicity
    
    # LFQ threshold
    
    LFQ_THRES <- input$infantum_LFQ
    
    # Data
    
    df <- infantum_FullDataTable1()
    
    rownames(df) <- df$my_label
    
    # Retrieve the mean columns & tidy
    
    df_means <- df[,grep(pattern = "mean",x = colnames(df))]
    
    if(length(input$infantum_ui_timepointID) > 0) {
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
    }else{
      infantum_time_df_ui <- infantum_time_df
    }
    
    colnames(df_means) <- unique(as.character(infantum_time_df_ui$hour))
    
    # Compute gini score
    
    avg.lfq     <- apply(df_means, 1, mean)
    max.lfq     <- apply(df_means, 1, max)
    df_means$gini <- apply(df_means, 1, gini)
    
    
    x <- avg.lfq
    y <- setNames(df_means$gini, rownames(df_means))
    
    x_df <- enframe(x = x)
    colnames(x_df)[2] <- "LFQ"
    y_df <- enframe(x = y)
    colnames(y_df)[2] <- "Dynamicity"
    
    plot_df <- left_join(x = x_df, y = y_df, by = "name")
    
    highlight_df <- data.frame(name = highlight)
    
    plot_df$Highlight <- ifelse(test = grepl(pattern = "ENSMUS",x = plot_df$name),
                                yes = gsub(pattern = "\\..*", replacement = "",plot_df$name),
                                no = plot_df$name) %in% highlight
    
    plot_df$DynamicityThreshold <- plot_df$Dynamicity > GINI_THRES
    
    plot_df$LFQThreshold <- plot_df$LFQ > LFQ_THRES
    
    plot_df$`Threshold (D.L)` <- interaction(plot_df$DynamicityThreshold, plot_df$LFQThreshold)
    
    plot_df <- plot_df[rev(order(plot_df$`Threshold (D.L)`)),]
    
    df_to_join <- df_infantum[,1:2]
    
    df_to_join <- unique(df_to_join)
    
    colnames(df_to_join) <- c("gene_name", "name")
    
    plot_df <- left_join(x = plot_df, y = df_to_join, by = "name")
    
    plotly_string <- paste0(sprintf('%s',
                                    sub('(.{30}).*', '\\1...', plot_df$gene_name)))
    
    plot_df$`Gene name` <- plotly_string
    
    infantum_scatter_plot <- ggplot(plot_df, aes(x=LFQ, y=Dynamicity, label =`Gene name`)) +
      geom_hline(yintercept = GINI_THRES, linetype = "longdash", color = exp_color, linewidth = 1)+
      geom_vline(xintercept = LFQ_THRES, linetype = "longdash", color = exp_color, linewidth = 1)+
      geom_pointdensity(data = subset(plot_df, Highlight == FALSE), adjust = .1,alpha = .3)+
      geom_point(data = subset(plot_df, Highlight == TRUE), color ="blue", size = 3)+
      scale_color_viridis(option = viridis_option, limits = c(0,1250),breaks = seq(0, 1250, by = 250))+
      scale_y_continuous(limits = c(0,0.6),breaks = seq(0, 0.5, by = 0.25))+
      scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
      xlab("Intensity")+
      annotate('text', label=paste(sum(y>GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
      annotate('text', label=paste(sum(y<GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=min(plot_df$Dynamicity), size = 7) +
      annotate('text', label=paste(sum(y>GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
      annotate('text', label=paste(sum(y<GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=min(plot_df$Dynamicity) +0.05, size = 7) +
      theme_minimal()+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
      theme(legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 14),
            legend.key.size = unit(1, 'cm'),
            legend.position = "none")
    
    infantum_scatter_plot <- ggplotly(infantum_scatter_plot)
    
    print(infantum_scatter_plot)
    
  }else{
    
    infantum_scatter_plot <- ggplot()+
      scale_y_continuous(limits = c(0,0.6),breaks = seq(0, 0.5, by = 0.25))+
      scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
      xlab("Intensity")+
      ylab("Dynamicity")+
      theme_minimal()+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
      theme(legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 14),
            legend.key.size = unit(1, 'cm'),
            legend.position = "none")
    
    print(infantum_scatter_plot)
    
  }

  })
  
  output$infantum_LinePlot <- renderPlot({
    
    if (length(input$infantum_ui_timepointID) > 1) {
      
      # Settings:
      
      exp_color <- "#A3C585"
      
      species <- "L. infantum"
      
      # Proteins to highlight:
      
      if(length(input$infantum_ui_DynamicControls) == 2) {
        
        highlight <- c("LINF_050008500", "LINF_230005400", "LINF_290017500", 
                       "LINF_040017600",
                       "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126",
                       "ENSMUSG00000004677")
      }
      
      if(length(input$infantum_ui_DynamicControls) == 1) {
        
        if(input$infantum_ui_DynamicControls == "Dynamic proteins") {
          
          highlight <- c("LINF_050008500", "LINF_230005400", "LINF_290017500",
                         "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126")
        }
        
        if(input$infantum_ui_DynamicControls == "Stable proteins") {
          
          highlight <- c("LINF_040017600",
                         "ENSMUSG00000004677")
        }
        
      }
      
      if(length(input$infantum_ui_DynamicControls) == 0) {
        
        highlight <- NULL
        
      }
      
      UI_ID <- input$infantum_ui_selectedIDs
      
      highlight <- c(highlight, unique(UI_ID))
      
      highlight <- ifelse(test = grepl(pattern = "ENSMUS",x = highlight),
                          yes = gsub(pattern = "\\..*", replacement = "",highlight),
                          no = highlight)
      
      if (length(highlight) > 0) {
        
        colors_lineplot <- c(rep(exp_color, times = length(which(!grepl(pattern = "ENSMUS", x = highlight)))),
                             rep("#b4b4b4", times = length(which(grepl(pattern = "ENSMUS", x = highlight)))))
        
        shapes_lineplot <- c(15,16,17,8,18)
        
        shapes_lineplot <- c(shapes_lineplot[1:length(which(!grepl(pattern = "ENSMUS", x = highlight)))],
                             shapes_lineplot[1:length(which(grepl(pattern = "ENSMUS", x = highlight)))])
        
        df <- infantum_FullDataTable1()
        
        rownames(df) <- df$my_label
        
        # Retrieve the mean columns & tidy
        
        df_means <- df[,grep(pattern = "mean",x = colnames(df))]
        
        if(length(input$infantum_ui_timepointID) > 0) {
          infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
        }else{
          infantum_time_df_ui <- infantum_time_df
        }
        
        colnames(df_means) <- unique(as.character(infantum_time_df_ui$hour))
        
        plot_df <- df_means
        
        plot_df$Highlight <- ifelse(test = grepl(pattern = "ENSMUS", x = rownames(plot_df)),
                                    yes = gsub(pattern = "\\..*", replacement = "", rownames(plot_df)),
                                    no = rownames(plot_df)) %in% highlight
        
        plot_df$name <- rownames(plot_df)
        
        df_to_join <- df_infantum[,1:2]
        
        df_to_join <- unique(df_to_join)
        
        colnames(df_to_join) <- c("gene_name", "name")
        
        plot_df <- left_join(x = plot_df, y = df_to_join, by = "name")
        
        plot_df <- plot_df[plot_df$Highlight == T,]
        
        plot_df <- plot_df %>% pivot_longer(-c(name, gene_name, Highlight), names_to = "Timepoint", values_to = "LFQ")
        
        plot_df$Species <- ifelse(test = grepl(pattern = "ENSMUS", x = plot_df$name),
                                  yes = "M. musculus",
                                  no = species)
        
        plot_df <- plot_df %>% arrange(Species,gene_name)
        
        plot_df$gene_name <- factor(x = plot_df$gene_name, levels = unique(plot_df$gene_name))
        
        plot_df$Timepoint <- factor(x = plot_df$Timepoint, levels = unique(as.character(infantum_time_df_ui$hour)))
        
        infantum_lineplot <- ggplot(plot_df, aes(x=Timepoint, y=LFQ, shape = gene_name, group = gene_name, color = gene_name))+
          facet_grid(cols = vars(Species))+
          geom_line(linewidth = 1)+
          geom_point(size = 4)+
          labs(shape = "Gene name")+
          guides(color="none")+
          # ylab(as.expression(bquote(bold(mean(log[2](LFQ)))))) +
          ylab("Intensity")+
          xlab("Hours post-infection") +
          scale_color_manual(labels = unique(plot_df$gene_name),
                             values = colors_lineplot)+
          scale_shape_manual(labels  = unique(plot_df$gene_name),
                             values = shapes_lineplot)+
          # scale_y_continuous(limits = c(24,32.5),breaks = seq(25, 32.5, by = 2.5))+
          # guides(color=FALSE)+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "top",
                strip.text = element_text(size = 20, face = "italic"))
        
        print(infantum_lineplot)
        
      }else{
        
        infantum_lineplot <- ggplot()+
          # scale_x_discrete(name = hours)+
          scale_y_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
          xlab("Hours post-infection")+
          ylab("Intensity")+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "none")
        
        print(infantum_lineplot)
        
        
      }
      
    }else{
      
      infantum_lineplot <- ggplot()+
        # scale_x_discrete(name = hours)+
        scale_y_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
        xlab("Hours post-infection")+
        ylab("Intensity")+
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      print(infantum_lineplot)
      
      showNotification(paste("Dynamicity analysis needs at least 2 timepoints!"), duration = 10,type = "error")
      
    }
    
  })
  
  #### Clustering analysis output 
  
  output$infantum_ClusterLinePlot <- renderPlot({
    
    # Fixed variables
    
    shortcut <- "inf"
    
    pg_LFQ <- infantum_FullDataTable1()
    
    hl_df <- infantum_FullDataTable2()
    
    # User input variables
    
    xdim <- input$infantum_xDIM
    
    ydim <- input$infantum_yDIM
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    if(length(input$infantum_ui_timepointID) > 0) {
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
    }else{
      infantum_time_df_ui <- infantum_time_df
    }
    
    colnames(pg_clust) <- c("Majority.Protein.IDs",infantum_time_df_ui$hour)
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- infantum_time_df_ui$hour
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_inf <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_inf$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    # This plot will always show, wether the user selects an input ID or not
    
    mycolors <- colorRampPalette(c("#000000", "#00A651"))(length(nclust))
    
    plot_list <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$TimePoint <- factor(data$TimePoint, levels = infantum_time_df_ui$hour)
      
      plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
        geom_line(color="lightgrey")+
        stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
        scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
        ylab("z-score(LFQ)")+
        ggtitle(paste0("Cluster ",i))+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(plot.title = element_text(size=22))+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      plot <- list(plot)
      
      plot_list <- c(plot_list, plot)
      
    }
    
    infantum_ClusterLinePlot <- ggarrange(plotlist = plot_list)
    
    # Input plot
    
    # User dependant plot. The IDs the user select will show up in this plot, when present.
    
    if (nrow(hl_df) > 0) {
      
      HID <- unique(hl_df$Majority.protein.IDs)
      
      plot_list <- c()
      
      for (i in 1:length(nclust)) {
        
        data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
        
        data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
        
        data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
        
        dataHL <- data[data$ID %in% HID,]
        
        if (nrow(dataHL) > 0) {
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          for (i in 1:length(HID)) {
            
            plot <- plot +
              geom_line(data=data[data$ID == HID[i],], aes(x=TimePoint, y=Value, group=ID,color = ID))+
              scale_color_brewer(name = "", palette = "Dark2")
              
            
          }
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }else{
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }
        
      }
      
      infantum_ClusterLinePlot <- ggarrange(plotlist = plot_list)
      
    }
    
    print(infantum_ClusterLinePlot)
    
  })
  
  output$infantum_ClusterBarPlot <- renderPlotly({
    
    # Fixed variables
    
    shortcut <- "inf"
    
    pg_LFQ <- infantum_FullDataTable1()
    
    hl_df <- infantum_FullDataTable2()
    
    # User input variables
    
    xdim <- input$infantum_xDIM
    
    ydim <- input$infantum_yDIM
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    if(length(input$infantum_ui_timepointID) > 0) {
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
    }else{
      infantum_time_df_ui <- infantum_time_df
    }
    
    colnames(pg_clust) <- c("Majority.Protein.IDs", infantum_time_df_ui$hour)
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- infantum_time_df_ui$hour
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_inf <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_inf$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    mycolors <- colorRampPalette(c("#000000", "#00A651"))(length(nclust))
    
    plot_data <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$Species <- "Unknown"
      
      data$Species <-
        ifelse(grepl("ENSMUS.*", data$ID),
               sub(".*", "M. musculus", data$Species),
               sub(".*", "Unknown", data$Species))
      
      data$Species <-
        ifelse(grepl("LINF.*", data$ID),
               sub(".*", "L. infantum", data$Species),
               sub("", "", data$Species))
      
      data$TimePoint <- factor(data$TimePoint, levels = infantum_time_df_ui$hour)
      
      data <- data[data$TimePoint == infantum_time_df_ui$hour[1],]
      
      data <- data %>%
        group_by(Species) %>%
        summarize(Count = n())%>% 
        mutate(perc = Count/sum(Count))
      
      data$Cluster <- as.character(i)
      
      plot_data <- rbind(plot_data, data)
      
    }
    
    infantum_ClusterBarPlot <- ggplot(plot_data, aes(fill=Species, y = Count , x= Cluster)) + 
      geom_bar(position="stack", stat="identity")+
      scale_fill_manual(values = c("L. infantum" = "#00A651",
                                   "M. musculus" = "#b4b4b4"))+
      theme_minimal()+
      ylab("Protein IDs")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(legend.title = element_text(size = 16))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    infantum_ClusterBarPlot <- ggplotly(infantum_ClusterBarPlot)
    
    print(infantum_ClusterBarPlot)
    
  })
  
  #### Pairwise analysis  time course output 
  
  output$infantum_Boxplot <- renderPlot({
    
    if(length(input$infantum_ui_timepointID) > 0) {
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
    }else{
      infantum_time_df_ui <- infantum_time_df
    }
    
    if (length(infantum_time_df_ui$hour)>1) {
      colors_infantum <- brewer.pal(n = 9,name = "Greys")[3:(length(infantum_time_df_ui$hour)+2)]
    }else{
      colors_infantum <- brewer.pal(n = 9,name = "Greys")[3:9]
    }
    
    colors_infantum_replica <- colorRampPalette(c("#FFFFFF", "#A3C585", "#000000"))(6)[2:5]
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_infantum <- infantum_MainDataTable1() 
    
    ## Selection data
    
    jitter_data_infantum <- infantum_MainDataTable2() 
    
    # Dummy plot
    
    infantum_Boxplot <-ggboxplot(data = plot_df_infantum, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
      scale_fill_manual(values = colors_infantum)+
      scale_color_brewer(palette="Dark2")+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))+
      theme(legend.title = element_text(size = 12))+
      theme(legend.text = element_text(size = 10))+
      theme(legend.key.size = unit(1, 'cm'))
    
    # Input plot
    
    hl_df <- infantum_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    if (nrow(jitter_data_infantum) > 0) {
      
      boxplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_infantum <- jitter_data_infantum[jitter_data_infantum$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        infantum_Boxplot <- ggboxplot(data = plot_df_infantum, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
          geom_jitter(data = jitter_data_infantum[jitter_data_infantum$Majority.protein.IDs == HID[i],],aes(x=TimePoint, y=Value, color = Replica),width = 0.20)+
          geom_line(data = line_data_infantum, aes(x=TimePoint, y=Value,group = 1))+
          scale_fill_manual(values = colors_infantum)+
          scale_color_manual(values = colors_infantum_replica)+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        infantum_Boxplot <- list(infantum_Boxplot)
        
        boxplot_list <- c(boxplot_list, infantum_Boxplot)
        
      }
      
      infantum_Boxplot <- ggarrange(plotlist = boxplot_list)
      
    }
    
    print(infantum_Boxplot)
  })
  
  output$infantum_Dotplot <- renderPlot({
    
    colors_infantum_replica <- colorRampPalette(c("#FFFFFF", "#A3C585", "#000000"))(6)[2:5]
    
    ## Retrieve highlight data
    
    hl_df <- infantum_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_infantum <- infantum_MainDataTable1()
    
    ## Selection data
    
    jitter_data_infantum <- infantum_MainDataTable2()
    
    test_data_infantum <- infantum_TestDataTable2()
    
    # Dummy plot
    
    infantum_Dotplot <- ggplot(data = jitter_data_infantum, aes(TimePoint, Value))+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    # Input plot
    
    if (nrow(jitter_data_infantum) > 0) {
      
      y_pos_dot <- max(jitter_data_infantum$Value)+0.2
      
      dotplot_list <- c()
      
      for (i in 1:length(HID)) {
      
      # Retrieve line data  
        
      line_data_infantum <- jitter_data_infantum[jitter_data_infantum$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
      
      # Format test data
      
      plot_test_data_infantum <- test_data_infantum[test_data_infantum$Majority.protein.IDs == HID[i],]
      plot_test_data_infantum <- plot_test_data_infantum[,c(4,5,3)]
      plot_test_data_infantum$p.adj <- format.pval(plot_test_data_infantum$p.adj, digits = 3)
        
      infantum_Dotplot <- ggdotplot(data = jitter_data_infantum[jitter_data_infantum$Majority.protein.IDs == HID[i],], x="TimePoint", y="Value", fill = "Replica",
                                   position = position_jitter(0.2))+
        geom_line(data = line_data_infantum, aes(x=TimePoint, y=Value,group = 1),alpha = 0.8)+
        stat_pvalue_manual(plot_test_data_infantum, y.position = y_pos_dot, step.increase = 0.1,label = "p.adj")+
        scale_fill_manual(values = colors_infantum_replica)+
        ylab("log2(LFQ)")+
        ggtitle(HID[i])+
        theme_minimal()+
        theme(plot.title = element_text(size=14,face = "bold"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 14))+
        theme(legend.text = element_text(size = 12))+
        theme(legend.key.size = unit(1, 'cm'))
      
      infantum_Dotplot <- list(infantum_Dotplot)
      
      dotplot_list <- c(dotplot_list, infantum_Dotplot)
      
      }
      
      infantum_Dotplot <- ggarrange(plotlist = dotplot_list)
      
    }
    
    print(infantum_Dotplot)
  })
  
  #### Pairwise analysis  reference output 
  
  output$infantum_Volcanoplot <- renderPlotly({
    
    #Input variables
    
    filter <- input$infantum_ui_datasetID
    exp_a <- input$infantum_TP
    exp_b <- input$infantum_referenceTP
    ID <- input$infantum_ui_selectedIDs
    enrich_s0 <- log2(2^as.numeric(input$infantum_LFC)) 
    enrich_pvalue <- as.numeric(input$infantum_Pvalue)
    
    # Fixed variables
    
    difference <- paste('difference', exp_a, exp_b, sep='_')
    pvalue <- paste('log10.pvalue', exp_a, exp_b, sep='_')
    val_count_a <- paste0('value_count_',exp_a) 
    val_count_b <- paste0('value_count_',exp_b)
    enriched_col <- paste('enriched', exp_a, exp_b, sep='_')
    
    min_quant_events <- 0 
    replicate_count <- 4
    enrich_c <- .05 
    volcano_threshold_linetype <- 2
    volcano_threshold_linesize <- 0.5
    difference_volcano_column <- 'mean2_naimp_meas_'
    volcano_column_type <- sub('([^[:digit:]]+).*', '\\1', difference_volcano_column)
    
    # Dummy plot
    
    ## Filter the data frame for data set user input.
    
    if (filter == paste0(unique(df_infantum$Species),collapse = " & ")) {
      pg_quant <- full_infantum
    }else{
      pg_quant <- full_infantum %>% filter(Species == filter)
    }
    
    ## Change the column pattern to user friendly reading (i. e., 0h instead of inft000).
    
    if(length(input$infantum_ui_timepointID) > 0) {
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
    }else{
      infantum_time_df_ui <- infantum_time_df
    }
    
    for (i in 1:length(infantum_time_df_ui$hour)) {
      
      colnames(pg_quant) <- gsub(pattern = infantum_time_df_ui$time_point[i], replacement = infantum_time_df_ui$hour[i], colnames(pg_quant))
      
    }
    
    ## Dummy plot to show until user selection.
    
    dummy_df <- data.frame(x = c(-5,5),
                           y = c(0,10))
    
    volcanoplot <-
      ggplot(dummy_df, aes(x=x, y=y))+
      theme_minimal()+
      ylab("-log<sub>10</sub>(p-value)") + 
      xlab("log<sub>2</sub>(Fold change)") +
      theme(legend.text = element_text(size = 12))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(plot.title = element_text(size=20))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
    
    infantum_Volcanoplot <-
      plotly::ggplotly(volcanoplot) %>%
      plotly::layout(showlegend = FALSE)
    
    # Input plot
    
    ## Upon proper user selection (if) produce the volcano plot.
    
    if (difference %in% colnames(pg_quant)) {
      
      ## Highlight user input, if any.
      
      pg_quant$highlight <- FALSE
      
      if (length(ID) != 0) {
        
        ids_searchstring <- paste(ID, collapse='|')
        
        matching_rows <- grep(ids_searchstring, as.character(pg_quant$Protein.IDs))
        
        pg_quant$highlight[matching_rows] <- TRUE
        
      }
     
      
      ## Filter count values. Standard set to 0, so they can check the dot all along the time course.
      ## In Mario's standard script is usually set to 2. But in this case, since it shown in the 
      ## box plot and in the dot plot, would be confusing not to show it here as well.
      ## Can be adjusted on the fixed variables section of this function.
      
      filter_min_quant_events <-
        pg_quant[val_count_a] >= min_quant_events |
        pg_quant[val_count_b] >= min_quant_events
      
      pg_quant <- pg_quant[filter_min_quant_events,]
      
      ## Select and highlight the enriched dots.
      
      my_enriched <- Enriched(pg_quant[c(difference, pvalue)],
                              c=enrich_c,
                              s0=enrich_s0,
                              pvalue=enrich_pvalue,
                              linetype=volcano_threshold_linetype,
                              size=volcano_threshold_linesize)
      
      pg_quant$enriched <- as.factor(my_enriched$enriched)
      
      
      my_points_nohighlight <- interaction(pg_quant[!pg_quant$highlight, c('enriched', 'highlight')])
      my_points_highlight <- interaction(pg_quant[pg_quant$highlight, c('enriched', 'highlight')])
      
      ## Generate a user friendly plotly string for the hoover.
      
      plotly_string <- paste0(sprintf('log2(foldchange): %.2f, -log10(pvalue): %.2f<br />Protein.IDs: %s<br />my_label: %s',
                                      pg_quant[[difference]], pg_quant[[pvalue]],
                                      sub('(.{30}).*', '\\1...', pg_quant$Protein.IDs),
                                      pg_quant$my_label),
                              switch('Gene.names' %in% names(pg_quant),
                                     sprintf('<br />Gene.names: %s', pg_quant$Gene.names),
                                     NULL),
                              switch('Protein.names' %in% names(pg_quant),
                                     sprintf('<br />Protein.names: %s', pg_quant$Protein.names),
                                     NULL),
                              switch('description' %in% names(pg_quant),
                                     sprintf('<br />description: %s', pg_quant$description),
                                     NULL))
      
      pg_quant$my_text <- plotly_string
      
      ## Subset the data for the volcano plot.
      
      pg_quant_subset <-
        pg_quant %>%
        rowwise() %>%
        mutate(imputation=ifelse(max(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count &
                                   min(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count,
                                 'measured',
                                 ifelse(min(range(.data[[val_count_a]], .data[[val_count_b]])) > 0,
                                        'some imputed',
                                        'ref imputed'))) %>%
        dplyr::select(Protein.IDs, my_label, matches(difference), matches(pvalue),
               my_text, highlight, enriched, imputation)
      
      ## Plot
      
      ## Select maximum x and y values of your volcano plot.
      
      max_x_value <- max(abs(pg_quant[grep('^difference', names(pg_quant))]),
                         na.rm=TRUE)
      max_y_value <- max(c(unlist(pg_quant[grep('^log10.pvalue', names(pg_quant))]),
                           10^enrich_pvalue + 0.5))
      
      volcanoplot <-
        ggplot(pg_quant_subset, aes_string(x=difference, y=pvalue, text='my_text')) +
        geom_hline(yintercept=0, color='grey80', na.rm=TRUE) +
        geom_vline(xintercept=0, color='grey80', na.rm=TRUE) +
        plot(my_enriched) +
        geom_point(data=subset(pg_quant_subset, !highlight),
                   aes(alpha=my_points_nohighlight,
                       color=my_points_nohighlight,
                       fill=my_points_nohighlight,
                       size=my_points_nohighlight),
                   na.rm=TRUE) +
        geom_point(data=subset(pg_quant_subset, highlight),
                   aes(alpha=my_points_highlight,
                       color=my_points_highlight,
                       fill=my_points_highlight,
                       size=my_points_highlight),
                   na.rm=TRUE) +
        scale_alpha_manual(values=c(`TRUE.FALSE`=.7, `FALSE.FALSE`=.4, `TRUE.TRUE`=.7, `FALSE.TRUE`=.7),
                           guide="none") +
        scale_color_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                           guide="none") +
        scale_fill_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                          guide="none") +
        scale_size_manual(values=c(`TRUE.FALSE`=1, `FALSE.FALSE`=.5, `TRUE.TRUE`=1.2, `FALSE.TRUE`=1),
                          guide="none")+
        coord_cartesian(xlim=c(-max_x_value,max_x_value), ylim=c(-.8,max_y_value)) +
        ylab("-log<sub>10</sub>(p-value)") + 
        xlab("log<sub>2</sub>(Fold change)") +
        annotate('text', label=exp_a, x=max_x_value / 2, y=-.4,size = 6) +
        annotate('text', label=exp_b, x=-max_x_value / 2, y=-.4,size = 6)+
        theme_minimal()+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      infantum_Volcanoplot <-
        plotly::ggplotly(volcanoplot,tooltip='text') %>%
        plotly::layout(showlegend = FALSE)
      
    }
    
    print(infantum_Volcanoplot)
    
  })
  
  #### L. major ####
  
  #### Main table output  
  
  output$major_MainDataTable <- renderDataTable({
    major_MainDataTable3() %>%
      datatable(options = list(pageLength = 7, autoWidth = TRUE), 
                rownames = F,escape = FALSE) %>%
      formatRound("mean(LFQ)", digits = 2)
  })
  
  major_MainDataTable_proxy <- dataTableProxy('major_MainDataTable')
  
  #### Exploratory analysis output 
  
  output$major_PCA <- renderPlot({
    
    # Load and format the data
    
    if(length(input$major_ui_timepointID) > 0) {
      major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
    }else{
      major_time_df_ui <- major_time_df
    }
    plot_df_major <- major_FullDataTable1()
    pca_columns <- grep(pattern = "^imputed.log2.LFQ.intensity.", colnames(plot_df_major))
    pca <- prcomp(t(na.omit(plot_df_major[,pca_columns])), scale.=TRUE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    scores <- data.frame(pca$x[,c("PC1","PC2")])
    scores$PC2 <- scores$PC2
    scores$TimePoint <- factor(rep(major_time_df_ui$hour, each = 4), levels = major_time_df_ui$hour)
    scores$Label <- gsub(".*_", replacement = "", x = gsub(pattern = "^imputed.log2.LFQ.intensity.",replacement = "",rownames(scores)))
    
    if (length(major_time_df_ui$hour)>1) {
      colors_major <- brewer.pal(n = 9,name = "Greys")[3:(length(major_time_df_ui$hour)+2)]
    }else{
      colors_major <- brewer.pal(n = 9,name = "Greys")[3:9]
    }
    
    # Input plot
    
    major_PCA <- ggplot(data=scores,aes(x=PC1, y=PC2, colour = TimePoint))+
      geom_hline(yintercept = 0,color = "#FFC9DE", linewidth = 1)+
      geom_vline(xintercept = 0, color = "#FFC9DE", linewidth = 1)+
      stat_ellipse(geom = "polygon", alpha = 0.15, aes(fill = TimePoint))+
      geom_point(size=3.5)+
      geom_text_repel(mapping = aes(label = Label), 
                      size = 5,
                      fontface = 'bold', 
                      color = 'black')+
      theme_minimal()+
      xlab(paste("PC1: ", round(percentVar[1] * 100,digits = 2), "%"))+
      ylab(paste("PC2: ", round(percentVar[2] * 100, digits = 2), "%"))+
      scale_y_continuous(limits = c(-150,175),breaks = seq(-150, 150, by = 50))+
      # scale_x_continuous(limits = c(-100,100),breaks = seq(-100, 100, by = 25))+
      scale_color_manual(name = "Hours\npost-infection", values = colors_major)+
      scale_fill_manual(name = "Hours\npost-infection", values = colors_major)+
      guides(colour = guide_legend(nrow = 1))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
      theme(legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 14),
            legend.key.size = unit(1, 'cm'),
            legend.position = "top")
    
    print(major_PCA)
    
  })
  
  output$major_HeatMap <- renderPlot({
    
    # Load and format the data
    
    if(length(input$major_ui_timepointID) > 0) {
      major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
    }else{
      major_time_df_ui <- major_time_df
    }
    plot_df_major <- major_FullDataTable1()
    lfq <- plot_df_major[, grep("^imputed.log2.LFQ.intensity.", names(plot_df_major))]
    colnames(lfq) <- gsub("^imputed.log2.LFQ.intensity.", "", colnames(lfq))
    pattern <- unique(gsub(pattern = "_.*", replacement = "", colnames(lfq)))
    
    for (i in 1:length(major_time_df_ui$hour)) {
      
      colnames(lfq) <- gsub(pattern = pattern[i], replacement = major_time_df_ui$hour[i], colnames(lfq))
      
    }
    
    group <- as.factor(gsub("_[1-4]$", "", colnames(lfq)))
    lfq.avg <- t(apply(lfq, 1, function(x) tapply(x, group, mean,na.rm = T)))
    v <- apply(lfq, 1, sd)
    
    corr <- apply(lfq[rev(order(v)), ], 2, function(x) {
      apply(lfq[rev(order(v)), ], 2, function(y) {
        cor(x, y, use = "na.or.complete")
      })
    })
    
    if (length(major_time_df_ui$hour)>0) {
      
      Exp <- rep("L. major", length(major_time_df_ui$hour), each = 4)
      HPI <- rep(as.character(major_time_df_ui$hour), each = 4)
      
    }else{
      
      Exp <- rep("L. major", length(as.character(hours)), each = 4)
      HPI <- rep(as.character(hours), each = 4)
      
    }
    
    df <- data.frame(HPI,Exp)
    
    rownames(df) <- colnames(corr)
    
    Exp <- "#FFC9DE"
    names(Exp) <- unique(df$Exp)
    
    HPI <- brewer.pal(n = 9,name = "Greys")[3:(length(major_time_df_ui$hour)+2)]
    names(HPI) <- unique(df$HPI)
    
    anno_colors <- list(HPI = HPI,
                        Exp = Exp)
    
    # Input plot
    
    hm_pal <- viridis(n = 15,option = "G",direction = -1)
    
    # major_HeatMap <- gplots::heatmap.2(corr,
    #                                      trace="none",
    #                                      Colv=T,
    #                                      Rowv=T,
    #                                      dendrogram="column",
    #                                      density.info="density",
    #                                      srtCol=45,
    #                                      margins=c(10, 10),
    #                                      key.xlab="Pearson's correlation coefficient",
    #                                      key.title="",
    #                                      keysize=1.5,
    #                                      col=hm_pal)
    
    major_heatmap <- pheatmap(corr,
                              col=hm_pal,
                              cluster_cols = F,
                              cluster_rows = F,
                              show_rownames = F,
                              show_colnames = F,
                              # trace = "none",
                              annotation_col = df,
                              annotation_row = df,
                              annotation_colors = anno_colors,
                              annotation_names_row = T,
                              annotation_legend = F,
                              border_color=NA,
                              fontsize = 14)
    
    dev.off()
    print(major_heatmap)
    
  })
  
  #### Dynamicity analysis output
  
  output$major_ScatterPlot <- renderPlotly({
    
    if (length(input$major_ui_timepointID) > 1) {
      
      # Settings:
      
      viridis_option <- "G"
      
      exp_color <- "#FFC9DE"
      
      # Proteins to highlight:
      
      if(length(input$major_ui_DynamicControls) == 2) {
        
        highlight <- c("LmjF.27.1870", "LmjF.29.2450", 
                       "LmjF.04.1230",
                       "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126",
                       "ENSMUSG00000004677")
      }
      
      if(length(input$major_ui_DynamicControls) == 1) {
        
        if(input$major_ui_DynamicControls == "Dynamic proteins") {
          
          highlight <- c("LmjF.27.1870", "LmjF.29.2450", 
                         "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126")
        }
        
        if(input$major_ui_DynamicControls == "Stable proteins") {
          
          highlight <- c("LmjF.04.1230",
                         "ENSMUSG00000004677")
        }
        
      }
      
      if(length(input$major_ui_DynamicControls) == 0) {
        
        highlight <- NULL
      }
      
      UI_ID <- input$major_ui_selectedIDs
      
      highlight <- c(highlight, unique(UI_ID))
      
      highlight <- ifelse(test = grepl(pattern = "ENSMUS",x = highlight),
                          yes = gsub(pattern = "\\..*", replacement = "",highlight),
                          no = highlight) 
      
      # gini function
      
      gini <- function(x) {
        x <- x - min(df_means)
        sum(sapply(1:length(x), function(i) sum(abs(x[i] - x)))) / (2 * length(x)^2 * mean(x))
      }
      
      # gini threshold
      
      GINI_THRES <- input$major_Dynamicity
      
      # LFQ threshold
      
      LFQ_THRES <- input$major_LFQ
      
      # Data
      
      df <- major_FullDataTable1()
      
      rownames(df) <- df$my_label
      
      # Retrieve the mean columns & tidy
      
      df_means <- df[,grep(pattern = "mean",x = colnames(df))]
      
      if(length(input$major_ui_timepointID) > 0) {
        major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
      }else{
        major_time_df_ui <- major_time_df
      }
      
      colnames(df_means) <- unique(as.character(major_time_df_ui$hour))
      
      # Compute gini score
      
      avg.lfq     <- apply(df_means, 1, mean)
      max.lfq     <- apply(df_means, 1, max)
      df_means$gini <- apply(df_means, 1, gini)
      
      
      x <- avg.lfq
      y <- setNames(df_means$gini, rownames(df_means))
      
      x_df <- enframe(x = x)
      colnames(x_df)[2] <- "LFQ"
      y_df <- enframe(x = y)
      colnames(y_df)[2] <- "Dynamicity"
      
      plot_df <- left_join(x = x_df, y = y_df, by = "name")
      
      highlight_df <- data.frame(name = highlight)
      
      plot_df$Highlight <- ifelse(test = grepl(pattern = "ENSMUS",x = plot_df$name),
                                  yes = gsub(pattern = "\\..*", replacement = "",plot_df$name),
                                  no = plot_df$name) %in% highlight
      
      plot_df$DynamicityThreshold <- plot_df$Dynamicity > GINI_THRES
      
      plot_df$LFQThreshold <- plot_df$LFQ > LFQ_THRES
      
      plot_df$`Threshold (D.L)` <- interaction(plot_df$DynamicityThreshold, plot_df$LFQThreshold)
      
      plot_df <- plot_df[rev(order(plot_df$`Threshold (D.L)`)),]
      
      df_to_join <- df_major[,1:2]
      
      df_to_join <- unique(df_to_join)
      
      colnames(df_to_join) <- c("gene_name", "name")
      
      plot_df <- left_join(x = plot_df, y = df_to_join, by = "name")
      
      plotly_string <- paste0(sprintf('%s',
                                      sub('(.{30}).*', '\\1...', plot_df$gene_name)))
      
      plot_df$`Gene name` <- plotly_string
      
      major_scatter_plot <- ggplot(plot_df, aes(x=LFQ, y=Dynamicity, label =`Gene name`)) +
        geom_hline(yintercept = GINI_THRES, linetype = "longdash", color = exp_color, linewidth = 1)+
        geom_vline(xintercept = LFQ_THRES, linetype = "longdash", color = exp_color, linewidth = 1)+
        geom_pointdensity(data = subset(plot_df, Highlight == FALSE), adjust = .1,alpha = .3)+
        geom_point(data = subset(plot_df, Highlight == TRUE), color ="blue", size = 3)+
        scale_color_viridis(option = viridis_option, limits = c(0,1250),breaks = seq(0, 1250, by = 250))+
        scale_y_continuous(limits = c(0,0.6),breaks = seq(0, 0.5, by = 0.25))+
        scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
        xlab("Intensity")+
        annotate('text', label=paste(sum(y>GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
        annotate('text', label=paste(sum(y<GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=min(plot_df$Dynamicity), size = 7) +
        annotate('text', label=paste(sum(y>GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
        annotate('text', label=paste(sum(y<GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=min(plot_df$Dynamicity) +0.05, size = 7) +
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      major_scatter_plot <- ggplotly(major_scatter_plot)
      
      print(major_scatter_plot)
      
    }else{
      
      major_scatter_plot <- ggplot()+
        scale_y_continuous(limits = c(0,0.6),breaks = seq(0, 0.5, by = 0.25))+
        scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
        xlab("Intensity")+
        ylab("Dynamicity")+
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      print(major_scatter_plot)
      
    }
    
  })
  
  output$major_LinePlot <- renderPlot({
    
    if (length(input$major_ui_timepointID) > 1) {
      
      # Settings:
      
      viridis_option <- "G"
      
      exp_color <- "#FFC9DE"
      
      species <- "L. major"
      
      # Proteins to highlight:
      
      if(length(input$major_ui_DynamicControls) == 2) {
        
        highlight <- c("LmjF.27.1870", "LmjF.29.2450", 
                       "LmjF.04.1230",
                       "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126",
                       "ENSMUSG00000004677")
      }
      
      if(length(input$major_ui_DynamicControls) == 1) {
        
        if(input$major_ui_DynamicControls == "Dynamic proteins") {
          
          highlight <- c("LmjF.27.1870", "LmjF.29.2450", 
                         "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126")
        }
        
        if(input$major_ui_DynamicControls == "Stable proteins") {
          
          highlight <- c("LmjF.04.1230",
                         "ENSMUSG00000004677")
        }
        
      }
      
      if(length(input$major_ui_DynamicControls) == 0) {
        
        highlight <- NULL
      }
      
      UI_ID <- input$major_ui_selectedIDs
      
      highlight <- c(highlight, unique(UI_ID))
      
      highlight <- ifelse(test = grepl(pattern = "ENSMUS",x = highlight),
                          yes = gsub(pattern = "\\..*", replacement = "",highlight),
                          no = highlight)
      
      if (length(highlight) > 0) {
        
        colors_lineplot <- c(rep(exp_color, times = length(which(!grepl(pattern = "ENSMUS", x = highlight)))),
                             rep("#b4b4b4", times = length(which(grepl(pattern = "ENSMUS", x = highlight)))))
        
        shapes_lineplot <- c(15,16,17,8,18)
        
        shapes_lineplot <- c(shapes_lineplot[1:length(which(!grepl(pattern = "ENSMUS", x = highlight)))],
                             shapes_lineplot[1:length(which(grepl(pattern = "ENSMUS", x = highlight)))])
        
        df <- major_FullDataTable1()
        
        rownames(df) <- df$my_label
        
        # Retrieve the mean columns & tidy
        
        df_means <- df[,grep(pattern = "mean",x = colnames(df))]
        
        if(length(input$major_ui_timepointID) > 0) {
          major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
        }else{
          major_time_df_ui <- major_time_df
        }
        
        colnames(df_means) <- unique(as.character(major_time_df_ui$hour))
        
        plot_df <- df_means
        
        plot_df$Highlight <- ifelse(test = grepl(pattern = "ENSMUS", x = rownames(plot_df)),
                                    yes = gsub(pattern = "\\..*", replacement = "", rownames(plot_df)),
                                    no = rownames(plot_df)) %in% highlight
        
        plot_df$name <- rownames(plot_df)
        
        df_to_join <- df_major[,1:2]
        
        df_to_join <- unique(df_to_join)
        
        colnames(df_to_join) <- c("gene_name", "name")
        
        plot_df <- left_join(x = plot_df, y = df_to_join, by = "name")
        
        plot_df <- plot_df[plot_df$Highlight == T,]
        
        plot_df <- plot_df %>% pivot_longer(-c(name, gene_name, Highlight), names_to = "Timepoint", values_to = "LFQ")
        
        plot_df$Species <- ifelse(test = grepl(pattern = "ENSMUS", x = plot_df$name),
                                  yes = "M. musculus",
                                  no = species)
        
        plot_df <- plot_df %>% arrange(Species,gene_name)
        
        plot_df$gene_name <- factor(x = plot_df$gene_name, levels = unique(plot_df$gene_name))
        
        plot_df$Timepoint <- factor(x = plot_df$Timepoint, levels = unique(as.character(major_time_df_ui$hour)))
        
        major_lineplot <- ggplot(plot_df, aes(x=Timepoint, y=LFQ, shape = gene_name, group = gene_name, color = gene_name))+
          facet_grid(cols = vars(Species))+
          geom_line(linewidth = 1)+
          geom_point(size = 4)+
          labs(shape = "Gene name")+
          guides(color="none")+
          # ylab(as.expression(bquote(bold(mean(log[2](LFQ)))))) +
          ylab("Intensity")+
          xlab("Hours post-infection") +
          scale_color_manual(labels = unique(plot_df$gene_name),
                             values = colors_lineplot)+
          scale_shape_manual(labels  = unique(plot_df$gene_name),
                             values = shapes_lineplot)+
          # scale_y_continuous(limits = c(24,32.5),breaks = seq(25, 32.5, by = 2.5))+
          # guides(color=FALSE)+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "top",
                strip.text = element_text(size = 20, face = "italic"))
        
        print(major_lineplot)
        
      }else{
        
        major_lineplot <- ggplot()+
          # scale_x_discrete(name = hours)+
          scale_y_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
          xlab("Hours post-infection")+
          ylab("Intensity")+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "none")
        
        print(major_lineplot)
        
        
      }
      
    }else{
      
      major_lineplot <- ggplot()+
        # scale_x_discrete(name = hours)+
        scale_y_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
        xlab("Hours post-infection")+
        ylab("Intensity")+
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      print(major_lineplot)
      
      showNotification(paste("Dynamicity analysis needs at least 2 timepoints!"), duration = 10,type = "error")
      
    }
    
  })
  
  #### Clustering analysis output 
  
  output$major_ClusterLinePlot <- renderPlot({
    
    # Fixed variables
    
    shortcut <- "major"
    
    pg_LFQ <- major_FullDataTable1()
    
    hl_df <- major_FullDataTable2()
    
    # User input variables
    
    xdim <- input$major_xDIM
    
    ydim <- input$major_yDIM
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    if(length(input$major_ui_timepointID) > 0) {
      major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
    }else{
      major_time_df_ui <- major_time_df
    }
    
    colnames(pg_clust) <- c("Majority.Protein.IDs",major_time_df_ui$hour)
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- major_time_df_ui$hour
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_inf <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_inf$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    # This plot will always show, wether the user selects an input ID or not
    
    mycolors <- colorRampPalette(c("#000000", "#EE2A7B"))(length(nclust))
    
    plot_list <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$TimePoint <- factor(data$TimePoint, levels = major_time_df_ui$hour)
      
      plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
        geom_line(color="lightgrey")+
        stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
        scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
        ylab("z-score(LFQ)")+
        ggtitle(paste0("Cluster ",i))+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(plot.title = element_text(size=22))+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      plot <- list(plot)
      
      plot_list <- c(plot_list, plot)
      
    }
    
    major_ClusterLinePlot <- ggarrange(plotlist = plot_list)
    
    # Input plot
    
    # User dependant plot. The IDs the user select will show up in this plot, when present.
    
    if (nrow(hl_df) > 0) {
      
      HID <- unique(hl_df$Majority.protein.IDs)
      
      plot_list <- c()
      
      for (i in 1:length(nclust)) {
        
        data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
        
        data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
        
        data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
        
        dataHL <- data[data$ID %in% HID,]
        
        if (nrow(dataHL) > 0) {
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          for (i in 1:length(HID)) {
            
            plot <- plot +
              geom_line(data=data[data$ID == HID[i],], aes(x=TimePoint, y=Value, group=ID,color = ID))+
              scale_color_brewer(name = "", palette = "Dark2")
            
            
          }
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }else{
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }
        
      }
      
      major_ClusterLinePlot <- ggarrange(plotlist = plot_list)
      
    }
    
    print(major_ClusterLinePlot)
    
  })
  
  output$major_ClusterBarPlot <- renderPlotly({
    
    # Fixed variables
    
    shortcut <- "major"
    
    pg_LFQ <- major_FullDataTable1()
    
    hl_df <- major_FullDataTable2()
    
    # User input variables
    
    xdim <- input$major_xDIM
    
    ydim <- input$major_yDIM
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    if(length(input$major_ui_timepointID) > 0) {
      major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
    }else{
      major_time_df_ui <- major_time_df
    }
    
    colnames(pg_clust) <- c("Majority.Protein.IDs", major_time_df_ui$hour)
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- major_time_df_ui$hour
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_inf <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_inf$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    mycolors <- colorRampPalette(c("#000000", "#EE2A7B"))(length(nclust))
    
    plot_data <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$Species <- "Unknown"
      
      data$Species <-
        ifelse(grepl("ENSMUS.*", data$ID),
               sub(".*", "M. musculus", data$Species),
               sub(".*", "Unknown", data$Species))
      
      data$Species <-
        ifelse(grepl("LmjF.*", data$ID),
               sub(".*", "L. major", data$Species),
               sub("", "", data$Species))
      
      data$TimePoint <- factor(data$TimePoint, levels = major_time_df_ui$hour)
      
      data <- data[data$TimePoint == major_time_df_ui$hour[1],]
      
      data <- data %>%
        group_by(Species) %>%
        summarize(Count = n())%>% 
        mutate(perc = Count/sum(Count))
      
      data$Cluster <- as.character(i)
      
      plot_data <- rbind(plot_data, data)
      
    }
    
    major_ClusterBarPlot <- ggplot(plot_data, aes(fill=Species, y = Count , x= Cluster)) + 
      geom_bar(position="stack", stat="identity")+
      scale_fill_manual(values = c("L. major" = "#EE2A7B",
                                   "M. musculus" = "#b4b4b4"))+
      theme_minimal()+
      ylab("Protein IDs")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(legend.title = element_text(size = 16))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    major_ClusterBarPlot <- ggplotly(major_ClusterBarPlot)
    
    print(major_ClusterBarPlot)
    
  })
  
  #### Pairwise analysis  time course output 
  
  output$major_Boxplot <- renderPlot({
    
    if(length(input$major_ui_timepointID) > 0) {
      major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
    }else{
      major_time_df_ui <- major_time_df
    }
    
    if (length(major_time_df_ui$hour)>1) {
      colors_major <- brewer.pal(n = 9,name = "Greys")[3:(length(major_time_df_ui$hour)+2)]
    }else{
      colors_major <- brewer.pal(n = 9,name = "Greys")[3:9]
    }
    
    colors_major_replica <- colorRampPalette(c("#FFFFFF", "#FFC9DE", "#000000"))(6)[2:5]
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_major <- major_MainDataTable1() 
    
    ## Selection data
    
    jitter_data_major <- major_MainDataTable2() 
    
    # Dummy plot
    
    major_Boxplot <-ggboxplot(data = plot_df_major, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
      scale_fill_manual(values = colors_major)+
      scale_color_brewer(palette="Dark2")+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))+
      theme(legend.title = element_text(size = 12))+
      theme(legend.text = element_text(size = 10))+
      theme(legend.key.size = unit(1, 'cm'))
    
    # Input plot
    
    hl_df <- major_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    if (nrow(jitter_data_major) > 0) {
      
      boxplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_major <- jitter_data_major[jitter_data_major$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        major_Boxplot <- ggboxplot(data = plot_df_major, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
          geom_jitter(data = jitter_data_major[jitter_data_major$Majority.protein.IDs == HID[i],],aes(x=TimePoint, y=Value, color = Replica),width = 0.20)+
          geom_line(data = line_data_major, aes(x=TimePoint, y=Value,group = 1))+
          scale_fill_manual(values = colors_major)+
          scale_color_manual(values = colors_major_replica)+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        major_Boxplot <- list(major_Boxplot)
        
        boxplot_list <- c(boxplot_list, major_Boxplot)
        
      }
      
      major_Boxplot <- ggarrange(plotlist = boxplot_list)
      
    }
    
    print(major_Boxplot)
  })
  
  output$major_Dotplot <- renderPlot({
    
    colors_major_replica <- colorRampPalette(c("#FFFFFF", "#FFC9DE", "#000000"))(6)[2:5]
    
    ## Retrieve highlight data
    
    hl_df <- major_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_major <- major_MainDataTable1()
    
    ## Selection data
    
    jitter_data_major <- major_MainDataTable2()
    
    test_data_major <- major_TestDataTable2()
    
    # Dummy plot
    
    major_Dotplot <- ggplot(data = jitter_data_major, aes(TimePoint, Value))+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    # Input plot
    
    if (nrow(jitter_data_major) > 0) {
      
      y_pos_dot <- max(jitter_data_major$Value)+0.2
      
      dotplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_major <- jitter_data_major[jitter_data_major$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        # Format test data
        
        plot_test_data_major <- test_data_major[test_data_major$Majority.protein.IDs == HID[i],]
        plot_test_data_major <- plot_test_data_major[,c(4,5,3)]
        plot_test_data_major$p.adj <- format.pval(plot_test_data_major$p.adj, digits = 3)
        
        major_Dotplot <- ggdotplot(data = jitter_data_major[jitter_data_major$Majority.protein.IDs == HID[i],], x="TimePoint", y="Value", fill = "Replica",
                                      position = position_jitter(0.2))+
          geom_line(data = line_data_major, aes(x=TimePoint, y=Value,group = 1),alpha = 0.8)+
          stat_pvalue_manual(plot_test_data_major, y.position = y_pos_dot, step.increase = 0.1,label = "p.adj")+
          scale_fill_manual(values= colors_major_replica)+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),
                axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        major_Dotplot <- list(major_Dotplot)
        
        dotplot_list <- c(dotplot_list, major_Dotplot)
        
      }
      
      major_Dotplot <- ggarrange(plotlist = dotplot_list)
      
    }
    
    print(major_Dotplot)
  })
  
  #### Pairwise analysis  reference output 
  
  output$major_Volcanoplot <- renderPlotly({
    
    #Input variables
    
    filter <- input$major_ui_datasetID
    exp_a <- input$major_TP
    exp_b <- input$major_referenceTP
    ID <- input$major_ui_selectedIDs
    enrich_s0 <- log2(2^as.numeric(input$major_LFC)) 
    enrich_pvalue <- as.numeric(input$major_Pvalue)
    
    # Fixed variables
    
    difference <- paste('difference', exp_a, exp_b, sep='_')
    pvalue <- paste('log10.pvalue', exp_a, exp_b, sep='_')
    val_count_a <- paste0('value_count_',exp_a) 
    val_count_b <- paste0('value_count_',exp_b)
    enriched_col <- paste('enriched', exp_a, exp_b, sep='_')
    
    min_quant_events <- 0 
    replicate_count <- 4
    enrich_c <- .05 
    volcano_threshold_linetype <- 2
    volcano_threshold_linesize <- 0.5
    difference_volcano_column <- 'mean2_naimp_meas_'
    volcano_column_type <- sub('([^[:digit:]]+).*', '\\1', difference_volcano_column)
    
    # Dummy plot
    
    ## Filter the data frame for data set user input.
    
    if (filter == paste0(unique(df_major$Species),collapse = " & ")) {
      pg_quant <- full_major
    }else{
      pg_quant <- full_major %>% filter(Species == filter)
    }
    
    ## Change the column pattern to user friendly reading (i. e., 0h instead of inft000).
    
    if(length(input$major_ui_timepointID) > 0) {
      major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
    }else{
      major_time_df_ui <- major_time_df
    }
    
    for (i in 1:length(major_time_df_ui$hour)) {
      
      colnames(pg_quant) <- gsub(pattern = major_time_df_ui$time_point[i], replacement = major_time_df_ui$hour[i], colnames(pg_quant))
      
    }
    
    ## Dummy plot to show until user selection.
    
    dummy_df <- data.frame(x = c(-5,5),
                           y = c(0,10))
    
    volcanoplot <-
      ggplot(dummy_df, aes(x=x, y=y))+
      theme_minimal()+
      ylab("-log<sub>10</sub>(p-value)") + 
      xlab("log<sub>2</sub>(Fold change)") +
      theme(legend.text = element_text(size = 12))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(plot.title = element_text(size=20))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
    
    major_Volcanoplot <-
      plotly::ggplotly(volcanoplot) %>%
      plotly::layout(showlegend = FALSE)
    
    # Input plot
    
    ## Upon proper user selection (if) produce the volcano plot.
    
    if (difference %in% colnames(pg_quant)) {
      
      ## Highlight user input, if any.
      
      pg_quant$highlight <- FALSE
      
      if (length(ID) != 0) {
        
        ids_searchstring <- paste(ID, collapse='|')
        
        matching_rows <- grep(ids_searchstring, as.character(pg_quant$Protein.IDs))
        
        pg_quant$highlight[matching_rows] <- TRUE
        
      }
      
      
      ## Filter count values. Standard set to 0, so they can check the dot all along the time course.
      ## In Mario's standard script is usually set to 2. But in this case, since it shown in the 
      ## box plot and in the dot plot, would be confusing not to show it here as well.
      ## Can be adjusted on the fixed variables section of this function.
      
      filter_min_quant_events <-
        pg_quant[val_count_a] >= min_quant_events |
        pg_quant[val_count_b] >= min_quant_events
      
      pg_quant <- pg_quant[filter_min_quant_events,]
      
      ## Select and highlight the enriched dots.
      
      my_enriched <- Enriched(pg_quant[c(difference, pvalue)],
                              c=enrich_c,
                              s0=enrich_s0,
                              pvalue=enrich_pvalue,
                              linetype=volcano_threshold_linetype,
                              size=volcano_threshold_linesize)
      
      pg_quant$enriched <- as.factor(my_enriched$enriched)
      
      
      my_points_nohighlight <- interaction(pg_quant[!pg_quant$highlight, c('enriched', 'highlight')])
      my_points_highlight <- interaction(pg_quant[pg_quant$highlight, c('enriched', 'highlight')])
      
      ## Generate a user friendly plotly string for the hoover.
      
      plotly_string <- paste0(sprintf('log2(foldchange): %.2f, -log10(pvalue): %.2f<br />Protein.IDs: %s<br />my_label: %s',
                                      pg_quant[[difference]], pg_quant[[pvalue]],
                                      sub('(.{30}).*', '\\1...', pg_quant$Protein.IDs),
                                      pg_quant$my_label),
                              switch('Gene.names' %in% names(pg_quant),
                                     sprintf('<br />Gene.names: %s', pg_quant$Gene.names),
                                     NULL),
                              switch('Protein.names' %in% names(pg_quant),
                                     sprintf('<br />Protein.names: %s', pg_quant$Protein.names),
                                     NULL),
                              switch('description' %in% names(pg_quant),
                                     sprintf('<br />description: %s', pg_quant$description),
                                     NULL))
      
      pg_quant$my_text <- plotly_string
      
      ## Subset the data for the volcano plot.
      
      pg_quant_subset <-
        pg_quant %>%
        rowwise() %>%
        mutate(imputation=ifelse(max(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count &
                                   min(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count,
                                 'measured',
                                 ifelse(min(range(.data[[val_count_a]], .data[[val_count_b]])) > 0,
                                        'some imputed',
                                        'ref imputed'))) %>%
        dplyr::select(Protein.IDs, my_label, matches(difference), matches(pvalue),
               my_text, highlight, enriched, imputation)
      
      ## Plot
      
      ## Select maximum x and y values of your volcano plot.
      
      max_x_value <- max(abs(pg_quant[grep('^difference', names(pg_quant))]),
                         na.rm=TRUE)
      max_y_value <- max(c(unlist(pg_quant[grep('^log10.pvalue', names(pg_quant))]),
                           10^enrich_pvalue + 0.5))
      
      volcanoplot <-
        ggplot(pg_quant_subset, aes_string(x=difference, y=pvalue, text='my_text')) +
        geom_hline(yintercept=0, color='grey80', na.rm=TRUE) +
        geom_vline(xintercept=0, color='grey80', na.rm=TRUE) +
        plot(my_enriched) +
        geom_point(data=subset(pg_quant_subset, !highlight),
                   aes(alpha=my_points_nohighlight,
                       color=my_points_nohighlight,
                       fill=my_points_nohighlight,
                       size=my_points_nohighlight),
                   na.rm=TRUE) +
        geom_point(data=subset(pg_quant_subset, highlight),
                   aes(alpha=my_points_highlight,
                       color=my_points_highlight,
                       fill=my_points_highlight,
                       size=my_points_highlight),
                   na.rm=TRUE) +
        scale_alpha_manual(values=c(`TRUE.FALSE`=.7, `FALSE.FALSE`=.4, `TRUE.TRUE`=.7, `FALSE.TRUE`=.7),
                           guide="none") +
        scale_color_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                           guide="none") +
        scale_fill_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                          guide="none") +
        scale_size_manual(values=c(`TRUE.FALSE`=1, `FALSE.FALSE`=.5, `TRUE.TRUE`=1.2, `FALSE.TRUE`=1),
                          guide="none")+
        coord_cartesian(xlim=c(-max_x_value,max_x_value), ylim=c(-.8,max_y_value)) +
        ylab("-log<sub>10</sub>(p-value)") + 
        xlab("log<sub>2</sub>(Fold change)") +
        annotate('text', label=exp_a, x=max_x_value / 2, y=-.4,size = 6) +
        annotate('text', label=exp_b, x=-max_x_value / 2, y=-.4,size = 6)+
        theme_minimal()+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      major_Volcanoplot <-
        plotly::ggplotly(volcanoplot,tooltip='text') %>%
        plotly::layout(showlegend = FALSE)
      
    }
    
    print(major_Volcanoplot)
    
  })
  
  #### L. mexicana ####
  
  #### Main table output  
  
  output$mexicana_MainDataTable <- renderDataTable({
    mexicana_MainDataTable3() %>%
      datatable(options = list(pageLength = 7, autoWidth = TRUE), 
                rownames = F,escape = FALSE) %>%
      formatRound("mean(LFQ)", digits = 2)
  })
  
  mexicana_MainDataTable_proxy <- dataTableProxy('mexicana_MainDataTable')
  
  #### Exploratory analysis output 
  
  output$mexicana_PCA <- renderPlot({
    
    # Load and format the data
    
    if(length(input$mexicana_ui_timepointID) > 0) {
      mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
    }else{
      mexicana_time_df_ui <- mexicana_time_df
    }
    plot_df_mexicana <- mexicana_FullDataTable1()
    pca_columns <- grep(pattern = "^imputed.log2.LFQ.intensity.", colnames(plot_df_mexicana))
    pca <- prcomp(t(na.omit(plot_df_mexicana[,pca_columns])), scale.=TRUE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    scores <- data.frame(pca$x[,c("PC1","PC2")])
    scores$PC2 <- scores$PC2 
    scores$TimePoint <- factor(rep(mexicana_time_df_ui$hour, each = 4), levels = mexicana_time_df_ui$hour)
    scores$Label <- gsub(".*_", replacement = "", x = gsub(pattern = "^imputed.log2.LFQ.intensity.",replacement = "",rownames(scores)))
    
    if (length(mexicana_time_df_ui$hour)>1) {
      colors_mexicana <- brewer.pal(n = 9,name = "Greys")[3:(length(mexicana_time_df_ui$hour)+2)]
    }else{
      colors_mexicana <- brewer.pal(n = 9,name = "Greys")[3:9]
    }
    
    # Input plot
    
    mexicana_PCA <- ggplot(data=scores,aes(x=PC1, y=PC2, colour = TimePoint))+
      geom_hline(yintercept = 0,color = "#6EB5FF", linewidth = 1)+
      geom_vline(xintercept = 0, color = "#6EB5FF", linewidth = 1)+
      stat_ellipse(geom = "polygon", alpha = 0.15, aes(fill = TimePoint))+
      geom_point(size=3.5)+
      geom_text_repel(mapping = aes(label = Label), 
                      size = 5,
                      fontface = 'bold', 
                      color = 'black')+
      theme_minimal()+
      xlab(paste("PC1: ", round(percentVar[1] * 100,digits = 2), "%"))+
      ylab(paste("PC2: ", round(percentVar[2] * 100, digits = 2), "%"))+
      scale_y_continuous(limits = c(-150,175),breaks = seq(-150, 150, by = 50))+
      # scale_x_continuous(limits = c(-100,100),breaks = seq(-100, 100, by = 25))+
      scale_color_manual(name = "Hours\npost-infection", values = colors_mexicana)+
      scale_fill_manual(name = "Hours\npost-infection", values = colors_mexicana)+
      guides(colour = guide_legend(nrow = 1))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
      theme(legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 14),
            legend.key.size = unit(1, 'cm'),
            legend.position = "top")
    
    print(mexicana_PCA)
    
  })
  
  output$mexicana_HeatMap <- renderPlot({
    
    # Load and format the data
    
    if(length(input$mexicana_ui_timepointID) > 0) {
      mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
    }else{
      mexicana_time_df_ui <- mexicana_time_df
    }
    plot_df_mexicana <- mexicana_FullDataTable1()
    lfq <- plot_df_mexicana[, grep("^imputed.log2.LFQ.intensity.", names(plot_df_mexicana))]
    colnames(lfq) <- gsub("^imputed.log2.LFQ.intensity.", "", colnames(lfq))
    pattern <- unique(gsub(pattern = "_.*", replacement = "", colnames(lfq)))
    
    for (i in 1:length(mexicana_time_df_ui$hour)) {
      
      colnames(lfq) <- gsub(pattern = pattern[i], replacement = mexicana_time_df_ui$hour[i], colnames(lfq))
      
    }
    
    group <- as.factor(gsub("_[1-4]$", "", colnames(lfq)))
    lfq.avg <- t(apply(lfq, 1, function(x) tapply(x, group, mean,na.rm = T)))
    v <- apply(lfq, 1, sd)
    
    corr <- apply(lfq[rev(order(v)), ], 2, function(x) {
      apply(lfq[rev(order(v)), ], 2, function(y) {
        cor(x, y, use = "na.or.complete")
      })
    })
    
    if (length(mexicana_time_df_ui$hour)>0) {
      
      Exp <- rep("L. mexicana", length(mexicana_time_df_ui$hour), each = 4)
      HPI <- rep(as.character(mexicana_time_df_ui$hour), each = 4)
      
    }else{
      
      Exp <- rep("L. mexicana", length(as.character(hours)), each = 4)
      HPI <- rep(as.character(hours), each = 4)
      
    }
    
    df <- data.frame(HPI,Exp)
    
    rownames(df) <- colnames(corr)
    
    Exp <- "#6EB5FF"
    names(Exp) <- unique(df$Exp)
    
    HPI <- brewer.pal(n = 9,name = "Greys")[3:(length(mexicana_time_df_ui$hour)+2)]
    names(HPI) <- unique(df$HPI)
    
    anno_colors <- list(HPI = HPI,
                        Exp = Exp)
    
    # Input plot
    
    hm_pal <- viridis(n = 15,option = "G",direction = -1)
    
    # mexicana_HeatMap <- gplots::heatmap.2(corr,
    #                                      trace="none",
    #                                      Colv=T,
    #                                      Rowv=T,
    #                                      dendrogram="column",
    #                                      density.info="density",
    #                                      srtCol=45,
    #                                      margins=c(10, 10),
    #                                      key.xlab="Pearson's correlation coefficient",
    #                                      key.title="",
    #                                      keysize=1.5,
    #                                      col=hm_pal)
    
    mexicana_heatmap <- pheatmap(corr,
                                 col=hm_pal,
                                 cluster_cols = F,
                                 cluster_rows = F,
                                 show_rownames = F,
                                 show_colnames = F,
                                 # trace = "none",
                                 annotation_col = df,
                                 annotation_row = df,
                                 annotation_colors = anno_colors,
                                 annotation_names_row = T,
                                 annotation_legend = F,
                                 border_color=NA,
                                 fontsize = 14)
    
    dev.off()
    print(mexicana_heatmap)
    
  })
  
  #### Dynamicity analysis output
  
  output$mexicana_ScatterPlot <- renderPlotly({
    
    if (length(input$mexicana_ui_timepointID) > 1) {
      
      # Settings:
      
      viridis_option <- "G"
      
      exp_color <- "#6EB5FF"
      
      # Proteins to highlight:
      
      if(length(input$mexicana_ui_DynamicControls) == 2) {
        
        highlight <- c("LmxM.34.1480", "LmxM.31.2260", "LmxM.27.1870", "LmxM.15.1160",
                       "LmxM.04.1230",
                       "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126",
                       "ENSMUSG00000004677")
      }
      
      if(length(input$mexicana_ui_DynamicControls) == 1) {
        
        if(input$mexicana_ui_DynamicControls == "Dynamic proteins") {
          
          highlight <- c("LmxM.34.1480", "LmxM.31.2260", "LmxM.27.1870", "LmxM.15.1160",
                         "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126")
        }
        
        if(input$mexicana_ui_DynamicControls == "Stable proteins") {
          
          highlight <- c("LmxM.04.1230",
                         "ENSMUSG00000004677")
        }
        
      }
      
      if(length(input$mexicana_ui_DynamicControls) == 0) {
        
        highlight <- NULL
      }
      
      UI_ID <- input$mexicana_ui_selectedIDs
      
      highlight <- c(highlight, unique(UI_ID))
      
      highlight <- ifelse(test = grepl(pattern = "ENSMUS",x = highlight),
                          yes = gsub(pattern = "\\..*", replacement = "",highlight),
                          no = highlight) 
      
      # gini function
      
      gini <- function(x) {
        x <- x - min(df_means)
        sum(sapply(1:length(x), function(i) sum(abs(x[i] - x)))) / (2 * length(x)^2 * mean(x))
      }
      
      # gini threshold
      
      GINI_THRES <- input$mexicana_Dynamicity
      
      # LFQ threshold
      
      LFQ_THRES <- input$mexicana_LFQ
      
      # Data
      
      df <- mexicana_FullDataTable1()
      
      rownames(df) <- df$my_label
      
      # Retrieve the mean columns & tidy
      
      df_means <- df[,grep(pattern = "mean",x = colnames(df))]
      
      if(length(input$mexicana_ui_timepointID) > 0) {
        mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
      }else{
        mexicana_time_df_ui <- mexicana_time_df
      }
      
      colnames(df_means) <- unique(as.character(mexicana_time_df_ui$hour))
      
      # Compute gini score
      
      avg.lfq     <- apply(df_means, 1, mean)
      max.lfq     <- apply(df_means, 1, max)
      df_means$gini <- apply(df_means, 1, gini)
      
      
      x <- avg.lfq
      y <- setNames(df_means$gini, rownames(df_means))
      
      x_df <- enframe(x = x)
      colnames(x_df)[2] <- "LFQ"
      y_df <- enframe(x = y)
      colnames(y_df)[2] <- "Dynamicity"
      
      plot_df <- left_join(x = x_df, y = y_df, by = "name")
      
      highlight_df <- data.frame(name = highlight)
      
      plot_df$Highlight <- ifelse(test = grepl(pattern = "ENSMUS",x = plot_df$name),
                                  yes = gsub(pattern = "\\..*", replacement = "",plot_df$name),
                                  no = plot_df$name) %in% highlight
      
      plot_df$DynamicityThreshold <- plot_df$Dynamicity > GINI_THRES
      
      plot_df$LFQThreshold <- plot_df$LFQ > LFQ_THRES
      
      plot_df$`Threshold (D.L)` <- interaction(plot_df$DynamicityThreshold, plot_df$LFQThreshold)
      
      plot_df <- plot_df[rev(order(plot_df$`Threshold (D.L)`)),]
      
      df_to_join <- df_mexicana[,1:2]
      
      df_to_join <- unique(df_to_join)
      
      colnames(df_to_join) <- c("gene_name", "name")
      
      plot_df <- left_join(x = plot_df, y = df_to_join, by = "name")
      
      plotly_string <- paste0(sprintf('%s',
                                      sub('(.{30}).*', '\\1...', plot_df$gene_name)))
      
      plot_df$`Gene name` <- plotly_string
      
      mexicana_scatter_plot <- ggplot(plot_df, aes(x=LFQ, y=Dynamicity, label =`Gene name`)) +
        geom_hline(yintercept = GINI_THRES, linetype = "longdash", color = exp_color, linewidth = 1)+
        geom_vline(xintercept = LFQ_THRES, linetype = "longdash", color = exp_color, linewidth = 1)+
        geom_pointdensity(data = subset(plot_df, Highlight == FALSE), adjust = .1,alpha = .3)+
        geom_point(data = subset(plot_df, Highlight == TRUE), color ="blue", size = 3)+
        scale_color_viridis(option = viridis_option, limits = c(0,1250),breaks = seq(0, 1250, by = 250))+
        scale_y_continuous(limits = c(0,0.6),breaks = seq(0, 0.5, by = 0.25))+
        scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
        xlab("Intensity")+
        annotate('text', label=paste(sum(y>GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
        annotate('text', label=paste(sum(y<GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=min(plot_df$Dynamicity), size = 7) +
        annotate('text', label=paste(sum(y>GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
        annotate('text', label=paste(sum(y<GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=min(plot_df$Dynamicity) +0.05, size = 7) +
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      mexicana_scatter_plot <- ggplotly(mexicana_scatter_plot)
      
      print(mexicana_scatter_plot)
      
    }else{
      
      mexicana_scatter_plot <- ggplot()+
        scale_y_continuous(limits = c(0,0.6),breaks = seq(0, 0.5, by = 0.25))+
        scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
        xlab("Intensity")+
        ylab("Dynamicity")+
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      print(mexicana_scatter_plot)
      
    }
    
  })
  
  output$mexicana_LinePlot <- renderPlot({
    
    if (length(input$mexicana_ui_timepointID) > 1) {
      
      # Settings:
      
      viridis_option <- "G"
      
      exp_color <- "#6EB5FF"
      
      species <- "L. mexicana"
      
      # Proteins to highlight:
      
      if(length(input$mexicana_ui_DynamicControls) == 2) {
        
        highlight <- c("LmxM.34.1480", "LmxM.31.2260", "LmxM.27.1870", "LmxM.15.1160",
                       "LmxM.04.1230",
                       "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126",
                       "ENSMUSG00000004677")
      }
      
      if(length(input$mexicana_ui_DynamicControls) == 1) {
        
        if(input$mexicana_ui_DynamicControls == "Dynamic proteins") {
          
          highlight <- c("LmxM.34.1480", "LmxM.31.2260", "LmxM.27.1870", "LmxM.15.1160",
                         "ENSMUSG00000027995", "ENSMUSG00000058818","ENSMUSG00000034459", "ENSMUSG00000022126")
        }
        
        if(input$mexicana_ui_DynamicControls == "Stable proteins") {
          
          highlight <- c("LmxM.04.1230",
                         "ENSMUSG00000004677")
        }
        
      }
      
      if(length(input$mexicana_ui_DynamicControls) == 0) {
        
        highlight <- NULL
      }
      
      UI_ID <- input$mexicana_ui_selectedIDs
      
      highlight <- c(highlight, unique(UI_ID))
      
      highlight <- ifelse(test = grepl(pattern = "ENSMUS",x = highlight),
                          yes = gsub(pattern = "\\..*", replacement = "",highlight),
                          no = highlight)
      
      if (length(highlight) > 0) {
        
        colors_lineplot <- c(rep(exp_color, times = length(which(!grepl(pattern = "ENSMUS", x = highlight)))),
                             rep("#b4b4b4", times = length(which(grepl(pattern = "ENSMUS", x = highlight)))))
        
        shapes_lineplot <- c(15,16,17,8,18)
        
        shapes_lineplot <- c(shapes_lineplot[1:length(which(!grepl(pattern = "ENSMUS", x = highlight)))],
                             shapes_lineplot[1:length(which(grepl(pattern = "ENSMUS", x = highlight)))])
        
        df <- mexicana_FullDataTable1()
        
        rownames(df) <- df$my_label
        
        # Retrieve the mean columns & tidy
        
        df_means <- df[,grep(pattern = "mean",x = colnames(df))]
        
        if(length(input$mexicana_ui_timepointID) > 0) {
          mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
        }else{
          mexicana_time_df_ui <- mexicana_time_df
        }
        
        colnames(df_means) <- unique(as.character(mexicana_time_df_ui$hour))
        
        plot_df <- df_means
        
        plot_df$Highlight <- ifelse(test = grepl(pattern = "ENSMUS", x = rownames(plot_df)),
                                    yes = gsub(pattern = "\\..*", replacement = "", rownames(plot_df)),
                                    no = rownames(plot_df)) %in% highlight
        
        plot_df$name <- rownames(plot_df)
        
        df_to_join <- df_mexicana[,1:2]
        
        df_to_join <- unique(df_to_join)
        
        colnames(df_to_join) <- c("gene_name", "name")
        
        plot_df <- left_join(x = plot_df, y = df_to_join, by = "name")
        
        plot_df <- plot_df[plot_df$Highlight == T,]
        
        plot_df <- plot_df %>% pivot_longer(-c(name, gene_name, Highlight), names_to = "Timepoint", values_to = "LFQ")
        
        plot_df$Species <- ifelse(test = grepl(pattern = "ENSMUS", x = plot_df$name),
                                  yes = "M. musculus",
                                  no = species)
        
        plot_df <- plot_df %>% arrange(Species,gene_name)
        
        plot_df$gene_name <- factor(x = plot_df$gene_name, levels = unique(plot_df$gene_name))
        
        plot_df$Timepoint <- factor(x = plot_df$Timepoint, levels = unique(as.character(mexicana_time_df_ui$hour)))
        
        mexicana_lineplot <- ggplot(plot_df, aes(x=Timepoint, y=LFQ, shape = gene_name, group = gene_name, color = gene_name))+
          facet_grid(cols = vars(Species))+
          geom_line(linewidth = 1)+
          geom_point(size = 4)+
          labs(shape = "Gene name")+
          guides(color="none")+
          # ylab(as.expression(bquote(bold(mean(log[2](LFQ)))))) +
          ylab("Intensity")+
          xlab("Hours post-infection") +
          scale_color_manual(labels = unique(plot_df$gene_name),
                             values = colors_lineplot)+
          scale_shape_manual(labels  = unique(plot_df$gene_name),
                             values = shapes_lineplot)+
          # scale_y_continuous(limits = c(24,32.5),breaks = seq(25, 32.5, by = 2.5))+
          # guides(color=FALSE)+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "top",
                strip.text = element_text(size = 20, face = "italic"))
        
        print(mexicana_lineplot)
        
      }else{
        
        mexicana_lineplot <- ggplot()+
          # scale_x_discrete(name = hours)+
          scale_y_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
          xlab("Hours post-infection")+
          ylab("Intensity")+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "none")
        
        print(mexicana_lineplot)
        
        
      }
      
    }else{
      
      mexicana_lineplot <- ggplot()+
        # scale_x_discrete(name = hours)+
        scale_y_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
        xlab("Hours post-infection")+
        ylab("Intensity")+
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      print(mexicana_lineplot)
      
      showNotification(paste("Dynamicity analysis needs at least 2 timepoints!"), duration = 10,type = "error")
      
    }
    
  })
  
  #### Clustering analysis output 
  
  output$mexicana_ClusterLinePlot <- renderPlot({
    
    # Fixed variables
    
    shortcut <- "mex"
    
    pg_LFQ <- mexicana_FullDataTable1()
    
    hl_df <- mexicana_FullDataTable2()
    
    # User input variables
    
    xdim <- input$mexicana_xDIM
    
    ydim <- input$mexicana_yDIM
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    if(length(input$mexicana_ui_timepointID) > 0) {
      mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
    }else{
      mexicana_time_df_ui <- mexicana_time_df
    }
    
    colnames(pg_clust) <- c("Majority.Protein.IDs",mexicana_time_df_ui$hour)
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- mexicana_time_df_ui$hour
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_inf <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_inf$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    # This plot will always show, wether the user selects an input ID or not
    
    mycolors <- colorRampPalette(c("#000000", "#2B3990"))(length(nclust))
    
    plot_list <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$TimePoint <- factor(data$TimePoint, levels = mexicana_time_df_ui$hour)
      
      plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
        geom_line(color="lightgrey")+
        stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
        scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
        ylab("z-score(LFQ)")+
        ggtitle(paste0("Cluster ",i))+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(plot.title = element_text(size=22))+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      plot <- list(plot)
      
      plot_list <- c(plot_list, plot)
      
    }
    
    mexicana_ClusterLinePlot <- ggarrange(plotlist = plot_list)
    
    # Input plot
    
    # User dependant plot. The IDs the user select will show up in this plot, when present.
    
    if (nrow(hl_df) > 0) {
      
      HID <- unique(hl_df$Majority.protein.IDs)
      
      plot_list <- c()
      
      for (i in 1:length(nclust)) {
        
        data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
        
        data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
        
        data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
        
        dataHL <- data[data$ID %in% HID,]
        
        if (nrow(dataHL) > 0) {
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          for (i in 1:length(HID)) {
            
            plot <- plot +
              geom_line(data=data[data$ID == HID[i],], aes(x=TimePoint, y=Value, group=ID,color = ID))+
              scale_color_brewer(name = "", palette = "Dark2")
            
            
          }
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }else{
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }
        
      }
      
      mexicana_ClusterLinePlot <- ggarrange(plotlist = plot_list)
      
    }
    
    print(mexicana_ClusterLinePlot)
    
  })
  
  output$mexicana_ClusterBarPlot <- renderPlotly({
    
    # Fixed variables
    
    shortcut <- "mex"
    
    pg_LFQ <- mexicana_FullDataTable1()
    
    hl_df <- mexicana_FullDataTable2()
    
    # User input variables
    
    xdim <- input$mexicana_xDIM
    
    ydim <- input$mexicana_yDIM
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    if(length(input$mexicana_ui_timepointID) > 0) {
      mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
    }else{
      mexicana_time_df_ui <- mexicana_time_df
    }
    
    colnames(pg_clust) <- c("Majority.Protein.IDs", mexicana_time_df_ui$hour)
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- mexicana_time_df_ui$hour
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_inf <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_inf$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    mycolors <- colorRampPalette(c("#000000", "#2B3990"))(length(nclust))
    
    plot_data <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$Species <- "Unknown"
      
      data$Species <-
        ifelse(grepl("ENSMUS.*", data$ID),
               sub(".*", "M. musculus", data$Species),
               sub(".*", "Unknown", data$Species))
      
      data$Species <-
        ifelse(grepl("LmxM.*", data$ID),
               sub(".*", "L. mexicana", data$Species),
               sub("", "", data$Species))
      
      data$TimePoint <- factor(data$TimePoint, levels = mexicana_time_df_ui$hour)
      
      data <- data[data$TimePoint == mexicana_time_df_ui$hour[1],]
      
      data <- data %>%
        group_by(Species) %>%
        summarize(Count = n())%>% 
        mutate(perc = Count/sum(Count))
      
      data$Cluster <- as.character(i)
      
      plot_data <- rbind(plot_data, data)
      
    }
    
    mexicana_ClusterBarPlot <- ggplot(plot_data, aes(fill=Species, y = Count , x= Cluster)) + 
      geom_bar(position="stack", stat="identity")+
      scale_fill_manual(values = c("L. mexicana" = "#2B3990",
                                   "M. musculus" = "#b4b4b4"))+
      theme_minimal()+
      ylab("Protein IDs")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(legend.title = element_text(size = 16))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    mexicana_ClusterBarPlot <- ggplotly(mexicana_ClusterBarPlot)
    
    print(mexicana_ClusterBarPlot)
    
  })
  
  #### Pairwise analysis  time course output 
  
  output$mexicana_Boxplot <- renderPlot({
    
    if(length(input$mexicana_ui_timepointID) > 0) {
      mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
    }else{
      mexicana_time_df_ui <- mexicana_time_df
    }
    
    if (length(mexicana_time_df_ui$hour)>1) {
      colors_mexicana <- brewer.pal(n = 9,name = "Greys")[3:(length(mexicana_time_df_ui$hour)+2)]
    }else{
      colors_mexicana <- brewer.pal(n = 9,name = "Greys")[3:9]
    }
    
    colors_mexicana_replica <- colorRampPalette(c("#FFFFFF", "#6EB5FF", "#000000"))(6)[2:5]
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_mexicana <- mexicana_MainDataTable1() 
    
    ## Selection data
    
    jitter_data_mexicana <- mexicana_MainDataTable2() 
    
    # Dummy plot
    
    mexicana_Boxplot <-ggboxplot(data = plot_df_mexicana, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
      scale_fill_manual(values = colors_mexicana)+
      scale_color_brewer(palette="Dark2")+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))+
      theme(legend.title = element_text(size = 12))+
      theme(legend.text = element_text(size = 10))+
      theme(legend.key.size = unit(1, 'cm'))
    
    # Input plot
    
    hl_df <- mexicana_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    if (nrow(jitter_data_mexicana) > 0) {
      
      boxplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_mexicana <- jitter_data_mexicana[jitter_data_mexicana$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        mexicana_Boxplot <- ggboxplot(data = plot_df_mexicana, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
          geom_jitter(data = jitter_data_mexicana[jitter_data_mexicana$Majority.protein.IDs == HID[i],],aes(x=TimePoint, y=Value, color = Replica),width = 0.20)+
          geom_line(data = line_data_mexicana, aes(x=TimePoint, y=Value,group = 1))+
          scale_fill_manual(values = colors_mexicana)+
          scale_color_manual(values = colors_mexicana_replica)+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        mexicana_Boxplot <- list(mexicana_Boxplot)
        
        boxplot_list <- c(boxplot_list, mexicana_Boxplot)
        
      }
      
      mexicana_Boxplot <- ggarrange(plotlist = boxplot_list)
      
    }
    
    print(mexicana_Boxplot)
  })
  
  output$mexicana_Dotplot <- renderPlot({
    
    colors_mexicana_replica <- colorRampPalette(c("#FFFFFF", "#6EB5FF", "#000000"))(6)[2:5]
    
    ## Retrieve highlight data
    
    hl_df <- mexicana_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_mexicana <- mexicana_MainDataTable1()
    
    ## Selection data
    
    jitter_data_mexicana <- mexicana_MainDataTable2()
    
    test_data_mexicana <- mexicana_TestDataTable2()
    
    # Dummy plot
    
    mexicana_Dotplot <- ggplot(data = jitter_data_mexicana, aes(TimePoint, Value))+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    # Input plot
    
    if (nrow(jitter_data_mexicana) > 0) {
      
      y_pos_dot <- max(jitter_data_mexicana$Value)+0.2
      
      dotplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_mexicana <- jitter_data_mexicana[jitter_data_mexicana$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        # Format test data
        
        plot_test_data_mexicana <- test_data_mexicana[test_data_mexicana$Majority.protein.IDs == HID[i],]
        plot_test_data_mexicana <- plot_test_data_mexicana[,c(4,5,3)]
        plot_test_data_mexicana$p.adj <- format.pval(plot_test_data_mexicana$p.adj, digits = 3)
        
        mexicana_Dotplot <- ggdotplot(data = jitter_data_mexicana[jitter_data_mexicana$Majority.protein.IDs == HID[i],], x="TimePoint", y="Value", fill = "Replica",
                                      position = position_jitter(0.2))+
          geom_line(data = line_data_mexicana, aes(x=TimePoint, y=Value,group = 1),alpha = 0.8)+
          stat_pvalue_manual(plot_test_data_mexicana, y.position = y_pos_dot, step.increase = 0.1,label = "p.adj")+
          scale_fill_manual(values = colors_mexicana_replica)+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),
                axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        mexicana_Dotplot <- list(mexicana_Dotplot)
        
        dotplot_list <- c(dotplot_list, mexicana_Dotplot)
        
      }
      
      mexicana_Dotplot <- ggarrange(plotlist = dotplot_list)
      
    }
    
    print(mexicana_Dotplot)
  })
  
  #### Pairwise analysis  reference output 
  
  output$mexicana_Volcanoplot <- renderPlotly({
    
    #Input variables
    
    filter <- input$mexicana_ui_datasetID
    exp_a <- input$mexicana_TP
    exp_b <- input$mexicana_referenceTP
    ID <- input$mexicana_ui_selectedIDs
    enrich_s0 <- log2(2^as.numeric(input$mexicana_LFC)) 
    enrich_pvalue <- as.numeric(input$mexicana_Pvalue)
    
    # Fixed variables
    
    difference <- paste('difference', exp_a, exp_b, sep='_')
    pvalue <- paste('log10.pvalue', exp_a, exp_b, sep='_')
    val_count_a <- paste0('value_count_',exp_a) 
    val_count_b <- paste0('value_count_',exp_b)
    enriched_col <- paste('enriched', exp_a, exp_b, sep='_')
    
    min_quant_events <- 0 
    replicate_count <- 4
    enrich_c <- .05 
    volcano_threshold_linetype <- 2
    volcano_threshold_linesize <- 0.5
    difference_volcano_column <- 'mean2_naimp_meas_'
    volcano_column_type <- sub('([^[:digit:]]+).*', '\\1', difference_volcano_column)
    
    # Dummy plot
    
    ## Filter the data frame for data set user input.
    
    if (filter == paste0(unique(df_mexicana$Species),collapse = " & ")) {
      pg_quant <- full_mexicana
    }else{
      pg_quant <- full_mexicana %>% filter(Species == filter)
    }
    
    ## Change the column pattern to user friendly reading (i. e., 0h instead of inft000).
    
    if(length(input$mexicana_ui_timepointID) > 0) {
      mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
    }else{
      mexicana_time_df_ui <- mexicana_time_df
    }
    
    for (i in 1:length(mexicana_time_df_ui$hour)) {
      
      colnames(pg_quant) <- gsub(pattern = mexicana_time_df_ui$time_point[i], replacement = mexicana_time_df_ui$hour[i], colnames(pg_quant))
      
    }
    
    ## Dummy plot to show until user selection.
    
    dummy_df <- data.frame(x = c(-5,5),
                           y = c(0,10))
    
    volcanoplot <-
      ggplot(dummy_df, aes(x=x, y=y))+
      theme_minimal()+
      ylab("-log<sub>10</sub>(p-value)") + 
      xlab("log<sub>2</sub>(Fold change)") +
      theme(legend.text = element_text(size = 12))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(plot.title = element_text(size=20))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
    
    mexicana_Volcanoplot <-
      plotly::ggplotly(volcanoplot) %>%
      plotly::layout(showlegend = FALSE)
    
    # Input plot
    
    ## Upon proper user selection (if) produce the volcano plot.
    
    if (difference %in% colnames(pg_quant)) {
      
      ## Highlight user input, if any.
      
      pg_quant$highlight <- FALSE
      
      if (length(ID) != 0) {
        
        ids_searchstring <- paste(ID, collapse='|')
        
        matching_rows <- grep(ids_searchstring, as.character(pg_quant$Protein.IDs))
        
        pg_quant$highlight[matching_rows] <- TRUE
        
      }
      
      
      ## Filter count values. Standard set to 0, so they can check the dot all along the time course.
      ## In Mario's standard script is usually set to 2. But in this case, since it shown in the 
      ## box plot and in the dot plot, would be confusing not to show it here as well.
      ## Can be adjusted on the fixed variables section of this function.
      
      filter_min_quant_events <-
        pg_quant[val_count_a] >= min_quant_events |
        pg_quant[val_count_b] >= min_quant_events
      
      pg_quant <- pg_quant[filter_min_quant_events,]
      
      ## Select and highlight the enriched dots.
      
      my_enriched <- Enriched(pg_quant[c(difference, pvalue)],
                              c=enrich_c,
                              s0=enrich_s0,
                              pvalue=enrich_pvalue,
                              linetype=volcano_threshold_linetype,
                              size=volcano_threshold_linesize)
      
      pg_quant$enriched <- as.factor(my_enriched$enriched)
      
      
      my_points_nohighlight <- interaction(pg_quant[!pg_quant$highlight, c('enriched', 'highlight')])
      my_points_highlight <- interaction(pg_quant[pg_quant$highlight, c('enriched', 'highlight')])
      
      ## Generate a user friendly plotly string for the hoover.
      
      plotly_string <- paste0(sprintf('log2(foldchange): %.2f, -log10(pvalue): %.2f<br />Protein.IDs: %s<br />my_label: %s',
                                      pg_quant[[difference]], pg_quant[[pvalue]],
                                      sub('(.{30}).*', '\\1...', pg_quant$Protein.IDs),
                                      pg_quant$my_label),
                              switch('Gene.names' %in% names(pg_quant),
                                     sprintf('<br />Gene.names: %s', pg_quant$Gene.names),
                                     NULL),
                              switch('Protein.names' %in% names(pg_quant),
                                     sprintf('<br />Protein.names: %s', pg_quant$Protein.names),
                                     NULL),
                              switch('description' %in% names(pg_quant),
                                     sprintf('<br />description: %s', pg_quant$description),
                                     NULL))
      
      pg_quant$my_text <- plotly_string
      
      ## Subset the data for the volcano plot.
      
      pg_quant_subset <-
        pg_quant %>%
        rowwise() %>%
        mutate(imputation=ifelse(max(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count &
                                   min(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count,
                                 'measured',
                                 ifelse(min(range(.data[[val_count_a]], .data[[val_count_b]])) > 0,
                                        'some imputed',
                                        'ref imputed'))) %>%
        dplyr::select(Protein.IDs, my_label, matches(difference), matches(pvalue),
               my_text, highlight, enriched, imputation)
      
      ## Plot
      
      ## Select maximum x and y values of your volcano plot.
      
      max_x_value <- max(abs(pg_quant[grep('^difference', names(pg_quant))]),
                         na.rm=TRUE)
      max_y_value <- max(c(unlist(pg_quant[grep('^log10.pvalue', names(pg_quant))]),
                           10^enrich_pvalue + 0.5))
      
      volcanoplot <-
        ggplot(pg_quant_subset, aes_string(x=difference, y=pvalue, text='my_text')) +
        geom_hline(yintercept=0, color='grey80', na.rm=TRUE) +
        geom_vline(xintercept=0, color='grey80', na.rm=TRUE) +
        plot(my_enriched) +
        geom_point(data=subset(pg_quant_subset, !highlight),
                   aes(alpha=my_points_nohighlight,
                       color=my_points_nohighlight,
                       fill=my_points_nohighlight,
                       size=my_points_nohighlight),
                   na.rm=TRUE) +
        geom_point(data=subset(pg_quant_subset, highlight),
                   aes(alpha=my_points_highlight,
                       color=my_points_highlight,
                       fill=my_points_highlight,
                       size=my_points_highlight),
                   na.rm=TRUE) +
        scale_alpha_manual(values=c(`TRUE.FALSE`=.7, `FALSE.FALSE`=.4, `TRUE.TRUE`=.7, `FALSE.TRUE`=.7),
                           guide="none") +
        scale_color_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                           guide="none") +
        scale_fill_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                          guide="none") +
        scale_size_manual(values=c(`TRUE.FALSE`=1, `FALSE.FALSE`=.5, `TRUE.TRUE`=1.2, `FALSE.TRUE`=1),
                          guide="none")+
        coord_cartesian(xlim=c(-max_x_value,max_x_value), ylim=c(-.8,max_y_value)) +
        ylab("-log<sub>10</sub>(p-value)") + 
        xlab("log<sub>2</sub>(Fold change)") +
        annotate('text', label=exp_a, x=max_x_value / 2, y=-.4,size = 6) +
        annotate('text', label=exp_b, x=-max_x_value / 2, y=-.4,size = 6)+
        theme_minimal()+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      mexicana_Volcanoplot <-
        plotly::ggplotly(volcanoplot,tooltip='text') %>%
        plotly::layout(showlegend = FALSE)
      
    }
    
    print(mexicana_Volcanoplot)
   
  })
  
  #### Clustering ####
  
  output$clustering_MainDataTable <- renderDataTable({
    
    clustering_MainDataTable1() %>%
      datatable(options = list(pageLength = 7, autoWidth = TRUE), 
                rownames = F,escape = FALSE) %>%
      formatRound("mean(LFQ)", digits = 2)
    
  })
  
  clustering_MainDataTable_proxy <- dataTableProxy('clustering_MainDataTable')
  
  #### Clustering analysis output 
  
  output$clustering_LinePlot <- renderPlot({
    
    # Load SOM model
    
    if (input$clustering_ui_datasetID == "Leishmania spp.") {
      
      som_model <- cluster_leishmania_model
      
    }
    
    if (input$clustering_ui_datasetID == "M. musculus") {
      
      som_model <- cluster_mouse_model
      
    }
    
    nclust <- unique(som_model$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    # This plot will always show, wether the user selects an input ID or not
    
    mycolors <- brewer.pal(n = 9, name = "Set3")
    
    plot_list <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model$data[[1]][which(som_model$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$TimePoint <- gsub(pattern = "h", replacement = "", x = data$TimePoint)
      
      data$TimePoint <- factor(data$TimePoint, levels = c("0.5", "2", "6", "12","24", "48", "72"))
      
      plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
        geom_line(color="lightgrey")+
        stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1, linewidth = 2)+
        scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
        ylab("z-score(LFQ)")+
        xlab("Hours post-infection")+
        ggtitle(paste0("Cluster ",i))+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(plot.title = element_text(size=22))+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      plot <- list(plot)
      
      plot_list <- c(plot_list, plot)
      
    }
    
    # Input plot
    
    # User dependant plot. The IDs the user select will show up in this plot, when present.
    
    hl_df <- clustering_MainDataTable2()
    
    if (nrow(hl_df) > 0) {
      
      if (input$clustering_ui_datasetID == "Leishmania spp.") {
        
        HID <- unique(hl_df$`Protein ID`)
        
      }
      
      
      if (input$clustering_ui_datasetID == "M. musculus") {
        
        HID <- unique(paste0(hl_df$`Protein ID`, "_", gsub(pattern = "\\. ", replacement = "", x = hl_df$Experiment)))
        
      }
      
      plot_list <- c()
      
      for (i in 1:length(nclust)) {
        
        data <- as.data.frame(som_model$data[[1]][which(som_model$unit.classif==i),])
        
        data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
        
        data$TimePoint <- gsub(pattern = "h", replacement = "", x = data$TimePoint)
        
        data$TimePoint <- factor(data$TimePoint, levels = c("0.5", "2", "6", "12","24", "48", "72"))
        
        dataHL <- data[data$ID %in% HID,]
        
        if (nrow(dataHL) > 0) {
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1, linewidth = 2)+
            scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
            ylab("z-score(LFQ)")+
            xlab("Hours post-infection")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          for (j in 1:length(HID)) {
            
            colors_ui <- c(mycolors[i], viridis(n = length(HID),option = "G"))
            
            names(colors_ui) <- c("Mean", HID)
            
            plot <- plot +
              geom_line(data=data[data$ID == HID[j],], aes(x=TimePoint, y=Value, group=ID,color = ID))+
              scale_color_manual(name = "", values = colors_ui)
            
            
          }
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }else{
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1, linewidth = 2)+
            scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
            ylab("z-score(LFQ)")+
            xlab("Hours post-infection")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }
        
      }
      
      ClusterLinePlot <- ggarrange(plotlist = plot_list)
      
    }
    
    ClusterLinePlot <- ggarrange(plotlist = plot_list)
    
    print(ClusterLinePlot)
    
  })
  
  output$clustering_BarPlot <- renderPlot({
    
    color_species <- c("#A3C585", "#FFC9DE", "#6EB5FF")
    
    if (input$clustering_ui_datasetID == "Leishmania spp.") {
      
      som_model <- cluster_leishmania_model
      
      rel_set <- 1:dim(som_model$codes[[1]])[1]
      
      df <- data.frame()
      for (i in rel_set) {
        names <- rownames(som_model$data[[1]])[which(som_model$unit.classif==i)]
        clust <- i
        loop_df <- as.data.frame(cbind(names, clust))
        df <- rbind(df, loop_df)
      }
      
      df$Species <- NA
      
      df[grepl(pattern = "LINF",x = df$names),]$Species <- "L. infantum"
      df[grepl(pattern = "LmjF",x = df$names),]$Species <- "L. major"
      df[grepl(pattern = "LmxM",x = df$names),]$Species <- "L. mexicana"
      
      colnames(df)[1] <- "Protein.ID"
      
      order  <- df %>% group_by(clust) %>% tally() %>% arrange(n)
      
      order <- order$clust
      
      barplot_df <- df %>% group_by(clust,Species) %>% tally()
      
      barplot_df$clust <- factor(barplot_df$clust, levels = rev(unique(order)))
      
      barplot <- ggplot(barplot_df, aes(fill=Species, y = n , x= clust, label = n))+
        geom_bar(position="stack", stat="identity")+
        geom_text(size = 6, position = position_stack(vjust = 0.5)
                  # color = c("black","black",
                  #           "black", "black",
                  #           "black", "black")
        )+
        scale_fill_manual(values = color_species)+
        ylab("Leishmania protein IDs")+
        guides(fill=guide_legend(title="Experiment"))+
        xlab("Cluster")+
        theme_minimal()+
        theme(legend.key.size = unit(1, 'cm'),
              legend.position = "top",
              legend.text = element_text(size = 14, face = "italic"),
              legend.title = element_text(size = 16, face = "bold"))+
        theme(axis.text = element_text(size = 18),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title =element_text(size=20,face="bold"))
      
    }
    
    if (input$clustering_ui_datasetID == "M. musculus") {
      
      som_model <- cluster_mouse_model
      
      rel_set <- 1:dim(som_model$codes[[1]])[1]
      
      df <- data.frame()
      for (i in rel_set) {
        names <- rownames(som_model$data[[1]])[which(som_model$unit.classif==i)]
        clust <- i
        loop_df <- as.data.frame(cbind(names, clust))
        df <- rbind(df, loop_df)
      }
      
      df$Species <- gsub(pattern = "L", 
                         replacement = "L. ",
                         x = gsub(pattern = ".*_",
                                  replacement = "",
                                  x = df$names))
      
      df$names <- gsub(pattern = "_.*",
                       replacement = "",
                       x = df$names)
      
      colnames(df)[1] <- "Protein.ID"
      
      order  <- df %>% group_by(clust) %>% tally() %>% arrange(n)
      
      order <- order$clust
      
      colnames(df)[3] <- "Experiment"
      
      df$Species <- "M. musculus"
      
      barplot_df <- df %>% group_by(clust, Experiment, Species) %>% tally()
      
      barplot_df$clust <- factor(barplot_df$clust, levels = rev(unique(order)))
      
      barplot <- ggplot(barplot_df, aes(fill=Experiment, y = n , x= clust, label = n))+
        geom_bar(position="stack", stat="identity")+
        geom_text(size = 6, position = position_stack(vjust = 0.5)
                  # color = c("black","black",
                  #           "black", "black",
                  #           "black", "black")
        )+
        scale_fill_manual(values = color_species)+
        ylab("Mouse protein IDs")+
        xlab("Cluster")+
        theme_minimal()+
        theme(legend.key.size = unit(1, 'cm'),
              legend.position = "top",
              legend.text = element_text(size = 14, face = "italic"),
              legend.title = element_text(size = 16, face = "bold"))+
        theme(axis.text = element_text(size = 18),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title =element_text(size=20,face="bold"))
    }
    
    print(barplot)
  })
  
  #### Orthology ####
  
  output$Ortho_table <- renderDataTable({orth_df})
  
  ortho_MainDataTable <- reactive({
    ortho_main_table <- orth_df 
    ortho_main_table
  })
  
  ortho_MainDataTable_proxy <- dataTableProxy('Ortho_table')
  
  output$OrthoPlot <- renderPlot({
    
    dummy <- df_infantum
    
    dummy$TimePoint <- gsub(pattern = "h", replacement = "", x = dummy$TimePoint)
    
    OrthoPlot <- ggplot(data = dummy, aes(TimePoint, Value))+
      ylab("log2(LFQ)")+
      xlab("Hours post-infection")+
      theme_minimal()+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    id_input <- input$ui_orthoSelectedIDs
    
    # id_input <- c("OG6_100012")
  
    df_all <- rbind(rbind(df_infantum, df_major), df_mexicana)
      
    id <- c()
    
    plot_df <- c()
    
    if (length(id_input) > 0) {
      
      for (i in 1:length(id_input)) {
        
        loop_id <- orth_df[orth_df$`OrthoMCL ID` == id_input[i],]
        
        df_loop <- df_all %>% filter(Majority.protein.IDs %in% loop_id$`Protein ID`)
        
        df_loop <- df_loop %>%
          group_by(TimePoint,Species, Majority.protein.IDs) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE))
        
        plot_df <- rbind(plot_df, df_loop)
        
      }
      
      infantum_id_number <- length(unique(plot_df[plot_df$Species == "L. infantum",]$Majority.protein.IDs))
      
      major_id_number <- length(unique(plot_df[plot_df$Species == "L. major",]$Majority.protein.IDs))
      
      mexicana_id_number <- length(unique(plot_df[plot_df$Species == "L. mexicana",]$Majority.protein.IDs))
      
      colors_infantum <- c()
      
      if (infantum_id_number == 0) {
        
        colors_infantum <- c()
        
      }else{
        
        if (infantum_id_number == 1) {
          
          colors_infantum <- "#A3C585"
          
        }else{
          
          colors_infantum <- colorRampPalette(c("#FFFFFF","#A3C585", "#000000"))(infantum_id_number+2)
          
          colors_infantum <- colors_infantum[-1]
          
          colors_infantum <- colors_infantum[-length(colors_infantum)]
          
        }
        
      }
      
      colors_major <- c()
      
      if (major_id_number == 0) {
        
        colors_major <- c()
        
      }else{
        
        if (major_id_number == 1) {
          
          colors_major <- "#FFC9DE"
          
        }else{
          
          colors_major <- colorRampPalette(c("#FFFFFF","#FFC9DE", "#000000"))(major_id_number+2)
          
          colors_major <- colors_major[-1]
          
          colors_major <- colors_major[-length(colors_major)]
          
        }
        
      }
      
      colors_mexicana <- c()
      
      if (mexicana_id_number == 0) {
        
        colors_mexicana <- c()
        
      }else{
        
        if (mexicana_id_number == 1) {
          
          colors_mexicana <- "#6EB5FF"
          
        }else{
          
          colors_mexicana <- colorRampPalette(c("#FFFFFF","#6EB5FF", "#000000"))(mexicana_id_number+2)
          
          colors_mexicana <- colors_mexicana[-1]
          
          colors_mexicana <- colors_mexicana[-length(colors_mexicana)]
          
        }
        
      }
      
      colors <- c(colors_infantum, colors_major, colors_mexicana)
      
      plot_df$TimePoint <- gsub(pattern = "h", replacement = "",x = plot_df$TimePoint)
      
      OrthoPlot <- ggplot(data=plot_df, aes(x=TimePoint, y=Value,shape =Majority.protein.IDs,  group=Majority.protein.IDs,color = Majority.protein.IDs)) +
        geom_line(linewidth = 1)+
        geom_point(size = 4)+
        guides(shape="none")+
        scale_color_manual("",values = colors)+
        ylab("log2(LFQ)")+
        xlab("Hours post-infection")+
        theme_minimal()+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=20,face="bold"))
      
    }
    
    print(OrthoPlot)
    
  })
  
  #### Observe functions ####
  
  ### L. infantum 
  
  observe({
    updateSelectizeInput(
      session,
      'infantum_ui_selectedIDs',
      choices=infantum_MainDataTable3()[['Protein ID']][as.numeric(input$infantum_MainDataTable_rows_selected)],
      selected=infantum_MainDataTable3()[['Protein ID']][as.numeric(input$infantum_MainDataTable_rows_selected)]
    )
  })
  
  observe({
    if(length(input$infantum_ui_timepointID) > 0) {
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
    }else{
      infantum_time_df_ui <- infantum_time_df
    }
    
    updateSelectInput(session = session,
                      inputId = 'infantum_referenceTP',
                      choices = infantum_time_df_ui$hour
    )
  })
  
  observe({
    if(length(input$infantum_ui_timepointID) > 0) {
      infantum_time_df_ui <- infantum_time_df %>% filter(hour %in% input$infantum_ui_timepointID)
    }else{
      infantum_time_df_ui <- infantum_time_df
    }
    
    updateSelectInput(session = session,
                      inputId = 'infantum_TP',
                      choices = infantum_time_df_ui$hour
    )
  })
  
  observeEvent(input$update_infantum_MainDataTable, {
    rows <-
      match(input$infantum_ui_selectedIDs,
            infantum_MainDataTable3()[['Protein ID']])
    selectRows(infantum_MainDataTable_proxy,
               as.numeric(rows))
  })
  
  ### L. major 
  
  observe({
    updateSelectizeInput(
      session,
      'major_ui_selectedIDs',
      choices=major_MainDataTable3()[['Protein ID']][as.numeric(input$major_MainDataTable_rows_selected)],
      selected=major_MainDataTable3()[['Protein ID']][as.numeric(input$major_MainDataTable_rows_selected)]
    )
  })
  
  observe({
    if(length(input$major_ui_timepointID) > 0) {
      major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
    }else{
      major_time_df_ui <- major_time_df
    }
    
    updateSelectInput(session = session,
                      inputId = 'major_referenceTP',
                      choices = major_time_df_ui$hour
    )
  })
  
  observe({
    if(length(input$major_ui_timepointID) > 0) {
      major_time_df_ui <- major_time_df %>% filter(hour %in% input$major_ui_timepointID)
    }else{
      major_time_df_ui <- major_time_df
    }
    
    updateSelectInput(session = session,
                      inputId = 'major_TP',
                      choices = major_time_df_ui$hour
    )
  })
  
  observeEvent(input$update_major_MainDataTable, {
    rows <-
      match(input$major_ui_selectedIDs,
            major_MainDataTable3()[['Protein ID']])
    selectRows(major_MainDataTable_proxy,
               as.numeric(rows))
  })
  
  ### L. mexicana 
  
  observe({
    updateSelectizeInput(
      session,
      'mexicana_ui_selectedIDs',
      choices=mexicana_MainDataTable3()[['Protein ID']][as.numeric(input$mexicana_MainDataTable_rows_selected)],
      selected=mexicana_MainDataTable3()[['Protein ID']][as.numeric(input$mexicana_MainDataTable_rows_selected)]
    )
  })
  
  observe({
    if(length(input$mexicana_ui_timepointID) > 0) {
      mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
    }else{
      mexicana_time_df_ui <- mexicana_time_df
    }
    
    updateSelectInput(session = session,
                      inputId = 'mexicana_referenceTP',
                      choices = mexicana_time_df_ui$hour
    )
  })
  
  observe({
    if(length(input$mexicana_ui_timepointID) > 0) {
      mexicana_time_df_ui <- mexicana_time_df %>% filter(hour %in% input$mexicana_ui_timepointID)
    }else{
      mexicana_time_df_ui <- mexicana_time_df
    }
    
    updateSelectInput(session = session,
                      inputId = 'mexicana_TP',
                      choices = mexicana_time_df_ui$hour
    )
  })
  
  observeEvent(input$update_mexicana_MainDataTable, {
    rows <-
      match(input$mexicana_ui_selectedIDs,
            mexicana_MainDataTable3()[['Protein ID']])
    selectRows(mexicana_MainDataTable_proxy,
               as.numeric(rows))
  })
  
  ### Clustering
  
  observe({
    updateSelectizeInput(
      session,
      'clustering_ui_selectedIDs',
      choices=clustering_MainDataTable1()[['Protein ID']][as.numeric(input$clustering_MainDataTable_rows_selected)],
      selected=clustering_MainDataTable1()[['Protein ID']][as.numeric(input$clustering_MainDataTable_rows_selected)]
    )
  })
  
  observeEvent(input$update_clustering_MainDataTable, {
    rows <-
      match(input$clustering_ui_selectedIDs,
            clustering_MainDataTable1()[['Protein ID']])
    selectRows(clustering_MainDataTable_proxy,
               as.numeric(rows))
  })
  
  ### Orthology
  
  observe({
    updateSelectizeInput(
      session,
      'ui_orthoSelectedIDs',
      choices=ortho_MainDataTable()[['OrthoMCL ID']][as.numeric(input$Ortho_table_rows_selected)],
      selected=ortho_MainDataTable()[['OrthoMCL ID']][as.numeric(input$Ortho_table_rows_selected)]
    )
  })
  
  observeEvent(input$update_ortho_MainDataTable, {
    rows <-
      match(input$ui_orthoSelectedIDs,
            ortho_MainDataTable()[['OrthoMCL ID']])
    selectRows(ortho_MainDataTable_proxy,
               as.numeric(rows))
  })
  
  
}

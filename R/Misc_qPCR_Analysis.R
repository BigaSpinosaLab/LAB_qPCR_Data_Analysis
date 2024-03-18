################################################################################
################################                ################################

#####################   Title:   Miscelaneous qPCR Analysis  ###################
#####################   Creator: Eric Canton Dominguez       ###################

################################                ################################
################################################################################


################################################################################
########################## 1. Install all Packages #############################
################################################################################

if(!require(tabulizer)){
  devtools::install_github("ropensci/tabulizer")
  library(tabulizer) # Extract the tables from the pdf file 
}

if(!require(pdftools)){
  install.packages("pdftools")
  library(pdftools) # Extract the text from the pdf file 
}

if(!require(stringr)){
  install.packages("stringr")
  library(stringr) # Recognize patterns
}

if(!require(openxlsx)){
  install.packages("openxlsx")
  library(openxlsx) # Save the results in an excel file
}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr) # Utility
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2) # Plots
}

if(!require(ggpubr)){
  install.packages("ggpubr")
  library(ggpubr) # Plot statistics
}

if(!require(ggsignif)){
  install.packages("ggsignif")
  library(ggsignif) # Plot statistics
}

if(!require(readxl)){
  install.packages("readxl")
  library(readxl) # Read excel files
}

if(!require(RCurl)){
  install.packages("RCurl")
  library(RCurl) # Read excel files
}

################################################################################
################################ 2. Functions ##################################
################################################################################

read_PDF <- function(pdf_path,
                     result_path,
                     exp_name,
                     wells,
                     analized_groups,
                     housekeeping_genes) {
  
  ############################################################################
  ############################## Sanity Checks ###############################
  ############################################################################
  
  # Sanity check 1: Paths  
  if (!file.exists(pdf_path)) {
    stop("This PDF file doesn't exist")
  }
  
  if (!file.exists(result_path)) {
    stop("This path doesn't exist")
  }
  
  # Sanity check 2: Mandatory parameters
  if (missing(pdf_path) || missing(result_path) || missing(exp_name) ||
      missing(wells) || missing(analized_groups) || missing(housekeeping_genes)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
  
  # Sanity check 3: Parameter values
  if (!wells == round(wells) || wells <= 0) {
    stop("The parameter *wells* must be a positive integer number.")
  }
  
  if (!is.character(analized_groups) || length(analized_groups) == 0) {
    stop("The parameter *analized_groups* must be a string vector.")
  }
  
  if (!is.character(housekeeping_genes) || length(housekeeping_genes) == 0) {
    stop("The parameter *housekeeping_genes* must be a string vector.")
  }
  
  
  ############################################################################
  ################################## Code ####################################
  ############################################################################
  
  
  # Extract the tables into a list of dataframes
  Dataframes <- extract_tables(pdf_path, method = "stream", 
                               output= "data.frame")
  
  # Extract the text from the pdf 
  pdf_text <- pdf_text(pdf_path)
  pdf_text_combined <- paste(pdf_text, collapse = "\n")
  
  # Use the pattern that always repeats before and after the gene name to 
  # extract that name
  
  pattern <- "Abs Quant/2nd Derivative Max for (.*?) \\(Abs Quant/2nd Derivative Max\\)"
  
  pattern_string <- str_extract_all(pdf_text_combined, pattern)
  
  # Keep the gene name without the pattern
  gene_names <- strsplit(pattern_string[[1]], 
                         split = " ", 
                         fixed= TRUE)
  
  gene_names <- unlist(lapply(gene_names, function(df){
    df[6]
  }))
  
  
  ## Keep the tables that contain CT values
  #Delete the unnecessary tables
  delete_dataframes <- function(list) {
    
    # Function that filters the list of dataframes
    filtered_list <- lapply(list, function(df) {
      if ("CP" %in% colnames(df)) {
        return(df)
      }
    })
    
    filtered_list <- filtered_list[!sapply(filtered_list, is.null)]
    
    return(filtered_list)
  }
  
  CT_List <- delete_dataframes(Dataframes)
  
  
  ## Join the dataframes that are split in two
  # Due to the placement of the tables in the file, we need to use the number of
  # wells to join the tables that were split
  
  combine_and_delete <- function(list, ind) {
    list[[ind]] <- rbind(list[[ind]], list[[ind + 1]])
    list[[ind + 1]] <- NULL
    list <- list[!sapply(list, is.null)]
    return(list)
  }
  
  
  # Verify that the lenght of the list is the same as the vector of gene names
  i <- 1
  while (i <= length(CT_List) - 1) {
    current_df <- CT_List[[i]]
    next_df <- CT_List[[i + 1]]
    
    # Verify if the dataframe has the same number of rows than wells
    if (nrow(current_df) != wells) {
      CT_List <- combine_and_delete(CT_List, i)
    } else {
      i <- i + 1
    }
  }
  
  # Delete possible 0 row dataframes
  CT_List <- Filter(function(df) {
    nrow(df) > 0
  }, CT_List)
  
  if (length(CT_List) == length(gene_names)) {
    
    # Asign gene names to the names of the dataframes
    names(CT_List) <- gene_names
  } else {
    # error handling (stop)
    stop("The number of dataframes is not the same as the amount of genes.")
  }
  
  ## List that contains the dataframes with the necessary values obtained from the pdf
  Result_List <- lapply(CT_List, function(df) {
    
    # Column with the position of the well of each sample
    pos_well <- c(1:wells)
    df <- cbind(df, Well = pos_well)
    
    # Column with the groups submitted by the user ( mirar which() )
    df$Group <- sapply(df$Name, function(name) {
      matching_group <- which(str_detect(name, analized_groups))
      if (length(matching_group) > 0) {
        return(analized_groups[matching_group[1]])
      } else {
        return("No match")
      }
    })
    
    
    # Column with the replicas of each group
    df$Rep <- ave(seq_along(df$Name), df$Name, df$Group, FUN = seq_along)
    
    return(df)
  })
  
  # Check if any of the groups introduced by the user is not in the pdf file
  missing <- analized_groups[!analized_groups %in% 
                               unlist(lapply(Result_List, function(df) unique(df$Group)))]
  if (length(missing) > 0) {
    warning(paste("Some groups were not found in the pdf file:", paste(missing, collapse = ", ")))
  }
  
  # Add a new column "HK". It contais information about the gene being a housekeeping.
  i = 1
  for (i in 1:(length(Result_List))) {
    
    is_hk <- housekeeping_genes[str_detect(gene_names[i], housekeeping_genes)]
    if (length(is_hk) > 0) {
      Result_List[[i]]$HK <- "Yes"
    }  else {
      Result_List[[i]]$HK <- "No"
    }
    
  }
  
  # Reorganize dataframes and delete H2O samples
  Result_List <- lapply(Result_List, function(df) {
    order <- c("Pos", "Group", "Name", "Rep", "CP", "Status", "HK")
    df <- df[, order]
    df <- subset(df, Group %in% analized_groups)
    
    return(df)
  })
  
  hs <- createStyle(textDecoration = "BOLD", 
                    fontColour = "white", 
                    fontSize = 12, 
                    fontName = "Arial Narrow",
                    fgFill = "slateblue")
  
  write.xlsx(x = Result_List, 
             file = paste0(result_path, "/", exp_name, "_qPCR_data.xlsx"),
             headerStyle = hs,
             rowNames = FALSE)
  
  save(Result_List, 
       result_path,
       exp_name,
       wells,
       analized_groups,
       housekeeping_genes, 
       file = paste0(result_path, "/", exp_name, ".RData"))
  
  return(Result_List)
  
}

read_excel <- function(Excel_path) { 
  
  sheets <- readxl::excel_sheets(Excel_path) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(Excel_path, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
  
  names(data_frame) <- sheets 
  
  return(data_frame)
} 

ddct_analysis <- function(data,
                          result_path,
                          exp_name,
                          housekeeping_genes,
                          control_variable){
  
  ############################################################################
  ############################## Sanity Checks ###############################
  ############################################################################
  
  # Sanity check 1: Paths
  if (!file.exists(result_path)) {
    stop("This path doesn't exist")
  }
  
  # Sanity check 2: Mandatory parameters
  if (missing(data) || missing(result_path) || missing(exp_name) ||
      missing(housekeeping_genes) || missing(control_variable)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
  
  # Sanity check 3: Parameter values
  if (!wells == round(wells) || wells <= 0) {
    stop("The parameter *wells* must be a positive integer number.")
  }
  
  if (!is.character(housekeeping_genes) || length(housekeeping_genes) == 0) {
    stop("The parameter *housekeeping_genes* must be a string vector.")
  }
  
  if (!is.character(control_variable) || length(control_variable) == 0) {
    stop("The parameter *control_variable* must be a string vector.")
  }
  
  ############################################################################
  ################################## Code ####################################
  ############################################################################
  
  
  sample_means <- lapply(Results, function(df) {
    means <- aggregate(CP ~ Name, df, FUN = "mean")
    return(means)
  })
  
  sample_means <-  sample_means[housekeeping_genes]
  hk_combined <- bind_rows(sample_means, .id = "Gene")
  hk_combined <- aggregate(CP ~ Name, hk_combined, FUN = "mean")
  
  results <- lapply(Results, function(df) {
    
    merged_df <- merge(df, hk_combined, by = "Name", all.x = TRUE)
    names(merged_df)[names(merged_df) == 'CP.x'] <- 'CP'
    names(merged_df)[names(merged_df) == 'CP.y'] <- 'CP.HK'
    
    merged_df$delta1 <- merged_df$CP - merged_df$CP.HK
    
    meanControl <- aggregate(delta1 ~ Group, merged_df, FUN= "mean")
    
    meanControl <- subset(meanControl, Group == control_variable)
    mean <- meanControl$delta1
    
    merged_df$delta2 <- merged_df$delta1 - mean
    
    merged_df$ddct <- 2^(-merged_df$delta2)
    
    return(merged_df)
  })
  
  final <- lapply(results, function(df){
    
    frecuencias <- table(df$Name)
    
    merged_df <- aggregate(cbind(delta1,delta2,ddct) ~ Name + Group, df, FUN = "mean") 
    
    sd <- aggregate(ddct ~ Name, df, FUN = "sd")
    
    merged_df$sd <- sd$ddct
    
    
    merged_df$ddct <- round(merged_df$ddct, 2)
    merged_df$sd <- round(merged_df$sd, 2)
    
    merged_df$sd[merged_df$Name %in% names(frecuencias[frecuencias < 3])] <- NA
    
    return(merged_df)
    
  })
  
  final <- final[!(names(final) %in% housekeeping_genes)]
  
  hs <- createStyle(textDecoration = "BOLD", 
                    fontColour = "white", 
                    fontSize = 14, 
                    fontName = "Arial Narrow", 
                    fgFill = "turquoise4")
  
  write.xlsx(x = final, 
             file = paste0(result_path, "/", exp_name, "_qPCR_ddct.xlsx"),
             headerStyle = hs,
             rowNames = FALSE)
  
  
  return(final)
}


ddct_plot <- function(data,
                      ddct_values, 
                      Genes_of_interest, 
                      title,
                      result_path,
                      Groups_of_interest,
                      analized_groups,
                      GroupTest,
                      PairsTest,
                      exp_name) {
  
  ############################################################################
  ############################## Sanity Checks ###############################
  ############################################################################
  
  # Sanity check 1: Paths
  if (!file.exists(result_path)) {
    stop("This path doesn't exist")
  }
  
  # Sanity check 2: Mandatory parameters
  if (missing(data) || missing(ddct_values) || missing(Genes_of_interest) ||
      missing(title) || missing(result_path)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
  
  # Sanity check 3: Parameter values
  if (!is.character(Genes_of_interest) || length(Genes_of_interest) == 0) {
    stop("The parameter *Genes_of_interest* must be a string vector.")
  }
  
  if (!is.character(title)) {
    stop("The parameter *title* must be a string vector.")
  }
  
  ############################################################################
  ################################## Code ####################################
  ############################################################################
  
  df_combined <- bind_rows(ddct_values, .id = "Gene")
  df_filtrado <- df_combined %>% filter(Gene %in% Genes_of_interest)
  
  df_filtrado$Group <- apply(df_filtrado, 1, function(row) {
    str_extract(row["Name"], paste(analized_groups, collapse = "|"))
  })
  
  df_filtrado <- df_filtrado[df_filtrado$Group %in% Groups_of_interest, ]
  
  res <- aggregate(ddct ~ Group + Gene, data = df_filtrado, FUN = "mean")
  res <- res%>% select (Gene, Group, ddct)
  
  sd <- aggregate(ddct ~ Group + Gene, data = df_filtrado, FUN = "sd")
  res$sd <- sd$ddct
  
  
  frecuencias <- table(df_filtrado$Gene, df_filtrado$Group)
  
  freq_values <- all(frecuencias >= 3)
  
  df_filtrado$Gene <- factor(df_filtrado$Gene, levels = Genes_of_interest)
  df_filtrado$Group <- factor(df_filtrado$Group, levels = Groups_of_interest)

  res$Gene <- factor(res$Gene, levels = Genes_of_interest)
  res$Group <- factor(res$Group, levels = Groups_of_interest)  
  
  plot <- ggplot(res, aes(x = Gene, y=ddct, fill=Group)) +
    ggtitle(title) +
    xlab("") +
    ylab("Relative mRNA levels") +
    
    scale_y_continuous(expand = expansion(mult = c(0, .15))) +
    
    geom_hline(yintercept=0) +
    
    theme_bw() +
    theme(plot.title = element_text(hjust = .5, vjust = 0.6, face = "bold", size = "15"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.title.y = element_text(face="bold"),
          axis.line = element_line(colour = "black", linewidth = 1),
          axis.text.x = element_text(face="bold", size=14, angle=45, 
                                     vjust = 0.6, colour="Black"),
          axis.text.y = element_text(colour = "black", face="bold", size=14))
  
  if (freq_values) {
    plot <- plot + geom_bar(stat="identity", position = position_dodge(.6), 
                            width = .5) +
      scale_fill_grey(start = .3, end = .7) +
      geom_errorbar(aes(ymin = ddct - sd, ymax = ddct + sd), 
                    position = position_dodge(.6), width = .2, color = "#2F2F2F")
    
    compare_means(ddct ~ Group, data = df_filtrado, group.by= "Gene", method = "anova")
    
    
    if (length(Groups_of_interest) > 2) {
      
      plot <- plot + stat_compare_means(aes(x = Gene, y = ddct), data = df_filtrado,
                                        method = GroupTest, vjust = -12, hjust = 0.5)      
      plot <- plot + geom_pwc(
        aes(x = Gene, y = ddct), data = df_filtrado, tip.length = .01,
        method = PairsTest, p.adjust.method = "BH", label = "p.adj.format",
        bracket.nudge.y = 0.02
      )
      
      plot <- facet(plot, facet.by="Gene", scales="free")
      
    } else {
      
      plot <- plot + geom_pwc(
        aes(x = Gene, y = ddct), data = df_filtrado, tip.length = .01,
        method = PairsTest, label = "p.adj.format",
        bracket.nudge.y = 0.02
      )
    }
    plot <- plot + geom_point(data = df_filtrado, aes(fill=Group), 
                              size = .8, color = "black",position = position_dodge(.6), show.legend = FALSE) +
                   scale_colour_grey(start = .2, end = .2) 
    
  } else {
    
    plot <- plot + geom_point(data = df_filtrado, aes(fill=Group, shape = Group), 
                              size = 1.5, color = "black",position = position_dodge(.6)) +
                   scale_colour_grey(start = .2, end = .2)
  }
  
  
  
  print(plot)
  
  ggsave(file = paste0(exp_name, "_ddct_plot.pdf"),
         plot= last_plot(),
         device = pdf,
         path=result_path,
         width = 6,
         height = 5)          
  
  return(plot)
}


complete_ddct_analysis <- function() {
  
  Results <- read_PDF(pdf_path = pdf_path, 
                      result_path = result_path, 
                      exp_name = exp_name, 
                      wells = wells, 
                      analized_groups = analized_groups, 
                      housekeeping_genes = housekeeping_genes)
  
  ddct_results <- ddct_analysis(data = Results,
                                result_path = result_path, 
                                exp_name = exp_name, 
                                housekeeping_genes = housekeeping_genes,
                                control_variable = control_variable)
  
  plot_result <- ddct_plot(data = Results,
                           ddct_values = ddct_results,
                           Genes_of_interest = Genes_of_interest,
                           title = title,
                           result_path = result_path,
                           Groups_of_interest = Groups_of_interest,
                           analized_groups = analized_groups,
                           GroupTest = GroupTest,
                           PairsTest = PairsTest,
                           exp_name = exp_name)
  
  qPCR_report(ct_data = Results,
              dd_data = ddct_results,
              pal = 1, # Number from 1 to 5
              output_dir = paste0(result_path),
              output_file = paste0(exp_name, "_qPCR_Analysis_Report"))
  
  Complete_list <- list(Results,
                        ddct_results)
  
  names(Complete_list) <- c("qPCR Results",
                            "DDCT results")
  return(Complete_list)
  
}


qPCR_experiments <- function(exp_paths, 
                             exp_names, 
                             Genes_of_interest, 
                             title,
                             result_path,
                             Groups_of_interest,
                             exp_name) {
  
  if (length(exp_paths) != length(exp_names)) {
    stop("'exp_paths' and 'exp_names' must be the same length.")
  }
  
  lista_dataframes <- lapply(exp_paths, function(archivo) {
    read_excel(archivo)
    
  })
  
  names(lista_dataframes) <- exp_names
  
  lista_reducida <- lapply(lista_dataframes, function(df){
    bind_rows(df, .id = "Gene")
  })
  
  lista_completa_reducida <- bind_rows(lista_reducida, .id = "Exp")
  
  lista_total_reducida <- aggregate(ddct ~ Group + Gene + Exp, data = lista_completa_reducida, FUN="mean")
  
  lista_total_reducida <- lista_total_reducida%>% select (Exp, Gene, Group, ddct)
  
  sd <- aggregate(ddct ~ Group + Gene, data = lista_completa_reducida, FUN = "sd")
  lista_total_reducida$sd <- sd$ddct
  
  
  lista_final <- lista_total_reducida %>%
    group_split(Gene) %>%
    setNames(unique(lista_total_reducida$Gene))
  
  lista_final <- lapply(lista_final, function(df) {
    dataframe <- as.data.frame(df)
    return(dataframe)
  })
  
  frecuencias <- table(lista_total_reducida$Gene, lista_total_reducida$Group)
  
  freq_values <- all(frecuencias >= 3)
  
  lista_total_reducida$Gene <- factor(lista_total_reducida$Gene, levels = Genes_of_interest)
  lista_total_reducida$Group <- factor(lista_total_reducida$Group, levels = Groups_of_interest)
  
  lista_total_reducida <- na.omit(lista_total_reducida)
  
  lista_mean <- aggregate(ddct ~ Group + Gene, data = lista_total_reducida, FUN = "mean")
  
  sd <- aggregate(ddct ~ Group + Gene, data = lista_total_reducida, FUN = "sd")
  lista_mean$sd <- sd$ddct
  
  ###### PLOT #####
  
  plot <- ggplot(lista_mean, aes(x = Gene, y=ddct, fill=Group)) +
    
    geom_point(data = lista_total_reducida, aes(color=Group), size = .8, position = 
                 position_dodge(.6), show.legend = FALSE) +
    scale_colour_grey(start = .2, end = .2) +
    
    ggtitle(title) +
    xlab("") +
    ylab("Relative mRNA levels") +
    
    scale_y_continuous(expand = expansion(mult = c(0, .15))) +
    
    geom_hline(yintercept=0) +
    
    theme_bw() 
  
  if (freq_values) {
    plot <- plot + geom_bar(stat="identity", position = position_dodge(.6), 
                            width = .5) +
      scale_fill_grey(start = .3, end = .7) +
      geom_errorbar(aes(ymin = ddct - sd, ymax = ddct + sd), 
                    position = position_dodge(.6), width = .2, color = "#2F2F2F")
    
    compare_means(ddct ~ Group, data = lista_total_reducida, group.by= "Gene", method = "anova")
    
    
    plot <- plot + geom_point(data = lista_total_reducida, aes(color=Group), size = .8, 
                              position = position_dodge(.6), show.legend = FALSE)
    
    
    if (length(Groups_of_interest) > 2) {
      
      plot <- plot + stat_compare_means(aes(x = Gene, y = ddct), data = lista_total_reducida,
                                        method = "anova", vjust = -10, hjust = 1.2)      
      plot <- plot + geom_pwc(
        aes(x = Gene, y = ddct), data = lista_total_reducida, tip.length = .01,
        method = "t.test", p.adjust.method = "bonferroni", label = "p.adj.format",
        bracket.nudge.y = 0.02
      )
    } else {
      
      plot <- plot + geom_pwc(
        aes(x = Gene, y = ddct), data = lista_total_reducida, tip.length = .01,
        method = "t.test", label = "p.adj.format",
        bracket.nudge.y = 0.02
      )
    }
    
    
  }
  
  if (length(Genes_of_interest) > 2) {
    plot <- facet(plot, facet.by="Gene", scales="free") +
      theme(plot.title = element_text(hjust = .5, vjust = 0.6, face = "bold", size = "15"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            axis.title.y = element_text(face="bold"),
            axis.line = element_line(colour = "black", linewidth = 1),
            axis.text.x = element_text(size=0),
            axis.text.y = element_text(colour = "black", face="bold", size=12))
  } else {
    plot <- plot +
      theme(plot.title = element_text(hjust = .5, vjust = 0.6, face = "bold", size = "15"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            axis.title.y = element_text(face="bold"),
            axis.line = element_line(colour = "black", linewidth = 1),
            axis.text.x = element_text(face="bold", size=14, angle=45, 
                                       vjust = 0.6, colour="Black"),
            axis.text.y = element_text(colour = "black", face="bold", size=14))
  }
  
  ggsave(file = paste0(exp_name, "_ddct_Experiments_plot.pdf"),
         plot= last_plot(),
         device = pdf,
         path = result_path,
         width = 6,
         height = 5)            
  
  hs <- createStyle(textDecoration = "BOLD", fontColour = "white", fontSize = 14, 
                    fontName = "Arial Narrow", fgFill = "turquoise4")
  
  write.xlsx(x = lista_final, 
             file = paste0(result_path, "/", exp_name, "_qPCR_ddct_Experiments.xlsx"),
             headerStyle = hs,
             rowNames = FALSE)
  
  
  return(lista_final)
  
}

qPCR_report <- function(ct_data,
                        dd_data,
                        output_format = rmarkdown::html_document(toc = TRUE, 
                                                                 theme = "flatly", 
                                                                 toc_float = TRUE),
                        pal = 1,
                        output_dir = getwd(),
                        output_file = "qPCR_Analysis_Report",
                        title = "qPCR Analysis Report") 
{
  
  ############################################################################
  ############################## Sanity Checks ###############################
  ############################################################################  
  
  if (!file.exists(output_dir)) {
    stop("This path doesn't exist")
  }
  
  if (missing(ct_data) || missing(dd_data)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
 
  if (!is.character(title)) {
    stop("The parameter *title* must be a string vector.")
  }
  
  ############################################################################
  ################################## Code ####################################
  ############################################################################
  
  URL <- "https://raw.githubusercontent.com/BigaSpinosaLab/LAB_qPCR_Data_Analysis/main/Rmd/qPCR_Analysis_Report.Rmd"

  knitrRmd <- paste(readLines(textConnection(getURL(URL))), collapse="\n")

  setwd(tempdir())
  fn=tempfile()
  fn.Rmd <- paste(fn, ".Rmd", sep="")
  cat(knitrRmd, file=fn.Rmd)
  tf <- fn.Rmd

  suppressWarnings(
    rmarkdown::render(
      input = tf,
      output_format = output_format,
      output_dir = output_dir,
      output_file = output_file,
      params = list(ct_data = ct_data, dd_data = dd_data, pal = pal)
    )

  )

  file.remove(tf)

  report_path <- path.expand(file.path(paste0(output_dir, "/", output_file, ".html")))
  utils::browseURL(report_path)
  

}

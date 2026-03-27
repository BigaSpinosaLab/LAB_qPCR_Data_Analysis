#------------------------------------------------------------------------------#
####          Title :   Miscelaneous qPCR Analysis                          ####
####    Description :   Functions used in the qPCR Analysis script          ####
####         Author :   Eric Canton Dominguez                               ####
####           Date :   Updated March 27th 2026                             ####
#------------------------------------------------------------------------------#


#-----------------------------------#
#####  1. Install all Packages  #####
#-----------------------------------#

required_packages <- c("pdftools", "stringr", "openxlsx", "dplyr", "tidyverse",
                       "ggplot2", "ggpubr", "ggsignif", "readxl", "RCurl")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

invisible(lapply(required_packages, library, character.only = TRUE))

#-----------------------------------#
#####        2. Functions       #####
#-----------------------------------#

#' Extract and Process qPCR Data from PDF
#'
#' This function parses specific qPCR PDF reports, extracts result tables using regular expressions, 
#' assigns experimental groups, and exports the processed data to an Excel file.
#'
#' @param pdf_path            String. The full file path to the source PDF.
#' @param result_path         String. The path where the generated Excel file will be saved.
#' @param exp_name            String. The experiment name (used as a prefix for the output file).
#' @param analyzed_groups     Character vector. Experimental group names to be identified in the 'Name' column.
#' @param housekeeping_genes  Character vector. Genes to be flagged as housekeeping.
#'
#' @return  A named list of data.frames (one per detected gene). 
#'          It also saves an `.xlsx` file to the specified `result_path`.
#' @export

read_PDF <- function(pdf_path, result_path, exp_name,  
                     analyzed_groups, housekeeping_genes) {
  
  #-----------------------#
  ##### Sanity Checks #####
  #-----------------------#
  
  # Sanity check 1: Paths  
  if (!file.exists(pdf_path)) stop("This PDF file doesn't exist")
  if (!file.exists(result_path)) stop("This path doesn't exist")
  
  # Sanity check 2: Mandatory parameters
  if (missing(pdf_path) || missing(result_path) || missing(exp_name) || missing(analyzed_groups) || missing(housekeeping_genes)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
  
  if (!is.character(analyzed_groups) || length(analyzed_groups) == 0) {
    stop("The parameter *analyzed_groups* must be a string vector.")
  }
  
  if (!is.character(housekeeping_genes) || length(housekeeping_genes) == 0) {
    stop("The parameter *housekeeping_genes* must be a string vector.")
  }
  
  #-----------------------#
  #####      Code     #####
  #-----------------------#
  
  analyzed_groups <- unique(analyzed_groups)
  analyzed_groups_sorted <- analyzed_groups[order(nchar(analyzed_groups), decreasing = TRUE)]
  
  pdf_text_raw <- pdftools::pdf_text(pdf_path)
  pdf_combined <- paste(pdf_text_raw, collapse = "\n")
  
  pattern_gene <- "Abs Quant/2nd Derivative Max for (.*?) \\(Abs Quant/2nd Derivative Max\\)"
  gene_matches <- str_match_all(pdf_combined, pattern_gene)
  gene_names <- gene_matches[[1]][, 2]
  
  content_blocks <- str_split(pdf_combined, pattern_gene)[[1]][-1]
  
  
  if (length(content_blocks) > 0) {
    content_blocks[length(content_blocks)] <- str_split(
      content_blocks[length(content_blocks)], 
      "Tm Calling"
    )[[1]][1]
  }
  
  if (length(content_blocks) != length(gene_names)) {
    stop("Error: The number of genes and tables differ.")
  }
  
  Result_List <- lapply(seq_along(content_blocks), function(i) {
    block <- content_blocks[i]
    lines <- str_split(block, "\n")[[1]]
    
    regex_pattern <- "^\\s*\\S?\\s+([A-Z]\\d{1,2})\\s+(.*?)\\s{2,}\\S+\\s+([\\d.]+)(?:\\s+(\\S+))?\\s*$"
    
    all_matches <- str_match(lines, regex_pattern)
    
    df_all <- as.data.frame(all_matches) %>%
      filter(!is.na(V2)) %>%
      select(Pos = V2, Name = V3, CP_raw = V4, Status_raw = V5) %>%
      mutate(
        Name = str_trim(Name),
        CP = as.numeric(CP_raw), 
        Status = str_trim(Status_raw)
      ) %>%
      select(Pos, Name, CP, Status)
    
    df_final <- df_all %>%
      mutate(
        Group = map_chr(Name, function(x) {
          matches <- analyzed_groups_sorted[str_detect(x, fixed(analyzed_groups_sorted))]
          if (length(matches) > 0) return(matches[1]) else return(NA_character_)
        })
      ) %>%
      filter(!is.na(Group)) %>%
      group_by(Name) %>% 
      mutate(
        Rep = row_number(),
        HK = if_else(gene_names[i] %in% housekeeping_genes, "Yes", "No")
      ) %>%
      ungroup()
    
    df_final <- df_final[c("Pos", "Group", "Name", "Rep", "CP", "Status", "HK")]
    return(df_final)
  })
  
  names(Result_List) <- gene_names
  
  hs <- createStyle(textDecoration = "BOLD", fontColour = "white", fgFill = "slateblue")
  write.xlsx(Result_List, file = paste0(result_path, "/", exp_name, "_qPCR_data.xlsx"), headerStyle = hs)
  
  return(Result_List)
}


#' The read_excel function reads a .xlsx file created by the read_pdf function. 
#' Used if data was already converted to xlsx
#' 
#' @param excel_path            Path containing the xlsx data
#' @returns                     Excel file with pdf information
#' @export

read_excel <- function(excel_path) { 
  
  sheets <- readxl::excel_sheets(excel_path) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(excel_path, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
  
  names(data_frame) <- sheets 
  
  return(data_frame)
} 


#' Perform Delta-Delta Ct (ddCt) Analysis
#'
#' This function calculates relative gene expression by normalizing target genes 
#' to housekeeping controls and comparing them against a reference group.
#'
#' @param data                List. The list of data.frames returned by `read_PDF`.
#' @param result_path         String. Path where the ddCt results will be saved.
#' @param exp_name            String. The experiment name (used as a prefix for the output file).
#' @param housekeeping_genes  Character vector. Genes used for normalization.
#' @param control_variable    String. The name of the group used as the control.
#'
#' @return A list of data.frames containing mean Ct values, dct, ddct, and the final 
#'         2^-ddCt (fold change) values.
#' @export


ddct_analysis <- function(data,
                          result_path,
                          exp_name,
                          housekeeping_genes,
                          control_variable){
  
  #-----------------------#
  ##### Sanity Checks #####
  #-----------------------#
  
  # Sanity check 1: Paths
  if (!file.exists(result_path)) {
    stop("This path doesn't exist")
  }
  
  # Sanity check 2: Mandatory parameters
  if (missing(data) || missing(result_path) || missing(exp_name) ||
      missing(housekeeping_genes) || missing(control_variable)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
  
  if (!is.character(housekeeping_genes) || length(housekeeping_genes) == 0) {
    stop("The parameter *housekeeping_genes* must be a string vector.")
  }
  
  if (!is.character(control_variable) || length(control_variable) == 0) {
    stop("The parameter *control_variable* must be a string vector.")
  }
  
  #-----------------------#
  #####      Code     #####
  #-----------------------#
  
  
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


#' Plot Relative mRNA Levels (ddCt)
#'
#' Creates a standardized plot (Bar or Dot plot) with 
#' automated statistical comparisons between experimental groups.
#'
#' @param data                List. Raw qPCR data list.
#' @param ddct_values         List. Processed data from `ddct_analysis`.
#' @param genes_of_interest   Character vector. Genes to be included in the plot.
#' @param title               String. The main title of the plot.
#' @param result_path         String. Path to save the resulting PDF plot.
#' @param groups_of_interest  Character vector. Specific experimental groups to include.
#' @param group_test           String. Method for global comparison (e.g., "anova").
#' @param pairs_test           String. Method for pairwise comparison (e.g., "t.test").
#'
#' @return A ggplot object. Also saves a PDF file to the `result_path`.
#' @export

ddct_plot <- function(data,
                      ddct_values, 
                      genes_of_interest, 
                      title,
                      result_path,
                      groups_of_interest,
                      analyzed_groups,
                      group_test,
                      pairs_test,
                      exp_name) {
  
  #-----------------------#
  ##### Sanity Checks #####
  #-----------------------#
  
  # Sanity check 1: Paths
  if (!file.exists(result_path)) {
    stop("This path doesn't exist")
  }
  
  # Sanity check 2: Mandatory parameters
  if (missing(data) || missing(ddct_values) || missing(genes_of_interest) ||
      missing(title) || missing(result_path)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
  
  # Sanity check 3: Parameter values
  if (!is.character(genes_of_interest) || length(genes_of_interest) == 0) {
    stop("The parameter *genes_of_interest* must be a string vector.")
  }
  
  if (!is.character(title)) {
    stop("The parameter *title* must be a string vector.")
  }
  
  #-----------------------#
  #####      Code     #####
  #-----------------------#
  
  df_combined <- bind_rows(ddct_values, .id = "Gene")
  df_filtrado <- df_combined %>% filter(Gene %in% genes_of_interest)
  
  df_filtrado$Group <- apply(df_filtrado, 1, function(row) {
    str_extract(row["Name"], paste(analyzed_groups, collapse = "|"))
  })
  
  df_filtrado <- df_filtrado[df_filtrado$Group %in% groups_of_interest, ]
  
  res <- aggregate(ddct ~ Group + Gene, data = df_filtrado, FUN = "mean")
  res <- res%>% select (Gene, Group, ddct)
  
  sd <- aggregate(ddct ~ Group + Gene, data = df_filtrado, FUN = "sd")
  res$sd <- sd$ddct
  
  
  frecuencias <- table(df_filtrado$Gene, df_filtrado$Group)
  
  freq_values <- all(frecuencias >= 3)
  
  df_filtrado$Gene <- factor(df_filtrado$Gene, levels = genes_of_interest)
  df_filtrado$Group <- factor(df_filtrado$Group, levels = groups_of_interest)

  res$Gene <- factor(res$Gene, levels = genes_of_interest)
  res$Group <- factor(res$Group, levels = groups_of_interest)  
  
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
    
    
    if (length(groups_of_interest) > 2) {
      
      plot <- plot + stat_compare_means(aes(x = Gene, y = ddct), data = df_filtrado,
                                        method = group_test, vjust = -12, hjust = 0.5)      
      plot <- plot + geom_pwc(
        aes(x = Gene, y = ddct), data = df_filtrado, tip.length = .01,
        method = pairs_test, p.adjust.method = "BH", label = "p.adj.format",
        bracket.nudge.y = 0.02
      )
      
      plot <- facet(plot, facet.by="Gene", scales="free")
      
    } else {
      
      plot <- plot + geom_pwc(
        aes(x = Gene, y = ddct), data = df_filtrado, tip.length = .01,
        method = pairs_test, label = "p.adj.format",
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


#' Run Complete qPCR Analysis Pipeline
#'
#' Full workflow: parses PDF, performs ddCt calculations, 
#' generates plots, and produces an HTML report.
#'
#' @return A named list containing:
#'   \item{qPCR Results}{Raw data extracted from the PDF.}
#'   \item{DDCT results}{Processed relative expression data.}
#' @export

complete_ddct_analysis <- function() {
  
  Results <- read_PDF(pdf_path = pdf_path, 
                      result_path = result_path, 
                      exp_name = exp_name, 
                      analyzed_groups = analyzed_groups, 
                      housekeeping_genes = housekeeping_genes)
  
  ddct_results <- ddct_analysis(data = Results,
                                result_path = result_path, 
                                exp_name = exp_name, 
                                housekeeping_genes = housekeeping_genes,
                                control_variable = control_variable)
  
  plot_result <- ddct_plot(data = Results,
                           ddct_values = ddct_results,
                           genes_of_interest = genes_of_interest,
                           title = title,
                           result_path = result_path,
                           groups_of_interest = groups_of_interest,
                           analyzed_groups = analyzed_groups,
                           group_test = group_test,
                           pairs_test = pairs_test,
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


#' Aggregate and Analyze Multiple qPCR Experiments
#'
#' Reads multiple Excel files (generated by read_PDF), merges the data, 
#' calculates means and SDs, and generates a comparative plot.
#'
#' @param exp_paths           Character vector. Paths to the .xlsx files to be aggregated.
#' @param exp_names           Character vector. Labels for each experiment (must match exp_paths length).
#' @param genes_of_interest   Character vector. Genes to be included in the final analysis.
#' @param title               String. Title for the comparative plot.
#' @param result_path         String. Directory to save the aggregated Excel and PDF plot.
#' @param groups_of_interest  Character vector. Subset of groups to display.
#' @param exp_name            String. Prefix for the output filenames.
#'
#' @return A list of data.frames grouped by gene. Also saves a merged .xlsx 
#'   file and a .pdf plot to the result_path.
#' @export

qPCR_experiments <- function(exp_paths, 
                             exp_names, 
                             genes_of_interest, 
                             title,
                             result_path,
                             groups_of_interest,
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
  
  lista_total_reducida$Gene <- factor(lista_total_reducida$Gene, levels = genes_of_interest)
  lista_total_reducida$Group <- factor(lista_total_reducida$Group, levels = groups_of_interest)
  
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
    
    
    if (length(groups_of_interest) > 2) {
      
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
  
  if (length(genes_of_interest) > 2) {
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


#' Generate HTML qPCR Analysis Report
#'
#' Fetches a remote RMarkdown template to generate a comprehensive HTML report 
#' including tables and visualizations. Automatically opens the report in the browser.
#'
#' @param ct_data       List. Raw Ct data from `read_PDF`.
#' @param dd_data       List. Processed ddCt data from `ddct_analysis`.
#' @param output_format RMarkdown output format.
#' @param pal           Integer. Color palette selection (1-5).
#' @param output_dir    String. Directory to save the HTML report.
#' @param output_file   String. Filename for the report (without extension).
#' @param title         String. Report header title.
#'
#' @return It creates an HTML file and opens it in the system's default web browser.
#' @export

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
  
  #-----------------------#
  ##### Sanity Checks #####
  #-----------------------#
  
  if (!file.exists(output_dir)) {
    stop("This path doesn't exist")
  }
  
  if (missing(ct_data) || missing(dd_data)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
 
  if (!is.character(title)) {
    stop("The parameter *title* must be a string vector.")
  }
  
  #-----------------------#
  #####      Code     #####
  #-----------------------#
  
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

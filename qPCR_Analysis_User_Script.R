#------------------------------------------------------------------------------#
####          Title :   qPCR Analysis User Script                           ####
####    Description :   Automated workflow for qPCR data extraction         ####
####                    (PDF/Excel), relative quantification using the      ####
####                    2^-ddCt method, plotting, and reporting.            ####
####         Author :   Eric Canton Dominguez                               ####
####           Date :   Updated March 27th 2026                             ####
####        Version :   2.0.0                                               ####
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Version 2.0 changes:  * Removed tabulizer dependency (deprecated)            #
#                       * Refactored variables and improved documentation      #
#                       * Optimized UI/Report layout                           #
#------------------------------------------------------------------------------#

# REMARK: MELTING TEMPERATURES ARE *NOT* CHECKED.                              #
# REMARK: OUTLIERS (IF ANY) ARE NOT IDENTIFIED FOR THE DDCT ANALYSIS.          #
#         CHECK HTML REPORT FOR EXPLORATORY ANALYSIS.                          #
#         (CP AND DDCT VALUES)                                                 #

#------------------------------------------------------------------------------#
####                          Summary of the script                         ####
#------------------------------------------------------------------------------#

### This script is divided in 3 parts:

#   1 - INTRODUCTION - Define all the parameters you want to use               #
#   2 - SOURCE       - Import all the functions                                #
#   3 - EXECUTION    - Execute the functions. A brief description of each      #
#                      function is shown below:                                #


#    1. Read_pdf: Use the read_pdf function to extract the qPCR data from a    #
#                 pdf file and export it to a .xlsx file. (JDK required)       #

#    2. Read_excel: Use the read_excel function to import the data directly    #  
#                   from a .xlsx file.                                         #

#    3. ddct_analysis: Use the ddct_analysis function to quickly perform       #
#                      the ddct method from an excel file or the result        # 
#                      from the pdf read by the read_pdf function.             #

#    4. qPCR_report: Use the qpcr_report function to generate a report of      #  
#                    the qPCR analysis.                                        #

#    5. ddct_plot: Use the ddct_plot to generate a plot comparing the values   #  
#                  of the different conditions.                                #

#    6. complete_ddct_analysis: Combination of the read_pdf, ddct_analysis     #
#                               and ddct_plot functions.                       #

#    7. qPCR_experiments: Combine 2 or more experiments and make a plot with   # 
#                         their ddct values.                                   #


#------------------------------------------------------------------------------#
#####                            1. INTRODUCTION                           #####
#                     This is the only section the user must fill              #
#------------------------------------------------------------------------------#

# Fill the lines that do not start with a "#".                                 #
# Do not change the name of the parameters, only change after the "=".         #
# If you change a parameter, remember to load the updated value into           #
# the environment.                                                             #

#------------------------------------------------------------------------------#
######    General variables  (Variables that you always need to fill)     ######
#------------------------------------------------------------------------------#
# Path where the results will be stored
# result_path <- "path_to_store_results"
result_path = "~/LAB_qPCR_Data_Analysis-main/Results"

# Name of the experiment (This will be used to name the files).
# exp_name <- "Name_of_the_experiment"
exp_name = "example_test"

# Housekeeping genes
# housekeeping_genes = c("ACTB", "...",...)
housekeeping_genes = c("GENE1")

# Names of all the conditions ("Groups") of the analysis
# If one group contains another ("SPP1" and "SPP12") put the longest one first.
# If there are special characters in a group ("SPP1+EPZ) put "\\" before the 
# special character ("SPP1\\+EPZ")
# analyzed_groups = c("UNT", "CAF", "...",...)
analyzed_groups = c("WT", "KO")


#------------------------------------------------------------------------------#
######                 FILE TO READ  (Choose PDF or excel)                 ######
#------------------------------------------------------------------------------#
# Path of the LightCycler 480 PDF Report 
# pdf_path <- "path_of_the_file"
pdf_path = "~/LAB_qPCR_Data_Analysis/Demo_Files/test_files/test.PDF"

# qPCR excel file
# Result = "path_of_the_excel_file"
excel_path = "~/LAB_qPCR_Data_Analysis/Demo_Files/Results/example_test_qPCR_data.xlsx"


#------------------------------------------------------------------------------#
######                       ddct_analysis variables                      ######
#------------------------------------------------------------------------------#
# Variable used as the control variable for the ddct analysis 
# control_variable = "UNT"
control_variable = "WT"


#------------------------------------------------------------------------------#
######                         ddct_plot variables                        ######
#------------------------------------------------------------------------------#
# Names of all the Genes to plot
# Genes_of_interest = c("NNMT1", "SERPINE1", "SNAI2", "THBS1")
Genes_of_interest = c("GENE1", "GENE2", "GENE3", "GENE4")

# Names of all the conditions ("Groups") to plot
# Groups_of_interest = c("WT", "KO", "...",...)
Groups_of_interest = c("WT", "KO")

# Test to perform if there are 3 or more groups to compare ("anova" or "kruskal.test")
group_test = "kruskal.test"

# Test to perform between groups ("t.test" or "wilcox.test")
pairs_test = "t.test"

# Title of the plot generated
# title = "Title"
title = "DDCT test results"


#------------------------------------------------------------------------------#
######                     qPCR_experiments variables                     ######
#------------------------------------------------------------------------------#
# List of all the experiments to mix
# exp_paths = list("~/Exp/Path/Exp1.xlsx", 
#                  "~/Exp/Path/Exp2.xlsx")
exp_paths = list("/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/Demo_Files/Results/xlsx/example_test_qPCR_ddct.xlsx",
                 "/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/Demo_Files/test_files/example_test2_qPCR_ddct.xlsx", 
                 "/home/user/Documents/Files/Projects/LAB_qPCR_Data_Analysis/Demo_Files/test_files/example_test3_qPCR_ddct.xlsx")

# List of all the experiments to mix
# exp_names = list("2023-11-10", "2023-12-10", "2023-12-20")
exp_names = c("Test1", "Test2", "Test3")



#------------------------------------------------------------------------------#
#####                        2. Source all functions                       #####
#------------------------------------------------------------------------------#

source("https://raw.githubusercontent.com/BigaSpinosaLab/LAB_qPCR_Data_Analysis/main/R/Misc_qPCR_Analysis.R") # All of the functions we will need 



#------------------------------------------------------------------------------#
#####                         3. Execute functions                         #####
#------------------------------------------------------------------------------#

# Read the pdf file of the qPCR
Results <- read_PDF(pdf_path = pdf_path, 
                    result_path = result_path, 
                    exp_name = exp_name, 
                    analyzed_groups = analyzed_groups, 
                    housekeeping_genes = housekeeping_genes)


# Read the excel file of the qPCR
Results <- read_excel(excel_path = excel_path)


# Perform a ddct analysis of the file
ddct_results <- ddct_analysis(data = Results,
                              result_path = result_path, 
                              exp_name = exp_name, 
                              housekeeping_genes = housekeeping_genes,
                              control_variable = control_variable)


# Create a report of the results
qPCR_report(ct_data = Results,
            dd_data = ddct_results,
            pal = 1, # Number from 1 to 5 (Different color palettes)
            output_dir = paste0(result_path),
            output_file = paste0(exp_name, "_qPCR_Analysis_Report"))


# Plot the results of the ddct analysis
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



# Perform the read_pdf, ddct_analysis and ddct_plot functions at one from the pdf file 
complete_ddct_analysis <- complete_ddct_analysis()

# Join multiple experiments and plot them
combine_qPCRs <- qPCR_experiments(exp_paths = exp_paths,
                              exp_names = exp_names,
                              genes_of_interest = genes_of_interest,
                              title = title,
                              result_path = result_path,
                              groups_of_interest = groups_of_interest,
                              exp_name = exp_name)

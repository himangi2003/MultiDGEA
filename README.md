# Multidgea

This package is specifically designed to carry out differential gene expression analysis using five different analysis methods that includes Limma Voom, Limma voom with duplicate correlation, DEseq, Kimma and DREM.


Author: Himangi Srivastava

The MultiDGEA package requires some of the packages dependency  that are used. The package specifically requires R version 4.2.

| Package | 
| -----------| 
| DESeq2   | 
|  dplyr | 
|doParallel|
| edgeR |
|ggplot2|
|kimma|
|limma|
|magrittr|
| readr|
| rlang|
|stringr|
|tibble|
|variancePartition|
| VennDiagram |
    
    

# Introduction
The package MultiDGEA is designed to carry out differetial gene expression analysis for the RNA seq data generated at different time points. The package includes the capacity to run and analyze RNASeq data to obtain Differentially expressed genes using different methods at the same time. It can also perform comparitive analysis of the results obtained from all different kinds of analysis with the capacity to store and summarize the results of each comparative analysis methods used.

The package is designed to carry out analysis by just extracting argument form an arguments file by using extract_input_data function which will automatically convert the data into the form required by the later function to run the analysis  or it can carry out analysis by dynamically assigning variables and inputs to the functions.
There are set of inputs used to run the analysis, some of it are specific to the particular analysis while others are common to all the anlysis types. The ones that are common are 
1. count matrix
2. sample meta data
3. contrast variable
4. patient id variable
5. model design parameters
   
In order to run the pacakage correctly the some of these inputs like Count matrix and sample meta data  must be preprocessed as a dataframe containing the right column names and row names while others must be correctly assigned as the string variable based on what kind of analysis the user wants to perform.

### Inputs
If the user is running the pacakage by extracting the inputs using the argument file then the sample argument file looks like this.
|argument	| value1 |
| -----------| ------------------- |
|data_file	|path to count matrix|
|sample_meta_data	|path to sample meta data|
|day_subset|	0,3|
|treatment_gp	|2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections),2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)|
|model design parameters	|~day+pubid|
|output_directory	|path to store the results|


The function extract_data_input will extract all the inputs from the arguments files which can be useful to run analysis in an iteration for all different kinds of analysis and many other different time points as a pipeline. However, to run the analysis as a pipeline or even separately there are few things that theuser has to keep in mind before running the analysis and that is the form of the input data structure.

This table explains the inputs used for the anlalyis and the structure and format that is required for the input to run the analysis.
|input data	| format and structure | analysis type that uses this input |
| -----------| ------------------- | ------------------- |
|count_matrix| R dataframe containing gene counts with genenames as rownames and the sample subset selected for the analysis as column names |all |
|sample_metadata| R dataframe containing information about each sample (day, participant id/ptid) with subset samples selected for the analysis as rownames and the information about samples like the day, participant ID as column names |all |
|contrast variable | string variable the is used to define the contrast variable which is the time point variable for the analysis, for example "day" | all|
|patient id variable  | string variable that is used to define the unique patient id variable used for the analysis, for example in the vaccine study it is the variable that stores participant id for each individual like "ptid" | all|
|model_design_parameter  | string variable that is used to define the model design parameter to be used to run the analysis. for example "~day+ptid",  note: it should be the same variable or combination of variable used to define the contrast variable and patient id variable for the analysis | all|
| form  | string variable that is used to define the variable to be tested for a fixed effect e.g. "~ day + (1/ptid) " | limma_voom with duplicate correlation |
| kimma_model | string variable that is used to define the model design parameter used to perform the kimma analysis e.g "day+ptid"| kimma|




## MultiDGEA R Package without using arguments file

The outline of this tutorial is as follows:
### Running the analysis 
1.  Install the package MultiDGEA from a zipped folder `` install.packages(file.choose(), repos = NULL, type="source") ``
2.  Use the library `` library(MultiDGEA) ``
3.  Load the  counts data 
   `` counts_data <- read.csv("counts.csv", row.names = 1)  ``
4. Load the sample meta data 
   `` sample_meta_data <- read.csv("sample_meta_data.csv", row.names = 1)  ``

6.  Find MultiDGEA analysis results using Kimma method

                Kimma_result <- fitting_data_kimma (count_matrix = counts_data,
                                                    meta_data = sample_meta_data,
                                                    model_design_parameter = "~day",
                                                    kimma_model=  "~ day + (1|pubid)",
                                                    contrast_var="day",
                                                    patient_id_var="pubid")
7.  Find MultiDGEA analysis results using DREM method

               Drem_results <- fitting_data_DREAM (count_matrix = counts_data,
                                                    meta_data = sample_meta_data,
                                                    model_design_parameter = "~day",
                                                    form = "~ day + (1|pubid) ",
                                                    contrast_var="day",
                                                    patient_id_var="pubid")
    
9.  Find MultiDGEA analysis results using Limma Voom method
   
                Limma_voom_results <- fitting_data_limma_voom (count_matrix = counts_data,
                                                    meta_data = sample_meta_data,
                                                    model_design_parameter = "~day",
                                                    contrast_var="day",
                                                    patient_id_var="pubid")

10. Find MultiDGEA analysis results using Limma Voom Duplicate Correlation method
   
                Limma_voom_duplicate_corr_results <- limma_voom_duplicate_correlation(count_matrix = counts_data,
                                                                                         meta_data = sample_meta_data,
                                                                                         model_design_parameter = "~day+pubid",
                                                                                         contrast_var="day",
                                                                                         patient_id_var="pubid")

11. Find MultiDGEA analysis results using Deseq method
   
                Deseq_result <- fitting_data_Deseq(count_matrix = counts_data,
                                        meta_data = sample_meta_data,
                                        model_design_parameter = "~day+pubid",
                                        contrast_var="day",
                                        patient_id_var="pubid")



### Visualization of the results

12.  Plot the normalized counts plot data set for different analysis:
  for  limma_voom


                        Normalized_data_plot(method="limma_voom",
                                            result=limma_voom_results,
                                                count_matrix = counts_data,
                                               meta_data = sample_meta_data,
                                                model_design_parameter ="~day")

For deseq


                        Normalized_data_plot(method="deseq",
                                            result=deseq_results,
                                                count_matrix = counts_data,
                                               meta_data = sample_meta_data,
                                                model_design_parameter ="~day") 


                    
14.  Plot the heatmap plot for differenr analysis

                     heatmap_top_100_MultiDGEA_edgeR(result=Limma_voom_results,
                                                 count_matrix = counts_data,
                                               meta_data = sample_meta_data,
                                                model_design_parameter ="~day")
                
16.  Plot the histogram plot for different analysis

                            histogram_plot_p_value(Deseq_results)
                          
18.  Plot the volcano plot data for different analysis

                        volcano_plot(Limma_voom_results,"volcano_plot_for_Limma_voom_duplicate_coorelation")

                        


### Comparision of different results
18. P_value scatter plot for comparing the p_values for different analysis

    for comparing limma_voom_duplicate_correlation and Deseq downregulated genes

                                    p_val_scatter_plot_comparison_down(res1=Limma_voom_duplicate_corr_results,
                                                                   res2=Deseq_results,
                                                                   p_val_xlabel="p_val_Limma_voom_duplicate_corr_results",
                                                                   p_val_ylabel="p_val_Deseq_result",
                                                                   p_val_title="Limma_voom_duplicate_corr_results_vs_Deseq_results",
                                                                   p_val=0.05,
                                                                   FDR_=0.2,
                                                                   LFC_=0.05)

                                   
    for comparing limma_voom_duplicate_correlation and Deseq upregulated genes

                                    p_val_scatter_plot_comparison_up(res1=Limma_voom_duplicate_corr_results,
                                                                   res2=Deseq_results,
                                                                   p_val_xlabel="p_val_Limma_voom_duplicate_corr_results",
                                                                   p_val_ylabel="p_val_Deseq_result",
                                                                   p_val_title="Limma_voom_duplicate_corr_results_vs_Deseq_results",
                                                                   p_val=0.05,
                                                                   FDR_=0.2,
                                                                   LFC_=0.05)
    

                               
20. LFC_value comparision for different comparision

for comparing limma_voom_duplicate_correlation and Deseq downregulated genes

    LFC_scatter_plot_comparison_down(res1=Limma_voom_duplicate_corr_results,
                                  res2=Deseq_result,
                                 LFC_xlabel="LFC_Limma_voom_duplicate_corr_results",
                                 LFC_ylabel="LFC_Deseq_result",
                                 LFC_title="Limma_voom_duplicate_corr_results_vs_Deseq_results",
                                 p_val=0.05,
                                 FDR_=0.2,
                                 LFC_=0.05)


for comparing limma_voom_duplicate_correlation and Deseq upregulated genes

    LFC_scatter_plot_comparison_up(res1=Limma_voom_duplicate_corr_results,
                                  res2=Deseq_results,
                                 LFC_xlabel="LFC_Limma_voom_duplicate_corr_results",
                                 LFC_ylabel="LFC_Deseq_result",
                                 LFC_title="Limma_voom_duplicate_corr_results_vs_Deseq_results",
                                 p_val=0.05,
                                 FDR_=0.2,
                                 LFC_=0.05)

                                 
22. Finding significant gene and plotting venn diagram representation for the comparing the number  of significant genes accross different analysis

Finding significant genes for Limma_voom_duplicate_correlation
                        
                    Limma_voom_duplicate_corr_results_MultiDGEA <- find_MultiDGEA_genes(result = Limma_voom_duplicate_corr_results,
                                                                             p_val=0.05,
                                                                              FDR_=0.2,
                                                                              LFC_=0.05)
Finding significant genes for Deseq
                        
                    Limma_voom_duplicate_corr_results_MultiDGEA <- find_MultiDGEA_genes(result = Limma_voom_duplicate_corr_results,
                                                         p_val=0.05,
                                                          FDR_=0.2,
                                                          LFC_=0.05)


                        Deseq_results_MultiDGEA <- find_MultiDGEA_genes(result = Deseq_results,
                                                            p_val=0.05,
                                                            FDR_=0.2,
                                                            LFC_=0.05)
                        
                        sig1=Deseq_results_MultiDGEA$up_regulated
                        sig2=Limma_voom_duplicate_corr_results_MultiDGEA$up_regulated
                        plot_venn(sig1,sig2,"Limma_voom_vs_Deseq")

                                                          



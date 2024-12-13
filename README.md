# MultiDGEA

This package is specifically designed to carry out differential gene expression analysis using five different analysis methods, including Limma Voom, Limma Voom with duplicate correlation, DESeq, Kimma, and DREM.

### Author: Himangi Srivastava

## Requirements
The MultiDGEA package requires R version 4.2 and several dependencies. The package is provided as `Multidgea_0.1.0.tar.gz`, along with sample input data:
- **Count matrix**: `count_matrix.csv`
- **Sample meta data**: `sample_meta_data.csv`

### Package Dependencies
| Package              |
|----------------------|
| DESeq2              |
| dplyr               |
| doParallel          |
| edgeR               |
| ggplot2             |
| kimma               |
| limma               |
| magrittr            |
| readr               |
| rlang               |
| stringr             |
| tibble              |
| variancePartition   |
| VennDiagram         |

---

## Introduction
MultiDGEA is designed for differential gene expression analysis of RNA-seq data across multiple time points. It can:
- Analyze RNA-seq data using multiple methods simultaneously.
- Perform comparative analyses of the results from different methods.
- Store and summarize the results.

### Features
1. **Automated Input Handling**: Extract arguments from an input file using `extract_input_data`.
2. **Dynamic Input Configuration**: Assign inputs dynamically for specific analyses.
3. **Comprehensive Output**: Summarized results for each method.

### Inputs
Common inputs include:
- **Count matrix**: Gene counts (rows: genes, columns: samples).
- **Sample metadata**: Sample details (e.g., day, participant ID).
- **Contrast variable**: Defines the time points (e.g., "day").
- **Patient ID variable**: Unique participant identifier (e.g., "ptid").
- **Model design parameters**: Defines the model (e.g., `~day+ptid`).

### Argument File Example
| Argument            | Value                          |
|---------------------|--------------------------------|
| data_file           | Path to count matrix          |
| sample_meta_data    | Path to sample meta data      |
| day_subset          | 0,3                           |
| treatment_gp        | Treatment group details       |
| model design parameters | ~day+pubid                 |
| output_directory    | Path to store results         |

---

## Running MultiDGEA Without Argument Files

### Steps to Run Analysis

1. **Install the package:**
   ```R
   install.packages(file.choose(), repos = NULL, type="source")
   ```

2. **Load the library:**
   ```R
   library(MultiDGEA)
   ```

3. **Load input data:**
   ```R
   counts_data <- read.csv("counts.csv", row.names = 1)
   sample_meta_data <- read.csv("sample_meta_data.csv", row.names = 1)
   ```

4. **Run analyses:**
   - **Kimma:**
     ```R
     Kimma_result <- fitting_data_kimma(count_matrix = counts_data,
                                        meta_data = sample_meta_data,
                                        model_design_parameter = "~day",
                                        kimma_model = "~ day + (1|pubid)",
                                        contrast_var = "day",
                                        patient_id_var = "pubid")
     ```
   - **DREM:**
     ```R
     Drem_results <- fitting_data_DREAM(count_matrix = counts_data,
                                         meta_data = sample_meta_data,
                                         model_design_parameter = "~day",
                                         form = "~ day + (1|pubid)",
                                         contrast_var = "day",
                                         patient_id_var = "pubid")
     ```
   - **Limma Voom:**
     ```R
     Limma_voom_results <- fitting_data_limma_voom(count_matrix = counts_data,
                                                   meta_data = sample_meta_data,
                                                   model_design_parameter = "~day",
                                                   contrast_var = "day",
                                                   patient_id_var = "pubid")
     ```
   - **Limma Voom Duplicate Correlation:**
     ```R
     Limma_voom_duplicate_corr_results <- limma_voom_duplicate_correlation(count_matrix = counts_data,
                                                                           meta_data = sample_meta_data,
                                                                           model_design_parameter = "~day+pubid",
                                                                           contrast_var = "day",
                                                                           patient_id_var = "pubid")
     ```
   - **DESeq:**
     ```R
     Deseq_result <- fitting_data_Deseq(count_matrix = counts_data,
                                        meta_data = sample_meta_data,
                                        model_design_parameter = "~day+pubid",
                                        contrast_var = "day",
                                        patient_id_var = "pubid")
     ```

---

## Visualization of Results

1. **Normalized Counts Plot:**
   - For Limma Voom:
     ```R
     Normalized_data_plot(method = "limma_voom",
                          result = limma_voom_results,
                          count_matrix = counts_data,
                          meta_data = sample_meta_data,
                          model_design_parameter = "~day")
     ```
   - For DESeq:
     ```R
     Normalized_data_plot(method = "deseq",
                          result = deseq_results,
                          count_matrix = counts_data,
                          meta_data = sample_meta_data,
                          model_design_parameter = "~day")
     ```

2. **Heatmap:**
   ```R
   heatmap_top_100_MultiDGEA_edgeR(result = Limma_voom_results,
                                   count_matrix = counts_data,
                                   meta_data = sample_meta_data,
                                   model_design_parameter = "~day")
   ```

3. **Histogram:**
   ```R
   histogram_plot_p_value(Deseq_results)
   ```

4. **Volcano Plot:**
   ```R
   volcano_plot(Limma_voom_results, "volcano_plot_for_Limma_voom_duplicate_correlation")
   ```

---

## Comparison of Results

1. **P-Value Scatter Plot:**
   - For downregulated genes:
     ```R
     p_val_scatter_plot_comparison_down(res1 = Limma_voom_duplicate_corr_results,
                                        res2 = Deseq_results,
                                        p_val_xlabel = "p_val_Limma_voom_duplicate_corr_results",
                                        p_val_ylabel = "p_val_Deseq_result",
                                        p_val_title = "Limma_voom_duplicate_corr_vs_Deseq_results",
                                        p_val = 0.05,
                                        FDR_ = 0.2,
                                        LFC_ = 0.05)
     ```
   - For upregulated genes:
     ```R
     p_val_scatter_plot_comparison_up(res1 = Limma_voom_duplicate_corr_results,
                                      res2 = Deseq_results,
                                      p_val_xlabel = "p_val_Limma_voom_duplicate_corr_results",
                                      p_val_ylabel = "p_val_Deseq_result",
                                      p_val_title = "Limma_voom_duplicate_corr_vs_Deseq_results",
                                      p_val = 0.05,
                                      FDR_ = 0.2,
                                      LFC_ = 0.05)
     ```

2. **LFC Comparison:**
   - For downregulated genes:
     ```R
     LFC_scatter_plot_comparison_down(res1 = Limma_voom_duplicate_corr_results,
                                       res2 = Deseq_result,
                                       LFC_xlabel = "LFC_Limma_voom_duplicate_corr_results",
                                       LFC_ylabel = "LFC_Deseq_result",
                                       LFC_title = "Limma_voom_duplicate_corr_vs_Deseq_results",
                                       p_val = 0.05,
                                       FDR_ = 0.2,
                                       LFC_ = 0.05)
     ```
   - For upregulated genes:
     ```R
     LFC_scatter_plot_comparison_up(res1 = Limma_voom_duplicate_corr_results,
                                     res2 = Deseq_results,
                                     LFC_xlabel = "LFC_Limma_voom_duplicate_corr_results",
                                     LFC_ylabel = "LFC_Deseq_result",
                                     LFC_title = "Limma_voom_duplicate_corr_vs_Deseq_results",
                                     p_val = 0.05,
                                     FDR_ = 0.2,
                                     LFC_ = 0.05)
     ```

3. **Venn Diagram:**
   - Finding significant genes:
     ```R
     Limma_voom_duplicate_corr_results_MultiDGEA <- find_MultiDGEA_genes(result = Limma_voom_duplicate_corr_results,
                                                                         p_val = 0.05,
                                                                         FDR_ = 0.2,
                                                                         LFC_ = 0.05)
     Deseq_results_MultiDGEA <- find_MultiDGEA_genes(result = Deseq_results,
                                                     p_val = 0.05,
                                                     FDR_ = 0.2,
                                                     LFC_ = 0.05)
     sig1 <- Deseq_results_MultiDGEA$up_regulated
     sig2 <- Limma_voom_duplicate_corr_results_MultiDGEA$up_regulated
     plot_venn(sig1, sig2, "Limma_voom_vs_Deseq")
     ```

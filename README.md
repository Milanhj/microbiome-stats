
# RCT AMR and microbiome analysis

Collected function documentation for comparing differences in microbiomes and resistomes between treatment arms.
Example code uses simulated count data. All scripts include commented functions with documentation. 

<br>

# DESeq2 with FDR

### Description

Using `DESeq2` for differential expression analysis.

**Usage**

`do_deseq(dat, stop_col, formula, alpha = 0.05, test = "Wald", sf_type = "custom", total_counts = NULL, ordered = FALSE, reduced = NULL)`:

**Arguments**

dat

-   Dataset with IDs, metadata for design matrix, and counts of all integer values.
requires an ID variable in the first column and all meta data before first count column

stop_col

-   Index for final column of metadata, must be directly before the count data starts

alpha

-   Significance cutoff used for optimizing the independent filtering. default = 0.05, but default is 0.1 for the package function. set to adjusted p-value cutoff

test

-   Either "Wald" for cross-sectional, or "LRT" for longitudinal

sf_type

- Indicates which method to use for estimating size factors. May be set to "median_ratio" for use of the DESeq2
package default method, or "custom" which estimates size factors with a vector of the total non-human reads in each sample
(default = "custom")

total_counts

- Optional vector of total counts for custom size factor estimation. Must input a vector of total reads if sf_type = "custom"

reduced

- for LRT, reduced formula without the final interaction term

ordered

- set to TRUE if the design formula contains an ordered categorical variable




<br>

### Description

Display DESeq output for a single independent variable. 
Returns a table with the log2foldchange, foldchange, standard error, and p-value
for each level of the dependent variable.  


**Usage**

`display_deseq_results(mod, vars = NULL, var_label = NULL, digits = NULL, p.adjust.method = NULL, na.rm = FALSE)`

**Arguments**

mod

-   Model output from my `do_deseq()` function or package `results(DESeq())`

vars

-   Optional vector with the names of which variables to include in the output table. 
For visualization purposes only. 

var_label

-   Optional string. What to name the column containing vars as rows (defaults to "outcome")

digits

- optional argument how many digits to round output values to. default = 4

p.adjust.method

- optional argument for p.adjust method, if NULL, no values are returned. 

na.rm

- option to remove NA rows. Do not set na.rm = TRUE if there are multiple independent variables in the design matrix 

<br>


### Description

Function wrapping `do_deseq()` and `display_deseq_results()` together.
Outputs a list of tables for each variable level (class, gene, etc..) with results for each variable. 
Used to output DESeq2 results for each variable in a design formula with one or more independent variables.


**Usage**

`fit_deseq2(dat, stop_col, formulas, alpha = 0.05, test = "Wald", sf_type = "custom", total_counts = NULL, ordered = FALSE, reduced = NULL, na.rm = FALSE, p.adjust.method = NULL, vars = NULL, var_label = NULL, digits = NULL)`

**Arguments**

Requires same arguments from `do_deseq()` and `display_deseq_results()`.

dat

-   Dataset with IDs, metadata for design matrix, and counts of all integer values.
requires an ID variable in the first column and all meta data before first count column

formulas

- a list of formulas to iterate through. 
DESeq2 returns results for the final variable in the formula. 
Must input formulas with condition of interest in the final position for each independent variable.

- `formulas = list(c(~ x1 + x2 + x3), c(~ x1 + x3 + x2), c(~ x2 + x3 + x1))`
- `formulas = list(c(~x1)`


p.adjust.method

- method for p.adjust(). If NULL, no values are returned


<br>

### Description

Function uses `fit_deseq2()` and a permutation loop to return
Westfall-Young minP adjusted p-values for all classes or a single class.


**Usage**

`deseq2_min_p(dat, B = 999, stop_col, formulas, tot_counts_col, time_col = NULL, alpha = 0.05, test = "Wald", sf_type = "custom", sep_time = TRUE, select_var = NULL, ordered = FALSE, reduced = NULL, digits = NULL)`

**Arguments**

Requires same arguments from `do_deseq()` and `display_deseq_results()`.

dat

-   Dataset with IDs, metadata for design matrix, and counts of all integer values.
requires an ID variable in the first column and all meta data before first count column


B

- number of permutations

stop_col

- index of final column of metadata

tot_reads_col 

- column index with total reads for computing size factors

group_col 

- permutation column index. will be renamed `tx`

select_var

- option to run min-p on a single drug class name (select_var = "MLS")
- must input a time column index when running min-P on a single variable, 
because min-p must be run at separate timepoints

sep_time

- option to obtain results for separate timepoints, but run all p-values through a single min-p procedure
- if sep_time = TRUE, you must provide the index for the time variable `time_col`
- If the model has a longitudinal component (i.e. ~ tx + time), sep_time = FALSE


formulas

- a list of formulas to iterate through. 
DESeq2 returns results for the final variable in the formula. 
Must input formulas with condition of interest in the final position for each independent variable.

- `formulas = list(c(~ x1 + x2 + x3), c(~ x1 + x3 + x2), c(~ x2 + x3 + x1))`
- `formulas = list(c(~x1)`

<br>

### Description

Function uses `fit_deseq2()` and a permutation loop to return
Westfall-Young minP adjusted p-values for all classes. Option to run 
select classes separate from all others in a longitudinal data set. 


**Usage**

`minp_deseq2_longitudinal(dat, B = 999, formulas, stop_col, tot_counts_col, time_col, group_col, alpha = 0.05, test = "Wald", select_vars = NULL, thresh = NULL, ordered = FALSE, reduced = NULL, digits = NULL)`

**Arguments**

Requires same arguments from `do_deseq()` and `display_deseq_results()`.

dat

-   Dataset with IDs, metadata for design matrix, and counts of all integer values.
requires an ID variable in the first column and all meta data before first count column

formulas

- a list of formulas to iterate through. 
DESeq2 returns results for the final variable in the formula. 
Must input formulas with condition of interest in the final position for each independent variable.

- `formulas = list(c(~ x1 + x2 + x3), c(~ x1 + x3 + x2), c(~ x2 + x3 + x1))`
- `formulas = list(c(~x1)`

B

- number of permutations

stop_col

- index of final column of metadata

tot_reads_col 

- column index with total reads for computing size factors

group_col 

- permutation column index. will be renamed `tx`

select_vars

- option to run seperate min-p on select variables (select_var = c("MLS", "betalactams"))

thresh

- threshold for proportion of non-zero values needed for removing a count column


<br>


# T-Tests

### Description

Pairwise comparison of treatment arms.
Option to log2 transform reads per million with $Log_2(rM + {\frac{1}{total reads}} \cdot 1,000,000)$

**Usage**

`ttest_log2rm(dat, group_col, tot_reads_col, start_col, time_col, list = FALSE, log2 = TRUE)`

**Arguments**

dat

- dataset with counts. All meta-data must precede count columns

start_col 

- the index of first column with counts

group_col 

- index of grouping column (tx)

total_reads_col 

- column with total non-human read numbers for log transforming

time_col 

- time (longitudinal) or measurment month/year

log2 

- should the data be log2 transformed?

list 

- list = TRUE will return results in list format
- list = FALSE will return results as a single table

<br>


### Description

Pairwise comparison of treatment arms 
with permutation based Westfall-Young Min-P adjustment for multiple comparisons.
Wrapping function that calls `ttest_log2rm()`.
Option to log2 transform reads per million with $Log_2(rM + {\frac{1}{total reads}} \cdot 1,000,000)$


**Usage**

`ttest_min_p(dat, B = 999, group_col, tot_reads_col, start_col, time_col, list = FALSE, log2 = TRUE)`

**Arguments**

dat

- dataset with counts. All meta-data must precede count columns

B

- number of permutations. Default is 9999

start_col 

- the index of first column with counts

group_col 

- index of grouping column (tx)

total_reads_col 

- column with total non-human read numbers for log transforming

time_col 

- time (longitudinal) or measurment month/year

log2 

- should the data be log2 transformed

list 

- list = TRUE will return results in list format
- list = FALSE will return results as a single table

<br>


# Wilcoxon Rank Sum


### Description

Wilcoxon rank sum test at a single time. 

**Usage**

`rank_sum(dat)`

**Arguments**

dat

- dataset with just rM or counts and grouping variable (treatment arm)
- grouping variable must be in 1st column, followed by count columns


<br>

### Description

Wilcoxon rank sum test for all timepoints in a study. 
Returns results for all classes in each timepoint in a single table

**Usage**

`rank_sum_longitudinal(dat, time_col, group_col, start_col)`

**Arguments**

dat

- dataset with counts. All meta-data must precede count columns

start_col 

- the index of first column with counts

group_col 

- index of grouping column (tx)

total_reads_col 

- column with total non-human read numbers for log transforming

time_col 

- time (longitudinal) or measurment month/year. 

<br>



### Description

Wilcoxon rank sum test for all timepoints in a study with Westfall-Young Min-p adjustment for multiple comparisons.
Calls `rank_sum_longitudinal()`

**Usage**

`rank_sum_min_p(dat, time_col, group_col, start_col, select_vars = NULL, thresh = NULL,)`

**Arguments**

dat

- dataset with counts. All meta-data must precede count columns

B

- number of permutations

start_col 

- the index of first column with counts

group_col 

- index of grouping column (tx)


time_col 

- time (longitudinal) or measurment month/year. 
If the study is not longitudinal, use seperate datasets for pre and post intervention min-p


select_vars 

- option to min-p select classes separately, 
- vector of class (outcome_var) names to analyze independent of the others

thresh

- the threshold for proportion of non-zero values needed to be included in analysis




<br>

# Diversity

### Description

$$L_{\alpha} = \left( \sum^n_{i \ = \ 1} P_i^{\ \ \alpha} \right) ^{\frac{1}{1 \ - \ \alpha}}$$

Calculate $L_1$ $(\alpha = 1)$ and $L_2$ $(\alpha = 2)$ diversity for a vector of counts or normalized reads. 

Shannon's $L_1$ diversity is calulated as $e^{-\sum^n_{i=1} P_i log(P_i)}$

Simpson's $L_2$ diversity returns a vector with *[1]* simpson's index, and *[2]* inverse simpson.



**Usage**

`shannon(v)`

`simpson(v)`

**Arguments**

v

-   Vector of counts or normalized read numbers

<br>

### Description

Calculate all diversity indices simultaneously.
Output table with treatment arm, shannon, simpson's index, and inverse simpson.

**Usage**

`calc_diversity(dat, start_col)`

**Arguments**

dat

-   Count data with treatment arm in column immediately before counts start

start_col

-   Index of first count column

<br>

### Description

Plot diversity values from `calc_diversity()` output, requires `patchwork` loaded.

**Usage**

`plot_diversity(dat, xlim, timepoint, groups, colors)`

**Arguments**

dat

-   Output table from `calc_diversity()`

xlim

-   Manually set limits of x-axis (if you don't do this, values get truncated)

timepoint

-   Character string for plot title (e.g. "baseline", "36 months"..)

groups

-   Character vector for naming groups in the legend

colors

-   Vector of colors matching groups

<br>

# PERMANOVA

### Description

Runs PERMANOVA on observed data or principal components from PCoA at all timepoints in a dataset.

**Usage**

`pcoa_permanova(dat, count_start_col, group_col, time_col, run_pcoa = TRUE, B = 9999, method, thresh_variance = 0.8)`

**Arguments**

dat

-   Count data with time and group variables

count_start_col

-   First column with counts

run_pcoa

- logical input. 
- if `run_pcoa == TRUE`, PCoA is computed and principal components are used for PERMANOVA 
- if `run_pcoa == FALSE`, the observed data is used for PERMANOVA

method

-   Distance measure for `vegan::adonis2()`
-   `"l1"` for manhattan or `"l2"` for euclidean

B

-   Number of permutations

thresh_variance

- threshold for determining how many principal components to use for PERMANOVA based on percent variance explained

<br>

# PCoA

### Description

Function returning Eigenvectors, Eigenvalues, and Centroids from PCoA.

**Usage**

`pcoa(dat, start_col, method)`

**Arguments**

dat

-   Count data with time variable and group column preceding count columns

start_col

-   First column with counts

method

-   Distance measure for vegan distance matrix (takes "l1" or "l2")

<br>


### Description

Plot Centroids from PCoA1 and PCoA2.

**Usage**

`plot_pcoa(vectors, centroids, label, colors = NULL, title_text = NULL, subtitle_text = NULL)`

**Arguments**

vectors

-   Eigenvectors from `pcoa()` output

centroids

-   Centroids from `pcoa()` output

label

-   Labels for groups in the legend

colors

-   Option to manually set group colors

title_text

-   Optional plot title

subtitle_text

-   Optional plot subtitle

<br>


# Correlation


### Description

Test for correlation/association between paired samples and correct for multiple comparisons.
Returns correlation coeffiecents, p-values, padjust values, a summary table with all values, 
and a vector of columns (if any) that were filtered out due to sparseness.


**Usage**

`corr_fdr(dat, method = "spearman", p.adjust.method = "fdr", nonzero_thres = round(0.1*nrow(dat)), var_lab = NULL, digits = NULL)`


**Arguments**

dat

- Dataframe with only count or normalized read data

method

- Which correlation coefficient to use (default = spearman)

p.adjust.method

- Method for p.adjust correction for multiple comparisons (default = "fdr")

nonzero_thres

- Threshold of non-zero values to use for filtering out columns with mostly zero counts.
Default removes columns with fewer than 10% non-zero values

var_lab

- Optional, what column names should be in result table (e.g "gene","class", etc). Defaults to "var"

digits

- Optional argument, how many digits to round final outputs to if desired


<br>





# Longitudinal Results


### Description

Combine all output values from different timepoints into a single table.

**Usage**

`combine_table(list_results, var_col, var_name, bind = FALSE)`

**Arguments**

grp_name

-   The name of the test group to pull out (for example, "MLS")


<br>




### Description

Combine p-value outputs for different timepoints in a single table.

**Usage**

`combine_pval_table(list_results, var_col, var_name, bind = FALSE)`

**Arguments**

list_results

-   `fit_permanova()` output for each timepoint collected in a list

var_col

-   The column of the result variable to be represented (e.g. column with `padj`, etc)

var_name

-   Name of the variable in var_col

time_name

- Manually set the name of timepoint labels appended at the end of var_name. 
Default is numerical order beginning at 0 (e.g. 0, 1, ... n - 1)

bind

-   If no key variable for joining exists, use bind = TRUE to bind columns instead of joining

<br>

### Description

Plot a 3D scatterplot for visualizing time centroids on 3 principal components.

**Usage**

`plotly_scatterplot(centroids, colors, labels, title_text = NULL)`

**Arguments**

centroids

-   Centroid coordinate output from `pcoa()`

colors

-   Vector of colors matching number of timepoints

labels

-   Text to label reach centroid with timepoint

title_text

-   Add an optional title

<br>


# RCT microbiome analysis

Collected function documentation for comparing differences in microbiomes and resistomes of treatment arms.
Example code uses simulated count data. All scripts include commented functions with documentation. 

<br>

# DESeq2 with FDR

### Description

Using `DESeq2` for differential expression analysis.

**Usage**

`do_deseq(dat, stop_col, formula, alpha = 0.05, test = "Wald", sf_type, total_counts = NULL, ordered = FALSE, reduced = NULL)`:

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

`display_deseq_results(mod, vars, var_label, digits = 4, fdr = FALSE, na.rm = FALSE)`

**Arguments**

mod

-   Model output from my `do_deseq()` function or package `results(DESeq())`

vars

-   The names of which variables to include in the output table

var_label

-   What to name the column containing vars as rows

digits

- how many digits to round output values to. default = 4

fdr

- if TRUE, q-values are calculated and returned

na.rm

- option to remove NA rows. Do not set na.rm = TRUE if there are multiple independent variables in the design matrix 

<br>


### Description

Function wrapping `do_deseq()` and `display_deseq_results()` together.
Outputs a list of tables for each variable level (class, gene, etc..) with results for each variable. 
Used to output DESeq2 results for each variable in a design formula with multiple independent variables.


**Usage**

`fit_deseq2(dat, stop_col, formulas, alpha = 0.05, test = "Wald", sf_type, total_counts = NULL, ordered = FALSE, reduced = NULL, na.rm = FALSE, fdr = FALSE, vars, var_label, digits = 4)`

**Arguments**

Same arguments from `do_deseq()` and `display_deseq_results()`.

formulas

- a list of formulas to iterate through. 
DESeq2 returns results for the final variable in the formula. 
Must input formulas with condition of interest in the final position for each independent variable.

- formulas = list(c(~ x1 + x2 + x3), c(~ x1 + x3 + x2), c(~ x2 + x3 + x1))


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

Combines runs PERMANOVA, comparing groups at a single timepoint.

**Usage**

`fit_permanova(dat, start_col, mthd, perm = 9999)`

**Arguments**

dat

-   Count data with time variable and group column preceding count columns

start_col

-   First column with counts

mthd

-   Distance measure for `vegan::adonis2()`

perm

-   Number of permutations

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

# Longitudinal Results

### Description

Combine p-value outputs for different timepoints in a single table.

**Usage**

`combine_table(list_results, var_col, var_name, bind = FALSE)`

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

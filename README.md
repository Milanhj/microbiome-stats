
# RCT microbiome analysis

Collected function documentation for comparing differences in microbiomes and resistomes of treatment arms.
Example code uses simulated count data. All scripts include commented functions with documentation. 

<br>

# DESeq2 with FDR

### Description

Using `DESeq2` for differential expression analysis.

**Usage**

`do_deseq(df, tx_col, alpha = 0.05, test = "Wald", time_col = NULL)`:

**Arguments**

df

-   Dataset with atleast IDs, treatment, and counts of all integer values.
requires an ID variable in the first column and treatment variable directly before first count column

tx_col

-   Treatment column index, must be directly before the count data starts

alpha

-   Significance cutoff used for optimizing the independent filtering. default = 0.05, but default is 0.1 for the package function. set to adjusted p-value cutoff

test

-   Either "Wald" for cross-sectional, or "LRT" for longitudinal

time_col

-   Optional column index of timepoint data (must be indexed between ID and treatment variables)

<br>

### Description

Display DESeq output in a cleaner, more easy to work with way.


**Usage**

`display_deseq_results(mod, vars, var_label)`

**Arguments**

mod

-   Model output from my `do_deseq()` function or package `results(DESeq())`

vars

-   The names of which variables to include in the output table

var_label

-   What to name the column containing vars as rows

<br>

# Diversity

### Description

$$L_{\alpha} = \left( \sum^n_{i=1} P_i^ \alpha \right) ^{\frac{1}{1-\alpha}}$$

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

`calc_diversity(df, start_col)`

**Arguments**

df

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

`fit_permanova(df, start_col, mthd, perm = 9999)`

**Arguments**

df

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

`pcoa(df, start_col, method)`

**Arguments**

df

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

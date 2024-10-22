[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/7ZSnxz5Q)
# module1-biomarker

This is the template repository for the module 1 group assignment conducting a proteomic analysis of serum data to identify possible early detection biomarkers of autism spectrum disorder (ASD).

Data for this assignment come from Hewitson *et al.* (2021). *Blood biomarker discovery for autism spectrum disorder: A proteomic analysis.* PLoS One. 2021 Feb 24;16(2):e0246581. doi: 10.1371/journal.pone.0246581. Data are publicly available and were accessed on August 19, 2022.

### Repository contents

-   `data` contains

    -   `biomarker-raw.csv` raw data as retrieved from PLoS One

    -   `biomarker-clean.RData` processed data

-   `scripts` contains

    -   `preprocessing.R` script used to generate the processed data

    -   `preliminary-analysis.R` concise version of in-class codes used to replicate published analysis

-   `results` contains a report template `report.qmd`

### Assignment tasks

Most of the lifting for this analysis was already done in class. As noted, the published analysis excluded 192 unidentified variables which are included in our analysis by default. In addition, we made some simplifications for ease of presentation -- specifically, repeated RF fitting, forward selection in the logistic regression model, and repeated estimation of AUROC for the final classifier(s) were omitted. As such, our results will not match those in the paper exactly.

Your group's task is to explore the sensitivity of the results to certain design choices in the methodology. Use the provided scripts `preprocessing.R` and `preliminary-analysis.R` as starting points to address the following questions:

1.  What do you imagine is the reason for log-transforming the protein levels in `biomarker-raw.csv`? (Hint: look at the distribution of raw values for a sample of proteins.)

2.  Temporarily remove the outlier trimming from preprocessing and do some exploratory analysis of outlying values. Are there specific *subjects* (not values) that seem to be outliers? If so, are outliers more frequent in one group or the other? (Hint: consider tabluating the number of outlying values per subject.)

3.  Experiment with the following modifications:

    -   repeat the analysis but carry out the entire selection procedure on a training partition -- in other words, set aside some testing data at the very beginning and don't use it until you are evaluating accuracy at the very end

    -   choose a larger number (more than ten) of top predictive proteins using each selection method

    -   use a fuzzy intersection instead of a hard intersection to combine the sets of top predictive proteins across selection methods

    How are results affected by each modification?

4.  Use any method to find either:

    -   a simpler panel that achieves comparable classification accuracy

    -   an alternative panel that achieves improved classification accuracy

    Benchmark your results against the in-class analysis.

### Assignment instructions

1.  Coordinate with your assigned group, divide tasks, and make a plan for how to proceed.
2.  Similar to last time, develop your scratch work via scripts in the repository. Store these in a subdirectory, *e.g.,* `scripts/drafts/FILENAME.R` .
3.  Pool your work and develop a primary script that contains the codes that produce your findings. Store this as `scripts/analysis-main.R` .
4.  Prepare a write-up using the template `results/report.qmd` .
5.  To submit, commit all changes to the main branch of your group repository by the assignment deadline.

*Remark*: ensure your contributions are recorded in the commit history.

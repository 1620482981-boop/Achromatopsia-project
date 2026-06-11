# Codes for RT-qPCR analysis

## Overview

Python codes instead of Excels are used to analyse and plot the summary plot since they produce better reproducibility
than simple excel-mediated workflows. This README file describes two versions of the code suited for routine qPCR analysis, plotting summary statistics(`delta_delta_ct`, `fold_change`) with respect to the control groups and last but not least,conducting statistical tests based on summary statistic values. 

## Input Files

| Files | Descriptions |
| ----- | ------------ |
| Chengjie_qPCR_Chrnb4.GFP_pde6c.csv | Chrnb4 mice with 3 timepoints using Pde6c probe |
| Pde6c wild type.csv | First run on Pde6c mice with 3 timepoints using Pde6c probe |
| Pde6c wild type rerun.csv | Second run on Pde6c mice with 3 timepoints using Pde6c probe |
| qPCR test_promoter.csv | Run on P2M Chrnb4 mice with testing promoter(BP400, BP700, BP1500) using mcherry probe |
| promoter control -  Quantification Summary.csv | Run on P2M Chrnb4 mice with control promoter(ProA1, PR1.7) using mcherry probe |
| pde promoter total.csv | Run on P2M Chrnb4 mice with all promoters using mcherry probe |
| P12 promoter control.csv | Run on P12 Chrnb4 mice with control promoters using mcherry probe |

## Software Requirement

This project can be run in python 3.

Required python packages

- `pandas`: data cleaning and manipulation
- `numpy`: numerical calculations
- `matplotlib`: plotting
- `seaborn`: statistical visualisation
- `pathlib`: results path creating
- `scipy`: statistical tests
- `statsmodels`: statistical tests

If using Conda, create and activate an environment using:

```bash
conda create -n qpcr_env python=3.10 pandas numpy matplotlib seaborn scipy statsmodels pathlib
conda activate qpcr_env
```

## Folder Structure

- promoter/: promoter .csv files with Ct values from qPCR machine
- wildtype/: wildtype .csv files with Ct values from qPCR machine
- scripts/: python scripts for generating outputs
- results/: summary statistics and plots
- README.md: This README

## How to run the workflow

**Step 1**: Open the computer terminal and type the command below

```bash
cd achromatopsia/
```

This would take you to the project directory

**Step 2**: Install the required packages indicated in session \#\# Softwave development 

```bash
conda create -n qpcr_env python=3.10 pandas numpy matplotlib seaborn scipy statsmodels pathlib
conda activate qpcr_env
```

**Step 3**: run python command in terminal

```bash
python script/[script name]
```
Pair-wise test statistics beteen delta delta ct groups and p value for levene test are shown in terminal after the command has been executed correctly.  

**Step 4**: Check results at results/

## Output files

Output files name have the same prefix as their respective input but with different suffix.

| Format | Suffix | Description |
|------| -----  | ----------- |
| .png | fold change | plots including delta delta ct and fold change |
| .png | linearity check | Check whether the linearity condition meets |
| .png | statistical test | boxplots of delta delta ct values among different groups |
| .csv | summary table | summary statistics |

## Notes

- Users must run the bash command under project directory.
- The results/ subdirectory will only appears after running **step 3**

## Author

Name: Chengjie Peng
Date 11/6/2026










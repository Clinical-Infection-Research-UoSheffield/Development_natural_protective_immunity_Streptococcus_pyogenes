# Early life serological profiles and the development of natural protective humoral immunity to Streptococcus pyogenes in a high burden setting 

Code to reproduce analyses from SpyCATS longitudinal study of immunity to Strep pyogenes in The Gambia

## Description

1. This repository contains the code to reproduce analyses from manuscript entitled: Early life serological profiles and the development of natural protective humoral immunity to Streptococcus pyogenes in a high burden setting  
2. Data including antibody measurements, clinical events, and demographics have been anonymized to protect participant confidentiality.
3. Analyses include mixed-effects logistic regression to model event probabilities and piecewise regression. Cox proportional hazards models are used to estimate risk across different outcomes, with results formatted into comprehensive tables and plots.

## Data anonymisation 

IDs are anonymized by mapping them to randomly generated codes, all dates are uniformly offset by a constant to preserve time intervals while concealing actual dates, and exact ages are replaced with pseudo-ages randomly generated within defined age groups.

Upon publication data to reproduce analyses will be made publically available, hosted on Zenodo at <insert DOI url>

## Repository Structure

```markdown
.
├── README.md
├── data
├── R_output
├── scripts


```

data: Contains processed data files.

R_output: Contains output from R scripts, such as figures and summaries.

scripts: Contains R scripts used for data analysis and visualization.


##  Prerequisites
Ensure you have R installed on your system. You will also need the librarian package to manage dependencies.

## Installing Dependencies
The following R packages are required for the analyses. The librarian package will ensure all dependencies are installed and loaded:

## License

This repository is licensed under the **Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)** license.

You are free to:

- **Share** — Copy and redistribute the material in any medium or format.

### Under the following terms

- **Attribution** — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
- **NonCommercial** — You may not use the material for commercial purposes.
- **NoDerivatives** — If you remix, transform, or build upon the material, you may not distribute the modified version.

For full license details, see: https://creativecommons.org/licenses/by-nc-nd/4.0/


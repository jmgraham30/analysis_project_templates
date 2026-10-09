# Analysis Project Templates

Template R projects for an undergraduate **neuroscience research methods** course. Each project is a small, self-contained example of one type of statistical analysis, built to model good habits for organizing, documenting, and reporting data analysis:

- every analysis lives in its own **RStudio project**,
- raw **data** are kept in a `Data/` folder and never edited by hand,
- **R scripts** are saved and documented in an `R/` folder, and
- the work is written up in a **Quarto lab notebook** that combines the question, methods, code, results, and conclusions in one reproducible document.

## The projects

| Project folder | Analysis topic | Status |
|:--|:--|:--|
| `data_organization_and_display` | Organization and display of data | Planned |
| `measures_tendency_variability` | Measures of central tendency and variability | Planned |
| `correlation_regression` | Correlation and regression | Planned |
| `one_sample_t_test` | One-sample t-test | **Ready for review** |
| `two_sample_t_test` | Two-sample t-test | **Ready for review** |
| `oneway_anova` | One-way analysis of variance | Planned |
| `twoway_anova` | Two-way analysis of variance | Planned |
| `chisquare_test` | Chi-square test | Planned |

Where appropriate, each project's notebook illustrates **null hypothesis testing** (including Type I and Type II errors), **effect size**, and **power**, and shows how to **write up and present results** in a publication-style format.

## What is inside each project

```
project_name/
├── project_name.Rproj      <- double-click this to open the project
├── project_name.qmd        <- the Quarto lab notebook (question, methods, code, results)
├── Data/
│   ├── data_file.csv       <- the data (raw data are never edited by hand)
│   └── README.md           <- describes the data set and every variable
└── R/
    ├── 00_setup.R          <- packages, plot theme, and helper functions
    ├── 01_simulate_data.R  <- how the (simulated) data were created
    └── 02_analysis.R       <- the full, commented analysis script
```

Running the code creates an `output/` folder with saved figures. It is regenerated every time and is not tracked by git.

## What is an RStudio project?

An **RStudio project** is a folder that RStudio treats as one self-contained piece of work, such as one analysis, one lab report, or one study. You recognize it by the `.Rproj` file inside. Opening that file does three useful things:

- **It sets the working directory** to the project folder. R looks for files relative to that folder, so a path like `Data/my_data.csv` works on any computer, with no `C:\Users\yourname\...` paths that break when the work moves.
- **It gives each analysis a clean workspace.** Objects, history, and open files from one project do not leak into another.
- **It keeps everything for one analysis together** (data, code, notebook, and output), so you can find it, share it, or return to it months later.

Think of it like a lab notebook with its own labeled shelf: everything for one experiment is in one place, and nothing gets mixed up with another experiment. Working this way is a core habit of reproducible research. A collaborator, or you in six months, should be able to open the project and run it without guessing where things belong.

**Rule of thumb:** open the `.Rproj` file first, every time. You can check that it worked by looking at the project name in the top-right corner of RStudio, or by running `getwd()` in the Console.

## Getting started (for students)

1. **Install the software** (once):
   - [R](https://cran.r-project.org/) and [RStudio Desktop](https://posit.co/download/rstudio-desktop/). Current RStudio releases include Quarto; otherwise install it from [quarto.org](https://quarto.org/docs/get-started/).
2. **Download this repository.** On GitHub, click the green **Code** button, then **Download ZIP**, and unzip it. (You do not need to use Git.)
3. **Open one project.** Open the folder for your topic and double-click the `.Rproj` file. *Always open the `.Rproj` file first.* This sets the working directory to the project folder so that file paths work.
4. **Open the lab notebook** (the `.qmd` file) and click **Render**. The first time, R will install any missing packages, which needs an internet connection and may take a few minutes.
5. **Explore.** Read the notebook, open the scripts in `R/`, and try the "Your Turn" exercises at the end of the notebook.

### Packages used

The projects use the **tidyverse** (`dplyr`, `ggplot2`, `readr`, `purrr`, and friends) and tidyverse-friendly packages: `here`, `broom`, `effectsize`, `pwr`, `knitr`, `kableExtra`, and `patchwork`. Other packages may be added in later projects (for example `car` and `emmeans` for ANOVA). Each project's `R/00_setup.R` installs anything missing.

## Conventions used throughout

- **Tidyverse style** wherever possible (pipes, `dplyr` verbs, `ggplot2`).
- **Colorblind-friendly palettes** (the Okabe-Ito palette) in all figures.
- **APA-style reporting** of statistics (test statistic, degrees of freedom, exact *p*, effect size with confidence interval).
- **American English** in text, code, and comments.
- **Relative paths** built with the `here` package, so projects run on any computer.
- **Simulated data** are clearly labeled. Their README documents how they were generated and what the "true" values were. Real data sets, when used, are documented with their source.

## Good workflow habits these templates model

1. One project per analysis; open the `.Rproj` file first.
2. Keep raw data raw. Do all cleaning and calculation in code.
3. State the question, hypotheses, and analysis plan (including exclusion rules and the significance level) *before* looking at results.
4. Look at your data (summaries and plots) before you test anything.
5. Check assumptions, then run the test.
6. Report effect sizes and confidence intervals, not just *p*-values.
7. Document your decisions in the notebook log so the analysis can be understood and repeated.
8. Re-run everything from the top (**Render**) before you trust the results.

## For instructors

The `.gitignore` excludes R/RStudio session files (`.Rproj.user/`, `.Rhistory`), Quarto caches, rendered HTML/PDF files, the `output/` folders, and operating system clutter. Rendered notebooks are therefore not stored in the repository; students create them by rendering. If you would like to share rendered examples, consider hosting them with GitHub Pages or posting the HTML files in your course site.

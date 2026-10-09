# Guides

Two companion documents that sit *beside* the analysis projects. Neither needs a data file; open the `.qmd` in RStudio and click **Render**.

| File | Purpose |
|:--|:--|
| `choosing_a_test.qmd` | **Decision guide.** Go from your research question and study design to an analysis: a flowchart, a design-to-test table (effect size, assumptions, R function, and the matching template project), rank-based alternatives, assumption checks, ten practice scenarios, and a list of common mistakes. Includes the "unit of analysis" (pseudoreplication) warning. |
| `analysis_plan_template.qmd` | **Analysis plan template.** Fill it in *before* collecting data: question, hypotheses and error rates, design, variables, sample size and power (edit a few settings and the numbers update), data-handling rules, the planned analysis, interpretation rules, a pre-collection checklist, and a deviations log. Ends with a completed example. |

**Suggested order for an independent project:** read the decision guide, complete the analysis plan, have your instructor review it, collect the data, and then use the matching template project as a model for your lab notebook.

Orange boxes in the decision-guide flowchart (repeated-measures designs and rank-based tests) do not yet have a template project; the guide gives the R functions to use.

These files install a few small packages if they are missing (`tibble`, `knitr`, and `kableExtra`; the plan template also uses `dplyr` and `pwr`).

library(gtsummary)
library(tidyverse)
library(flextable)
# summarize the data with our package
table1 <-
  trial |> 
  tbl_summary(include = c(age, grade, response))


table2 <-
  tbl_summary(
    trial,
    include = c(age, grade, response),
    by = trt, # split table by group
    missing = "no" # don't list missing data separately
  ) |> 
  add_n() |> # add column with total number of non-missing observations
  add_p() |> # test for a difference between groups
  modify_header(label = "**Variable**") |> # update the column header
  bold_labels()


data(trial)
trial %>% 
  tbl_summary(
    by = trt, # stratify by treatment group
    statistic = list(all_continuous() ~ "{mean} ({sd})", 
                     all_categorical() ~ "{n} ({p}%)"), 
    digits = all_continuous() ~ 2
  )


as_gt(table2) %>% gt::gtsave("summary_table.html")
as_flex_table(table2) %>% flextable::save_as_docx(path = "summary_table.docx")

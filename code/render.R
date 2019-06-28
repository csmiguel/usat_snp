#render article
rmarkdown::render("code/rmarkdown/article.Rmd",
  output_file = "usat_snp.docx", output_dir = "article")

rmarkdown::render("code/rmarkdown/tables.Rmd",
  output_file = "Tables.docx", output_dir = "article")

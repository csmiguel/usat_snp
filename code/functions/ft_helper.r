#helper functions flextable

## helper functions
# 1. ft2ggplot function
ft2ggplot <- function(x){
  library(ggplot2)
  ggplot() +
    theme_void() +
    annotation_custom(grid::rasterGrob(x),
                      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}

# 2. export flextable to Word table
ft2word <- function(ft = NULL, ftcaption = NULL, name = "test.docx", ...){
  #ft, is a flextable
  #name, is the export path
  #ftcaption, is a table caption
  h <- officer::read_docx(...)
  if (!is.null(ftcaption))
    h <- officer::body_add_par(h, ftcaption)
  h <- flextable::body_add_flextable(h, value = ft)
  print(h, target = name)
}

#' ---
#' title: "test markdown"
#' author: "luke durell"
#' date: "10/1/2020"
#' output: github_document
#' ---
#' 
#' Hello this is from Rmarkdown.

#' r
#' - here's some code that does something
r <- 5
r + 5
for (i in 4:(r+4)) {
  message(paste0("This is the ", i+4, "th iteration!"))
}


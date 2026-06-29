## Make RStudio's bundled pandoc available outside the RStudio IDE.
rstudio_pandoc <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
if (dir.exists(rstudio_pandoc) && !nzchar(Sys.getenv("RSTUDIO_PANDOC"))) {
  Sys.setenv(RSTUDIO_PANDOC = rstudio_pandoc)
}

## This makes sure that R loads the workflowr package automatically when the
## project is loaded, without making all R startup fail if a dependency is
## temporarily broken.
if (requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  tryCatch(
    library("workflowr"),
    error = function(e) {
      message("workflowr could not be loaded: ", conditionMessage(e))
    }
  )
} else {
  message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
}

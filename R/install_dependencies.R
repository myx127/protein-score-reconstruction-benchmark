# Run in an isolated R library; no participant data are accessed.
options(repos=c(CRAN="https://cloud.r-project.org"))
if(getRversion()<"4.1.0")stop("R >= 4.1.0 required; manuscript used R 4.4.3")
install_if_missing<-function(pkg){if(!requireNamespace(pkg,quietly=TRUE))install.packages(pkg)}
install_if_missing("BiocManager");install_if_missing("remotes")
bioc_packages<-c("impute","preprocessCore","GO.db","AnnotationDbi")
missing_bioc<-bioc_packages[!vapply(bioc_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_bioc)) {
  # Select the Bioconductor release compatible with the caller's installed R.
  BiocManager::install(missing_bioc,ask=FALSE,update=FALSE)
}
versions<-c(WGCNA="1.73",data.table="1.18.2.1",jsonlite="2.0.0",microbenchmark="1.5.0")
for(pkg in names(versions)) {
  if(!requireNamespace(pkg,quietly=TRUE)||as.character(packageVersion(pkg))!=versions[[pkg]]) {
    remotes::install_version(pkg,version=versions[[pkg]],upgrade="never",dependencies=NA)
  }
}
cat("R dependencies installed. CRNP source is loaded directly from vendor/crnp/R.\n")
cat("Manuscript environment: R 4.4.3, WGCNA 1.73; see environment/manuscript_versions.json.\n")

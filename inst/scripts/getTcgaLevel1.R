#Run this script from the directory where you want to store the downloaded files
#Get all Level 1 files for a tumor type from TCGA portal open https directory
#You must know how many batches there are
#You must also know the version you want

URL <- "https://tcga-data.nci.nih.gov"
path <- "/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor"
target <- paste0(URL, path)
domain <- "cgcc/jhu-usc.edu/humanmethylation450/methylation"
center <- "jhu-usc.edu"
platform <- "HumanMethylation450"
dir <- "Level_1"
batch <- 1:11 ## or whatever is correct for your tumor type
disease <- "BLCA" ## replace with the name of your tumor type
version <- "4" ## replace with the correct version (could default to "latest")
url <- paste(target,tolower(disease),  domain, paste(paste(center, disease, sep="_"), platform, dir, 1, version, 0, "tar", "gz", sep="."), sep="/")
cmd <- paste("wget", url, sep=" ")
sapply(cmd, system, USE.NAMES=F)  ## Moiz this is a little horrifying ;-)

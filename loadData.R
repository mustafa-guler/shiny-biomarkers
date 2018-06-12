loadTCGA <- function (root_dir, force = FALSE) {
    # remove trailing slash
    root_dir <- .cleanRootDir(root_dir)
    sample_sheet <- .loadSampleSheet(root_dir)
    if(file_test('-f', file.path(root_dir, "data.rds")) && !force) {
        result <- readRDS(file.path(root_dir, "data.rds"))
        return(result)
    }
    
    # column names are case ids from the summary file
    columns <- dir(root_dir)
    columns <- columns[columns != "summary"]
    columns <- columns[!grepl("\\.", columns)]
    columns <- sample_sheet[match(columns, sample_sheet$File.ID),]$Sample.ID
    columns <- c("gene_id", columns)

    result <- data.frame()
    files <- dir(root_dir)
    num_files <- length(files)
    
    # add each sample to data frame 
    for(i in files) {
        if(i == "summary" || file_test("-f", file.path(root_dir, i))) {
            next
        }
        
        f <- file.path(root_dir, i)
        f <- file.path(f, dir(f, pattern = ".*FPKM-UQ.*"))
        curr <- data.frame(read.delim(gzfile(f), header = FALSE, stringsAsFactors = FALSE))
        
        new_name <- sample_sheet[sample_sheet$File.ID == i,]$Sample.ID
        
        # initialize data frame if it doesn't exist yet
        if(nrow(result) == 0) {
            result <- data.frame(matrix(ncol = length(columns), nrow = nrow(curr)))
            colnames(result) <- columns
            result$gene_id <- curr[1]
        }
        
        result[new_name] <- curr[2]
        # incProgress(1 / num_files)
    }
    saveRDS(result, file=file.path(root_dir, "data.rds"))
    return(result)
}

loadClinical <- function(root_dir) {
    root_dir <- .cleanRootDir(root_dir)
    clinical_file <- file.path(root_dir, "summary", "clinical.tsv")
    clinical <- read.delim(clinical_file, stringsAsFactors = FALSE)
    clinical <- clinical[clinical$age_at_diagnosis != "--",]
    clinical <- clinical[clinical$year_of_birth != '--',]
    clinical$year_of_birth <- as.numeric(clinical$year_of_birth)
    clinical$age_at_diagnosis <- as.numeric(clinical$age_at_diagnosis)
    filt <- clinical$days_to_death == '--'
    clinical$days_to_death[filt] <- 2016.5 - clinical$year_of_birth[filt] - clinical$age_at_diagnosis[filt] / 365
    clinical$days_to_death <- as.numeric(clinical$days_to_death)
    clinical$days_to_death[!filt] <- clinical$days_to_death[!filt] / 365
    names(clinical)[names(clinical) == 'days_to_death'] <- 'years_to_death'
    return(data.frame(clinical))
}

.loadSampleSheet <- function(root_dir) {
    sample_file <- file.path(root_dir, "summary")
    sample_file <- file.path(sample_file, dir(sample_file, pattern = "gdc.*tsv"))
    return(data.frame(read.delim(sample_file, as.is = TRUE)))
}

.cleanRootDir <- function(root_dir) {
    return(sub(paste(.Platform$file.sep, "$", sep = ""), "", root_dir))
}
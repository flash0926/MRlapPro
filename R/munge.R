.munge <- function (data, hm3, trait.names = NULL, N = NULL, info.filter = 0.9,
          maf.filter = 0.01, log.name = NULL, column.names = list(),
          parallel = FALSE, cores = NULL, overwrite = TRUE)
{
  # Annotate useless code
  # if (is.list(files)) {
  #   wrn <- paste0("DeprecationWarning: In future versions a list of filenames will no longer be accepted.\n",
  #                 "                    Please change files to a vector to ensure future compatibility.")
  #   warning(wrn)
  #   files_ <- c()
  #   for (i in 1:length(files)) {
  #     files_ <- c(files_, files[[i]])
  #   }
  #   files <- files_
  # }

  if (is.null(N)) {
    # Annotate useless code
    #N <- rep(NA, length(files))
    N = 1 #only one file
  }
  # Annotate useless code
  # .check_file_exists(files)
  .check_file_exists(hm3)
  # Annotate useless code
  # .check_equal_length(files, trait.names)
  # .check_equal_length(files, N)
  .check_range(N, min = 0, max = Inf, allowNA = TRUE)
  .check_range(info.filter)
  .check_range(maf.filter)
  if (any(!(names(column.names) %in% c("SNP", "A1", "A2",
                                       "effect", "INFO", "P", "N", "MAF", "Z")))) {
    stop(paste0("Names in column.names not recognized. Please use the following keys:\n        ",
                paste(c("SNP", "A1", "A2", "effect", "INFO", "P",
                        "N", "MAF", "Z"), collapse = ", ")))
  }
  .check_boolean(parallel)
  if (!is.null(cores))
    .check_range(cores, min = 0, max = Inf)
  .check_boolean(overwrite)
  filenames <- as.vector(trait.names)
  if (is.null(log.name)) {
    log2 <- paste(trait.names, collapse = "_")
    if (nchar(log2) > 200) {
      log2 <- substr(log2, 1, 100)
    }
    log.file <- file(paste0(log2, "_munge.log"), open = "wt")
  }
  if (!is.null(log.name)) {
    log.file <- file(paste0(log.name, "_munge.log"), open = "wt")
  }
  # Annotate useless code
  # if (parallel & (length(files) == 1)) {
  #   .LOG("Parallel munging requested for a single file.\nParallel munging only has benefits for munging multiple files.\nParallel disabled",
  #        file = log.file)
  #   parallel <- FALSE
  # }
  begin.time <- Sys.time()
  .LOG("The munging of ", length(trait.names), " summary statistics started at ",
       begin.time, file = log.file)
  if (overwrite) {
    existing_files <- c()
    for (trait.name in trait.names) {
      if (file.exists(paste0(trait.name, ".sumstats")))
        existing_files <- c(existing_files, paste0(trait.name,
                                                   ".sumstats"))
    }
    if (length(existing_files) > 0)
      .LOG("File(s) ", paste0(existing_files, collapse = ", "),
           " already exist and will be overwritten", file = log.file)
  }
  .LOG("Reading in reference file", file = log.file)
  ref <- data.table::fread(hm3, header = T, data.table = F)
  if (!parallel) {
    .LOG("Reading summary statistics for ", paste(trait.names,
                                                  collapse = " "), ". Please note that this step usually takes a few minutes due to the size of summary statistic files.",
         file = log.file)
    # Annotate useless code
    # files <- lapply(files, data.table::fread, header = T, quote = "\"",
    #                 fill = T, na.string = c(".", NA, "NA", ""))
    .LOG("All files loaded into R!", file = log.file)
    # Annotate useless code
    # for (i in 1:length(files)) {
    #   .munge_main(i, NULL, files[[i]], filenames[i], trait.names[i],
    #               N[i], ref, hm3, info.filter, maf.filter, column.names,
    #               overwrite, log.file)
    # }
    #add code: only need one data
    .munge_main(1, NULL, data, filenames[1], trait.names[1],
                N[1], ref, hm3, info.filter, maf.filter, column.names,
                overwrite, log.file)
  }
  else {
    # Annotate useless code
    # if (is.null(cores)) {
    #   int <- detectCores() - 1
    # }
    # else {
    #   int <- cores
    # }
    # if (int > length(filenames)) {
    #   .LOG("Number of requested cores(", int, ") greater than the number of files (",
    #        length(filenames), "). Deferring to the lowest number",
    #        file = log.file)
    #   int <- length(filenames)
    # }
    # cl <- makeCluster(int, type = "PSOCK")
    # registerDoParallel(cl)
    # on.exit(stopCluster(cl))
    # utilfuncs <- list()
    # utilfuncs[[".get_renamed_colnames"]] <- .get_renamed_colnames
    # utilfuncs[[".LOG"]] <- .LOG
    # utilfuncs[["gzip"]] <- gzip
    # .LOG("As parallel munging was requested, logs of each sumstats file will be saved separately",
    #      file = log.file)
    # foreach(i = 1:length(filenames), .export = c(".munge_main"),
    #         .packages = c("stringr")) %dopar% {
    #           .munge_main(i, utilfuncs, NULL, filenames[i], trait.names[i],
    #                       N[i], ref, hm3, info.filter, maf.filter, column.names,
    #                       overwrite, NULL)
    #         }
  }
  end.time <- Sys.time()
  total.time <- difftime(time1 = end.time, time2 = begin.time,
                         units = "sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time - mins * 60
  .LOG("     ", file = log.file)
  .LOG("Munging was completed at ", end.time, file = log.file)
  .LOG("The munging of all files took ", mins, " minutes and ",
       secs, " seconds", file = log.file)
  .LOG("Please check the .log file(s) to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files",
       file = log.file)
  flush(log.file)
  close(log.file)
}



.check_file_exists <- function (path, name = NULL)
{
  if (is.null(name))
    name <- deparse(substitute(path))
  if (length(path) > 1) {
    for (x in path) {
      .check_file_exists(x, name = name)
    }
  }
  else {
    if (!(file.exists(path))) {
      stop(paste("File", path, "passed to", name, "does not exist"),
           call. = FALSE)
    }
  }
}


.check_range <- function (val, min = 0, max = 1, allowNA = FALSE, inclusive = FALSE,
          name = NULL)
{
  if (is.null(name))
    name <- deparse(substitute(val))
  if (length(val) > 1) {
    for (x in val) {
      .check_range(x, min = min, max = max, allowNA = allowNA,
                   inclusive = inclusive, name = name)
    }
  }
  else {
    if (is.na(val)) {
      if (!(allowNA))
        stop(paste("Value(s) of", name, "should not be NA"),
             call. = FALSE)
    }
    else {
      if (!(is.numeric(val)))
        stop(paste("Value(s) of", name, "should be numeric"),
             call. = FALSE)
      if (inclusive) {
        if (val <= min) {
          stop(paste("Value(s) of", name, "should be above",
                     min), call. = FALSE)
        }
        else if (val >= max) {
          stop(paste("Value(s) of", name, "should be below",
                     max), call. = FALSE)
        }
      }
      else {
        if (val < min) {
          stop(paste("Value(s) of", name, "should be above",
                     min), call. = FALSE)
        }
        else if (val > max) {
          stop(paste("Value(s) of", name, "should be below",
                     max), call. = FALSE)
        }
      }
    }
  }
}

.check_boolean <- function (val, name = NULL)
{
  if (is.null(name))
    name <- deparse(substitute(val))
  if (length(val) > 1) {
    for (x in val) {
      .check_boolean(x, name = name)
    }
  }
  else {
    if (!(val %in% c(TRUE, FALSE))) {
      stop(paste(name, "should be TRUE or FALSE"))
    }
  }
}


#change all file var to data
#only need data to operate
#no need to read and write file again an again
.munge_main <- function (i, utilfuncs, data, filename, trait.name, N, ref,
          hm3, info.filter, maf.filter, column.names, overwrite, log.file = NULL)
{
  if (is.null(log.file)) {
    log.file <- file(paste0(trait.name, "_munge.log"), open = "wt")
    on.exit(flush(log.file))
    on.exit(close(log.file))
  }
  else {
    .LOG("\n\n", file = log.file, print = FALSE)
  }
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
      assign(j, utilfuncs[[j]], envir = environment())
    }
  }
  if (is.null(data)) {
    data <- data.table::fread(filename, header = T, quote = "\"",
                       fill = T, na.string = c(".", NA, "NA", ""))
  }
  .LOG("Munging file: ", filename, file = log.file, print = TRUE)
  N_provided <- (!is.na(N))
  colnames(data) <- toupper(colnames(data))
  if ("NEFFDIV2" %in% colnames(data)) {
    .LOG("Found an NEFFDIV2 column for sample size. \n\nPlease note that this is likely effective sample size cut in half. The function is automatically doubling this value. This should only be used for liability h^2 conversion for binary traits and that it should reflect the sum of effective sample sizes across cohorts.",
         file = log.file)
    data$NEFF <- file$NEFFDIV2 * 2
    data$NEFFDIV2 <- NULL
  }
  if ("NEFF_HALF" %in% colnames(data)) {
    .LOG("Found an NEFF_HALF column for sample size. \n\nPlease note that this is likely effective sample size cut in half. The function is automatically doubling this value. This should only be used for liability h^2 conversion for binary traits and that it should reflect the sum of effective sample sizes across cohorts.",
         file = log.file)
    data$NEFF <- file$NEFF_HALF * 2
    data$NEFF_HALF <- NULL
  }
  hold_names <- .get_renamed_colnames(toupper(names(data)),
                                      column.names, c("P", "A1", "A2", "effect", "SNP"), filename,
                                      N_provided = N_provided, log.file = log.file, warn_for_missing = c("P",
                                                                                                         "A1", "A2", "effect", "SNP", "N"), utilfuncs = utilfuncs)
  colnames(data) <- hold_names
  if (N_provided) {
    data$N <- N
    .LOG("Using provided N (", N, ") for file:", filename,
         file = log.file)
  }
  if ("MAF" %in% colnames(data)) {
    data$MAF <- ifelse(data$MAF <= 0.5, data$MAF, (1 - data$MAF))
  }
  data$A1 <- factor(toupper(data$A1), c("A", "C", "G", "T"))
  data$A2 <- factor(toupper(data$A2), c("A", "C", "G", "T"))
  .LOG("Merging file:", filename, " with the reference file:",
       hm3, file = log.file)
  b <- nrow(data)
  .LOG(b, " rows present in the full ", filename, " summary statistics file.",
       file = log.file)
  data <- merge(ref, data, by = "SNP", all.x = F, all.y = F)
  .LOG((b - nrow(data)), " rows were removed from the ", filename,
       " summary statistics file as the rs-ids for these rows were not present in the reference file.",
       file = log.file)
  b <- nrow(data)
  if ("P" %in% colnames(data)) {
    data <- subset(data, !(is.na(data$P)))
  }
  if (b - nrow(data) > 0)
    .LOG(b - nrow(data), " rows were removed from the ",
         filename, " summary statistics file due to missing values in the P-value column",
         file = log.file)
  b <- nrow(data)
  if ("effect" %in% colnames(data)) {
    data <- subset(data, !(is.na(data$effect)))
  }
  if (b - nrow(data) > 0)
    .LOG(b - nrow(data), " rows were removed from the ",
         filename, " summary statistics file due to missing values in the effect column",
         file = log.file)
  a1 <- data$effect[[1]]
  data$effect <- ifelse(rep(round(median(data$effect, na.rm = T)) ==
                              1, nrow(data)), log(data$effect), data$effect)
  a2 <- data$effect[[1]]
  if (a1 != a2)
    .LOG("The effect column was determined to be coded as an odds ratio (OR) for the ",
         filename, " summary statistics file. Please ensure this is correct.",
         file = log.file)
  data$effect <- ifelse(data$A1.x != (data$A1.y) & data$A1.x ==
                          (data$A2.y), data$effect * -1, data$effect)
  b <- nrow(data)
  data <- subset(data, !(data$A1.x != (data$A1.y) & data$A1.x !=
                           (data$A2.y)))
  if (b - nrow(data) > 0)
    .LOG(b - nrow(data), " row(s) were removed from the ",
         filename, " summary statistics file due to the effect allele (A1) column not matching A1 or A2 in the reference file.",
         file = log.file)
  b <- nrow(data)
  data <- subset(data, !(data$A2.x != (data$A2.y) & data$A2.x !=
                           (data$A1.y)))
  if (b - nrow(data) > 0)
    .LOG(b - nrow(data), " row(s) were removed from the ",
         filename, " summary statistics file due to the other allele (A2) column not matching A1 or A2 in the reference file.",
         file = log.file)
  if ((sum(data$P > 1) + sum(data$P < 0)) > 100) {
    .LOG("In excess of 100 SNPs have P val above 1 or below 0. The P column may be mislabled!",
         file = log.file)
  }
  data$Z <- sign(data$effect) * sqrt(qchisq(data$P, 1, lower = F))
  if ("INFO" %in% colnames(data)) {
    b <- nrow(data)
    data <- data[data$INFO >= info.filter, ]
    .LOG(b - nrow(data), " rows were removed from the ",
         filename, " summary statistics file due to INFO values below the designated threshold of",
         info.filter, file = log.file)
  }
  else {
    .LOG("No INFO column, cannot filter on INFO, which may influence results",
         file = log.file)
  }
  if ("MAF" %in% colnames(data)) {
    data$MAF <- as.numeric(as.character(data$MAF))
    b <- nrow(data)
    data <- data[data$MAF >= maf.filter, ]
    data <- subset(data, !(is.na(data$MAF)))
    .LOG(b - nrow(data), " rows were removed from the ",
         filename, " summary statistics file due to missing MAF information or MAFs below the designated threshold of",
         maf.filter, file = log.file)
  }
  else {
    .LOG("No MAF column, cannot filter on MAF, which may influence results",
         file = log.file)
  }
  if ("N" %in% colnames(data)) {
    output <- cbind.data.frame(data$SNP, data$N, data$Z,
                               data$A1.x, data$A2.x)
  }
  else {
    output <- cbind.data.frame(data$SNP, N, data$Z, data$A1.x,
                               data$A2.x)
  }
  if (!("N" %in% names(data)) & (exists("N") == FALSE))
    .LOG("Cannot find sample size column for", filename,
         " and a sample size was not provided for the N argument. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N.",
         file = log.file)
  colnames(output) <- c("SNP", "N", "Z", "A1", "A2")
  .LOG(nrow(output), "SNPs are left in the summary statistics file ",
       filename, " after QC.", file = log.file)
  trait.name <- stringr::str_replace_all(trait.name, stringr::fixed(" "), "")
  data.table::fwrite(x = output, file = paste0(trait.name, ".sumstats.gz"),
              sep = "\t", quote = FALSE, row.names = F)
  #gzip(paste0(trait.name, ".sumstats"), overwrite = overwrite)
  .LOG("I am done munging file: ", filename, file = log.file)
  .LOG("The file is saved as ", paste0(trait.name, ".sumstats.gz"),
       " in the current working directory.", file = log.file)
  return()
}



.LOG <- function (..., file, print = TRUE)
{
  msg <- paste0(..., "\n")
  if (print)
    cat(msg)
  cat(msg, file = file, append = TRUE)
}


.get_renamed_colnames <- function (hold_names, userprovided, checkforsingle = c(), filename,
          N_provided, log.file, warnz = FALSE, warn_for_missing = c(),
          stop_on_missing = c(), utilfuncs = NULL)
{
  interpreted_names <- list(SNP = c("SNP", "SNPID", "RSID",
                                    "RS_NUMBER", "RS_NUMBERS", "MARKERNAME", "ID", "PREDICTOR",
                                    "SNP_ID", "VARIANTID", "VARIANT_ID", "RSIDS"), A1 = c("A1",
                                                                                          "ALLELE1", "EFFECT_ALLELE", "INC_ALLELE", "REFERENCE_ALLELE",
                                                                                          "EA", "REF"), A2 = c("A2", "ALLELE2", "ALLELE0", "OTHER_ALLELE",
                                                                                                               "NON_EFFECT_ALLELE", "DEC_ALLELE", "OA", "NEA", "ALT",
                                                                                                               "A0"), effect = c("OR", "B", "BETA", "LOG_ODDS", "EFFECTS",
                                                                                                                                 "EFFECT", "SIGNED_SUMSTAT", "EST", "BETA1", "LOGOR"),
                            INFO = c("INFO", "IMPINFO"), P = c("P", "PVALUE", "PVAL",
                                                               "P_VALUE", "P-VALUE", "P.VALUE", "P_VAL", "GC_PVALUE",
                                                               "WALD_P"), N = c("N", "WEIGHT", "NCOMPLETESAMPLES",
                                                                                "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N", "N_COMPLETE_SAMPLES",
                                                                                "SAMPLESIZE", "NEFF", "N_EFF", "N_EFFECTIVE", "SUMNEFF"),
                            MAF = c("MAF", "CEUAF", "FREQ1", "EAF", "FREQ1.HAPMAP",
                                    "FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU",
                                    "EFFECT_ALLELE_FREQ", "FREQ.A1", "A1FREQ", "ALLELEFREQ"),
                            Z = c("Z", "ZSCORE", "Z-SCORE", "ZSTATISTIC", "ZSTAT",
                                  "Z-STATISTIC"), SE = c("STDERR", "SE", "STDERRLOGOR",
                                                         "SEBETA", "STANDARDERROR"), DIRECTION = c("DIRECTION",
                                                                                                   "DIREC", "DIRE", "SIGN"))
  full_names <- list(P = "P-value", A1 = "effect allele",
                     A2 = "other allele", effect = "beta or effect", SNP = "rs-id",
                     SE = "standard error", DIRECTION = "direction")
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
      assign(j, utilfuncs[[j]], envir = environment())
    }
  }
  if (all(c("ALT", "REF") %in% hold_names)) {
    .LOG(paste0("Found REF and ALT columns in the summary statistic file ",
                filename, ". Please note that REF will be interpreted as A1 (effect allele) and ALT as A2 (other allele)"),
         print = TRUE, file = log.file)
  }
  if (N_provided) {
    interpreted_names[["N"]] <- NULL
  }
  else {
    if ("NEFF" %in% hold_names | "N_EFF" %in% hold_names |
        "N_EFFECTIVE" %in% hold_names | "SUMNEFF" %in% hold_names) {
      .LOG("Found an NEFF column for sample size. \n\nPlease note that this is likely effective sample size and should only be used for liability h^2 conversion for binary traits and that it should reflect the sum of effective sample sizes across cohorts.\n\nBe aware that some NEFF columns reflect half of the effective sample size; the function will automatically double the column names if recognized [check above in .log file to determine if this is the case].\nIf the Neff value is halved in the summary stats, but not recognized by the munge function, this should be manually doubled prior to running munge.",
           file = log.file)
    }
  }
  for (col in names(interpreted_names)) {
    if (col %in% names(userprovided)) {
      .LOG("Interpreting the ", userprovided[[col]], " column as the ",
           col, " column, as requested", file = log.file)
      hold_names[hold_names == toupper(userprovided[[col]])] <- col
    }
    else if (col %in% hold_names) {
      .LOG("Interpreting the ", col, " column as the ",
           col, " column.", file = log.file)
    }
    else if (any(interpreted_names[[col]] %in% hold_names)) {
      .LOG("Interpreting the ", hold_names[hold_names %in%
                                             interpreted_names[[col]]], " column as the ",
           col, " column.", file = log.file)
      hold_names[hold_names %in% interpreted_names[[col]]] <- col
    }
    else if ((col == "effect")) {
      if (any(interpreted_names[["Z"]] %in% hold_names)) {
        if (!warnz) {
          .LOG("Interpreting the ", hold_names[hold_names %in%
                                                 interpreted_names[["Z"]]], " column as the ",
               col, " column.", file = log.file)
          hold_names[hold_names %in% interpreted_names[["Z"]]] <- col
        }
        else {
          .LOG("There appears to be a Z-statistic column in the summary statistic file ",
               filename, ". Please set linprob to TRUE for binary traits or OLS to true for continuous traits in order to back out the betas or if betas are already available remove this column.",
               print = FALSE, file = log.file)
          warning(paste0("There appears to be a Z-statistic column in the summary statistic file ",
                         filename, ". Please set linprob to TRUE for binary traits or OLS to true for continuous traits in order to back out the betas or if betas are already available remove this column."))
        }
      }
    }
    else {
      if (col %in% warn_for_missing) {
        .LOG("Cannot find ", col, " column, try renaming it to ",
             col, " in the summary statistics file for:",
             filename, file = log.file)
      }
      else if (col %in% stop_on_missing) {
        stop(paste0("Cannot find ", col, " column, try renaming it to ",
                    col, " in the summary statistics file for:",
                    filename))
      }
    }
  }
  if (length(checkforsingle) > 0) {
    for (col in checkforsingle) {
      if (sum(hold_names == col) == 0) {
        .LOG("Cannot find ", full_names[[col]], " column, try renaming it ",
             col, " in the summary statistics file for:",
             filename, file = log.file)
        warning(paste0("Cannot find ", full_names[[col]],
                       " column, try renaming it ", col, " in the summary statistics file for:",
                       filename))
      }
      if (sum(hold_names == col) > 1) {
        .LOG("Multiple columns are being interpreted as the ",
             full_names[[col]], " column, try renaming the column you dont want interpreted to ",
             col, "2 in the summary statistics file for:",
             filename, file = log.file)
        warning(paste0("Multiple columns are being interpreted as the ",
                       full_names[[col]], " column, try renaming the column you dont want interpreted to ",
                       col, "2 in the summary statistics file for:",
                       filename))
      }
    }
  }
  return(hold_names)
}


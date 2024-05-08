rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
launch_dir <- args[1]

if (grepl("rds$", args[2])) {
  assemblies <- readRDS(args[2])
}
if (grepl("csv$", args[2])) {
  assemblies <- read.csv(args[2], sep = ",")
}
if (grepl("tsv$", args[2])) {
  assemblies <- read.csv(args[2], sep = "\t")
}

if ("assembly" %in% names(assemblies) == FALSE) {
  stop("Input table must contain a column called 'assembly'.")
}

if ("relative_path" %in% names(assemblies) == FALSE) {
  stop("Input table must contain a column called 'relative_path'.")
}

if (any(grepl("root", assemblies$assembly, ignore.case = TRUE))) {
   stop("Assembly name cannot be 'root'. Specify another name.")
}

if (any(grepl("^Node_", assemblies$assembly, ignore.case = TRUE))) {
  stop("Assembly name cannot start with 'Node_'. Specify another name.")
}

if (any(!is.na(suppressWarnings(as.numeric(assemblies$assembly))))) {
  stop("Assembly name cannot be a number. Specify another name.")
}

if (any(length(unique(assemblies$assembly)) != length(assemblies$assembly))) {
  stop("Assembly names must be unique.")
}

if (any(grepl("\\.fna|fasta$",basename(assemblies$relative_path))) == FALSE) {
  stop("Assembly files must have either '.fna' or '.fasta' extension.")
}

assemblies$absolute_path <- paste0(launch_dir, "/", assemblies$relative_path)

if (any(file.exists(assemblies$absolute_path) == FALSE)) {
  index <- which(!file.exists(assemblies$absolute_path))

  df <- data.frame(
    assembly = assemblies$assembly[index],
    relative_path = assemblies$relative_path[index],
    absolute_path = assemblies$absolute_path[index]
  )

  write.table(
    df, 
    file = "missing_assemblies.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  missing <- paste(df$assembly, collapse = ", ")

  stop("One or more assembly files could not be found: ", missing)
} else {
  file.create("missing_assemblies.tsv")
}

q50 <- quantile(file.size(assemblies$absolute_path), 0.5)

if (any(file.size(assemblies$absolute_path) < 0.1 * q50)) {
  index <- which(file.size(assemblies$absolute_path) < 0.1 * q50)

  smallfiles <- data.frame(
    assembly = assemblies$assembly[index],
    relative_path = assemblies$relative_path[index],
    absolute_path = assemblies$absolute_path[index],
    file_size = file.size(assemblies$absolute_path[index])
  )

  write.table(
    smallfiles,
    file = "smallfiles.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  smallfiles <- paste(assemblies$assembly[index], collapse = ", ")

  stop("One or more assembly files are very small: ", smallfiles)
} else {
  file.create("smallfiles.tsv")
}

if (any(file.size(assemblies$absolute_path) > 10 * q50)) {
  index <- which(file.size(assemblies$absolute_path) > 10 * q50)

  bigfiles <- data.frame(
    assembly = assemblies$assembly[index],
    relative_path = assemblies$relative_path[index],
    absolute_path = assemblies$absolute_path[index],
    file_size = file.size(assemblies$absolute_path[index])
  )

  write.table(
    bigfiles, 
    file = "bigfiles.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  bigfiles <- paste(assemblies$assembly[index], collapse = ", ")
  stop("One or more assembly files are very large: ", bigfiles)
} else {
  file.create("bigfiles.tsv")
}

# read_error <- sapply(assemblies$absolute_path, function(x) {
#   f <- try(seqinr::read.fasta(x, forceDNAtolower = FALSE), silent = TRUE)
#   out <- ifelse(inherits(f, "try-error"), TRUE, FALSE)
#   return(out)
# })

# if (any(read_error == TRUE)) {

#   index <- which(read_error == TRUE)

#   df <- data.frame(
#     assembly = assemblies$assembly[index],
#     relative_path = assemblies$relative_path[index],
#     absolute_path = assemblies$absolute_path[index],
#     file_size = file.size(assemblies$absolute_path[index])
#   )

#   write.table(
#     df, 
#     file = "read_error.tsv",
#     sep = "\t",
#     row.names = FALSE,
#     quote = FALSE
#   )

#   read_error <- paste(read_error, collapse = ", ")

#   stop("One or more assemblies could not be read: ", read_error)
# } else {
#   file.create("read_error.tsv")
# }

assemblies <- assemblies[, c("assembly", "absolute_path")]

write.table(
    assemblies,
    file = "validated_input.csv",
    sep = ",",
    row.names = FALSE,
    quote = FALSE
)

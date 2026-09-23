#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(monocle3)
  library(cicero)
  library(Matrix)
})


# -------------------------------------------------------------------------
# Command-line interface
# -------------------------------------------------------------------------

option_list <- list(
  make_option(
    "--matrix-dir",
    dest = "matrix_dir",
    type = "character",
    help = paste(
      "Directory containing matrix.mtx.gz, barcodes.tsv.gz,",
      "features.tsv.gz, and peaks.bed."
    )
  ),
  make_option(
    "--output-dir",
    dest = "output_dir",
    type = "character",
    default = "results",
    help = "Output directory [default: %default]."
  ),
  make_option(
    "--target-gene",
    dest = "target_gene",
    type = "character",
    help = "Target gene symbol, for example Macrod2."
  ),
  make_option(
    "--chromosome",
    dest = "chromosome",
    type = "character",
    help = "Target chromosome, for example chr2."
  ),
  make_option(
    "--gene-correlation-threshold",
    dest = "gene_correlation_threshold",
    type = "double",
    default = 0.15,
    help = "Minimum gene–peak correlation [default: %default]."
  ),
  make_option(
    "--coaccessibility-threshold",
    dest = "coaccessibility_threshold",
    type = "double",
    default = 0.20,
    help = "Minimum Cicero co-accessibility score [default: %default]."
  ),
  make_option(
    "--sample-num",
    dest = "sample_num",
    type = "integer",
    default = 100,
    help = "Number of Cicero sample windows [default: %default]."
  ),
  make_option(
    "--seed",
    dest = "seed",
    type = "integer",
    default = 2017,
    help = "Random seed [default: %default]."
  ),
  make_option(
    "--genome",
    dest = "genome",
    type = "character",
    help = paste(
      "Genome specification. Use 'mouse.mm9' for Cicero's built-in",
      "mouse genome or provide a chromosome-sizes file."
    )
  ),
  make_option(
    "--gtf",
    dest = "gtf",
    type = "character",
    default = NULL,
    help = paste(
      "Optional local GTF or GTF.GZ file for plotting connections.",
      "No annotation file is downloaded automatically."
    )
  ),
  make_option(
    "--plot-start",
    dest = "plot_start",
    type = "integer",
    default = NA_integer_,
    help = "Optional start coordinate for a connection plot."
  ),
  make_option(
    "--plot-end",
    dest = "plot_end",
    type = "integer",
    default = NA_integer_,
    help = "Optional end coordinate for a connection plot."
  )
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = paste(
    "Calculate gene–peak correlations and Cicero peak",
    "co-accessibility for a target gene and chromosome."
  )
)

opt <- parse_args(parser)


# -------------------------------------------------------------------------
# Validate parameters
# -------------------------------------------------------------------------

required_arguments <- c(
  "matrix_dir",
  "target_gene",
  "chromosome",
  "genome"
)

missing_arguments <- required_arguments[
  vapply(
    required_arguments,
    function(argument) {
      value <- opt[[argument]]
      is.null(value) || !nzchar(value)
    },
    logical(1)
  )
]

if (length(missing_arguments) > 0) {
  print_help(parser)

  stop(
    paste(
      "Missing required argument(s):",
      paste(
        paste0("--", gsub("_", "-", missing_arguments)),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

if (
  opt$gene_correlation_threshold < -1 ||
    opt$gene_correlation_threshold > 1
) {
  stop(
    "--gene-correlation-threshold must be between -1 and 1.",
    call. = FALSE
  )
}

if (
  opt$coaccessibility_threshold < -1 ||
    opt$coaccessibility_threshold > 1
) {
  stop(
    "--coaccessibility-threshold must be between -1 and 1.",
    call. = FALSE
  )
}

if (opt$sample_num < 1) {
  stop(
    "--sample-num must be at least 1.",
    call. = FALSE
  )
}

matrix_dir <- normalizePath(
  opt$matrix_dir,
  mustWork = TRUE
)

# Use this sanitized name in logs and saved run metadata.
matrix_dir_label <- basename(matrix_dir)

output_dir <- opt$output_dir

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------------------
# Input and output paths
# -------------------------------------------------------------------------

matrix_file <- file.path(
  matrix_dir,
  "matrix.mtx.gz"
)

barcode_file <- file.path(
  matrix_dir,
  "barcodes.tsv.gz"
)

feature_file <- file.path(
  matrix_dir,
  "features.tsv.gz"
)

peak_file <- file.path(
  matrix_dir,
  "peaks.bed"
)

required_files <- c(
  matrix_file,
  barcode_file,
  feature_file,
  peak_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Required input files are missing from",
      matrix_dir_label,
      ":",
      paste(basename(missing_files), collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

gene_peak_output <- file.path(
  output_dir,
  "gene_peak_correlations.csv"
)

all_connections_output <- file.path(
  output_dir,
  "peak_coaccessibility_all.csv"
)

filtered_connections_output <- file.path(
  output_dir,
  "peak_coaccessibility_filtered.csv"
)

gene_matrix_output <- file.path(
  output_dir,
  "gene_expression_matrix.mtx"
)

gene_names_output <- file.path(
  output_dir,
  "genes.tsv"
)

barcode_names_output <- file.path(
  output_dir,
  "barcodes.tsv"
)

umap_output <- file.path(
  output_dir,
  "cicero_umap.pdf"
)

connection_plot_output <- file.path(
  output_dir,
  "cicero_connections.pdf"
)

connection_plot_data_output <- file.path(
  output_dir,
  "cicero_connection_plot_data.rds"
)


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

load_genome_coordinates <- function(genome_argument) {
  if (genome_argument == "mouse.mm9") {
    data(
      "mouse.mm9.genome",
      package = "cicero",
      envir = environment()
    )

    genome_coordinates <- get(
      "mouse.mm9.genome",
      envir = environment()
    )
  } else {
    genome_file <- normalizePath(
      genome_argument,
      mustWork = TRUE
    )

    genome_coordinates <- read.table(
      genome_file,
      header = FALSE,
      stringsAsFactors = FALSE
    )

    if (ncol(genome_coordinates) < 2) {
      stop(
        paste(
          "The genome file must contain at least two columns:",
          "chromosome and chromosome length."
        ),
        call. = FALSE
      )
    }

    genome_coordinates <- genome_coordinates[, 1:2]
    names(genome_coordinates) <- c("V1", "V2")
  }

  genome_coordinates$V1 <- as.character(
    genome_coordinates$V1
  )

  genome_coordinates$V2 <- as.numeric(
    genome_coordinates$V2
  )

  if (anyNA(genome_coordinates$V2)) {
    stop(
      "Genome chromosome lengths must be numeric.",
      call. = FALSE
    )
  }

  genome_coordinates
}


sparse_row_correlations <- function(
  feature_matrix,
  response
) {
  response <- as.numeric(response)

  if (ncol(feature_matrix) != length(response)) {
    stop(
      paste(
        "The feature matrix and response vector",
        "contain different numbers of cells."
      ),
      call. = FALSE
    )
  }

  n <- length(response)

  sum_x <- Matrix::rowSums(feature_matrix)
  sum_x_squared <- Matrix::rowSums(
    feature_matrix ^ 2
  )

  sum_y <- sum(response)
  sum_y_squared <- sum(response ^ 2)

  sum_xy <- as.numeric(
    feature_matrix %*% response
  )

  numerator <- n * sum_xy - sum_x * sum_y

  denominator <- sqrt(
    (n * sum_x_squared - sum_x ^ 2) *
      (n * sum_y_squared - sum_y ^ 2)
  )

  correlations <- numerator / denominator
  correlations[denominator == 0] <- NA_real_

  correlations
}


gene_enhancer_corr <- function(
  matrix_data,
  feature_data,
  target_gene,
  chromosome,
  threshold
) {
  if (ncol(feature_data) < 4) {
    stop(
      paste(
        "features.tsv.gz must contain a fourth column",
        "specifying the chromosome for each feature."
      ),
      call. = FALSE
    )
  }

  gene_indices <- which(
    feature_data$V3 == "Gene Expression" &
      feature_data$V2 == target_gene
  )

  if (length(gene_indices) == 0) {
    stop(
      paste(
        "Target gene was not found:",
        target_gene
      ),
      call. = FALSE
    )
  }

  if (length(gene_indices) > 1) {
    stop(
      paste(
        "More than one feature matched target gene",
        shQuote(target_gene),
        "Use a unique feature identifier."
      ),
      call. = FALSE
    )
  }

  peak_indices <- which(
    feature_data$V3 == "Peaks" &
      feature_data$V4 == chromosome
  )

  if (length(peak_indices) == 0) {
    stop(
      paste(
        "No peaks were found on chromosome",
        chromosome
      ),
      call. = FALSE
    )
  }

  gene_expression <- as.numeric(
    matrix_data[gene_indices, ]
  )

  peak_accessibility <- matrix_data[
    peak_indices,
    ,
    drop = FALSE
  ]

  # Binarize peak accessibility, but preserve gene-expression counts.
  peak_accessibility@x[
    peak_accessibility@x > 0
  ] <- 1

  correlations <- sparse_row_correlations(
    feature_matrix = peak_accessibility,
    response = gene_expression
  )

  results <- data.frame(
    peak = rownames(feature_data)[peak_indices],
    coefficient = correlations,
    stringsAsFactors = FALSE
  )

  results <- results[
    complete.cases(results) &
      results$coefficient >= threshold,
    ,
    drop = FALSE
  ]

  results <- results[
    order(
      results$coefficient,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]

  rownames(results) <- NULL

  results
}


# -------------------------------------------------------------------------
# Load 10x multiome data
# -------------------------------------------------------------------------

message(
  "Reading matrix: ",
  file.path(
    matrix_dir_label,
    basename(matrix_file)
  )
)

# Convert the Matrix Market dgTMatrix into a compressed sparse-column matrix.
indata <- as(
  readMM(matrix_file),
  "CsparseMatrix"
)

cellinfo <- read.table(
  barcode_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(cellinfo) < 1) {
  stop(
    "The barcode file does not contain any columns.",
    call. = FALSE
  )
}

cellinfo <- data.frame(
  cells = cellinfo[[1]],
  row.names = cellinfo[[1]],
  stringsAsFactors = FALSE
)

features <- read.table(
  feature_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

if (ncol(features) < 3) {
  stop(
    paste(
      "The feature file must contain at least",
      "feature ID, feature name, and feature type."
    ),
    call. = FALSE
  )
}

if (nrow(indata) != nrow(features)) {
  stop(
    paste(
      "Matrix rows do not match feature rows:",
      nrow(indata),
      "matrix rows versus",
      nrow(features),
      "feature rows."
    ),
    call. = FALSE
  )
}

if (ncol(indata) != nrow(cellinfo)) {
  stop(
    paste(
      "Matrix columns do not match barcode rows:",
      ncol(indata),
      "matrix columns versus",
      nrow(cellinfo),
      "barcodes."
    ),
    call. = FALSE
  )
}

rownames(features) <- features$V1
rownames(indata) <- rownames(features)
colnames(indata) <- rownames(cellinfo)

peak_indices <- which(
  features$V3 == "Peaks"
)

gene_indices <- which(
  features$V3 == "Gene Expression"
)

if (length(peak_indices) == 0) {
  stop(
    "No features with type 'Peaks' were found.",
    call. = FALSE
  )
}

if (length(gene_indices) == 0) {
  stop(
    "No features with type 'Gene Expression' were found.",
    call. = FALSE
  )
}

peak_names <- rownames(features)[peak_indices]
gene_names <- rownames(features)[gene_indices]


# -------------------------------------------------------------------------
# Calculate gene–peak correlations
# -------------------------------------------------------------------------

message(
  "Calculating peak correlations for ",
  opt$target_gene,
  " on ",
  opt$chromosome,
  "."
)

gene_peak_results <- gene_enhancer_corr(
  matrix_data = indata,
  feature_data = features,
  target_gene = opt$target_gene,
  chromosome = opt$chromosome,
  threshold = opt$gene_correlation_threshold
)

write.csv(
  gene_peak_results,
  gene_peak_output,
  row.names = FALSE
)

message(
  "Retained ",
  nrow(gene_peak_results),
  " gene-correlated peaks."
)


# -------------------------------------------------------------------------
# Save the gene-expression component
# -------------------------------------------------------------------------

gene_expression_matrix <- indata[
  gene_indices,
  ,
  drop = FALSE
]

# writeMM() returns NULL; invisible() prevents it from being printed.
invisible(
  writeMM(
    obj = gene_expression_matrix,
    file = gene_matrix_output
  )
)

writeLines(
  gene_names,
  con = gene_names_output
)

writeLines(
  colnames(gene_expression_matrix),
  con = barcode_names_output
)


# -------------------------------------------------------------------------
# Prepare peak metadata
# -------------------------------------------------------------------------

peakinfo <- read.table(
  peak_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(peakinfo) < 3) {
  stop(
    "peaks.bed must contain at least three columns.",
    call. = FALSE
  )
}

peakinfo <- peakinfo[, 1:3]
names(peakinfo) <- c(
  "chr",
  "bp1",
  "bp2"
)

peakinfo$site_name <- paste0(
  peakinfo$chr,
  ":",
  peakinfo$bp1,
  "-",
  peakinfo$bp2
)

rownames(peakinfo) <- peakinfo$site_name

if (nrow(peakinfo) != length(peak_names)) {
  stop(
    paste(
      "The number of rows in peaks.bed does not match",
      "the number of peak features:",
      nrow(peakinfo),
      "versus",
      length(peak_names)
    ),
    call. = FALSE
  )
}

if (all(peak_names %in% rownames(peakinfo))) {
  peakinfo <- peakinfo[
    peak_names,
    ,
    drop = FALSE
  ]
} else {
  warning(
    paste(
      "Peak feature IDs do not exactly match the coordinates",
      "constructed from peaks.bed. Assuming identical row order."
    )
  )

  rownames(peakinfo) <- peak_names
}


# -------------------------------------------------------------------------
# Construct the Monocle 3/Cicero cell-data set
# -------------------------------------------------------------------------

peak_matrix <- indata[
  peak_indices,
  ,
  drop = FALSE
]

# Binarize only peak accessibility.
peak_matrix@x[
  peak_matrix@x > 0
] <- 1

input_cds <- suppressWarnings(
  new_cell_data_set(
    peak_matrix,
    cell_metadata = cellinfo,
    gene_metadata = peakinfo
  )
)

nonzero_peaks <- Matrix::rowSums(
  SingleCellExperiment::counts(input_cds)
) != 0

input_cds <- input_cds[
  nonzero_peaks,
]

set.seed(opt$seed)

input_cds <- monocle3::detect_genes(
  input_cds
)

input_cds <- estimate_size_factors(
  input_cds
)

input_cds <- preprocess_cds(
  input_cds,
  method = "LSI"
)

input_cds <- reduce_dimension(
  input_cds,
  reduction_method = "UMAP",
  preprocess_method = "LSI"
)


# -------------------------------------------------------------------------
# Save a UMAP plot without requiring clustering or trajectory inference
# -------------------------------------------------------------------------

umap_coordinates <- reducedDims(input_cds)$UMAP

pdf(
  umap_output,
  width = 7,
  height = 6
)

plot(
  umap_coordinates[, 1],
  umap_coordinates[, 2],
  pch = 16,
  cex = 0.4,
  col = "#3366AA80",
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  main = "Cicero input cells"
)

invisible(
  dev.off()
)

cicero_cds <- make_cicero_cds(
  input_cds,
  reduced_coordinates = umap_coordinates
)


# -------------------------------------------------------------------------
# Load genome coordinates
# -------------------------------------------------------------------------

genome_coordinates <- load_genome_coordinates(
  opt$genome
)

target_genome <- genome_coordinates[
  genome_coordinates$V1 == opt$chromosome,
  ,
  drop = FALSE
]

if (nrow(target_genome) == 0) {
  stop(
    paste(
      "Chromosome",
      opt$chromosome,
      "was not found in the genome specification."
    ),
    call. = FALSE
  )
}


# -------------------------------------------------------------------------
# Calculate Cicero co-accessibility
# -------------------------------------------------------------------------

message(
  "Running Cicero on ",
  opt$chromosome,
  " using sample_num = ",
  opt$sample_num,
  "."
)

connections <- run_cicero(
  cicero_cds,
  target_genome,
  sample_num = opt$sample_num
)

write.csv(
  connections,
  all_connections_output,
  row.names = FALSE
)

filtered_connections <- connections[
  complete.cases(connections) &
    connections$coaccess >= opt$coaccessibility_threshold,
  ,
  drop = FALSE
]

filtered_connections <- filtered_connections[
  order(
    filtered_connections$coaccess,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

write.csv(
  filtered_connections,
  filtered_connections_output,
  row.names = FALSE
)

message(
  "Retained ",
  nrow(filtered_connections),
  " connections with co-accessibility >= ",
  opt$coaccessibility_threshold,
  "."
)


# -------------------------------------------------------------------------
# Optional connection plot
# -------------------------------------------------------------------------

plot_region_supplied <- (
  !is.na(opt$plot_start) &&
    !is.na(opt$plot_end)
)

if (
  !is.null(opt$gtf) &&
    plot_region_supplied
) {
  if (
    !requireNamespace(
      "rtracklayer",
      quietly = TRUE
    )
  ) {
    stop(
      paste(
        "Package 'rtracklayer' is required when",
        "--gtf is supplied."
      ),
      call. = FALSE
    )
  }

  gtf_file <- normalizePath(
    opt$gtf,
    mustWork = TRUE
  )

  gene_annotation <- rtracklayer::readGFF(
    gtf_file
  )

  gene_annotation$chromosome <- ifelse(
    grepl(
      "^chr",
      gene_annotation$seqid
    ),
    gene_annotation$seqid,
    paste0(
      "chr",
      gene_annotation$seqid
    )
  )

  gene_annotation$gene <- gene_annotation$gene_id
  gene_annotation$transcript <- gene_annotation$transcript_id
  gene_annotation$symbol <- gene_annotation$gene_name

  pdf(
    connection_plot_output,
    width = 12,
    height = 6
  )

  plot_connections(
    connections,
    opt$chromosome,
    opt$plot_start,
    opt$plot_end,
    gene_model = gene_annotation,
    coaccess_cutoff = opt$coaccessibility_threshold,
    connection_width = 0.5,
    peak_color = "purple",
    gene_model_color = "green",
    collapseTranscripts = "longest"
  )

  invisible(
    dev.off()
  )

  connection_plot_data <- plot_connections(
    connections,
    opt$chromosome,
    opt$plot_start,
    opt$plot_end,
    gene_model = gene_annotation,
    coaccess_cutoff = opt$coaccessibility_threshold,
    connection_width = 0.5,
    peak_color = "purple",
    gene_model_color = "green",
    collapseTranscripts = "gene",
    return_as_list = TRUE
  )

  saveRDS(
    connection_plot_data,
    connection_plot_data_output
  )
} else if (
  !is.null(opt$gtf) ||
    plot_region_supplied
) {
  warning(
    paste(
      "Connection plotting requires --gtf,",
      "--plot-start, and --plot-end together.",
      "The plot was skipped."
    )
  )
}


# -------------------------------------------------------------------------
# Save sanitized provenance
# -------------------------------------------------------------------------

run_parameters <- data.frame(
  parameter = c(
    "matrix_dir",
    "target_gene",
    "chromosome",
    "genome",
    "gene_correlation_threshold",
    "coaccessibility_threshold",
    "sample_num",
    "seed"
  ),
  value = c(
    matrix_dir_label,
    opt$target_gene,
    opt$chromosome,
    basename(opt$genome),
    opt$gene_correlation_threshold,
    opt$coaccessibility_threshold,
    opt$sample_num,
    opt$seed
  ),
  stringsAsFactors = FALSE
)

write.table(
  run_parameters,
  file = file.path(
    output_dir,
    "run_parameters.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

message("Pipeline stage completed successfully.")
message(
  "Outputs written to: ",
  basename(output_dir)
)
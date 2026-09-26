#!/usr/bin/env Rscript
# Differential expression (limma or msqrob2) + ROAST experiment validity.

suppressMessages({
  library(limma)
  library(jsonlite)
})

parse_args <- function(argv) {
  parsed <- list()
  i <- 1
  while (i <= length(argv)) {
    arg <- argv[i]
    if (startsWith(arg, "--")) {
      key <- substring(arg, 3)
      if (i < length(argv) && !startsWith(argv[i + 1], "--")) {
        parsed[[key]] <- argv[i + 1]
        i <- i + 2
      } else {
        parsed[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  parsed
}

params <- parse_args(commandArgs(trailingOnly = TRUE))

args <- list(
  abundance_long_file = params$abundance_long_file,
  samples_file = params$samples_file,
  peptide_abundance_file = params$peptide_abundance_file,
  impute_method = ifelse(is.null(params$impute_method), "MinDet", params$impute_method),
  min_completeness = as.numeric(ifelse(is.null(params$min_completeness), 0.3, params$min_completeness)),
  n_rotations = as.integer(ifelse(is.null(params$n_rotations), 1999, params$n_rotations)),
  min_proteins_for_validity = as.integer(ifelse(is.null(params$min_proteins_for_validity), 20, params$min_proteins_for_validity)),
  validity_alpha = as.numeric(ifelse(is.null(params$validity_alpha), 0.05, params$validity_alpha)),
  output_folder = params$output_folder
)

dir.create(args$output_folder, showWarnings = FALSE, recursive = TRUE)

MIN_REPLICATES <- 2

# @step[id=load]: Loading samples and abundance data
samples <- read.delim(args$samples_file, stringsAsFactors = FALSE)
abundance <- read.delim(args$abundance_long_file, stringsAsFactors = FALSE)
abundance <- merge(abundance, samples[, c("sample", "group")], by = "sample")

gene_by_protein <- unique(abundance[, c("protein", "gene")])
gene_lookup <- setNames(gene_by_protein$gene, gene_by_protein$protein)

write_validity <- function(n_tested, p_value) {
  is_valid <- NA
  if (!is.na(n_tested) && n_tested >= args$min_proteins_for_validity && !is.na(p_value)) {
    is_valid <- p_value <= args$validity_alpha
  }
  write(
    toJSON(list(n_tested = n_tested, omnibus_pvalue = p_value, is_valid = is_valid),
           auto_unbox = TRUE, na = "null"),
    file.path(args$output_folder, "experiment_validity.json")
  )
}

run_limma_path <- function() {
  # @step[id=limma_gate,from=check:no]: Gating proteins by per-group replicate count
  ip_samples <- samples$sample[samples$group == "ip"]
  wcl_samples <- samples$sample[samples$group == "wcl"]
  qualifying <- character(0)
  for (protein in unique(abundance$protein)) {
    protein_rows <- abundance[abundance$protein == protein, ]
    ip_count <- sum(protein_rows$sample %in% ip_samples)
    wcl_count <- sum(protein_rows$sample %in% wcl_samples)
    if (ip_count >= MIN_REPLICATES && wcl_count >= MIN_REPLICATES) {
      qualifying <- c(qualifying, protein)
    }
  }

  if (length(qualifying) == 0) {
    write.table(
      data.frame(protein = character(0), gene = character(0), log2_fold_change = numeric(0),
                 p_value = numeric(0), p_adjusted = numeric(0)),
      file.path(args$output_folder, "differential_expression.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    write_validity(0, NA)
    return(invisible(NULL))
  }
  abundance <- abundance[abundance$protein %in% qualifying, ]

  # @step[id=limma,from=limma_gate]: Running limma differential expression
  abundance$log2_value <- log2(abundance$value + 1.0)

  wide <- reshape(abundance[, c("protein", "sample", "log2_value")],
                   idvar = "protein", timevar = "sample", direction = "wide")
  proteins <- wide$protein
  value_cols <- grep("^log2_value\\.", names(wide), value = TRUE)
  expr <- as.matrix(wide[, value_cols, drop = FALSE])
  colnames(expr) <- sub("^log2_value\\.", "", value_cols)
  rownames(expr) <- proteins

  sample_group <- unique(abundance[, c("sample", "group")])
  group <- factor(sample_group$group[match(colnames(expr), sample_group$sample)], levels = c("wcl", "ip"))

  design <- model.matrix(~group)
  fit <- eBayes(lmFit(expr, design))

  result <- topTable(fit, coef = "groupip", number = Inf, sort.by = "none")
  output <- data.frame(
    protein = proteins,
    gene = gene_lookup[proteins],
    log2_fold_change = result$logFC,
    p_value = result$P.Value,
    p_adjusted = result$adj.P.Val
  )
  write.table(output, file.path(args$output_folder, "differential_expression.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # @step[id=limma_roast,from=limma]: Running ROAST experiment validity test
  roast_expr <- expr[complete.cases(expr), , drop = FALSE]
  set.seed(20260101)
  roast_result <- roast(
    roast_expr, index = seq_len(nrow(roast_expr)), design = design,
    contrast = which(colnames(design) == "groupip"), nrot = args$n_rotations
  )
  write_validity(nrow(roast_expr), roast_result$p.value["Mixed", "P.Value"])
}

run_msqrob2_path <- function() {
  suppressMessages({
    library(QFeatures)
    library(msqrob2)
  })

  # @step[id=msqrob2_load,from=check:yes]: Loading peptide-level abundance
  peptides <- read.delim(args$peptide_abundance_file, stringsAsFactors = FALSE)
  peptides <- merge(peptides, samples[, c("sample", "group")], by = "sample")

  ip_samples <- samples$sample[samples$group == "ip"]
  wcl_samples <- samples$sample[samples$group == "wcl"]

  # @step-if[id=gate,from=msqrob2_load]: Gating peptides by per-group completeness
  qualifying <- character(0)
  if (length(ip_samples) >= MIN_REPLICATES && length(wcl_samples) >= MIN_REPLICATES) {
    all_peptides <- unique(peptides$peptide)
    for (pep in all_peptides) {
      pep_rows <- peptides[peptides$peptide == pep, ]
      ip_present <- sum(pep_rows$sample %in% ip_samples & !is.na(pep_rows$value))
      wcl_present <- sum(pep_rows$sample %in% wcl_samples & !is.na(pep_rows$value))
      if (ip_present / length(ip_samples) >= args$min_completeness &&
          wcl_present / length(wcl_samples) >= args$min_completeness) {
        qualifying <- c(qualifying, pep)
      }
    }
  }

  if (length(qualifying) == 0) {
    write.table(
      data.frame(protein = character(0), gene = character(0), log2_fold_change = numeric(0),
                 p_value = numeric(0), p_adjusted = numeric(0)),
      file.path(args$output_folder, "differential_expression.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    write_validity(0, NA)
    return(invisible(NULL))
  }

  long <- peptides[peptides$peptide %in% qualifying, c("peptide", "protein", "sample", "value")]

  # @step[id=msqrob2,from=gate]: Running msqrob2 differential expression
  wide <- reshape(long[, c("peptide", "sample", "value")], idvar = "peptide", timevar = "sample", direction = "wide")
  value_cols <- grep("^value\\.", names(wide), value = TRUE)
  sample_keys <- sub("^value\\.", "", value_cols)
  names(wide)[match(value_cols, names(wide))] <- sample_keys

  protein_by_peptide <- unique(long[, c("peptide", "protein")])
  wide$protein <- protein_by_peptide$protein[match(wide$peptide, protein_by_peptide$peptide)]

  pe <- readQFeatures(assayData = wide, fnames = "peptide", quantCols = sample_keys, name = "peptideRaw")

  sample_names <- colnames(pe[["peptideRaw"]])
  sample_group <- unique(long[, c("sample")])
  colData(pe)$group <- factor(samples$group[match(sample_names, samples$sample)], levels = c("wcl", "ip"))

  pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
  pe <- normalize(pe, i = "peptideLog", name = "peptideNorm", method = "center.median")
  if (args$impute_method == "none") {
    pe <- addAssay(pe, pe[["peptideNorm"]], name = "peptideImputed")
  } else {
    pe <- impute(pe, method = args$impute_method, i = "peptideNorm", name = "peptideImputed")
  }
  pe <- aggregateFeatures(pe, i = "peptideImputed", fcol = "protein", na.rm = TRUE, name = "protein")

  # @step[id=msqrob2_roast,from=msqrob2]: Running ROAST experiment validity test
  protein_expr <- assay(pe[["protein"]])
  protein_expr <- protein_expr[complete.cases(protein_expr), , drop = FALSE]
  roast_group <- colData(pe)$group
  protein_design <- model.matrix(~roast_group)
  set.seed(20260101)
  roast_result <- roast(
    protein_expr, index = seq_len(nrow(protein_expr)), design = protein_design,
    contrast = which(colnames(protein_design) == "roast_groupip"), nrot = args$n_rotations
  )
  write_validity(nrow(protein_expr), roast_result$p.value["Mixed", "P.Value"])

  pe <- msqrob(object = pe, i = "protein", formula = ~group)
  L <- makeContrast("groupip = 0", parameterNames = c("groupip"))
  pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

  result <- rowData(pe[["protein"]])[["groupip"]]
  output <- data.frame(
    protein = rownames(result),
    gene = gene_lookup[rownames(result)],
    log2_fold_change = result$logFC,
    p_value = result$pval,
    p_adjusted = result$adjPval
  )
  write.table(output, file.path(args$output_folder, "differential_expression.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# @step-if[id=check,from=load]: Has peptide-level data?
if (!is.null(args$peptide_abundance_file) && nzchar(args$peptide_abundance_file)) {
  run_msqrob2_path()
} else {
  run_limma_path()
}

# @step[from=limma_roast+msqrob2_roast]: Differential expression complete
cat("Differential expression complete.\n", file = stderr())

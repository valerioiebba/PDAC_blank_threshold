#!/usr/bin/env Rscript

# Microbiome Contamination Threshold Analysis
# Command line version with optparse
# Created by Iebba Valerio, 6th-7th September 2025
# valerio.iebba@gmail.com
# 
# Usage: Rscript blank_threshold.R \
#   --pos DF_sp_prev20_SELECTED.csv \
#   --blank DF_sp_BLANK.csv \
#   --md MD_samples_blanks.csv \
#   --threshold 99 \
#   --output PDAC_intratumoral

# Usage: Rscript blank_threshold.R --pos DF_sp_prev20_SELECTED.csv --blank DF_sp_BLANK.csv --md MD_samples_blanks.csv --threshold 99 --output PDAC_intratumoral

# Load required libraries
required_packages <- c("data.table", "openxlsx", "ggplot2", "optparse", 
                      "dplyr", "tidyr", "cowplot", "effsize", "stringr", 
                      "scales", "RColorBrewer", "gridExtra", "grid", "parallel", "foreach", "doParallel")

# Function to install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    cat(paste("Installing package:", pkg, "\n"))
    install.packages(pkg, dependencies = TRUE, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# Install and load packages
cat("Installing and loading required packages...\n")
sapply(required_packages, install_if_missing)

# Setup parallel processing
n_cores <- max(1, parallel::detectCores() - 1)
cat(sprintf("Setting up parallel processing with %d cores\n", n_cores))
registerDoParallel(cores = n_cores)

# Parse command line arguments
option_list <- list(
  make_option(c("--pos"), type="character", default=NULL, 
              help="Path to positive samples CSV file", metavar="character"),
  make_option(c("--blank"), type="character", default=NULL,
              help="Path to blank samples CSV file", metavar="character"),
  make_option(c("--md"), type="character", default=NULL,
              help="Path to metadata CSV file", metavar="character"),
  make_option(c("--threshold"), type="numeric", default=99,
              help="Percentile for blank threshold (0-99) [default= %default]", metavar="number"),
  make_option(c("--output"), type="character", default="microbiome_analysis",
              help="Output folder name (used for all output files) [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$pos) || is.null(opt$blank) || is.null(opt$md)) {
  print_help(opt_parser)
  stop("Must specify --pos, --blank, and --md arguments.", call.=FALSE)
}

# Check if files exist
if (!file.exists(opt$pos)) stop(paste("Positive samples file not found:", opt$pos))
if (!file.exists(opt$blank)) stop(paste("Blank samples file not found:", opt$blank))
if (!file.exists(opt$md)) stop(paste("Metadata file not found:", opt$md))

# Create output directory
output_dir <- opt$output
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define output file paths
excel_file <- file.path(output_dir, paste0(opt$output, ".xlsx"))
pdf_file <- file.path(output_dir, paste0(opt$output, ".pdf"))
log_file <- file.path(output_dir, paste0(opt$output, ".log"))

# Start logging
log_conn <- file(log_file, open = "w")
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste0("[", timestamp, "] ", msg)
  cat(log_entry, file = log_conn)
  cat("\n", file = log_conn)
  cat(log_entry, "\n")
  flush(log_conn)
}

log_message("=== MICROBIOME CONTAMINATION THRESHOLD ANALYSIS ===")
log_message("Starting analysis...")
log_message(paste("Command line arguments:"))
log_message(paste("  --pos:", opt$pos))
log_message(paste("  --blank:", opt$blank))
log_message(paste("  --md:", opt$md))
log_message(paste("  --threshold:", opt$threshold))
log_message(paste("  --output:", opt$output))
log_message("")
log_message(paste("Output directory:", output_dir))
log_message(paste("Excel output:", excel_file))
log_message(paste("PDF output:", pdf_file))
log_message(paste("Log output:", log_file))
log_message("")

log_message("STATISTICAL METHODS:")
log_message("- Contamination threshold: Specified percentile of blank samples")
log_message("- Statistical test: Mann-Whitney U test (wilcox.test) for pairwise comparisons")
log_message("- Multiple testing correction: Benjamini-Hochberg FDR correction (method='fdr' in p.adjust)")
log_message("- Effect size: Cliff's delta for pairwise comparisons")
log_message("- Data transformation: log1p for visualizations")
log_message("- Minimum sample size: 3 samples per species-cohort group")
log_message("")

# 1. DATA LOADING AND PREPROCESSING
log_message("1. Loading and preprocessing data...")

# Load data using data.table for speed
blank_data <- fread(opt$blank)
positive_data <- fread(opt$pos)
metadata <- fread(opt$md)

log_message(sprintf("Loaded %d blank samples, %d positive samples, %d metadata entries", 
            nrow(blank_data), nrow(positive_data), nrow(metadata)))

# Convert to data.frame for consistency with dplyr
blank_data <- as.data.frame(blank_data)
positive_data <- as.data.frame(positive_data)
metadata <- as.data.frame(metadata)

# Standardize metadata column names
if ("Status_222" %in% colnames(metadata)) {
  metadata <- metadata %>% rename(Cohort = Status_222)
  log_message("Using 'Status_222' column as Cohort")
} else {
  # Use first non-Group column as cohort
  cohort_col <- setdiff(colnames(metadata), "Group")[1]
  if (length(cohort_col) > 0) {
    metadata <- metadata %>% rename(Cohort = all_of(cohort_col))
    log_message(paste("Using", cohort_col, "column as Cohort"))
  } else {
    stop("No suitable cohort column found in metadata")
  }
}

# Add blank classification
metadata <- metadata %>%
  mutate(
    BlankClass = case_when(
      str_detect(tolower(Cohort), "^nfw$") ~ "nfW",
      str_detect(tolower(Cohort), "^eb$") ~ "EB", 
      str_detect(tolower(Cohort), "^ffpeb$") ~ "FFPEb", 
      str_detect(tolower(Cohort), "^pcrmix$") ~ "PCRmix",
      str_detect(tolower(Cohort), "^pcrmix_w_nfw$") ~ "PCRmix_w_nfW",
      str_detect(tolower(Cohort), "^pcr_nan$") ~ "PCR_NaN",
      TRUE ~ NA_character_
    ),
    IsBlank = !is.na(BlankClass)
  )

blank_cohorts <- unique(metadata$Cohort[metadata$IsBlank == TRUE])
positive_cohorts <- unique(metadata$Cohort[metadata$IsBlank == FALSE | is.na(metadata$IsBlank)])
log_message(paste("Blank cohorts identified:", paste(blank_cohorts, collapse = ", ")))
log_message(paste("Positive cohorts identified:", paste(positive_cohorts, collapse = ", ")))

# Filter datasets to only include samples present in metadata
blank_data <- blank_data %>% filter(Group %in% metadata$Group)
positive_data <- positive_data %>% filter(Group %in% metadata$Group)

# Get species columns (all except Group)
species_cols <- setdiff(colnames(positive_data), "Group")

# Ensure both datasets have the same species columns
missing_in_blank <- setdiff(species_cols, colnames(blank_data))
if (length(missing_in_blank) > 0) {
  for (col in missing_in_blank) {
    blank_data[[col]] <- 0
  }
  log_message(paste("Added", length(missing_in_blank), "missing species columns to blank data"))
}

# Reorder columns to match
blank_data <- blank_data %>% select(Group, all_of(species_cols))
positive_data <- positive_data %>% select(Group, all_of(species_cols))

log_message(sprintf("Processing %d species across datasets", length(species_cols)))

# 2. CONVERT TO LONG FORMAT AND MERGE WITH METADATA
log_message("2. Converting to long format and merging with metadata...")

# Convert to long format
blank_long <- blank_data %>%
  pivot_longer(cols = -Group, names_to = "Species", values_to = "RelAbund") %>%
  left_join(metadata, by = "Group") %>%
  filter(!is.na(Cohort), IsBlank == TRUE) %>%
  mutate(RelAbund = as.numeric(RelAbund),
         RelAbund = pmax(RelAbund, 0)) # Ensure non-negative

positive_long <- positive_data %>%
  pivot_longer(cols = -Group, names_to = "Species", values_to = "RelAbund") %>%
  left_join(metadata, by = "Group") %>%
  filter(!is.na(Cohort), IsBlank == FALSE | is.na(IsBlank)) %>%
  mutate(RelAbund = as.numeric(RelAbund),
         RelAbund = pmax(RelAbund, 0)) # Ensure non-negative

# Apply minimum sample size filter (3 samples per species-cohort)
min_sample_size <- 3

# Filter for groups with sufficient sample size
pos_counts <- positive_long %>%
  count(Species, Cohort) %>%
  filter(n >= min_sample_size)

blank_counts <- blank_long %>%
  count(Species, Cohort) %>%
  filter(n >= min_sample_size)

positive_long <- positive_long %>%
  semi_join(pos_counts, by = c("Species", "Cohort"))

blank_long <- blank_long %>%
  semi_join(blank_counts, by = c("Species", "Cohort"))

log_message(sprintf("After filtering: %d unique blank samples, %d unique positive samples",
            length(unique(blank_long$Group)), length(unique(positive_long$Group))))

# 3. CALCULATE CONTAMINATION THRESHOLDS
log_message("3. Calculating contamination thresholds...")

# Calculate threshold percentile for each species using pooled blanks (parallel)
threshold_percentile <- opt$threshold / 100

# Parallel threshold calculation
compute_thresholds_parallel <- function(blank_data, percentile) {
  species_list <- unique(blank_data$Species)
  
  results <- foreach(species = species_list, .combine = rbind, .packages = c("dplyr")) %dopar% {
    species_data <- blank_data %>% filter(Species == species)
    
    data.frame(
      Species = species,
      Blank_threshold = quantile(species_data$RelAbund, percentile, na.rm = TRUE),
      n_blank_samples = nrow(species_data),
      stringsAsFactors = FALSE
    )
  }
  
  return(results %>% filter(n_blank_samples >= 3))
}

thresholds <- compute_thresholds_parallel(blank_long, threshold_percentile)

log_message(sprintf("Calculated %d%% thresholds for %d species", opt$threshold, nrow(thresholds)))

# 4. STATISTICAL ANALYSIS
log_message("4. Performing statistical analysis...")

# Calculate counts above threshold (parallel by species)
compute_counts_above_threshold_parallel <- function(pos_data, threshold_data) {
  # Join threshold data with positive data
  merged_data <- pos_data %>%
    left_join(threshold_data %>% select(Species, Blank_threshold), by = "Species") %>%
    filter(!is.na(Blank_threshold))
  
  species_list <- unique(merged_data$Species)
  
  results <- foreach(species = species_list, .combine = rbind, .packages = c("dplyr")) %dopar% {
    species_data <- merged_data %>% filter(Species == species)
    
    species_data %>%
      group_by(Species, Cohort) %>%
      summarise(
        n_above_threshold = sum(RelAbund > Blank_threshold, na.rm = TRUE),
        cohort_n = n(),
        frac_above_threshold = mean(RelAbund > Blank_threshold, na.rm = TRUE),
        n_above_or_equal = sum(RelAbund >= Blank_threshold, na.rm = TRUE),
        frac_above_or_equal = mean(RelAbund >= Blank_threshold, na.rm = TRUE),
        # Add strict version that includes zeros when threshold is zero
        strict_frac_above_threshold = ifelse(Blank_threshold[1] == 0, 
                                           mean(RelAbund >= Blank_threshold, na.rm = TRUE),
                                           mean(RelAbund > Blank_threshold, na.rm = TRUE)),
        .groups = "drop"
      )
  }
  
  return(results)
}

counts_above_threshold <- compute_counts_above_threshold_parallel(positive_long, thresholds)

# Descriptive statistics using parallel processing by species
log_message("Computing descriptive statistics using parallel processing...")

# Parallel descriptive statistics function
compute_descriptive_stats_parallel <- function(data_long, sample_type) {
  species_list <- unique(data_long$Species)
  
  results <- foreach(species = species_list, .combine = rbind, .packages = c("dplyr")) %dopar% {
    species_data <- data_long %>% filter(Species == species)
    
    species_data %>%
      group_by(Species, Cohort) %>%
      summarise(
        n = n(),
        mean = mean(RelAbund, na.rm = TRUE),
        median = median(RelAbund, na.rm = TRUE),
        sd = sd(RelAbund, na.rm = TRUE),
        min = min(RelAbund, na.rm = TRUE),
        max = max(RelAbund, na.rm = TRUE),
        q25 = quantile(RelAbund, 0.25, na.rm = TRUE),
        q75 = quantile(RelAbund, 0.75, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  return(results)
}

# Compute stats in parallel
descriptive_stats_pos <- compute_descriptive_stats_parallel(positive_long, "positive")
descriptive_stats_blank <- compute_descriptive_stats_parallel(blank_long, "blank")

# Pairwise comparisons with Mann-Whitney U tests
log_message("5. Performing pairwise comparisons with FDR correction...")
log_message("Using Benjamini-Hochberg FDR correction (method='fdr' in R p.adjust function)")
log_message(sprintf("Using parallel processing with %d cores for faster computation", n_cores))

perform_pairwise_tests_parallel <- function(data_long) {
  species_list <- unique(data_long$Species)
  
  # Parallel processing for each species
  results_list <- foreach(species = species_list, .combine = rbind, .packages = c("dplyr", "effsize")) %dopar% {
    species_data <- data_long %>% filter(Species == species)
    cohorts <- unique(species_data$Cohort)
    
    if (length(cohorts) < 2) return(data.frame())
    
    species_results <- data.frame()
    
    # All pairwise combinations for this species
    for (i in 1:(length(cohorts)-1)) {
      for (j in (i+1):length(cohorts)) {
        cohort1 <- cohorts[i]
        cohort2 <- cohorts[j]
        
        data1 <- species_data %>% filter(Cohort == cohort1) %>% pull(RelAbund)
        data2 <- species_data %>% filter(Cohort == cohort2) %>% pull(RelAbund)
        
        if (length(data1) >= 3 && length(data2) >= 3) {
          # Mann-Whitney U test
          test_result <- wilcox.test(data1, data2, exact = FALSE)
          
          # Calculate Cliff's delta (effect size)
          cliffs_delta <- tryCatch({
            effsize::cliff.delta(data1, data2)$estimate
          }, error = function(e) NA)
          
          species_results <- rbind(species_results, data.frame(
            Species = species,
            Cohort1 = cohort1,
            Cohort2 = cohort2,
            n1 = length(data1),
            n2 = length(data2),
            p_value = test_result$p.value,
            statistic = test_result$statistic,
            cliffs_delta = cliffs_delta,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    return(species_results)
  }
  
  # FDR correction using Benjamini-Hochberg method
  if (nrow(results_list) > 0) {
    results_list$p_adj_fdr <- p.adjust(results_list$p_value, method = "fdr")
    results_list$significant_fdr <- results_list$p_adj_fdr < 0.05
  }
  
  return(results_list)
}

pairwise_results <- perform_pairwise_tests_parallel(positive_long)

log_message(sprintf("Performed %d pairwise comparisons, %d significant after FDR correction",
            nrow(pairwise_results), sum(pairwise_results$significant_fdr, na.rm = TRUE)))

# 6. CREATE VISUALIZATIONS FOR PDF - ALL SPECIES
log_message("6. Creating ggplot2 visualizations for multi-page PDF...")
log_message("Including ALL species present in positive samples dataset")

# Create species plot function with improved P-value positioning
create_species_plot_for_pdf <- function(species_name, page_number = 1) {
  # Filter data for this species
  pos_species <- positive_long %>% filter(Species == species_name)
  blank_species <- blank_long %>% filter(Species == species_name)
  
  if (nrow(pos_species) == 0 && nrow(blank_species) == 0) {
    return(NULL)
  }
  
  # Get threshold
  threshold_val <- thresholds %>% 
    filter(Species == species_name) %>% 
    pull(Blank_threshold)
  
  if (length(threshold_val) == 0) threshold_val <- 0
  
  # Transform data (log1p) - ensure no negative values
  pos_species <- pos_species %>% mutate(RelAbund_log1p = log1p(pmax(RelAbund, 0)))
  blank_species <- blank_species %>% mutate(RelAbund_log1p = log1p(pmax(RelAbund, 0)))
  
  # Get pairwise results for this species (all, not just significant)
  species_pairwise <- pairwise_results %>% 
    filter(Species == species_name, p_adj_fdr < 0.1) %>%  # Show if FDR < 0.1
    arrange(p_adj_fdr) %>%
    slice_head(n = 3)  # Top 3 most significant
  
  # Create plots
  plots <- list()
  
  # Positive samples plot (top half)
  if (nrow(pos_species) > 0) {
    p1 <- ggplot(pos_species, aes(x = Cohort, y = RelAbund_log1p, fill = Cohort)) +
      geom_violin(alpha = 0.7, scale = "width", trim = FALSE) +
      geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 0.8) +
      scale_fill_manual(values = c("H_Adjuvant" = "#E31A1C", "M_Autophagy_deficiency " = "#FF7F00",
                                   "H_Chronic_Pancreatitis" = "#1F78B4", "H_Donor" = "#33A02C",
                                   "H_IPMN" = "#FFFF99",
                                   "H_NeoAdjuvant" = "#B15928", "M_Tumor_precursor" = "#A6CEE3",
                                   "M_Wildtype" = "#FB9A99")) +
      labs(
        title = str_replace_all(species_name, "_", " "),
        y = "Relative abundance (log1p)",
        x = ""
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        panel.background = element_rect(fill = "grey95"),
        axis.title.y = element_text(size = 10)
      )
    
    # Add threshold line with descriptive label (always show, even if threshold is 0)
    if (!is.na(threshold_val)) {
      threshold_label <- sprintf("%dth percentile pooled blanks", opt$threshold)
      
      p1 <- p1 + 
        geom_hline(yintercept = log1p(threshold_val), 
                   color = "red", linetype = "dashed", size = 0.8) +
        annotate("text", 
                 x = Inf, y = log1p(threshold_val), 
                 label = threshold_label,
                 hjust = 1.02, vjust = -0.3,
                 color = "red", size = 3, fontface = "italic")
    }
    
    # Calculate proper y-axis limits to include all violin/box plot shapes
    # Get the full range including violin plot extensions
    violin_data <- pos_species$RelAbund_log1p
    y_min <- min(violin_data, na.rm = TRUE)
    y_max <- max(violin_data, na.rm = TRUE)
    y_span <- y_max - y_min
    
    # Add extra padding for violin plot shapes (25% above, limited below to prevent negatives)
    y_plot_min <- max(0, y_min - y_span * 0.1)  # Don't go below 0 for log1p data
    y_plot_max <- y_max + y_span * 0.25
    
    # Add P-value annotations ABOVE the plots (without Cliff's delta)
    if (nrow(species_pairwise) > 0) {
      # Start annotations well above the violin plot area
      y_start <- y_plot_max + y_span * 0.15
      
      for (i in 1:nrow(species_pairwise)) {
        row <- species_pairwise[i, ]
        
        # P-value formatting only (no Cliff's delta)
        p_text <- if (row$p_adj_fdr < 0.001) {
          "FDR<0.001"
        } else if (row$p_adj_fdr < 0.01) {
          sprintf("FDR=%.3f", row$p_adj_fdr)
        } else {
          sprintf("FDR=%.3f", row$p_adj_fdr)
        }
        
        # Find cohort positions on x-axis
        cohorts <- unique(pos_species$Cohort)
        if (row$Cohort1 %in% cohorts && row$Cohort2 %in% cohorts) {
          x1 <- which(cohorts == row$Cohort1)
          x2 <- which(cohorts == row$Cohort2)
          x_center <- (x1 + x2) / 2
          
          # Current y position for this annotation (well above the data)
          y_current <- y_start + (i - 1) * y_span * 0.12
          
          # Add horizontal line connecting the two cohorts
          p1 <- p1 + 
            annotate("segment", x = x1, xend = x2, 
                    y = y_current, yend = y_current,
                    color = "black", size = 0.5) +
            # Add vertical lines at the ends
            annotate("segment", x = x1, xend = x1, 
                    y = y_current - y_span * 0.02, yend = y_current + y_span * 0.02,
                    color = "black", size = 0.5) +
            annotate("segment", x = x2, xend = x2, 
                    y = y_current - y_span * 0.02, yend = y_current + y_span * 0.02,
                    color = "black", size = 0.5) +
            # Add P-value text above the line
            annotate("text", x = x_center, y = y_current + y_span * 0.04, 
                    label = p_text, hjust = 0.5, vjust = 0, 
                    size = 3, fontface = "bold")
        }
      }
      
      # Final y-axis limits including annotations - use scale_y_continuous instead of coord_cartesian
      final_y_max <- y_start + nrow(species_pairwise) * y_span * 0.12 + y_span * 0.08
      p1 <- p1 + scale_y_continuous(limits = c(y_plot_min, final_y_max), expand = c(0, 0))
    } else {
      # No annotations, just ensure violin plots fit properly - use scale_y_continuous
      p1 <- p1 + scale_y_continuous(limits = c(y_plot_min, y_plot_max), expand = c(0, 0))
    }
    
    plots$positive <- p1
  }
  
  # Blank samples plot (bottom half) 
  if (nrow(blank_species) > 0) {
    p2 <- ggplot(blank_species, aes(x = Cohort, y = RelAbund_log1p, fill = Cohort)) +
      geom_violin(alpha = 0.7, scale = "width", trim = FALSE, adjust = 2) +  # Increase smoothing for better violin shapes
      geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA) +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.6, size = 0.8) +
      scale_fill_manual(values = c("EB" = "#E41A1C", "nfW" = "#377EB8", 
                                   "PCRmix" = "#4DAF4A", "PCRmix_w_nfW" = "#984EA3",
                                   "FFPEb" = "#FF7F00")) +
      labs(
        y = "Relative abundance (log1p)",
        x = ""
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        panel.background = element_rect(fill = "grey95"),
        axis.title.y = element_text(size = 10)
      )
    
    # Calculate proper y-axis limits for blank plots - fix spacing issue
    blank_violin_data <- blank_species$RelAbund_log1p
    blank_y_min <- min(blank_violin_data, na.rm = TRUE)
    blank_y_max <- max(blank_violin_data, na.rm = TRUE)
    blank_y_span <- blank_y_max - blank_y_min
    
    # Handle case when all values are the same (zero span)
    if (blank_y_span == 0 || is.na(blank_y_span)) {
      blank_y_span <- max(0.1, abs(blank_y_max) * 0.1)  # Use 10% of the value or 0.1 minimum
    }
    
    # Add reasonable padding for violin plot shapes (limited below to prevent negatives)
    blank_y_plot_min <- max(0, blank_y_min - blank_y_span * 0.1)  # Don't go below 0 for log1p data
    blank_y_plot_max <- blank_y_max + blank_y_span * 0.25
    
    # Use scale_y_continuous with proper limits to prevent over-expansion
    p2 <- p2 + scale_y_continuous(limits = c(blank_y_plot_min, blank_y_plot_max), expand = c(0, 0))
    
    plots$blank <- p2
  }
  
  # Combine plots in vertical layout exactly like PDF
  if (length(plots) == 2) {
    combined_plot <- plot_grid(plots$positive, plots$blank, 
                              ncol = 1, rel_heights = c(2, 1),
                              align = "v", axis = "lr")
  } else if ("positive" %in% names(plots)) {
    combined_plot <- plots$positive
  } else {
    combined_plot <- plots$blank
  }
  
  # Add page number at bottom center (not bold)
  combined_plot <- combined_plot + 
    draw_label(as.character(page_number), x = 0.5, y = 0.02, 
               hjust = 0.5, vjust = 0, size = 14, fontface = "plain")
  
  return(combined_plot)
}

# Get ALL species present in positive samples for PDF
all_species <- unique(positive_long$Species)
log_message(sprintf("Creating PDF with ALL %d species from positive samples", length(all_species)))
log_message("Generating plots in parallel for faster PDF creation...")

# Generate all plots in parallel first
generate_plots_parallel <- function(species_list) {
  results <- foreach(i = seq_along(species_list), .packages = c("ggplot2", "dplyr", "cowplot", "stringr", "scales", "grid")) %dopar% {
    species <- species_list[i]
    tryCatch({
      plot_result <- create_species_plot_for_pdf(species, page_number = i)
      list(species = species, plot = plot_result, page = i, success = TRUE)
    }, error = function(e) {
      list(species = species, plot = NULL, page = i, success = FALSE, error = e$message)
    })
  }
  return(results)
}

# Generate all plots in parallel
all_plots <- generate_plots_parallel(all_species)

# Create PDF by writing pre-generated plots
pdf(pdf_file, width = 8.5, height = 11, paper = "letter")

for (plot_info in all_plots) {
  if (plot_info$success && !is.null(plot_info$plot)) {
    print(plot_info$plot)
    log_message(sprintf("Added plot for species %d/%d: %s", plot_info$page, length(all_species), plot_info$species))
  } else {
    if (plot_info$success) {
      log_message(sprintf("Skipped species %d/%d (no data): %s", plot_info$page, length(all_species), plot_info$species))
    } else {
      log_message(sprintf("Error with species %d/%d: %s - %s", plot_info$page, length(all_species), plot_info$species, plot_info$error))
    }
  }
}

dev.off()
log_message(sprintf("Multi-page PDF created with %d species: %s", length(all_species), pdf_file))

# 7. CREATE COMPREHENSIVE EXCEL OUTPUT
log_message("7. Creating Excel report...")

# Create comprehensive summary statistics matching the sample Excel format
summary_stats_combined <- descriptive_stats_pos %>%
  left_join(
    counts_above_threshold %>% 
      select(Species, Cohort, n_above_threshold, frac_above_threshold, strict_frac_above_threshold),
    by = c("Species", "Cohort")
  ) %>%
  left_join(
    thresholds %>% select(Species, Blank_threshold),
    by = "Species"
  ) %>%
  mutate(
    threshold_percentile = Blank_threshold,
    samples_above_threshold = n_above_threshold,
    fraction_above_threshold = round(frac_above_threshold, 4),
    strict_fraction_above_threshold = round(strict_frac_above_threshold, 4),
    mean_abundance = round(mean, 6),
    median_abundance = round(median, 6),
    std_dev = round(sd, 6),
    min_abundance = round(min, 6),
    max_abundance = round(max, 6),
    q25_abundance = round(q25, 6),
    q75_abundance = round(q75, 6),
    sample_count = n
  ) %>%
  select(Species, Cohort, sample_count, mean_abundance, median_abundance, 
         std_dev, min_abundance, max_abundance, q25_abundance, q75_abundance,
         threshold_percentile, samples_above_threshold, fraction_above_threshold, strict_fraction_above_threshold)

# Pairwise comparison results with proper formatting
pairwise_formatted <- pairwise_results %>%
  mutate(
    comparison = paste(Cohort1, "vs", Cohort2),
    p_value_formatted = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
    p_adj_fdr_formatted = ifelse(p_adj_fdr < 0.001, "<0.001", sprintf("%.3f", p_adj_fdr)),
    cliffs_delta_formatted = ifelse(is.na(cliffs_delta), "NA", sprintf("%.3f", cliffs_delta)),
    mann_whitney_u = as.numeric(statistic),
    significance = ifelse(significant_fdr, "Significant", "Not Significant")
  ) %>%
  select(Species, comparison, Cohort1, Cohort2, n1, n2, mann_whitney_u, 
         p_value_formatted, p_adj_fdr_formatted, cliffs_delta_formatted, significance)

# Create workbook with multiple sheets
wb <- createWorkbook()

# Sheet 1: Analysis Overview
analysis_overview <- data.frame(
  Parameter = c("Analysis Date", "Analysis Type", "Contamination Threshold", 
                "Statistical Test", "Multiple Testing Correction", "FDR Method Details",
                "Minimum Sample Size", "Total Species Analyzed", "Species with Thresholds", 
                "Unique Positive Samples", "Unique Blank Samples",
                "Total Pairwise Comparisons", "Significant Comparisons (FDR < 0.05)",
                "Input Files", "Blank Samples File", "Positive Samples File", "Metadata File"),
  Value = c(as.character(Sys.time()), "Microbiome Contamination Threshold Analysis",
            paste0(opt$threshold, "th percentile of blank samples"), "Mann-Whitney U test", 
            "Benjamini-Hochberg FDR correction", "method='fdr' in R p.adjust function (standard Benjamini-Hochberg procedure)",
            as.character(min_sample_size), as.character(length(species_cols)), as.character(nrow(thresholds)),
            as.character(length(unique(positive_long$Group))), 
            as.character(length(unique(blank_long$Group))),
            as.character(nrow(pairwise_results)),
            as.character(sum(pairwise_results$significant_fdr, na.rm = TRUE)),
            "3 CSV files", basename(opt$blank), basename(opt$pos), basename(opt$md)),
  stringsAsFactors = FALSE
)

addWorksheet(wb, "Analysis_Overview")
writeData(wb, "Analysis_Overview", analysis_overview)

# Sheet 2: Summary Statistics (main results)
addWorksheet(wb, "Summary_Statistics")
writeData(wb, "Summary_Statistics", summary_stats_combined)

# Sheet 3: Contamination Thresholds
addWorksheet(wb, "Contamination_Thresholds") 
writeData(wb, "Contamination_Thresholds", thresholds)

# Sheet 4: Pairwise Comparisons
addWorksheet(wb, "Pairwise_Comparisons")
writeData(wb, "Pairwise_Comparisons", pairwise_formatted)

# Sheet 5: Significant Results Only
significant_only <- pairwise_formatted %>% filter(significance == "Significant")
addWorksheet(wb, "Significant_Results")
writeData(wb, "Significant_Results", significant_only)

# Sheet 6: Descriptive Stats Positive
addWorksheet(wb, "Positive_Sample_Stats")
writeData(wb, "Positive_Sample_Stats", descriptive_stats_pos)

# Sheet 7: Descriptive Stats Blank
addWorksheet(wb, "Blank_Sample_Stats")
writeData(wb, "Blank_Sample_Stats", descriptive_stats_blank)

# Sheet 8: Raw Data for Violin Plots
# Combine positive and blank data for all species in PDF
violin_plot_data <- bind_rows(
  positive_long %>% 
    mutate(Sample_Type = "Positive") %>%
    select(Group, Species, Cohort, RelAbund, Sample_Type),
  blank_long %>% 
    mutate(Sample_Type = "Blank") %>%
    select(Group, Species, Cohort, RelAbund, Sample_Type)
) %>%
  arrange(Species, Sample_Type, Cohort, Group)

addWorksheet(wb, "Violin_Plot_Raw_Data")
writeData(wb, "Violin_Plot_Raw_Data", violin_plot_data)

# Apply professional formatting
headerStyle <- createStyle(fontSize = 12, fontColour = "white", 
                          fgFill = "#4F81BD", textDecoration = "bold",
                          border = "TopBottomLeftRight", borderColour = "white")

# Format all sheets with headers
sheet_names <- names(wb)
for (sheet_name in sheet_names) {
  # Apply header style to first row of each sheet
  addStyle(wb, sheet_name, headerStyle, rows = 1, 
           cols = 1:10, gridExpand = TRUE)  # Apply to first 10 columns
}

# Save Excel file
saveWorkbook(wb, excel_file, overwrite = TRUE)
log_message(sprintf("Excel report saved: %s", excel_file))

# 8. GENERATE EXCEL DOCUMENTATION
log_message("8. Generating Excel file documentation...")

# Create documentation file path
doc_file <- file.path(output_dir, paste0(opt$output, "_Excel_Documentation.txt"))

# Generate comprehensive documentation
doc_content <- c(
  "===============================================================================",
  "EXCEL FILE DOCUMENTATION: MICROBIOME CONTAMINATION THRESHOLD ANALYSIS",
  "===============================================================================",
  "",
  sprintf("Analysis Date: %s", Sys.time()),
  sprintf("Excel File: %s", basename(excel_file)),
  sprintf("Analysis Parameters: %d%% threshold, Mann-Whitney U test with FDR correction", opt$threshold),
  "",
  "===============================================================================",
  "TAB 1: ANALYSIS_OVERVIEW",
  "===============================================================================",
  "PURPOSE: Summary of analysis parameters and key results",
  "",
  "COLUMNS:",
  "- Parameter: Analysis setting or result type",
  "- Value: Corresponding value or description",
  "",
  "KEY INFORMATION:",
  "- Analysis type and statistical methods used",
  "- Input files and contamination threshold settings",
  "- Summary counts of species, samples, and significant results",
  "",
  "===============================================================================",
  "TAB 2: SUMMARY_STATISTICS", 
  "===============================================================================",
  "PURPOSE: Comprehensive statistical summary for positive sample species-cohort combinations",
  "",
  "COLUMNS:",
  "- Species: Bacterial species name",
  "- Cohort: Sample group/category (positive samples only)",
  "- sample_count: Number of samples in this species-cohort combination",
  "- mean_abundance: Mean relative abundance value (6 decimal precision)",
  "- median_abundance: Median relative abundance value (6 decimal precision)",
  "- std_dev: Standard deviation of relative abundance (6 decimal precision)",
  "- min_abundance: Minimum relative abundance value (6 decimal precision)",
  "- max_abundance: Maximum relative abundance value (6 decimal precision)",
  "- q25_abundance: 25th percentile/first quartile (6 decimal precision)",
  "- q75_abundance: 75th percentile/third quartile (6 decimal precision)",
  "- threshold_percentile: Contamination threshold for this species",
  "- samples_above_threshold: Number of samples exceeding contamination threshold",
  "- fraction_above_threshold: Proportion of samples above threshold (4 decimal precision)",
  "- strict_fraction_above_threshold: Adjusted fraction handling zero thresholds (4 decimal precision)",
  "",
  "NOTES:",
  "- Only includes POSITIVE samples (tissue/biological samples)",
  "- Only species-cohort combinations with ≥3 samples included",
  "- Main results table combining descriptive stats with contamination assessment",
  "",
  "===============================================================================",
  "TAB 3: CONTAMINATION_THRESHOLDS",
  "===============================================================================", 
  "PURPOSE: Contamination threshold values calculated from blank samples",
  "",
  "COLUMNS:",
  "- Species: Bacterial species name",
  sprintf("- Blank_threshold: %d%% percentile of pooled blank samples", opt$threshold),
  "- n_blank_samples: Number of blank samples used for threshold calculation",
  sprintf("- Threshold_Percentile: The percentile used (%d%%)", opt$threshold),
  "",
  "NOTES:",
  "- Thresholds calculated using ALL blank samples pooled together",
  "- Critical reference values for contamination assessment",
  "- Used to determine which positive samples exceed background levels",
  "",
  "===============================================================================",
  "TAB 4: PAIRWISE_COMPARISONS",
  "===============================================================================",
  "PURPOSE: Statistical significance testing between all positive sample groups",
  "",
  "COLUMNS:",
  "- Species: Bacterial species name",
  "- comparison: Description of comparison ('Cohort1 vs Cohort2' format)",
  "- Cohort1: First comparison group name",
  "- Cohort2: Second comparison group name",
  "- n1: Sample size in first group",
  "- n2: Sample size in second group", 
  "- mann_whitney_u: Mann-Whitney U test statistic (numeric)",
  "- p_value_formatted: Raw p-value ('<0.001' or 3 decimal places)",
  "- p_adj_fdr_formatted: FDR-adjusted p-value ('<0.001' or 3 decimal places)",
  "- cliffs_delta_formatted: Effect size measure ('NA' or 3 decimal places)",
  "- significance: 'Significant' or 'Not Significant' based on FDR < 0.05",
  "",
  "STATISTICAL METHODS:",
  "- Test: Mann-Whitney U test (non-parametric, wilcox.test in R)",
  "- Multiple testing correction: Benjamini-Hochberg FDR (p.adjust method='fdr')",
  "- Significance threshold: FDR-adjusted p-value < 0.05",
  "- Effect size: Cliff's delta (-1 to +1 scale, measures non-parametric effect)",
  "",
  "===============================================================================",
  "TAB 5: SIGNIFICANT_RESULTS",
  "===============================================================================",
  "PURPOSE: Filtered view of statistically significant comparisons only",
  "",
  "CONTENT: Subset of Pairwise_Comparisons tab showing only significant results",
  "FILTERING: Only comparisons with FDR-adjusted p < 0.05",
  "",
  "UTILITY:",
  "- Quick identification of meaningful differences",
  "- Reduced dataset for focused analysis",
  "- Priority results for biological interpretation",
  "",
  "===============================================================================",
  "TAB 6: POSITIVE_SAMPLE_STATS",
  "===============================================================================",
  "PURPOSE: Detailed descriptive statistics for positive samples only",
  "",
  "COLUMNS:",
  "- Species: Bacterial species name", 
  "- Cohort: Positive sample group (tissue/biological samples)",
  "- n: Number of samples in this species-cohort combination",
  "- mean: Mean relative abundance (raw values, not rounded)",
  "- median: Median relative abundance (raw values, not rounded)",
  "- sd: Standard deviation of relative abundance (raw values, not rounded)",
  "- min: Minimum relative abundance value (raw values, not rounded)",
  "- max: Maximum relative abundance value (raw values, not rounded)",
  "- q25: 25th percentile/first quartile (raw values, not rounded)",
  "- q75: 75th percentile/third quartile (raw values, not rounded)",
  "",
  "SAMPLE TYPES INCLUDED:",
  paste("- Positive cohorts:", paste(positive_cohorts, collapse = ", ")),
  "",
  "NOTES:",
  "- Only species-cohort combinations with ≥3 samples included",
  "- Raw statistical values (higher precision than Summary_Statistics tab)",
  "",
  "===============================================================================",
  "TAB 7: BLANK_SAMPLE_STATS", 
  "===============================================================================",
  "PURPOSE: Detailed descriptive statistics for blank control samples only",
  "",
  "COLUMNS:",
  "- Species: Bacterial species name", 
  "- Cohort: Blank sample group (control samples)",
  "- n: Number of samples in this species-cohort combination",
  "- mean: Mean relative abundance (raw values, not rounded)",
  "- median: Median relative abundance (raw values, not rounded)",
  "- sd: Standard deviation of relative abundance (raw values, not rounded)",
  "- min: Minimum relative abundance value (raw values, not rounded)",
  "- max: Maximum relative abundance value (raw values, not rounded)",
  "- q25: 25th percentile/first quartile (raw values, not rounded)",
  "- q75: 75th percentile/third quartile (raw values, not rounded)",
  "",
  "SAMPLE TYPES INCLUDED:",
  paste("- Blank cohorts:", paste(blank_cohorts, collapse = ", ")),
  "",
  "UTILITY:",
  "- Assessment of background contamination levels in control samples",
  "- Quality control evaluation of blank performance",
  "- Raw data used for contamination threshold calculations",
  "- Only species-cohort combinations with ≥3 samples included",
  "",
  "===============================================================================",
  "TAB 8: VIOLIN_PLOT_RAW_DATA",
  "===============================================================================",
  "PURPOSE: Complete raw data used for PDF violin plot generation",
  "",
  "COLUMNS:",
  "- Group: Individual sample identifier (unique sample ID)",
  "- Species: Bacterial species name",
  "- Cohort: Sample group/category (e.g., H_Donor, nfW, EB, etc.)",
  "- RelAbund: Raw relative abundance value (exact values used in plots)",
  "- Sample_Type: 'Positive' or 'Blank' classification based on source file",
  "",
  "CONTENT:",
  "- ALL individual data points displayed in PDF violin plots",
  "- Combined positive and blank samples in single table",
  "- Sorted by Species, Sample_Type, Cohort, Group for easy navigation",
  "- Exact data used for log₁₀(RelAbund + 1) transformation in visualizations",
  "",
  "UTILITY:",
  "- Complete reproducibility of PDF plots",
  "- Custom visualization and analysis capability", 
  "- Data verification and quality control",
  "- Full traceability of all plotted points",
  "",
  "===============================================================================",
  "DATA PROCESSING AND QUALITY CONTROL",
  "===============================================================================",
  "",
  "SAMPLE CLASSIFICATION:",
  sprintf("- Blank samples: Determined by source file (%s)", basename(opt$blank)),
  sprintf("- Positive samples: Determined by source file (%s)", basename(opt$pos)),
  sprintf("- Metadata: Used only for cohort assignment (%s)", basename(opt$md)),
  "",
  "FILTERING CRITERIA:",
  sprintf("- Minimum sample size: %d samples per species-cohort group", min_sample_size),
  "- Negative abundance values: Set to zero",
  "- Missing species: Added with zero values",
  "",
  "STATISTICAL METHODS:",
  sprintf("- Contamination threshold: %d%% percentile of pooled blank samples", opt$threshold),
  "- Statistical test: Mann-Whitney U test (wilcox.test in R)",
  "- Multiple testing correction: Benjamini-Hochberg FDR (p.adjust method='fdr')",
  "- Effect size: Cliff's delta",
  "- Significance level: FDR-adjusted p < 0.05",
  "",
  "DATA TRANSFORMATIONS:",
  "- Visualization: log₁₀(relative abundance + 1) for PDF plots only",
  "- Analysis: Raw relative abundance values used for statistics",
  "- Thresholds: Applied to raw values, not transformed data",
  "",
  "===============================================================================",
  "FILE INTEGRATION",
  "===============================================================================",
  "",
  sprintf("RELATED OUTPUT FILES:"),
  sprintf("- Multi-page PDF: %s", basename(pdf_file)),
  sprintf("- Analysis log: %s", basename(log_file)),
  sprintf("- R script: %s", basename(opt$pos)), # This will show the script name used
  "",
  "REPRODUCIBILITY:",
  "- All parameters logged and documented",
  "- Complete data lineage maintained",
  "- Statistical methods explicitly defined",
  "- Quality control steps documented",
  "",
  "===============================================================================",
  sprintf("Documentation generated: %s", Sys.time()),
  "==============================================================================="
)

# Write documentation to file
writeLines(doc_content, doc_file)
log_message(sprintf("Excel documentation saved: %s", doc_file))

# Final summary
log_message("")
log_message("=== ANALYSIS COMPLETE ===")
log_message(sprintf("Analysis completed successfully on: %s", Sys.time()))
log_message(sprintf("- Species analyzed: %d", length(species_cols)))
log_message(sprintf("- Species with thresholds: %d", nrow(thresholds)))
log_message(sprintf("- Total pairwise comparisons: %d", nrow(pairwise_results)))
log_message(sprintf("- Significant comparisons (FDR < 0.05): %d", sum(pairwise_results$significant_fdr, na.rm = TRUE)))
log_message(sprintf("- All species included in PDF: %d", length(all_species)))
log_message("")
log_message("OUTPUT FILES:")
log_message(sprintf("- Multi-page PDF: %s", pdf_file))
log_message(sprintf("- Excel report: %s", excel_file))
log_message(sprintf("- Excel documentation: %s", doc_file))
log_message(sprintf("- Log file: %s", log_file))
log_message("")
log_message("PDF MODIFICATIONS:")
log_message("- Removed Cliff's delta from plot annotations")
log_message("- Positioned P-values above violin/box plots to avoid superimposition")
log_message("- Included ALL species from positive samples dataset")
log_message("")
log_message("Analysis completed successfully!")

# Close log file
close(log_conn)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("PDF: %s\n", pdf_file))
cat(sprintf("Excel: %s\n", excel_file))
cat(sprintf("Log: %s\n", log_file))
cat("Analysis completed successfully!\n")

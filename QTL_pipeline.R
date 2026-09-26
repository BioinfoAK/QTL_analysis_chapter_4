# ==============================================================================
# 1. ENVIRONMENT SETUP & DEPENDENCIES
# ==============================================================================
library(CMplot)
library(qtl2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(GenomicRanges)
library(regioneR)
library(regioneReloaded)
library(gridExtra)

# Custom physical chromosome lengths (bp)
custom_lengths <- c(
  "CM001764.2" = 54367034, 
  "CM001765.2" = 52723162, 
  "CM001766.2" = 76723173,
  "CM001767.2" = 65803235,
  "CM001768.2" = 79327772,
  "CM001769.2" = 70817592,
  "CM001770.2" = 62532141,
  "CM001771.2" = 28684778
)

custom_genome <- data.frame(
  chr   = names(custom_lengths),
  start = 0,
  end   = as.numeric(custom_lengths)
)

custom_genome_fixed <- data.frame(
  chr    = custom_genome$chr,
  length = custom_genome$end
)

# Chromosome label mapping for plots
chromosome_map <- c(
  "CM001764.2" = "Chromosome 1",
  "CM001765.2" = "Chromosome 2",
  "CM001766.2" = "Chromosome 3",
  "CM001767.2" = "Chromosome 4",
  "CM001768.2" = "Chromosome 5",
  "CM001769.2" = "Chromosome 6",
  "CM001770.2" = "Chromosome 7",
  "CM001771.2" = "Chromosome 8"
)

# ==============================================================================
# 2. MARKER DENSITY VISUALIZATION (CMPLOT)
# ==============================================================================
library(CMplot)

# 1. Base Chromosome Mapping & Labels
chromosome_map <- c(
  "CM001764.2" = 1,
  "CM001765.2" = 2,
  "CM001766.2" = 3,
  "CM001767.2" = 4,
  "CM001768.2" = 5,
  "CM001769.2" = 6,
  "CM001770.2" = 7,
  "CM001771.2" = 8
)

custom_chr_labels <- paste("Chromosome", 1:8)

# 2. Build lengthchrom mapped directly to numeric chr IDs (1-8)
lengthchrom <- data.frame(
  marker = c("CM001764.2:54367034", "CM001765.2:52723162", "CM001766.2:76723173", "CM001767.2:65803235",
             "CM001768.2:79327772", "CM001769.2:70817592", "CM001770.2:62532141", "CM001771.2:28684778"), 
  chr    = unname(chromosome_map[c("CM001764.2", "CM001765.2", "CM001766.2", "CM001767.2",
                                   "CM001768.2", "CM001769.2", "CM001770.2", "CM001771.2")]),
  pos    = c(54367034, 52723162, 76723173, 65803235, 79327772, 70817592, 62532141, 28684778)
)

# 3. Robust map prep function
prep_map_data <- function(map_path, length_df, map_vec) {
  df <- read.csv(map_path)
  
  if (any(df$chr %in% names(map_vec))) {
    df$chr <- map_vec[as.character(df$chr)]
  } else {
    df$chr <- as.numeric(gsub("[^0-9]", "", as.character(df$chr)))
  }
  
  full_df <- rbind(df, length_df)
  full_df <- full_df[order(full_df$chr, full_df$pos), ]
  return(full_df)
}

# --- SPop Density Map ---
gmap_spop_full <- prep_map_data("/Users/nastya/Downloads/spop_qtl/map.csv", lengthchrom, chromosome_map)

CMplot(
  gmap_spop_full, 
  plot.type   = "d", 
  bin.size    = 1e6, 
  chr.den.col = c("darkgreen", "yellow", "red"),
  chr.labels  = rev(custom_chr_labels), # REVERSED vector to fix CMplot's top-to-bottom axis assignment
  file        = "jpg", 
  file.name   = "spop_density", 
  dpi         = 600,
  main        = "ProgxDom: Number of SNPs within 1Mb window size",
  file.output = TRUE, 
  verbose     = TRUE, 
  width       = 10, 
  height      = 7
)

# --- KPop Density Map ---
gmap_kpop_full <- prep_map_data("/Users/nastya/Downloads/kpop_qtl/map.csv", lengthchrom, chromosome_map)

CMplot(
  gmap_kpop_full, 
  plot.type   = "d", 
  bin.size    = 1e6, 
  chr.den.col = c("darkgreen", "yellow", "red"),
  chr.labels  = rev(custom_chr_labels), # REVERSED vector to fix CMplot's top-to-bottom axis assignment
  file        = "jpg", 
  file.name   = "kpop_density", 
  dpi         = 600,
  main        = "NDWxDom: Number of SNPs within 1Mb window size",
  file.output = TRUE, 
  verbose     = TRUE, 
  width       = 10, 
  height      = 7
)
# ==============================================================================
# PHENOTYPE HISTOGRAMS, BINARY FREQUENCIES & STATISTICAL TESTING (SPOP vs KPOP)
# ==============================================================================
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")

library(cowplot)
library(dplyr)
library(ggplot2)
library(tools)   # For automatic title capitalization
library(stringr) # For str_wrap

# Defined High-Contrast Palette: SPop (Orange), KPop (Blue)
pop_colors <- c(
  "Progenitor" = "#E66101", 
  "NDW"  = "#0571B0"
)

# Load CSV files (check.names = FALSE keeps exact column names with spaces)
pheno_spop_raw <- read.csv("/Users/nastya/Downloads/spop_qtl/pheno.csv", check.names = FALSE)
pheno_kpop_raw <- read.csv("/Users/nastya/Downloads/kpop_qtl/pheno.csv", check.names = FALSE)

# ------------------------------------------------------------------------------
# CATEGORY TRAIT DEFINITIONS
# ------------------------------------------------------------------------------

# 1. Phenological Traits (Timing / Life Stage)
pheno_traits_cont   <- c("Days to germination","Days to flower","days to pod","days to harvest")

# 2. Morphological Traits (Plant Architecture / Vegetative Features)
morph_traits_cont   <- c("total biomass", "stem thickness", "root weight")
morph_traits_binary <- c("Flower colour", "thick side roots", "Big leaves", "Reddening post harvest", "angle of primary branches", "shed leaves")

# 3. Seed / Yield Traits (Reproduction & Yield Components)
yield_traits_cont   <- c("seed mass", "number undeveloped pods", "number seeds", "shell mass")
yield_traits_binary <- c("seed colour beige", "seed colour orange/red", 
                         "seed colour grey", "seed colour brown", "seed colour yellow", 
                         "seed colour green", "seed pattern", "shape", "smoothness")
color_traits_all  <- c("seed colour beige", "seed colour grey", "seed colour brown", "seed colour yellow",
                       "seed colour orange/red", "seed pattern", "shape", "smoothness" )
yield_seed_traits_all <- c("number seeds","seed mass", "number undeveloped pods", "shell mass")
# Combine per category
pheno_traits_all <- c("Days to germination","Days to flower","days to pod","days to harvest")
morph_traits_all <- c("total biomass", "stem thickness", "thick side roots", "root weight", "Flower colour",
                      "Reddening post harvest", "shed leaves", "Big leaves", "angle of primary branches")
yield_traits_all <- c(yield_traits_cont, yield_traits_binary)

# All Continuous and Binary lists for testing logic
all_cont_traits   <- c(pheno_traits_cont, morph_traits_cont, yield_traits_cont)
all_binary_traits <- c(morph_traits_binary, yield_traits_binary)

plot_list <- list()
stats_results <- data.frame(
  Trait = character(),
  Category = character(),
  Type = character(),
  SPop_N = integer(),
  KPop_N = integer(),
  Test_Statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Function to tag trait category
get_trait_category <- function(trait_name) {
  if (trait_name %in% pheno_traits_all) return("Phenological")
  if (trait_name %in% morph_traits_all) return("Morphological")
  if (trait_name %in% yield_traits_all) return("Seed/Yield")
  return("Uncategorized")
}
# Make sure required package for striped patterns is loaded
if (!requireNamespace("ggpattern", quietly = TRUE)) {
  install.packages("ggpattern")
}
library(ggpattern)
library(ggplot2)

# ------------------------------------------------------------------------------
# PART 1: CONTINUOUS TRAIT ANALYSIS (Wilcoxon Rank-Sum Test & Histograms)
# ------------------------------------------------------------------------------
for (trait in all_cont_traits) {
  if (trait %in% colnames(pheno_spop_raw) && trait %in% colnames(pheno_kpop_raw)) {
    vals_spop <- suppressWarnings(as.numeric(pheno_spop_raw[[trait]]))
    vals_kpop <- suppressWarnings(as.numeric(pheno_kpop_raw[[trait]]))
    
    vals_spop <- vals_spop[!is.na(vals_spop)]
    vals_kpop <- vals_kpop[!is.na(vals_kpop)]
    
    trait_display <- stringr::str_wrap(tools::toTitleCase(tolower(trait)), width = 25)
    wt <- wilcox.test(vals_spop, vals_kpop)
    
    stats_results <- rbind(stats_results, data.frame(
      Trait = trait_display,
      Category = get_trait_category(trait),
      Type = "Continuous",
      SPop_N = length(vals_spop),
      KPop_N = length(vals_kpop),
      Test_Statistic = as.numeric(wt$statistic),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    ))
    
    df_plot <- rbind(
      data.frame(Value = vals_spop, Population = "Progenitor"),
      data.frame(Value = vals_kpop, Population = "NDW")
    )
    
    p_str <- ifelse(wt$p.value < 0.001, "p < 0.001", paste0("p = ", round(wt$p.value, 4)))
    
    p <- ggplot(df_plot, aes(x = Value, fill = Population)) +
      geom_histogram(alpha = 0.55, position = "identity", bins = 20, color = "black", linewidth = 0.2) +
      scale_fill_manual(values = pop_colors) +
      theme_bw() +
      labs(
        title = trait_display, 
        subtitle = paste("Wilcoxon:", p_str), 
        x = "Raw Trait Value", 
        y = "Frequency"
      ) +
      theme(
        plot.title = element_text(size = 10, face = "bold", lineheight = 1.1),
        plot.subtitle = element_text(size = 8.5),
        axis.title = element_text(size = 8.5),
        axis.text = element_text(size = 7.5, color = "black"),
        legend.position = "none"
      )
    
    plot_list[[trait_display]] <- p
  } else {
    message(paste("Warning: Continuous trait", trait, "was not found in one or both pheno files."))
  }
}

# ------------------------------------------------------------------------------
# PART 2: BINARY TRAIT ANALYSIS (Fisher's Exact Test & % Stacked Barplots)
# ------------------------------------------------------------------------------
for (trait in all_binary_traits) {
  if (trait %in% colnames(pheno_spop_raw) && trait %in% colnames(pheno_kpop_raw)) {
    vals_spop <- suppressWarnings(as.numeric(pheno_spop_raw[[trait]]))
    vals_kpop <- suppressWarnings(as.numeric(pheno_kpop_raw[[trait]]))
    
    vals_spop <- vals_spop[!is.na(vals_spop)]
    vals_kpop <- vals_kpop[!is.na(vals_kpop)]
    
    trait_display <- stringr::str_wrap(tools::toTitleCase(tolower(trait)), width = 25)
    
    spop_1s <- sum(vals_spop == 1); spop_0s <- sum(vals_spop == 0)
    kpop_1s <- sum(vals_kpop == 1); kpop_0s <- sum(vals_kpop == 0)
    
    spop_tot <- length(vals_spop)
    kpop_tot <- length(vals_kpop)
    
    contingency_tab <- matrix(
      c(spop_1s, spop_0s, kpop_1s, kpop_0s), 
      nrow = 2, byrow = TRUE
    )
    
    ft <- fisher.test(contingency_tab)
    
    stats_results <- rbind(stats_results, data.frame(
      Trait = trait_display,
      Category = get_trait_category(trait),
      Type = "Binary",
      SPop_N = spop_tot,
      KPop_N = kpop_tot,
      Test_Statistic = as.numeric(ft$estimate), # Odds Ratio
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    ))
    
    # Calculate Percentage per Population
    df_bar <- data.frame(
      Population = factor(rep(c("Progenitor", "NDW"), each = 2), levels = c("Progenitor", "NDW")),
      Observation = factor(rep(c("0", "1"), times = 2), levels = c("0", "1")),
      Percentage = c(
        if (spop_tot > 0) (spop_0s / spop_tot) * 100 else 0,
        if (spop_tot > 0) (spop_1s / spop_tot) * 100 else 0,
        if (kpop_tot > 0) (kpop_0s / kpop_tot) * 100 else 0,
        if (kpop_tot > 0) (kpop_1s / kpop_tot) * 100 else 0
      )
    )
    
    p_str <- ifelse(ft$p.value < 0.001, "p < 0.001", paste0("p = ", round(ft$p.value, 3)))
    
    p_bin <- ggplot(df_bar, aes(x = Population, y = Percentage, fill = Population, pattern = Observation)) +
      geom_col_pattern(
        position = "stack",
        color = "black",
        linewidth = 0.2,
        width = 0.5,
        alpha = 0.8,
        pattern_color = "black",
        pattern_fill = "white",
        pattern_density = 0.15,
        pattern_spacing = 0.05
      ) +
      scale_fill_manual(values = pop_colors) +
      scale_pattern_manual(values = c("0" = "stripe", "1" = "none")) +
      scale_y_continuous(limits = c(0, 100.01), expand = c(0, 0)) +
      theme_bw() +
      labs(
        title = trait_display, 
        subtitle = paste("Fisher's exact:", p_str), 
        x = NULL, 
        y = "Percentage (%)"
      ) +
      theme(
        plot.title = element_text(size = 10, face = "bold", lineheight = 1.1),
        plot.subtitle = element_text(size = 8.5),
        axis.title = element_text(size = 8.5),
        axis.text.y = element_text(size = 7.5, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        legend.position = "none"
      )
    
    plot_list[[trait_display]] <- p_bin
  } else {
    message(paste("Warning: Binary trait", trait, "was not found in one or both pheno files."))
  }
}

# ------------------------------------------------------------------------------
# CHI-SQUARE TEST: OVERALL SEED COLOUR DISTRIBUTION BETWEEN POPULATIONS
# ------------------------------------------------------------------------------

# Define the specific seed colour presence/absence trait columns
seed_color_cols <- c(
  "seed colour beige", 
  "seed colour orange/red", 
  "seed colour grey", 
  "seed colour brown", 
  "seed colour yellow", 
  "seed colour green"
)

# Sum total presences (1s) for each colour trait in Progenitor
spop_color_counts <- sapply(seed_color_cols, function(trait) {
  if (trait %in% colnames(pheno_spop_raw)) {
    sum(suppressWarnings(as.numeric(pheno_spop_raw[[trait]])) == 1, na.rm = TRUE)
  } else {
    0
  }
})

# Sum total presences (1s) for each colour trait in NDW
kpop_color_counts <- sapply(seed_color_cols, function(trait) {
  if (trait %in% colnames(pheno_kpop_raw)) {
    sum(suppressWarnings(as.numeric(pheno_kpop_raw[[trait]])) == 1, na.rm = TRUE)
  } else {
    0
  }
})

# Combine counts into a contingency table (Rows = Colours, Columns = Populations)
color_contingency_table <- cbind(
  Progenitor = spop_color_counts, 
  NDW        = kpop_color_counts
)

# Run Chi-Square Test of Independence
chisq_overall_color <- chisq.test(color_contingency_table)

# Display Contingency Table and Results
cat("=== Seed Colour Contingency Table ===\n")
print(color_contingency_table)

cat("\n=== Chi-Square Test Results ===\n")
print(chisq_overall_color)

# Store raw p-values for each color trait
p_values <- numeric(length(seed_color_cols))
names(p_values) <- seed_color_cols

for (trait in seed_color_cols) {
  spop_1s <- sum(suppressWarnings(as.numeric(pheno_spop_raw[[trait]])) == 1, na.rm = TRUE)
  spop_0s <- sum(suppressWarnings(as.numeric(pheno_spop_raw[[trait]])) == 0, na.rm = TRUE)
  
  kpop_1s <- sum(suppressWarnings(as.numeric(pheno_kpop_raw[[trait]])) == 1, na.rm = TRUE)
  kpop_0s <- sum(suppressWarnings(as.numeric(pheno_kpop_raw[[trait]])) == 0, na.rm = TRUE)
  
  tab <- matrix(c(spop_1s, spop_0s, kpop_1s, kpop_0s), nrow = 2)
  
  # Fisher's exact test (or chisq.test) per color
  p_values[trait] <- fisher.test(tab)$p.value
}

# Adjust p-values for multiple comparisons (Benjamini-Hochberg / FDR)
p_adj_fdr <- p.adjust(p_values, method = "BH")
p_adj_bonf <- p.adjust(p_values, method = "bonferroni")

# Combine into a results data frame
posthoc_summary <- data.frame(
  Color = seed_color_cols,
  Raw_p = round(p_values, 4),
  FDR_p = round(p_adj_fdr, 4),
  Bonferroni_p = round(p_adj_bonf, 4),
  Significant_FDR = p_adj_fdr < 0.05
)

print(posthoc_summary)
# Sum directly from the Chi-Square test object
N <- sum(chisq_overall_color$observed)

print(N)
# ------------------------------------------------------------------------------
# PART 3: SAVE & PRINT STATISTICAL RESULTS, BUILD SHARED LEGEND
# ------------------------------------------------------------------------------
write.csv(stats_results, "Phenotype_Distribution_All_Tests.csv", row.names = FALSE)

# Print results directly to console
cat("\n=== STATISTICAL TEST RESULTS ===\n")
print(stats_results)
cat("================================\n\n")

legend_plot <- ggplot(
  data.frame(Value = 1, Population = c("Progenitor", "NDW")), 
  aes(x = Value, fill = Population)
) +
  geom_histogram(bins = 20) +
  scale_fill_manual(values = pop_colors) +
  theme_bw() + 
  theme(
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold")
  )

shared_legend <- cowplot::get_legend(legend_plot)

# ------------------------------------------------------------------------------
# PART 4: SAVE THE CATEGORY A4 PANELS (WITH ORDERING FUNCTIONALITY)
# ------------------------------------------------------------------------------
save_a4_portrait_panel <- function(plots, filename_prefix, n_cols = 2) {
  if (length(plots) == 0) {
    message(paste("No plots found for", filename_prefix))
    return(NULL)
  }
  
  plots_grid <- cowplot::plot_grid(
    plotlist = plots, 
    ncol = n_cols,
    labels = "AUTO",
    label_size = 11,
    label_fontface = "bold"
  )
  
  final_panel <- cowplot::plot_grid(
    plots_grid, 
    shared_legend, 
    ncol = 1, 
    rel_heights = c(1, 0.04)
  )
  
  filename <- paste0(filename_prefix, "_A4_Portrait_Panel.png")
  ggsave(
    filename, 
    plot = final_panel, 
    width = 8,        # A4 Width (inches) 
    height = 10,      # A4 Height (inches) make 6 for pheno 
    dpi = 300
  )
  message(paste("Saved panel to:", filename))
}

# Helper to filter AND strictly order plot_list based on trait list order
get_plots_for_traits <- function(trait_vec) {
  titles <- stringr::str_wrap(tools::toTitleCase(tolower(trait_vec)), width = 25)
  valid_titles <- titles[titles %in% names(plot_list)]
  # Retains the exact order specified in trait_vec
  plot_list[valid_titles]
}

# ------------------------------------------------------------------------------
# EXPORT PANELS
# To reorder plots on any panel, simply change the order of trait names inside 
# the vector (e.g., morph_traits_all) or pass a custom ordered vector directly.
# ------------------------------------------------------------------------------

# Panel 1: Morphological Traits
save_a4_portrait_panel(
  plots = get_plots_for_traits(morph_traits_all), 
  filename_prefix = "Phenotype_Distributions_Morphological_Traits"
)

# Panel 2: Phenological Traits
save_a4_portrait_panel(
  plots = get_plots_for_traits(pheno_traits_all), 
  filename_prefix = "Phenotype_Distributions_Phenological_Traits"
)

# Panel 3: Seed and Yield Traits (Excluding Color)
save_a4_portrait_panel(
  plots = get_plots_for_traits(yield_seed_traits_all), 
  filename_prefix = "Phenotype_Distributions_Seed_Yield_Traits"
)

# Panel 4: Color Traits
save_a4_portrait_panel(
  plots = get_plots_for_traits(color_traits_all), 
  filename_prefix = "Phenotype_Distributions_Color_Traits"
)

 # ==============================================================================
# 4. GENOTYPE BINNING (PER-MARKER AA FREQUENCY AVERAGED IN 20-SNP BINS)
# ==============================================================================
plot_genotype_bins <- function(geno_path, map_path, pop_label, output_file) {
  # 1. Read files preserving raw strings and special characters
  geno_raw <- read.csv(geno_path, check.names = FALSE, stringsAsFactors = FALSE)
  gmap_df  <- read.csv(map_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  gmap_df$marker <- trimws(as.character(gmap_df$marker))
  
  # 2. Reshape matrix so Rows = Individuals, Columns = Markers
  id_col <- grep("^id$", colnames(geno_raw), ignore.case = TRUE, value = TRUE)
  if (length(id_col) > 0) {
    marker_names <- trimws(as.character(geno_raw[[id_col[1]]]))
    geno_mat <- t(geno_raw[, colnames(geno_raw) != id_col[1], drop = FALSE])
    colnames(geno_mat) <- marker_names
    geno <- as.data.frame(geno_mat, stringsAsFactors = FALSE)
  } else {
    geno <- geno_raw
    colnames(geno) <- trimws(colnames(geno))
  }
  
  # Clean missing values ("-", "", "NA", "N/A") to real R NA
  geno[] <- lapply(geno, function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("-", "", "NA", "N/A", "missing")] <- NA
    return(x)
  })
  
  common_snps <- intersect(colnames(geno), gmap_df$marker)
  
  if (length(common_snps) == 0) {
    stop("Error: No matching markers between files.\n",
         "Geno markers sample: ", paste(head(colnames(geno), 3), collapse = ", "), "\n",
         "Map markers sample: ", paste(head(gmap_df$marker, 3), collapse = ", "))
  }
  
  geno <- geno[, common_snps, drop = FALSE]
  gmap_df <- gmap_df %>% filter(marker %in% common_snps) %>% arrange(match(marker, common_snps))
  
  # 3. Calculate % AA per individual marker using non-NA calls
  marker_valid <- colSums(!is.na(geno)) 
  marker_aa    <- colSums(geno == "AA", na.rm = TRUE) 
  
  # Marker frequency (0.0 to 1.0)
  gmap_df$marker_AA_freq <- marker_aa / marker_valid
  
  # 4. Group markers into 20-SNP bins per chromosome and average frequencies
  bin_size <- 20
  gmap_df <- gmap_df %>%
    group_by(chr) %>%
    mutate(bin = ceiling(row_number() / bin_size)) %>%
    ungroup()
  
  bin_summary <- gmap_df %>%
    group_by(chr, bin) %>%
    summarise(
      start_pos    = min(pos),
      end_pos      = max(pos),
      mid_pos      = mean(pos),
      mean_AA_freq = mean(marker_AA_freq, na.rm = TRUE),
      snp_count    = n(),
      .groups      = "drop"
    )
  
  bin_summary$chr <- factor(bin_summary$chr, levels = names(custom_lengths))
  
  chr_limits <- data.frame(
    chr = factor(names(custom_lengths), levels = names(custom_lengths)), 
    pos = as.numeric(custom_lengths)
  )
  
  # 5. Plot genomic bin segment frequencies
  p <- ggplot(bin_summary, aes(x = mid_pos, y = mean_AA_freq)) +
    geom_blank(data = chr_limits, aes(x = pos, y = 0)) +
    geom_segment(aes(x = start_pos, xend = end_pos, y = mean_AA_freq, yend = mean_AA_freq, color = mean_AA_freq), linewidth = 1.2) +
    scale_color_gradient(low = "#E66101", high = "#0571B0", name = "Mean % AA", limits = c(0, 1)) +
    facet_wrap(~chr, scales = "free_x", ncol = 2, labeller = as_labeller(chromosome_map)) +
    theme_bw() +
    labs(
      title = paste0(pop_label, ": Domesticated Allele (AA) Frequency Across Genome"),
      subtitle = "Calculated per marker and averaged across 20-SNP genomic bins",
      x = "Genomic Position (bp)",
      y = "Domesticated Allele Frequency (% AA)"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  ggsave(output_file, plot = p, width = 10, height = 7, dpi = 300)
}
plot_genotype_bins("/Users/nastya/Downloads/spop_qtl/geno_snps_Q20DP6Miss10MAC2.csv", "/Users/nastya/Downloads/spop_qtl/map_snps_Q20DP6Miss10MAC2.csv", "SPop (Prog x Dom)", "SPop_20SNP_AA_Proportion.png")
plot_genotype_bins("/Users/nastya/Downloads/kpop_qtl/geno_snps_Q20DP6Miss10MAC2.csv", "/Users/nastya/Downloads/kpop_qtl/map_snps_Q20DP6Miss10MAC2.csv", "KPop (NDW x Dom)", "KPop_20SNP_AA_Proportion.png")

# ==============================================================================
# 5. CORE QTL PIPELINE FUNCTION
# ==============================================================================
run_qtl_pipeline <- function(cross_obj, pop_name, bin_list, cont_list) {
  message(paste("Starting analysis for:", pop_name))
  
  plot_dir <- paste0(pop_name, "_QTL_plots")
  if(!dir.exists(plot_dir)) dir.create(plot_dir)
  
  # 1. Map Prep
  gmap <- insert_pseudomarkers(cross_obj$gmap, step=1000000)
  pr <- calc_genoprob(cross_obj, gmap, error_prob=0.002)
  pr <- clean_genoprob(pr)
  
  # 2. Run Scans
  out_bin <- scan1(pr, cross_obj$pheno[, bin_list, drop=FALSE], model="binary", maxit=1000, tol=1e-4)
  out_cont <- scan1(pr, cross_obj$pheno[, cont_list, drop=FALSE], model="normal", maxit=1000, tol=1e-4)
  out_all <- cbind(out_bin, out_cont)
  
  # 3. Permutations
  # NOTE FOR PRODUCTION: Switch n_perm from 100L to 1000L after testing is complete.
  message("Running binary permutations...")
  operm_bin <- scan1perm(genoprobs = pr, pheno = cross_obj$pheno[, bin_list, drop=FALSE], 
                         model = "binary", n_perm = 100L, cores = parallel::detectCores())
  
  message("Running continuous permutations...")
  operm_cont <- scan1perm(genoprobs = pr, pheno = cross_obj$pheno[, cont_list, drop=FALSE], 
                          model = "normal", n_perm = 100L, cores = parallel::detectCores())
  
  thr_bin <- summary(operm_bin, alpha=0.05)
  thr_cont <- summary(operm_cont, alpha=0.05)
  
  # 4. Find All Peaks
  p_bin <- find_peaks(out_bin, gmap, threshold=2.5, drop=1, peakdrop = 2.5)
  p_cont <- find_peaks(out_cont, gmap, threshold=2.5, drop=1, peakdrop = 2.5)
  all_peaks <- rbind(p_bin, p_cont)
  
  if(nrow(all_peaks) == 0) {
    message(paste("No peaks found for", pop_name))
    return(NULL)
  }
  
  # 5. Effect Sizes, PVE Loop, and Plotting
  n_samples <- nrow(cross_obj$pheno)
  results_list <- list()
  
  for(i in 1:nrow(all_peaks)) {
    trait  <- all_peaks$lodcolumn[i]
    chr    <- all_peaks$chr[i]
    pos    <- all_peaks$pos[i]
    m_type <- ifelse(trait %in% bin_list, "binary", "normal")
    
    plot_filename <- file.path(plot_dir, paste0(pop_name, "_Chr", chr, "_", trait, "_LOD.png"))
    png(plot_filename, width = 800, height = 500, res = 120)
    plot(out_all, gmap, lodcolumn = trait, chr = chr, 
         main = paste(pop_name, "-", trait, "- Chromosome", chr),
         col = "slateblue", lwd = 2)
    abline(h = all_peaks$lod[i], lty = 3, col = "gray50")
    dev.off()
    
    coefs <- scan1coef(pr[, chr], cross_obj$pheno[, trait, drop=FALSE], model=m_type)
    target_marker <- find_marker(gmap, chr, pos)
    
    if(target_marker %in% rownames(coefs)) {
      est_effects <- coefs[target_marker, , drop=FALSE]
    } else {
      pos_num <- as.numeric(rownames(coefs))
      est_effects <- coefs[which.min(abs(pos_num - pos)), , drop=FALSE]
    }
    
    current_thr <- if(m_type == "binary") thr_bin["0.05", trait] else thr_cont["0.05", trait]
    is_sig <- all_peaks$lod[i] >= current_thr
    pve_val <- (1 - 10^(-2 * all_peaks$lod[i] / n_samples)) * 100
    
    eff_df <- as.data.frame(est_effects)
    if("AA" %in% colnames(eff_df) & "BB" %in% colnames(eff_df) & "intercept" %in% colnames(eff_df)) {
      eff_df$additive_effect  <- (eff_df$BB - eff_df$AA) / 2
      eff_df$dominance_effect <- eff_df$intercept - ((eff_df$AA + eff_df$BB) / 2)
    } else {
      eff_df$additive_effect  <- NA
      eff_df$dominance_effect <- NA
    }
    
    rownames(eff_df) <- NULL
    results_list[[i]] <- cbind(all_peaks[i, ], eff_df, PVE = pve_val, significant = is_sig, population = pop_name)
  }
  
  return(list(
    summary_table = do.call(rbind, results_list),
    scan_effects_probs = pr,
    gmap = gmap,
    scan_outputs = out_all
  ))
}

# ==============================================================================
# RUN PIPELINE & DATA CLEANING
# ==============================================================================
binary_traits <- c("thick side roots", "shed leaves", "Big leaves", "Reddening post harvest", "Flower colour", 
                   "angle of primary branches", "seed colour beige", "seed colour orange red", "seed colour grey",
                   "seed colour brown", "seed colour yellow", "seed pattern") 

binary_traits_kpop <- c("thick side roots", "shed leaves", "Big leaves", "Reddening post harvest", "Flower colour", 
                   "angle of primary branches", "seed colour beige", "seed colour orange red", "seed colour grey",
                   "seed colour brown", "seed colour yellow", "seed colour green", "seed pattern") 

cont_traits <- c("Days to flower", "days to harvest","shell mass", "total biomass", "root weight",
                 "stem thickness", "seed mass", "number seeds", "Days to germination",  "days to pod" )

cross_spop <- read_cross2("/Users/nastya/Downloads/spop_qtl/cross.yaml")
cross_kpop <- read_cross2("/Users/nastya/Downloads/kpop_qtl/cross.yaml")

# Replace dots with spaces in the cross phenotype matrices
colnames(cross_spop$pheno) <- gsub("\\.", " ", colnames(cross_spop$pheno))
colnames(cross_kpop$pheno) <- gsub("\\.", " ", colnames(cross_kpop$pheno))

summary_spop <- run_qtl_pipeline(cross_spop, "SPop", binary_traits, cont_traits)
summary_kpop <- run_qtl_pipeline(cross_kpop, "KPop", binary_traits_kpop, cont_traits)

peaks_spop <- summary_spop$summary_table
peaks_kpop <- summary_kpop$summary_table

write.csv(peaks_spop, "Summary_SPop_Comprehensive.csv", row.names=FALSE)
write.csv(peaks_kpop, "Summary_KPop_Comprehensive.csv", row.names=FALSE)

peaks_spop_clean <- peaks_spop %>%
  mutate(
    chr = as.character(chr),
    ci_lo = as.numeric(gsub("[^0-9.]", "", as.character(ci_lo))),
    ci_hi = as.numeric(gsub("[^0-9.]", "", as.character(ci_hi)))
  ) %>%
  filter(!is.na(chr) & chr != "" & !is.na(ci_lo) & !is.na(ci_hi))

peaks_kpop_clean <- peaks_kpop %>%
  mutate(
    chr = as.character(chr),
    ci_lo = as.numeric(gsub("[^0-9.]", "", as.character(ci_lo))),
    ci_hi = as.numeric(gsub("[^0-9.]", "", as.character(ci_hi)))
  ) %>%
  filter(!is.na(chr) & chr != "" & !is.na(ci_lo) & !is.na(ci_hi))

# Master LOD Profiles
map_df_spop <- data.frame(
  chr = rep(names(summary_spop$gmap), lengths(summary_spop$gmap)),
  pos = unlist(summary_spop$gmap),
  row.names = NULL
)
lod_profiles_spop <- cbind(chr = map_df_spop$chr, pos = map_df_spop$pos, as.data.frame(summary_spop$scan_outputs))

map_df_kpop <- data.frame(
  chr = rep(names(summary_kpop$gmap), lengths(summary_kpop$gmap)),
  pos = unlist(summary_kpop$gmap),
  row.names = NULL
)
lod_profiles_kpop <- cbind(chr = map_df_kpop$chr, pos = map_df_kpop$pos, as.data.frame(summary_kpop$scan_outputs))

# ==============================================================================
# 6. WITHIN-POPULATION REGIONER & REGIONERELOADED TESTS
# ==============================================================================
spop_gr <- makeGRangesFromDataFrame(peaks_spop_clean, keep.extra.columns=TRUE, ignore.strand=TRUE,
                                    seqnames.field="chr", start.field="ci_lo", end.field="ci_hi")

kpop_gr <- makeGRangesFromDataFrame(peaks_kpop_clean, keep.extra.columns=TRUE, ignore.strand=TRUE,
                                    seqnames.field="chr", start.field="ci_lo", end.field="ci_hi")

# ==============================================================================
# WITHIN-POPULATION REGIONER TESTS (WITH PLOTS & SEED)
# ==============================================================================
run_within_pop_regioner <- function(gr_obj, genome_df, pop_label, seed = 42) {
  message(paste("Running within-population regioneR tests for:", pop_label))
  
  set.seed(seed)
  
  # 1. Circular Randomization (Chromosome-constrained)
  pt_dist_circ <- permTest(
    A = gr_obj, B = gr_obj, ntimes = 1000, 
    randomize.function = circularRandomizeRegions,
    evaluate.function = meanDistance, genome = genome_df, verbose = FALSE
  )
  
  pt_ovlp_circ <- permTest(
    A = gr_obj, B = gr_obj, ntimes = 1000, 
    randomize.function = circularRandomizeRegions,
    evaluate.function = numOverlaps, genome = genome_df, count.once = TRUE, verbose = FALSE
  )
  
  # 2. Per-Chromosome Randomization
  pt_dist_rand <- permTest(
    A = gr_obj, B = gr_obj, ntimes = 1000, 
    randomize.function = randomizeRegions,
    ran.par = list(per.chromosome = TRUE), evaluate.function = meanDistance, 
    genome = genome_df, verbose = FALSE
  )
  
  pt_ovlp_rand <- permTest(
    A = gr_obj, B = gr_obj, ntimes = 1000, 
    randomize.function = randomizeRegions,
    ran.par = list(per.chromosome = TRUE), evaluate.function = numOverlaps, 
    genome = genome_df, count.once = TRUE, verbose = FALSE
  )
  
  # Save plots to 2x2 PDF
  pdf_filename <- paste0(pop_label, "_regioneR_Permutation_Plots.pdf")
  pdf(pdf_filename, width = 10, height = 8)
  par(mfrow = c(2, 2))
  
  plot(pt_dist_circ, main = paste(pop_label, "- Circular: Mean Distance"))
  plot(pt_ovlp_circ, main = paste(pop_label, "- Circular: Num Overlaps"))
  plot(pt_dist_rand, main = paste(pop_label, "- Per-Chr: Mean Distance"))
  plot(pt_ovlp_rand, main = paste(pop_label, "- Per-Chr: Num Overlaps"))
  
  par(mfrow = c(1, 1))
  dev.off()
  
  # Render interactively
  par(mfrow = c(2, 2))
  plot(pt_dist_circ, main = paste(pop_label, "- Circular: Mean Distance"))
  plot(pt_ovlp_circ, main = paste(pop_label, "- Circular: Num Overlaps"))
  plot(pt_dist_rand, main = paste(pop_label, "- Per-Chr: Mean Distance"))
  plot(pt_ovlp_rand, main = paste(pop_label, "- Per-Chr: Num Overlaps"))
  par(mfrow = c(1, 1))
  
  res_table <- data.frame(
    Population = pop_label,
    Randomization = c("Circular", "Circular", "Per-Chromosome", "Per-Chromosome"),
    Metric = c("Mean Distance", "Num Overlaps", "Mean Distance", "Num Overlaps"),
    Observed = c(
      pt_dist_circ$meanDistance$observed, 
      pt_ovlp_circ$numOverlaps$observed,
      pt_dist_rand$meanDistance$observed, 
      pt_ovlp_rand$numOverlaps$observed
    ),
    P_Value = c(
      pt_dist_circ$meanDistance$pval, 
      pt_ovlp_circ$numOverlaps$pval,
      pt_dist_rand$meanDistance$pval, 
      pt_ovlp_rand$numOverlaps$pval
    ),
    Z_Score = c(
      pt_dist_circ$meanDistance$zscore,
      pt_ovlp_circ$numOverlaps$zscore,
      pt_dist_rand$meanDistance$zscore,
      pt_ovlp_rand$numOverlaps$zscore
    ),
    stringsAsFactors = FALSE
  )
  
  return(res_table)
}

within_spop_regioner <- run_within_pop_regioner(spop_gr, custom_genome, "SPop")
within_kpop_regioner <- run_within_pop_regioner(kpop_gr, custom_genome, "KPop")
within_regioner_summary <- rbind(within_spop_regioner, within_kpop_regioner)
write.csv(within_regioner_summary, "Within_Population_regioneR_Global.csv", row.names = FALSE)

# ==============================================================================
# PAIRWISE REGIONER PERMUTATION TEST (CROSSWISEPERMTEST)
# ==============================================================================
# Wrapper to convert S4 meanDistance output into a clean double vector
meanDist_numeric <- function(A, B, ...) {
  res <- meanDistance(A, B, ...)
  return(as.numeric(res$meanDistance))
}

run_pairwise_reloaded <- function(gr_list, genome_df, pop_label, min_regions = 2, seed = 42) {
  message(paste("Running pairwise crosswisePermTest for:", pop_label))
  
  set.seed(seed)
  
  # Auto-split single GRanges into GRangesList by trait
  if (inherits(gr_list, "GRanges")) {
    if ("lodcolumn" %in% colnames(mcols(gr_list))) {
      gr_list <- split(gr_list, gr_list$lodcolumn)
    } else if ("trait" %in% colnames(mcols(gr_list))) {
      gr_list <- split(gr_list, gr_list$trait)
    } else {
      stop("Input GRanges object lacks a 'lodcolumn' or 'trait' metadata column.")
    }
  }
  
  # Filter out trait sets with fewer than `min_regions` to avoid zero-variance warnings
  gr_list <- gr_list[elementNROWS(gr_list) >= min_regions]
  
  pw_res <- crosswisePermTest(
    A = gr_list, 
    ntimes = 1000, 
    randomize.function = circularRandomizeRegions,
    evaluate.function = meanDist_numeric, # Use numeric wrapper to avoid S4 error
    genome = genome_df, 
    verbose = FALSE
  )
  
  pdf_filename <- paste0(pop_label, "_Pairwise_Crosswise_Plots.pdf")
  pdf(pdf_filename, width = 10, height = 8)
  plot(pw_res)
  dev.off()
  
  plot(pw_res)
  
  return(pw_res)
}

pw_spop_reloaded <- run_pairwise_reloaded(spop_gr, custom_genome, "SPop")
pw_kpop_reloaded <- run_pairwise_reloaded(kpop_gr, custom_genome, "KPop")

# ==============================================================================
# 7. MANUAL OVERLAP & INTER-QTL DISTANCE CALCULATIONS (WITHIN POPULATION)
# ==============================================================================
analyze_qtl_geometry <- function(peaks_clean, lod_prof, pop_label) {
  num_peaks <- nrow(peaks_clean)
  overlap_rows <- list()
  distance_rows <- list()
  
  if (num_peaks < 2) return(list(overlaps = data.frame(), distances = data.frame()))
  
  for (i in 1:(num_peaks - 1)) {
    for (j in (i + 1):num_peaks) {
      p1 <- peaks_clean[i, ]
      p2 <- peaks_clean[j, ]
      
      if (p1$chr == p2$chr) {
        dist_bp <- abs(p1$pos - p2$pos)
        distance_rows[[length(distance_rows) + 1]] <- data.frame(
          Population = pop_label,
          Chromosome = p1$chr,
          Trait_1 = p1$lodcolumn,
          Trait_1_Pos = p1$pos,
          Trait_2 = p2$lodcolumn,
          Trait_2_Pos = p2$pos,
          Peak_Distance_bp = dist_bp
        )
        
        start_ovlp <- max(p1$ci_lo, p2$ci_lo)
        end_ovlp   <- min(p1$ci_hi, p2$ci_hi)
        
        if (start_ovlp < end_ovlp) {
          ovlp_len_bp <- end_ovlp - start_ovlp
          
          window_lod <- lod_prof %>%
            filter(chr == p1$chr, pos >= start_ovlp, pos <= end_ovlp)
          
          if (nrow(window_lod) > 0) {
            lod_1 <- window_lod[[p1$lodcolumn]]
            lod_2 <- window_lod[[p2$lodcolumn]]
            
            valid_pts <- which(lod_1 >= 2.5 & lod_2 >= 2.5)
            if (length(valid_pts) > 0) {
              overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
                Population = pop_label,
                Chromosome = p1$chr,
                Trait_1 = p1$lodcolumn,
                Trait_2 = p2$lodcolumn,
                Overlap_Start_bp = start_ovlp,
                Overlap_End_bp = end_ovlp,
                Overlap_Size_bp = ovlp_len_bp,
                Min_LOD_Trait_1 = min(lod_1[valid_pts]),
                Min_LOD_Trait_2 = min(lod_2[valid_pts])
              )
            }
          }
        }
      }
    }
  }
  
  return(list(
    overlaps = do.call(rbind, overlap_rows),
    distances = do.call(rbind, distance_rows)
  ))
}

spop_geom <- analyze_qtl_geometry(peaks_spop_clean, lod_profiles_spop, "SPop")
kpop_geom <- analyze_qtl_geometry(peaks_kpop_clean, lod_profiles_kpop, "KPop")

write.csv(rbind(spop_geom$overlaps, kpop_geom$overlaps), "Within_Population_Validated_Overlaps.csv", row.names = FALSE)
write.csv(rbind(spop_geom$distances, kpop_geom$distances), "Within_Population_InterQTL_Distances.csv", row.names = FALSE)

# ==============================================================================
# 8. INTER-POPULATION COMPARISONS (SPOP vs KPOP)
# ==============================================================================
message("Running inter-population regioneR tests (SPop vs KPop)...")

pt_inter_dist_circ <- permTest(A = spop_gr, B = kpop_gr, ntimes = 1000, 
                               randomize.function = circularRandomizeRegions,
                               evaluate.function = meanDistance, genome = custom_genome, verbose = FALSE)

pt_inter_ovlp_circ <- permTest(A = spop_gr, B = kpop_gr, ntimes = 1000, 
                               randomize.function = circularRandomizeRegions,
                               evaluate.function = numOverlaps, genome = custom_genome, verbose = FALSE)

pt_inter_dist_rand <- permTest(A = spop_gr, B = kpop_gr, ntimes = 1000, 
                               randomize.function = randomizeRegions,
                               ran.par = list(per.chromosome = TRUE),
                               evaluate.function = meanDistance, genome = custom_genome, verbose = FALSE)

pt_inter_ovlp_rand <- permTest(A = spop_gr, B = kpop_gr, ntimes = 1000, 
                               randomize.function = randomizeRegions,
                               ran.par = list(per.chromosome = TRUE),
                               evaluate.function = numOverlaps, genome = custom_genome, verbose = FALSE)

inter_regioner_summary <- data.frame(
  Comparison = "SPop vs KPop",
  Randomization = c("Circular", "Circular", "Per-Chromosome", "Per-Chromosome"),
  Metric = c("Mean Distance", "Num Overlaps", "Mean Distance", "Num Overlaps"),
  Observed = c(pt_inter_dist_circ$meanDistance$observed, pt_inter_ovlp_circ$numOverlaps$observed,
               pt_inter_dist_rand$meanDistance$observed, pt_inter_ovlp_rand$numOverlaps$observed),
  P_Value = c(pt_inter_dist_circ$meanDistance$pval, pt_inter_ovlp_circ$numOverlaps$pval,
              pt_inter_dist_rand$meanDistance$pval, pt_inter_ovlp_rand$numOverlaps$pval)
)
write.csv(inter_regioner_summary, "Inter_Population_regioneR_Global.csv", row.names = FALSE)

# Inter-Population Manual Overlaps & Distances
analyze_inter_qtl_geometry <- function(peaks_A, lod_A, peaks_B, lod_B) {
  overlap_rows <- list()
  distance_rows <- list()
  
  for (i in 1:nrow(peaks_A)) {
    for (j in 1:nrow(peaks_B)) {
      p1 <- peaks_A[i, ]
      p2 <- peaks_B[j, ]
      
      if (p1$chr == p2$chr) {
        dist_bp <- abs(p1$pos - p2$pos)
        distance_rows[[length(distance_rows) + 1]] <- data.frame(
          Chromosome = p1$chr,
          SPop_Trait = p1$lodcolumn,
          SPop_Pos = p1$pos,
          KPop_Trait = p2$lodcolumn,
          KPop_Pos = p2$pos,
          Peak_Distance_bp = dist_bp
        )
        
        start_ovlp <- max(p1$ci_lo, p2$ci_lo)
        end_ovlp   <- min(p1$ci_hi, p2$ci_hi)
        
        if (start_ovlp < end_ovlp) {
          ovlp_len_bp <- end_ovlp - start_ovlp
          
          win_lod_A <- lod_A %>% filter(chr == p1$chr, pos >= start_ovlp, pos <= end_ovlp)
          win_lod_B <- lod_B %>% filter(chr == p2$chr, pos >= start_ovlp, pos <= end_ovlp)
          
          common_positions <- intersect(win_lod_A$pos, win_lod_B$pos)
          if (length(common_positions) > 0) {
            sub_A <- win_lod_A %>% filter(pos %in% common_positions)
            sub_B <- win_lod_B %>% filter(pos %in% common_positions)
            
            lA <- sub_A[[p1$lodcolumn]]
            lB <- sub_B[[p2$lodcolumn]]
            
            valid_pts <- which(lA >= 2.5 & lB >= 2.5)
            if (length(valid_pts) > 0) {
              overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
                Chromosome = p1$chr,
                SPop_Trait = p1$lodcolumn,
                KPop_Trait = p2$lodcolumn,
                Overlap_Start_bp = start_ovlp,
                Overlap_End_bp = end_ovlp,
                Overlap_Size_bp = ovlp_len_bp,
                Min_LOD_SPop = min(lA[valid_pts]),
                Min_LOD_KPop = min(lB[valid_pts])
              )
            }
          }
        }
      }
    }
  }
  return(list(
    overlaps = do.call(rbind, overlap_rows),
    distances = do.call(rbind, distance_rows)
  ))
}

inter_geom <- analyze_inter_qtl_geometry(peaks_spop_clean, lod_profiles_spop, peaks_kpop_clean, lod_profiles_kpop)
write.csv(inter_geom$overlaps, "Inter_Population_Validated_Overlaps.csv", row.names = FALSE)
write.csv(inter_geom$distances, "Inter_Population_InterQTL_Distances.csv", row.names = FALSE)

# Chromosome Map Plotting
plot_data <- bind_rows(peaks_kpop_clean, peaks_spop_clean) %>%
  mutate(
    chr = factor(chr, levels = names(custom_lengths)), 
    qtl_track = paste(population, lodcolumn, sep = " | ")
  )

chr_limits <- data.frame(chr = factor(names(custom_lengths), levels = names(custom_lengths)), pos = as.numeric(custom_lengths))
strip_pop_prefix <- function(x) gsub("^(SPop|KPop) \\| ", "", x)

qtl_map <- ggplot(plot_data, aes(x = pos, y = qtl_track, color = population)) +
  geom_blank(data = chr_limits, aes(x = pos, y = NULL, color = NULL)) +
  geom_blank(aes(x = 0)) + 
  geom_segment(
    aes(x = ci_lo, xend = ci_hi, y = qtl_track, yend = qtl_track), 
    linewidth = 1.2, 
    alpha = 0.7
  ) +
  geom_point(aes(size = lod)) +
  facet_wrap(~chr, scales = "free_y", ncol = 2, labeller = as_labeller(chromosome_map)) + 
  scale_color_manual(
    name = "Population",
    labels = c("KPop" = "NDW x Dom", "SPop" = "Prog x Dom"),
    values = c("KPop" = "#4294F7", "SPop" = "#196B24")
  ) +
  scale_y_discrete(labels = strip_pop_prefix) +
  theme_bw(base_size = 11) + # Set global base text size for A4 scaling
  labs(
    title = "Chromosome-Specific QTL Map", 
    x = "Genomic Position (bp)", 
    size = "LOD Score"
  ) +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    title = element_text(size = 13, face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  )
ggsave(
  filename = "QTL_map_Plot_A4.png", 
  plot = qtl_map, 
  width = 11.69, 
  height = 10.0, 
  units = "in", 
  dpi = 300
)
# ==============================================================================
# SPOP vs KPOP: QTL COUNT PER TRAIT STATISTICAL COMPARISON
# ==============================================================================
library(dplyr)
library(ggplot2)

# Normalization helper preserving single spaces between words
clean_trait_names <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00a0", " ", x)        # Remove non-breaking spaces
  x <- gsub("[_.-]", " ", x)         # Replace underscores/dots/hyphens with space
  x <- gsub("[[:space:]]+", " ", x)  # Collapse multiple spaces
  x <- trimws(tolower(x))            # Lowercase and trim edges
  return(x)
}

pop_colors <- c(
  "SPop (Prog x Dom)" = "#E66101", 
  "KPop (NDW x Dom)"  = "#0571B0"
)

# 1. Summarize counts using safe trait names
qtl_counts_spop <- peaks_spop_clean %>% 
  mutate(trait_clean = clean_trait_names(lodcolumn)) %>% 
  group_by(trait_clean) %>% 
  summarise(SPop_QTLs = n(), .groups = "drop")

qtl_counts_kpop <- peaks_kpop_clean %>% 
  mutate(trait_clean = clean_trait_names(lodcolumn)) %>% 
  group_by(trait_clean) %>% 
  summarise(KPop_QTLs = n(), .groups = "drop")

# 2. Build master traits from UNION to ensure ZERO missing phenotypes
all_unique_traits <- unique(c(
  clean_trait_names(cont_traits),
  qtl_counts_spop$trait_clean,
  qtl_counts_kpop$trait_clean
))

all_traits <- data.frame(
  trait = all_unique_traits,
  stringsAsFactors = FALSE
)

# 3. Complete Join
qtl_comp <- all_traits %>%
  left_join(qtl_counts_spop, by = c("trait" = "trait_clean")) %>%
  left_join(qtl_counts_kpop, by = c("trait" = "trait_clean")) %>%
  mutate(
    SPop_QTLs = ifelse(is.na(SPop_QTLs), 0, SPop_QTLs),
    KPop_QTLs = ifelse(is.na(KPop_QTLs), 0, KPop_QTLs),
    Diff = SPop_QTLs - KPop_QTLs
  )

# ------------------------------------------------------------------------------
# STATISTICAL TESTING & VISUALIZATION
# ------------------------------------------------------------------------------
wilcox_res <- wilcox.test(qtl_comp$SPop_QTLs, qtl_comp$KPop_QTLs, paired = TRUE)
ttest_res  <- t.test(qtl_comp$SPop_QTLs, qtl_comp$KPop_QTLs, paired = TRUE)

cat("\n=======================================================\n")
cat("          QTL COUNT COMPARISON (SPop vs KPop)          \n")
cat("=======================================================\n")
cat(sprintf("Mean QTLs per trait  -> SPop: %.2f | KPop: %.2f\n", mean(qtl_comp$SPop_QTLs), mean(qtl_comp$KPop_QTLs)))
cat(sprintf("Median QTLs per trait -> SPop: %.1f | KPop: %.1f\n", median(qtl_comp$SPop_QTLs), median(qtl_comp$KPop_QTLs)))
cat("-------------------------------------------------------\n")
cat(sprintf("Paired Wilcoxon Signed-Rank Test: V = %g, p-value = %.4f\n", wilcox_res$statistic, wilcox_res$p.value))
cat(sprintf("Paired t-test: t = %.3f, df = %d, p-value = %.4f\n", ttest_res$statistic, ttest_res$parameter, ttest_res$p.value))
cat("=======================================================\n\n")

write.csv(qtl_comp, "QTL_Counts_Per_Trait_Comparison.csv", row.names = FALSE)

qtl_long <- qtl_comp %>%
  tidyr::pivot_longer(
    cols = c("SPop_QTLs", "KPop_QTLs"),
    names_to = "Population",
    values_to = "QTL_Count"
  ) %>%
  mutate(
    Population = ifelse(Population == "SPop_QTLs", "SPop (Prog x Dom)", "KPop (NDW x Dom)"),
    trait_display = tools::toTitleCase(trait)
  )

p_val_text <- ifelse(wilcox_res$p.value < 0.001, "p < 0.001", sprintf("p = %.4f", wilcox_res$p.value))

p_qtl_count <- ggplot(qtl_long, aes(x = Population, y = QTL_Count, fill = Population)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.4) +
  geom_line(aes(group = trait_display), color = "grey60", linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(color = Population), size = 3, alpha = 0.9) +
  scale_fill_manual(values = pop_colors) +
  scale_color_manual(values = pop_colors) +
  scale_y_continuous(breaks = function(x) seq(0, max(x) + 1, by = 1)) +
  theme_bw() +
  labs(
    title = "Comparison of Detected QTL Counts Per Trait",
    subtitle = paste("Paired Wilcoxon Signed-Rank Test:", p_val_text),
    x = NULL,
    y = "Number of Detected QTLs"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black")
  )

plot(p_qtl_count)
ggsave("QTL_Counts_Comparison_Plot.png", plot = p_qtl_count, width = 6, height = 5, dpi = 300)

#MANUAL CHECK FOR DIFFERENCE IN QTL COUNT PER CHROMOSOME

# 1. Manually enter your QTL counts per chromosome (Chr 1 to Chr 8)
pop1_counts <- c(3, 1, 3, 8, 13, 6, 3, 0) # Replace with Pop 1 counts
pop2_counts <- c(2, 5, 8, 11, 7, 9, 1, 3) # Replace with Pop 2 counts

# Define chromosome names
chr_names <- paste0("Chr", 1:8)

# 2. Build the contingency table
qtl_table <- matrix(
  c(pop1_counts, pop2_counts),
  ncol = 2,
  dimnames = list(Chromosome = chr_names, Population = c("Prog x Dom", "NDW x Dom"))
)

# View your matrix to double-check numbers
print(qtl_table)

# 3. Run Fisher's Exact Test
# simulate.p.value = TRUE uses Monte Carlo simulation, which is standard practice
# for tables larger than 2x2 to ensure accuracy and prevent computational memory errors.
fisher_result <- fisher.test(qtl_table, simulate.p.value = TRUE, B = 100000)

# Print test results
print(fisher_result)

# 1. Recreate the dataset directly from your image
qtl_data <- data.frame(
  trait = c(
    "days to harvest", "days to pod", "seed mass", "days to flower", 
    "days to germination", "number undeveloped pods", "number seeds", 
    "shell mass", "stem thickness", "root weight", "total biomass", 
    "angle of primary branches", "big leaves", "flower colour", 
    "reddening post harvest", "seed colour beige", "seed colour grey", 
    "seed colour orange red", "seed pattern", "shed leaves", 
    "thick side roots", "seed colour brown", "seed colour green", 
    "seed colour yellow"
  ),
  SPop_QTLs = c(2, 4, 0, 2, 0, 2, 2, 0, 1, 2, 5, 3, 1, 2, 2, 2, 1, 3, 1, 1, 1, 0, 0, 0),
  KPop_QTLs = c(2, 0, 0, 3, 2, 2, 0, 0, 2, 0, 2, 12, 0, 2, 2, 1, 3, 0, 0, 2, 5, 2, 2, 2)
)

# 2. Extract matrix for testing (excluding rows where both populations have 0 QTLs)
qtl_matrix <- as.matrix(qtl_data[, c("SPop_QTLs", "KPop_QTLs")])
rownames(qtl_matrix) <- qtl_data$trait

# Filter out empty rows (like 'seed mass' and 'shell mass')
qtl_matrix_clean <- qtl_matrix[rowSums(qtl_matrix) > 0, ]

# View the matrix
print(qtl_matrix_clean)

# 3. Run Fisher's Exact Test with Monte Carlo simulation
# (Essential here since the matrix is 22 x 2)
fisher_result <- fisher.test(qtl_matrix_clean, simulate.p.value = TRUE, B = 100000)

# Print test results
print(fisher_result)

library(regioneReloaded)
library(GenomicRanges)
library(ggplot2)

# Set seed for reproducible permutation tests
set.seed(123)

# ==============================================================================
# HELPER FUNCTION: Convert crosswise result object to a single combined Table
# ==============================================================================
extract_results_table <- function(cw_obj, test_label = "") {
  # cw_obj stores individual results inside the @multiOverlaps slot
  res_list <- cw_obj@multiOverlaps
  
  df_list <- lapply(names(res_list), function(regionA_name) {
    df <- res_list[[regionA_name]]
    df$RegionA <- regionA_name
    df$Comparison <- paste(df$RegionA, "vs", df$name)
    df$Test_Type <- test_label
    
    # Reorder and select relevant columns
    cols <- c("Test_Type", "Comparison", "RegionA", "name", 
              "n_regionA", "n_regionB", "z_score", "p_value", 
              "n_overlaps", "mean_perm_test", "norm_zscore", "adj.p_value")
    
    # Keep columns that exist
    cols_present <- intersect(cols, colnames(df))
    return(df[, cols_present])
  })
  
  do.call(rbind, df_list)
}


# ==============================================================================
# SECTION 1: Within and Between Population QTL Clustering
# ==============================================================================

# 1. Define list of populations
qtl_list <- list(
  Pop1_QTLs = spop_gr,
  Pop2_QTLs = kpop_gr
)

# 1A. Number of Overlaps + Randomize Regions
cw_ovlp_rand <- crosswisePermTest(
  Alist = qtl_list, Blist = qtl_list,
  ntimes = 1000,
  ranFUN = "randomizeRegions",
  ran.par = list(per.chromosome = TRUE),
  evFUN = "numOverlaps",
  genome = custom_genome,
  mc.cores = 4
)

# 1B. Number of Overlaps + Circular Randomize Regions
cw_ovlp_circ <- crosswisePermTest(
  Alist = qtl_list, Blist = qtl_list,
  ntimes = 1000,
  ranFUN = "circularRandomizeRegions",
  evFUN = "numOverlaps",
  genome = custom_genome,
  mc.cores = 4
)

# 1C. Mean Distance + Randomize Regions
cw_dist_rand <- crosswisePermTest(
  Alist = qtl_list, Blist = qtl_list,
  ntimes = 1000,
  ranFUN = "randomizeRegions",
  ran.par = list(per.chromosome = TRUE),
  evFUN = "meanDistance",
  genome = custom_genome,
  mc.cores = 4
)

# 1D. Mean Distance + Circular Randomize Regions
cw_dist_circ <- crosswisePermTest(
  Alist = qtl_list, Blist = qtl_list,
  ntimes = 1000,
  ranFUN = "circularRandomizeRegions",
  evFUN = "meanDistance",
  genome = custom_genome,
  mc.cores = 4
)

# Combine all 4 tests into one summary table
pop_summary_table <- rbind(
  extract_results_table(cw_ovlp_rand, "Overlaps (Randomize)"),
  extract_results_table(cw_ovlp_circ, "Overlaps (Circular)"),
  extract_results_table(cw_dist_rand, "MeanDistance (Randomize)"),
  extract_results_table(cw_dist_circ, "MeanDistance (Circular)")
)

print(pop_summary_table)
write.csv(pop_summary_table, "Population_QTL_Clustering_Summary.csv", row.names = FALSE)

# ==============================================================================
# SECTION 1.5: Within and Between Population QTL Clustering not per chromosome
# ==============================================================================

# 1. Define list of populations
qtl_list <- list(
  Pop1_QTLs = spop_gr,
  Pop2_QTLs = kpop_gr
)

# 1A. Number of Overlaps + Randomize Regions
cw_ovlp_rand <- crosswisePermTest(
  Alist = qtl_list, Blist = qtl_list,
  ntimes = 1000,
  ranFUN = "randomizeRegions",
  ran.par = list(per.chromosome = FALSE),
  evFUN = "numOverlaps",
  genome = custom_genome,
  mc.cores = 4
)

# 1B. Number of Overlaps + Circular Randomize Regions
cw_ovlp_circ <- crosswisePermTest(
  Alist = qtl_list, Blist = qtl_list,
  ntimes = 1000,
  ranFUN = "circularRandomizeRegions",
  ran.par = list(per.chromosome = FALSE),
  evFUN = "numOverlaps",
  genome = custom_genome,
  mc.cores = 4
)

# 1C. Mean Distance + Randomize Regions
cw_dist_rand <- crosswisePermTest(
  Alist = qtl_list, Blist = qtl_list,
  ntimes = 1000,
  ranFUN = "randomizeRegions",
  ran.par = list(per.chromosome = FALSE),
  evFUN = "meanDistance",
  genome = custom_genome,
  mc.cores = 4
)

# 1D. Mean Distance + Circular Randomize Regions
cw_dist_circ <- crosswisePermTest(
  Alist = qtl_list, Blist = qtl_list,
  ntimes = 1000,
  ranFUN = "circularRandomizeRegions",
  ran.par = list(per.chromosome = FALSE),
  evFUN = "meanDistance",
  genome = custom_genome,
  mc.cores = 4
)

# Combine all 4 tests into one summary table
pop_summary_table <- rbind(
  extract_results_table(cw_ovlp_rand, "Overlaps (Randomize)"),
  extract_results_table(cw_ovlp_circ, "Overlaps (Circular)"),
  extract_results_table(cw_dist_rand, "MeanDistance (Randomize)"),
  extract_results_table(cw_dist_circ, "MeanDistance (Circular)")
)

print(pop_summary_table)
write.csv(pop_summary_table, "Population_QTL_Clustering_Summary_noperchrom.csv", row.names = FALSE)





# ==============================================================================
# SECTION 2: Trait vs Trait Comparisons Within Pop1 and Within Pop2
# ==============================================================================

# Split GRanges by trait (e.g., 'lodcolumn' column)
pop1_by_trait <- split(spop_gr, spop_gr$lodcolumn)
pop2_by_trait <- split(kpop_gr, kpop_gr$lodcolumn)

# Helper function to run all 4 permutation combinations, save heatmaps & generate table
run_trait_analysis <- function(trait_list, pop_name) {
  
  cat("\n--- Running Trait Analysis for:", pop_name, "---\n")
  
  # 1. Overlaps + Randomize
  cw_tr_ovlp_rand <- crosswisePermTest(
    Alist = trait_list, Blist = trait_list, ntimes = 1000,
    ranFUN = "randomizeRegions", ran.par = list(per.chromosome = TRUE),
    evFUN = "numOverlaps", genome = custom_genome, mc.cores = 4
  )
  
  # 2. Overlaps + Circular
  cw_tr_ovlp_circ <- crosswisePermTest(
    Alist = trait_list, Blist = trait_list, ntimes = 1000,
    ranFUN = "circularRandomizeRegions",
    evFUN = "numOverlaps", genome = custom_genome, mc.cores = 4
  )
  
  # --- Heatmap Generation (Safely handling S4 object) ---
  # Setting hc.method = NULL prevents hclust evaluation failures
  mat_ovlp_circ <- makeCrosswiseMatrix(
    cw_tr_ovlp_circ, 
    value = "zscore", 
    hc.method = NULL
  )
  
  # Safe PDF generation
  pdf_file <- paste0(pop_name, "_Trait_vs_Trait_Circular_Overlaps_Heatmap.pdf")
  pdf(pdf_file, width = 8, height = 8)
  
  # Try printing plot; fallback to unordered if clustering fails
  tryCatch({
    p <- plotCrosswiseMatrix(mat_ovlp_circ)
    print(p)
  }, error = function(e) {
    message("Clustered heatmap failed. Plotting without clustering...")
    p <- plotCrosswiseMatrix(mat_ovlp_circ, ord0 = FALSE)
    print(p)
  })
  
  dev.off()
  
  # --- Summary Table Compilation ---
  trait_table <- rbind(
    extract_results_table(cw_tr_ovlp_rand, paste0(pop_name, " - Overlaps (Randomize)")),
    extract_results_table(cw_tr_ovlp_circ, paste0(pop_name, " - Overlaps (Circular)"))
  )
  
  write.csv(trait_table, paste0(pop_name, "_Trait_vs_Trait_Summary.csv"), row.names = FALSE)
  return(trait_table)
}
# Run pipeline for Pop1
pop1_trait_summary <- run_trait_analysis(pop1_by_trait, "Pop1")

# Run pipeline for Pop2
pop2_trait_summary <- run_trait_analysis(pop2_by_trait, "Pop2")

run_trait_analysis_distance <- function(trait_list, pop_name) {
  
  cat("\n--- Running Trait Analysis for:", pop_name, "---\n")
  
  # 1. Overlaps + Randomize
  cw_tr_ovlp_rand <- crosswisePermTest(
    Alist = trait_list, Blist = trait_list, ntimes = 1000,
    ranFUN = "randomizeRegions", ran.par = list(per.chromosome = TRUE),
    evFUN = "meanDistance", genome = custom_genome, mc.cores = 4
  )
  
  # 2. Overlaps + Circular
  cw_tr_ovlp_circ <- crosswisePermTest(
    Alist = trait_list, Blist = trait_list, ntimes = 1000,
    ranFUN = "circularRandomizeRegions",
    evFUN = "meanDistance", genome = custom_genome, mc.cores = 4
  )
  
  # --- Heatmap Generation (Safely handling S4 object) ---
  # Setting hc.method = NULL prevents hclust evaluation failures
  mat_ovlp_circ <- makeCrosswiseMatrix(
    cw_tr_ovlp_circ, 
    value = "zscore", 
    hc.method = NULL
  )
  
  # Safe PDF generation
  pdf_file <- paste0(pop_name, "_Trait_vs_Trait_Circular_Distance_Heatmap.pdf")
  pdf(pdf_file, width = 8, height = 8)
  
  # Try printing plot; fallback to unordered if clustering fails
  tryCatch({
    p <- plotCrosswiseMatrix(mat_ovlp_circ)
    print(p)
  }, error = function(e) {
    message("Clustered heatmap failed. Plotting without clustering...")
    p <- plotCrosswiseMatrix(mat_ovlp_circ, ord0 = FALSE)
    print(p)
  })
  
  dev.off()
  
  # --- Summary Table Compilation ---
  trait_table <- rbind(
    extract_results_table(cw_tr_ovlp_rand, paste0(pop_name, " - Overlaps (Randomize)")),
    extract_results_table(cw_tr_ovlp_circ, paste0(pop_name, " - Overlaps (Circular)"))
  )
  
  write.csv(trait_table, paste0(pop_name, "_Trait_vs_Trait_Summary.csv"), row.names = FALSE)
  return(trait_table)
}
# Run pipeline for Pop1
pop1_trait_summary_distance <- run_trait_analysis_distance(pop1_by_trait, "Pop1")

# Run pipeline for Pop2
pop2_trait_summary_distance <- run_trait_analysis_distance(pop2_by_trait, "Pop2")

library(dplyr)

###Compare average pve across qtl

obj1<- peaks_spop_clean
obj2<- peaks_kpop_clean
# Summary for Object 1
avg_pve_1 <- obj1 %>%
  group_by(lodcolumn) %>%
  summarise(
    mean_PVE_1 = mean(PVE, na.rm = TRUE),
    qtl_count_1 = n(),
    .groups = "drop"
  )

# Summary for Object 2
avg_pve_2 <- obj2 %>%
  group_by(lodcolumn) %>%
  summarise(
    mean_PVE_2 = mean(PVE, na.rm = TRUE),
    qtl_count_2 = n(),
    .groups = "drop"
  )

# Combine summaries
trait_pve_comp <- inner_join(avg_pve_1, avg_pve_2, by = "lodcolumn")

# View the merged data
head(trait_pve_comp)

wilcox_res <- wilcox.test(
  trait_pve_comp$mean_PVE_1, 
  trait_pve_comp$mean_PVE_2, 
  paired = TRUE
)

print(wilcox_res)

# Define target traits
target_traits <- c("Days to flower", "stem thickness", "angle of primary branches", "shed leaves", "thick side roots",
                   "seed colour grey")

# Filter merged summary
subset_pve_comp <- trait_pve_comp %>%
  filter(lodcolumn %in% target_traits)

# Paired test on the subset
wilcox.test(subset_pve_comp$mean_PVE_1, subset_pve_comp$mean_PVE_2, paired = TRUE)

library(ggplot2)

# 1. Define your custom color palette
pop_colors <- c(
  "Progenitor" = "#E66101", 
  "NDW"        = "#0571B0"
)

# 2. Input your data as c() vectors
chromosome <- c(1, 2, 3, 4, 5, 6, 7, 8, 
                1, 2, 3, 4, 5, 6, 7, 8)         # Replace with your chromosome numbers
population <- c("Progenitor","Progenitor","Progenitor","Progenitor","Progenitor","Progenitor","Progenitor","Progenitor",
                "NDW","NDW","NDW","NDW","NDW","NDW","NDW","NDW" ) # Replace with populations
qtl_count  <- c(3,1,3,8,13,6,3,0,
                2,5,8,11,7,9,1,3)       # Replace with your QTL counts

# 3. Combine into a data frame
df <- data.frame(
  chromosome = factor(chromosome), # Converted to factor for discrete axis ticks
  population = population,
  qtl_count  = qtl_count
)

ggplot(df, aes(x = chromosome, y = qtl_count, fill = population)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = pop_colors) +
  labs(
    x = "Chromosome",
    y = "Number of QTL",
    fill = "Population",
    title = "QTL Count per Chromosome by Population"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # Set axis lines, ticks, labels, and titles to black
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    # Layout adjustments
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )


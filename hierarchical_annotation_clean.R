########################################################################
# COMPLETE HIERARCHICAL CELL-TYPE ANNOTATION PIPELINE v3.0
# CONSENSUS ANNOTATION ACROSS 5 RUNS - MODIFIED VERSION
########################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(openai)
  library(GPTCelltype)
})

## -------- CONFIGURATION ------------------------------------------------
Sys.setenv(OPENAI_API_KEY = "your-api-key-here")

TISSUE_CONTEXT <- "human PBMC peripheral blood mononuclear cells"
TENX_DIR <- "~/Downloads/pbmc3k_filtered_gene_bc_matrices.tar/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19"
N_RUNS <- 5  # Number of annotation runs

## -------- FUNCTIONS ----------------------------------------------------

create_coarse_prompt <- function(tissue_context) {
  paste0(
    "You are analyzing single-cell RNA-seq data from ", tissue_context, ". ",
    "Based on the marker genes provided, identify the MAIN CELL TYPE for each cluster.\n\n",
    "CRITICAL: Identify the actual biological cell type, NOT the genes being expressed.\n\n",
    "Example expected major cell types in immune cells are as follow, it can be different in solid tissue (this is just an example):\n",
    "- T cells (CD3D, CD3E, CD3G, CD2, CD5)\n",
    "- B cells (CD79A, CD79B, MS4A1, CD19)\n", 
    "- NK cells (NKG7, GNLY, KLRD1, KLRF1)\n",
    "- Monocytes (CD14, LYZ, S100A8, S100A9)\n",
    "- Macrophages (CD68, MSR1, MRC1)\n",
    "- Dendritic cells (FCER1A, CLEC10A)\n\n",
    "RULES:\n",
    "1. Output ONLY the cell type name (2-3 words maximum)\n",
    "2. NEVER ever annotate cells as describing based on gene expression for example this is unaccepatble annotation: 'Ribosomal protein')\n",
    "3. Focus on cell identity match with ", tissue_context, ", not cell state\n",
    "4. If unclear, use broad categories: 'Lymphocytes', 'Myeloid cells'\n"
  )
}

create_subtype_prompt <- function(parent_celltype, tissue_context) {
  paste0(
    "CRITICAL CONTEXT: ALL clusters below are ", parent_celltype, " from ", tissue_context, ".\n",
    "There are NO other cell types - ONLY ", parent_celltype, ".\n\n",
    "TASK: Identify the specific SUBTYPE of ", parent_celltype, " for each cluster.\n\n",
    "You are analyzing the markers for each cluster as a subtype of ", parent_celltype, " and based on the marker you should say which cell subtype are they", 
    "RULES:\n",
    "1. Every cluster is a known subtype of ", parent_celltype, "\n",
    "2. Use marker genes to distinguish ", parent_celltype, " subtypes\n",
    "3. Keep annotations concise (≤4 words)\n",
    "4. Must include or imply '", parent_celltype, "' in the name\n",
    "5. If thre is more than one cluster with unclear annoation, use: '", parent_celltype, " 1', '", parent_celltype, " 2', etc. NEVER use number for unique annoation\n",
    "6. DO NOT use sepeicifc symbols, use routine alphabetic symbols and text\n",
    "7. NEVER identify a different cell type - these are ALL ", parent_celltype, "\n"
  )
}

extract_top_markers <- function(seurat_obj, cluster_id, n_markers = 30) {
  tryCatch({
    markers <- FindMarkers(seurat_obj, 
                           ident.1 = cluster_id,
                           only.pos = TRUE, 
                           min.pct = 0.25, 
                           logfc.threshold = 0.25)
    
    if(nrow(markers) == 0) return(NULL)
    
    # Keep the full marker information like FindAllMarkers() returns
    top_markers <- markers %>% 
      arrange(desc(avg_log2FC)) %>% 
      head(n_markers)
    
    # Add cluster and gene columns to match FindAllMarkers() format
    top_markers$cluster <- cluster_id
    top_markers$gene <- rownames(top_markers)
    
    return(top_markers)
  }, error = function(e) {
    cat("Error finding markers for cluster", cluster_id, ":", e$message, "\n")
    return(NULL)
  })
}

## -------- ROBUST ANNOTATION ASSIGNMENT FUNCTION ----------------------

assign_annotations_robust <- function(seurat_obj, annotations, column_name = "celltype_coarse") {
  cat("🔧 Assigning", column_name, "annotations...\n")
  
  # Initialize the metadata column
  seurat_obj@meta.data[[column_name]] <- NA
  
  # Get all unique cluster IDs from Seurat object
  seurat_clusters <- unique(as.character(Idents(seurat_obj)))
  
  # Assign annotations cluster by cluster
  for(cluster_id in seurat_clusters) {
    # Find cells in this cluster
    cells_in_cluster <- which(as.character(Idents(seurat_obj)) == cluster_id)
    
    # Get annotation for this cluster
    if(cluster_id %in% names(annotations)) {
      annotation <- annotations[cluster_id]
    } else {
      annotation <- paste("Unknown", cluster_id)
      cat("⚠️ No annotation found for cluster", cluster_id, ", using:", annotation, "\n")
    }
    
    # Assign annotation to these cells
    seurat_obj@meta.data[[column_name]][cells_in_cluster] <- annotation
    
    cat("✓ Cluster", cluster_id, ":", length(cells_in_cluster), "cells →", annotation, "\n")
  }
  
  # Convert to factor
  seurat_obj@meta.data[[column_name]] <- factor(seurat_obj@meta.data[[column_name]])
  
  cat("✅ Successfully assigned", column_name, "!\n")
  return(seurat_obj)
}

## -------- SINGLE ANNOTATION RUN FUNCTION ------------------------------

run_single_annotation <- function(pbmc, run_id) {
  cat("\n🔄 ====== ANNOTATION RUN", run_id, "======\n")
  
  # Create a copy for this run
  pbmc_run <- pbmc
  
  ## -------- COARSE ANNOTATION --------------------------------------------
  
  cat("\n--- STEP 1: COARSE CELL TYPE ANNOTATION ---\n")
  
  # Find markers (add some randomness by slightly varying parameters)
  min_pct_val <- runif(1, 0.2, 0.3)  # Add some variability
  logfc_val <- runif(1, 0.2, 0.3)
  
  cat("Finding markers for all clusters (min.pct =", round(min_pct_val, 2), ")...\n")
  all_markers <- FindAllMarkers(pbmc_run, only.pos = TRUE, 
                                min.pct = min_pct_val, 
                                logfc.threshold = logfc_val)
  
  # Get coarse annotations
  cat("Getting coarse annotations from GPT...\n")
  coarse_prompt <- create_coarse_prompt(TISSUE_CONTEXT)
  
  coarse_annotations <- gptcelltype(all_markers, 
                                    model = 'gpt-4.1',
                                    tissuename = coarse_prompt)
  
  # Clean annotations
  coarse_clean <- sapply(coarse_annotations, function(x) {
    clean <- trimws(gsub("^[-\\s]*", "", x))
    if(nchar(clean) > 30) {
      clean <- substr(clean, 1, 27)
      clean <- paste0(clean, "...")
    }
    return(clean)
  })
  
  # Ensure names are character
  if(is.null(names(coarse_clean))) {
    names(coarse_clean) <- as.character(0:(length(coarse_clean)-1))
  }
  
  # Use robust assignment
  pbmc_run <- assign_annotations_robust(pbmc_run, coarse_clean, "celltype_coarse")
  
  ## -------- FINE ANNOTATION ----------------------------------------------
  
  cat("\n--- STEP 2: FINE-GRAINED HIERARCHICAL ANNOTATION ---\n")
  
  # Initialize fine annotation
  pbmc_run$celltype_fine <- as.character(pbmc_run$celltype_coarse)
  
  # Get unique coarse types
  unique_coarse_types <- unique(pbmc_run$celltype_coarse)
  cat("Processing", length(unique_coarse_types), "unique cell types\n")
  
  # Process each coarse cell type
  for(celltype in unique_coarse_types) {
    
    cat("\n━━━ PROCESSING:", celltype, "━━━\n")
    
    # Skip problematic types
    if(grepl("Ribosomal|Unknown", celltype, ignore.case = TRUE)) {
      cat("⚠️ Skipping:", celltype, "\n")
      next
    }
    
    # Get cells of this type
    cells_of_type <- which(pbmc_run$celltype_coarse == celltype)
    
    if(length(cells_of_type) < 30) {
      cat("❌ Too few cells (", length(cells_of_type), ") for", celltype, "\n")
      next
    }
    
    # Subset and re-cluster
    subset_seurat <- pbmc_run[, cells_of_type]
    
    cat("🔄 Re-clustering", length(cells_of_type), celltype, "cells...\n")
    
    # Add some randomness to clustering resolution
    resolution_val <- runif(1, 0.6, 1.0)
    
    subset_seurat <- subset_seurat %>%
      FindNeighbors(dims = 1:15) %>%
      FindClusters(resolution = resolution_val)
    
    n_subclusters <- length(unique(Idents(subset_seurat)))
    cat("Found", n_subclusters, "subclusters (resolution =", round(resolution_val, 2), ")\n")
    
    if(n_subclusters == 1) {
      cat("❌ No sub-clustering for", celltype, "\n")
      next
    }
    
    # Process each subcluster
    subtype_prompt <- create_subtype_prompt(celltype, TISSUE_CONTEXT)
    subcluster_ids <- unique(Idents(subset_seurat))
    
    for(sc_id in subcluster_ids) {
      cat("  Processing subcluster", sc_id, "...\n")
      
      # Get markers for this subcluster
      sc_markers <- extract_top_markers(subset_seurat, sc_id, n_markers = 30)
      
      if(is.null(sc_markers)) {
        cat("    ❌ No markers for subcluster", sc_id, "\n")
        next
      }
      
      # Get annotation
      tryCatch({
        cat("    🤖 Calling GPT for subcluster", sc_id, "...\n")
        sc_annotation <- gptcelltype(sc_markers, 
                                     model = 'gpt-4.1',
                                     tissuename = subtype_prompt)
        
        clean_annotation <- trimws(gsub("^[-\\s]*", "", unlist(sc_annotation)[1]))
        
        # Find cells in this subcluster
        cells_in_subcluster <- colnames(subset_seurat)[Idents(subset_seurat) == sc_id]
        
        # Update main object
        cell_indices <- which(colnames(pbmc_run) %in% cells_in_subcluster)
        pbmc_run$celltype_fine[cell_indices] <- clean_annotation
        
        cat("    ✅ Subcluster", sc_id, "→", clean_annotation, "\n")
        
      }, error = function(e) {
        cat("    ❌ Error for subcluster", sc_id, ":", e$message, "\n")
        cat("    📝 Using fallback annotation:", paste(celltype, sc_id), "\n")
        # Fallback annotation
        cell_indices <- which(colnames(pbmc_run) %in% colnames(subset_seurat)[Idents(subset_seurat) == sc_id])
        pbmc_run$celltype_fine[cell_indices] <- paste(celltype, sc_id)
      })
    }
  }
  
  # Convert to factor
  pbmc_run$celltype_fine <- factor(pbmc_run$celltype_fine)
  
  ## -------- POST-PROCESSING: MERGE NUMBERED ANNOTATIONS ----------------
  
  cat("\n--- STEP 3: POST-PROCESSING AND MERGING ---\n")
  
  # Function to merge numbered cell types
  merge_numbered_annotations <- function(annotations) {
    cat("🔄 Merging numbered annotations...\n")
    
    # Get unique annotations
    unique_annots <- unique(as.character(annotations))
    
    # Find patterns like "NK cell 1", "NK cell 2", etc.
    merge_groups <- list()
    
    for(annot in unique_annots) {
      # Remove trailing numbers and spaces to get base name
      base_name <- trimws(gsub("\\s+[0-9]+$", "", annot))
      
      # If the base name is different from original, it had a number
      if(base_name != annot && nchar(base_name) > 0) {
        # Group by base name
        if(!base_name %in% names(merge_groups)) {
          merge_groups[[base_name]] <- c()
        }
        merge_groups[[base_name]] <- c(merge_groups[[base_name]], annot)
      }
    }
    
    # Create mapping of old -> new annotations
    new_annotations <- as.character(annotations)
    
    # Apply merging
    for(base_name in names(merge_groups)) {
      numbered_variants <- merge_groups[[base_name]]
      
      if(length(numbered_variants) > 1) {
        cat("🔗 Merging:", paste(numbered_variants, collapse = ", "), "→", base_name, "\n")
        
        # Replace all numbered variants with base name
        for(variant in numbered_variants) {
          new_annotations[new_annotations == variant] <- base_name
        }
      }
    }
    
    return(new_annotations)
  }
  
  # Apply merging to fine annotations
  original_fine <- pbmc_run$celltype_fine
  merged_fine <- merge_numbered_annotations(pbmc_run$celltype_fine)
  pbmc_run$celltype_fine <- factor(merged_fine)
  
  # Show merging results
  cat("📊 MERGING SUMMARY FOR RUN", run_id, ":\n")
  cat("Before merging:", length(unique(original_fine)), "unique types\n")
  cat("After merging:", length(unique(pbmc_run$celltype_fine)), "unique types\n")
  
  # Return the annotations for this run
  return(list(
    coarse = as.character(pbmc_run$celltype_coarse),
    fine = as.character(pbmc_run$celltype_fine),
    cell_barcodes = colnames(pbmc_run)
  ))
}

## -------- DATA LOADING -------------------------------------------------

cat("=== LOADING AND PREPROCESSING DATA ===\n")
cat("Tissue context:", TISSUE_CONTEXT, "\n")

# Load data
pbmc.data <- Read10X(data.dir = TENX_DIR)
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "hierarchical_annotation", 
                           min.cells = 3, 
                           min.features = 200)

cat("Loaded:", ncol(pbmc), "cells,", nrow(pbmc), "genes\n")

# Quality control
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

cat("After filtering:", ncol(pbmc), "cells\n")

# Standard preprocessing
pbmc <- pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(object = pbmc)) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:15)

n_clusters <- length(unique(Idents(pbmc)))
cat("Initial clustering complete:", n_clusters, "clusters\n")

## -------- MULTIPLE ANNOTATION RUNS ------------------------------------

cat("\n🚀 === RUNNING", N_RUNS, "ANNOTATION ITERATIONS ===\n")

# Initialize results storage
all_results <- list()
cell_barcodes <- colnames(pbmc)

# Run annotation N_RUNS times
for(run_i in 1:N_RUNS) {
  
  # Set seed for some reproducibility within randomness
  set.seed(run_i * 12345)
  
  cat("\n" %+% paste(rep("=", 50), collapse="") %+% "\n")
  cat("🏃 STARTING RUN", run_i, "of", N_RUNS, "\n")
  cat(paste(rep("=", 50), collapse="") %+% "\n")
  
  # Run annotation
  run_result <- run_single_annotation(pbmc, run_i)
  
  # Store results
  all_results[[paste0("Run_", run_i)]] <- run_result
  
  cat("✅ COMPLETED RUN", run_i, "\n")
  cat("   - Coarse types:", length(unique(run_result$coarse)), "\n")
  cat("   - Fine types:", length(unique(run_result$fine)), "\n")
  
  # Optional: save intermediate results
  saveRDS(run_result, paste0("annotation_run_", run_i, "_", Sys.Date(), ".rds"))
}

## -------- CREATE CONSENSUS CSV -----------------------------------------

cat("\n📊 === CREATING CONSENSUS ANNOTATION CSV ===\n")

# Create consensus dataframe
consensus_df <- data.frame(Cell_Barcode = cell_barcodes)

# Add coarse annotations from each run
for(run_i in 1:N_RUNS) {
  run_name <- paste0("Run_", run_i)
  consensus_df[[paste0("Coarse_", run_name)]] <- all_results[[run_name]]$coarse
}

# Add fine annotations from each run
for(run_i in 1:N_RUNS) {
  run_name <- paste0("Run_", run_i)
  consensus_df[[paste0("Fine_", run_name)]] <- all_results[[run_name]]$fine
}

# Save consensus CSV
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
csv_filename <- paste0("consensus_annotations_", timestamp, ".csv")
write.csv(consensus_df, csv_filename, row.names = FALSE)

cat("✅ Consensus CSV created:", csv_filename, "\n")
cat("📏 Dimensions:", nrow(consensus_df), "cells ×", ncol(consensus_df), "columns\n")

## -------- SUMMARY -------------------------------------------------------

cat("\n🎯 === FINAL SUMMARY ===\n")
cat("✅ Completed", N_RUNS, "annotation runs successfully!\n")
cat("📁 Files created:\n")
cat("   - Consensus CSV:", csv_filename, "\n")
for(run_i in 1:N_RUNS) {
  cat("   - Run", run_i, "RDS: annotation_run_", run_i, "_", Sys.Date(), ".rds\n")
}

cat("\n📊 ANNOTATION DIVERSITY ACROSS RUNS:\n")
for(run_i in 1:N_RUNS) {
  run_name <- paste0("Run_", run_i)
  cat("Run", run_i, "- Coarse:", length(unique(all_results[[run_name]]$coarse)), 
      "types, Fine:", length(unique(all_results[[run_name]]$fine)), "types\n")
}

cat("\n🔍 Next step: Use the ChatGPT prompt with", csv_filename, "to determine final consensus annotations!\n") 

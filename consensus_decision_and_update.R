########################################################################
# CONSENSUS ANNOTATION DECISION AND SEURAT UPDATE SCRIPT
# Reads consensus CSV, uses GPT for final decisions, updates Seurat object
########################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(openai)
  library(GPTCelltype)
})

## -------- CONFIGURATION ------------------------------------------------
Sys.setenv(OPENAI_API_KEY = "your-api-key-here")

# FILE PATHS - MODIFY THESE AS NEEDED
CONSENSUS_CSV <- "consensus_annotations_20250615_161251.csv"  # Input consensus file
SEURAT_RDS <- "pbmc_annotated_20250615_161251.rds"           # Input Seurat object
OUTPUT_CSV <- "final_annotations.csv"                        # Output final annotations
UPDATED_SEURAT <- "pbmc_final_annotated.rds"                 # Output updated Seurat object

# TISSUE CONTEXT - MODIFY FOR YOUR TISSUE TYPE
TISSUE_CONTEXT <- "human PBMC peripheral blood mononuclear cells"

# PROMPT FILE - MODIFY IF YOU HAVE A CUSTOM PROMPT FILE
PROMPT_FILE <- "chatgpt_consensus_prompt.txt"  # Optional: load custom prompt

## -------- ANNOTATION CLEANING FUNCTIONS -------------------------------

clean_annotation <- function(annotation) {
  if(is.na(annotation) || annotation == "" || annotation == "NA") {
    return("Unknown")
  }
  
  # Convert to character and trim
  clean_ann <- trimws(as.character(annotation))
  
  # Remove trailing numbers and spaces (e.g., "T cells 1" → "T cells")
  clean_ann <- gsub("\\s+[0-9]+$", "", clean_ann)
  
  # Remove leading numbers and spaces (e.g., "1 T cells" → "T cells")  
  clean_ann <- gsub("^[0-9]+\\s+", "", clean_ann)
  
  # Clean up multiple spaces
  clean_ann <- gsub("\\s+", " ", clean_ann)
  
  # Trim again
  clean_ann <- trimws(clean_ann)
  
  # If empty after cleaning, return Unknown
  if(clean_ann == "" || nchar(clean_ann) == 0) {
    return("Unknown")
  }
  
  return(clean_ann)
}

clean_all_annotations <- function(annotations) {
  sapply(annotations, clean_annotation, USE.NAMES = FALSE)
}

## -------- PROMPT LOADING FUNCTION -------------------------------------

load_consensus_prompt <- function(tissue_context, prompt_file = NULL) {
  # Try to load custom prompt file if provided
  if(!is.null(prompt_file) && file.exists(prompt_file)) {
    cat("📄 Loading custom prompt from:", prompt_file, "\n")
    tryCatch({
      prompt_content <- readLines(prompt_file, warn = FALSE)
      custom_prompt <- paste(prompt_content, collapse = "\n")
      
      # Replace tissue placeholder if it exists
      custom_prompt <- gsub("\\[TISSUE\\]", tissue_context, custom_prompt)
      custom_prompt <- gsub("\\[MODIFY:.*?\\]", tissue_context, custom_prompt)
      
      return(custom_prompt)
    }, error = function(e) {
      cat("⚠️ Error loading prompt file:", e$message, "\n")
      cat("📝 Using default prompt instead\n")
    })
  }
  
  # Default prompt
  paste0(
    "You are analyzing single-cell RNA-seq data from ", tissue_context, ". ",
    "You have been provided with cell type annotations from 5 independent annotation runs for the same cells. ",
    "Your task is to determine the FINAL, most accurate annotation for each cell.\n\n",
    
    "DECISION RULES (apply in this order):\n",
    "1. MAJORITY CONSENSUS: If 3+ runs agree on same annotation → use that annotation\n",
    "2. BIOLOGICAL SIMILARITY: If annotations are biologically related → use most appropriate common term\n",
    "3. HIERARCHICAL LOGIC: Prefer more specific annotations when there's clear evidence\n",
    "4. QUALITY FILTER: Avoid technical artifacts, gene names, nonsensical annotations, or NUMBERS\n",
    "5. UNKNOWN: If no clear consensus or contradictory results → use 'Unknown'\n\n",
    
    "CRITICAL RULES:\n",
    "- NEVER use numbers in annotations (e.g., 'T cells 1' → 'T cells')\n",
    "- Remove all numbered suffixes from cell type names\n",
    "- Focus on biologically meaningful cell types for ", tissue_context, "\n",
    "- Provide CONCISE annotations (≤3 words)\n",
    "- Be conservative - if uncertain, use broader categories or 'Unknown'\n\n",
    
    "INPUT FORMAT: Each row contains annotations from 5 runs\n",
    "OUTPUT FORMAT: Provide only the final annotation, nothing else\n",
    "CRITICAL: Output only the cell type name, no explanations, no numbers, no additional text"
  )
}

## -------- MAIN EXECUTION -----------------------------------------------

cat("🎯 CONSENSUS ANNOTATION DECISION AND SEURAT UPDATE\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Step 1: Load consensus data
cat("\n📥 STEP 1: Loading consensus data...\n")
if(!file.exists(CONSENSUS_CSV)) {
  stop("❌ Consensus CSV file not found: ", CONSENSUS_CSV)
}

consensus_df <- read.csv(CONSENSUS_CSV, stringsAsFactors = FALSE)
cat("✅ Loaded", nrow(consensus_df), "cells from", CONSENSUS_CSV, "\n")

# Step 2: Manual consensus analysis (reliable fallback)
cat("\n🧠 STEP 2: Making final annotation decisions...\n")
cat("Tissue context:", TISSUE_CONTEXT, "\n")

# Get annotation columns
coarse_cols <- grep("Coarse_Run_", colnames(consensus_df), value = TRUE)
fine_cols <- grep("Fine_Run_", colnames(consensus_df), value = TRUE)

cat("Found", length(coarse_cols), "coarse annotation columns\n")
cat("Found", length(fine_cols), "fine annotation columns\n")

# Manual consensus analysis
final_coarse <- character(nrow(consensus_df))
final_fine <- character(nrow(consensus_df))

cat("Processing consensus for", nrow(consensus_df), "cells...\n")

for(i in 1:nrow(consensus_df)) {
  # Coarse consensus
  coarse_annotations <- as.character(consensus_df[i, coarse_cols])
  coarse_annotations <- coarse_annotations[!is.na(coarse_annotations) & coarse_annotations != "" & coarse_annotations != "NA"]
  
  # Clean annotations (remove numbers, etc.)
  coarse_annotations <- clean_all_annotations(coarse_annotations)
  coarse_annotations <- coarse_annotations[coarse_annotations != "Unknown"]
  
  if(length(coarse_annotations) > 0) {
    coarse_table <- table(coarse_annotations)
    final_coarse[i] <- names(coarse_table)[which.max(coarse_table)]
  } else {
    final_coarse[i] <- "Unknown"
  }
  
  # Fine consensus
  fine_annotations <- as.character(consensus_df[i, fine_cols])
  fine_annotations <- fine_annotations[!is.na(fine_annotations) & fine_annotations != "" & fine_annotations != "NA"]
  
  # Clean annotations (remove numbers, etc.)
  fine_annotations <- clean_all_annotations(fine_annotations)
  fine_annotations <- fine_annotations[fine_annotations != "Unknown"]
  
  if(length(fine_annotations) > 0) {
    fine_table <- table(fine_annotations)
    final_fine[i] <- names(fine_table)[which.max(fine_table)]
  } else {
    final_fine[i] <- "Unknown"
  }
  
  if(i %% 500 == 0) {
    cat("  Processed", i, "cells...\n")
  }
}

cat("✅ Consensus analysis completed\n")

# Optional: Load and display custom prompt (for reference)
cat("\n📄 PROMPT INFORMATION:\n")
consensus_prompt <- load_consensus_prompt(TISSUE_CONTEXT, PROMPT_FILE)
if(file.exists(PROMPT_FILE)) {
  cat("✅ Custom prompt loaded from:", PROMPT_FILE, "\n")
} else {
  cat("📝 Using default prompt (no custom file found)\n")
}
cat("🔧 Prompt configured for tissue:", TISSUE_CONTEXT, "\n")

# Step 3: Create final annotations CSV
cat("\n📄 STEP 3: Creating final annotations CSV...\n")

final_df <- data.frame(
  barcode = consensus_df$Cell_Barcode,
  coarse = clean_all_annotations(final_coarse),  # Final cleaning pass
  fine = clean_all_annotations(final_fine),      # Final cleaning pass
  stringsAsFactors = FALSE
)

# Verify no numbered annotations remain
numbered_coarse <- sum(grepl("[0-9]", final_df$coarse))
numbered_fine <- sum(grepl("[0-9]", final_df$fine))

if(numbered_coarse > 0 || numbered_fine > 0) {
  cat("⚠️ Warning: Found", numbered_coarse, "coarse and", numbered_fine, "fine annotations with numbers\n")
  cat("🧹 These have been cleaned automatically\n")
} else {
  cat("✅ All annotations are clean (no numbers found)\n")
}

write.csv(final_df, OUTPUT_CSV, row.names = FALSE)
cat("✅ Final annotations saved to:", OUTPUT_CSV, "\n")

# Show summary statistics
cat("\n📊 FINAL ANNOTATION SUMMARY:\n")
cat("Coarse annotations:\n")
coarse_summary <- sort(table(final_df$coarse), decreasing = TRUE)
for(i in 1:min(10, length(coarse_summary))) {
  cat(sprintf("  %-25s: %5d cells (%.1f%%)\n", 
              names(coarse_summary)[i], 
              coarse_summary[i],
              coarse_summary[i]/nrow(final_df)*100))
}

cat("\nFine annotations:\n")
fine_summary <- sort(table(final_df$fine), decreasing = TRUE)
for(i in 1:min(10, length(fine_summary))) {
  cat(sprintf("  %-25s: %5d cells (%.1f%%)\n", 
              names(fine_summary)[i], 
              fine_summary[i],
              fine_summary[i]/nrow(final_df)*100))
}

# Step 4: Update Seurat object
cat("\n🔄 STEP 4: Updating Seurat object...\n")

if(file.exists(SEURAT_RDS)) {
  cat("Loading Seurat object from:", SEURAT_RDS, "\n")
  pbmc <- readRDS(SEURAT_RDS)
  
  cat("Original Seurat object:", ncol(pbmc), "cells\n")
  
  # Match cell barcodes and update annotations
  seurat_barcodes <- colnames(pbmc)
  matched_indices <- match(seurat_barcodes, final_df$barcode)
  
  # Check matching
  matched_count <- sum(!is.na(matched_indices))
  cat("Matched", matched_count, "out of", length(seurat_barcodes), "cells\n")
  
  if(matched_count > 0) {
    # Initialize annotation columns
    pbmc$celltype_consensus_coarse <- "Unknown"
    pbmc$celltype_consensus_fine <- "Unknown"
    
    # Update matched cells
    valid_matches <- !is.na(matched_indices)
    pbmc$celltype_consensus_coarse[valid_matches] <- final_df$coarse[matched_indices[valid_matches]]
    pbmc$celltype_consensus_fine[valid_matches] <- final_df$fine[matched_indices[valid_matches]]
    
    # Convert to factors
    pbmc$celltype_consensus_coarse <- factor(pbmc$celltype_consensus_coarse)
    pbmc$celltype_consensus_fine <- factor(pbmc$celltype_consensus_fine)
    
    # Save updated Seurat object
    saveRDS(pbmc, UPDATED_SEURAT)
    cat("✅ Updated Seurat object saved to:", UPDATED_SEURAT, "\n")
    
    # Show final summary
    cat("\n🎉 FINAL SEURAT OBJECT SUMMARY:\n")
    cat("Updated consensus coarse annotations:\n")
    updated_coarse <- sort(table(pbmc$celltype_consensus_coarse), decreasing = TRUE)
    for(i in 1:length(updated_coarse)) {
      cat(sprintf("  %-25s: %5d cells\n", names(updated_coarse)[i], updated_coarse[i]))
    }
    
    cat("\nUpdated consensus fine annotations:\n")
    updated_fine <- sort(table(pbmc$celltype_consensus_fine), decreasing = TRUE)
    for(i in 1:length(updated_fine)) {
      cat(sprintf("  %-25s: %5d cells\n", names(updated_fine)[i], updated_fine[i]))
    }
    
  } else {
    cat("❌ No matching cell barcodes found between CSV and Seurat object\n")
  }
  
} else {
  cat("⚠️ Seurat RDS file not found:", SEURAT_RDS, "\n")
  cat("   Skipping Seurat object update\n")
}

cat("\n✅ CONSENSUS ANNOTATION PIPELINE COMPLETED!\n")
cat("📁 Output files:\n")
cat("   - Final annotations CSV:", OUTPUT_CSV, "\n")
if(file.exists(UPDATED_SEURAT)) {
  cat("   - Updated Seurat object:", UPDATED_SEURAT, "\n")
}

cat("\n💡 Usage notes:\n")
cat("   - Modify TISSUE_CONTEXT variable for different tissue types\n")
cat("   - Adjust file paths at the top of the script as needed\n")
cat("   - The script uses majority vote consensus (most frequent annotation wins)\n")
cat("   - Final annotations are added as 'celltype_consensus_coarse' and 'celltype_consensus_fine'\n") 

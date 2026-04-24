# ==============================================================================
# Chunked Bootstrap for swa-oBN
# Runs bootstrap in blocks of `chunk_size`, saving after each block.
# Can resume from where it left off if the job crashes.
#
# Usage (after step3 has loaded data and fit_point):
#   source("2026_03_11_bootstrap_chunked.R")
#   result <- run_chunked_bootstrap(
#     y_data, Z_data, w_data, 
#     blacklist = blacklist,
#     B_total = 500, chunk_size = 50,
#     save_dir = file.path(here::here(), "analyses", "finals", "boot_chunks")
#   )
# ==============================================================================

run_chunked_bootstrap <- function(y_data, Z_data, w_data,
                                  blacklist = NULL, whitelist = NULL,
                                  B_total = 500, chunk_size = 50,
                                  nstart = 5, maxit = 50,
                                  ic = "bic", link = "probit",
                                  save_dir = "boot_chunks",
                                  seed = 42) {
  
  # --- Setup output directory ---
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  n <- nrow(y_data)
  q <- ncol(y_data)
  var_names <- colnames(y_data)
  
  # --- Check for existing chunks (resume support) ---
  existing_files <- sort(list.files(save_dir, pattern = "^chunk_\\d+\\.rds$",
                                    full.names = TRUE))
  
  if (length(existing_files) > 0) {
    cat(sprintf("Found %d existing chunk file(s). Resuming...\n",
                length(existing_files)))
    
    # Load all existing chunks
    all_gams <- list()
    for (f in existing_files) {
      chunk <- readRDS(f)
      all_gams <- c(all_gams, chunk$boot_gams)
    }
    b_done <- length(all_gams)
    cat(sprintf("  %d / %d bootstrap resamples already completed.\n", b_done, B_total))
    
    # Recover the RNG state so new resamples are reproducible
    # We do this by burning through the first b_done resamples
    set.seed(seed)
    for (i in seq_len(b_done)) sample(n, replace = TRUE)
    
  } else {
    all_gams <- list()
    b_done <- 0
    set.seed(seed)
  }
  
  if (b_done >= B_total) {
    cat("All bootstrap resamples already completed.\n")
  } else {
    b_remaining <- B_total - b_done
    n_chunks <- ceiling(b_remaining / chunk_size)
    
    cat(sprintf("Running %d remaining resamples in %d chunk(s) of up to %d...\n",
                b_remaining, n_chunks, chunk_size))
    
    # --- Normalize weights (same as swa_obn does internally) ---
    w_norm <- as.numeric(w_data) * (n / sum(as.numeric(w_data)))
    
    for (ch in seq_len(n_chunks)) {
      b_start <- b_done + (ch - 1) * chunk_size + 1
      b_end   <- min(b_done + ch * chunk_size, B_total)
      b_this  <- b_end - b_start + 1
      
      cat(sprintf("\n=== Chunk %d/%d: resamples %d–%d ===\n",
                  ch, n_chunks, b_start, b_end))
      
      chunk_gams <- vector("list", b_this)
      chunk_start_time <- Sys.time()
      
      for (b in seq_len(b_this)) {
        b_global <- b_start + b - 1
        
        # Draw bootstrap resample (RNG state is sequential from seed)
        idx <- sample(n, replace = TRUE)
        y_b <- y_data[idx, , drop = FALSE]
        w_b <- w_norm[idx]
        Z_b <- if (!is.null(Z_data)) Z_data[idx, , drop = FALSE] else NULL
        
        t0 <- Sys.time()
        
        boot_fit <- tryCatch({
          swa_obn(
            y = y_b, Z = Z_b, weights = w_b,
            search = "greedy", ic = ic, link = link,
            blacklist = blacklist, whitelist = whitelist,
            nstart = nstart, boot = NULL,
            verbose = FALSE, maxit = maxit
          )
        }, error = function(e) {
          cat(sprintf("  ERROR on resample %d: %s\n", b_global, conditionMessage(e)))
          list(gam = matrix(NA, q, q))
        })
        
        chunk_gams[[b]] <- boot_fit$gam
        
        elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        if (b %% 10 == 0 || b == 1) {
          cat(sprintf("  Resample %d/%d (global %d) — %.1f sec\n",
                      b, b_this, b_global, elapsed))
        }
      }
      
      chunk_elapsed <- as.numeric(difftime(Sys.time(), chunk_start_time, units = "mins"))
      
      # --- Save this chunk ---
      chunk_file <- file.path(save_dir, sprintf("chunk_%03d.rds", b_start))
      saveRDS(
        list(
          b_start = b_start,
          b_end = b_end,
          boot_gams = chunk_gams,
          elapsed_minutes = chunk_elapsed,
          timestamp = Sys.time()
        ),
        file = chunk_file
      )
      
      cat(sprintf("  Saved %s (%.1f min for %d resamples)\n",
                  basename(chunk_file), chunk_elapsed, b_this))
      
      # Running estimate of time remaining
      avg_per_resample <- chunk_elapsed / b_this
      remaining <- B_total - b_end
      cat(sprintf("  Est. remaining: %.0f min (%d resamples @ %.1f min each)\n",
                  remaining * avg_per_resample, remaining, avg_per_resample))
      
      all_gams <- c(all_gams, chunk_gams)
    }
  }
  
  # --- Assemble final results ---
  cat(sprintf("\nAssembling %d bootstrap results...\n", length(all_gams)))
  
  adj_sum <- matrix(0, q, q)
  n_valid <- 0
  for (g in all_gams) {
    if (!any(is.na(g))) {
      adj_sum <- adj_sum + g
      n_valid <- n_valid + 1
    }
  }
  
  edge_probs <- adj_sum / n_valid
  colnames(edge_probs) <- var_names
  rownames(edge_probs) <- var_names
  
  n_failed <- length(all_gams) - n_valid
  if (n_failed > 0) {
    cat(sprintf("  WARNING: %d / %d resamples failed and were excluded.\n",
                n_failed, length(all_gams)))
  }
  
  # --- Save final assembled result ---
  final_file <- file.path(save_dir, "bootstrap_final.rds")
  final_result <- list(
    edge_probs = edge_probs,
    boot_gams = all_gams,
    B_total = B_total,
    B_valid = n_valid,
    B_failed = n_failed
  )
  saveRDS(final_result, file = final_file)
  cat(sprintf("Saved final result to %s\n", final_file))
  
  return(final_result)
}


# ==============================================================================
# Convenience: load and assemble from saved chunks (no re-running)
# ==============================================================================

load_bootstrap_chunks <- function(save_dir = "boot_chunks", var_names = NULL) {
  
  chunk_files <- sort(list.files(save_dir, pattern = "^chunk_\\d+\\.rds$",
                                 full.names = TRUE))
  
  if (length(chunk_files) == 0) stop("No chunk files found in ", save_dir)
  
  all_gams <- list()
  total_time <- 0
  
  for (f in chunk_files) {
    chunk <- readRDS(f)
    all_gams <- c(all_gams, chunk$boot_gams)
    total_time <- total_time + chunk$elapsed_minutes
    cat(sprintf("  Loaded %s: resamples %d–%d (%.1f min)\n",
                basename(f), chunk$b_start, chunk$b_end, chunk$elapsed_minutes))
  }
  
  cat(sprintf("Total: %d resamples, %.1f min compute time\n",
              length(all_gams), total_time))
  
  q <- nrow(all_gams[[1]])
  adj_sum <- matrix(0, q, q)
  n_valid <- 0
  for (g in all_gams) {
    if (!any(is.na(g))) {
      adj_sum <- adj_sum + g
      n_valid <- n_valid + 1
    }
  }
  
  edge_probs <- adj_sum / n_valid
  if (!is.null(var_names)) {
    colnames(edge_probs) <- var_names
    rownames(edge_probs) <- var_names
  }
  
  return(list(
    edge_probs = edge_probs,
    boot_gams = all_gams,
    B_valid = n_valid,
    B_failed = length(all_gams) - n_valid
  ))
}
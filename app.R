### PACKAGE INSTALLATION AND SET-UP ### 

# Set CRAN mirror for automated installation
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Define the list of required packages (dependencies)
packages <- c("geomorph", "mvMORPH", "phytools", "vegan", "Rphylopars", "viridis", "plotly", "shiny", "ggplot2", "dplyr")

# Check for missing packages and install them
missing_pkgs <- packages[!(packages %in% installed.packages()[, "Package"])]

# Install missing packages if any exist
if(length(missing_pkgs)) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, dependencies = TRUE)
} else {
  message("All packages are already installed.")
}

# Bioconductor specific installation for ggtree
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ggtree", quietly = TRUE)) {
  BiocManager::install("ggtree")
}

# Load all libraries
lapply(packages, library, character.only = TRUE)


### HELPER FUNCTION ###

# Convert a row-wise shape matrix into a geomorph-compatible 3D array.
reshape_to_array <- function(sim_data, n_land) {
  n_sp <- nrow(sim_data)
  # Transpose and set dimensions: [Landmarks, Dimensions, Specimens]
  array(t(sim_data), dim = c(n_land, 3, n_sp))
}



### SIMULATE THE SKULLS IN THE TIPS AND THE NODES ###

# Purpose:
# Simulate tip trait scores in PC space under BM or OU dynamics.
# BM: X(t) = X(0) + epsilon with Var[e] proportional to branch length.
# OU: dX_t = alpha(theta - X_t)dt + Sigma dW_t.

simulate_tip_scores <- function(tree, root_state, sigma_matrix, target_state = NULL, model = "BM", alpha_val = NULL, pc_labels = NULL) {
  if (model == "BM") {
    # BM: Shape changes randomly; variance increases linearly with time
    tip_scores <- mvSIM(
      tree,
      model = "BM1",
      param = list(theta = root_state, sigma = sigma_matrix)
    )
  } else {
    # OU: Shape is pulled toward a 'target_state' (theta) with strength 'alpha'
    n_pcs <- length(root_state)
    n_pcs <- length(root_state)
    alpha_matrix <- if (is.matrix(alpha_val)) alpha_val else diag(alpha_val, n_pcs)
    tip_scores <- mvSIM(
      tree,
      model = "OU1",
      param = list(theta = target_state, sigma = sigma_matrix, alpha = alpha_matrix)
    )

    # Shift all tips so the mean is anchored to the ancestor consensus (root_state)
    # rather than the stationary optimum.  This preserves relative variation while
    # keeping the root at the correct starting position in PC space
    
    root_offset <- root_state - colMeans(tip_scores)
    tip_scores <- sweep(tip_scores, 2, root_offset, "+")
  }

  # Assigns names of species to the rows and PC names to columns
  rownames(tip_scores) <- tree$tip.label
  if (!is.null(pc_labels) && length(pc_labels) == ncol(tip_scores)) {
    colnames(tip_scores) <- pc_labels
  }

  tip_scores
}

# Purpose: 
# Reconstruct internal node states with maximum likelihood using phylopars.

estimate_ancestral_scores <- function(tree, tip_scores, model = "BM") {
  trait_data <- data.frame(species = tree$tip.label, tip_scores, check.names = FALSE)
  # phylopars handles the estimation of ancestral states across the whole tree topology
  # Computes the most likely values for every internal split in the tree
  recon_fit <- phylopars(
    trait_data = trait_data,
    tree = tree,
    model = if (model == "OU") "OU" else "BM",
    pheno_error = FALSE,
    phylo_correlated = TRUE,
    pheno_correlated = TRUE
  )

  # Organize data into a single matrix combining both tips and reconstructed nodes
  reconstructed_scores <- recon_fit$anc_recon
  n_tips <- Ntip(tree)
  full_history <- matrix(NA, nrow = n_tips + tree$Nnode, ncol = ncol(reconstructed_scores))
  rownames(full_history) <- as.character(seq_len(n_tips + tree$Nnode))
  colnames(full_history) <- colnames(reconstructed_scores)
  # Fills the matrix with tips
  full_history[as.character(seq_len(n_tips)), ] <- reconstructed_scores[tree$tip.label, , drop = FALSE]
  internal_node_ids <- as.character((n_tips + 1):(n_tips + tree$Nnode))
  # and then fills it with internal nodes
  full_history[internal_node_ids, ] <- reconstructed_scores[internal_node_ids, , drop = FALSE]

  list(
    fit = recon_fit,
    full_history = full_history
  )
}


## Purpose: 
# Run simulation and ancestral reconstruction and return full history.

simulate_and_reconstruct_history <- function(tree, root_state, sigma_matrix, target_state = NULL, model = "BM", alpha_val = NULL, pc_labels = NULL) {
  tip_scores <- simulate_tip_scores(
    tree,
    root_state,
    sigma_matrix,
    target_state = target_state,
    model = model,
    alpha_val = alpha_val,
    pc_labels = pc_labels
  )

  ancestral_estimate <- estimate_ancestral_scores(tree, tip_scores, model = model)
  ancestral_estimate$full_history[as.character(Ntip(tree) + 1), ] <- root_state

  list(
    tip_scores = tip_scores,
    reconstruction = ancestral_estimate$fit,
    full_history = ancestral_estimate$full_history
  )
}


### GET LANDMARK CONFIGURATIONS AT SPECIFIC TIMESTAMPS ###

# Purpose: 
# Extract branch-wise shape states at one requested geological time.
# Sample branch states at a geological time slice (Ma) by interpolating
# between parent and child node states in PC space.

sample_tree_at_time <- function(tree, full_data, target_time, alpha_val = NULL, model = NULL) {
  edges <- tree$edge
  heights <- nodeHeights(tree) # Start/End times for every branch

  slice_coords <- matrix(NA, nrow = 0, ncol = ncol(full_data))
  lineage_ids <- c()

  for (e in seq_len(nrow(edges))) {
    start_h <- heights[e, 1]
    end_h <- heights[e, 2]

    # Check if the branch exists at the target_time
    if (target_time >= start_h && target_time <= end_h) {
      parent_node <- edges[e, 1]
      child_node <- edges[e, 2]
      p_coords <- full_data[as.character(parent_node), ]
      c_coords <- full_data[as.character(child_node), ]
      dt_slice <- target_time - start_h
      branch_duration <- end_h - start_h

      if (model == "OU") {
        # OU interpolation: nonlinear pull toward the endpoint state.
        local_pull <- if (branch_duration == 0) 1 else (1 - exp(-alpha_val * dt_slice)) / (1 - exp(-alpha_val * branch_duration))
        interp <- p_coords + local_pull * (c_coords - p_coords)
      } else {
        # BM interpolation: linear interpolation along the branch segment.
        prop <- dt_slice / branch_duration
        interp <- p_coords + prop * (c_coords - p_coords)
      }

      slice_coords <- rbind(slice_coords, interp)
      lineage_ids <- c(lineage_ids, child_node)
    }
  }

  rownames(slice_coords) <- lineage_ids
  return(slice_coords)
}


## Purpose: 
# Generate time-ordered snapshot matrices across all requested intervals.

generate_all_snapshots <- function(tree, intervals, full_history, alpha_val = NULL, model = NULL) {
  snapshot_list <- list()
  for (i in seq_along(intervals)) {
    target_time <- intervals[i]
    snapshot_list[[paste0("Time_", target_time)]] <- sample_tree_at_time(
      tree,
      full_history,
      target_time,
      alpha_val = alpha_val,
      model = model
    )
    message(paste("Processed time slice:", target_time, "Ma"))
  }
  return(snapshot_list)
}

## Purpose: 
# After back-projection, align all temporal snapshots again so
# any coordinate-space noise is removed before diagnostics or morphospace plots.
# Apply GPA jointly to all snapshots and return aligned snapshot matrices.

perform_final_gpa_on_snapshots <- function(snapshot_list, n_land) {
  snapshot_sizes <- vapply(snapshot_list, nrow, integer(1))
  all_snapshots <- do.call(rbind, snapshot_list)
  all_array <- array(t(all_snapshots), dim = c(n_land, 3, nrow(all_snapshots)))
  gpa_result <- gpagen(all_array, print.progress = FALSE)
  aligned_matrix <- two.d.array(gpa_result$coords)

  aligned_snapshots <- list()
  start_idx <- 1
  for (snapshot_name in names(snapshot_list)) {
    end_idx <- start_idx + snapshot_sizes[[snapshot_name]] - 1
    aligned_snapshots[[snapshot_name]] <- aligned_matrix[start_idx:end_idx, , drop = FALSE]
    rownames(aligned_snapshots[[snapshot_name]]) <- rownames(snapshot_list[[snapshot_name]])
    start_idx <- end_idx + 1
  }

  list(
    aligned_snapshots = aligned_snapshots,
    gpa_result = gpa_result
  )
}

## Purpose: 
# GPA-align shapes and compute PCA coordinates for diagnostics.

align_shape_matrix_for_diagnostics <- function(shape_matrix, n_land, target_config = NULL) {
  all_points_matrix <- shape_matrix
  target_index <- NULL

  if (!is.null(target_config)) {
    all_points_matrix <- rbind(all_points_matrix, as.numeric(target_config))
    target_index <- nrow(all_points_matrix)
  }

  all_array <- array(t(all_points_matrix), dim = c(n_land, 3, nrow(all_points_matrix)))
  gpa_results <- gpagen(all_array, print.progress = FALSE)
  aligned_2d <- two.d.array(gpa_results$coords)
  pca_res <- prcomp(aligned_2d)

  list(
    aligned_matrix = aligned_2d,
    pca = pca_res,
    target_index = target_index
  )
}

## Purpose: 
# Plot tree-connected node/tip trajectories in PCA morphospace.
# show the reconstructed tips and internal nodes as a tree in
# morphospace rather than only as disconnected snapshot point clouds.

plot_phylomorphospace <- function(tree, full_history, n_land, title = "Phylomorphospace", target_config = NULL) {
  diag_data <- align_shape_matrix_for_diagnostics(full_history, n_land, target_config = target_config)
  pca_scores <- diag_data$pca$x[seq_len(nrow(full_history)), 1:2, drop = FALSE]

  plot(
    pca_scores[, 1],
    pca_scores[, 2],
    type = "n",
    xlab = paste0("PC1 (", round(summary(diag_data$pca)$importance[2, 1] * 100, 1), "%)"),
    ylab = paste0("PC2 (", round(summary(diag_data$pca)$importance[2, 2] * 100, 1), "%)"),
    main = title
  )

  for (edge_index in seq_len(nrow(tree$edge))) {
    parent_node <- tree$edge[edge_index, 1]
    child_node <- tree$edge[edge_index, 2]
    segments(
      pca_scores[parent_node, 1], pca_scores[parent_node, 2],
      pca_scores[child_node, 1], pca_scores[child_node, 2],
      col = rgb(0.3, 0.3, 0.3, 0.6), lwd = 1.2
    )
  }

  tip_ids <- seq_len(Ntip(tree))
  node_ids <- (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)
  points(pca_scores[tip_ids, 1], pca_scores[tip_ids, 2], pch = 16, col = "black", cex = 0.8)
  points(pca_scores[node_ids, 1], pca_scores[node_ids, 2], pch = 21, bg = "steelblue", col = "white", cex = 1)

  if (!is.null(diag_data$target_index)) {
    target_pca <- diag_data$pca$x[diag_data$target_index, 1:2, drop = FALSE]
    points(target_pca[1, 1], target_pca[1, 2], col = "red", pch = 8, cex = 1.2, lwd = 2)
    text(target_pca[1, 1], target_pca[1, 2], labels = "optimum", pos = 3, col = "red", cex = 0.6)
  }
}

## Purpose: 
# Build an interactive 3D landmark + wireframe Plotly object.

build_wireframe_plotly <- function(coords, wireframe_links, title_text, point_color = "steelblue") {
  fig <- plot_ly(type = "scatter3d") |>
    add_trace(
      x = coords[, 1],
      y = coords[, 2],
      z = coords[, 3],
      mode = "markers+text",
      type = "scatter3d",
      marker = list(size = 3, color = point_color),
      text = as.character(seq_len(nrow(coords))),
      textposition = "top center",
      textfont = list(size = 10, color = "black"),
      hovertemplate = "Landmark %{text}<extra></extra>"
    )

  if (!is.null(wireframe_links) && nrow(wireframe_links) > 0) {
    for (link_index in seq_len(nrow(wireframe_links))) {
      link_pair <- wireframe_links[link_index, ]
      fig <- fig |>
        add_trace(
          x = coords[link_pair, 1],
          y = coords[link_pair, 2],
          z = coords[link_pair, 3],
          type = "scatter3d",
          mode = "lines",
          line = list(color = "rgba(40, 40, 40, 0.7)", width = 4),
          hoverinfo = "skip",
          showlegend = FALSE
        )
    }
  }

  plotly::layout(
    fig,
    title = list(text = title_text),
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z"),
      aspectmode = "data"
    ),
    margin = list(l = 0, r = 0, b = 0, t = 40)
  )
}

# Purpose: 
# Visualize temporal dispersion of simulated shapes in PCA morphospace.
# Perform GPA -> PCA on all snapshots and render the temporal morphospace.

analyze_rigorous_morphospace <- function(snapshot_list, n_land, title = "Shape Space Expansion", target_config = NULL, start_time_ma = 1, end_time_ma = 50) {
  all_points_matrix <- do.call(rbind, snapshot_list)

  if (!is.null(target_config)) {
    target_vector <- as.numeric(target_config)
    all_points_matrix <- rbind(all_points_matrix, target_vector)
  }

  n_total_observations <- nrow(all_points_matrix)
  all_array <- array(t(all_points_matrix), dim = c(n_land, 3, n_total_observations))
  message("Performing Procrustes Alignment...")
  gpa_results <- gpagen(all_array, print.progress = FALSE)
  aligned_2d <- two.d.array(gpa_results$coords)
  pca_res <- prcomp(aligned_2d)

  if (!is.null(target_config)) {
    target_pca <- pca_res$x[n_total_observations, , drop = FALSE]
  }

  n_steps <- length(snapshot_list)
  colors_vec <- rev(viridis(n_steps))
  par(mar = c(5, 4, 4, 7))

  plot(
    pca_res$x[, 1],
    pca_res$x[, 2],
    type = "n",
    xlab = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
    ylab = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)"),
    main = title
  )

  curr_row <- 1
  for (i in seq_along(snapshot_list)) {
    n_taxa <- nrow(snapshot_list[[i]])
    rows <- curr_row:(curr_row + n_taxa - 1)
    points(pca_res$x[rows, 1], pca_res$x[rows, 2], col = colors_vec[i], pch = 16, cex = 0.8)
    curr_row <- curr_row + n_taxa
  }

  if (!is.null(target_config)) {
    points(target_pca[1, 1], target_pca[1, 2], col = "red", pch = 8, cex = 1, lwd = 2)
    text(target_pca[1, 1], target_pca[1, 2], labels = "optimum", pos = 3, col = "red", cex = 0.4, font = 2)
  }

  u <- par("usr")
  bar_width <- (u[2] - u[1]) * 0.05
  left_edge <- u[2] + (u[2] - u[1]) * 0.02
  right_edge <- left_edge + bar_width
  y_steps <- seq(u[3], u[4], length.out = n_steps + 1)
  for (i in seq_len(n_steps)) {
    rect(left_edge, y_steps[i], right_edge, y_steps[i + 1], col = colors_vec[i], border = NA, xpd = TRUE)
  }
  rect(left_edge, u[3], right_edge, u[4], border = "black", xpd = TRUE)

  xr <- (max(pca_res$x[,1]) + (max(pca_res$x[,1]) - min(pca_res$x[,1])) * 0.15)
  yb <- min(pca_res$x[,2])
  yt <- max(pca_res$x[,2])
  xl <- (max(pca_res$x[,1]) + (max(pca_res$x[,1]) - min(pca_res$x[,1])) * 0.05)
  text(x = xr, y = yb, labels = paste0(start_time_ma, " Ma"), pos = 4, cex = 0.8, xpd = TRUE)
  text(x = xr, y = yt, labels = paste0(end_time_ma, " Ma"), pos = 4, cex = 0.8, xpd = TRUE)
  text(x = (xl + xr) / 2, y = yt + (yt - yb) * 0.03, labels = "Time", pos = 3, cex = 0.9, font = 2, xpd = TRUE)
}

### RUN THE SCRIPT  ###


## Purpose: 
# Read Morphologika input and return coordinates plus metadata.

load_morphologika_data <- function(file_path) {
  all_skulls_local <- read.morphologika(file_path)
  coords_data_local <- all_skulls_local$coords

  reference_lines <- readLines(file_path, warn = FALSE)
  trimmed_reference_lines <- trimws(reference_lines)
  names_start <- match("[names]", trimmed_reference_lines)
  labels_start <- match("[labels]", trimmed_reference_lines)

  name_lines <- trimmed_reference_lines[(names_start + 1):(labels_start - 1)]
  name_lines <- name_lines[name_lines != ""]
  specimen_taxa_local <- vapply(name_lines, function(line) {
    parts <- strsplit(line, "\\s+")[[1]]
    paste(parts[1:min(2, length(parts))], collapse = " ")
  }, character(1))

  list(
    all_skulls = all_skulls_local,
    coords_data = coords_data_local,
    specimen_taxa = specimen_taxa_local,
    unique_taxa = sort(unique(specimen_taxa_local)),
    n_land = dim(coords_data_local)[1],
    p_dim = dim(coords_data_local)[2],
    n_vars = dim(coords_data_local)[1] * dim(coords_data_local)[2]
  )
}

## Purpose: 
#Compute consensus landmark configuration for one selected taxon.
# for a user-selected taxon, align its specimens with gpagen
# and use the resulting consensus as the representative shape.

compute_taxon_consensus <- function(morph_data, taxon_name) {
  selected_idx <- which(morph_data$specimen_taxa == taxon_name)
  gpagen(morph_data$coords_data[, , selected_idx, drop = FALSE], print.progress = FALSE)$consensus
}

## Purpose: 
# Apply optional node-age overrides and rebuild a calibrated branch-time tree.

build_calibrated_tree <- function(tree, calibration_text = NULL) {
  node_ids <- as.character((Ntip(tree) + 1):(Ntip(tree) + tree$Nnode))
  default_ages <- branching.times(tree)
  node_ages <- as.numeric(default_ages[node_ids])
  names(node_ages) <- node_ids

  overrides <- numeric(0)
  if (!is.null(calibration_text) && nzchar(trimws(calibration_text))) {
    lines <- trimws(unlist(strsplit(calibration_text, "\n")))
    lines <- lines[lines != ""]
    for (line in lines) {
      parts <- trimws(unlist(strsplit(line, "[:=,]")))
      parts <- parts[parts != ""]
      if (length(parts) == 2) {
        node_id <- as.character(as.integer(parts[1]))
        node_ages[node_id] <- as.numeric(parts[2])
        overrides[node_id] <- as.numeric(parts[2])
      }
    }
  }

  calibrated_tree <- compute.brtime(tree, method = node_ages, force.positive = FALSE)

  list(
    tree = calibrated_tree,
    node_ages = node_ages,
    applied_overrides = overrides
  )
}

## Purpose: 
# Build PCA simulation space, module assignments, and Sigma structure.

build_pc_simulation_space <- function(morph_data, ancestor_taxon, target_taxon, pc_mode, variance_threshold, requested_pcs, taxa_limit, module_text, sigma_diag, sigma_offdiag) {
  aligned_coords <- gpagen(morph_data$coords_data, print.progress = FALSE)$coords # Run GPA
  aligned_matrix <- two.d.array(aligned_coords) # 3D coordinates to 2D coordinates
  pca_result <- prcomp(aligned_matrix, center = TRUE, scale. = FALSE) # Run PCA

  explained_variance <- (pca_result$sdev ^ 2) / sum(pca_result$sdev ^ 2) # calculates how much information about the shape each PC has
  cumulative_variance <- cumsum(explained_variance)
  threshold_pcs <- which(cumulative_variance >= variance_threshold)[1] # number of PCs to reach the variance threshold
  max_pcs <- min(ncol(aligned_matrix), nrow(aligned_matrix) - 1, taxa_limit - 1)
  retained_pcs <- if (pc_mode == "manual") requested_pcs else threshold_pcs # checks between manual PCs input or the percentage threshold
  if (is.null(retained_pcs) || is.na(retained_pcs)) retained_pcs <- threshold_pcs
  retained_pcs <- min(retained_pcs, max_pcs)
  retained_indices <- seq_len(retained_pcs)
  retained_rotation <- pca_result$rotation[, retained_indices, drop = FALSE]

  
  # Purpose:
  # Build biologically motivated modules from landmark sets.
  # PCs are later grouped by dominant module loadings to induce block structure in Sigma.
  
  module_definitions <- if (is.null(module_text) || !nzchar(trimws(module_text))) {
    list(Global = seq_len(morph_data$n_land))
  } else {
    lines <- trimws(unlist(strsplit(module_text, "\n")))
    lines <- lines[lines != ""]
    modules <- list()
    assigned_landmarks <- integer(0)
    for (line in lines) {
      parts <- strsplit(line, ":", fixed = TRUE)[[1]]
      module_name <- trimws(parts[1])
      landmark_tokens <- trimws(unlist(strsplit(parts[2], "[,[:space:]]+")))
      landmark_tokens <- landmark_tokens[landmark_tokens != ""]
      landmark_ids <- suppressWarnings(as.integer(landmark_tokens))
      modules[[module_name]] <- sort(unique(landmark_ids))
      assigned_landmarks <- c(assigned_landmarks, landmark_ids)
    }
    unassigned <- setdiff(seq_len(morph_data$n_land), assigned_landmarks)
    if (length(unassigned) > 0) modules[["Unassigned"]] <- unassigned
    modules
  }

  module_names <- names(module_definitions)
  module_scores <- matrix(0, nrow = retained_pcs, ncol = length(module_definitions))
  for (module_index in seq_along(module_definitions)) {
    coord_ids <- unlist(lapply(module_definitions[[module_index]], function(landmark_id) {
      ((landmark_id - 1) * 3 + 1):(landmark_id * 3)
    }))
    module_scores[, module_index] <- colSums(abs(retained_rotation[coord_ids, , drop = FALSE])) # how much each landmark contributes to each PC
  }

  dominant_modules <- max.col(module_scores, ties.method = "first") # Assigns each PC to the module that uses it the most
  pc_order <- order(dominant_modules, seq_len(retained_pcs))
  ordered_rotation <- retained_rotation[, pc_order, drop = FALSE]
  ordered_pc_indices <- retained_indices[pc_order]
  ordered_module_names <- module_names[dominant_modules[pc_order]]

  # Sigma is diagonal by default; within-module off-diagonal entries encode correlated evolution.
  sigma_matrix <- diag(sigma_diag, retained_pcs) # Creates a matrix where every trait evolves at a baseline rate
  for (i in seq_len(retained_pcs)) {
    for (j in seq_len(retained_pcs)) {
      if (i != j && ordered_module_names[i] == ordered_module_names[j]) sigma_matrix[i, j] <- sigma_offdiag # e.g., if PC1 and PC2 both belong to the "Snout", this adds a covariance term
    }
  }
  pc_labels <- paste0("PC", ordered_pc_indices)
  rownames(sigma_matrix) <- pc_labels
  colnames(sigma_matrix) <- pc_labels

  ancestor_idx <- which(morph_data$specimen_taxa == ancestor_taxon)
  target_idx <- which(morph_data$specimen_taxa == target_taxon)
  ancestor_consensus <- apply(aligned_coords[, , ancestor_idx, drop = FALSE], c(1, 2), mean)
  target_consensus <- apply(aligned_coords[, , target_idx, drop = FALSE], c(1, 2), mean)
  ancestor_vector <- as.vector(t(ancestor_consensus))
  target_vector <- as.vector(t(target_consensus))

  # projects their 3D coordinates into the PC space
  pca_model <- list(
    aligned_coords = aligned_coords,
    aligned_matrix = aligned_matrix,
    pca_result = pca_result,
    explained_variance = explained_variance,
    cumulative_variance = cumulative_variance,
    threshold_pcs = threshold_pcs,
    max_pcs = max_pcs,
    retained_pcs = retained_pcs,
    retained_indices = retained_indices,
    center = pca_result$center,
    retained_rotation = retained_rotation
  )

  sigma_info <- list(
    sigma_matrix = sigma_matrix,
    ordered_rotation = ordered_rotation,
    ordered_pc_indices = ordered_pc_indices,
    ordered_module_names = ordered_module_names,
    pc_labels = pc_labels
  )

# Returns root/target states in PC space plus back-projection objects.
# rotation maps from retained PC space back to landmark coordinate space.
  
  list(
    root_state = as.numeric((ancestor_vector - pca_model$center) %*% sigma_info$ordered_rotation),
    theta_state = as.numeric((target_vector - pca_model$center) %*% sigma_info$ordered_rotation),
    ancestor_vector = ancestor_vector,
    target_vector = target_vector,
    center = pca_model$center,
    rotation = sigma_info$ordered_rotation,
    sigma_matrix = sigma_info$sigma_matrix,
    pca_model = pca_model,
    module_definitions = module_definitions,
    sigma_info = sigma_info
  )
}

## Purpose: 
# Format calibrated node ages into a readable reference table.

format_node_age_reference <- function(tree, node_ages) {
  node_ids <- as.integer(names(node_ages))
  parent_lookup <- setNames(tree$edge[, 1], tree$edge[, 2])
  child_counts <- table(tree$edge[, 1])

  lines <- vapply(node_ids, function(node_id) {
    parent_id <- parent_lookup[as.character(node_id)]
    descendant_count <- if (as.character(node_id) %in% names(child_counts)) child_counts[[as.character(node_id)]] else 0
    parent_label <- if (is.na(parent_id)) "root" else as.character(parent_id)
    paste0(
      "Node ", node_id,
      " | age = ", round(node_ages[as.character(node_id)], 3),
      " | parent = ", parent_label,
      " | child branches = ", descendant_count
    )
  }, character(1))

  paste(lines, collapse = "\n")
}

# The app starts empty and is populated after the user uploads a Morphologika file.
unique_taxa <- character(0)
default_ancestor_taxon <- NULL


### INTERACTIVE PLOT WITH NODES, TIPS, AND TIMESTAMP INTERSECTIONS ###

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      h4("Simulation Parameters"),
      fileInput("morph_file", "Morphologika file (.txt):", accept = ".txt"),
      numericInput("n_tips", "Tree tip count (n):", value = 30, min = 2, step = 1),
      numericInput("tree_scale", "Tree scale:", value = 50, min = 0.01, step = 1),
      numericInput("interval_by", "Interval step (by):", value = 1, min = 0.01, step = 0.1),
      selectInput("pc_mode", "PC retention mode:", choices = c("Variance threshold" = "threshold", "Manual count" = "manual"), selected = "threshold"),
      sliderInput("pc_variance", "Variance retained:", min = 0.95, max = 0.99, value = 0.95, step = 0.01),
      numericInput("pc_count", "PC count (manual):", value = 5, min = 1, step = 1),
      textAreaInput(
        "landmark_modules",
        "Landmark modules:",
        value = "",
        rows = 6,
        placeholder = "Snout: 1,2,3,4\nBraincase: 5,6,7,8"
      ),
      helpText("Define landmark modules as 'ModuleName: landmark_ids'. PCs are assigned to modules by their landmark loadings."),
      textAreaInput(
        "node_age_overrides",
        "Internal node age overrides:",
        value = "",
        rows = 6,
        placeholder = "Example:\n31 = 45\n32 = 28.5"
      ),
      helpText("Use 'node_id = age' to calibrate internal nodes after the topology is generated."),
      selectInput("ancestor_taxon", "Ancestor taxon:", choices = unique_taxa, selected = default_ancestor_taxon),
      selectInput("target_taxon", "Target taxon:", choices = unique_taxa, selected = default_ancestor_taxon),
      selectInput("evol_model", "Evolution model:", choices = c("OU", "BM"), selected = "OU"),
      numericInput("alpha_val", "Alpha (OU):", value = 0.3, min = 0.001, step = 0.05),
      numericInput("sigma_offdiag", "Sigma off-diagonal:", value = 0.02, min = 0, step = 0.001),
      numericInput("sigma_diag", "Sigma diagonal:", value = 0.04, min = 0, step = 0.001),
      actionButton("run_sim", "Run Simulation")
    ),
    mainPanel(
      fluidRow(
        column(7,
               h3("Tree Explorer"),
               plotOutput("tree_plot", click = "tree_click", height = "800px")
        ),
        column(5,
               tabsetPanel(
                 tabPanel("Previews", plotlyOutput("ancestor_consensus_preview", height = "400px"), plotlyOutput("ancestor_wireframe_preview", height = "400px"), hr(), plotlyOutput("target_consensus_preview", height = "400px"), plotlyOutput("target_wireframe_preview", height = "400px")),
                 tabPanel("Tree", plotlyOutput("skull_3d", height = "500px"), verbatimTextOutput("selection_info"), verbatimTextOutput("node_age_reference")),
                 tabPanel("Diagnostics", verbatimTextOutput("pca_info"), verbatimTextOutput("module_info"), verbatimTextOutput("mantel_info"), plotOutput("shape_time_plot_bm", height = "360px"), plotOutput("shape_time_plot_ou", height = "360px"), plotOutput("morphospace_plot_bm", height = "480px"), plotOutput("morphospace_plot_ou", height = "480px"))
               )
        )
      )
    )
  )
)

## Purpose: 
# Define Shiny callbacks that execute and visualize the simulation pipeline.

server <- function(input, output, session) {
  # Pipeline stage 1: persistent reactive container for tree topology.
  # Updated by topology controls; consumed when the simulation is triggered.
  tree_topology <- reactiveVal(NULL)

  # Pipeline stage 2 (callback): load and parse morphometric input file.
  # Trigger: input$morph_file changes.
  # Output: coordinates, taxa labels, and metadata used in all downstream stages.
  morph_data_reactive <- reactive({
    req(input$morph_file)
    selected_file <- input$morph_file$datapath
    load_morphologika_data(selected_file)
  })

  # Pipeline stage 3 (callback): generate/update random tree topology.
  # Trigger: n_tips or tree_scale changes.
  observeEvent(list(input$n_tips, input$tree_scale), {
    req(input$n_tips, input$tree_scale)
    tree_topology(pbtree(n = as.integer(input$n_tips), scale = as.numeric(input$tree_scale)))
  }, ignoreInit = FALSE)

  # Pipeline stage 4 (callback): keep ancestor/target selectors synchronized
  # with taxa available in the uploaded Morphologika file.
  # Trigger: morph_data_reactive() changes.
  observeEvent(morph_data_reactive(), {
    md <- morph_data_reactive()
    updateSelectInput(session, "ancestor_taxon", choices = md$unique_taxa, selected = md$unique_taxa[1])
    updateSelectInput(session, "target_taxon", choices = md$unique_taxa, selected = md$unique_taxa[1])
  }, ignoreInit = FALSE)

  output$ancestor_consensus_preview <- renderPlotly({
    req(input$ancestor_taxon)
    md <- morph_data_reactive()
    consensus <- compute_taxon_consensus(md, input$ancestor_taxon)
    
    plot_ly(x = ~consensus[,1], y = ~consensus[,2], 
            text = ~seq_len(nrow(consensus)),
            type = 'scatter', mode = 'markers+text',
            textposition = "top center",
            marker = list(color = 'black', size = 8)) %>%
      layout(title = paste("Consensus:", input$ancestor_taxon),
             xaxis = list(title = "X", scaleanchor = "y", scaleratio = 1),
             yaxis = list(title = "Y"))
  })

  output$ancestor_wireframe_preview <- renderPlotly({
    req(input$ancestor_taxon)
    md <- morph_data_reactive()
    consensus <- compute_taxon_consensus(md, input$ancestor_taxon)
    wireframe_links <- matrix(md$all_skulls$wireframe, ncol = 2)
    
    p <- plot_ly()
    
    # Add Wireframe segments
    if (nrow(wireframe_links) > 0) {
      for (i in seq_len(nrow(wireframe_links))) {
        p <- p %>% add_segments(
          x = consensus[wireframe_links[i,1], 1], xend = consensus[wireframe_links[i,2], 1],
          y = consensus[wireframe_links[i,1], 2], yend = consensus[wireframe_links[i,2], 2],
          line = list(color = 'blue', width = 1.5), showlegend = FALSE, hoverinfo = "none"
        )
      }
    }
    
    # Add Points and Numbers
    p %>% add_trace(x = consensus[,1], y = consensus[,2], 
                    text = seq_len(nrow(consensus)),
                    type = 'scatter', mode = 'markers+text',
                    textposition = "top center",
                    marker = list(color = 'red', size = 7, line = list(color = 'black', width = 1))) %>%
      layout(title = paste("2D Wireframe:", input$ancestor_taxon),
             xaxis = list(title = "X", scaleanchor = "y", scaleratio = 1),
             yaxis = list(title = "Y"))
  })

  output$target_consensus_preview <- renderPlotly({
    req(input$target_taxon)
    md <- morph_data_reactive()
    consensus <- compute_taxon_consensus(md, input$target_taxon)
    
    plot_ly(x = ~consensus[,1], y = ~consensus[,2], 
            text = ~seq_len(nrow(consensus)),
            type = 'scatter', mode = 'markers+text',
            textposition = "top center",
            marker = list(color = 'black', size = 8)) %>%
      layout(title = paste("Consensus:", input$target_taxon),
             xaxis = list(title = "X", scaleanchor = "y", scaleratio = 1),
             yaxis = list(title = "Y"))
  })

  output$target_wireframe_preview <- renderPlotly({
    req(input$target_taxon)
    md <- morph_data_reactive()
    consensus <- compute_taxon_consensus(md, input$target_taxon)
    wireframe_links <- matrix(md$all_skulls$wireframe, ncol = 2)
    
    p <- plot_ly()
    
    # Add Wireframe segments
    if (nrow(wireframe_links) > 0) {
      for (i in seq_len(nrow(wireframe_links))) {
        p <- p %>% add_segments(
          x = consensus[wireframe_links[i,1], 1], xend = consensus[wireframe_links[i,2], 1],
          y = consensus[wireframe_links[i,1], 2], yend = consensus[wireframe_links[i,2], 2],
          line = list(color = 'darkgreen', width = 1.5), showlegend = FALSE, hoverinfo = "none"
        )
      }
    }
    
    # Add Points and Numbers
    p %>% add_trace(x = consensus[,1], y = consensus[,2], 
                    text = seq_len(nrow(consensus)),
                    type = 'scatter', mode = 'markers+text',
                    textposition = "top center",
                    marker = list(color = 'orange', size = 7, line = list(color = 'black', width = 1))) %>%
      layout(title = paste("2D Wireframe:", input$target_taxon),
             xaxis = list(title = "X", scaleanchor = "y", scaleratio = 1),
             yaxis = list(title = "Y"))
  })

  # Core pipeline callback:
  # Trigger: input$run_sim button click.
  # Stages inside this callback:
  # 1) Build PCA simulation space and Sigma structure from modules.
  # 2) Calibrate tree node ages.
  # 3) Simulate/reconstruct histories in PC space (active model + BM + OU controls).
  # 4) Slice trajectories through time and back-project to landmark coordinates.
  # 5) Align snapshots (GPA) and compute Mantel diagnostics.
  # Output: single state object consumed by all diagnostic/tree/3D outputs.
  sim_state <- eventReactive(input$run_sim, {
    req(input$n_tips, input$tree_scale, input$interval_by, input$ancestor_taxon, input$target_taxon, input$evol_model, input$sigma_offdiag, input$sigma_diag)
    md <- morph_data_reactive()

    pca_space_local <- build_pc_simulation_space(
      md,
      input$ancestor_taxon,
      input$target_taxon,
      pc_mode = input$pc_mode,
      variance_threshold = as.numeric(input$pc_variance),
      requested_pcs = as.integer(input$pc_count),
      taxa_limit = as.integer(input$n_tips),
      module_text = input$landmark_modules,
      sigma_diag = as.numeric(input$sigma_diag),
      sigma_offdiag = as.numeric(input$sigma_offdiag)
    )

    root_state_local <- pca_space_local$root_state
    target_local <- pca_space_local$theta_state
    sigma_matrix_local <- pca_space_local$sigma_matrix

    base_tree_local <- tree_topology()
    req(base_tree_local)
    calibrated_tree_info <- build_calibrated_tree(base_tree_local, input$node_age_overrides)
    tree_local <- calibrated_tree_info$tree
    tree_depth_local <- max(nodeHeights(tree_local))
    intervals_local <- seq(1, tree_depth_local, by = as.numeric(input$interval_by))
    if (length(intervals_local) == 0 || tail(intervals_local, 1) < tree_depth_local) {
      intervals_local <- sort(unique(c(intervals_local, tree_depth_local)))
    }

    active_history_pc_local <- simulate_and_reconstruct_history(
      tree_local,
      root_state_local,
      sigma_matrix_local,
      target_state = target_local,
      model = input$evol_model,
      alpha_val = as.numeric(input$alpha_val),
      pc_labels = pca_space_local$sigma_info$pc_labels
    )
    active_snapshots_pc_local <- generate_all_snapshots(
      tree_local,
      intervals_local,
      active_history_pc_local$full_history,
      alpha_val = if (input$evol_model == "OU") as.numeric(input$alpha_val) else NULL,
      model = input$evol_model
    )

    bm_history_pc_local <- simulate_and_reconstruct_history(
      tree_local,
      root_state_local,
      sigma_matrix_local,
      model = "BM",
      pc_labels = pca_space_local$sigma_info$pc_labels
    )
    ou_history_pc_local <- simulate_and_reconstruct_history(
      tree_local,
      root_state_local,
      sigma_matrix_local,
      target_state = target_local,
      model = "OU",
      alpha_val = as.numeric(input$alpha_val),
      pc_labels = pca_space_local$sigma_info$pc_labels
    )
    bm_snapshots_pc_local <- generate_all_snapshots(tree_local, intervals_local, bm_history_pc_local$full_history, model = "BM")
    ou_snapshots_pc_local <- generate_all_snapshots(tree_local, intervals_local, ou_history_pc_local$full_history, alpha_val = as.numeric(input$alpha_val), model = "OU")

    histories_pc <- list(
      active = active_history_pc_local$full_history,
      bm = bm_history_pc_local$full_history,
      ou = ou_history_pc_local$full_history
    )
    histories_local <- lapply(histories_pc, function(history_matrix) {
      reconstructed <- t(apply(history_matrix, 1, function(pc_scores) as.numeric(pc_scores %*% t(pca_space_local$rotation) + pca_space_local$center)))
      rownames(reconstructed) <- rownames(history_matrix)
      reconstructed
    })

    snapshots_pc <- list(
      active = active_snapshots_pc_local,
      bm = bm_snapshots_pc_local,
      ou = ou_snapshots_pc_local
    )
    snapshots_local <- lapply(snapshots_pc, function(snapshot_list) {
      lapply(snapshot_list, function(snapshot_matrix) {
        reconstructed <- t(apply(snapshot_matrix, 1, function(pc_scores) as.numeric(pc_scores %*% t(pca_space_local$rotation) + pca_space_local$center)))
        rownames(reconstructed) <- rownames(snapshot_matrix)
        reconstructed
      })
    })

    active_results_local <- histories_local$active
    bm_results_local <- histories_local$bm
    ou_results_local <- histories_local$ou
    active_snapshots_local <- snapshots_local$active
    bm_snapshots_local <- snapshots_local$bm
    ou_snapshots_local <- snapshots_local$ou

    active_snapshots_gpa_local <- perform_final_gpa_on_snapshots(active_snapshots_local, md$n_land)
    bm_snapshots_gpa_local <- perform_final_gpa_on_snapshots(bm_snapshots_local, md$n_land)
    ou_snapshots_gpa_local <- perform_final_gpa_on_snapshots(ou_snapshots_local, md$n_land)

    bm_tips_local <- bm_results_local[seq_len(Ntip(tree_local)), , drop = FALSE]
    bm_array_local <- reshape_to_array(bm_tips_local, md$n_land)
    gpa_bm_local <- gpagen(bm_array_local, print.progress = FALSE)

    ou_tips_local <- ou_results_local[seq_len(Ntip(tree_local)), , drop = FALSE]
    ou_array_local <- reshape_to_array(ou_tips_local, md$n_land)
    gpa_ou_local <- gpagen(ou_array_local, print.progress = FALSE)

    # Mantel diagnostics: compare Procrustes shape distances vs phylogenetic path distances.
    morph_dist_bm_local <- dist(two.d.array(gpa_bm_local$coords))
    morph_dist_ou_local <- dist(two.d.array(gpa_ou_local$coords))
    phylo_dist_local <- as.dist(cophenetic(tree_local))
    mantel_bm_local <- mantel(morph_dist_bm_local, phylo_dist_local, method = "pearson", permutations = 999)
    mantel_ou_local <- mantel(morph_dist_ou_local, phylo_dist_local, method = "pearson", permutations = 999)

    list(
      tree = tree_local,
      intervals = intervals_local,
      active_model = input$evol_model,
      active_results = active_results_local,
      active_results_pc = active_history_pc_local$full_history,
      active_reconstruction = active_history_pc_local$reconstruction,
      active_snapshots = active_snapshots_local,
      active_snapshots_aligned = active_snapshots_gpa_local$aligned_snapshots,
      bm_results = bm_results_local,
      ou_results = ou_results_local,
      bm_snapshots = bm_snapshots_local,
      ou_snapshots = ou_snapshots_local,
      bm_snapshots_aligned = bm_snapshots_gpa_local$aligned_snapshots,
      ou_snapshots_aligned = ou_snapshots_gpa_local$aligned_snapshots,
      pca_space = pca_space_local,
      node_ages = calibrated_tree_info$node_ages,
      applied_overrides = calibrated_tree_info$applied_overrides,
      target = pca_space_local$target_vector,
      n_land = md$n_land,
      phylo_dist = phylo_dist_local,
      morph_dist_bm = morph_dist_bm_local,
      morph_dist_ou = morph_dist_ou_local,
      mantel_bm = mantel_bm_local,
      mantel_ou = mantel_ou_local,
      wireframe_links = matrix(md$all_skulls$wireframe, ncol = 2)
    )
  }, ignoreInit = FALSE)

  # Diagnostics callbacks (all depend on sim_state()).
  # These render model-fit summaries and shape-time relationships.
  output$mantel_info <- renderPrint({
    sim <- sim_state()
    cat("BM Mantel Result:\n")
    print(sim$mantel_bm)
    cat("\nOU Mantel Result:\n")
    print(sim$mantel_ou)
  })

  output$pca_info <- renderText({
    sim <- sim_state()
    pca_model <- sim$pca_space$pca_model
    retained_variance <- pca_model$cumulative_variance[pca_model$retained_pcs] * 100
    threshold_variance <- pca_model$cumulative_variance[pca_model$threshold_pcs] * 100

    paste(
      "PCA Reduction:",
      paste0("Mode: ", if (input$pc_mode == "manual") "Manual count" else "Variance threshold"),
      paste0("Retained PCs: ", pca_model$retained_pcs, " of ", length(pca_model$explained_variance)),
      paste0("Variance captured by retained PCs: ", round(retained_variance, 2), "%"),
      paste0("PCs needed for ", round(as.numeric(input$pc_variance) * 100, 1), "% variance: ", pca_model$threshold_pcs, " (", round(threshold_variance, 2), "%)"),
      paste0("Max PCs allowed by taxa constraint: ", pca_model$max_pcs),
      sep = "\n"
    )
  })

  output$module_info <- renderText({
    sim <- sim_state()
    module_lines <- vapply(names(sim$pca_space$module_definitions), function(module_name) {
      landmark_ids <- paste(sim$pca_space$module_definitions[[module_name]], collapse = ",")
      paste0(module_name, ": ", landmark_ids)
    }, character(1))

    pc_lines <- paste0(sim$pca_space$sigma_info$pc_labels, " -> ", sim$pca_space$sigma_info$ordered_module_names)

    paste(
      "Modules:",
      paste(module_lines, collapse = "\n"),
      "",
      "Retained PC assignments:",
      paste(pc_lines, collapse = "\n"),
      sep = "\n"
    )
  })

  output$shape_time_plot_bm <- renderPlot({
    sim <- sim_state()
    plot(sim$phylo_dist, sim$morph_dist_bm, pch = 16, col = rgb(0, 0, 1, 0.5),
         main = "BM: Shape vs Time", xlab = "Phylo Distance", ylab = "Procrustes Distance")
    abline(lm(sim$morph_dist_bm ~ sim$phylo_dist), col = "red")
  })

  output$shape_time_plot_ou <- renderPlot({
    sim <- sim_state()
    plot(sim$phylo_dist, sim$morph_dist_ou, pch = 16, col = rgb(1, 0, 0, 0.5),
         main = "OU: Shape vs Time", xlab = "Phylo Distance", ylab = "Procrustes Distance")
    abline(lm(sim$morph_dist_ou ~ sim$phylo_dist), col = "blue")
  })

  output$morphospace_plot_bm <- renderPlot({
    sim <- sim_state()
    analyze_rigorous_morphospace(sim$bm_snapshots, sim$n_land, title = "BM: Shapes under drift", start_time_ma = min(sim$intervals), end_time_ma = max(sim$intervals))
  })

  output$morphospace_plot_ou <- renderPlot({
    sim <- sim_state()
    analyze_rigorous_morphospace(sim$ou_snapshots, sim$n_land, title = "OU: Shapes under selection", target_config = sim$target, start_time_ma = min(sim$intervals), end_time_ma = max(sim$intervals))
  })

  # Purpose: 
  # Access cached coordinates from the most recent phylo base plot.
  
  last_plot <- function() get("last_plot.phylo", envir = .PlotPhyloEnv)

  # Tree callback: renders the clickable topology view used for node/slice selection.
  output$tree_plot <- renderPlot({
    sim <- sim_state()
    plot(sim$tree, show.tip.label = TRUE, edge.width = 2)
    slice_times <- sim$intervals
    abline(v = slice_times, col = rgb(1, 0, 0, 0.2), lty = 2)
    internal_node_ids <- (Ntip(sim$tree) + 1):(Ntip(sim$tree) + sim$tree$Nnode)
    nodelabels(text = internal_node_ids, frame = "circle", bg = "steelblue", cex = 0.8)
    tiplabels(pch = 21, bg = "black", cex = 1)
  })

  output$node_age_reference <- renderText({
    sim <- sim_state()
    header <- "Current internal node ages (after calibration):"

    if (length(sim$applied_overrides) == 0) {
      paste(header, format_node_age_reference(sim$tree, sim$node_ages), sep = "\n")
    } else {
      paste(
        header,
        "Applied overrides:",
        paste0(names(sim$applied_overrides), " = ", round(sim$applied_overrides, 3), collapse = ", "),
        "",
        format_node_age_reference(sim$tree, sim$node_ages),
        sep = "\n"
      )
    }
  })

  # Selection callback:
  # Maps a click on the tree panel to either (a) nearest node or (b) active time-slice
  # branch intersection. This acts as the interaction bridge from tree -> 3D skull preview.
  selection_data <- reactive({
    sim <- sim_state()
    req(input$tree_click)
    lpt <- last_plot()

    node_dists <- sqrt((lpt$xx - input$tree_click$x)^2 + (lpt$yy - input$tree_click$y)^2)
    min_node_idx <- which.min(node_dists)
    if (node_dists[min_node_idx] < 0.6) return(list(source = "node", id = min_node_idx))

    slice_times <- sim$intervals
    slice_idx <- which.min(abs(slice_times - input$tree_click$x))
    edge_mat <- sim$tree$edge
    branch_y <- lpt$yy[edge_mat[,2]]
    valid_edges <- which(lpt$xx[edge_mat[,1]] <= slice_times[slice_idx] & lpt$xx[edge_mat[,2]] >= slice_times[slice_idx])

    if (length(valid_edges) > 0) {
      best_edge <- valid_edges[which.min(abs(branch_y[valid_edges] - input$tree_click$y))]
      if (abs(branch_y[best_edge] - input$tree_click$y) < 0.5) {
        return(list(source = "slice", time_idx = slice_idx, lineage_id = edge_mat[best_edge, 2]))
      }
    }

    return(NULL)
  })

  # Selection summary callback: human-readable description of current tree selection.
  output$selection_info <- renderPrint({
    sim <- sim_state()
    sel <- selection_data()
    if (is.null(sel)) return("No valid position selected.")
    if (sel$source == "node") {
      type <- if (sel$id <= Ntip(sim$tree)) "Tip" else "Internal Node"
      if (sel$id > Ntip(sim$tree)) {
        cat(type, "ID:", sel$id, "\nAge:", round(sim$node_ages[as.character(sel$id)], 3), "\nModel:", sim$active_model)
      } else {
        cat(type, "ID:", sel$id, "\nModel:", sim$active_model)
      }
    } else {
      cat("Intersection Found!\nLineage:", sel$lineage_id, "\nTime Slice:", sel$time_idx, "\nLocal target node:", sel$lineage_id, "\nModel:", sim$active_model)
    }
  })

  # 3D callback:
  # Uses selection_data() + sim_state() to fetch either node/tip coordinates
  # or a lineage snapshot at a selected time slice, then draws a wireframe skull.
  output$skull_3d <- renderPlotly({
    sim <- sim_state()
    sel <- selection_data()
    md <- morph_data_reactive()
    req(sel, md)

    if (sel$source == "node") {
      row_data <- as.numeric(sim$active_results[sel$id, ])
      coords <- matrix(row_data, ncol = 3, byrow = TRUE)
      p_color <- if (sim$active_model == "OU") "firebrick" else "steelblue"
      p_title <- paste(sim$active_model, "Node/Tip Configuration:", sel$id)
    } else {
      snap_full <- sim$active_snapshots[[sel$time_idx]]
      lineage_row <- which(rownames(snap_full) == as.character(sel$lineage_id))
      lineage_data <- snap_full[lineage_row, , drop = FALSE]
      coords <- matrix(as.numeric(lineage_data), ncol = 3, byrow = TRUE)
      p_color <- if (sim$active_model == "OU") "forestgreen" else "dodgerblue3"
      p_title <- paste(sim$active_model, "Lineage", sel$lineage_id, "at Slice", sel$time_idx)
    }

    build_wireframe_plotly(coords, sim$wireframe_links, p_title, point_color = p_color)
  })
}

shinyApp(ui, server)

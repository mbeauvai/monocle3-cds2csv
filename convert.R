# -----------------------------------------------------------------------------
# Title:        Monocle3 CDS to CSV Converter
# Description:  Converts a Monocle3 CellDataSet (CDS) object into a DynWrap 
#               output and exports multiple CSV files for trajectory analysis.
#
# Repository:   https://github.com/mbeauvai/monocle3-cds2csv
# Author:       mbeauvai
#
# Required Libraries:
# - dplyr               (License: MIT)            - Data manipulation
# - igraph              (License: GPL-2 | GPL-3)  - Graph operations
# - SingleCellExperiment (License: Artistic-2.0)  - Single-cell experiment structure
# - dynwrap             (License: MIT)            - Trajectory inference wrapper
#
# -----------------------------------------------------------------------------


# Required libraries
library(dplyr)               # MIT License
library(igraph)              # GPL-2 | GPL-3
library(SingleCellExperiment) # Artistic-2.0
library(dynwrap)             # MIT License

# ----------------------------------------------------------------------------
# Function: convertMonocle3CDSToDynWrap
# Description:
# - Converts a Monocle3 cell_data_set (CDS) object into a DynWrap output object.
# - Saves the milestone percentages, milestone network, progressions, 
#   dimensional reductions, and metadata to CSV files in the provided output directory.
# ----------------------------------------------------------------------------
convertMonocle3CDSToDynWrap <- function(cds, reduction_method, output_dir) {
  
  # Step 1: Extract the graph for the reduction method's trajectory
  g <- cds@principal_graph@listData[[reduction_method]]
  
  # Define dimensions (1 and 2) for reduced space visualization
  dims <- c(1, 2)
  x <- dims[[1]]
  y <- dims[[2]]
  
  # Step 2: Extract coordinates of the cells in the reduction_method space
  S_matrix <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[, c(x, y)])  # Use first two dimensions
  
  # Rename columns and add sample names
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- rownames(data_df)
  
  # Combine the extracted coordinates with the column data (metadata)
  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  
  # Step 3: Extract coordinates of vertices in the trajectory
  ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
    as.data.frame() %>%
    dplyr::select(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
  
  # Extract the minimum spanning tree (MST) from the graph
  dp_mst <- cds@principal_graph[[reduction_method]]
  
  # Step 4: Extract the edges between the vertices in the trajectory
  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    dplyr::select(source = "from", target = "to") %>%
    dplyr::left_join(ica_space_df %>%
                       dplyr::select(
                         source = "sample_name",
                         source_prin_graph_dim_1 = "prin_graph_dim_1",
                         source_prin_graph_dim_2 = "prin_graph_dim_2"),
                     by = "source") %>%
    dplyr::left_join(ica_space_df %>%
                       dplyr::select(
                         target = "sample_name",
                         target_prin_graph_dim_1 = "prin_graph_dim_1",
                         target_prin_graph_dim_2 = "prin_graph_dim_2"),
                     by = "target")
  
  # Step 5: Create the trajectory graph using the edge data
  g <- graph_from_data_frame(edge_df, directed = FALSE)
  
  # Step 6: Extract unique vertices and their coordinates from the edge data
  unique_vertices <- unique(c(edge_df$source, edge_df$target))
  x_coords <- setNames(c(edge_df$source_prin_graph_dim_1, edge_df$target_prin_graph_dim_1),
                       c(edge_df$source, edge_df$target))
  y_coords <- setNames(c(edge_df$source_prin_graph_dim_2, edge_df$target_prin_graph_dim_2),
                       c(edge_df$source, edge_df$target))
  
  # Assign coordinates to the graph's vertices
  V(g)$x <- x_coords[V(g)$name]
  V(g)$y <- y_coords[V(g)$name]
  
  # Step 7: Create a vertex dataframe (source and target vertices) with their coordinates
  source_df <- edge_df[, c("source", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
  names(source_df) <- c("vertex", "dim_1", "dim_2")
  
  target_df <- edge_df[, c("target", "target_prin_graph_dim_1", "target_prin_graph_dim_2")]
  names(target_df) <- c("vertex", "dim_1", "dim_2")
  
  vertices_df <- rbind(source_df, target_df)  # Combine source and target vertices
  
  unique_vertices_df <- unique(vertices_df)   # Remove duplicates
  
  # Convert vertices data to a matrix for milestone calculations
  unique_vertices_matrix <- as.matrix(unique_vertices_df[, c("dim_1", "dim_2")])
  rownames(unique_vertices_matrix) <- unique_vertices_df$vertex
  
  # Step 8: Calculate degrees of the graph vertices
  degs <- igraph::degree(g)
  
  # Optionally simplify the graph if desired
  simplify <- FALSE
  if (simplify) {
    milestone_ids <- names(which(degs != 2))
    mst <- dynwrap::simplify_igraph_network(g)
  } else {
    milestone_ids <- names(degs)
    mst <- dynwrap::simplify_igraph_network(g, force_keep = milestone_ids)
  }
  
  # Convert the simplified graph to a milestone network data frame
  milestone_network <- mst %>%
    igraph::as_data_frame()
  
  # Step 9: Add length to the milestone network and remove weight
  milestone_network$length <- milestone_network$weight
  milestone_network$weight <- NULL
  
  # Step 10: Construct dimensional reduction (dimred) output
  dimred <- matrix(
    c(data_df$data_dim_1, data_df$data_dim_2),
    ncol = 2,
    byrow = FALSE,
    dimnames = list(data_df$sample_name, c('dim_1', 'dim_2'))
  )
  
  # Construct milestone dimred output
  dimred_milestones <- unique_vertices_matrix[milestone_ids, , drop = FALSE]
  
  # Rename milestone IDs
  milestone_id_map <- setNames(paste0("M", seq_along(milestone_ids)), milestone_ids)
  rownames(dimred_milestones) <- milestone_id_map[rownames(dimred_milestones)]
  
  # Update milestone network with new milestone IDs
  milestone_ids <- milestone_id_map[milestone_ids] %>% setNames(NULL)
  milestone_network <- milestone_network %>%
    mutate(from = milestone_id_map[from] %>% setNames(NULL),
           to = milestone_id_map[to] %>% setNames(NULL))
  
  # Remove unneeded columns from the milestone network
  milestone_network$source_prin_graph_dim_1 <- NULL
  milestone_network$source_prin_graph_dim_2 <- NULL
  milestone_network$target_prin_graph_dim_1 <- NULL
  milestone_network$target_prin_graph_dim_2 <- NULL
  
  # Step 11: Wrap the data and save it as a dynwrap object
  cell_ids <- rownames(dimred)
  
  # Create milestone network and dimred projection using dynwrap
  milestone_network2 <- milestone_network %>% distinct(from, to, .keep_all = TRUE)
  
  # Create the dynwrap output object
  output_dynwrap <- dynwrap::wrap_data(cell_ids = cell_ids) %>%
    dynwrap::add_dimred_projection(
      milestone_network = milestone_network2,
      dimred = dimred,
      dimred_milestones = dimred_milestones
    )
  
  # Step 12: Save all output data to CSV files
  write.csv(output_dynwrap$milestone_percentages, file.path(output_dir, "milestone_percentages.csv"), row.names = FALSE)
  write.csv(output_dynwrap$milestone_network, file.path(output_dir, "trajectory_edges.csv"), row.names = FALSE)
  write.csv(output_dynwrap$progressions, file.path(output_dir, "progressions.csv"), row.names = FALSE)
  
  # Save the dimensional reduction and milestones as CSV
  dimred_milestones_df <- as.data.frame(output_dynwrap$dimred_milestones)
  dimred_milestones_df$milestone_id <- rownames(dimred_milestones_df)
  write.csv(dimred_milestones_df, file.path(output_dir, "dimred_milestone.csv"), row.names = FALSE)
  
  dimred_df <- as.data.frame(output_dynwrap$dimred)
  dimred_df$cell_id <- rownames(dimred_df)
  write.csv(dimred_df, file.path(output_dir, "dimred.csv"), row.names = FALSE)
  
  # Save cell metadata as CSV
  cell_metadata_df <- as.data.frame(colData(cds))
  cell_metadata_df$cell_id <- rownames(cell_metadata_df)
  write.csv(cell_metadata_df, file.path(output_dir, "cell_metadata.csv"), row.names = FALSE)
}

# Example usage for various datasets
# Pancreas Dataset
output_dir_pancreas <- "/home/user1/Pancreas/"
cds_pancreas <- readRDS(paste0(output_dir_pancreas, "study1-trajectory-method_monocle3-default_cds.rds"))
convertMonocle3CDSToDynWrap(cds_pancreas, reduction_method = "UMAP", output_dir = output_dir_pancreas)

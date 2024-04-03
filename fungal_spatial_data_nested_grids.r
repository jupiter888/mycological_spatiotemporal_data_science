library(sf)
library(dplyr)
library(ggplot2)
library(parallel)
library(pryr)
library(pbapply)

bbox_to_sf <- function(bbox, crs) {
  # create polygon from bbox coordinates
  poly <- st_polygon(list(rbind(
    c(bbox["xmin"], bbox["ymin"]),
    c(bbox["xmax"], bbox["ymin"]),
    c(bbox["xmax"], bbox["ymax"]),
    c(bbox["xmin"], bbox["ymax"]),
    c(bbox["xmin"], bbox["ymin"])
  )))
  
  # create sf object from polygon w/ specified CRS
  sf_object <- st_sf(geometry = st_sfc(poly, crs = crs))
  return(sf_object)
}

create_grids_based_on_data <- function(data, cell_sizes) {
  if (is.na(st_crs(data))) {
    stop("Data does not have a CRS defined.")
  }
  
  bbox <- st_bbox(data)
  bbox_sf <- bbox_to_sf(bbox, st_crs(data))
  
  grids <- lapply(cell_sizes, function(size) {
    cell_size_in_degrees <- size / 110.574
    grid <- st_make_grid(bbox_sf, cellsize = c(cell_size_in_degrees, cell_size_in_degrees), what = "polygons") %>%
      st_sf() %>%
      mutate(id = 1:nrow(.))
    return(grid)
  })
  print("created grids")
  return(grids)
}
perform_spatial_join <- function(grid, data) {
  joined_data <- st_join(data, grid, join = st_within)
  print("Joined data....")
  return(joined_data)
}

visualize_data <- function(data, grid_list) {
  plot <- ggplot() + theme_minimal()
  
  grid_colors <- c("red", "green", "blue", "orange")
  for (i in seq_along(grid_list)) {
    plot <- plot + geom_sf(data = grid_list[[i]], fill = NA, color = grid_colors[i], size = 0.5)
  }
  plot1 <- plot + geom_sf(data = data, aes(color = species), size = 2) +
    labs(title = "Fungi Distribution Across Grids")
  show(plot1())
  return(plot1)
}

visualize_2 <- function(data, grids) {
  plot <- ggplot(data) +
    geom_sf() +
    geom_point(aes(x = decimalLongitude, y = decimalLatitude, color = species))
  
  for (grid in grids) {
    plot <- plot + geom_sf(data = grid, fill = NA, color = "black")
  }
  
  plot <- plot + labs(title = "Species Distribution with Grids", x = "Longitude", y = "Latitude") + theme_bw()
  return(plot)
}

(bbox_to_sf, create_grids_based_on_data, perform_spatial_join, visualize_data) remain unchanged

process_and_visualize_fungi_data <- function(fungi_data_path, output_path, cell_sizes, original_crs) {
  # Start profiling
  Rprof("my_profile_output.out")
  
  start_time <- Sys.time()
  
  
  print("reading fungi data...")
  fungi_data <- read.csv(fungi_data_path, header = TRUE, sep = "\t")
  
  if (!all(sapply(fungi_data[,c("decimalLongitude", "decimalLatitude")], is.numeric))) {
    stop("Long & Lat cols must must be numeric.")
  }
  
  if(any(is.na(fungi_data$decimalLongitude)) || any(is.na(fungi_data$decimalLatitude))) {
    stop("Long & Lat cols must !contain NA values.")
  }
  
  fungi_sf <- st_as_sf(fungi_data, coords = c("decimalLongitude", "decimalLatitude"), crs = original_crs, agr = "constant")
  
  print("creating grids...")     
  grids <- create_grids_based_on_data(fungi_sf, cell_sizes)
  print("Printing class of grids var: \n")
  print(class(grids))
  # To check the class of each element in the 'grids' list
  print("Printing class of each element in the grids list: \n")
  print(lapply(grids, class))
  print("Printing grids info completed. \n Now beginning spatial joins in parallel")
  
  # performing spatial joins (parallel processing)
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, varlist = c("perform_spatial_join", "fungi_sf"))
  joined_results <- parLapply(cl, grids, function(grid) perform_spatial_join(grid, fungi_sf))
  stopCluster(cl)
  
  print("visualizing data (1/2)...")
  plot <- visualize_data(fungi_sf, grids)
  print(plot)
  show(plot)
  
  print("visualizing data (2/2)...")
  plot2 <- visualize_2(fungi_sf, grids)
  show(plot2)                                 
  print(plot2)
  
  # parallel processing results
  cl <- makeCluster(num_cores)
  clusterExport(cl, varlist = c("_joined_results_", "output_path", "cell_sizes"))
  parLapply(cl, seq_along(joined_results), function(i) {
    write.csv(as.data.frame(joined_results[[i]]), paste0(output_path, "joined_data_grid_", cell_sizes[i], ".csv"), row.names = FALSE)
  })
  stopCluster(cl)
  
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  print(paste("SUCCESSFUL COMPLETION! Total process duration using parallel processing:", round(duration, 2), "minutes"))
  # stop profiling
  Rprof(NULL)
  
  # analyzing profiling data
  print(summaryRprof("my_profile_output.out"))
}

# set cell sizes parameters here
cell_sizes <- c(10, 25, 100 )
original_crs <- "EPSG:4326"
result <- process_and_visualize_fungi_data("/Desktop/myco_1/Fungi_1970_cleaned.csv", "/Desktop/myco_1/results/10_25_100/", cell_sizes, original_crs)
# load results back into the global environment
saveRDS(results, "results_3.rds")

results_read_in <- readRDS("results_3.rds")

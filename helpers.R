# Convert degrees to radians
convertDegrees <- function(deg) {
  rad <- deg * pi / 180
  return(rad)
}

# Haversine 
haverFunction <- function(lat1, lon1, lat2, lon2){
  '
  Calculate the straightline distance between a pair of points.

  Inputs: Point1 lat, Point1 lon, Point2 lat, Point2 lon
  Output: Straightline distance in km between two points.

  #haversine = sin^2(theta/2)
  #d = [2r][arcsin(sqrt(hav(lat2 - lat1)+cos(lat1)*cos(lat2)*hav(long2 - long1))]
  #d = [2r][arcsin(sqrt(sin^2((lat2 - lat1)/2)+cos(lat1)*cos(lat2)*sin^2((long2 - long1)/2)))]
  '
  
  earthRadius <- 6378.1 #in km
  lat1Rad <- convertDegrees(lat1)
  lat2Rad <- convertDegrees(lat2)
  lon1Rad <- convertDegrees(lon1)
  lon2Rad <- convertDegrees(lon2)
  dLat <- lat2Rad - lat1Rad
  dLon <- lon2Rad - lon1Rad
  
  bracketCalc <- sin(dLat/2)^2 + cos(lat1Rad) * cos(lat2Rad) * (sin(dLon/2))^2
  d <- 2 * earthRadius * asin(bracketCalc^0.5)
  return(d)
}

# Two dimensional balancing matrix
balance <- function(matrix, tot, axis) {
  if (axis == 1) {
    sum <- rowSums(matrix)
  } else if (axis == 2) {
    sum <- colSums(matrix)
  }
  sc <- tot / sum 
  sc[is.nan(sc)] <- 0
  
  # MARGIN = 1 indicates rows, MARGIN = 2 indicates columns
  matrix2 <- sweep(matrix, MARGIN = axis, sc, `*`)
  return(matrix2)
}

calc_error <- function(matrix, a, b) {
  row_sum <- sum(abs(a - rowSums(matrix)))
  col_sum <- sum(abs(b - colSums(matrix)))
  return(row_sum + col_sum)
}

matrix_balancing <- function(matrix, a, b, totals_to_use = "raise", max_iterations = 10000, rel_error = 0.0001) {
  valid_totals_to_use = c("rows", "columns", "average", "raise")
  if (!(totals_to_use %in% valid_totals_to_use)) {
    stop("totals_to_use is invalid, not one of ('rows', 'columns', 'average', 'raise')")
  }
  
  # Scale the column and row totals, if specified
  a_sum = sum(a)
  b_sum = sum(b)
  print(paste0("Sum of a is: ", sum(a), ". Sum of b is: ", sum(b), "."))
  
  if (!(a_sum == b_sum)) {
    if (totals_to_use == "rows") {
      b = b * (a_sum / b_sum)
      print("Scaled b to the row totals")
    } else if (totals_to_use == "columns") {
      a = a * (b_sum / a_sum)
      print("Scaled a to the column totals")
    } else if (totals_to_use == "average") {
      avg_sum = 0.5 * (a_sum + b_sum)
      a = a * (avg_sum / a_sum)
      b = b * (avg_sum / b_sum)
      print("Scaled a and b to the average totals")
    } else {
      stop("a and b vector totals do not match")
    }
  } 
  
  print(paste0("Sum of scaled a is: ", sum(a), ". Sum of scaled b is: ", sum(b), "."))
  
  error <- 1.0
  i <- 0
  init_error <- calc_error(matrix, a, b)
  matrix2 <- matrix
  
  while (error > rel_error) {
    if (i > max_iterations) {
      print("Matrix balancing did not converge within iteration limit")
      break
    }
    matrix2 <- balance(matrix2, a, 1)
    matrix2 <- balance(matrix2, b, 2)
    error <- calc_error(matrix2, a, b) / init_error
    i <- i + 1
  }
  return(matrix2)
}

# Bucket rounding
# [https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum]
smart_round <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  # Get the biggest residuals that will add up the difference
  # from the floor round
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

# Apply sample() row by row
sample_by_row <- function(row) {
  x <- row["sfis.list"][[1]]
  size <- row["enrolment.rounded"][[1]]
  
  print(paste0("Number of schools in selection is: ", length(x), ". Size of students: ", size))
  
  if (length(x) == 1) {
    # Quirk 
    prob = c(rep(0, x - 1), row["school.weight.prob.list"][[1]])
  } else {
    prob = row["school.weight.prob.list"][[1]]
  }
  list(sample(x, size=size, replace=TRUE, prob=prob))
}

# Gravity model 
furness <- function(cmat, observed_matrix){
  
  '
  This function conducts furnessing or two dimensional balancing
  
  Inputs: cost function, observed matrix
  Output: balanced matrix
  '
  
  #' row sum of cost function matrix
  skim_trow <- as.data.frame(rowSums(cmat))
  #' row sum of observed matrix
  obs_tatt <- as.data.frame(rowSums(observed_matrix))
  
  #' row balancing factors
  rowbal <- obs_tatt/skim_trow
  colnames(rowbal) <- 'ratio'
  # rowbal1 <- do.call('cbind', replicate(ncol(observed_matrix), rowbal, 
  #                                       simplify = FALSE))
  
  #first iteration
  cfunc1 <- cmat*rowbal$ratio
  
  #' col sum of cost function matrix
  skim_tcol <- as.data.frame(colSums(cfunc1))
  #' col sum of observed matrix
  obs_tprod <- as.data.frame(colSums(observed_matrix))
  
  #' column balancing factors
  colbal <- obs_tprod/skim_tcol
  colnames(colbal) <- 'ratio'
  # colbal1 <- do.call('rbind', replicate(nrow(observed_matrix), colbal, 
  #                                       simplify = FALSE))
  
  # next iteration. Transpose the matrix to allow multiplying by the column ratios.
  # once done retranspose the matrix before running the row balancing
  cfunc1 <- t(cfunc1)*colbal$ratio
  cfunc1 <- t(cfunc1)
  
  return(cfunc1)
} 

# Plot TLFD with OD skim 
generate_tlfd <- function(observed_trips, simulated_trips, dist_matrix) {
  
  # Cuts the the length into 50 bins
  dist_matrix <- transform(dist_matrix, km_bin = cut(dist_matrix$value, 50))
  
  # Merge with dsitance matrix
  obs_tlfd <- merge(observed_trips, dist_matrix, by = c('orig', 'dest')) %>%
    group_by(km_bin) %>%
    summarise(trips = sum(trips)) %>%
    transform(flag = 'obs')
  
  # Prepare the simulated trips to merge with distance
  sim_tlfd <- as.data.frame(simulated_trips) %>%
    transform(orig = rownames()) %>%
    # collapse the dataframe
    melt(id.vars = c('orig')) %>%
    transform(variable = substring(.$variable, 2)) %>%
    filter(value > 0)
  colnames(sim_tlfd) <- c('orig', 'dest', 'trips')
  
  # Merge with distance matrix
  sim_tlfd <- merge(sim_tlfd, dist_matrix, by = c('orig', 'dest')) %>%
    group_by(km_bin) %>%
    summarise(trips = sum(trips)) %>%
    transform(flag = 'model')
  
  # Combine the observed and simulated tlfd into one datframe
  combined_tlfd <- rbind(obs_tlfd, sim_tlfd)
  
  return(combined_tlfd)
}

# Plot TLFD with segmentation
plot_straight_line_tlfd <- function(tlfd) {
  g1 <- ggplot(tlfd, aes(dist, color = as.factor(board_type_name))) +
    geom_freqpoly(stat = 'bin', binwidth = 0.1, size = 1) + 
    facet_wrap(vars(board_type_name), scales = 'free_y') +
    labs(title = 'Euclidean Distance', x = 'distance (km)', color = 'Board Type') +
    theme_grey() +
    coord_cartesian(xlim = c(-1, 10))
  
  g2 <- ggplot(tlfd, aes(dist, color = as.factor(mof_region))) +
    geom_freqpoly(stat = 'bin', binwidth = 0.1, size = 1) + 
    facet_wrap(vars(board_type_name), scales = 'free_y') +
    labs(title = 'Euclidean Distance', x = 'distance (km)', color = 'Geographic Location') +
    theme_grey() +
    coord_cartesian(xlim = c(-1, 10))
  
  g3 <- tlfd %>%
    filter(startsWith(board_type_name, 'English')) %>%
    ggplot(aes(dist, color = as.factor(mof_region))) +
    # geom_area(stat = 'bin', binwidth = 0.1) +
    geom_freqpoly(stat = 'bin', binwidth = 0.1, size = 1) + 
    facet_grid(cols = vars(board_type_name), rows = vars(as.factor(panel)), scales = 'free_y') +
    labs(title = 'Euclidean Distance', x = 'distance (km)', color = 'Geographic Location') +
    theme_grey() +
    coord_cartesian(xlim = c(-1, 10))
  
  g4 <- tlfd %>%
    filter(startsWith(board_type_name, 'French')) %>%
    ggplot(aes(dist, color = as.factor(mof_region))) +
    # geom_area(stat = 'bin', binwidth = 0.1) +
    geom_freqpoly(stat = 'bin', binwidth = 0.1, size = 1) + 
    facet_grid(cols = vars(board_type_name), rows = vars(as.factor(panel)), scales = 'free_y') +
    labs(title = 'Euclidean Distance', x = 'distance (km)', color = 'Geographic Location') +
    theme_grey() +
    coord_cartesian(xlim = c(-1, 10))
  
  grid.arrange(g1, g2, g3, g4, nrow = 4)
}

# Plot TLFD with segmentation
plot_travel_time_tlfd <- function(tlfd) {
  g1 <- ggplot(tlfd, aes(travel.time, color = as.factor(board_type_name))) +
    geom_freqpoly(stat = 'bin', binwidth = 2, size = 1) + 
    facet_wrap(vars(board_type_name), scales = 'free_y') +
    labs(title = 'Travel Time', x = 'time (min)', color = 'Board Type') +
    theme_grey() +
    coord_cartesian(xlim = c(-1, 30))
  
  g2 <- ggplot(tlfd, aes(travel.time, color = as.factor(mof_region))) +
    geom_freqpoly(stat = 'bin', binwidth = 2, size = 1) +
    facet_wrap(vars(board_type_name), scales = 'free_y') +
    labs(title = 'Travel Time', x = 'time (min)', color = 'Geographic Location') +
    theme_grey() +
    coord_cartesian(xlim = c(-1, 30))
  
  g3 <- tlfd %>%
    filter(startsWith(board_type_name, 'English')) %>%
    ggplot(aes(travel.time, color = as.factor(mof_region))) +
    geom_freqpoly(stat = 'bin', binwidth = 2, size = 1) +
    facet_grid(cols = vars(board_type_name), rows = vars(as.factor(panel)), scales = 'free_y') +
    labs(title = 'Travel Time', x = 'time (min)', color = 'Geographic Location') +
    theme_grey() +
    coord_cartesian(xlim = c(-1, 30))
  
  g4 <- tlfd %>%
    filter(startsWith(board_type_name, 'French')) %>%
    ggplot(aes(travel.time, color = as.factor(mof_region))) +
    geom_freqpoly(stat = 'bin', binwidth = 2, size = 1) +
    facet_grid(cols = vars(board_type_name), rows = vars(as.factor(panel)), scales = 'free_y') +
    labs(title = 'Travel Time', x = 'time (min)', color = 'Geographic Location') +
    theme_grey() +
    coord_cartesian(xlim = c(-1, 30))
  
  g <- arrangeGrob(g1, g2, g3, g4, nrow = 4)
}

# Obtain school and student xy from `student_travel_##` where ## is the catchment distance
create_student_xy <- function(student_travel) {
  '
  This function filters the individual student information and transform the projection to the same
  as the TRESO shapefiles

  inputs: Dataframe of the students within a certain catchment distance of the school
  output: SpatialPointsDataFrame of each student
  '
  # Convert the student dataframe into SpatialPointsDataframe
  student_spdf <- student_travel %>%
    ungroup() %>%
    select(student.lat, student.long, dist, school.name, sfis, dsb.index, panel, enrolment, student.postal.code) %>%
    rename(
      lat = student.lat,
      long = student.long,
      euclidean.dist = dist
    ) %>%
    mutate(id = row_number())
  
  coordinates(student_spdf) <- c('long', 'lat')
  
  # Project the student and school points from lat/long to TRESO's LCC specification
  proj4string(student_spdf) <- CRS('+proj=longlat +datum=WGS84')
  treso_projarg = treso_shp@proj4string@projargs
  
  student_xy <- spTransform(student_spdf, CRS(treso_projarg))
  
  return(student_xy)
}

create_school_xy <- function(student_travel) {
  '
  This function filters the individual school information and transform the projection to the same
  as the TRESO shapefiles

  inputs: Dataframe of the students within a certain catchment distance of the school
  output: SpatialPointsDataFrame of each school
  '
  # Convert the school dataframe into SpatialPointsDataframe
  school_spdf <- student_travel %>%
    ungroup() %>%
    select(school.name, school.lat, school.long, dsb.index, bsid, sfis, panel, perc.dist) %>%
    group_by(sfis) %>%
    summarise_all(funs(first)) %>%
    rename(
      lat = school.lat,
      long = school.long,
      catchment.dist = perc.dist
    ) %>%
    mutate(id = row_number())
  coordinates(school_spdf) <- c('long', 'lat')
  
  # Project the student and school points from lat/long to TRESO's LCC specification
  proj4string(school_spdf) <- CRS('+proj=longlat +datum=WGS84')
  treso_projarg = treso_shp@proj4string@projargs
  
  school_xy <- spTransform(school_spdf, CRS(treso_projarg))
  
  return(school_xy)
}

create_school_xy_from_school <- function(school_sfis) {
  '
  This function differs fro `create_school_xy` in that this takes in the school dataframe
  without the catchment distnace
  
  input: Dataframe of school with lat long
  output: SpatialPointsDataFrame of each school
  '
  # Convert the school dataframe into SpatialPointsDataframe
  school_spdf <- school_sfis %>%
    ungroup() %>%
    select(school.name, school.lat, school.long, dsb.index, bsid, sfis) %>%
    rename(
      lat = school.lat,
      long = school.long
    ) %>%
    mutate(id = row_number())
  coordinates(school_spdf) <- c('long', 'lat')
  
  # Project the student and school points from lat/long to TRESO's LCC specification
  proj4string(school_spdf) <- CRS('+proj=longlat +datum=WGS84')
  treso_projarg = treso_shp@proj4string@projargs
  
  school_xy <- spTransform(school_spdf, CRS(treso_projarg))
  
  return(school_xy)
}

# Spatial gymnastics to find the matching TRESO zone for student or school XY locations
create_overlay <- function(xy_location, treso_shp, type = 'student') {
  '
  This function takes the SpatialPointsDataFrame of either school or students and maps the
  XY location of the object to the TRESO zone and returns the appropriate TRESO zone ID.

  inputs: SpatialPointsDataFrame of school/student, TRESO shapefile, string indicating what is the object
  output: Dataframe of the school or student with the appropriate TRESO zone ID
  '
  if (type == 'student'){
    # Find the treso zones which the student points layover
    overlay <- over(xy_location, treso_shp, returnList = FALSE) %>%
      cbind(euclidean.dist = xy_location@data$euclidean.dist,
            student.postal.code = xy_location@data$student.postal.code,
            school.name = xy_location@data$school.name,
            sfis = xy_location@data$sfis,
            enrolment = xy_location@data$enrolment) %>%
      as_tibble() %>%
      select(Treso_ID, euclidean.dist, student.postal.code, enrolment, school.name, sfis) %>%
      rename(
        treso.id.por = Treso_ID
      )
  }
  if (type == 'school'){
    # Find the treso zones which the school points layover
    overlay <- over(xy_location, treso_shp, returnList = FALSE) %>%
      cbind(catchment.dist = xy_location@data$catchment.dist,
            dsb.index= xy_location@data$dsb.index,
            school.name = xy_location@data$school.name,
            sfis = xy_location@data$sfis,
            bsid = xy_location@data$bsid) %>%
      as_tibble() %>%
      select(Treso_ID, catchment.dist, dsb.index, school.name, sfis, bsid) %>%
      rename(
        treso.id.pos = Treso_ID
      )
  }
  return(overlay)
}

# Buffer the TRESO zones for each school based on the catchment distance
buffer_zones <- function(school_xy, treso_shp) {
  "
  Performs spatial buffers to select the TRESO zones within the catchment distance of a school.

  Inputs: School's XY location, TRESO shapefile
  Output: Dataframe of TRESO zones within each School
  "
  buffered_points <- gBuffer(school_xy, width = school_xy@data$catchment.dist * 1000, byid = TRUE)
  
  tic('Multi-Core Buffering')
  # Start multi-core clustering for performance
  cores = parallel::detectCores()
  cl <- parallel::makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  
  # foreach() loop to find the TRESO zones touching the school buffers
  datalist <- list()
  datalist <- foreach(i = 1:nrow(buffered_points@data), .packages=c('sp', 'tidyverse', 'rgeos')) %dopar% {
    
    # Select the buffer
    buffer <- buffered_points[buffered_points@data$id == i, ]
    
    # Get all TRESO zones that touches the buffer with over(), gIntersect() could be another solution
    treso_overlay <- over(treso_shp, buffer, returnList = FALSE) %>%
      cbind(treso.id = treso_shp@data$Treso_ID,
            households = treso_shp@data$TotHH,
            persons = treso_shp@data$TotPers,
            shape.area = treso_shp@data$Area) %>%
      as_tibble() %>%
      select(school.name, sfis, dsb.index, catchment.dist, treso.id, households, persons, shape.area) %>%
      drop_na()
    
    # Save each tibble in a list of tibbles
    datalist[[i]] <- treso_overlay
  }
  
  stopCluster(cl)
  toc()
  
  # Write the datalist to csv
  datalist_df <- rbindlist(datalist)
  return(datalist_df)
}

# Condense the socio-economic information of the buffered TRESO zones for each school
summarize_buffered_zones <- function(buffered_df, treso_tb, school_ade, school_board_def, treso_zone_def) {
  '
  This function takes the Dataframe of TRESO zones associated with each school and summarize the socio-economic
  data for each school.

  inputs: SpatialPointsDataFrame of school/student, TRESO shapefile, string indicating what is the object
  output: Dataframe of the school or student with the appropriate TRESO zone ID
  '
  # First summarise data that is calculated with `mean()` or `first()`
  # I think using purr, one could summarise different columns with different functions in one go
  
  school_tb_temp <- as_tibble(buffered_df) %>%
    left_join(select(treso_tb, treso_zone, mean_income, mean_age), by = c('treso.id' = 'treso_zone')) %>%
    left_join(select(school_board_def, dsb, board_type_name), by = c('dsb.index' = 'dsb')) %>%
    left_join(select(treso_zone_def, treso_id, area, mof_region), by = c('treso.id' = 'treso_id')) %>%
    replace(is.na(.), 0) %>%
    group_by(sfis) %>%
    summarise(
      mean.income = mean(mean_income),
      mean.age = mean(mean_age),
      dsb.index = first(dsb.index),
      school.name = first(school.name),
      catchment.dist = first(catchment.dist),
      board.type.name = first(board_type_name),
      area = names(which.max(table(area))),
      mof.region = names(which.max(table(mof_region)))
    )
  
  # Then, combine with data that is calculated with sum
  school_tb <- as_tibble(buffered_df) %>%
    left_join(treso_tb, by = c('treso.id' = 'treso_zone')) %>%
    replace(is.na(.), 0) %>%
    group_by(sfis) %>%
    summarise_at(
      .vars = vars(starts_with('n_'), starts_with('occu_'), starts_with('deg_'), 'attend_school', 'shape.area'),
      .funs = c(sum = 'sum')
    ) %>%
    left_join(school_tb_temp, by = 'sfis')
  
  # Combine School's ADE and utilization info with School's buffered zones
  school_tb <- school_tb %>%
    left_join(school_ade, by = c('sfis'))
  
  # Clean up the global environments
  rm(school_tb_temp)
  
  return(school_tb)
}

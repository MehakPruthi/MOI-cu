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

matrix_balancing_1d <- function(matrix, a, weight, axis=1, constrained=TRUE) {
  
  #' One dimensional balances a matrix
  #' 
  #' One dimensional matrix balancing with optional scaling. If there are known scale weights that needs to be
  #' applied to the propensity matrix, it can be supplied in argument `weight`. The `axis` argument determine which
  #' direction to balance the matrix, 1 indicates rows (origin) and 2 indicates columns (destination). The 
  #' `constrained` flag determines if the weighting vector is used as a hard limit or as a boolean value to 
  #' prevent trips to be sent. 
  #' 
  #' @param matrix The propensity matrix
  #' @param a The totals to balance against
  #' @param weight The weight vector to scale against
  #' @param axis The direction to perform the balance
  #' @param constrianed A flag to determine if the weight is to be used as a constraint or not
  #' @return One dimentionally balanced matrix
  
  # Check if `axis` is 1 or 2
  if (!(axis %in% c(1, 2))) {
    stop("`axis` value is invalid, not one of (1, 2)")
  }
  
  # Check if weigth vector is supplied
  if (missing(weight)) {
    print("`weight` was not supplied, 1D balancing done without any scaling")
    tot = a
    matrix2 <- balance(matrix, tot, axis)
  } else {
    print("`weight` was supplied, 1D balancing done with scaling")
    tot = a
    # If `axis` is 1, balance the matrix against the rows but scale against the columns
    if (axis == 1) {
      sum = colSums(matrix)
      axis_to_scale = 2
    } 
    # If `axis` is 2, balance the matrix against the columns but scale against the rows
    else if (axis == 2) {
      sum = rowSums(matrix)
      axis_to_scale = 1
    }
    # Normalize the propensity matrix to ensure the scaled probailitiy does not exceed 1
    matrix_norm = sweep(matrix, MARGIN = axis_to_scale, sum, `/`)
    
    # If constrained, scale the matrix to the values of the weight vector
    # It not constrained, apply a boolean value to the matrix
    if (!(constrained)) {
      weight = ifelse(weight == 0.0, 0, 1)
    }
    # Scale the normalized propensity matrix to the weight vector
    matrix_scaled = sweep(matrix_norm, MARGIN = axis_to_scale, weight, `*`)
    
    # 1D balance against the tot vector on the specified axis
    matrix2 <- balance(matrix_scaled, tot, axis)
  }
  
  return(matrix2)
}

matrix_balancing_2d <- function(matrix, a, b, totals_to_use = "raise", max_iterations = 10000, rel_error = 0.0001) {
  '
  Two dimensional matrix balancing.
  
  Inputs: Propensity matrix, Origin totals to balance against, Destination totals to balance against, 
          The method of matching the totals, Maximum iterations to be used, Relative error to be required
  Ouput: Balanced matrix
  '
  valid_totals_to_use = c("rows", "columns", "average", "raise")
  if (!(totals_to_use %in% valid_totals_to_use)) {
    stop("totals_to_use is invalid, not one of ('rows', 'columns', 'average', 'raise')")
  }
  
  # Match the column and row totals, if not matching already and specified
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

# Calculate the simulated trips based on alpha and beta values
calculate_simulated_trips <- function(observed_trips, cost, alpha, beta) {
  cfunc <- cost %>%
    mutate(value = value^alpha * exp(beta*value)) %>%
    replace_na(value = 0)
  
  cfunc_rowsums <- cfunc %>%
    group_by(treso.id.por) %>%
    summarise(rowsum = sum(value))
  
  t <- left_join(cfunc, cfunc_rowsums, by = "treso.id.por") %>%
    mutate(prob_scaled = value / rowsum) %>%
    select(treso.id.por, treso.id.pos, prob_scaled)
  
  simulated_trips <- select(observed_trips, treso.id.por, enrolment, value) %>%
    group_by(treso.id.por) %>%
    summarise(enrolment = sum(enrolment)) %>%
    left_join(t, by = "treso.id.por") %>%
    mutate(enrolment = enrolment * prob_scaled) %>%
    left_join(cost, by = c("treso.id.por", "treso.id.pos"))
  
  return(simulated_trips)
}

# Read observed trips and travel time skim with the calculated intra-zonal travel times
read_observed_trips <- function(filepath, school_board_def, treso_zone_def, school_sfis_2017, travel_time_skim,
                                panel_id = "Elementary", board_id = "English Public") {
  #'
  #'
  #'
  #'
  
  # Check if panel and board_type_name are valid
  if (!(panel_id %in% c("Elementary", "Secondary")) | !(board_id %in% c("English Public", "English Catholic", "French Public", "French Catholic"))) {
    stop("`panel` or `board_type_name` is not valid!")
  }
  
  observed_trips <- readRDS(filepath) %>%
    # Drop NAs in POR, it didn't overlay on the TRESO shapefile
    drop_na(treso.id.por) %>%
    left_join(select(school_board_def, dsb, board_type_name), by = c("dsb.index" = "dsb")) %>%
    left_join(select(treso_zone_def, treso_id, area, mof_region), by = c("treso.id.por" = "treso_id")) %>%
    left_join(select(school_sfis_2017, sfis, panel), by = "sfis") %>%
    # Filter the user selected panel and board type name
    filter(panel == panel_id, board_type_name == board_id) %>%
    # Join with travel_time_skim
    left_join(travel_time_skim, by = c("treso.id.por", "treso.id.pos")) %>%
    group_by(treso.id.por, treso.id.pos) %>%
    summarise(value = weighted.mean(value, enrolment),
              euclidean.dist = weighted.mean(euclidean.dist, enrolment),
              enrolment = sum(enrolment))
  
  return(observed_trips)
}

generate_tlfd <- function(observed_trips, simulated_trips, max_value=85, bin_size=1) {
  obs_tlfd <- select(observed_trips, treso.id.por, treso.id.pos, enrolment, value) %>%
    mutate(bin = cut(value, seq(0, max_value, bin_size), labels = seq(bin_size, max_value, bin_size))) %>%
    group_by(bin) %>%
    summarise(enrolment = sum(enrolment)) %>%
    transform(., flag = "obs")
  
  sim_tlfd <- select(simulated_trips, treso.id.por, treso.id.pos, enrolment, value) %>%
    mutate(bin = cut(value, seq(0, max_value, bin_size), labels = seq(bin_size, max_value, bin_size))) %>%
    group_by(bin) %>%
    summarise(enrolment = sum(enrolment)) %>%
    transform(., flag = "model")
  
  # Combine into one dataframe and plot TLFD
  combined_tlfd <- rbind(obs_tlfd, sim_tlfd)
  return(combined_tlfd)
}

calculate_school_weight_forecasting <- function(trip_list, school_list_master, eqao_2017, year_id, panel_id, board_id) {
  
  #filter school list based on panel and board type
  school_list_master <- school_list_master %>%
    filter(year == year_id, panel == panel_id, board_type_name == board_id)
  
  print(head(trip_list))
  print(head(school_list_master))
  
  # Calculate school weighting for TRESO zones with multiple schools and export it for
  pos_school <- trip_list %>%
    group_by(treso.id.pos) %>%
    summarise(trips = sum(trips)) %>%
    right_join(select(school_list_master, otg, sfis, treso.id.pos, school.name, dsb.index), by = c("treso.id.pos") ) %>% 
    left_join(select(eqao_2017, eqao.standardized, sfis), by = "sfis") %>%
    mutate(eqao.standardized = replace_na(eqao.standardized, mean(.$eqao.standardized, na.rm=TRUE)))
  
  print(head(pos_school))
  
  # get treso zone otg totals
  pos_otg_total <- pos_school %>%
    group_by(treso.id.pos) %>%
    summarise(otg.total = sum(otg))
  
  # calculate a combined weight between eqao and otg ratio
  pos_school_weight <- left_join(pos_school, pos_otg_total, by = "treso.id.pos") %>%
    mutate(school.weight = eqao.standardized * (otg / otg.total)) 
  
  # get treso zone school.weight totals
  pos_school_weight_total <- pos_school_weight %>%
    group_by(treso.id.pos) %>%
    summarise(school.weight.total = sum(school.weight))
  
  # school weight
  pos_school_weight_new <- left_join(pos_school_weight, pos_school_weight_total, by = "treso.id.pos") %>%
    mutate(school.weight.prob = school.weight / school.weight.total) %>%
    select(treso.id.pos, sfis, school.name, dsb.index, school.weight.prob) %>%
    group_by(treso.id.pos) %>%
    summarise(
      sfis.list = paste(sfis, collapse = ","),
      school.name.list = paste(school.name, collapse = ","),
      dsb.index.ist = paste(dsb.index, collapse = ","),
      school.weight.prob.list = paste(school.weight.prob, collapse = ",")
    ) %>%
    rowwise() %>%
    mutate(school.weight.prob.list = list(as.numeric(unlist(strsplit(school.weight.prob.list, ","))))) %>%
    mutate(sfis.list = list(as.numeric(unlist(strsplit(sfis.list, ",")))))
  
  return(pos_school_weight_new)
}

calculate_school_weight <- function(observed_trips, school_board_def, school_sfis_2017,
                                    eqao_2017, panel_id = "Elementary", board_id = "English Public") {
  
  # Calculate school weighting for TRESO zones with multiple schools and export it for
  pos_school_EPE <- select(observed_trips, school.name, sfis, treso.id.pos, dsb.index) %>%
    left_join(select(school_board_def, dsb, board_type_name), by = c("dsb.index" = "dsb")) %>%
    left_join(select(school_sfis_2017, sfis, panel), by = "sfis") %>%
    filter(panel == panel_id, board_type_name == board_id) %>%
    group_by(sfis, school.name) %>%
    summarise(treso.id.pos = first(treso.id.pos), dsb.index = first(dsb.index)) %>%
    left_join(select(school_sfis_2017, sfis, otg), by = "sfis") %>%
    left_join(select(eqao_2017, eqao.standardized, sfis), by = "sfis") %>%
    mutate(eqao.standardized = replace_na(eqao.standardized, mean(.$eqao.standardized, na.rm=TRUE))) 
  
  # get treso zone otg totals
  pos_otg_total <- pos_school_EPE %>%
    group_by(treso.id.pos) %>%
    summarise(otg.total = sum(otg))
  
  # calculate a combined weight between eqao and otg ratio
  pos_school_weight <- left_join(pos_school_EPE, pos_otg_total, by = "treso.id.pos") %>%
    mutate(school.weight = eqao.standardized * (otg / otg.total)) 
  
  # get treso zone school.weight totals
  pos_school_weight_total <- pos_school_weight %>%
    group_by(treso.id.pos) %>%
    summarise(school.weight.total = sum(school.weight))
  
  # school weight
  pos_school_weight <- left_join(pos_school_weight, pos_school_weight_total, by = "treso.id.pos") %>%
    mutate(school.weight.prob = school.weight/school.weight.total) %>%
    select(treso.id.pos, sfis, school.name, dsb.index, school.weight.prob) %>%
    group_by(treso.id.pos) %>%
    summarise(
      sfis.list = paste(sfis, collapse = ","),
      school.name.list = paste(school.name, collapse = ","),
      dsb.index.ist = paste(dsb.index, collapse = ","),
      school.weight.prob.list = paste(school.weight.prob, collapse = ",")
    ) %>%
    rowwise() %>%
    mutate(school.weight.prob.list = list(as.numeric(unlist(strsplit(school.weight.prob.list, ","))))) %>%
    mutate(sfis.list = list(as.numeric(unlist(strsplit(sfis.list, ",")))))
  
  return(pos_school_weight)
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
  
  return(school_tb)
}


convertDegrees <- function(deg) {
  #' Convert degrees to radians
  #' 
  #' @param deg The degree value to be converted
  #' @return radian equivalent of the input
  return(deg * pi / 180)
}

haverFunctionEuclidean <- function(lat1, lon1, lat2, lon2) {
  #' Haversine / Euclidean Distance
  #' 
  #' Calculates the euclidean (or crow's fly) distance between two points
  #' 
  #' @param lat1 Latitude value of the first point
  #' @param lon1 Longitude value of the first point
  #' @param lat2 Latitude value of the second point
  #' @param lon2 Longitude value of the second point
  #' @return Euclidean distance in KM between the two points

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

haverFunctionManhattan <- function(lat1, lon1, lat2, lon2) {
  #' Haversine / Manhattan Distance
  #' 
  #' Calculates the manhattan (or grid) distance between two points
  #' 
  #' @param lat1 Latitude value of the first point
  #' @param lon1 Longitude value of the first point
  #' @param lat2 Latitude value of the second point
  #' @param lon2 Longitude value of the second point
  #' @return Manhattan distance in KM between the two points
  
  # Determine 'corner' point in Pythagorean triangle by assiging lat, lon from other points
  lat3 <- lat1
  lon3 <- lon2
  
  d1 <- haverFunction(lat1, lon1, lat3, lon3)
  d2 <- haverFunction(lat2, lon2, lat3, lon3)
  
  # Total Distance
  d = d1 + d2
  
  return(d)
}

balance <- function(matrix, tot, axis) {
  #' Balance function
  #' 
  #' @param matrix The matrix to be balanced
  #' @param tot The vector to be balanced against
  #' @param axis The direction of the balance, 1 indicates rows, 2 indicates columns
  #' @return A balanced matrix

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
  #' Error calculation function
  #' 
  #' @param matrix The matrix
  #' @param a The row vector
  #' @param b The column vector
  #' @return The total difference between the row sum of the matrix and the row vector and the column sum of the matrix and column vector

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
    print("`weight` was not supplied, 1D balancing done without any destination constraints")
    tot = a
    matrix2 <- balance(matrix, tot, axis)
  } else {
    print("`weight` was supplied, 1D balancing done with destination constraints")
    tot = a
    
    # If constrained, scale the matrix to the values of the weight vector
    # It not constrained, apply a boolean value to the matrix
    if (!(constrained)) {
      print("Weight is used as a boolean control for destination")
      weight = ifelse(weight == 0.0, 0, 1)
    } else {
      print("Weight is used directly")
    }
    
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
    
    # Scale the propensity matrix to the weight vector
    matrix_scaled = sweep(matrix, MARGIN = axis_to_scale, weight, `*`)
    # Normalize the propensity matrix to ensure the scaled probailitiy does not exceed 1
    matrix_norm = sweep(matrix_scaled, MARGIN = axis_to_scale, sum, `/`)
    
    # 1D balance against the tot vector on the specified axis
    matrix2 <- balance(matrix_scaled, tot, axis)
  }
  
  return(matrix2)
}

matrix_balancing_2d <- function(matrix, a, b, totals_to_use="raise", max_iterations=10000, rel_error=0.0001) {
  #' Two dimensional balances a matrix
  #' 
  #' Two dimensional matrix balancing with the option to scale the rows or columns to match one another or the average of the two.
  #' The user can control the number of iterations and the relative error this function uses to terminate the iterative procedure. 
  #' 
  #' @param matrix The matrix to be balanced
  #' @param a The row vector to be balanced against
  #' @param b The column vector to be balanced against
  #' @param totals_to_use A flag to determine which totals to use, it could be "row", "column" or "average"
  #' @param max_iterations The maximum number of iterations this procedure should run for
  #' @param rel_error The stopping threshold for this procedure
  #' @return Two dimentionally balanced matrix

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
    
    print(paste0("Iteration: ", i))
  }
  return(matrix2)
}

calculate_simulated_trips <- function(observed_trips, travel_time, alpha, beta) {
  #' Calculate the simulated trips based on alpha and beta values
  #' 
  #' The cost function between every origin and destination pair is computed based on the equation: \eqn{cost = t^{\alpha} * e^{\beta * t|}
  #' 
  #' @param observed_trips The observed trip list 
  #' @param travel_time The travel time matrix
  #' @param alpha The alpha value to compute the cost function
  #' @param beta The beta value to compute the cost function
  #' @return Simulated trips
  
  cfunc <- travel_time %>%
    mutate(value = value^alpha * exp(beta*value)) %>%
    replace_na(value = 0.001)
  
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
    left_join(travel_time, by = c("treso.id.por", "treso.id.pos"))
  
  return(simulated_trips)
}

# TODO check that renaming school_panel_def didn't mess anything up
read_observed_trips <- function(filepath, school_board_def, treso_zone_def, school_panel_def, travel_time_skim,
                                panel_id = "", board_id = "") {
  #' Read the observed trip list and return the trips of the specified panel and board type with travel time, travel distance and enrolment
  #' 
  #' @param filepath The filepath to the observed trip list
  #' @param school_board_def The school board definition file to convert DSB index to board types
  #' @param treso_zone_def The treso zone definition file
  #' @param school_panel_def The school panel defintion file based on sfis
  #' @param travel_time_skim The travel time skim 
  #' @param panel_id The panel id string
  #' @param board_id The board id string
  #' @return Observed trip list with 
  
  # Check if panel and board_type_name are valid
  if (!(panel_id %in% c("Elementary", "Secondary")) | !(board_id %in% c("English Public", "English Catholic", "French Public", "French Catholic"))) {
    stop("`panel` or `board_type_name` is not valid!")
  }
  
  observed_trips <- readRDS(filepath) %>%
    # Drop NAs in POR, it didn't overlay on the TRESO shapefile
    drop_na(treso.id.por) %>%
    left_join(select(school_board_def, dsb, board_type_name), by = c("dsb.index" = "dsb")) %>%
    left_join(select(treso_zone_def, treso_id, area, mof_region), by = c("treso.id.por" = "treso_id")) %>%
    left_join(select(school_panel_def, sfis, panel), by = "sfis") %>%
    # Filter the user selected panel and board type name
    filter(panel == panel_id, board_type_name == board_id) %>%
    # Join with travel_time_skim
    left_join(travel_time_skim, by = c("treso.id.por", "treso.id.pos")) %>%
    group_by(treso.id.por, treso.id.pos) %>%
    summarise(value = weighted.mean(value, enrolment),
              manhattan.dist = weighted.mean(manhattan.dist, enrolment),
              euclidean.dist = weighted.mean(euclidean.dist, enrolment),
              enrolment = sum(enrolment))
  
  return(observed_trips)
}

trip_mean <- function(zone_trips, calc_type) {
  '
  Calculate mean travel time or distance for students based on TRESO zone trips and associated travel times / distances
  
  Input: Dataframe describing, for all relevant zone pairs: number of trips; enrolment; travel time or distance
  Ouput: Mean student travel time or distance
  '
  if(tolower(calc_type) == 'time')
  {
    meanVal <- zone_trips %>%
      ungroup() %>%
      summarise(mean_time = sum(value * enrolment) / sum(enrolment))
  }
  
  if(tolower(calc_type) == 'distance')
  {
    meanVal <- zone_trips %>%
      ungroup() %>%
      summarise(mean_distance = sum(distance * enrolment) / sum(enrolment))
  }
  
  return(meanVal)
}

trip_percentile <- function(zone_trips, calc_type, percentile) {
  '
  Calculate the x percentile travel time or distance for students based on TRESO zone trips and associated travel times / distance
  
  Input: Dataframe describing, for all relevant zone pairs: number of trips; enrolment; travel time or distance
  Ouput: x Percentile student travel time or distance
  
  '
  if(tolower(calc_type) == 'time')
  {
    percentileVal <- zone_trips %>% 
      ungroup() %>%
      arrange(value) %>% 
      mutate(cumSumEnrol = cumsum(enrolment)) %>%
      mutate(enrol90 = (0.9 * max(cumSumEnrol))) %>% 
      filter(cumSumEnrol >= enrol90) %>%
      summarise(timePercent = first(value))
  }
  
  if(tolower(calc_type) == 'distance')
  {
    percentileVal <- zone_trips %>% 
      ungroup() %>%
      arrange(distance) %>% 
      mutate(cumSumEnrol = cumsum(enrolment)) %>%
      mutate(enrolPercent = (percentile * max(cumSumEnrol))) %>% 
      filter(cumSumEnrol >= enrolPercent) %>%
      summarise(distPercent = first(distance))
  }
  
  return(percentileVal)
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

calculate_school_weight_forecasting <- function(trip_list, school_list_master, eqao_2017, new_school=NULL,
                                                year_id, panel_id, board_id) {
  
  # Include the new_school in the master school list
  if(!is.null(new_school)){
    school_list_master <- school_list_master %>% 
      bind_rows(new_school)
  }

  # Filter school list based on panel and board type
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
    select(student.lat, student.long, dist, man.dist, school.name, sfis, dsb.index, panel, enrolment, student.postal.code) %>%
    rename(
      lat = student.lat,
      long = student.long,
      euclidean.dist = dist,
      manhattan.dist = man.dist
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
    summarise_all(list(~first(.))) %>%
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
  This function differs from `create_school_xy` in that this takes in the school dataframe
  without the catchment distance
  
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

create_school_xy_simple <- function(school_sfis) {
  '
  This function differs from `create_school_xy` in that this takes in the school dataframe
  without the catchment distnace. This function differs from `create_school_xy_from_schools` 
  in that it pulls different metadata fields.
  
  input: Dataframe of school with lat long
  output: SpatialPointsDataFrame of each school
  '
  # Convert the school dataframe into SpatialPointsDataframe
  school_spdf <- school_sfis %>%
    ungroup() %>%
    select(sfis, school.lat, school.long) %>%
    mutate(schoolLat = school.lat, schoolLong = school.long) %>%
    mutate(id = row_number())
  coordinates(school_spdf) <- c('schoolLong', 'schoolLat')
  
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
  if (type == 'student') {
    # Find the treso zones which the student points layover
    overlay <- over(xy_location, treso_shp, returnList = FALSE) %>%
      cbind(manhattan.dist = xy_location@data$manhattan.dist,
            euclidean.dist = xy_location@data$euclidean.dist,
            student.postal.code = xy_location@data$student.postal.code,
            school.name = xy_location@data$school.name,
            sfis = xy_location@data$sfis,
            enrolment = xy_location@data$enrolment) %>%
      as_tibble() %>%
      select(Treso_ID, manhattan.dist, euclidean.dist, student.postal.code, enrolment, school.name, sfis) %>%
      rename(
        treso.id.por = Treso_ID
      )
  }
  else if (type == 'school') {
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
  else if (type == 'schoolSimple') {
    # Find the treso zones which the school points layover
    overlay <- over(xy_location, treso_shp, returnList = FALSE) %>%
      cbind(school.lat = xy_location@data$school.lat,
            school.long = xy_location@data$school.long,
            sfis = xy_location@data$sfis) %>%
      as_tibble() %>%
      select(Treso_ID, sfis, school.lat, school.long) %>%
      rename(
        treso.id.pos = Treso_ID
      )
  }
  else if (type == 'marker') {
    overlay <- over(xy_location, treso_shp, returnList = FALSE) %>% 
      cbind(., school.name = xy_location@data$school.name,
            year = xy_location@data$year,
            sfis = xy_location@data$sfis,
            board_type_name = xy_location@data$board_type_name,
            panel = xy_location@data$panel,
            otg = xy_location@data$otg) %>% 
      as_tibble() %>%
      mutate(school.name = as.character(school.name), board_type_name = as.character(board_type_name), panel = as.character(panel)) %>% 
      select(Treso_ID, school.name, year, sfis, board_type_name, panel, otg) %>% 
      rename(treso.id.pos = Treso_ID)
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
  
  # Then, combine with data  that is calculated with sum
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

# Calculate the Geographic Adjustment Factor of a new school based on its enrolment, panel, and location
school_gaf <- function(row, gaf_lookup) {
  '
  This function calculates the geographic adjustment factor for a school based on its location.

  inputs: list of new schools containing postal code by school, and list of geographic adjustment factors by postal code
  output: GAF of school for all schools in list
  '
  # Using function to calculate GAF on a per-row basis
  user_input_pcode <- row['user_input_pcode']
  
  # Determine number of digits for GAF lookup in postal_code lookup table
  for (pcode_length in 3:6) {
    user_pcode <- substr(toupper(user_input_pcode), 1, pcode_length)
    colname <- noquote(paste0('pcode_',pcode_length,'_digit'))
    gafnum <- gaf_lookup %>%
      filter(grepl(user_pcode, !!as.symbol(colname))) %>%
      pull()
    
    if (length(gafnum) > 0) {
      i <- pcode_length
    }
  }
  
  # Finding appropriate GAF based on correct number of postal code digits
  user_pcode <- substr(toupper(user_input_pcode), 1, i)
  colname <- noquote(paste0('pcode_',i,'_digit'))
  
  gaf <- gaf_lookup %>%
    filter(!!as.symbol(colname) == user_pcode) %>% 
    select(gaf) %>% 
    pull()
  
return(gaf)
}  

# Calculate the inflated cost of construction of one or more facilities based on time horizon
construction_cost <- function(user_input_inflation, user_input_scenario_year, currency_year, current_year, new_facility_list) {
  '
  This function calculates the annualized, inflated construction cost for a list of new facilities, and sums the 
  inflated costs together. The function assumes total capital cost is split equally over construction years, and
  that inflation is constant throughout the time horizon. The total cost presented is a sum of nominal dollars from
  each year, in keeping with Treasury Board budget estimates, though in reality the unit of currency is therefore an
  amalgam of dollars belonging to each year in the time horizon.

  inputs: user provided inflation rate; year of scenario being tested; base year of cost estimates; current year upon
    which construction timelines are set; list of new facilities to be built
  output: list of spending in nominal dollars per year, and in total, to build the new facilities
  '

  # Determine min and max years for building inflation index
  min_year = min(user_input_scenario_year, current_year) 
  max_year = max(user_input_scenario_year, current_year)
  
  # "current_year + 1" is used assuming construction wouldn't begin until next year for any FUTURE scenarios
  if (user_input_scenario_year > current_year) {
    min_year = min(user_input_scenario_year, current_year+1) 
    max_year = max(user_input_scenario_year, current_year+1)
  }
  
  # Number of years to use in 'amortization' calc below
  year_count = max_year - min_year + 1
  
  # Calculate inflation index, cost factors
  compounded_inflation <- tibble(index_year = c(seq(min_year,max_year,1))) %>% 
    mutate(inflation_factor = (1 + user_input_inflation)^(index_year - currency_year)) %>% 
    mutate(annual_cost_index_uninflated = 1 / year_count) %>% 
    mutate(annual_cost_index_inflated = inflation_factor / year_count) %>%
    mutate(total_cost_index_inflated = sum(annual_cost_index_inflated))
  
  # Assign cost index to each new facility for each year in the scenario
  scen_cost <- crossing(new_facility_list, compounded_inflation)
  
  # Calculate actual annualized inflated costs, and sum of same
  annualized_inflated_cost <- scen_cost %>% 
    mutate(annual_cost_inflated = cost_2018 * annual_cost_index_inflated) %>% 
    group_by(index_year) %>% 
    mutate(total_yearly_cost_inflated = sum(annual_cost_inflated)) %>%
    ungroup() %>% 
    mutate(total_cost = sum(annual_cost_inflated))
  
  return(annualized_inflated_cost)
} 

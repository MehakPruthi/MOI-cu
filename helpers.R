# unnest in tidyr 1.0.0 currently has slow-down issues
# https://github.com/tidyverse/tidyr/issues/694
# So using the legacy unnest for now.
if (exists("unnest_legacy", where="package:tidyr", mode="function")) {
  unnest <- unnest_legacy
}

if (exists("nest_legacy", where="package:tidyr", mode="function")) {
  nest <- nest_legacy
}

# Common Helpers -----
write_excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  #' Copies the entire dataframe into clipboard
  write.table(x, "clipboard", sep="\t", row.names=row.names, col.names=col.names,...)
}

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
  
  d1 <- haverFunctionEuclidean(lat1, lon1, lat3, lon3)
  d2 <- haverFunctionEuclidean(lat2, lon2, lat3, lon3)
  
  # Total Distance
  d = d1 + d2
  
  return(d)
}

smart_round <- function(x, digits = 0) {
  #' Bucket rounding
  #' [https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum]
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  # Get the biggest residuals that will add up the difference
  # from the floor round
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

sample_by_row <- function(row) {
  #' Apply `sample()`` row by row
  #' 
  #' @param row A row of data
  
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

matrix_balancing_1d <- function(matrix, a, weight, axis=1, is_weight_continuous=TRUE) {
  #' One dimensional balances a matrix
  #' 
  #' One dimensional matrix balancing with optional scaling. If there are known scale weights that needs to be
  #' applied to the propensity matrix, it can be supplied in argument `weight`. The `axis` argument determine which
  #' direction to balance the matrix, 1 indicates rows (origin) and 2 indicates columns (destination). The 
  #' `is_weight_continuous` flag determines if the weighting vector is used as a hard limit or as a boolean value to 
  #' prevent trips to be sent. 
  #' 
  #' @param matrix The propensity matrix
  #' @param a The totals to balance against
  #' @param weight The weight vector to scale against
  #' @param axis The direction to perform the balance
  #' @param is_weight_continuous A flag to determine if the weight is to be used as a constraint or not
  #' @return One dimentionally balanced matrix
  
  # Check if `axis` is 1 or 2
  if (!(axis %in% c(1, 2))) {
    stop("`axis` value is invalid, not one of (1, 2)")
  }
  
  # Check if weight vector is supplied
  if (missing(weight)) {
    print("`weight` was not supplied, 1D balancing done without any destination constraints")
    tot = a
    matrix2 <- balance(matrix, tot, axis)
  } else {
    print("`weight` was supplied, 1D balancing done with destination constraints")
    tot = a
    
    # If is_weight_continuous is TRUE, scale the matrix to the values of the weight vector
    # It is_weight_continuous is FALSE, apply a boolean value to the matrix
    if (!(is_weight_continuous)) {
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

matrix_balancing_2d <- function(matrix, a, b, totals_to_use="raise", max_iterations=10000, rel_error=0.007) {
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
    print(paste0("Error: ", error))
  }
  return(matrix2)
}

# Transportation Related Helpers -----

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

trip_mean <- function(trip_df, calc_type="time") {
  #' Calculate mean travel time or distance for trips
  #' 
  #' @param trip_df A dataframe for all zone pairs the number trips, enrolment, travel time or distance
  #' @param calc_type A string of either 'time' or 'distance'
  #' @return The mean value
  
  if (tolower(calc_type) == 'time') {
    meanVal <- trip_df %>%
      ungroup() %>%
      summarise(mean.time = sum(value * enrolment) / sum(enrolment))
  } else if (tolower(calc_type) == 'distance') {
    meanVal <- trip_df %>%
      ungroup() %>%
      summarise(mean.distance = sum(distance * enrolment) / sum(enrolment))
  } else {
    stop("Expecting `calc_type` to be a string of 'time' or 'distance'")
  }
  
  return(meanVal)
}

trip_percentile <- function(trip_df, calc_type="time", percentile) {
  #' Calculate percentile travel time or distance for trips
  #' 
  #' @param trip_df A dataframe for all zone pairs the number trips, enrolment, travel time or distance
  #' @param calc_type A string of either 'time' or 'distance'
  #' @param percentile A double that represent the percentile the user wish to calculate
  #' @return The percentile value

  if (tolower(calc_type) == 'time') {
    percentileVal <- trip_df %>% 
      ungroup() %>%
      arrange(value) %>% 
      mutate(cumSumEnrol = cumsum(enrolment)) %>%
      mutate(enrol90 = (0.9 * max(cumSumEnrol))) %>% 
      filter(cumSumEnrol >= enrol90) %>%
      summarise(timePercent = first(value))
  } else if (tolower(calc_type) == 'distance') {
    percentileVal <- trip_df %>% 
      ungroup() %>%
      arrange(distance) %>% 
      mutate(cumSumEnrol = cumsum(enrolment)) %>%
      mutate(enrolPercent = (percentile * max(cumSumEnrol))) %>% 
      filter(cumSumEnrol >= enrolPercent) %>%
      summarise(distPercent = first(distance))
  } else {
    stop("Expecting `calc_type` to be a string of 'time' or 'distance'")
  }
  
  return(percentileVal)
}

generate_tlfd <- function(observed_trips, simulated_trips, max_value=85, bin_size=1) {
  #' Combine the observed trips and simulated trips into a binned dataframe. Used to plot trip length frequency diagrams
  #' 
  #' @param observed_trips A dataframe of observed trips
  #' @param simulated_trips A dataframe of simulated trips
  #' @param max_value An integer specifying the maximum binned value 
  #' @param bin_size An numeric specifying the bin size
  #' @return A dataframe with both types of trips and binned
  
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

plot_straight_line_tlfd <- function(tlfd) {
  #' Plot TLFD with segmentation
  
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

plot_travel_time_tlfd <- function(tlfd) {
  #' Plot TLFD with segmentation

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

create_overlay <- function(xy_location, treso_shp, type='student') {
  #' This function takes the SpatialPointsDataFrame and maps the XY location of the object on the 
  #' TRESO shapefile and returns the appropriate TRESO zone ID the objects fall on top of.
  #' 
  #' @param xy_location A SpatialPointsDataFrame of the object
  #' @param treso_shp A Shapefile of the TRESO zones
  #' @param type A string indicating the type of the object in order to pull the correct metadata fields
  #' @return A database with the TRESO zone ID for each object
  
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
  } else if (type == 'school') {
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
  } else if (type == 'schoolSimple') {
    # Find the treso zones which the school points layover
    overlay <- over(xy_location, treso_shp, returnList = FALSE) %>%
      cbind(school.lat = xy_location@data$school.lat,
            school.long = xy_location@data$school.long,
            sfis = xy_location@data$sfis) %>%
      as_tibble() %>%
      select(Treso_ID, sfis, school.lat, school.long) %>%
      rename(treso.id.pos = Treso_ID)
  } else if (type == 'court') {
    # Find the treso zones which the school points layover
    overlay <- over(xy_location, treso_shp, returnList = FALSE) %>%
      cbind(lat = xy_location@data$courthouse.lat,
            long = xy_location@data$courthouse.long,
            bid = xy_location@data$bid) %>%
      as_tibble() %>%
      mutate(bid = as.character(bid)) %>%  
      select(Treso_ID, bid, lat, long) %>%
      rename(treso.id.pos = Treso_ID)
  } else if (type == 'hospital') {
    # Find the treso zones which the hospital points layover
    overlay <- over(xy_location, treso_shp, returnList = FALSE) %>%
      cbind(lat = xy_location@data$hospital.lat,
            long = xy_location@data$hospital.long,
            id = xy_location@data$id) %>%
      as_tibble() %>%
      mutate(bid = as.character(id)) %>%  
      select(Treso_ID, id, lat, long) %>%
      rename(treso.id.pos = Treso_ID)
  } else if (type == 'marker') {
    overlay <- over(xy_location, treso_shp, returnList = FALSE) %>% 
      cbind(., school.name = xy_location@data$school.name,
            year = xy_location@data$year,
            sfis = xy_location@data$sfis,
            dsb.index = xy_location@data$dsb.index,
            school.lat = xy_location@data$school.lat,
            school.long = xy_location@data$school.long,
            board_type_name = xy_location@data$board_type_name,
            board.name = xy_location@data$board.name,
            panel = xy_location@data$panel,
            otg = xy_location@data$otg) %>% 
      as_tibble() %>%
      mutate(year = as.integer(as.character(year)), school.name = as.character(school.name),
             board.name = as.character(board.name),
             board_type_name = as.character(board_type_name), panel = as.character(panel)) %>% 
      select(Treso_ID, year, dsb.index, school.name, school.lat, school.long, sfis, board.name, board_type_name, panel, otg) %>% 
      rename(treso.id.pos = Treso_ID)
  } else {
    stop("Expecting `type` to be either 'student', 'school', 'schoolSimple', 'court', 'hospital' or 'marker'.")
  }
  return(overlay)
}

# EDU -----
## EDU Cleaning ----

create_student_xy <- function(student_travel, treso_shp) {
  #' This function transforms each student location into the same projection as the TRESO shapefiles
  #' 
  #' @param student_travel A dataframe of the students and schools within a certain catchment distance of the school
  #' @param treso_shp A Shapefile of TRESO zones
  #' @return A SpatialPointsDataFrame 

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

create_school_xy <- function(student_travel, treso_shp) {
  #' This function transforms each school location into the same projection as the TRESO shapefiles
  #' 
  #' @param student_travel A dataframe of the students and schools within a certain catchment distance of the school
  #' @param treso_shp A Shapefile of TRESO zones
  #' @return A SpatialPointsDataFrame 
  
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

buffer_zones <- function(school_xy, treso_shp) {
  #' Performs spatial buffers to select the TRESO zones within the catchment distance of a school.
  #' 
  #' @param school_xy A SpatialPointsDataFrame of the schools
  #' @param treso_shp A Shapefile of TRESO zones
  #' @return A Dataframe of TRESO zones in close proximity to each school
  
  buffered_points <- gBuffer(school_xy, width = school_xy@data$catchment.dist * 1000, byid = TRUE)
  
  tic('Multi-Core Buffering')
  # Start multi-core clustering for performance
  # if you do not remove at lease one core the machine will lock until model has run
  cores = parallel::detectCores()
  cl <- parallel::makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  
  # `foreach()` loop to find the TRESO zones touching the school buffers
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

summarize_buffered_zones <- function(buffered_df, treso_tb, school_ade, school_board_def, treso_zone_def) {
  #' Summarize the socio-economic information of the buffered TRESO zones for each school.
  #' 
  #' @param buffered_df A Dataframe of the buffered TRESO zones
  #' @param treso_tb A Dataframe of TRESO socio-economic information
  #' @param school_ade A Dataframe of the school's ADE information
  #' @param school_board_def A Dataframe with the school board definitions
  #' @param treso_zone_def A Dataframe with the TRESO zone definitions
  #' @return A Dataframe of the summarized socio-economic information for each school

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

school_gaf <- function(row, gaf_lookup) {
  #' Calculate the Geographic Adjustment Factor of a new school based on its enrolment, panel, and location
  #' 
  #' @param row A row in a Dataframe that contains the postal code by school
  #' @param gaf_lookup A Dataframe of geographic adjustment factors by postal code
  #' @return GAF of school for all schools in list
  
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

construction_cost <- function(user_input_inflation, user_input_scenario_year, currency_year, current_year, new_facility_list) {
  #' Calculates the annualized, inflated construction cost for a list of new facilities, and sums the 
  #' inflated costs together. The function assumes total capital cost is split equally over construction years, and
  #' that inflation is constant throughout the time horizon. The total cost presented is a sum of nominal dollars from
  #' each year, in keeping with Treasury Board budget estimates, though in reality the unit of currency is therefore an
  #' amalgam of dollars belonging to each year in the time horizon.
  #' 
  #' @param user_input_inflation A numeric value that represents the inflation rate
  #' @param user_input_scenario_year An integer that represents the forecast scenario year
  #' @param currency_year An integer that represents the base year of cost esimates
  #' @param current_year An integer that represents the current year upon which construction timelines are set
  #' @param new_facility_list A list of new facilities to be built
  #' @return A list of spending in nominal dollars per year, and in total, to build the new facilities

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

# EDU Model -----
create_pos_vector <- function(df, new_school, full_vector, year_id, panel_id, board_id) {
  #' Creates a POS vector from the input dataframe
  #' 
  #' @param df The input dataframe, expects the master school list
  #' @param new_school Input dataframe of new schools added by user/decision-making layer
  #' @param full_vector The full list of TRESO zones in order
  #' @param year_id Integer indicating the modelling year
  #' @param panel_id String indicating the modelling panel
  #' @param board_id String indicating the modelling board type
  #' @return a data matrix
  #' 

  # Include the new_school in the master school list
  if(!is.null(new_school)){
    df <- df %>% 
      bind_rows(new_school)
  }
  
  pos <- df %>% 
    filter(year == year_id, panel == panel_id, board_type_name == board_id) %>%
    select(treso.id.pos, otg) %>%
    group_by(treso.id.pos) %>%
    summarise(otg = sum(otg)) %>%
    right_join(full_vector, by=c("treso.id.pos" = "dest")) %>%
    replace_na(list(otg = 0)) %>%
    arrange(treso.id.pos) %>%
    column_to_rownames(var = "treso.id.pos") %>%
    data.matrix()

  return(pos)
}

create_por_forecast_vector <- function(df, full_vector, panel_id, board_id) {
  #' Creates a POR vector from the input dataframe
  #' 
  #' @param df The input dataframe, expects the forecasted treso population
  #' @param full_vector The full list of TRESo zones in order
  #' @param panel_id String indicating the modelling panel
  #' @param board_id String indicating the modelling board type
  #' @return a data matrix
  #' 
  por <- df %>% 
    filter(panel == panel_id, board_type_name == board_id) %>% 
    select(treso.id.por, potential.enrolment) %>% 
    right_join(full_vector, by=c("treso.id.por" = "orig")) %>% 
    replace_na(list(potential.enrolment = 0)) %>% 
    arrange(treso.id.por) %>% 
    column_to_rownames(var = "treso.id.por") %>% 
    data.matrix()
  
  return(por)
}

apply_sampling_to_population <- function(forecast_population, board_type_sample) {
  #' Apply the board type sample to the forecasted population dataframe
  #' 
  #' @param forecast_population The dataframe with forecasted population information
  #' @param board_type_sample The dataframe with the board type sampling probabilities for each CSD and Panel
  #' @return A dataframe with population segmented by panel and board type
  #' 
  set.seed(42)
  forecast_population_by_board <- forecast_population %>%
    left_join(board_type_sample, by = c("panel", "cduid")) %>%
    filter(!is.na(cduid)) %>%
    # Create probability list for sample
    unite(prob, `English Catholic`, `English Public`, `French Catholic`, `French Public`, sep = ",") %>%
    rowwise() %>%
    mutate(prob = list(as.double(unlist(strsplit(prob, ","))))) %>%
    # Sample the different board types 
    mutate(sample.result = list(sample(c("EC", "EP", "FC", "FP"), size=potential.enrolment, replace=TRUE, prob=prob))) %>%
    mutate(`English Catholic` = sum(sample.result == "EC"),
           `English Public` = sum(sample.result == "EP"),
           `French Catholic` = sum(sample.result == "FC"),
           `French Public` = sum(sample.result == "FP")) %>% 
    select(treso.id.por, panel, cduid, `English Catholic`:`French Public`) %>%
    gather(key="board_type_name", value="potential.enrolment", `English Catholic`:`French Public`)
  
  return(forecast_population_by_board)
}

calculate_school_weight_forecasting <- function(trip_list, school_list_master, eqao_2017, new_school,
                                                year_id, panel_id, board_id) {
  #' Calculate weight factor for distributing students between two or more schools within a single TRESO zone
  #' 
  #' @param trip_list
  #' @param school_list_master
  #' @param eqao_2017
  #' @param new_school
  #' @param year_id
  #' @param panel_id
  #' @param board_id
  #' @return
  
  # Include the new_school in the master school list
  if(!is.null(new_school)){
    school_list_master <- school_list_master %>% 
      bind_rows(new_school)
  }
  
  # Filter school list based on panel and board type
  school_list_master <- school_list_master %>%
    filter(year == year_id, panel == panel_id, board_type_name == board_id)
  
  # Calculate school weighting for TRESO zones with multiple schools and export it for
  pos_school <- trip_list %>%
    group_by(treso.id.pos) %>%
    summarise(trips = sum(trips)) %>%
    right_join(select(school_list_master, otg, sfis, treso.id.pos, school.name, dsb.index), by = c("treso.id.pos") ) %>% 
    left_join(select(eqao_2017, eqao.standardized, sfis), by = "sfis") %>%
    mutate(eqao.standardized = replace_na(eqao.standardized, mean(.$eqao.standardized, na.rm=TRUE)))
  
  # calculate a combined weight between eqao and otg ratio
  pos_school_weight <- pos_school %>%
    group_by(treso.id.pos) %>% 
    mutate(school.weight = eqao.standardized * (otg / sum(otg)), 
           school.weight.prob = school.weight / sum(school.weight))
  
  # school weight
  pos_school_weight_new <- pos_school_weight %>%
    select(treso.id.pos, sfis, school.name, dsb.index, school.weight.prob) %>%
    group_by(treso.id.pos) %>%
    summarise(
      sfis.list = paste(sfis, collapse = ","),
      school.name.list = paste(school.name, collapse = ","),
      dsb.index.list = paste(dsb.index, collapse = ","),
      school.weight.prob.list = paste(school.weight.prob, collapse = ",")
    ) %>%
    rowwise() %>%
    mutate(school.weight.prob.list = list(as.numeric(unlist(strsplit(school.weight.prob.list, ","))))) %>%
    mutate(sfis.list = list(as.numeric(unlist(strsplit(sfis.list, ",")))))
  
  return(pos_school_weight_new)
}

distribute_students_to_schools <- function(school_master, new_school_df, por, pos_full, prop_matrix, year_id, panel_id, board_id, is_weight_continuous) {
  #' Creates a new POS vector with new schools and distribute the students to the schools
  #' 
  #' @param school_df
  #' @param new_school_df
  #' @param por
  #' @param pos_full
  #' @param prop_matrix
  #' @param panel_id
  #' @param board_id
  #' @param school_diff
  #' @param schools_consolidated_closed
  #' @return A dataframe of schools with simulated ADE 
  
  pos <- create_pos_vector(school_master, new_school_df, pos_full, year_id=year_id, panel_id=panel_id, board_id=board_id)
  
  print(paste0("ADE is: ", sum(por), ". OTG is: ", sum(pos), "."))
  
  prop_matrix_balanced <- matrix_balancing_1d(prop_matrix, por, pos, axis=1, is_weight_continuous=is_weight_continuous)
  trip_list <- reshape2::melt(prop_matrix_balanced) %>%
    arrange(Var1, Var2)
  colnames(trip_list) <- c("treso.id.por", "treso.id.pos", "trips")
  
  # Distribute the ADE at each zone to the individual schools of the zone
  school_summary_df_list <- forecast_school_ade(prop_matrix, trip_list, school_master, eqao_2017, new_school_df,
                                                year_id, panel_id, board_id)
  
  school_summary <- school_summary_df_list[[1]]
  
  return(school_summary)
}

distribute_students_within_zones <- function(school_forecast_df) {
  #' Create the forecasted school list with OTG, OTG Threshold, ADE and Simulated ADE
  #' Also redistribute any overfilled schools in a zone to 
  #' 
  #' @param school_forecast_df
  #' @return A dataframe of schools with ADE distributed among underfilled schools in the TRESO zone 
  
  school_forecast_distributed_df <- school_forecast_df %>% 
    mutate(ade.diff = simulated.ade - otg.threshold,
           overfill = pmax(ade.diff, 0),
           underfill = pmin(ade.diff, 0)) %>% 
    # Create flag to determine if there is capacity remaining in each school
    mutate(capacity.flag = 1 - pmax(pmin(ade.diff, 1), 0)) %>%
    group_by(treso.id.pos) %>% 
    mutate(treso.capacity = sum(underfill),
           treso.overfill = sum(overfill)) %>%
    # Calculate likelihood of going each of the other school(s) in the TRESO zone if redirected from full school
    mutate(otg.weight.underfill = otg * capacity.flag,
           underfill.share = otg.weight.underfill / sum(otg.weight.underfill),
           underfill.share = replace_na(underfill.share, 0)) %>%
    # Calculate likelihood of coming from each overfilled school in a TRESO zone
    mutate(otg.weight.overfill = overfill * (1 - capacity.flag),
           overfill.share = otg.weight.overfill / sum(otg.weight.overfill),
           overfill.share = replace_na(overfill.share, 0)) %>%
    # Determine how many students will be redistributed within zone
    mutate(treso.redist = min(abs(treso.overfill), abs(treso.capacity)),
           school.redist.to = treso.redist * underfill.share,
           school.redist.from = treso.redist * overfill.share) %>% 
    # Adjust simulated ADE to account for shift from one school to another within zones
    mutate(simulated.ade.rev = simulated.ade + school.redist.to - school.redist.from) %>% 
    ungroup() %>% 
    select(treso.id.pos, sfis, school.name, otg, otg.threshold, ade, simulated.ade = simulated.ade.rev)
  
  
  return(school_forecast_distributed_df)
}

create_forecast_school_list <- function(school_df, new_school_df, panel_id, board_id, school_diff, schools_consolidated_closed) {
  # TODO remove schools_consolidated_closed when it is in the original school_df
  
  #' Create the forecasted school list with OTG, OTG Threshold, ADE and Simulated ADE
  #' Also redistribute any overfilled schools in a zone to 
  #' 
  #' @param school_df
  #' @param new_school_df
  #' @param panel_id
  #' @param board_id
  #' @param school_diff
  #' @param schools_consolidated_closed
  #' @return A dataframe of schools with forecasted/simulated ADE and 
  
  school_forecast_df <- school_df %>% 
    filter((status == "Open" | is.na(status)), otg != 0, ade != 0) %>% 
    filter(panel == panel_id, board_type_name == board_id) %>% 
    # Add in new schools
    bind_rows(new_school_df) %>% 
    full_join(select(school_diff, sfis, simulated.ade.20xx, simulated.ade.base, change.ade), by="sfis") %>%
    # Anti-join to remove any schools which existed in 2017 but no longer exist in future year
    anti_join(select(schools_consolidated_closed, sfis), by = c('sfis')) %>% 
    replace_na(list(ade = 0, change.ade = 0)) %>%                 
    # 'Actual' forecast ADE = existing ADE (2017 actual historical data) plus change in ADE estimated in Step 3
    mutate(simulated.ade.raw = ade + change.ade) %>% 
    # If simulated.ade.raw < 0 for any school, this value is rounded up to 0 to prevent having a negative number of students at each school.
    # However, by rounding up to 0, the model could be adding 'phantom' students to the system
    # So the total ADE as estimated by the distribution model will exceed total ADE estimated in population forecasts.
    # However, in test runs, this did not occur at any schools, so risk appears low. If it were to happen, the result would be a slight increase in the total number of ADE across the province, which would almost certainly be negligible. 
    mutate(simulated.ade = ifelse(ade + change.ade < 0, 0, ade + change.ade)) %>%
    # Create OTG threshold based on user input
    mutate(otg.threshold = otg * USER_OTG_THRESHOLD) %>% 
    select(treso.id.pos, sfis, school.name, otg, otg.threshold, ade, simulated.ade)
  
  # Redistribute students from overfilled schools to underfilled schools in the same zone
  school_forecast_distributed_df <- distribute_students_within_zones(school_forecast_df)
  
  return(school_forecast_distributed_df)
}
  
forecast_school_ade <- function(prop_matrix, trip_list, school_master, eqao_2017, new_school, year_id, panel_id, board_id) {
  #' Produce a dataframe with the summary of the school's forecasted ADE
  #'
  #' @param prop_matrix
  #' @param trip_list The trip list from the balanced matrix
  #' @param school_master
  #' @param eqao_2017
  #' @param new_school
  #' @param year_id
  #' @param panel_id
  #' @param board_id
  #' @return A dataframe of schools with forecasted ADE
  #' 
  # Calculate the school weight
  pos_school_weight <- calculate_school_weight_forecasting(trip_list, school_master, eqao_2017, new_school = new_school,
                                                           year_id, panel_id, board_id)
  print(paste0("For ", panel_id, "-", board_id, ", there are ",
               nrow(filter(pos_school_weight, length(school.weight.prob.list) == 1)),
               " TRESO zones with a single school, and ",
               nrow(filter(pos_school_weight, length(school.weight.prob.list) > 1)),
               " TRESO zones with multiple schools."))
  
  # Apply bucket rounding to chunks of data by TRESO POS
  results <- by(trip_list$trips, trip_list[c("treso.id.pos")], smart_round, simplify = TRUE)
  
  # Convert the output of by() to a dataframe
  results2 <- sapply(results, I)
  colnames(results2) <- colnames(prop_matrix)
  rownames(results2) <- rownames(prop_matrix)
  
  df <- reshape2::melt(results2) %>%
    arrange(Var1, Var2)
  colnames(df) <- c('treso.id.por', 'treso.id.pos', 'enrolment.rounded')
  
  # Save a list of schools with 0 students assigned
  df_0 <- df %>% 
    group_by(treso.id.pos) %>% 
    summarise(enrolment.rounded = sum(enrolment.rounded)) %>% 
    filter(enrolment.rounded == 0)
  
  # After combining with travel time, trip list can be shortened
  df <- filter(df, enrolment.rounded != 0)
  
  # Combine the school weight and sample using the weight
  trip_list <- left_join(df, pos_school_weight, by=c("treso.id.pos"))
  
  # Apply sampling procedure to assign students to schools in the same TRESO zone
  plan(multiprocess)
  results <- future_apply(trip_list, 1, sample_by_row)
  results_tb <- t(as.data.table(results))
  trip_list_schools <- cbind(trip_list, results_tb)
  
  schools_summary <- trip_list_schools %>% 
    unnest(results_tb) %>%
    mutate(enrol.value = 1) %>%
    rename(sfis = results_tb) %>% 
    group_by(sfis) %>% 
    summarise(simulated.ade = sum(enrol.value))
  
  df_list <- list(schools_summary, df_0)
  
  return(df_list)
  
}

edu_dm <- function(treso_travel_time, trip_list, treso_zone_def, por_additional, travel_time_threshold_factor, zone_proximity_threshold, min_tt_threshold) {
  #' Produce a list of 2 dataframes containing: 
  #' 1. Dataframe of TRESO zones in which to consider building schools
  #' 2. Dataframe of TRESO zones ruled out from building schools due to proximity to higher-ranked location for building school
  #'
  #' @param treso_travel_time A dataframe containing travel times to and from all origin-destination TRESO pairs
  #' @param trip_list The trip list from the balanced matrix
  #' @param treso_zone_def Details of treso zones, e.g., CSDUID/CDUID
  #' @param por_additional A list of TRESO origins which have residual students due to school overfills --> candidate zones for building a new school
  #' @param travel_time_threshold_factor A user-set factor for selecting how much travel time is considered 'excess' time
  #' @param zone_proximity_threshold A user-set threshold to preclude construction of schools in each of two TRESO zones too close together
  #' @param min_tt_threshold A user-set threshold to preclude construction of schools in TRESO zones without sufficient 'excess' travel time
  #' @return A list of 2 dataframes with TRESO zones in which to build, and TRESO zones ruled out due to proximity
  
  # Create full matrix of travel times for students across the province
  potential_zones <- cbind(trip_list, select(treso_travel_time, value)) %>%
    rename(travel.time = value) %>%
    left_join(select(treso_zone_def, treso_id, treso.por.csduid = csduid, treso.por.cduid = cduid), by = c('treso.id.por' = 'treso_id')) %>% 
    left_join(select(treso_zone_def, treso_id, treso.pos.csduid = csduid, treso.pos.cduid = cduid), by = c('treso.id.pos' = 'treso_id')) 

  # Calculate mean travel time across CD and filter out the shorter trips
  threshold_tt_zones <- potential_zones %>%
    group_by(treso.por.cduid) %>%
    mutate(mean.travel.time = sum(trips * travel.time) / sum(trips)) %>%
    # Filtered to remove trips below the travel time threshold set
    mutate(travel.time.threshold = travel_time_threshold_factor * mean.travel.time) %>%
    filter(travel.time >= travel.time.threshold) %>%
    mutate(excess.travel.time = (travel.time - travel.time.threshold) * trips) %>%
    # Sum excess travel time for all destinations (i.e., school locations) for each TRESO zone to get total excess travel time for students living in that TRESO zone
    group_by(treso.id.por) %>%
    mutate(excess.travel.time = sum(excess.travel.time)) %>%
    ungroup() %>%
    arrange(desc(excess.travel.time))
  
  # Calculate zone pairs for consideration due to school overfill
  additional_pairs <- threshold_tt_zones %>%
    right_join(select(por_additional, treso.id.por), by = c('treso.id.por'))

  # Determine shortlist of zones within which construction should be considered by ruling out zones too close to superior zones
  shortlist_zones <- threshold_tt_zones %>%
    # Cut down list of zones under consideration based on minimum bar for 'excess' travel time
    filter(excess.travel.time > min_tt_threshold) %>%
    select(treso.id.por) %>% 
    rbind(select(por_additional, treso.id.por)) %>%
    distinct(treso.id.por)

  # Retain the original copy of the shortlisted zones
  shortlist_zones_dup <- shortlist_zones

  shortlist_pairs <- potential_zones %>%
    select(treso.id.por, treso.id.pos, travel.time) %>%
    inner_join(shortlist_zones, by = c('treso.id.por')) %>%
    inner_join(shortlist_zones, by = c('treso.id.pos' = 'treso.id.por'))

  print(paste0("The number of zones in shortlist is: ", nrow(shortlist_zones)))
  print(paste0("The number of zone pairs in shortlist is: ", nrow(shortlist_pairs)))
  
  # Initialize dataframes and while-loop check
  num_shortlist_zones <- nrow(shortlist_zones)
  build_df <- tibble(zone = integer())
  proximity_df <- tibble(zone = integer())

  start.time <- Sys.time()
  print(Sys.time())

  # Loop through the shortlisted TRESO locations to choose best zones in which to build
  while (num_shortlist_zones > 0) {
    # Append the highest prioirty TRESO zone to the build_df
    build_zone_id = first(shortlist_zones$treso.id.por)
    build_df <- build_df %>%
      add_row(., zone = build_zone_id)
  
    print(paste0('Build Zone: ', build_zone_id))

    # Create a list of zones to anti_join with shortlist_zones to remove close proximity zones
    origin_zones_removed <- shortlist_pairs %>%
      filter(treso.id.por != build_zone_id) %>%
      filter(treso.id.pos == build_zone_id) %>%
      filter(travel.time <= zone_proximity_threshold) %>%
      distinct(treso.id.por)

    # Remove TRESO zones that have been selected or too close to the selected from shortlist_zones
    shortlist_zones <- shortlist_zones %>%
      filter(treso.id.por != build_zone_id) %>%
      anti_join(origin_zones_removed, by = c('treso.id.por'))

    print(paste0("The number of zones in shortlist after proximity check is: ", nrow(shortlist_zones)))

    # Reduce short_list pairs to exlude previously eliminated zones from future consideration
    shortlist_pairs <- shortlist_pairs %>%
      # Keep OD pairs where the origin zone is not the same as build_zone
      filter(treso.id.por != build_zone_id) %>%
      # Keep OD pairs where the destination zone is not the same as build_zone 
      # OR pairs where the destination is the same as build_zone BUT outside of the proximity threshold
      filter(treso.id.pos != build_zone_id | (treso.id.pos == build_zone_id & travel.time > zone_proximity_threshold))

    print(paste0('The number of OD pairs in shortlist after proximity check is: ', nrow(shortlist_pairs)))

    # Update the number of zones in shortlist_zones
    num_shortlist_zones <- nrow(shortlist_zones)
  }
  
  end.time <- Sys.time()
  time.diff <- end.time - start.time
  print(time.diff)
  
  proximity_df <- shortlist_zones_dup %>%
    anti_join(build_df, by = c('treso.id.por' = 'zone')) %>%
    distinct(treso.id.por)

  return(list(build_df, proximity_df))
}

# MOH ----
## MOH Cleaning ----

create_hospital_xy <- function(hospital_master) {
  #' This function transforms each hospital location into the same projection as the TRESO shapefiles
  #' 
  #' @param hospital_master A dataframe of the hospitals
  #' @param treso_shp A Shapefile of TRESO zones
  #' @return A SpatialPointsDataFrame 

  hospital_spdf <- hospital_master %>%
    ungroup() %>%
    select(id, hospital.lat, hospital.long) %>%
    mutate(lat = hospital.lat, long = hospital.long) %>%
    mutate(ref = row_number())
  coordinates(hospital_spdf) <- c('long', 'lat')
  
  # Project the hospital lat/long to TRESO's LCC specification
  proj4string(hospital_spdf) <- CRS('+proj=longlat +datum=WGS84')
  treso_projarg = treso_shp@proj4string@projargs
  
  hospital_xy <- spTransform(hospital_spdf, CRS(treso_projarg))
  
  return(hospital_xy)
}

# MAG -----
## MAG Cleaning ----

create_court_xy <- function(court_master) {
  #' This function transforms each court location into the same projection as the TRESO shapefiles
  #' 
  #' @param court_master A dataframe of the courts
  #' @param treso_shp A Shapefile of TRESO zones
  #' @return A SpatialPointsDataFrame 
  
  court_spdf <- court_master %>%
    ungroup() %>%
    select(bid, courthouse.lat, courthouse.long) %>%
    mutate(lat = courthouse.lat, long = courthouse.long) %>%
    mutate(id = row_number())
  coordinates(court_spdf) <- c('long', 'lat')
  
  # Project the student and school points from lat/long to TRESO's LCC specification
  proj4string(court_spdf) <- CRS('+proj=longlat +datum=WGS84')
  treso_projarg = treso_shp@proj4string@projargs
  
  court_xy <- spTransform(court_spdf, CRS(treso_projarg))
  
  return(court_xy)
}

courtroom_size <- function(courtrooms) {
  #' Determine the appropriate area per courtroom as a function of number of courtrooms in a courthouse
  #' 
  #' @param courtrooms An integer for the number of courtrooms
  #' @return An numeric value of the required area per courtroom
  #' 
  
  # Gross Area for smallest courthouses: 1,600 sqm per courtroom
  courtroom_count_min = 1
  area_max_sqm = 1600
  
  # Gross Area for largest courthouses: 1,200 sqm per courtroom
  courtroom_count_max = 50
  area_min_sqm = 1200
  
  # Gross Area function calculation
  area_diff = area_max_sqm - area_min_sqm
  courtroom_diff = courtroom_count_max - courtroom_count_min
  area_change_per_courtroom = area_diff / courtroom_diff
  
  courtroom_count_scaled = max(min(courtrooms, 50), 1) # Setting area standard to have a minimum courtroom count of 1 and a max of 50 in linear sizing scale
  required_area_per_courtroom = area_max_sqm - courtroom_count_scaled * area_change_per_courtroom
  
  required_area_per_courtroom_sqft = required_area_per_courtroom * 3.28^2
  
  return(required_area_per_courtroom_sqft)
}

# Shiny Helpers ----
getUtilizationColor <- function(value) {
  sapply(value, function(value) {
    if (value <= 0.75) {
      "green"
    } else if (value > 0.75 & value < 1.0) {
      "orange"
    } else {
      "red"
    }
  })
}

utility_rename <- function(x) {
  #' Rename columns by taking in x and adding "courtrooms.needed." as a prefix to x
  #' 
  #' @param x The string to be renamed
  #' @return Renamed string
  name = paste0("courtrooms.needed.", x)
}




















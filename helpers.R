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
  #' Copies the entire Dataframe into clipboard
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
  #' @param trip_df A Dataframe for all zone pairs the number trips, enrolment, travel time or distance
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
  #' @param trip_df A Dataframe for all zone pairs the number trips, enrolment, travel time or distance
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
  #' Combine the observed trips and simulated trips into a binned Dataframe. Used to plot trip length frequency diagrams
  #' 
  #' @param observed_trips A Dataframe of observed trips
  #' @param simulated_trips A Dataframe of simulated trips
  #' @param max_value An integer specifying the maximum binned value 
  #' @param bin_size An numeric specifying the bin size
  #' @return A Dataframe with both types of trips and binned
  
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
  
  # Combine into one Dataframe and plot TLFD
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
  #' This function takes the SpatialPointsDataframe and maps the XY location of the object on the 
  #' TRESO shapefile and returns the appropriate TRESO zone ID the objects fall on top of.
  #' 
  #' @param xy_location A SpatialPointsDataframe of the object
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
  #' @param student_travel A Dataframe of the students and schools within a certain catchment distance of the school
  #' @param treso_shp A Shapefile of TRESO zones
  #' @return A SpatialPointsDataframe 
  
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
  #' @param student_travel A Dataframe of the students and schools within a certain catchment distance of the school
  #' @param treso_shp A Shapefile of TRESO zones
  #' @return A SpatialPointsDataframe 
  
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
  #' @param school_xy A SpatialPointsDataframe of the schools
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
  
  user_input_pcode <- row['postal.code']
  # Initialize i
  i = 0
  
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
  
  # If postal code has no GAF associated, assign GAF = 1.0
  if (i == 0) {
    gaf <- 1
    return(gaf)
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

construction_cost <- function(inflation_rate, scenario_year, currency_year, current_year, new_facility_list) {
  #' Calculates the annualized, inflated construction cost for a list of new facilities, and sums the 
  #' inflated costs together. The function assumes total capital cost is split equally over construction years, and
  #' that inflation is constant throughout the time horizon. The total cost presented is a sum of nominal dollars from
  #' each year, in keeping with Treasury Board budget estimates, though in reality the unit of currency is therefore an
  #' amalgam of dollars belonging to each year in the time horizon.
  #' 
  #' @param inflation_rate A numeric value that represents the inflation rate
  #' @param scenario_year An integer that represents the forecast scenario year
  #' @param currency_year An integer that represents the base year of cost esimates
  #' @param current_year An integer that represents the current year upon which construction timelines are set
  #' @param new_facility_list A list of new facilities to be built
  #' @return A list of spending in nominal dollars per year, and in total, to build the new facilities
  
  # Determine min and max years for building inflation index
  min_year = min(scenario_year, current_year) 
  max_year = max(scenario_year, current_year)
  
  # "current_year + 1" is used assuming construction wouldn't begin until next year for any FUTURE scenarios
  if (scenario_year > current_year) {
    min_year = min(scenario_year, current_year+1) 
    max_year = max(scenario_year, current_year+1)
  }
  
  # Number of years to use in 'amortization' calc below
  year_count = max_year - min_year + 1
  
  # Calculate inflation index, cost factors
  compounded_inflation <- tibble(index_year = c(seq(min_year,max_year,1))) %>% 
    mutate(inflation_factor = (1 + inflation_rate)^(index_year - currency_year)) %>% 
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
    ungroup()
  
  return(annualized_inflated_cost)
} 

# EDU Model -----
create_pos_vector <- function(school_df, new_school, full_vector) {
  #' Creates a POS vector from the input Dataframe
  #' 
  #' @param school_df The input Dataframe, expects a list of school
  #' @param new_school Input Dataframe of new schools added by user/decision-making layer
  #' @param full_vector The full list of TRESO zones in order
  #' @return a data matrix
  #' 
  
  # Include the new_school in the master school list
  if(!is.null(new_school)){
    school_df <- school_df %>% 
      bind_rows(new_school)
  }
  
  pos <- school_df %>% 
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
  #' Creates a POR vector from the input Dataframe
  #' 
  #' @param df The input Dataframe, expects the forecasted treso population
  #' @param full_vector The full list of TRESo zones in order
  #' @param panel_id String indicating the modelling panel
  #' @param board_id String indicating the modelling board type
  #' @return a data matrix
  #' 
  por <- df %>% 
    filter(panel == panel_id, board_type_name == board_id) %>% 
    select(treso.id.por, enrolment) %>% 
    right_join(full_vector, by=c("treso.id.por" = "orig")) %>% 
    replace_na(list(enrolment = 0)) %>% 
    arrange(treso.id.por) %>% 
    column_to_rownames(var = "treso.id.por") %>% 
    data.matrix()
  
  return(por)
}

apply_sampling_to_population <- function(forecast_population, board_type_sample) {
  #' Apply the board type sample to the forecasted population Dataframe
  #' 
  #' @param forecast_population The Dataframe with forecasted population information
  #' @param board_type_sample The Dataframe with the board type sampling probabilities for each CSD and Panel
  #' @return A Dataframe with population segmented by panel and board type
  #' 
  set.seed(42)
  forecast_population_by_board <- forecast_population %>%
    left_join(board_type_sample, by = c("panel", "cduid")) %>%
    filter(!is.na(cduid)) %>%
    # Calculate actual enrolment by multiplying potential.enrolment with the total participations rate
    mutate(enrolment = potential.enrolment * sample.total) %>% 
    # Create probability list for sample
    unite(prob, `English Catholic`, `English Public`, `French Catholic`, `French Public`, sep = ",") %>%
    rowwise() %>%
    mutate(prob = list(as.double(unlist(strsplit(prob, ","))))) %>%
    # Sample the different board types 
    mutate(sample.result = list(sample(c("EC", "EP", "FC", "FP"), size=enrolment, replace=TRUE, prob=prob))) %>%
    mutate(`English Catholic` = sum(sample.result == "EC"),
           `English Public` = sum(sample.result == "EP"),
           `French Catholic` = sum(sample.result == "FC"),
           `French Public` = sum(sample.result == "FP")) %>% 
    select(treso.id.por, panel, cduid, `English Catholic`:`French Public`) %>%
    gather(key="board_type_name", value="enrolment", `English Catholic`:`French Public`)
  
  return(forecast_population_by_board)
}

calculate_school_weight_forecasting <- function(trip_list, school_df, eqao_2017, new_school) {
  #' Calculate weight factor for distributing students between two or more schools within a single TRESO zone
  #' 
  #' @param trip_list
  #' @param school_df
  #' @param eqao_2017
  #' @param new_school
  #' @return
  
  # Include the new_school in the master school list
  if(!is.null(new_school)){
    school_df <- school_df %>% 
      bind_rows(new_school)
  }
  
  # Calculate school weighting for TRESO zones with multiple schools and export it for
  pos_school <- trip_list %>%
    group_by(treso.id.pos) %>%
    summarise(trips = sum(trips)) %>%
    right_join(select(school_df, otg, sfis, treso.id.pos, school.name, dsb.index), by = c("treso.id.pos") ) %>% 
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

distribute_students_to_schools <- function(school_df, new_school_df, por, pos_full, prop_matrix, 
                                           treso_travel_time, treso_distance,
                                           is_weight_continuous) {
  #' Creates a new POS vector with new schools and distribute the students to the schools
  #' 
  #' @param school_df
  #' @param new_school_df
  #' @param por
  #' @param pos_full
  #' @param prop_matrix
  #' @param treso_travel_time
  #' @param treso_distance
  #' @param is_weight_continous
  #' @return A Dataframe of schools with simulated ADE 
  
  pos <- create_pos_vector(school_df, new_school_df, pos_full)
  
  print(paste0("ADE is: ", sum(por), ". OTG is: ", sum(pos), "."))
  
  prop_matrix_balanced <- matrix_balancing_1d(prop_matrix, por, pos, axis=1, is_weight_continuous=is_weight_continuous)
  trip_list <- reshape2::melt(prop_matrix_balanced) %>%
    arrange(Var1, Var2)
  colnames(trip_list) <- c("treso.id.por", "treso.id.pos", "trips")
  
  # Distribute the ADE at each zone to the individual schools of the zone
  school_summary_df_list <- forecast_school_ade(prop_matrix, trip_list, school_df, eqao_2017, new_school_df,
                                                treso_travel_time, treso_distance)
  
  school_summary <- school_summary_df_list[[1]]
  pos_travel_time <- school_summary_df_list[[2]]
  pos_travel_distance <- school_summary_df_list[[3]]
  
  return(list(school_summary, prop_matrix_balanced, pos_travel_time, pos_travel_distance))
}

distribute_students_within_zones <- function(school_forecast_df) {
  #' Create the forecasted school list with OTG, OTG Threshold, ADE and Simulated ADE
  #' Also redistribute any overfilled schools in a zone to 
  #' 
  #' @param school_forecast_df
  #' @return A Dataframe of schools with ADE distributed among underfilled schools in the TRESO zone 
  
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

create_forecast_school_list <- function(school_base, school_20xx, new_school_df, school_diff, otg_threshold) {
  
  #' Create the forecasted school list with OTG, OTG Threshold, ADE and Simulated ADE
  #' Also redistribute any overfilled schools in a zone to 
  #' 
  #' @param school_base
  #' @param school_20xx
  #' @param new_school_df
  #' @param school_diff
  #' @param otg_threshold
  #' @return A Dataframe of schools with forecasted/simulated ADE and 
  
  school_forecast_df <- filter(school_base, is.closed == 0) %>% 
    # Join in the is.consolidated flag from school_20xx
    left_join(select(school_20xx, sfis, is.consolidated), by="sfis") %>% 
    # Add in new schools
    bind_rows(new_school_df) %>%
    left_join(select(school_diff, sfis, simulated.ade.20xx, simulated.ade.base, change.ade), by = "sfis") %>%
    replace_na(list(ade = 0, change.ade = 0, is.consolidated = 0)) %>%
    # 'Actual' forecast ADE = existing ADE (2017 actual historical data) plus change in ADE estimated in Step 3
    mutate(simulated.ade = ade + change.ade) %>% 
    # If simulated.ade.raw < 0 for any school, this value is rounded up to 0 to prevent having a negative number of students at each school.
    # However, by rounding up to 0, the model could be adding 'phantom' students to the system
    # So the total ADE as estimated by the distribution model will exceed total ADE estimated in population forecasts.
    # However, in test runs, this did not occur at any schools, so risk appears low. If it were to happen, the result would be a slight increase in the total number of ADE across the province, which would almost certainly be negligible. 
    mutate(simulated.ade = ifelse(simulated.ade < 0, 0, simulated.ade),
           simulated.ade = simulated.ade * (1 - is.consolidated)) %>%
    # Create OTG threshold based on user input
    mutate(otg.threshold = otg * otg_threshold) %>% 
    select(treso.id.pos, sfis, school.name, otg, otg.threshold, ade, simulated.ade)
  
  print('debug create_forecast_school_list')
  print(school_forecast_df %>% filter(sfis > 99000))
  
  # Redistribute students from overfilled schools to underfilled schools in the same zone
  school_forecast_distributed_df <- distribute_students_within_zones(school_forecast_df)
  
  return(school_forecast_distributed_df)
}

forecast_school_ade <- function(prop_matrix, trip_list, school_df, eqao_2017, new_school, 
                                treso_travel_time = NULL, treso_distance = NULL) {
  #' Produce a Dataframe with the summary of the school's forecasted ADE
  #'
  #' @param prop_matrix
  #' @param trip_list The trip list from the balanced matrix
  #' @param school_df
  #' @param eqao_2017
  #' @param new_school
  #' @param treso_travel_time
  #' @param treso_distance
  #' @return A Dataframe of schools with forecasted ADE
  #' 
  # Calculate the school weight
  pos_school_weight <- calculate_school_weight_forecasting(trip_list, school_df, eqao_2017, new_school = new_school)
  print(paste0("There are ",
               nrow(filter(pos_school_weight, length(school.weight.prob.list) == 1)),
               " TRESO zones with a single school, and ",
               nrow(filter(pos_school_weight, length(school.weight.prob.list) > 1)),
               " TRESO zones with multiple schools."))
  
  # Apply bucket rounding to chunks of data by TRESO POS
  results <- by(trip_list$trips, trip_list[c("treso.id.pos")], smart_round, simplify = TRUE)
  
  # Convert the output of `by()` to a Dataframe
  results2 <- sapply(results, I)
  colnames(results2) <- colnames(prop_matrix)
  rownames(results2) <- rownames(prop_matrix)
  
  df <- reshape2::melt(results2) %>%
    arrange(Var1, Var2)
  colnames(df) <- c('treso.id.por', 'treso.id.pos', 'enrolment.rounded')
  
  if (!is.null(treso_travel_time)) {
    # Get weighted mean travel time for each treso POS by cbind() the full trip list with travel time
    pos_travel_time <- bind_cols(df, select(treso_travel_time, value)) %>%
      group_by(treso.id.pos) %>%
      summarise(sim.avg.travel.time = weighted.mean(value, enrolment.rounded)) %>% 
      replace_na(list(sim.avg.travel.time = 0))
  } else {
    pos_travel_time <- NULL
  }
  
  if (!is.null(treso_distance)) {
    # Get wegithed mean travel distance for each TRESO POS by cbind() the full trip list with travel distnace
    pos_travel_distance <- bind_cols(df, select(treso_distance, value)) %>% 
      group_by(treso.id.pos) %>% 
      summarise(sim.avg.travel.distance = weighted.mean(value, enrolment.rounded)) %>% 
      replace_na(list(sim.avg.travel.distance = 0))
  } else {
    pos_travel_distance <- NULL
  }
  
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
  
  df_list <- list(schools_summary, pos_travel_time, pos_travel_distance)
  
  return(df_list)
}

edu_dm <- function(treso_travel_time, trip_list, treso_zone_def, por_additional,
                   travel_time_threshold_factor, zone_proximity_threshold, min_tt_threshold) {
  #' Produce a list of 2 Dataframes containing: 
  #' 1. Dataframe of TRESO zones in which to consider building schools
  #' 2. Dataframe of TRESO zones ruled out from building schools due to proximity to higher-ranked location for building school
  #'
  #' @param treso_travel_time A Dataframe containing travel times to and from all origin-destination TRESO pairs
  #' @param trip_list The trip list from the balanced matrix
  #' @param treso_zone_def Details of treso zones, e.g., CSDUID/CDUID
  #' @param por_additional A list of TRESO origins which have residual students due to school overfills --> candidate zones for building a new school
  #' @param travel_time_threshold_factor A user-set factor for selecting how much travel time is considered 'excess' time
  #' @param zone_proximity_threshold A user-set threshold to preclude construction of schools in each of two TRESO zones too close together
  #' @param min_tt_threshold A user-set threshold to preclude construction of schools in TRESO zones without sufficient 'excess' travel time
  #' @return A list of 2 Dataframes with TRESO zones in which to build, and TRESO zones ruled out due to proximity
  
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
  
  if (!is.null(por_additional)) {
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
  } else {
    additional_pairs <- threshold_tt_zones
    
    shortlist_zones <- threshold_tt_zones %>%
      # Cut down list of zones under consideration based on minimum bar for 'excess' travel time
      filter(excess.travel.time > min_tt_threshold) %>%
      select(treso.id.por) %>% 
      distinct(treso.id.por)
  }
  
  # Retain the original copy of the shortlisted zones
  shortlist_zones_dup <- shortlist_zones
  
  shortlist_pairs <- potential_zones %>%
    select(treso.id.por, treso.id.pos, travel.time) %>%
    inner_join(shortlist_zones, by = c('treso.id.por')) %>%
    inner_join(shortlist_zones, by = c('treso.id.pos' = 'treso.id.por'))
  
  # Initialize Dataframes and while-loop check
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

calculate_delta_between_simulated_scenarios <- function(school_base, school_20xx, school_summary_base, school_summary_forecast, 
                                                        new_school, user_otg_threshold) {
  
  #' Calculate the delta between the two simulated scenarios. Apply this delta to the `school_base` Dataframe in order to obtain the
  #' forecasted ADE of the schools.
  #'
  #' @param school_base A Dataframe of the base year list of schools
  #' @param school_20xx A Dataframe of the 20xx list of schools
  #' @param school_summary_base A Dataframe of simulated ADE for schools in the base scenario
  #' @param school_summary_forecast A Dataframe of simulated ADE for schools in the forecast scenario
  #' @param new_school A Dataframe of new schools
  #' @param user_otg_threshold A numeric value set by the user
  #' @return A Dataframe of schools with the forecasted ADE
  
  # Take the difference between the simulated ADEs
  trip_list_schools_summary_difference <- full_join(school_summary_base, school_summary_forecast, by="sfis", suffix=c(".base", ".20xx")) %>% 
    replace_na(list(simulated.ade.base = 0, simulated.ade.20xx = 0)) %>% 
    mutate(change.ade = simulated.ade.20xx - simulated.ade.base)
  
  # Create the forecasted school list with OTG, OTG Threshold, ADE and Simulated ADE
  school_forecast <- create_forecast_school_list(school_base, school_20xx, new_school, trip_list_schools_summary_difference, 
                                                 user_otg_threshold)
  
  return(school_forecast)
}

distribution_model <- function(school_base, school_20xx, school_summary_2017, school_summary_20xx, prop_matrix,
                               prop_matrix_balanced, new_school, pos_full, overfill_threshold = 100, max_iter = 20,
                               user_otg_threshold, capacity_constrained = FALSE) {
  
  # Calculate the delta between simulated base and forecast scenarios and apply to the base scenario
  school_forecast <- calculate_delta_between_simulated_scenarios(school_base, school_20xx, school_summary_2017, school_summary_20xx,
                                                                 new_school, user_otg_threshold)
  
  # Calculate the total overfilled ADE
  overfill_ade <- school_forecast %>% 
    mutate(overfill.flag = ifelse(simulated.ade >= otg.threshold, 1, 0),
           overfill.ade = simulated.ade - otg.threshold) %>% 
    filter(overfill.flag == 1) %>% 
    select(overfill.ade) %>% 
    colSums()
  
  if (capacity_constrained & overfill_ade[[1]] >= overfill_threshold) {
    print(paste0("Running distribution capacity constrained"))
    i <- 1
    pos <- create_pos_vector(filter(school_20xx, is.consolidated == 0), new_school, pos_full)
    
    # Loop through until the overfill ade is assigned
    while(overfill_ade[[1]] >= overfill_threshold) {
      print(paste0("Iteration: ", i, ". Overfill ADE: ", round(overfill_ade[[1]], 0)))
      if (i <= max_iter) {
        # Create a POS vector with schools over OTG Threshold after redistributing to other schools within the TRESO zone
        pos_overfill <- school_forecast %>%
          mutate(overfill.flag = ifelse(simulated.ade >= otg.threshold, 1, 0)) %>% 
          filter(overfill.flag == 1) %>% 
          select(treso.id.pos, overfill.flag) %>% 
          group_by(treso.id.pos) %>% 
          summarise(overfill.flag = max(overfill.flag)) %>% 
          right_join(pos_full, by=c("treso.id.pos" = "dest")) %>% 
          replace_na(list(overfill.flag = 0)) %>% 
          arrange(treso.id.pos) %>% 
          column_to_rownames(var = "treso.id.pos") %>% 
          data.matrix()
        
        # Calculate the total overfilled ADE
        overfill_ade <- school_forecast %>% 
          mutate(overfill.flag = ifelse(simulated.ade >= otg.threshold, 1, 0),
                 overfill.ade = simulated.ade - otg.threshold) %>% 
          filter(overfill.flag == 1) %>% 
          select(overfill.ade) %>% 
          colSums()
        
        # Create a POS vector with schools under OTG Threshold - pos_underfill
        # Use the difference between OTG Threshold and Simulated ADE as the weight
        pos_underfill <- school_forecast %>% 
          mutate(underfill.flag = ifelse(simulated.ade < otg.threshold, 1, 0)) %>% 
          filter(underfill.flag == 1) %>%
          mutate(otg.leftover = otg.threshold - simulated.ade) %>% 
          select(treso.id.pos, otg.leftover) %>% 
          group_by(treso.id.pos) %>% 
          summarise(otg = sum(otg.leftover)) %>% 
          right_join(pos_full, by=c("treso.id.pos" = "dest")) %>% 
          replace_na(list(otg = 0)) %>% 
          arrange(treso.id.pos) %>%
          column_to_rownames(var = "treso.id.pos") %>% 
          data.matrix()
        
        # Apply pos_overfill to balanced matrix to obtain origin zones of the overflown schools
        por_overfill <- prop_matrix_balanced %*% pos_overfill %>% 
          as_tibble(rownames = NA) %>% 
          rename(ade = overfill.flag) %>% 
          rownames_to_column(var = "treso.id.por") %>% 
          # Scale the zone origin ADE to match the overfill ADE total previously calculated (overfill_ade)
          mutate(ade.overfill = ade * overfill_ade[[1]] / sum(ade)) %>% 
          select(treso.id.por, ade.overfill)
        
        # Store a list of TRESO zones where the students are assigned to overfilled schools
        por_additional <- por_overfill %>% 
          filter(ade.overfill >= 15) %>% 
          mutate(treso.id.por = as.numeric(treso.id.por)) %>% 
          select(treso.id.por, ade.overfill)
        
        # Add in any origin zones sending students to school zones still overfilled after several rounds of balancing
        por_overfill <- por_overfill %>% 
          column_to_rownames(var = "treso.id.por") %>% 
          data.matrix()
        
        # Sanity check with the dot product operation
        print(paste0("The first element of por_overfill is: ", (prop_matrix_balanced %*% pos_overfill)[1,1]))
        print(paste0("Compare that with the sum of first column of prop_matrix_balanced: ", sum(prop_matrix_balanced[1,as.logical(pos_overfill)])))
        
        print(paste0("The scaled POR vector sums to: ", sum(por_overfill), ". Which should be the same as the overfill_ade: ", overfill_ade[[1]]))
        
        # Balance the por_overfill against pos_underfill schools to redistribute the overfilled ADE
        print(paste0("ADE overfill total is: ", sum(por_overfill), ". OTG remaining is: ", sum(pos_underfill), "."))
        
        # Use total OTG to balance with the POR_OVERFILL 
        prop_matrix_balanced <- matrix_balancing_1d(prop_matrix, por_overfill, pos, axis=1, is_weight_continuous=TRUE)
        trip_list <- reshape2::melt(prop_matrix_balanced) %>%
          arrange(Var1, Var2)
        colnames(trip_list) <- c("treso.id.por", "treso.id.pos", "trips")
        
        # Produce the simulated ADE for each school
        school_summary_20xx_iterated <- forecast_school_ade(prop_matrix, trip_list, filter(school_20xx, is.consolidated == 0), 
                                                            eqao_2017, new_school, treso_travel_time = NULL, treso_distance = NULL)[[1]]
        
        # Calculate the ADE overfill for each school, so the overfill ADE can be removed in the summary
        school_summary_ade_overfill <- school_forecast %>% 
          mutate(overfill.flag = ifelse(simulated.ade >= otg.threshold, 1, 0),
                 overfill.ade = simulated.ade - otg.threshold) %>% 
          filter(overfill.flag == 1) %>% 
          select(sfis, overfill.ade)
        
        school_summary_20xx <- full_join(school_summary_20xx, school_summary_20xx_iterated,
                                         by="sfis", suffix=c("", ".iteration")) %>% 
          replace_na(list(simulated.ade = 0, simulated.ade.iteration = 0)) %>% 
          mutate(simulated.ade = simulated.ade + simulated.ade.iteration) %>% 
          select(-simulated.ade.iteration) %>% 
          left_join(school_summary_ade_overfill, by="sfis") %>% 
          replace_na(list(overfill.ade = 0)) %>% 
          mutate(simulated.ade = simulated.ade - overfill.ade) %>% 
          select(-overfill.ade)
        
        # Calculate the delta between simulated base and forecast scenarios and apply to the base scenario
        school_forecast <- calculate_delta_between_simulated_scenarios(school_base, school_20xx, school_summary_2017, 
                                                                       school_summary_20xx, new_school, user_otg_threshold)
        
        # Increment the counter
        i <- i + 1
      } else {
        print("Reached Max Iterations")
        return(list(school_forecast, por_additional))
      }
    }
    return(list(school_forecast, por_additional))
    
  } else {
    print(paste0("Running distribution capacity unconstrained"))
    # Calculate the delta between simulated base and forecast scenarios and apply to the base scenario
    school_forecast <- calculate_delta_between_simulated_scenarios(school_base, school_20xx, school_summary_2017, 
                                                                   school_summary_20xx, new_school, user_otg_threshold)
    return(list(school_forecast, NULL))
  }
}

# MOH ----
## MOH Cleaning ----

create_hospital_xy <- function(hospital_master) {
  #' This function transforms each hospital location into the same projection as the TRESO shapefiles
  #' 
  #' @param hospital_master A Dataframe of the hospitals
  #' @param treso_shp A Shapefile of TRESO zones
  #' @return A SpatialPointsDataframe 
  
  hospital_spdf <- hospital_master %>%
    ungroup() %>%
    select(id, hospital.lat, hospital.long) %>%
    mutate(lat = hospital.lat, long = hospital.long) %>%
    mutate(ref = row_number())
  coordinates(hospital_spdf) <- c('long', 'lat')
  
  # Project the hospital lat/long to TRESO's LCC specification
  proj4string(hospital_spdf) <- CRS('+proj=longlat +datum=WGS84')
  
  hospital_xy <- spTransform(hospital_spdf, CRS(treso_projarg))
  
  return(hospital_xy)
}

## MOH Model ----
base_scenario_alc_calculation <- function(hospital_lookup_ALC, ALC_hosp_vol_disag, ltc_homes, treso_population_moh_agecluster, ltc_op_days, LTC_FLAG) {
  # Calculate ALC days by hospital by caretype for 2016; this will be used in average length of stay calculations to distribute ALC based on caretype
  ALC_hosp_vol <- ALC_hosp_vol_disag %>% 
    left_join(hospital_lookup_ALC, by = c('siteid' = 'ALCsite_id')) %>%
    filter(!is.na(csduid)) %>% # Double-check in case any ALC sites don't match to Hospital PAI or master list
    group_by(year, id, agecluster) %>% 
    mutate(hosp.alc.patients = sum(patients),
           hosp.alc.days = sum(adj.hosp.days.tot),
           avg.los = hosp.alc.days / hosp.alc.patients) %>%
    group_by(year, csduid) %>% 
    mutate(hosp.alc.mrkt = hosp.alc.days / sum(hosp.alc.days)) %>% #MARKET-SHARE OF HOSPITAL ALC DAYS TO CSD TOTAL DAYS
    ungroup() %>% 
    distinct(year, id, name = name.x, csduid, csdname, agecluster, hosp.alc.patients, hosp.alc.days, avg.los, hosp.alc.mrkt)
  
  # SUMMARIZE BED INFO BY CSD
  # This assumes that LTC facilities are full, ie., all beds are taken, i.e., capacity (patients) = demand
  ltc_homes_csd <- ltc_homes %>% 
    group_by(csduid, csdname) %>% 
    summarise(csd.ltc.patients = sum(beds) * LTC_FLAG) %>% 
    mutate(ltc.bed.days.cap = csd.ltc.patients * ltc_op_days,
           csdname = toupper(csdname)) %>% 
    ungroup()
  
  # Calculate the percentage of ALC days used by agecluster == adult vs. agecluster == senior; JOINING LTC CAPACITY DATA FROM LTC-HOMES DATASET
  ALC_csd_vol_join <- ALC_hosp_vol %>%
    group_by(year, csduid, agecluster, csdname) %>% 
    summarise_at(vars(hosp.alc.patients:hosp.alc.days), sum) %>%
    rename(csd.alc.patients = hosp.alc.patients,
           csd.alc.days = hosp.alc.days) %>%
    group_by(year, csduid) %>% 
    mutate(csd.alc.days.sum.ages = sum(csd.alc.days), csd.alc.patients.sum.ages = sum(csd.alc.patients)) %>% 
    mutate(csd.alc.days.percent.ages = csd.alc.days / csd.alc.days.sum.ages, csd.alc.patients.percent.ages = csd.alc.patients / csd.alc.patients.sum.ages) %>% 
    mutate(csdname = toupper(csdname)) %>% 
    # Full join to ltc_homes_csd since there may be CSDs with hospitals but no LTC facilities, and vice versa
    full_join(ltc_homes_csd, by = c('csduid')) %>%
    # Split existing LTC patients along hospital observed proportions in terms of ageclusters, then sum total by CSD
    # Replace NA values for: CSDs w/ LTC facilities but no hospital ALC days; CSDs w/ LTC but no hospital ALC patients; CSDs w/ no hospital ALC; CSDs w/ no hosp. ALC
    replace_na(list(csd.ltc.patients = 0, ltc.bed.days.cap = 0, csd.alc.days = 0, csd.alc.patients = 0, csd.alc.days.percent.ages = 0,
                    csd.alc.patients.percent.ages = 0)) %>% 
    ungroup() %>% 
    filter(is.na(csduid) == FALSE) %>%
    mutate(csdname = ifelse(is.na(csdname.x), csdname.y, csdname.x)) %>% 
    select(-csdname.x, -csdname.y)
  
  # Calculate average % breakdown for agecluster = senior vs. agecluster = adult for CSDs where there ARE hospitals with ALC patients
  ALC_LTC_age_breakdown <- ALC_csd_vol_join %>%
    filter(is.na(agecluster) == FALSE) %>% 
    group_by(year, agecluster) %>%  
    mutate(sum.alc.days = sum(csd.alc.days), sum.alc.patients = sum(csd.alc.patients)) %>% 
    group_by(year) %>% 
    mutate(prov.alc.days.percent = sum.alc.days / sum(csd.alc.days), 
           prov.alc.patients.percent = sum.alc.patients / sum(csd.alc.patients)) %>% 
    select(year, agecluster, prov.alc.days.percent, prov.alc.patients.percent) %>% 
    distinct()
  
  # CSDs with LTC capacity but not hospital ALC capacity do not have an inherent split into age clusters; simulate this split by creating a dataframe with CSDs, LTC capacity, and ageclusters for use in next query
  ALC_LTC_csd_ageclusters <- ALC_csd_vol_join %>% 
    filter(is.na(agecluster))
  
  agelist <- c('ad', 'sr')
  yearlist <- c(2013, 2014, 2015, 2016) # List of historical years during which LTC beds in CSDs with no hospital ALC 
  
  ALC_LTC_csd_age_cross <- crossing(ALC_LTC_csd_ageclusters, agelist, yearlist) %>% 
    select(csduid, agecluster = agelist, year = yearlist)
  
  # Clean up ALC_csd_vol_join and calculate ALC demand at the CSD level by summing ALC beds/bed-days plus LTC beds/bed-days for each age cluster
  ALC_csd_vol <- ALC_csd_vol_join %>% 
    left_join(ALC_LTC_csd_age_cross, by = c('csduid')) %>% 
    mutate(agecluster.x = ifelse(is.na(agecluster.x), agecluster.y, agecluster.x)) %>%
    rename(agecluster = agecluster.x) %>% 
    select(-agecluster.y) %>% 
    mutate(year.x = ifelse(is.na(year.x), year.y, year.x)) %>%
    rename(year = year.x) %>% 
    select(-year.y) %>% 
    # Join in average % breakdown for agecluster = senior vs. agecluster = adult for use in CSDs where there are no hospitals with ALC patients and therefore cannot calculate a CSD-specific breakdown
    left_join(ALC_LTC_age_breakdown, by = c('year', 'agecluster')) %>% 
    # Calculate total ALC demand (ALC days observed plus LTC days)
    group_by(year, csduid) %>% 
    mutate(csd.alc.days.percent.ages = ifelse(csd.alc.days.percent.ages > 0, csd.alc.days.percent.ages, prov.alc.days.percent),
           csd.alc.patients.percent.ages = ifelse(csd.alc.patients.percent.ages > 0, csd.alc.patients.percent.ages, prov.alc.patients.percent)) %>% 
    mutate(csd.demand.days.ages = round(ltc.bed.days.cap * csd.alc.days.percent.ages + csd.alc.days, 0),
           csd.demand.patients.ages = round(csd.ltc.patients * csd.alc.patients.percent.ages + csd.alc.patients, 0),
           csd.demand.days.total = sum(csd.demand.days.ages),
           csd.demand.patients.total = sum(csd.demand.patients.ages)) %>% 
    ungroup() 
  
  # CSD/CD rates not particularly meaningful given the huge variation from one CSD/CD to the next; use provincial average values, disaggregated by agecluster, for ALC rates per population.
  treso_population_long <- gather(treso_population_moh_agecluster, 'year', 'population', population.2041:population.2040, -age.group, -agecluster) %>% 
    separate(year, into = c('pop', 'year'), sep = "n.") %>% 
    mutate(year = as.numeric(year)) %>%
    select(-pop) %>% 
    group_by(year, treso_zone, agecluster) %>% 
    mutate(treso.population.agecluster = sum(population))
  
  ALC_csd_rates <- treso_population_long %>%
    filter(agecluster != 'ped', year > 2012, year < 2017) %>% 
    mutate(cduid = substr(csduid, 1, 4)) %>% 
    group_by(year, cduid, csduid, agecluster) %>% 
    summarise(csd.population.ages = sum(population)) %>% 
    left_join(select(ALC_csd_vol, year, csduid, agecluster, csd.demand.patients.ages, csd.demand.days.ages), by = c('year', "csduid", 'agecluster')) %>% 
    mutate(csd.demand.patients.ages = replace_na(csd.demand.patients.ages, 0),
           csd.demand.days.ages = replace_na(csd.demand.days.ages, 0),
           csd.alc.rate.ages = csd.demand.days.ages / csd.population.ages) %>% 
    group_by(year, cduid, agecluster) %>% 
    mutate(cd.demand.days.ages = sum(csd.demand.days.ages),
           cd.alc.rate.ages = cd.demand.days.ages / sum(csd.population.ages)) %>% 
    ungroup()
  
  # Calculate provincial rates averaging across all CDs, all years, since huge variability observed in terms of ALC days per population by CSD/CD. E.g., for adults, some CDs have as little as 0.01 ALC days per person, whereas some have 2.6 days per person. For seniors, some CDs have as little as 0.5 days per person, whereas others have as many as 56 days per person (when ignoring CDs with no facilities and therefore a rate of 0). Clearly, because hospitals and ALC facilities are not evenly distributed throughout all CDs, 'demand' by CD swings wildly based on availability (i.e., capacity), with out-of-CD users certainly skewing the ALC rates when compared to a CD's population. 
  ALC_prov_rate <- ALC_csd_rates %>% 
    filter(year %in% (2016:2016)) %>% # Adjust this filter to choose specific years to include in base rate calculations
    group_by(agecluster) %>%  # Not grouped by 'year' so populations in calculations below appear high, but provincial ALC rate calculations are accurate, they simply represent an average of historical years
    summarise(total.pop.ages = sum(csd.population.ages),
              total.dem.days.ages = sum(csd.demand.days.ages)) %>% 
    mutate(prov.alc.rate.ages = total.dem.days.ages / total.pop.ages)
  
  return(ALC_prov_rate)
}

base_scenario_crm_calculation <- function(HBAMhistorical_master_tz, treso_population_moh_agecluster) {
  # Calculate total cases by hospital, and the percentage of cases belonging to each CSD for each hospital siteid
  HBAMhistorical_origcsd_cases <- HBAMhistorical_master_tz %>%
    select(-agegroup) %>% 
    group_by(year, id, agecluster, caretype, orig_csd, treso_zone, caretype, hosp_csduid) %>%
    summarise(hosp.csd.agecluster.caretype.cases = sum(cases)) %>%
    ungroup() %>% 
    group_by(year, id, agecluster, caretype) %>% 
    mutate(hosp.agecluster.caretype.cases = sum(hosp.csd.agecluster.caretype.cases),
           csd.agecluster.caretype.mrkt = hosp.csd.agecluster.caretype.cases / hosp.agecluster.caretype.cases) %>%
    rename(hosp.treso.zone = treso_zone) %>% 
    ungroup()
  
  # CSD/CD rates not particularly meaningful given the huge variation from one CSD/CD to the next; use provincial average values, disaggregated by agecluster, for ALC rates per population.
  treso_population_long <- gather(treso_population_moh_agecluster, 'year', 'population', population.2041:population.2040, -age.group, -agecluster) %>% 
    separate(year, into = c('pop', 'year'), sep = "n.") %>% 
    mutate(year = as.numeric(year)) %>%
    select(-pop) %>% 
    group_by(year, treso_zone, agecluster) %>% 
    mutate(treso.population.agecluster = sum(population))
  
  # Join historical HBAM cases by CSD/hospital pair to TRESO information to get TRESO populations
  # Then calculate case rates for TRESO zones at the CSD-Hospital pair, segmented by caretype, year, and agecluster
  treso_origin_cases <- HBAMhistorical_origcsd_cases %>% 
    filter(year == 2016) %>% # select most recent historical year
    left_join(distinct(treso_population_long, treso_zone, agecluster, csduid, year, treso.population.agecluster), by = c('year', 'orig_csd' = 'csduid', 'agecluster')) %>% 
    rename(orig.treso.zone = treso_zone,
           orig.treso.pop.agecluster = treso.population.agecluster) %>%
    # A few CSDs don't exist in TRESO, but small case volumes (~2,700 rows out of 2.64M, or 73,000 cases out of 8.6M) so removed for simplicity
    drop_na(orig.treso.zone) %>% 
    # Divide market share from CSDs into TRESO zones
    group_by(id, caretype, agecluster, orig_csd) %>% 
    mutate(orig.csd.pop.agecluster = sum(orig.treso.pop.agecluster)) %>% 
    mutate(treso.caserate = hosp.csd.agecluster.caretype.cases / orig.csd.pop.agecluster) %>% 
    # Check to see how many cases remain after dropping N/A TRESO zones
    mutate(testCases = treso.caserate * orig.treso.pop.agecluster) %>% 
    mutate(hosp.CSD.age.care.cases = sum(testCases)) %>% 
    ungroup() %>% 
    # Select key fields to bring forwards
    select(id,
           hosp.treso.zone,
           hosp_csduid,
           caretype,
           agecluster,
           orig_csd,
           orig.treso.zone,
           treso.caserate)
  
  return(treso_origin_cases)
}

forecast_scenario_alc_calculation <- function(treso_population_moh, ltc_turnover_raw, new_ltc_cap_csd, hospitallocations_tz, 
                                              ALC_hosp_vol_disag, ALC_prov_rate, hospital_lookup_ALC, ltc_homes, ltc_op_days, scenario_year, 
                                              LTC_FLAG, USER_LTC_RATE_FACTOR) {
  
  # Converts annual turnover into length of stay (days) WITHIN a given calendar year --> this does not imply the average length of stay overall is less than a year. However, WITHIN a given year, the length of stay tops out at 365 days.
  ltc_los_per_year <- ltc_turnover_raw %>% 
    mutate(avg.annual.los = 1 / (1 + system.turnover.2018) * LTC_OP_DAYS) %>% 
    left_join(distinct(hospitallocations_tz, lhin, csduid), by = c('lhin')) %>% 
    # This function fails for Toronto's CSD (3520005) which is associated with several LHINs; assume associated with LHIN 7 for purposes of length of stay evaluation
    mutate(csd.lhin.concat = paste(lhin, csduid, sep = '_')) %>% 
    filter(csd.lhin.concat != '5_3520005', csd.lhin.concat != '6_3520005', csd.lhin.concat != '8_3520005', csd.lhin.concat != '9_3520005', csd.lhin.concat != '10_3520005') %>% 
    select(-csd.lhin.concat)
  
  ALC_hosp_vol <- ALC_hosp_vol_disag %>% 
    left_join(hospital_lookup_ALC, by = c('siteid' = 'ALCsite_id')) %>%
    filter(!is.na(csduid)) %>% # Double-check in case any ALC sites don't match to Hospital PAI or master list
    group_by(year, id, agecluster) %>% 
    mutate(hosp.alc.patients = sum(patients),
           hosp.alc.days = sum(adj.hosp.days.tot),
           avg.los = hosp.alc.days / hosp.alc.patients) %>%
    group_by(year, csduid) %>% 
    mutate(hosp.alc.mrkt = hosp.alc.days / sum(hosp.alc.days)) %>% #MARKET-SHARE OF HOSPITAL ALC DAYS TO CSD TOTAL DAYS
    ungroup() %>% 
    distinct(year, id, name = name.x, csduid, csdname, agecluster, hosp.alc.patients, hosp.alc.days, avg.los, hosp.alc.mrkt)
  
  # This query pulls actual historical ALC days for 2016
  ALC_proj_hosp_2016 <- ALC_hosp_vol %>% 
    group_by(year, id, name, agecluster, csduid, csdname) %>% 
    filter(year == 2016) %>% 
    summarise_at(vars(hosp.alc.patients:hosp.alc.days), sum) %>% 
    ungroup() %>% 
    mutate(total.LTCalc.days = sum(hosp.alc.days))
  
  # Calculate capacity of LTC by CSD, including any additional beds added by visualization tool user
  ltc_proj_homes <- ltc_homes %>% 
    group_by(csduid) %>%
    mutate(existing.csd.beds = sum(beds)) %>% 
    left_join(select(new_ltc_cap_csd, csduid, added.ltc.cap), by = c('csduid')) %>% 
    mutate(added.ltc.cap = replace_na(added.ltc.cap, 0)) %>% 
    mutate(total.csd.beds = (existing.csd.beds + added.ltc.cap) * LTC_FLAG,
           ltc.bed.days.cap = total.csd.beds * ltc_op_days) %>% 
    ungroup() %>% 
    mutate(csdname = toupper(csdname)) %>%
    mutate(sumBeds = sum(beds))
  
  # Calculate population by TRESO zone / CSD for projection year
  proj_pop_control <- treso_population_moh %>%
    filter(age.group != 1) %>% # Filter out ages 0-17 as this age group is not represented in ALC data
    mutate(agecluster = ifelse(age.group < 5, 'ad', 'sr')) %>% 
    group_by(treso_zone, agecluster, csduid) %>%
    summarise_at(vars(paste0("population.", scenario_year)), sum)
  
  # Calculate ALC demand as a function of population growth, join to capacity; because of the join to capacity, agecluster is phased out (since capacity is not broken down by agecluster)
  ALC_proj_csd <- proj_pop_control %>%
    left_join(select(ALC_prov_rate, agecluster, prov.alc.rate.ages), by = c('agecluster')) %>% 
    mutate(ltc.treso.demand.days = prov.alc.rate.ages * get(paste0("population.", scenario_year))) %>% 
    # Adjust LTC rate using user-input factor; when this function is run for non-LTC ALC, ignore this factor change
    mutate(ltc.treso.demand.days = ltc.treso.demand.days + ltc.treso.demand.days * LTC_FLAG * (USER_LTC_RATE_FACTOR - 1)) %>% 
    group_by(csduid) %>%
    mutate(ltc.csd.demand.days = sum(ltc.treso.demand.days)) %>%
    # Join to capacity data
    left_join(distinct(ltc_proj_homes, csduid, total.csd.beds, ltc.bed.days.cap), by = 'csduid') %>% 
    # Join to length-of-stay data
    left_join(select(ltc_los_per_year, avg.annual.los, csduid), by = c('csduid')) %>%
    group_by(csduid) %>% 
    mutate(total.csd.beds = replace_na(total.csd.beds, 0), 
           ltc.bed.days.cap = replace_na(ltc.bed.days.cap, 0),
           csd.residual.alc.days = ltc.csd.demand.days - ltc.bed.days.cap,
           csd.patients.demand = ltc.csd.demand.days / avg.annual.los,
           csd.patients.cap = ltc.bed.days.cap / avg.annual.los,
           csd.residual.patients = csd.patients.demand - csd.patients.cap) %>%
    mutate(cduid = substr(csduid, 1, 4)) %>% 
    ungroup() %>% 
    distinct(cduid, csduid, total.csd.beds, ltc.bed.days.cap, avg.annual.los, csd.patients.demand, csd.patients.cap, ltc.csd.demand.days, csd.residual.alc.days, csd.residual.patients) %>%
    # Use the code below to sum to the CD level if desired
    # group_by(cduid) %>% 
    # mutate(total.cd.beds = sum(total.csd.beds), ltc.bed.days.cap.cd = sum(ltc.bed.days.cap), cd.patients.demand = sum(csd.patients.demand), cd.patients.cap = sum(csd.patients.cap), ltc.cd.demand.days = sum(ltc.csd.demand.days)) %>% 
    # mutate(cd.residual.patients = cd.patients.demand - cd.patients.cap,  
    #        cd.residual.alc.days = ltc.cd.demand.days - ltc.bed.days.cap.cd) %>% 
    # distinct(total.cd.beds, ltc.bed.days.cap.cd, cd.patients.demand, cd.patients.cap, ltc.cd.demand.days, cd.residual.patients, cd.residual.alc.days) %>% 
    # ungroup() %>%
    # Set CD Length of Stay to be equal to CSD length of stay since length of stay is calculated at the LHIN level. This is used so that CDs/CSDs with no LTC capacity can have patient demand estimated as a function of their bed-days demand
    group_by(cduid) %>%
    mutate(avg.annual.los = replace_na(avg.annual.los, 0)) %>% 
    mutate(avg.annual.los = max(avg.annual.los)) %>% 
    mutate(csd.patients.cap = replace_na(csd.patients.cap, 0),
           csd.patients.demand = replace_na(csd.patients.demand, 0)) %>%
    mutate(csd.patients.cap = ifelse(csd.patients.cap == 0, ltc.bed.days.cap / avg.annual.los, csd.patients.cap),
           csd.patients.demand = ifelse(csd.patients.demand == 0, ltc.csd.demand.days / avg.annual.los, csd.patients.demand),
           csd.residual.patients = ifelse(is.na(csd.residual.patients), csd.patients.demand - csd.patients.cap, csd.residual.patients)) %>% 
    ungroup() %>% 
    mutate(prov.demand.days = sum(ltc.csd.demand.days), prov.cap.days = sum(ltc.bed.days.cap), 
           prov.residual.days = sum(csd.residual.alc.days))
  
  # Roll up ALC_proj_csd to the CD level since 'residual' ALC demand is hard to reconcile with capacity at the CSD level
  ALC_proj_cd <- ALC_proj_csd %>% 
    group_by(cduid) %>% 
    summarise(ltc.cd.demand.days = sum(ltc.csd.demand.days), 
              ltc.cd.bed.days.cap = sum(ltc.bed.days.cap), 
              cd.patients.demand = sum(csd.patients.demand),
              cd.patients.cap = sum(csd.patients.cap),
              cd.residual.alc.days = sum(csd.residual.alc.days),
              cd.residual.patients = sum(csd.residual.patients),
              prov.demand.days = first(prov.demand.days),
              prov.cap.days = first(prov.cap.days),
              prov.residual.days = first(prov.residual.days)) %>%
    # Calculate 'adjusted' residual ALC by CD if any excess LTC capacity is shared with other CDs
    #mutate(cd.extra.ltc.days = ifelse(cd.residual.alc.days < 0, -cd.residual.alc.days, 0)) %>% 
    ungroup() %>%
    mutate(sum.cd.extra.ltc.days = sum(cd.residual.alc.days[cd.residual.alc.days < 0]),
           sum.cd.shortfall.ltc.days = sum(cd.residual.alc.days[cd.residual.alc.days >= 0]),
           percent.cd.shortfall.ltc.days = ifelse(cd.residual.alc.days >= 0, cd.residual.alc.days / sum.cd.shortfall.ltc.days, 0),
           borrowed.cd.ltc.days = percent.cd.shortfall.ltc.days * sum.cd.extra.ltc.days * -1,
           adj.cd.residual.alc.days = ifelse(cd.residual.alc.days < 0, 0, pmax(cd.residual.alc.days - borrowed.cd.ltc.days,0)),
           sumBorrowed = sum(borrowed.cd.ltc.days), sumResidual = sum(adj.cd.residual.alc.days))
  
  # Similar to ALC_proj_cd but maintains agecluster disaggregation - used to reduce demand days in scenario forecast years
  ALC_proj_cd_agecluster <- proj_pop_control %>%
    left_join(select(ALC_prov_rate, agecluster, prov.alc.rate.ages), by = c('agecluster')) %>% 
    mutate(ltc.treso.demand.days = prov.alc.rate.ages * get(paste0("population.", scenario_year))) %>%
    # Adjust LTC rate using user-input factor; when this function is run for non-LTC ALC, ignore this factor change
    mutate(ltc.treso.demand.days = ltc.treso.demand.days + ltc.treso.demand.days * LTC_FLAG * (USER_LTC_RATE_FACTOR - 1)) %>% 
    mutate(cduid = substr(csduid, 1, 4)) %>% 
    group_by(cduid, agecluster) %>%
    mutate(ltc.cd.demand.days = sum(ltc.treso.demand.days)) %>% 
    # Determine what percentage of csd ltc demand days are applied to each age cluster
    select(agecluster, cduid, ltc.cd.demand.days) %>%
    distinct() %>% 
    group_by(cduid) %>% 
    mutate(ltc.cd.demand.days.agepercent = ltc.cd.demand.days / sum(ltc.cd.demand.days)) %>% 
    ungroup() %>% 
    select(-ltc.cd.demand.days)
  
  return(list(ALC_proj_cd_agecluster, ALC_proj_cd))
}

forecast_scenario_hbam_calculation <- function(HBAMprojected_master_tz, num_beds, year_id) {
  HBAM_control <- HBAMprojected_master_tz %>%
    mutate(los = replace_na(los, 0),
           beddays = replace_na(beddays, 0),
           beds.forecasted = beddays / HOSP_OP_DAYS) %>%
    filter(year == year_id) %>%
    rename(demand.days.total = beddays) %>% 
    group_by(id, caretype, agecluster) %>%
    summarise_at(vars(cases:beds.forecasted), sum, na.rm=TRUE) %>%
    group_by(id) %>% 
    mutate(sum.cases = sum(cases),
           sum.demand.days.total = sum(demand.days.total),
           sum.beds.forecasted = sum(beds.forecasted)) %>%
    ungroup() %>% 
    # Pivot table from long to wide for agecluster for demand columns (cases, demand-days)
    pivot_wider(names_from = c(caretype, agecluster), values_from = c(cases, demand.days.total, beds.forecasted), 
                values_fill = c(cases = 0, demand.days.total = 0, beds.forecasted = 0)) %>% 
    mutate(sum.demand.days.total_AM = demand.days.total_AM_ad + demand.days.total_AM_ped + demand.days.total_AM_sr,
           sum.demand.days.total_AT = demand.days.total_AT_ad + demand.days.total_AT_ped + demand.days.total_AT_sr,
           sum.demand.days.total_CR = demand.days.total_CR_ad + demand.days.total_CR_ped + demand.days.total_CR_sr,
           sum.demand.days.total_GR = demand.days.total_GR_ad + demand.days.total_GR_ped + demand.days.total_GR_sr,
           sum.demand.days.total_MH = demand.days.total_MH_ad + demand.days.total_MH_ped + demand.days.total_MH_sr,
           sum.demand.days = sum.demand.days.total_AM + sum.demand.days.total_AT + sum.demand.days.total_CR + 
             sum.demand.days.total_GR + sum.demand.days.total_MH)
  
  # Convert num_beds to be in wide form (i.e., one row per asset)
  num_beds_wide <- num_beds %>%
    rename(beds.existing = staffedbeds) %>% 
    mutate(caretype.beds = paste0('beds.existing.', caretype)) %>%
    select(-year, -caretype, -count) %>% 
    pivot_wider(names_from = c(caretype.beds), values_from = c(beds.existing), values_fill = c(beds.existing = 0)) %>% 
    rename(beds.existing = total.hosp.beds)
  
  hospital_wide_final_hbam <- HBAM_control %>% 
    left_join(num_beds_wide, by = c('id')) %>% 
    mutate(beds.forecasted.AT = beds.forecasted_AT_ad + beds.forecasted_AT_ped + beds.forecasted_AT_sr,
           beds.forecasted.CR = beds.forecasted_CR_ad + beds.forecasted_CR_ped + beds.forecasted_CR_sr,
           beds.forecasted.GR = beds.forecasted_GR_ad + beds.forecasted_GR_ped + beds.forecasted_GR_sr,
           beds.forecasted.MH = beds.forecasted_MH_ad + beds.forecasted_MH_ped + beds.forecasted_MH_sr,
           beds.forecasted = beds.forecasted.AT + beds.forecasted.CR + beds.forecasted.GR + beds.forecasted.MH,
           utilization.AM = 0,
           utilization.AT = sum.demand.days.total_AT / beds.forecasted.AT / HOSP_OP_DAYS,
           utilization.CR = sum.demand.days.total_CR / beds.forecasted.CR / HOSP_OP_DAYS,
           utilization.GR = sum.demand.days.total_GR / beds.forecasted.GR / HOSP_OP_DAYS,
           utilization.MH = sum.demand.days.total_MH / beds.forecasted.MH / HOSP_OP_DAYS,
           utilization = sum.demand.days / beds.forecasted / HOSP_OP_DAYS,
           utilization.AT.existing = sum.demand.days.total_AM / beds.existing.AT / HOSP_OP_DAYS,
           utilization.CR.existing = sum.demand.days.total_CR / beds.existing.CR / HOSP_OP_DAYS,
           utilization.GR.existing = sum.demand.days.total_GR / beds.existing.GR / HOSP_OP_DAYS,
           utilization.MH.existing = sum.demand.days.total_MH / beds.existing.MH / HOSP_OP_DAYS,
           utilization.existing = sum.demand.days / beds.existing / HOSP_OP_DAYS) %>% 
    # Replace 'Inf' utilizations into NA
    mutate_at(vars(utilization.existing, utilization.AT.existing, utilization.CR.existing, utilization.GR.existing,
                   utilization.MH.existing), ~replace(., is.infinite(.), NA))
  
  return(hospital_wide_final_hbam)
}


forecast_scenario_crm_calculation <- function(treso_population_moh_agecluster, treso_origin_cases, hospital_lookup_master, 
                                              HBAMprojected_alos_adj, crm_treso, LTC_proj_cd_agecluster, LTC_proj_cd, 
                                              ALC_proj_cd_agecluster_ex_LTC, ALC_proj_cd_ex_LTC, caretype_list, 
                                              age_cluseter_list, hosp_op_days, scenario_year) {
  # Calculate TRESO populations at the agecluster level instead of agegroup
  colname <- noquote(paste0('population.', scenario_year))
  treso_population_moh_agecluster_scenario <- treso_population_moh_agecluster %>% 
    select(-age.group) %>% 
    group_by(treso_zone, agecluster, csduid) %>% 
    summarise_all(sum) %>% 
    select(csduid, treso_zone, agecluster, !!as.symbol(colname)) %>% 
    rename(scen.proj.pop.treso = !!as.symbol(colname)) %>% 
    ungroup()
  
  # CRM Output table calculates the number of cases, bed-days and beds-needed (bed-days/operating-days) for each hospital-treso pair and then summarizes it to each hospital, calculating the totals for each unique hospital-segment for the given projection year
  # treso_origin_cases <- readRDS('cache/moh/treso_origin_cases.rds') # Use this query if reading in from stored data instead of keeping in RAM
  crm_treso <- treso_origin_cases %>%
    left_join(treso_population_moh_agecluster_scenario, by = c("orig.treso.zone" = "treso_zone", "agecluster")) %>%
    select(-csduid) %>% 
    left_join(select(HBAMprojected_alos_adj, id, agecluster, caretype, adj.caretype.agecluster.los), by = c('id', "agecluster", "caretype")) %>% 
    # Filter caretypes/ageclusters down to those selected by the user
    filter(caretype %in% caretype_list, agecluster %in% age_cluseter_list) %>% 
    # Calculate the number of cases, beddays, bedsneeded depending on the input year
    mutate(hosp.treso.cases = treso.caserate * scen.proj.pop.treso,
           hosp.treso.beddays = hosp.treso.cases * adj.caretype.agecluster.los,
           hosp.treso.bedsneeded = hosp.treso.beddays / hosp_op_days) %>%
    # Assign 0 bed days and 0 beds needed for Emergency (AM)
    mutate(adj.caretype.agecluster.los = replace_na(adj.caretype.agecluster.los, 0),
           hosp.treso.beddays = replace_na(hosp.treso.beddays, 0),
           hosp.treso.bedsneeded = replace_na(hosp.treso.bedsneeded,0)) %>% 
    # Remove any TRESO/hospital pairs which do not historically have cases
    filter(treso.caserate > 0)
  
  # Calculate breakdown by caretype for hospitals for use in next step (since ALC affects 'unknown' caretypes)
  crm_treso_caretypes <- crm_treso %>%
    group_by(id, agecluster, caretype) %>%
    summarise(hosp.treso.caretype.beddays = sum(hosp.treso.beddays),
              hosp.treso.caretype.cases = sum(hosp.treso.cases)) %>%
    group_by(id, agecluster) %>% 
    mutate(hosp.treso.caretype.beddays.percent = hosp.treso.caretype.beddays / sum(hosp.treso.caretype.beddays),
           hosp.treso.caretype.cases.percent = hosp.treso.caretype.cases / sum(hosp.treso.caretype.cases)) %>% 
    mutate(hosp.treso.caretype.beddays.percent = replace_na(hosp.treso.caretype.beddays.percent, 0),
           hosp.treso.caretype.cases.percent = replace_na(hosp.treso.caretype.cases.percent, 0)) %>% 
    select(id, agecluster, caretype, hosp.treso.caretype.beddays.percent, hosp.treso.caretype.cases.percent)
  
  # Calculate CRM at the CD level for scenario year in order to add ALC days to forecasted hospital non-ALC demand
  # Start by collecting demand data at the hospital to treso-origin level
  crm_alc <- select(hospital_lookup_master, year, id, hosp.csd = csduid) %>%
    filter(year == 2016) %>%
    select(-year) %>% 
    group_by(id) %>% 
    mutate(hosp.cd = substr(hosp.csd, 1, 4)) %>%
    left_join(select(crm_treso, id, agecluster, caretype, orig_csd, orig.treso.zone, hosp.treso.cases, hosp.treso.beddays, hosp.treso.bedsneeded), by = c('id')) %>% 
    mutate(hosp.cd = substr(hosp.csd, 1, 4)) %>%
    group_by(id, agecluster) %>%
    ungroup() %>%
    # Join in LTC ALC and non-LTC ALC by CD, and then join in agecluster data
    left_join(select(LTC_proj_cd, cduid, ltc.cd.residual.alc.days = adj.cd.residual.alc.days), by = c('hosp.cd' = 'cduid')) %>%
    left_join(select(ALC_proj_cd_ex_LTC, cduid, ex.ltc.cd.residual.alc.days = adj.cd.residual.alc.days), by = c('hosp.cd' = 'cduid')) %>%
    left_join(ALC_proj_cd_agecluster_ex_LTC, by = c('agecluster', 'hosp.cd' = 'cduid')) %>%
    rename(alc.ex.ltc.cd.demand.days.agepercent = ltc.cd.demand.days.agepercent) %>% 
    left_join(LTC_proj_cd_agecluster, by = c('agecluster', 'hosp.cd' = 'cduid')) %>%
    replace_na(list(alc.ex.ltc.cd.demand.days.agepercent = 0, ltc.cd.demand.days.agepercent = 0)) %>% 
    mutate(cd.ltc.days.agecluster = ltc.cd.residual.alc.days * ltc.cd.demand.days.agepercent,
           alc.ex.ltc.cd.days.agecluster = ex.ltc.cd.residual.alc.days * alc.ex.ltc.cd.demand.days.agepercent) %>% 
    group_by(hosp.cd, agecluster) %>% 
    mutate(hosp.treso.beddays.percent = hosp.treso.beddays / sum(hosp.treso.beddays),
           hosp.treso.beddays.percent = replace_na(hosp.treso.beddays.percent, 0),
           hosp.treso.beddays.ltc = hosp.treso.beddays.percent * cd.ltc.days.agecluster,
           hosp.treso.beddays.alc.ex.ltc = hosp.treso.beddays.percent * alc.ex.ltc.cd.days.agecluster) %>%
    ungroup() %>% 
    select(id, hosp.csd, hosp.cd, agecluster, caretype, orig.csd = orig_csd, orig.treso.zone, hosp.treso.cases, 
           hosp.treso.beddays, hosp.treso.beddays.ltc, hosp.treso.beddays.alc.ex.ltc) %>% 
    mutate(hosp.treso.days.total = hosp.treso.beddays + hosp.treso.beddays.ltc + hosp.treso.beddays.alc.ex.ltc,
           hosp.treso.bedsneeded = hosp.treso.days.total / hosp_op_days) %>%
    # Finalize demand totals by caretype, agecluster, treso zone, hospital
    mutate(demand.days.total = hosp.treso.days.total,
           cases = hosp.treso.cases,
           demand.days.nonALC = hosp.treso.beddays,
           demand.days.LTC = hosp.treso.beddays.ltc,
           demand.days.ALC.ex.LTC = hosp.treso.beddays.alc.ex.ltc,
           beds.needed = demand.days.total/hosp_op_days) %>% 
    select(id, hosp.csd, hosp.cd, caretype, agecluster, orig.csd, orig.treso.zone, cases, demand.days.total, demand.days.nonALC, demand.days.LTC, 
           demand.days.ALC.ex.LTC, beds.needed)
  
  # Calculate number of cases, demand days, and beds needed by hospital, caretype, and agecluster
  crm_hosp_agecluster <- crm_alc %>%
    group_by(id, agecluster, caretype) %>%
    summarise_at(vars(cases:beds.needed), sum) %>%
    ungroup()
  
  # Calculate number of cases, demand days, and beds needed by CSD, caretype, and agecluster
  crm_csd_agecluster <- crm_alc %>%
    group_by(orig.csd, agecluster, caretype) %>%
    summarise_at(vars(cases:beds.needed), sum) %>%
    ungroup()
  
  # Calculate number of cases, demand days, and beds needed by hospital and caretype
  crm_hosp <- crm_alc %>% 
    group_by(id, caretype) %>%
    summarise_at(vars(cases:beds.needed), sum) %>%
    ungroup()
  
  # Calculate number of cases, demand days, and beds needed by origin CSD, caretype
  crm_orig <- crm_alc %>% 
    group_by(orig.csd, id, hosp.csd, caretype, agecluster) %>%
    summarise_at(vars(cases:beds.needed), sum) %>%
    ungroup()
  
  return(list(crm_alc, crm_orig))
}

moh_decision_making <- function(num_beds, hospitallocations_tz, crm_alc, crm_orig, utilization_targets) {
  # Determine total existing bed count by CSD by caretype
  csd_beds_existing <- num_beds %>% 
    left_join(select(hospitallocations_tz, id, hosp.csd = csduid), by = ('id')) %>% 
    rename(beds.existing = staffedbeds) %>% 
    mutate(hosp.cd = substr(hosp.csd, 1, 4)) %>%
    # Analyse by CSD and CD
    group_by(hosp.cd, hosp.csd, caretype) %>%
    summarise(csd.beds.existing = sum(beds.existing)) %>%
    ungroup()
  
  # Determine total existing bed count by CD by caretype
  cd_beds_existing <- csd_beds_existing %>%   
    group_by(hosp.cd, caretype) %>% 
    summarise(cd.beds.existing = sum(csd.beds.existing)) %>% 
    ungroup()
  
  # Calculate number of cases, demand days, and beds needed by CSD, caretype, and agecluster
  crm_csd_agecluster <- crm_alc %>%
    group_by(orig.csd, agecluster, caretype) %>%
    summarise_at(vars(cases:beds.needed), sum) %>%
    ungroup()
  
  # Join demand to capacity data and assess demand and capacity by CSD and CD
  demand_capacity_difference <- crm_csd_agecluster %>% 
    mutate(orig.cd = substr(orig.csd, 1, 4)) %>%
    left_join(csd_beds_existing, by = c('orig.csd' = 'hosp.csd', 'caretype')) %>% 
    left_join(cd_beds_existing, by = c('orig.cd' = 'hosp.cd', 'caretype')) %>% 
    replace_na(list(csd.beds.existing = 0, cd.beds.existing = 0)) %>%
    # Join in utilization targets provided by user
    left_join(utilization_targets, by = c('caretype')) %>% 
    # Calculate total beds needed across all ageclusters
    group_by(orig.csd, caretype) %>% 
    mutate(csd.beds.needed = sum(beds.needed),
           adjusted.csd.beds.needed = round(csd.beds.needed / target, 0),
           csd.beds.to.build = adjusted.csd.beds.needed - csd.beds.existing) %>% 
    group_by(orig.cd, caretype) %>% 
    mutate(cd.beds.needed = sum(beds.needed),
           adjusted.cd.beds.needed = round(cd.beds.needed / target, 0),
           cd.beds.to.build = adjusted.cd.beds.needed - cd.beds.existing) %>%
    ungroup()
  
  # Test number of beds required by CSD/CD if maintaining existing travel patterns
  demand_capacity_crm_patterns <- crm_orig %>% 
    mutate(hosp.cd = substr(hosp.csd, 1, 4)) %>%
    left_join(select(num_beds, caretype, id, beds.existing = staffedbeds), by = c('id', 'caretype')) %>% 
    replace_na(list(beds.existing = 0)) %>% 
    # Calculate the number of beds needed by caretype by hospital
    group_by(hosp.cd, hosp.csd, id, caretype, beds.existing) %>% 
    # mutate(hosp.beds.needed = sum(beds.needed)) %>% 
    summarise(hosp.beds.needed = sum(beds.needed)) %>% 
    ungroup() %>% 
    left_join(utilization_targets, by = c('caretype')) %>% 
    mutate(adjusted.hosp.beds.needed = round(hosp.beds.needed / target, 0),
           hosp.beds.to.build = adjusted.hosp.beds.needed - beds.existing)
  
  # Calculate overburden by CSD by hospital
  overburden_csd_hosp <- demand_capacity_crm_patterns %>% 
    select(id, caretype, hosp.beds.to.build) %>% 
    filter(hosp.beds.to.build > 0) %>% 
    left_join(select(crm_orig, orig.csd, id, hosp.csd, caretype, agecluster, cases:beds.needed), by = c('id', 'caretype')) %>% 
    # Sum cases, demand-days across ageclusters to find total demand by orig.csd, caretype
    mutate(sumCases = sum(cases)) %>% 
    group_by(orig.csd, id, caretype) %>% 
    mutate(cases = sum(cases), demand.days.total = sum(demand.days.total), demand.days.nonALC = sum(demand.days.nonALC), 
           demand.days.LTC = sum(demand.days.LTC), demand.days.ALC.ex.LTC = sum(demand.days.ALC.ex.LTC), beds.needed = sum(beds.needed)) %>% 
    select(-agecluster) %>%
    ungroup() %>% 
    distinct() %>% 
    mutate(tot.cases = sum(cases)) %>% 
    # filter(orig.csd != hosp.csd) %>% 
    group_by(orig.csd, caretype) %>% 
    # Caclulate sum of beds.needed by origin CSD 
    mutate(csd.caretype.beds.needed = sum(beds.needed)) %>% 
    ungroup()
  
  overburden_csd_hosp_new <- overburden_csd_hosp %>% 
    anti_join(select(hospitallocations_tz, csduid), by = c('orig.csd' = 'csduid'))
  
  return(list(demand_capacity_difference, overburden_csd_hosp_new))
}


moh_new_hospitals <- function(overburden_csd_hosp_new, treso_zone_system, treso_population_moh_agecluster, treso_centroids,
                              new_hospital_user = NULL, min_beds_toggle, min_beds, scenario_year) {
  
  # Initialize new_hospital dataframe
  new_hospital <- tibble(id = integer(), hospital.lat = double(), hospital.long = double(), name = character(), 
                         caretype = character(), beds = integer(), specialityfactor = integer(), total.hosp.beds = integer())
  
  # Calculate TRESO populations at the agecluster level instead of agegroup
  colname <- noquote(paste0('population.', scenario_year))
  treso_population_moh_agecluster_scenario <- treso_population_moh_agecluster %>% 
    select(-age.group) %>% 
    group_by(treso_zone, agecluster, csduid) %>% 
    summarise_all(sum) %>% 
    select(csduid, treso_zone, agecluster, !!as.symbol(colname)) %>% 
    rename(scen.proj.pop.treso = !!as.symbol(colname)) %>% 
    ungroup()
  
  new_hospital_dm <- overburden_csd_hosp_new %>% 
    distinct(orig.csd, caretype, csd.caretype.beds.needed) %>%
    mutate(csd.caretype.beds.needed = ceiling(csd.caretype.beds.needed)) %>%
    # Calculate sum of beds needed by CSD
    group_by(orig.csd) %>% 
    mutate(csd.beds.needed = sum(csd.caretype.beds.needed)) %>% 
    ungroup() %>% 
    # Limit new hospital beds to be built only if minimum size threshold is crossed
    mutate(beds.threshold = pmax(min_beds_toggle * csd.beds.needed, csd.caretype.beds.needed)) %>% 
    filter(beds.threshold >= min_beds) %>% 
    rename(beds.to.build = csd.caretype.beds.needed) %>% 
    mutate(specialityfactor = 0) %>% 
    # # Assign ID number to new hospitals; assume max of 1 new hospital per cSD
    mutate(row = group_indices(., orig.csd),
           id = 1111 * row) %>% 
    select(-row) %>%
    # Determine most populous TRESO zone in CSD
    left_join(select(treso_zone_system, csduid, treso_id), by = c('orig.csd' = 'csduid')) %>% 
    left_join(select(treso_population_moh_agecluster_scenario, treso_zone, agecluster, scen.proj.pop.treso), by = c('treso_id' = 'treso_zone')) %>% 
    rename(pop.treso = scen.proj.pop.treso) %>% 
    group_by(treso_id) %>% 
    mutate(pop.treso.all.ages = sum(pop.treso)) %>% 
    select(-agecluster, -pop.treso) %>% 
    distinct() 
  
  # Check to see if any new hospitals are required
  if (nrow(new_hospital_dm) > 0) {
    # Build new hospital(s) in most populous TRESO zone of CSDs needing a new hospital
    new_hospital_dm_treso <- new_hospital_dm %>% 
      arrange(desc(pop.treso.all.ages)) %>% 
      group_by(orig.csd, caretype) %>% 
      mutate(row = 1:n()) %>% 
      filter(row == 1) %>% 
      select(-row) %>% 
      ungroup() %>% 
      # Join in treso centroid to get lat/long for new hospitals
      left_join(select(treso_centroids, treso_id, lat, long), by = c('treso_id')) %>% 
      # Generate hospital name based on ID
      mutate(name = paste0('Hospital ', id)) %>% 
      # Rename and remove fields to match inputs for redistribution
      rename(beds = beds.to.build, total.hosp.beds = csd.beds.needed, hospital.lat = lat, hospital.long = long) %>% 
      select(id, hospital.lat, hospital.long, name, caretype, beds, specialityfactor, total.hosp.beds) 
    
    # Add in row for Emergency cases (AM) with 0 beds for each new hospital
    new_hospital <- new_hospital_dm_treso %>% 
      # Start by creating list of new hospitals and adding a column for 'AM' caretype
      distinct(id) %>% 
      # Set number of beds for Emerg (AM): 0 = no emerg, -1 = emerg to be built
      mutate(caretype = 'AM', beds = -1) %>%
      # Fill in blanks for AM caretype based on other caretypes
      bind_rows(., new_hospital_dm_treso) %>% 
      arrange(id, caretype) %>% 
      group_by(id) %>% 
      mutate(hospital.lat = last(hospital.lat), hospital.long = last(hospital.long), name = last(name), 
             specialityfactor = 0, total.hosp.beds = last(total.hosp.beds)) %>% 
      ungroup()
    
    if (!is.null(new_hospital_user)) {
      # Combine user-set new hospitals plus decision-making new hospitals (if running decision-making)
      new_hospital <- bind_rows(new_hospital, new_hospital_user)
    }
    
  }
  return(new_hospital)
}

redistribute_demand_leaving_for_new_hospitals <- function(new_hospital, speciality_hospitals, num_beds, crm_orig,
                                                          treso_shp, treso_zone_system, treso_population_moh_agecluster, 
                                                          age_cluster_list, stay_in_csd_factor, scenario_year) {
  new_hospital <- crossing(new_hospital, age_cluster_list) %>% 
    rename(agecluster = age_cluster_list)
  
  # Overlay with TRESO Shapefile to get the TRESO/CSD Information for new hospitals
  new_hospital_xy <- create_hospital_xy(new_hospital) 
  
  new_hospital_overlay <- create_overlay(new_hospital_xy, treso_shp, type = "hospital") %>% 
    left_join(treso_zone_system, by = c("treso.id.pos" = "treso_id")) %>%
    distinct() %>% 
    left_join(new_hospital, by = "id") %>% 
    mutate_if(is.factor, as.character) %>% 
    select(-lat, -long, -csdname, -area, -mof_region, -hospital.lat, -hospital.long, -cdname) %>% 
    rename(treso_zone = treso.id.pos, hosp.csd = csduid, hospitalname = name, hosp.cd = cduid) %>% 
    mutate(id = as.numeric(id), hosp.cd=as.character(hosp.cd)) %>% 
    mutate(new.hospital = 1)
  
  # Calculate TRESO populations at the agecluster level instead of agegroup
  colname <- noquote(paste0('population.', scenario_year))
  treso_population_moh_agecluster_scenario <- treso_population_moh_agecluster %>% 
    select(-age.group) %>% 
    group_by(treso_zone, agecluster, csduid) %>% 
    summarise_all(sum) %>% 
    select(csduid, treso_zone, agecluster, !!as.symbol(colname)) %>% 
    rename(scen.proj.pop.treso = !!as.symbol(colname)) %>% 
    ungroup()
  
  # Identify the TRESO zones within the CSD of each new hospital, and join in populations
  new_hospital_CSD_treso <- distinct(new_hospital_overlay, id, hosp.csd) %>% 
    left_join(select(treso_zone_system, csduid, treso_id), by = c('hosp.csd' = 'csduid')) %>% 
    left_join(select(treso_population_moh_agecluster_scenario, treso_zone, agecluster, scen.proj.pop.treso), by = c('treso_id' = 'treso_zone')) %>% 
    rename(pop.treso = scen.proj.pop.treso, csduid = hosp.csd) %>% 
    # Determine percentage of population by agecluster by treso zone per CSD
    group_by(csduid, agecluster) %>% 
    mutate(pop.percent = pop.treso / sum(pop.treso)) %>% 
    ungroup()
  
  # For CSDs with existing hospitals: redistribute based on weight, eg:
  # Existing hospitals in CSD A = 1,000 beds; other CSDs attracting cases from CSD A = 10,000 beds;
  # New hospital in CSD A = 100 beds
  # Total cases = 5,000 at existing hospitals in CSD A, 5,000 at hospitals outside CSD A
  # Cases per bed = 5,000 / 1,000 at existing hosp. CSD A = 5 cases/bed; 5,000/10,000 = 0.5 cases/bed outside CSD A
  # Assume new beds in CSD A = 5 cases/bed as per existing
  # New cases per bed = 1,100 beds x 5 cases = 5,500 cases CSD A; 5,000 cases non-CSD-A; total cases should = 10,000
  # Therefore: (1,100 x 5) * (10,000 / (1,110 * 5 + 10,000 * 0.5)) = 5,238 in CSD A
  # Therefore: 10,000 - 5,238 = 4,762 cases leaving CSD A for other CSDs
  # Now redistribute cases in CSD A based on new hospital
  # Then redistribute cases at hospitals outside CSD A (coming from CSD A) based on the above
  
  # Determine the number of cases/demand-days staying within CSDs vs. the number with origin CSD != hospital CSD, 
  # and then produce a rate at which patients choose to stay or leave their origin CSD based on the number of beds per hospital
  crm_staying <- crm_orig %>% 
    # Filter CRM outputs to CSDs with new hospitals being built
    filter(orig.csd %in% distinct(new_hospital_overlay, hosp.csd)[[1]]) %>% 
    # Identify existing hospitals which are speciality hospitals, and join in TRESO/CSD info
    left_join(select(speciality_hospitals, id, specialityfactor), by = c('id')) %>% 
    #left_join(select(hospitallocations_tz, id, hosp.csd = csduid), by = c('id')) %>% 
    mutate(new.hospital = 0) %>% 
    # Join in bed counts for existing hospitals
    left_join(select(num_beds, id, beds = staffedbeds, caretype), by = c('id', 'caretype')) %>%
    # Assign a representative number of beds to AM (Emergency) caretype to weight distribution of AM cases in later steps
    group_by(id) %>% 
    mutate(beds = replace_na(beds, 0), 
           beds = ifelse(beds == 0, sum(beds), beds)) %>% 
    # Calculate number of beds in CSDs with new hospitals, vs. number of beds in other CSDs which attracted at least some number of  residents from the CSD with new hospitals
    group_by(orig.csd, caretype, agecluster) %>% 
    mutate(beds = replace_na(beds, 0),
           beds.in = sum(beds[hosp.csd == orig.csd]),
           beds.out = sum(beds) - beds.in) %>% 
    # Calculate number of cases, demand days in CSDs with new hospitals vs. other CSDs
    mutate(cases.in = sum(cases[hosp.csd == orig.csd]),
           cases.out = sum(cases) - cases.in,
           demand.days.total.in = sum(demand.days.total[hosp.csd == orig.csd]),
           demand.days.total.out = sum(demand.days.total) - demand.days.total.in,
           demand.days.nonALC.in = sum(demand.days.nonALC[hosp.csd == orig.csd]),
           demand.days.nonALC.out = sum(demand.days.nonALC) - demand.days.nonALC.in,
           demand.days.LTC.in = sum(demand.days.LTC[hosp.csd == orig.csd]),
           demand.days.LTC.out = sum(demand.days.LTC) - demand.days.LTC.in,
           demand.days.ALC.ex.LTC.in = sum(demand.days.ALC.ex.LTC[hosp.csd == orig.csd]),
           demand.days.ALC.ex.LTC.out = sum(demand.days.ALC.ex.LTC) - demand.days.ALC.ex.LTC.in) %>%
    # Tested shorter code but not working so stuck with long form instead: mutate_at(vars(cases:demand.days.ALC), .funs = list(inin = ~sum(vars(cases:demand.days.ALC)[hosp.csd == orig.csd])))
    # Calculate rates of cases/days per bed in and bed out of CSDs with new hospitals
    mutate(cases.per.bed.in = cases.in / beds.in,
           cases.per.bed.out = cases.out / beds.out,
           demand.days.total.per.bed.in = demand.days.total.in / beds.in,
           demand.days.total.per.bed.out = demand.days.total.out / beds.out,
           demand.days.nonALC.per.bed.in = demand.days.nonALC.in / beds.in,
           demand.days.nonALC.per.bed.out = demand.days.nonALC.out / beds.out,
           demand.days.LTC.per.bed.in = demand.days.LTC.in / beds.in,
           demand.days.LTC.per.bed.out = demand.days.LTC.out / beds.out,
           demand.days.ALC.ex.LTC.per.bed.in = demand.days.ALC.ex.LTC.in / beds.in,
           demand.days.ALC.ex.LTC.per.bed.out = demand.days.ALC.ex.LTC.out / beds.out) %>% 
    # Fix N/A values if there are no existing hospitals in CSDs with a new hospital being added
    replace_na(list(cases.per.bed.in = 0, cases.per.bed.out = 0, demand.days.total.per.bed.in = 0, demand.days.total.per.bed.out = 0,
                    demand.days.nonALC.per.bed.in = 0, demand.days.nonALC.per.bed.out = 0, demand.days.LTC.per.bed.in = 0, demand.days.ALC.ex.LTC.per.bed.in = 0,
                    demand.days.LTC.per.bed.out = 0, demand.days.ALC.ex.LTC.per.bed.out = 0)) %>% 
    # Calculate total cases, demand-days by origin CSD, caretype
    mutate(orig.cases = sum(cases),
           orig.demand.days.total = sum(demand.days.total),
           orig.demand.days.nonALC = sum(demand.days.nonALC),
           orig.demand.days.LTC = sum(demand.days.LTC),
           orig.demand.days.ALC.ex.LTC = sum(demand.days.ALC.ex.LTC)) %>% 
    ungroup()
  
  # Apply rates to new hospitals
  new_hospital_rates <- new_hospital_overlay %>% 
    select(-hospitalname, -treso_zone, -hosp.cd) %>% 
    distinct() %>% 
    left_join(distinct(crm_staying, orig.csd, caretype, agecluster, cases.per.bed.in, demand.days.total.per.bed.in, 
                       demand.days.nonALC.per.bed.in, demand.days.LTC.per.bed.in, demand.days.ALC.ex.LTC.per.bed.in, cases.per.bed.out, 
                       demand.days.total.per.bed.out, demand.days.nonALC.per.bed.out, demand.days.LTC.per.bed.out, demand.days.ALC.ex.LTC.per.bed.out, 
                       orig.cases, orig.demand.days.total, orig.demand.days.nonALC, orig.demand.days.LTC, orig.demand.days.ALC.ex.LTC), 
              by = c('hosp.csd' = 'orig.csd', 'caretype', 'agecluster')) %>% 
    # Solve any N/A values based on lack of historical demand
    replace_na(list(demand.days.total.per.bed.in = 0, demand.days.nonALC.per.bed.in = 0, demand.days.LTC.per.bed.in = 0, demand.days.ALC.ex.LTC.per.bed.in = 0,
                    cases.per.bed.out = 0, demand.days.total.per.bed.out = 0, demand.days.nonALC.per.bed.out = 0, demand.days.LTC.per.bed.out = 0, 
                    demand.days.ALC.ex.LTC.per.bed.out = 0, orig.cases = 0, orig.demand.days.total = 0, orig.demand.days.nonALC = 0, orig.demand.days.LTC = 0, 
                    orig.demand.days.ALC.ex.LTC = 0)) %>% 
    # For caretypes with no existing beds in the CSD, assume likelihood of staying in CSD is a function of likelihood of going to external CSDs
    mutate(cases.per.bed.in = ifelse(cases.per.bed.in == 0, stay_in_csd_factor * cases.per.bed.out, cases.per.bed.in),
           demand.days.total.per.bed.in = ifelse(demand.days.total.per.bed.in == 0, stay_in_csd_factor * demand.days.total.per.bed.out, demand.days.total.per.bed.in),
           demand.days.nonALC.per.bed.in = ifelse(demand.days.nonALC.per.bed.in == 0, stay_in_csd_factor * demand.days.nonALC.per.bed.out, demand.days.nonALC.per.bed.in),
           demand.days.LTC.per.bed.in = ifelse(demand.days.LTC.per.bed.in == 0, stay_in_csd_factor * demand.days.LTC.per.bed.out, demand.days.LTC.per.bed.in),
           demand.days.ALC.ex.LTC.per.bed.in = ifelse(demand.days.ALC.ex.LTC.per.bed.in == 0, stay_in_csd_factor * demand.days.ALC.ex.LTC.per.bed.out, demand.days.ALC.ex.LTC.per.bed.in)) %>%
    select(-cases.per.bed.out, -demand.days.total.per.bed.out, -demand.days.nonALC.per.bed.out, -demand.days.LTC.per.bed.out, -demand.days.ALC.ex.LTC.per.bed.out) %>% 
    # Assign a representative number of beds to AM (Emergency) caretype to weight distribution of AM cases in later steps
    group_by(id) %>% 
    mutate(beds = ifelse(beds == -1 & caretype == 'AM', sum(beds), beds)) %>% 
    ungroup() %>% 
    # Estimate share of cases/demand-days going to new hospitals
    mutate(cases = cases.per.bed.in * beds,
           demand.days.total = demand.days.total.per.bed.in * beds,
           demand.days.nonALC = demand.days.nonALC.per.bed.in * beds,
           demand.days.LTC = demand.days.LTC.per.bed.in * beds,
           demand.days.ALC.ex.LTC = demand.days.ALC.ex.LTC.per.bed.in * beds)
  
  # Bind in new hospitals to estimate how many people will be converted to stay in the CSD
  crm_staying_new <- crm_staying %>% 
    # Bind in new hospital data
    bind_rows(., new_hospital_rates) %>% 
    mutate(orig.csd = ifelse(is.na(orig.csd), hosp.csd, orig.csd)) %>% 
    group_by(orig.csd, caretype, agecluster) %>% 
    # Determine market share per hospital
    mutate(orig.cases.percent = cases / sum(cases),
           orig.demand.days.total.percent = demand.days.total / sum(demand.days.total),
           orig.demand.days.nonALC.percent = demand.days.nonALC / sum(demand.days.nonALC),
           orig.demand.days.LTC.percent = demand.days.LTC / sum(demand.days.LTC),
           orig.demand.days.ALC.ex.LTC.percent = demand.days.ALC.ex.LTC / sum(demand.days.ALC.ex.LTC)) %>% 
    ungroup() %>%
    # Fix N/A values for groupings with no cases/demand
    mutate(orig.cases.percent = replace_na(orig.cases.percent, 0),
           orig.demand.days.total.percent = replace_na(orig.demand.days.total.percent, 0),
           orig.demand.days.nonALC.percent = replace_na(orig.demand.days.nonALC.percent, 0),
           orig.demand.days.LTC.percent = replace_na(orig.demand.days.LTC.percent, 0),
           orig.demand.days.ALC.ex.LTC.percent = replace_na(orig.demand.days.ALC.ex.LTC.percent, 0)) %>% 
    # Calculate revised number of cases based on percentages above
    mutate(cases.revised = orig.cases.percent * orig.cases,
           demand.days.total.revised = orig.demand.days.total.percent * orig.demand.days.total,
           demand.days.nonALC.revised = orig.demand.days.nonALC.percent * orig.demand.days.nonALC,
           demand.days.LTC.revised = orig.demand.days.LTC.percent * orig.demand.days.LTC,
           demand.days.ALC.ex.LTC.revised = orig.demand.days.ALC.ex.LTC.percent * orig.demand.days.ALC.ex.LTC)
  
  # Summarise results of redistribution of cases/demand-days from other CSDs to CSDs with new hospitals, and divide into TRESO zones  
  crm_orig_revised <- crm_staying_new %>% 
    select(orig.csd, id, hosp.csd, caretype, agecluster, specialityfactor, new.hospital, beds, cases.revised, demand.days.total.revised, 
           demand.days.nonALC.revised, demand.days.LTC.revised, demand.days.ALC.ex.LTC.revised) %>% 
    # Replace N/A values for AM cases (since no demand-days for AM, aka Emergency)
    mutate(demand.days.total.revised = replace_na(demand.days.total.revised, 0), 
           demand.days.nonALC.revised = replace_na(demand.days.nonALC.revised, 0), 
           demand.days.LTC.revised = replace_na(demand.days.LTC.revised, 0),
           demand.days.ALC.ex.LTC.revised = replace_na(demand.days.ALC.ex.LTC.revised, 0)) %>% 
    # Divide demand from CSDs down to individual TRESO zones
    left_join(select(new_hospital_CSD_treso, csduid, agecluster, pop.percent, treso_id), by = c('orig.csd' = 'csduid', 'agecluster')) %>% 
    mutate(cases.revised = cases.revised * pop.percent,
           demand.days.total.revised = demand.days.total.revised * pop.percent,
           demand.days.nonALC.revised = demand.days.nonALC.revised * pop.percent,
           demand.days.LTC.revised = demand.days.LTC.revised * pop.percent,
           demand.days.ALC.ex.LTC.revised = demand.days.ALC.ex.LTC.revised * pop.percent) %>% 
    rename(orig.treso.zone = treso_id)
  
  return(crm_orig_revised)
}
redistribute_demand_within_for_new_hospitals <- function(new_hospital, speciality_hospitals, num_beds, crm_alc, crm_orig_revised) {
  # Revise cases and demand-days in crm_alc based on redistribution from other CSDs due to new hospitals
  crm_alc_redist <- crm_alc %>% 
    full_join(select(crm_orig_revised, orig.csd, orig.treso.zone, hosp.csd, id, caretype, agecluster, cases.revised, 
                     demand.days.total.revised, demand.days.nonALC.revised, demand.days.LTC.revised, demand.days.ALC.ex.LTC.revised), 
              by = c('orig.treso.zone', 'id', 'orig.csd', 'hosp.csd', 'caretype', 'agecluster')) %>% 
    # Reuse existing cases, demand-days where cases.revised and demand.days.revised are empty (i.e., no changes)
    mutate(cases.revised = ifelse(is.na(cases.revised), cases, cases.revised),
           demand.days.total.revised = ifelse(is.na(demand.days.total.revised), demand.days.total, demand.days.total.revised),
           demand.days.nonALC.revised = ifelse(is.na(demand.days.nonALC.revised), demand.days.nonALC, demand.days.nonALC.revised),
           demand.days.LTC.revised = ifelse(is.na(demand.days.LTC.revised), demand.days.LTC, demand.days.LTC.revised),
           demand.days.ALC.ex.LTC.revised = ifelse(is.na(demand.days.ALC.ex.LTC.revised), demand.days.ALC.ex.LTC, demand.days.ALC.ex.LTC.revised)) %>% 
    select(id, hosp.csd, orig.csd, orig.treso.zone, agecluster, caretype, cases.revised, demand.days.total.revised, 
           demand.days.nonALC.revised, demand.days.LTC.revised, demand.days.ALC.ex.LTC.revised) %>% 
    rename(cases = cases.revised, demand.days.total = demand.days.total.revised, demand.days.nonALC = demand.days.nonALC.revised, 
           demand.days.LTC = demand.days.LTC.revised, demand.days.ALC.ex.LTC = demand.days.ALC.ex.LTC.revised)
  
  # Summarise revised crm_alc data based on redistributions within CSDs with new hospitals
  crm_hosp_agecluster_redist <- crm_alc_redist %>%
    group_by(id, orig.csd, hosp.csd, agecluster, caretype) %>%
    summarise_at(vars(cases:demand.days.ALC.ex.LTC), sum) %>%
    ungroup()
  
  crm_hosp_agecluster_all <- crm_hosp_agecluster_redist %>%
    # Identify existing hospitals which are speciality hospitals
    left_join(select(speciality_hospitals, id, specialityfactor), by = c('id')) %>%
    mutate(new.hospital = 0) %>% 
    # Join in bed counts for existing hospitals
    left_join(select(num_beds, id, beds = staffedbeds, caretype), by = c('id', 'caretype')) %>%
    # Join in total bed counts for existing hospitals
    left_join(distinct(num_beds, id, total.hosp.beds), by = c('id')) %>% 
    # Address existing hospitals with missing bed counts by estimating beds based on cases relative to other hospitals in CSD
    mutate(beds = ifelse(caretype == 'AM' & id < 999, -1, beds)) %>% 
    # Join in bed counts and speciality hospital flag for new hospitals
    left_join(distinct(new_hospital, id, caretype, beds, specialityfactor), by = c('id', 'caretype')) %>% 
    # Join in total bed counts for new hospitals
    left_join(distinct(new_hospital, id, total.hosp.beds), by = c('id')) %>% 
    rename(specialityfactor = specialityfactor.x, beds = beds.x, total.hosp.beds = total.hosp.beds.x) %>% 
    mutate(specialityfactor = ifelse(is.na(specialityfactor), specialityfactor.y, specialityfactor),
           beds = ifelse(is.na(beds), beds.y, beds),
           beds = replace_na(beds, 0),
           total.hosp.beds = ifelse(is.na(total.hosp.beds), total.hosp.beds.y, total.hosp.beds),
           total.hosp.beds = replace_na(total.hosp.beds, 0)) %>% 
    select(-specialityfactor.y, -beds.y, -total.hosp.beds.y) %>%
    # Flag new hospitals as new
    mutate(new.hospital = replace_na(new.hospital, 1)) %>% 
    # Estimate total.hosp.beds for hospitals with no bed data provided, in order to set weighting for AM (emerg) beds
    group_by(id) %>% 
    mutate(AM.cases.no.bed.data = sum(cases[caretype == 'AM' & id < 999])) %>% 
    group_by(hosp.csd) %>% 
    mutate(AM.cases.bed.data = sum(cases[caretype == 'AM' & id < 999]) - AM.cases.no.bed.data,
           AM.cases.weight = ifelse(AM.cases.bed.data == 0, 1, AM.cases.no.bed.data / AM.cases.bed.data),
           AM.cases.weight = replace_na(AM.cases.weight, 0),
           sum.tot.beds = sum(total.hosp.beds[id < 999]),
           total.hosp.beds = ifelse(total.hosp.beds == 0, AM.cases.weight * sum(total.hosp.beds[id < 999]), total.hosp.beds),
           total.hosp.beds = ifelse(total.hosp.beds == 0, 100, total.hosp.beds)) %>% 
    ungroup() %>% 
    # Set relative size of each Emergency department relative to Emergency Dept's at other hospitals by using sum of inpatient beds (i.e., total hospital beds) as a proxy. Note that this is only used as a weighting factor relative to other Emergency Dept's, so it is the relative size that counts, not the total number of beds
    group_by(id) %>% 
    mutate(beds = ifelse(beds == -1, total.hosp.beds, beds)) %>% 
    # Adjust number of beds for existing hospitals where cases/demand was nonzero in historical years but where the bedcount = 0
    mutate(beds = ifelse(beds == 0 & id < 300, 1, beds)) %>% 
    ungroup() %>% 
    # Remove existing hospitals which are speciality hospitals from redistribution
    # filter(new.hospital == 1 | specialityfactor == 0) %>%
    # Calculate percentage of beds by hospital by caretype by agecluster
    group_by(orig.csd, hosp.csd, caretype, agecluster) %>%
    mutate(beds.percent = beds / sum(beds)) %>% 
    mutate(beds.percent = replace_na(beds.percent, 0)) %>% 
    # Sum cases, bed-days, etc. by CSD
    mutate(cases.csd.caretype.age = sum(cases),
           demand.days.total.csd.caretype.age = sum(demand.days.total),
           demand.days.LTC.csd.caretype.age = sum(demand.days.LTC),
           demand.days.ALC.ex.LTC.csd.caretype.age = sum(demand.days.ALC.ex.LTC),
           demand.days.nonALC.csd.caretype.age = sum(demand.days.nonALC)) %>% 
    mutate(cases.revised = beds.percent * cases.csd.caretype.age,
           demand.days.total.revised = beds.percent * demand.days.total.csd.caretype.age,
           demand.days.LTC.revised = beds.percent * demand.days.LTC.csd.caretype.age,
           demand.days.ALC.ex.LTC.revised = beds.percent * demand.days.ALC.ex.LTC.csd.caretype.age,
           demand.days.nonALC.revised = beds.percent * demand.days.nonALC.csd.caretype.age) %>% 
    ungroup()
}

format_demand_output <- function(treso_shp, treso_zone_system, treso_population_moh_agecluster, travel_time_skim,
                                 travel_distance_skim, hospitallocations_tz, new_hospital=NULL, crm_hosp_agecluster_all,
                                 utilization_targets, trips_per_case, scenario_year) {
  
  # Calculate TRESO populations at the agecluster level instead of agegroup
  colname <- noquote(paste0('population.', scenario_year))
  treso_population_moh_agecluster_scenario <- treso_population_moh_agecluster %>% 
    select(-age.group) %>% 
    group_by(treso_zone, agecluster, csduid) %>% 
    summarise_all(sum) %>% 
    select(csduid, treso_zone, agecluster, !!as.symbol(colname)) %>% 
    rename(scen.proj.pop.treso = !!as.symbol(colname)) %>% 
    ungroup()
  
  # Simplify demand outputs
  crm_demand <- crm_hosp_agecluster_all %>% 
    select(-cases, -demand.days.total, -demand.days.nonALC, -demand.days.LTC, -demand.days.ALC.ex.LTC) %>% 
    rename(cases = cases.revised, demand.days.total = demand.days.total.revised, demand.days.nonALC = demand.days.nonALC.revised,
           demand.days.LTC = demand.days.LTC.revised, demand.days.ALC.ex.LTC = demand.days.ALC.ex.LTC.revised) %>% 
    # Join in TRESO population data to distribute cases from CSDs to TRESO zones
    left_join(select(treso_zone_system, csduid, treso_id), by = c('orig.csd' = 'csduid')) %>% 
    left_join(select(treso_population_moh_agecluster_scenario, treso_zone, agecluster, scen.proj.pop.treso), by = c('treso_id' = 'treso_zone', 'agecluster')) %>% 
    rename(pop.treso = scen.proj.pop.treso, orig.treso = treso_id) %>% 
    # Determine percentage of population by agecluster by treso zone per CSD
    replace_na(list(pop.treso = 0)) %>% 
    group_by(orig.csd, id, caretype, agecluster) %>% 
    mutate(pop.percent = pop.treso / sum(pop.treso)) %>% 
    ungroup() %>%
    replace_na(list(pop.percent = 0)) %>% 
    # Distribute Origin CSD cases/demand-days into TRESO zones
    mutate(cases = cases * pop.percent, demand.days.total = demand.days.total * pop.percent, 
           demand.days.nonALC = demand.days.nonALC * pop.percent, demand.days.LTC = demand.days.LTC * pop.percent, 
           demand.days.ALC.ex.LTC = demand.days.ALC.ex.LTC * pop.percent) %>% 
    mutate(hosp.cd = substr(hosp.csd, 1, 4)) %>%
    select(id, hosp.cd, hosp.csd, specialityfactor, new.hospital, orig.csd, orig.treso, caretype, agecluster, pop.treso, cases, 
           demand.days.total, demand.days.nonALC, demand.days.LTC, demand.days.ALC.ex.LTC) %>% 
    # Join in utilization targets
    left_join(utilization_targets, by = c('caretype')) %>% 
    # Calculate beds needed by hospital, caretype
    group_by(id, caretype) %>% 
    mutate(beds.needed = ceiling(sum(demand.days.total) / HOSP_OP_DAYS / target)) %>% 
    ungroup()
  
  if (is.null(new_hospital)) {
    new_hospital_overlay <- tibble(treso_zone = double(), id = double(), hosp.csd = double(), hosp.cd = double(), 
                                   caretype = character(), beds = double(), hospitalname = character(),
                                   specialtyfactor = double(), total.hosp.beds = double(), agecluster = character(),
                                   new.hospital = double())
  } else {
    # Overlay with TRESO Shapefile to get the TRESO/CSD Information for new hospitals
    new_hospital_xy <- create_hospital_xy(new_hospital) 
    
    new_hospital_overlay <- create_overlay(new_hospital_xy, treso_shp, type = "hospital") %>% 
      left_join(treso_zone_system, by = c("treso.id.pos" = "treso_id")) %>%
      distinct() %>% 
      left_join(new_hospital, by = "id") %>% 
      mutate_if(is.factor, as.character) %>% 
      select(-lat, -long, -csdname, -area, -mof_region, -hospital.lat, -hospital.long, -cdname) %>% 
      rename(treso_zone = treso.id.pos, hosp.csd = csduid, hospitalname = name, hosp.cd = cduid) %>% 
      mutate(id = as.numeric(id), hosp.cd=as.character(hosp.cd)) %>% 
      mutate(new.hospital = 1)
  }
  
  # Add travel times and distances to demand
  crm_demand_travel <- crm_demand %>% 
    # Join in treso.zone info for existing hospitals
    left_join(select(hospitallocations_tz, id, hosp.treso = treso_zone), by = c('id')) %>%
    # Join in treso.zone info for new hospitals
    left_join(distinct(new_hospital_overlay, id, treso_zone), by = c('id')) %>% 
    # Merge hosp.treso column (containing treso zones for existing hospitals) with treso_id column (containing treso zones for new hospitals)
    replace_na(list(hosp.treso = 0, treso_zone = 0)) %>% 
    mutate(hosp.treso = pmax(hosp.treso, treso_zone)) %>% 
    select(-treso_zone) %>% 
    # Join in travel times for treso zone to treso zone travel
    left_join(travel_time_skim, by = c('orig.treso' = 'treso.id.por', 'hosp.treso' = 'treso.id.pos')) %>%
    rename(trip.time = value) %>%
    # Join in travel distan for treso zone to treso zone travel
    left_join(travel_distance_skim, by=c('orig.treso' = 'treso.id.por', 'hosp.treso' = 'treso.id.pos')) %>% 
    rename(trip.distance = value) %>% 
    # Calculate travel time by origin treso zone to hospital treso zone
    mutate(cases.travel.time = cases * trip.time * trips_per_case,
           cases.travel.distance = cases * trip.distance * trips_per_case)%>%
    # Find total and average travel time by CSD origin
    group_by(orig.csd, caretype) %>%
    mutate(orig.csd.caretype.travel.time.total = sum(cases.travel.time),
           orig.csd.caretype.travel.time.avg = orig.csd.caretype.travel.time.total / sum(cases) / trips_per_case,
           orig.csd.caretype.travel.distance.total = sum(cases.travel.distance),
           orig.csd.caretype.travel.distance.avg = orig.csd.caretype.travel.distance.total / sum(cases) / trips_per_case) %>%
    ungroup() %>% 
    # Calculate total and average travel times at the hospital/caretype level
    group_by(id, caretype) %>% 
    mutate(hosp.caretype.travel.time.total = sum(cases.travel.time),
           hosp.caretype.travel.time.avg = hosp.caretype.travel.time.total / sum(cases) / trips_per_case,
           hosp.caretype.travel.distance.total = sum(cases.travel.distance),
           hosp.caretype.travel.distance.avg = hosp.caretype.travel.distance.total / sum(cases) / trips_per_case) %>% 
    ungroup() %>% 
    # Calculate total and average travel times at the hospital level
    group_by(id) %>% 
    mutate(hosp.travel.time.total = sum(cases.travel.time),
           hosp.travel.time.avg = hosp.travel.time.total / sum(cases) / trips_per_case,
           hosp.travel.distance.total = sum(cases.travel.distance),
           hosp.travel.distance.avg = hosp.travel.distance.total / sum(cases) / trips_per_case) %>% 
    ungroup()
  
  return(list(crm_demand, crm_demand_travel))
}

format_hospital_asset_output <- function(num_beds, crm_demand, crm_demand_travel, dm) {
  # Summarise total and average travel times by hospital by caretype, convert to wide format
  hosp_travel_wide <- crm_demand_travel %>% 
    # Remove the temoporary hospital if DM is off
    filter(id != 99999) %>% 
    distinct(id, caretype,
             hosp.caretype.travel.time.total, hosp.caretype.travel.time.avg,
             hosp.caretype.travel.distance.total, hosp.caretype.travel.distance.avg) %>% 
    # Convert to wide format
    pivot_wider(names_from = c(caretype), 
                values_from = c(hosp.caretype.travel.time.total, hosp.caretype.travel.time.avg, 
                                hosp.caretype.travel.distance.total, hosp.caretype.travel.distance.avg)) %>% 
    left_join(distinct(crm_demand_travel, id, hosp.travel.time.total, hosp.travel.time.avg,
                       hosp.travel.distance.total, hosp.travel.distance.avg), by="id")
  
  # Convert num_beds to be in wide form (i.e., one row per asset)
  num_beds_wide <- num_beds %>%
    rename(beds.existing = staffedbeds) %>% 
    mutate(caretype.beds = paste0('beds.existing.', caretype)) %>%
    select(-year, -caretype, -count) %>% 
    pivot_wider(names_from = c(caretype.beds), values_from = c(beds.existing), values_fill = c(beds.existing = 0)) %>% 
    rename(beds.existing = total.hosp.beds)
  
  # Temporary dataframe to store the beds needed columns
  beds_needed <- crm_demand %>%   
    # Remove the temoporary hospital if DM is off
    filter(id != 99999) %>% 
    # Calculate beds needed and demand, summarised at the hospital level
    group_by(id, caretype) %>% 
    summarise(beds.needed = first(beds.needed)) %>% 
    ungroup() %>% 
    select(id, caretype, beds.needed) %>%
    # Pivot table from long to wide for caretype and agecluster for demand columns (cases, demand-days)
    pivot_wider(names_from = c(caretype), values_from = c(beds.needed), names_prefix = "beds.needed_") %>% 
    mutate_if(is.numeric, replace_na, replace=0)
  
  # Calculate demand and capacity data by asset (i.e., by hospital)
  hospital_wide <- crm_demand %>%   
    # Remove the temoporary hospital if DM is off
    filter(id != 99999) %>% 
    # Calculate beds needed and demand, summarised at the hospital level
    group_by(id, hosp.csd, caretype, agecluster) %>% 
    summarise(cases = sum(cases), demand.days.total = sum(demand.days.total), 
              demand.days.nonALC = sum(demand.days.nonALC), demand.days.LTC = sum(demand.days.LTC), 
              demand.days.ALC.ex.LTC = sum(demand.days.ALC.ex.LTC)) %>% 
    ungroup() %>% 
    select(id, csduid = hosp.csd, caretype, agecluster, cases, demand.days.total, demand.days.nonALC, demand.days.LTC, demand.days.ALC.ex.LTC) %>%
    # Pivot table from long to wide for caretype and agecluster for demand columns (cases, demand-days)
    pivot_wider(names_from = c(caretype, agecluster), 
                values_from = c(cases, demand.days.total, demand.days.nonALC, demand.days.LTC, demand.days.ALC.ex.LTC)) %>% 
    mutate_if(is.numeric, replace_na, replace=0) %>% 
    # Calculate sum of cases, demand-days as a double-check, and make sure there are exactly 0 beds for AM cases (Emergency)
    mutate(sum.cases = rowSums(select(., starts_with("cases_")))) %>%
    mutate(sum.demand.days.total_AM = rowSums(select(., starts_with("demand.days.total_AM"))),
           sum.demand.days.total_AT = rowSums(select(., starts_with("demand.days.total_AT"))),
           sum.demand.days.total_CR = rowSums(select(., starts_with("demand.days.total_CR"))),
           sum.demand.days.total_GR = rowSums(select(., starts_with("demand.days.total_GR"))),
           sum.demand.days.total_MH = rowSums(select(., starts_with("demand.days.total_MH"))),
           sum.demand.days = sum.demand.days.total_AM + sum.demand.days.total_AT + sum.demand.days.total_CR + 
             sum.demand.days.total_GR + sum.demand.days.total_MH) %>% 
    mutate(sum.LTC.days_AM = rowSums(select(., starts_with("demand.days.LTC_AM"))),
           sum.LTC.days_AT = rowSums(select(., starts_with("demand.days.LTC_AT"))),
           sum.LTC.days_CR = rowSums(select(., starts_with("demand.days.LTC_CR"))),
           sum.LTC.days_GR = rowSums(select(., starts_with("demand.days.LTC_GR"))),
           sum.LTC.days_MH = rowSums(select(., starts_with("demand.days.LTC_MH"))),
           sum.LTC.days = sum.LTC.days_AM + sum.LTC.days_AT + sum.LTC.days_CR + sum.LTC.days_GR + sum.LTC.days_MH) %>% 
    mutate(sum.ALC.ex.LTC.days_AM = rowSums(select(., starts_with("demand.days.ALC.ex.LTC_AM"))),
           sum.ALC.ex.LTC.days_AT = rowSums(select(., starts_with("demand.days.ALC.ex.LTC_AT"))),
           sum.ALC.ex.LTC.days_CR = rowSums(select(., starts_with("demand.days.ALC.ex.LTC_CR"))),
           sum.ALC.ex.LTC.days_GR = rowSums(select(., starts_with("demand.days.ALC.ex.LTC_GR"))),
           sum.ALC.ex.LTC.days_MH = rowSums(select(., starts_with("demand.days.ALC.ex.LTC_MH"))),
           sum.ALC.ex.LTC.days = sum.ALC.ex.LTC.days_AM + sum.ALC.ex.LTC.days_AT + sum.ALC.ex.LTC.days_CR + sum.ALC.ex.LTC.days_GR + sum.ALC.ex.LTC.days_MH) %>% 
    left_join(beds_needed, by="id") %>% 
    mutate(beds.needed_AM = 0) %>% 
    mutate(total.hosp.beds.needed = beds.needed_AM + beds.needed_AT + beds.needed_CR + beds.needed_GR + beds.needed_MH) 
  
  # Join in capacity data to hospital_wide
  hospital_wide_cap <- hospital_wide %>% 
    # Join in existing capacity (i.e., beds.existing) data
    left_join(num_beds_wide, by = c('id')) %>% 
    replace_na(list(total.hosp.beds.existing = 0, beds.existing = 0, beds.existing.AM = 0, beds.existing.AT = 0, beds.existing.CR = 0, 
                    beds.existing.GR = 0, beds.existing.MH = 0)) %>% 
    # Add one bed to hospital/caretype pairs which have historical demand but 0 beds shown
    mutate(beds.existing.AT = ifelse((sum.demand.days.total_AT > 0) & (beds.existing.AT == 0), 1, beds.existing.AT),
           beds.existing.CR = ifelse((sum.demand.days.total_CR > 0) & (beds.existing.CR == 0), 1, beds.existing.CR),
           beds.existing.GR = ifelse((sum.demand.days.total_GR > 0) & (beds.existing.GR == 0), 1, beds.existing.GR),
           beds.existing.MH = ifelse((sum.demand.days.total_MH > 0) & (beds.existing.MH == 0), 1, beds.existing.MH),
           beds.existing = beds.existing.AT + beds.existing.CR + beds.existing.GR + beds.existing.MH) %>% 
    # Replace beds.needed columns (by caretype) with the pmax of (beds.needed, beds.existing)
    mutate(dm = dm,
           beds.forecasted.AT = ifelse(dm, pmax(beds.needed_AT, beds.existing.AT), beds.existing.AT),
           beds.forecasted.CR = ifelse(dm, pmax(beds.needed_CR, beds.existing.CR), beds.existing.CR),
           beds.forecasted.GR = ifelse(dm, pmax(beds.needed_GR, beds.existing.GR), beds.existing.GR),
           beds.forecasted.MH = ifelse(dm, pmax(beds.needed_MH, beds.existing.MH), beds.existing.MH),
           beds.forecasted = beds.forecasted.AT + beds.forecasted.CR + beds.forecasted.GR + beds.forecasted.MH) %>% 
    select(-dm)
    
  # Calculate utilization by hospital, caretype for forecasted beds, and utilizations using existing beds
  hospital_wide_util <- hospital_wide_cap %>% 
    # Utilization calculations use the sum of demand-days rather than beds.needed because beds.needed 
    # incorporates user-input utilization targets 
    mutate(utilization.AM = 0,
           utilization.AT = sum.demand.days.total_AT / beds.forecasted.AT / HOSP_OP_DAYS,
           utilization.CR = sum.demand.days.total_CR / beds.forecasted.CR / HOSP_OP_DAYS,
           utilization.GR = sum.demand.days.total_GR / beds.forecasted.GR / HOSP_OP_DAYS,
           utilization.MH = sum.demand.days.total_MH / beds.forecasted.MH / HOSP_OP_DAYS,
           utilization = sum.demand.days / beds.forecasted / HOSP_OP_DAYS,
           utilization.AT.existing = sum.demand.days.total_AM / beds.existing.AT / HOSP_OP_DAYS,
           utilization.CR.existing = sum.demand.days.total_CR / beds.existing.CR / HOSP_OP_DAYS,
           utilization.GR.existing = sum.demand.days.total_GR / beds.existing.GR / HOSP_OP_DAYS,
           utilization.MH.existing = sum.demand.days.total_MH / beds.existing.MH / HOSP_OP_DAYS,
           utilization.existing = sum.demand.days / beds.existing / HOSP_OP_DAYS) %>%
    # Replace 'Inf' utilizations into NA
    mutate_at(vars(utilization.existing, utilization.AT.existing, utilization.CR.existing, utilization.GR.existing, utilization.MH.existing), ~replace(., is.infinite(.), NA))
  
  # Add in travel times, distances to final hospital_wide table
  hospital_wide_final <- hospital_wide_util %>% 
    left_join(hosp_travel_wide, by = c('id'))
  
  return(hospital_wide_final)
  
}

# MAG -----
## MAG Cleaning ----

create_court_xy <- function(court_master) {
  #' This function transforms each court location into the same projection as the TRESO shapefiles
  #' 
  #' @param court_master A Dataframe of the courts
  #' @param treso_shp A Shapefile of TRESO zones
  #' @return A SpatialPointsDataframe 
  
  court_spdf <- court_master %>%
    ungroup() %>%
    select(bid, courthouse.lat, courthouse.long) %>%
    mutate(lat = courthouse.lat, long = courthouse.long) %>%
    mutate(id = row_number())
  coordinates(court_spdf) <- c('long', 'lat')
  
  # Project the student and school points from lat/long to TRESO's LCC specification
  proj4string(court_spdf) <- CRS('+proj=longlat +datum=WGS84')
  
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
  
  courtroom_count_scaled = pmax(pmin(courtrooms, 50), 1) # Setting area standard to have a minimum courtroom count of 1 and a max of 50 in linear sizing scale
  required_area_per_courtroom = area_max_sqm - courtroom_count_scaled * area_change_per_courtroom
  
  required_area_per_courtroom_sqft = required_area_per_courtroom * 3.28^2
  
  return(required_area_per_courtroom_sqft)
}

## MAG Model ----
prepare_treso_population <- function(treso_population, treso_zone_system) {
  # Gather the TRESO Population table into long format to prepare data for business line, year, and TRESO Zone grouping
  treso_population_gather <- treso_population %>% 
    select(-csduid) %>% 
    rename_at(vars(contains("population")), list(~substr(., 12, 15))) %>% 
    gather(key="year", value = "population", -treso_zone, -age.group)
  
  # Transform age ranges to business line populations by zone and projection year to sum zone population by matched business line age ranges
  treso_population_bl <- treso_population_gather %>% 
    mutate(OCJcrim = ifelse(age.group %in% c("12-17", "18", "19", "20-49"), population, 0)) %>% 
    mutate(SCJcrim = ifelse(age.group %in% c("18", "19", "20-49"), population, 0)) %>% 
    mutate(OCJfam = ifelse(age.group %in% c("0-11", "12-17"), population, 0)) %>% 
    mutate(SCJfam = ifelse(age.group %in% c("0-11", "12-17", "18", "19", "20-49", "50-59"), population, 0)) %>% 
    mutate(civil = ifelse(age.group %in% c("20-49", "50-59", "60-69"), population, 0)) %>%
    mutate(smac = ifelse(age.group %in% c("20-49", "50-59", "60-69"), population, 0)) %>%
    group_by(treso_zone, year) %>% 
    summarise_if(is.numeric, sum) %>% 
    left_join(select(treso_zone_system, treso_id, cduid, cdname), by = c("treso_zone" = "treso_id"))
  
  # Output the population by business line as long format for historical court hour matching
  treso_population_bl_long <- treso_population_bl %>%
    gather(casetypes, population, OCJcrim:smac) %>%
    arrange(treso_zone, casetypes, year) %>% 
    ungroup()
  
  return(treso_population_bl_long)
}

combine_hours_and_population <- function(treso_population_bl_long, court_hours) {
  # Summarize to CD level for population due to court hours being at the base court CD level
  treso_population_bl_cd <- treso_population_bl_long %>% 
    group_by(year, cduid, cdname, casetypes) %>%
    summarise(population = sum(population)) %>% 
    ungroup() %>% 
    mutate(cdname = as.character(cdname))
  
  # Combine CD populations with court hours
  court_hours_population <- court_hours %>% 
    left_join(treso_population_bl_cd, by = c('year', 'casetypes', 'cduid', 'cdname')) %>% 
    mutate_at("year", as.numeric) %>% 
    rename(hours = value, treso_id = treso.id.pos, courthouse.year.total.hours = total) %>% 
    select(-lat, -long, -csduid, -csdname, - courthouse.lat, - courthouse.long)
  
  # CD 3546 (Haliburton) is missing court hours - requested from MAG but not received as of Dec 2019
  # CD 3553 (Greater Sudbury) is within CD 3552 (Sudbury)
  # The courts are located in Greater Sudbury, so the two CDs are combined and treated as one CD
  sudbury_population_bl_cd <- treso_population_bl_cd %>% 
    mutate_at("year", as.numeric) %>% 
    filter(cduid == 3552) %>% 
    mutate(cduid = 3553) %>% 
    rename(population.sudbury = population) %>% 
    select(year, cduid, casetypes, population.sudbury)
  
  # Combine Subury population/hours with remainder of populations/hours  
  court_hours_population <- court_hours_population %>% 
    left_join(sudbury_population_bl_cd, by = c("cduid", "year", "casetypes")) %>% 
    replace_na(list(population.sudbury = 0)) %>% 
    mutate(population = population + population.sudbury)
  
  cat(paste0("Quick check that the population in", "\n",
             "CD 3553 and 3552 is: ", 
             sum(filter(treso_population_bl_cd, year == 2018, cduid %in% c(3552, 3553))$population), "\n",
             "Combined is: ",
             (filter(court_hours_population, year == 2018, cduid == 3553) %>% group_by(bid, cduid) %>% summarise(population = sum(population)))$population[1]))
  
  return(court_hours_population)
}

calculate_baseyear_rates <- function(court_hours_population, weight_rate) {
  # Determine rate of court hours per person (relevant population) for each courthouse; weighting determined as set by user in inputs
  court_rate <- court_hours_population %>% 
    left_join(WEIGHT_RATE, by = "year") %>% 
    filter(!is.na(weight)) %>%
    mutate(pop.rate = hours / population) %>% 
    # Average the years, as weighted by user
    group_by(bid, name, cduid, cdname, casetypes) %>% 
    summarise(pop.rate = weighted.mean(pop.rate, weight),
              hours = sum(hours)) %>% 
    ungroup()
  
  # Add courthouses together, grouped by CD, to get total CD rates, since relevant population for courthoues are equal (by casetype) within any CD 
  baseyear_rates <- court_rate %>% 
    group_by(cduid, cdname, casetypes) %>% 
    summarise(pop.rate = sum(pop.rate)) %>% 
    replace_na(list(pop.rate = 0)) %>% 
    ungroup()
  
  return(baseyear_rates)
}

calculate_forecast_cases <- function(baseyear_rates, treso_population_bl_long) {
  # Summarize to CD level for population due to court hours being at the base court CD level
  treso_population_bl_cd <- treso_population_bl_long %>% 
    group_by(year, cduid, cdname, casetypes) %>%
    summarise(population = sum(population)) %>% 
    ungroup() %>% 
    mutate(cdname = as.character(cdname))
  
  future_cases_bl <- baseyear_rates %>%
    left_join(treso_population_bl_cd, by = c("cduid", "cdname", "casetypes")) %>%
    spread(key = "year", value = "population") %>%
    # Remove historical data
    select(-`2011`:-`2017`) %>%
    # Multiply populations by hours rate
    mutate_at(vars(starts_with("20")), list(~(. * pop.rate))) %>% 
    gather(year, court.hours, `2018`:`2041`)
  
  # Calculate the cases at TRESO level
  # CD 3552 (Sudbury) will share rates from CD 3553 (Greater Sudbury)
  # The case rates were calculated using the combined population
  sudbury_rates <- baseyear_rates %>% 
    filter(cduid == 3553) %>% 
    select(cduid, casetypes, pop.rate) %>% 
    mutate(cduid = 3552) %>% 
    rename(pop.rate.sudbury = pop.rate)
  
  # Project future cases at the TRESO level based on CD rates and TRESO populations
  future_cases_bl_treso <- treso_population_bl_long %>% 
    mutate(cdname = as.character(cdname)) %>% 
    left_join(baseyear_rates, by = c("cduid", "cdname", "casetypes")) %>% 
    left_join(sudbury_rates, by = c("cduid", "casetypes")) %>% 
    mutate(pop.rate = ifelse(!is.na(pop.rate.sudbury), pop.rate.sudbury, pop.rate)) %>% 
    select(-pop.rate.sudbury) %>% 
    mutate(court.hours = population * pop.rate) %>% 
    # Replace N/A values where no hours are available - this should not be necessary once data for Haliburton is obtained
    replace_na(list(court.hours = 0, pop.rate = 0)) 
  
  return(future_cases_bl_treso)
}

calculate_courtrooms_required <- function(future_cases_bl_treso, utilization_rate, courthouse_asset, op_days, year_id) {
  
  # Calculate future court hours needed based on the business-line-specific utilization characteristics
  future_needs <- future_cases_bl_treso %>%
    filter(year == year_id) %>% 
    left_join(utilization_rate, by = "casetypes") %>% 
    mutate(op.utilization = 0) %>% 
    # If courtroom hours are <= inflection hours, use (inflection hours - 0) to calculate where on the capacity slope each courthouse falls
    mutate(op.utilization = base.hours.dy + (inflect.hours.dy - base.hours.dy) / (inflect.hours - 0 ) * court.hours) %>% 
    # If courtroom hours are > inflection hours, replace calculation above with (max.hours - inflection hours) to calculate where on capacity slope courthouses fall
    mutate(op.utilization = ifelse(court.hours > inflect.hours & court.hours <= max.hours,
                                   inflect.hours.dy + ((max.hours.dy - inflect.hours.dy) * (court.hours - inflect.hours) / (max.hours - inflect.hours)),
                                   op.utilization)) %>% 
    mutate(op.utilization = ifelse(court.hours > max.hours,
                                   max.hours.dy,
                                   op.utilization)) %>% 
    mutate(courtrooms.needed = round(court.hours / (op.utilization * op_days), digits = 2))
  
  # Summarise to CD level
  future_needs_cd <- future_needs %>% 
    group_by(cduid, cdname) %>% 
    summarise(courthours.needed = sum(court.hours),
              courtrooms.needed = sum(courtrooms.needed)) 
  
  # Summarise courtrooms needed with the existing 2018 courtrooms dataset, including GFA required
  courthouse_asset_with_projected_cd <- courthouse_asset %>%
    group_by(cduid, cdname) %>% 
    summarise(number.of.courtrooms.2018 = sum(courtrooms),
              rentable.square.feet.2018 = sum(rentable.square.feet)) %>% 
    left_join(future_needs_cd, by = c("cduid", "cdname")) %>% 
    mutate(rentable.square.feet.needed = courtrooms.needed * courtroom_size(courtrooms.needed)) %>% 
    ungroup()
  
  return(courthouse_asset_with_projected_cd)
  
}

distribute_mag_demand <- function(future_cases_bl_treso, courthouse_asset, new_courthouse, 
                                  treso_shp, travel_time_skim, travel_distance_skim,
                                  modernization_factor, appearance_factor, distribution_percentage, u_avg, cost_per_sqft,
                                  cost_per_courthouse, op_days, year_id) {
  # Prepare the origin demand vector for selected year
  future_cases_bl_treso <- future_cases_bl_treso %>%
    rename(treso.id.por = treso_zone) %>% 
    filter(year == year_id) %>%
    select(-year, -population, -pop.rate) %>%
    left_join(modernization_factor, by = "casetypes") %>%
    mutate(court.hours = court.hours * mod.factor)
  
  # Overlay with TRESO Shapefile to get the TRESO/CD Information
  if (!is.null(new_courthouse)) {
    new_courthouse_xy <- create_court_xy(new_courthouse)
    new_courthouse_overlay <- create_overlay(new_courthouse_xy, treso_shp, type = "court") %>% 
      left_join(treso_zone_system, by = c("treso.id.pos" = "treso_id")) %>% 
      left_join(new_courthouse, by = "bid") %>% 
      mutate_if(is.factor, as.character)
    
    # Find out the treso zone and CD region
    courthouse_asset <- bind_rows(courthouse_asset, new_courthouse_overlay)
  }
  
  # Create a list of courthouses available in each CD
  courthouse_available <- courthouse_asset %>%
    unnest(casetypes = strsplit(casetypes.serviced, ",")) %>% 
    group_by(cduid, cdname, casetypes) %>% 
    summarise(treso.id.pos = paste(treso.id.pos, collapse = ","),
              bid = paste(bid, collapse = ","),
              name = paste(name, collapse = ","),
              building.type = paste(building.type, collapse = ","),
              courtrooms = paste(courtrooms, collapse = ",")
    )
  
  # Join the list of courthouses available by CD to the origin demand
  # NOTE: Sudbury (CD 3552) needs to be matched with Greater Sudbury (CD 3553) courts
  # Distribute a percentage of TRESO demand to closest courthouse within the CD
  # Distribute another percentage of TRESO demand weight by available courhouses within the CD
  courthours_od <- future_cases_bl_treso %>%
    filter(court.hours != 0) %>% 
    mutate(cduid = ifelse(cduid == 3552, 3553, cduid),
           cdname = ifelse(cdname == "Sudbury", "Greater Sudbury / Grand Sudbury", cdname)) %>% 
    left_join(courthouse_available, by = c("cduid", "cdname", "casetypes")) %>% 
    unnest(treso.id.pos = strsplit(treso.id.pos, ","),
           bid = strsplit(bid, ","),
           name = strsplit(name, ","),
           building.type = strsplit(building.type, ","),
           courtrooms = strsplit(courtrooms, ",")) %>% 
    mutate(treso.id.por = as.numeric(treso.id.por),
           treso.id.pos = as.numeric(treso.id.pos),
           courtrooms = as.numeric(courtrooms)) %>% 
    # Merge with TRESO travel times
    left_join(travel_time_skim, by = c("treso.id.por", "treso.id.pos")) %>%
    rename(travel.time = value) %>% 
    left_join(travel_distance_skim, by=c("treso.id.por", "treso.id.pos")) %>% 
    rename(travel.distance = value) %>% 
    # Identify the closest travel time pair by calculating a travel utility (e^(1/tt))
    mutate(travel.utility = exp(1/travel.time)) %>%
    arrange(treso.id.por, casetypes, desc(travel.utility)) %>% 
    group_by(treso.id.por, casetypes) %>% 
    # Give a utility ranking for each combination of POR and casetypes
    mutate(travel.utility.ranking = row_number(),
           courtroom.weight = courtrooms / sum(courtrooms)) %>% 
    ungroup() %>%
    # Flag the courthouses in the same TRESO zone servicing the same casetype that has the closest travel utility
    group_by(treso.id.por, casetypes, travel.utility) %>% 
    mutate(identical = ifelse(nchar(paste(unique(bid), collapse = ",")) > 6 && travel.utility.ranking == 1,
                              ceiling(nchar(paste(unique(bid), collapse = "")) / 6), 0)) %>% 
    ungroup() %>% 
    # Join in closest travel and random assignment scaling factors
    left_join(distribution_percentage, by = c('casetypes')) %>% 
    # If the courthouse is the closest but also have a identical flag - change the travel utility ranking to 2
    group_by(treso.id.por, casetypes, identical) %>% 
    mutate(travel.utility.ranking = ifelse(identical > 0, 2, travel.utility.ranking)) %>% 
    mutate(court.hours.to.closest = ifelse(travel.utility.ranking == 1, court.hours * closest.percentage, 0),
           court.hours.to.random = court.hours * random.percentage * courtroom.weight,
           court.hours.to.identical = ifelse(identical > 0, court.hours * closest.percentage * courtroom.weight / sum(courtroom.weight), 0)) %>% 
    ungroup()
  
  # Summarize the distributed demand to each courthouses and calculate the utilization for each courthouse
  courthouse_asset_updated <- courthours_od %>% 
    mutate(court.hours.total = court.hours.to.closest + court.hours.to.random + court.hours.to.identical) %>% 
    # Add in factors for converting court-hours to cases, appearances, and trips
    left_join(appearance_factor, by = "casetypes") %>% 
    mutate(trips.per.court.hour = case.factor * app.factor * trip.factor,
           trips = trips.per.court.hour * court.hours.total) %>% 
    # Calculate travel times, distances
    group_by(bid, name, building.type, courtrooms, cduid, cdname, casetypes) %>% 
    summarise(court.hours.to.closest = sum(court.hours.to.closest),
              court.hours.to.random = sum(court.hours.to.random),
              court.hours.to.identical = sum(court.hours.to.identical),
              travel.time.avg = weighted.mean(travel.time, trips),
              travel.distance.avg = weighted.mean(travel.distance, trips),
              travel.time = sum(travel.time * trips),
              travel.distance = sum(travel.distance * trips),
              trips = sum(trips)) %>% 
    ungroup() %>% 
    # Join in rentable square feet to determine area shortage
    left_join(select(courthouse_asset, bid, treso.id.pos, courthouse.lat, courthouse.long, rentable.square.feet), by = c('bid')) %>% 
    rename(existing.rentable.sqft = rentable.square.feet) %>% 
    replace_na(list(existing.rentable.sqft = 0)) %>% 
    # Join in target utilizations (or neutral targets) based on decision-making choice
    left_join(u_avg, by = c('casetypes')) %>% 
    mutate(court.hours.distributed.prescaled = court.hours.to.closest + court.hours.to.random + court.hours.to.identical,
           court.hours.distributed = court.hours.distributed.prescaled / target.utilization) %>% 
    # Join in daily utilizable hours per courtroom and calculate utilization by courthouse
    left_join(utilization_rate, by = "casetypes") %>% 
    mutate(op.utilization = 0) %>% 
    mutate(op.utilization = ifelse(court.hours.distributed <= inflect.hours,
                                   base.hours.dy + (inflect.hours.dy - base.hours.dy) / (inflect.hours - 0 ) * court.hours.distributed,
                                   op.utilization)) %>% 
    mutate(op.utilization = ifelse(court.hours.distributed > inflect.hours & court.hours.distributed <= max.hours,
                                   inflect.hours.dy + ((max.hours.dy - inflect.hours.dy) * (court.hours.distributed - inflect.hours)                                                       / (max.hours - inflect.hours)),
                                   op.utilization)) %>% 
    mutate(op.utilization = ifelse(court.hours.distributed > max.hours,
                                   max.hours.dy,
                                   op.utilization)) %>% 
    # Summarise courtrooms needed, courtrooms existing, utilization, etc.
    mutate(courtrooms.needed = round(court.hours.distributed / (op.utilization * op_days), digits = 2)) %>% 
    mutate(utilization = court.hours.distributed.prescaled / courtrooms.needed / op.utilization / op_days) %>% 
    group_by(bid, treso.id.pos, courthouse.lat, courthouse.long, cduid, cdname, name, building.type) %>% 
    summarise(existing.rentable.sqft = first(existing.rentable.sqft),
              travel.time = sum(travel.time),
              travel.time.avg = weighted.mean(travel.time.avg, trips),
              travel.distance = sum(travel.distance),
              travel.distance.avg = weighted.mean(travel.distance.avg, trips),
              courtrooms = first(courtrooms),
              courtrooms.needed = sum(courtrooms.needed),
              utilization = weighted.mean(utilization, court.hours.distributed.prescaled),
              court.hours.distributed.prescaled = sum(court.hours.distributed.prescaled),
              trips = sum(trips)) %>% 
    ungroup()
  
  if(!is.null(new_courthouse)) {
    courthouse_asset_updated <- courthouse_asset_updated %>%
      # Join in list of new courthouses input by user to flag new courthouses
      left_join(select(new_courthouse, bid, new.courthouse), by = c('bid')) %>% 
      replace_na(list(new.courthouse = 0))
  } else {
    courthouse_asset_updated <- courthouse_asset_updated %>% 
      mutate(new.courthouse = 0)
  }
  
  courthouse_asset_updated <- courthouse_asset_updated %>% 
    # Calculate the number of courtrooms to use for area calculations; use new.courthouse flag to calculate area for ALL new courtrooms, not just courtrooms.needed, at new courthouses
    mutate(courtrooms.for.std.area.calc = courtrooms.needed + new.courthouse * (courtrooms - courtrooms.needed),
           net.new.courtrooms = pmax(courtrooms.needed - courtrooms, 0) + new.courthouse * courtrooms) %>% 
    # Calculate standard areas for courtrooms and area shortfalls
    mutate(standard.sqft.per.courtroom = courtroom_size(courtrooms.for.std.area.calc),
           built.area = existing.rentable.sqft + standard.sqft.per.courtroom * (net.new.courtrooms),
           area.shortfall = standard.sqft.per.courtroom * courtrooms.for.std.area.calc - built.area,
           courtrooms.shortfall = courtrooms.needed - courtrooms) %>% 
    # Select final display fields
    select(bid, treso.id.pos, cdname, cduid, courthouse.lat, courthouse.long, 
           courtrooms, courtrooms.needed, new.courthouse, courtrooms.shortfall, court.hours.distributed.prescaled, 
           utilization, existing.rentable.sqft, standard.sqft.per.courtroom, area.shortfall, built.area, trips,
           travel.time, travel.time.avg, travel.distance, travel.distance.avg)
  
  # Calculate 2018 costs
  courthouse_asset_updated <- courthouse_asset_updated %>% 
    mutate(cost_2018 = round((built.area - existing.rentable.sqft) * cost_per_sqft + new.courthouse * cost_per_courthouse, 0))
  
  return(list(courthouse_asset_updated, courthours_od))
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

utilizationColorLegend <- data.frame(x = c("green", "orange", "red"), y = c("<= 0.75", "0.75 to 1.0", ">= 1.0"))


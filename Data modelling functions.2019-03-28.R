
##### Convert degrees to radians
convertDegrees <- function(deg) {
  rad <- deg * pi / 180
  return(rad)
}

##### Haversine function

haverFunction <- function(lat1, lon1, lat2, lon2){
  "
  Calculate the straightline distance between a pair of points.

  Inputs: Point1 lat, Point1 lon, Point2 lat, Point2 lon
  Output: Straightline distance in km between two points.

  #haversine = sin^2(theta/2)
  #d = [2r][arcsin(sqrt(hav(lat2 - lat1)+cos(lat1)*cos(lat2)*hav(long2 - long1))]
  #d = [2r][arcsin(sqrt(sin^2((lat2 - lat1)/2)+cos(lat1)*cos(lat2)*sin^2((long2 - long1)/2)))]
  "

  earthRadius <- 6378.1 #in km
  lat1Rad <- convertDegrees(lat1)
  lat2Rad <- convertDegrees(lat2)
  long1Rad <- convertDegrees(long1)
  long2Rad <- convertDegrees(long2)
  dLat <- lat2Rad - lat1Rad
  dLong <- long2Rad - long1Rad
    
  
  bracketCalc <- sin(dLat/2)^2 + cos(lat1Rad) * cos(lat2Rad) * (sin(dLong/2))^2
  d <- 2 * earthRadius * asin(bracketCalc^0.5)
  return(d)
}

##### Gravity model function

furness <- function(cmat){
  
  "
  This function conducts furnessing or two dimensional balancing
  
  Inputs: observed truck matrix
  
  Arguments: cost function
  
  Results: balanced matrix
  "
  
  #' row sum of cost function matrix
  skim_trow <- as.data.frame(rowSums(cmat))
  #' row sum of observed matrix
  obs_tatt <- as.data.frame(rowSums(proflinkages_zero_cars1))
  
  #' row balancing factors
  rowbal <- obs_tatt/skim_trow
  colnames(rowbal) <- "ratio"
  # rowbal1 <- do.call("cbind", replicate(ncol(proflinkages_zero_cars1), rowbal, 
  #                                       simplify = FALSE))
  
  #first iteration
  cfunc1 <- cmat*rowbal$ratio
  
  #' col sum of cost function matrix
  skim_tcol <- as.data.frame(colSums(cfunc1))
  #' col sum of observed matrix
  obs_tprod <- as.data.frame(colSums(proflinkages_zero_cars1))
  
  #' column balancing factors
  colbal <- obs_tprod/skim_tcol
  colnames(colbal) <- "ratio"
  # colbal1 <- do.call("rbind", replicate(nrow(proflinkages_zero_cars1), colbal, 
  #                                       simplify = FALSE))
  
  # next iteration. Transpose the matrix to allow multiplying by the column ratios.
  # once done retranspose the matrix before running the row balancing
  cfunc1 <- t(cfunc1)*colbal$ratio
  cfunc1 <- t(cfunc1)
  
  return(cfunc1)
} 

##### Importing Schools Data into SQL Server

#Loading school postal code data from CSV to SQL Server
#schoolPost <- read.csv('C:/Users/azk/Documents/Working away files/EDU_enrol_SCH_PCode_WSP_Temp.csv')

#schoolPostSlim <- schoolPost %>%
#  mutate(rowNum = 1:n()) %>%
#  select(rowNum, Sch_YR, BoardNumber, dsbindex, BSID, SchoolName=School.Name, PostalCode=Student.Postal.Code, DataID, TableID)

#schoolPost1 <- schoolPostSlim %>%
#filter(rowNum %in% (1010372:1999999))

# schoolPost2 <- schoolPostSlim %>%
#   filter(rowNum %in% (2000000:2999999))
# 
# schoolPost3 <- schoolPostSlim %>%
#   filter(rowNum %in% (2348080:2999999))
# 
# schoolPost4 <- schoolPostSlim %>%
#   filter(rowNum %in% (3000000:333384299))
# 
# schoolPost5 <- schoolPostSlim %>%
#   filter(rowNum %in% (3338430:3999999))
# 
# schoolPost6 <- schoolPostSlim %>%
#   filter(rowNum %in% (4000000:4999999))
# 
# schoolPost7 <- schoolPostSlim %>%
#   filter(rowNum %in% (5000000:5607520))
# 
# schoolPost8 <- schoolPostSlim %>%
#   filter(rowNum %in% (5607521:6702802))
# 
# schoolPost9 <- schoolPostSlim %>%
#   filter(rowNum %in% (6702802:7025441))

# schoolPost10 <- schoolPostSlim %>%
#   filter(rowNum > 7025441)

#sqlSave(schoolRodbc,schoolPost10,"cap.EDU_enrol_SCH_PCodeData",safer=FALSE, fast=TRUE)
#sqlUpdate(schoolRodbc,schoolPostSlim,"cap.EDU_enrol_SCH_PCodeData",safer=FALSE, fast=TRUE, index=c(rownames))


  
#' Landform Classification
#'
#' This function classifies a landscape into different categories based on the Topographic Position Index (TPI) and slope.
#' It offers two types of classifications: Slope Position Classification (Weiss 2001) and Landform Classification (Weiss 2001; Jenness 2003).\cr
#' Visit this \href{https://drive.google.com/file/d/1vjy8_HPtuyKkzSkYwKdQbHJLrVtO79ft/view?usp=sharing}{LINK} to access the package's vignette.\cr
#'
#' @param x A RasterLayer or a SpatRaster object representing the landscape.
#' @param scale The size of the neighbourhood for calculating the TPI. Default is 3.
#' @param sn The size of the small neighbourhood for calculating the TPI in the Landform Classification. Default is 3.
#' @param ln The size of the large neighbourhood for calculating the TPI in the Landform Classification. Default is 7.
#' @param class.type The type of classification to be performed. Options are "slope.position" for Slope Position Classification
#' and "landform.classification" for Landform Classification. Default is "slope.position".
#' @param descriptive A logical parameter (default is FALSE). If set to TRUE, the function will calculate and return
#' additional descriptive statistics for each class. These statistics include the count of pixels, total area (in Km2), and percentage
#' of the total area for each class. Additionally, a bar plot showing the area of each class will be generated.
#' @param leg.pos The position of the legend in the plot. Default is "topright".
#' @param leg.cex The magnification to be used for sizing the legend relative to the current setting of 'cex'. Default is 0.65.
#'
#' @return The package plots the raster of all the classes combined. If the \code{descriptive} parameter is set to \code{TRUE},
#' a bar plot showing the area of each class is generated. A list of SpatRaster objects
#' representing the different classes of the landscape is also returned. The list will contain an item named
#' \code{Descriptive Statistics} if the \code{descriptive} parameter is set to \code{TRUE}. For each class, the item reports the
#' cell counts, the corresponding area (in Km2), and the corresponding percentage.\cr
#'
#' For the "slope.position" classification, the returned layers are:
#' \itemize{
#'   \item "all": all classes combined
#'   \item "Valley": TPI <= -1
#'   \item "Lower Slope": -1 < TPI <= -0.5
#'   \item "Flat Slope": -0.5 < TPI < 0.5 and slope <= 5
#'   \item "Middle Slope": -0.5 < TPI < 0.5 and slope > 5
#'   \item "Upper Slope": 0.5 < TPI <= 1
#'   \item "Ridge": TPI > 1
#' }
#'
#' For the "landform.classification", the returned layers are:
#' \itemize{
#'   \item "all": all classes combined
#'   \item "Canyons-Deeply Incised Streams": sn <= -1 and ln <= -1
#'   \item "Midslope Drainage-Shallow Valleys": sn <= -1 and -1 < ln < 1
#'   \item "Upland Drainages-Headwaters": sn <= -1 and ln >= 1
#'   \item "U-shaped Valleys": -1 < sn < 1 and ln <= -1
#'   \item "Plains": -1 < sn < 1 and -1 < ln < 1 and slope <= 5
#'   \item "Open Slopes": -1 < sn < 1 and -1 < ln < 1 and slope > 5
#'   \item "Upper Slopes-Mesas": -1 < sn < 1 and ln >= 1
#'   \item "Local Ridges-Hills in Valley": sn >= 1 and ln <= -1
#'   \item "Midslopes Ridges-Small Hills in Plains": sn >= 1 and -1 < ln < 1
#'   \item "Mountain Tops-High Ridges": sn >= 1 and ln >= 1
#' }
#'
#' @details
#' The function internally calculates the Standardised Topographic Position Index (TPI) for the given landscape.
#' The TPI is calculated as the difference between the elevation of each cell in the landscape and the mean elevation of
#' its neighbouring cells. The TPI is then standardised by dividing the difference between each cell's TPI and the mean TPI
#' of its neighbouring cells by the standard deviation of the TPI values within the neighbourhood.
#'
#' The Slope Position Classification uses the standardised TPI and slope to classify the landscape into
#' six categories: Valley, Lower Slope, Flat Slope, Middle Slope, Upper Slope, and Ridge.
#' This classification is useful for identifying the position of a location on a slope (Weiss 2001).
#'
#' The Landform Classification uses two standardised TPI grids at different scales to classify the landscape
#' into ten categories: Canyons, Midslope Drainage, Upland Drainage, U-shaped Valleys, Plains, Open Slopes,
#' Upper Slopes, Local Ridges, Midslope Ridges, and Mountain Tops. This classification is useful for identifying
#' broader landform types (Weiss 2001; Jenness 2003).
#'
#' As for the descriptive statistics returned by the function, please note that if the input raster layer is not in
#' a projected coordinate system, the area calculation for each class (which is based on the resolution of the raster and the
#' count of cells in each class) may not accurately represent the true area on the ground. Therefore, it's recommended to use
#' a raster layer in a suitable projected coordinate system for your area of interest to ensure accurate area calculations.
#'
#' @references
#' Weiss, A. (2001). Topographic Position and Landforms Analysis. Poster presentation, ESRI User Conference, San Diego, CA.
#'
#' Jenness, J. (2003). TPI ArcView Extension. Jenness Enterprises. Available at: http://www.jennessent.com/arcview/tpi.htm
#'
#' @seealso
#' \code{\link[raster]{raster}}, \code{\link[terra]{rast}}, \code{\link[terra]{focal}}, \code{\link[terra]{terrain}}, \code{\link[terra]{plot}}
#'
#' @examples
#'
#' # Create a toy elevation dataset
#' # Define the raster dimensions
#' nrows <- 100
#' ncols <- 100
#'
#' # Create a matrix with random values
#' set.seed(123) # to make the example reproducible
#' mat <- matrix(runif(nrows*ncols), nrow=nrows, ncol=ncols)
#'
#' # Convert the matrix to a raster object
#' library(terra)
#' elev <- rast(mat)
#'
#' # EXAMPLE 1
#'
#' # Run the analysis
#' result <- landform(elev, scale = 3, class.type = "slope.position")
#'
#' # EXAMPLE 2
#'
#' # Run the analysis
#' result <- landform(elev, class.type = "landform.classification", sn=3, ln=11)
#'
#'
#' @importFrom terra rast values terrain plot
#' @importFrom graphics legend barplot
#' @importFrom stats sd
#' @import grDevices
#'
#' @export
#'
#'
landform <- function (x, scale = 3, sn=3, ln=7, class.type="slope.position", descriptive = FALSE, leg.pos="topright", leg.cex=0.65) {

  # Check if x is a Raster object and, if it is, convert to a SpatRaster object
  if("RasterLayer" %in% class(x)) {
    x <- terra::rast(x)
  }

  # Function to calculate and standardise TPI
  calculate_standardised_tpi <- function(raster_layer, neighbourhood_size) {

    # Define the matrix for the moving window
    w <- matrix(1, nrow=neighbourhood_size, ncol=neighbourhood_size)

    # Calculate mean elevation in the neighbourhood
    mean_elevation <- terra::focal(raster_layer, w, fun=mean, pad=TRUE, na.rm=TRUE)

    # Calculate TPI
    tpi <- raster_layer - mean_elevation

    # Calculate standard deviation of the neighbourhood values
    sd_neighbourhood <- terra::focal(raster_layer, w, fun=sd, pad=TRUE, na.rm=TRUE)

    # Standardise TPI
    tpi_standardised <- tpi / sd_neighbourhood

    return(list("tpi_standardised" = tpi_standardised))
  }

  # Put the function to work and calculate the standardised TPI
  tp <- calculate_standardised_tpi(x, neighbourhood_size = scale)$tpi_standardised

  # Calculate the slope from the input DTM, to be used for either the six or ten class slope position
  slp <- terra::terrain(x, v="slope", unit="degrees", neighbors=8)

  if (class.type == "slope.position") {

    # Define the six classes on the basis of thresholds of tp and slope
    valley <- (tp <= -1)*1
    lower.slp <- (tp > -1 & tp <= -0.5)*2
    flat.slp <- ((tp > -0.5 & tp < 0.5) & (slp <= 5))*3
    middle.slp <- ((tp > -0.5 & tp < 0.5) & (slp > 5))*4
    upper.slp <- (tp > 0.5 & tp <= 1)*5
    ridge <- (tp > 1)*6

    result <- valley + lower.slp + flat.slp + middle.slp + upper.slp + ridge

    # Plot the result
    # Define the number of colors matching the number of classes
    colors <- rainbow(6)

    # Plot the landform raster
    terra::plot(result,
                col = colors,
                legend=FALSE)

    labels <- c("1: Valley", "2: Lower Slope", "3: Flat Slope", "4: Middle Slope", "5: Upper Slope", "6: Ridge")

    legend(leg.pos,
           legend = labels,
           fill = colors,
           title = "Slope Position (Weiss)",
           cex = leg.cex)

    if(descriptive) {
      # Calculate the count, area and percentage of each class
      class_count <- table(factor(terra::values(result), levels = 1:6)) # Ensure all classes are included
      class_area <- round(((class_count * terra::res(x)[1] * terra::res(x)[2]) / 1000000),3)
      class_percentage <- round((class_count / sum(class_count) * 100),3)

      # Create a data frame
      df <- data.frame(
        Class = names(class_count),
        Count = as.numeric(class_count),
        Area = as.numeric(class_area),
        Percentage = as.numeric(class_percentage)
      )

      # Define class names
      class_names <- c("Valley",
                       "Lower Slope",
                       "Flat Slope",
                       "Middle Slope",
                       "Upper Slope",
                       "Ridge")

      # Map class numbers to names
      df$Class <- class_names[as.numeric(df$Class)]

      # Generate a bar plot
      barplot(df$Area,
              names.arg = 1:nrow(df),
              main = "Area (Km2)",
              xlab = "Class",
              col = rainbow(nrow(df)),
              cex.names=0.7,
              cex.main=0.85,
              cex.axis = 0.85)

      # Add a legend
      legend("topright",
             legend = df$Class,
             fill = rainbow(nrow(df)),
             title = "Classes",
             cex = 0.7)

      return(list(
        "all"=result,
        "Valley" = valley,
        "Lower Slope" = lower.slp,
        "Flat Slope" = flat.slp,
        "Middle Slope" = middle.slp,
        "Upper Slope" = upper.slp,
        "Ridge" = ridge,
        "Descriptive Statistics" = df))
    }

    return(list(
      "all"=result,
      "Valley" = valley,
      "Lower Slope" = lower.slp,
      "Flat Slope" = flat.slp,
      "Middle Slope" = middle.slp,
      "Upper Slope" = upper.slp,
      "Ridge" = ridge))
  }

  if (class.type == "landform.classification") {

    # Calculate two standardized tpi, one with small neighbour, one with large neighbour
    sn <- calculate_standardised_tpi(x, neighbourhood_size = sn)$tpi_standardised
    ln <- calculate_standardised_tpi(x, neighbourhood_size = ln)$tpi_standardised

    # Define the ten classes on the basis of thresholds of sn, sl, and slope
    canyons <- ((sn <= -1) & (ln <= -1))*1
    midslope.dr <- ((sn <= -1) & (ln > -1 & ln < 1))*2
    upland.dr <-  ((sn <= -1) & (ln >= 1))*3
    us.valley <-  ((sn > -1 & sn < 1) & (ln <=-1))*4
    plains <- ((sn > -1 & sn < 1) & (ln > -1 & ln < 1) & (slp <= 5))*5
    open.slp <-  ((sn > -1 & sn < 1) & (ln > -1 & ln < 1) & (slp > 5))*6
    upper.slp <- ((sn > -1 & sn < 1) & (ln >= 1))*7
    local.rdg <- ((sn >= 1) & (ln <= -1))*8
    midslp.rdg <- ((sn >= 1) & (ln > -1 & ln < 1))*9
    mount.top <- ((sn >= 1) & (ln >=1))*10

    result <- canyons + midslope.dr + upland.dr + us.valley + plains + open.slp + upper.slp + local.rdg + midslp.rdg + mount.top

    # Plot the result
    colors <- rainbow(10)

    labels <- c("1: Canyons-Deeply Incised Streams",
                "2: Midslope Drainage-Shallow Valleys",
                "3: Upland Drainages-Headwaters",
                "4: U-shaped Valleys",
                "5: Plains",
                "6: Open Slopes",
                "7: Upper Slopes-Mesas",
                "8: Local Ridges-Hills in Valley",
                "9: Midslopes Ridges-Small Hills in Plains",
                "10: Mountain Tops-High Ridges")

    terra::plot(result,
                col = colors,
                legend=FALSE)

    legend(leg.pos,
           legend = labels,
           fill = colors,
           title = "Landform Classification (Jenness)",
           cex = leg.cex)

    if(descriptive) {
      # Calculate the count, area and percentage of each class
      class_count <- table(factor(terra::values(result), levels = 1:10)) # Ensure all classes are included
      class_area <- round(((class_count * terra::res(x)[1] * terra::res(x)[2]) / 1000000),3)
      class_percentage <- round((class_count / sum(class_count) * 100),3)

      # Create a data frame
      df <- data.frame(
        Class = names(class_count),
        Count = as.numeric(class_count),
        Area = as.numeric(class_area),
        Percentage = as.numeric(class_percentage)
      )

      class_names <- c("Canyons-Deeply Incised Streams",
                       "Midslope Drainage-Shallow Valleys",
                       "Upland Drainages-Headwaters",
                       "U-shaped Valleys",
                       "Plains",
                       "Open Slopes",
                       "Upper Slopes-Mesas",
                       "Local Ridges-Hills in Valley",
                       "Midslopes Ridges-Small Hills in Plains",
                       "Mountain Tops-High Ridges")

      # Map class numbers to names
      df$Class <- class_names[as.numeric(df$Class)]

      # Generate a bar plot
      barplot(df$Area,
              names.arg = 1:nrow(df),
              main = "Area (Km2)",
              xlab = "Class",
              col = rainbow(nrow(df)),
              cex.names=0.7,
              cex.main=0.85,
              cex.axis = 0.85)

      # Add a legend
      legend("topright",
             legend = df$Class,
             fill = rainbow(nrow(df)),
             title = "Classes",
             cex = 0.7)


      return(list(
        "all"=result,
        "Canyons-Deeply Incised Streams" = canyons,
        "Midslope Drainage-Shallow Valleys" = midslope.dr,
        "Upland Drainages-Headwaters" = upland.dr,
        "U-shaped Valleys" = us.valley,
        "Plains" = plains,
        "Open Slopes" = open.slp,
        "Upper Slopes-Mesas" = upper.slp,
        "Local Ridges-Hills in Valley" = local.rdg,
        "Midslopes Ridges-Small Hills in Plains" = midslp.rdg,
        "Mountain Tops-High Ridges" = mount.top,
        "Descriptive Statistics" = df))
    }

    return(list(
      "all"=result,
      "Canyons-Deeply Incised Streams" = canyons,
      "Midslope Drainage-Shallow Valleys" = midslope.dr,
      "Upland Drainages-Headwaters" = upland.dr,
      "U-shaped Valleys" = us.valley,
      "Plains" = plains,
      "Open Slopes" = open.slp,
      "Upper Slopes-Mesas" = upper.slp,
      "Local Ridges-Hills in Valley" = local.rdg,
      "Midslopes Ridges-Small Hills in Plains" = midslp.rdg,
      "Mountain Tops-High Ridges" = mount.top))
  }
}

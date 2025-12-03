#' Custom Geom for Encircling Points in ggplot2
#' @importFrom ggplot2 ggproto Geom aes layer
#' @importFrom grid unit convertUnit get.gpar gpar xsplineGrob grobTree rectGrob
#' @importFrom grDevices chull
#' @importFrom scales alpha
#' @format NULL
#' @usage NULL
#' @author Jared Andrews, heavily based on ggalt code from Ben Bolker
#'   (\url{https://github.com/hrbrmstr/ggalt/blob/master/R/geom_encircle.r})
#' @export
GeomEncircle <- ggproto("GeomEncircle", Geom,
  required_aes = c("x", "y"),
  default_aes = aes(
    colour = "black",
    fill = NA,
    alpha = 1,
    linetype = 1,
    size = 1,
    s_shape = 0.5,
    s_open = FALSE,
    expand = 0.05,
    spread = 0.1
  ),
  
  draw_group = function(data, panel_scales, coord) {
    # Apply coordinate transformation
    transformed_coords <- coord$transform(data, panel_scales)
    reference_row <- transformed_coords[1, , drop = FALSE]
    rownames(reference_row) <- NULL
    
    # Find center point of all coordinates
    center_point <- lapply(transformed_coords[, c("x", "y")], mean, na.rm = TRUE)
    
    # Get convex hull vertex indices
    hull_vertices <- grDevices::chull(transformed_coords[c("x", "y")])
    
    # Factory function for coordinate data frames
    build_coord_frame <- function(x_coords, y_coords) {
      non_xy_cols <- reference_row[!names(reference_row) %in% c("x", "y")]
      data.frame(x = x_coords, y = y_coords, non_xy_cols)
    }
    
    transformed_coords <- transformed_coords[hull_vertices, ]
    
    # Coordinate conversion utilities
    native_to_millimeters <- function(value, dimension = "x") {
      grid::convertUnit(
        grid::unit(value, "native"), "mm",
        typeFrom = "dimension", axisFrom = dimension, valueOnly = TRUE
      )
    }
    
    location_to_native <- function(value, dimension = "x") {
      grid::convertUnit(
        value, "native",
        typeFrom = "location", axisFrom = dimension, valueOnly = TRUE
      )
    }
    
    merge_native_snpc <- function(native_component, snpc_component, dimension = "x") {
      location_to_native(
        unit(native_component, "native") + unit(snpc_component, "snpc"),
        dir = dimension
      )
    }
    
    # Calculate unit vector from point2 to point1
    unit_vector <- function(point1, point2) {
      x_diff <- native_to_millimeters(point1$x - point2$x)
      y_diff <- native_to_millimeters(point1$y - point2$y)
      vector_length <- sqrt(x_diff^2 + y_diff^2)
      list(x = x_diff / vector_length, y = y_diff / vector_length)
    }
    
    # Handle edge cases for small point sets
    if (nrow(transformed_coords) == 1) {
      # Single point: expand into diamond
      transformed_coords <- with(transformed_coords, build_coord_frame(
        c(x, x + spread, x, x - spread),
        c(y + spread, y, y - spread, y)
      ))
    } else if (nrow(transformed_coords) == 2) {
      # Two points: create perpendicular diamond
      perpendicular_rotation <- matrix(c(0, 1, -1, 0), nrow = 2)
      direction <- unit_vector(transformed_coords[1, ], transformed_coords[2, ])
      offset <- c(perpendicular_rotation %*% unlist(direction)) * 
                transformed_coords$spread
      
      transformed_coords <- with(transformed_coords, {
        new_x <- c(x[1], center_point$x + offset[1], x[2], center_point$x - offset[1])
        new_y <- c(y[1], center_point$y + offset[2], y[2], center_point$y - offset[2])
        build_coord_frame(new_x, new_y)
      })
    }
    
    # Calculate outward directions from center
    outward_vectors <- unit_vector(transformed_coords, center_point)
    
    # Configure graphics parameters
    gpar_config <- grid::get.gpar()
    aesthetic_inputs <- c("colour", "linetype", "alpha", "fill", "size")
    gpar_outputs <- c("col", "lty", "alpha", "fill", "lwd")
    gpar_config[gpar_outputs] <- reference_row[aesthetic_inputs]
    
    # Generate the encircling spline
    grid::xsplineGrob(
      with(transformed_coords, 
           unit(x, "npc") + outward_vectors$x * unit(expand, "snpc")),
      with(transformed_coords, 
           unit(y, "npc") + outward_vectors$y * unit(expand, "snpc")),
      shape = transformed_coords$s_shape - 1,
      open = reference_row$s_open,
      gp = gpar_config
    )
  }
)

#' Automatically enclose points in a polygon
#'
#' Creates a smooth encircling polygon around a set of points using convex hull
#' calculation and xspline smoothing. Useful for highlighting groups of points
#' in scatter plots.
#'
#' @param mapping Set of aesthetic mappings created by \code{\link[ggplot2]{aes}}. 
#'   If specified and \code{inherit.aes = TRUE} (the default), it is combined with 
#'   the default mapping at the top level of the plot.
#' @param data The data to be displayed in this layer. If \code{NULL}, the default,
#'   the data is inherited from the plot data as specified in the call to 
#'   \code{\link[ggplot2]{ggplot}}.
#' @param stat The statistical transformation to use on the data for this layer,
#'   as a string.
#' @param position Position adjustment, either as a string, or the result of a call
#'   to a position adjustment function.
#' @param na.rm If \code{FALSE}, the default, missing values are removed with a warning.
#'   If \code{TRUE}, missing values are silently removed.
#' @param show.legend Logical. Should this layer be included in the legends?
#'   \code{NA}, the default, includes if any aesthetics are mapped.
#' @param inherit.aes If \code{FALSE}, overrides the default aesthetics, rather
#'   than combining with them.
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}. These are
#'   often aesthetics, used to set an aesthetic to a fixed value, like 
#'   \code{colour = "red"} or \code{size = 3}. They may also be parameters to the
#'   paired geom/stat. Additional parameters include:
#'   \describe{
#'     \item{s_shape}{Controls the shape of the spline (default = 0.5).}
#'     \item{s_open}{Logical indicating whether the spline should be open (default = FALSE).}
#'     \item{expand}{Amount to expand the encircling polygon outward (default = 0.05).}
#'     \item{spread}{Spread factor for single or double point sets (default = 0.1).}
#'   }
#'
#' @return A ggplot2 layer that can be added to a plot.
#'
#' @author Jared Andrews, heavily based on ggalt code from Ben Bolker
#'   (\url{https://github.com/hrbrmstr/ggalt/blob/master/R/geom_encircle.r})
#'
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' d <- data.frame(x=c(1,1,2),y=c(1,2,2)*100)
#'
#' gg <- ggplot(d,aes(x,y))
#' gg <- gg + scale_x_continuous(expand=c(0.5,1))
#' gg <- gg + scale_y_continuous(expand=c(0.5,1))
#'
#' gg + geom_encircle(s_shape=1, expand=0) + geom_point()
#'
#' gg + geom_encircle(s_shape=1, expand=0.1, colour="red") + geom_point()
#'
#' gg + geom_encircle(s_shape=0.5, expand=0.1, colour="purple") + geom_point()
#'
#' gg + geom_encircle(data=subset(d, x==1), colour="blue", spread=0.02) +
#'   geom_point()
#'
#' gg + geom_encircle(data=subset(d, x==2), colour="cyan", spread=0.04) +
#'   geom_point()
#'
#' gg <- ggplot(mpg, aes(displ, hwy))
#' gg + geom_encircle(data=subset(mpg, hwy>40)) + geom_point()
#' gg + geom_encircle(aes(group=manufacturer)) + geom_point()
#' gg + geom_encircle(aes(group=manufacturer,fill=manufacturer),alpha=0.4)+
#'        geom_point()
#' gg + geom_encircle(aes(group=manufacturer,colour=manufacturer))+
#'        geom_point()
#'
#' ss <- subset(mpg,hwy>31 & displ<2)
#'
#' gg + geom_encircle(data=ss, colour="blue", s_shape=0.9, expand=0.07) +
#'   geom_point() + geom_point(data=ss, colour="blue")
#' }
geom_encircle <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity", na.rm = FALSE, show.legend = NA,
                          inherit.aes = TRUE, ...) {
  layer(
    geom = GeomEncircle, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

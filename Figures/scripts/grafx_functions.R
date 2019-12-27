library(ggplot2)
library(igraph)

# graphix helper functions (aging paper)

### ggplot functions
# convex hull (source: http://docs.ggplot2.org/dev/vignettes/extending-ggplot2.html )
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     required_aes = c("x", "y")
)
stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 0) xminv else xmaxv), if (grp %% 2 == 0) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               aesthetics$colour <- rep("black", nrow(quantiles))
                               if (grp %% 2 == 0) {
                                 quantiles$x[2] <- data$xminv[1]
                               } else {
                                 quantiles$x[1] <- data$xminv[1]
                               }
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

### igraph functions
# triangle vertex shape (igraph)
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/150 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,fg = vertex.frame.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

## Function for plotting an elliptical node
myellipse <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/100 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  draw.ellipse(x=coords[,1], y=coords[,2],
               a = vertex.size, b=vertex.size/2, col=vertex.color)
}
add_shape("ellipse", clip=shapes("circle")$clip,
          plot=myellipse)

## Function to draw diamond shape
myrhombus <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color,
          stars = cbind(1.2*vertex.size, vertex.size, 1.2*vertex.size, vertex.size),
          add = TRUE, inches = FALSE)
}
# clips as a circle
add_shape("rhombus", clip = shapes("circle")$clip,
          plot = myrhombus)


### Data conversion functions
get_density <- function(x, y, n = 250) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
vec2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x,na.rm = T)
  pal[findInterval(x,seq(-max(abs(limits)),max(abs(limits)),length.out=length(pal)+1), all.inside=TRUE)]
}

# Z-score normalization
znorm <- function(x) (x-mean(x))/sd(x)

detach("package:igraph", unload=TRUE)
#' Plot Functional Time Series (funts) with Plotly
#'
#' Visualize univariate or multivariate Functional Time Series (funts) using Plotly-based plots.
#'
#' @param x An object of class \code{\link{funts}}.
#' @param vars Numeric specifying which variables in the FTS to plot (default: all).
#' @param types Tuple of strings specifying plot types for each variable.
#' @param subplot Logical for subplotting line plots.
#' @param main Titles for each plot.
#' @param ylab Y-axis titles.
#' @param xlab X-axis titles.
#' @param tlab Time-axis titles.
#' @param zlab Z-axis titles.
#' @param xticklabels Tick labels for the domain of the functions.
#' @param xticklocs Positions of tick labels for the domain of the functions.
#' @param yticklabels Tick labels for the domain of the functions.
#' @param yticklocs Positions of tick labels for the domain of the functions.
#' @param color_palette Color palette for two-dimensional FTS plots.
#' @param reverse_color_palette Reverse the color palette scale.
#' @param ... Additional arguments to pass to Plotly methods.
#'
#' @details
#' Supported plot types for one-dimensional domain variables:
#'   - "line": Line plots (default).
#'   - "heatmap": Heatmaps.
#'   - "3Dsurface": 3D surface plots.
#'   - "3Dline": 3D line plots.
#'
#' Supported plot type for two-dimensional domain variables:
#'   - "heatmap"
#'
#' Each variable can be plotted multiple times with different types.
#'
#' @importFrom plotly plot_ly add_lines layout subplot add_surface hide_colorbar ggplotly
#' @importFrom ggplot2 ggplot aes_string unit geom_tile scale_fill_distiller xlab ylab labs ggtitle scale_y_continuous scale_x_continuous waiver theme element_line element_text element_blank
#' @importFrom tibble as_tibble
#'
#' @seealso \code{\link{funts}}, \code{\link{Callcenter}}, \code{\link{Montana}}
#'
#' @examples
#' data("Callcenter") # Univariate FTS example
#' plotly_funts(Callcenter, main = "Call Center Data Line Plot",
#'              xlab = "Time (6 minutes aggregated)",
#'              ylab = "Sqrt of Call Numbers", type = "line",
#'              xticklabels = list(c("00:00","06:00","12:00","18:00","24:00")),
#'              xticklocs = list(c(1,60,120,180,240)))
#'
#' data("Montana") # Multivariate FTS example
#' plotly_funts(Montana[1:100],
#'              xlab = c("Time", "Longitude"),
#'              ylab = c("Normalized Temperature (\u00B0C)", "Latitude"),
#'              zlab = c("", "NDVI"),
#'              main = c("Temperature Curves", "NDVI Images"),
#'              color_palette = "RdYlGn",
#'              xticklabels = list(c("00:00","06:00","12:00","18:00","24:00"),
#'                                 c("113.40\u00B0 W", "113.30\u00B0 W")),
#'              xticklocs = list(c(1,6,12,18,24),c(1,33)),
#'              yticklabels = list(NA,c("48.70\u00B0 N", "48.77\u00B0 N")),
#'              yticklocs = list(NA,c(1,33))
#' )
#'
#' @export
plotly_funts <- function(x, vars = NULL, types = NULL, subplot = TRUE, main = NULL, ylab = NULL, xlab = NULL, tlab = NULL,
                         zlab = NULL, xticklabels = NULL, xticklocs = NULL, yticklabels = NULL, yticklocs = NULL,
                         color_palette = "RdYlBu", reverse_color_palette = FALSE, ...) {
  p <- length(x$dimSupp)
  N <- x$N
  time <- x$time
  count_twod <- 0
  Pl <- list()
  if (is.null(types)) types <- rep(NA, p)
  if (is.null(ylab)) ylab <- rep(NA, p)
  if (is.null(xlab)) xlab <- rep(NA, p)
  if (is.null(tlab)) tlab <- rep(NA, p)
  if (is.null(zlab)) zlab <- rep(NA, p)
  if (is.null(xticklabels)) xticklabels <- as.list(rep(NA, p))
  if (is.null(xticklocs)) xticklocs <- as.list(rep(NA, p))
  if (is.null(yticklabels)) yticklabels <- as.list(rep(NA, p))
  if (is.null(yticklocs)) yticklocs <- as.list(rep(NA, p))
  if (is.null(zlab)) zlab <- rep(NA, p)
  if (is.null(zlab)) zlab <- rep(NA, p)
  if (is.null(main)) main <- rep(NA, p)
  if (!is.null(vars) && length(types) != length(vars)) warning("\"vars\" and \"types\" are not the same length. Some plots might not appear as expected.")
  if (is.null(vars)) vars <- 1:p
  cat("Plotting, please wait...\n")
  for (j in 1:length(vars)) {
    if (j > length(xticklocs) || j > length(xticklabels)) xticklocs[[j]] <- xticklabels[[j]] <- NA
    if (j > length(yticklocs) || j > length(yticklabels)) yticklocs[[j]] <- yticklabels[[j]] <- NA
    if ((is.na(xticklocs[[j]][1]) && !is.na(xticklabels[[j]][1])) || (!is.na(xticklocs[[j]][1]) && is.na(xticklabels[[j]][1]))) warning(paste0("Please provide horizontal axis labels and label locations for variable ", as.character(j), "."))
    if ((is.na(yticklocs[[j]][1]) && !is.na(yticklabels[[j]][1])) || (!is.na(yticklocs[[j]][1]) && is.na(yticklabels[[j]][1]))) warning(paste0("Please provide vertical axis labels and label locations for variable ", as.character(j), "."))
    if (is.numeric(x$argval[[vars[j]]]) && !is.matrix(x$argval[[vars[j]]])) {
      if (is.na(types[j]) || types[j] == "line") {
        if (is.na(ylab[j])) ylab[j] <- "y"
        if (is.na(xlab[j])) xlab[j] <- "x"
        if (is.na(main[j])) main[j] <- paste("Variable", vars[j])
        y <- tibble::as_tibble(data.frame(y = c(x$B_mat[[vars[j]]] %*% x$coefs[[vars[j]]])))
        y$time <- as.factor(rep(time, each = length(x$argval[[vars[j]]])))
        y$x <- rep(1:length(x$argval[[vars[j]]]), ncol(x$coefs[[vars[j]]]))
        if ((!is.na(xticklabels[[j]][1]) && !is.na(xticklocs[[j]][1])) && (!is.na(yticklabels[[j]][1]) && !is.na(yticklocs[[j]][1]))) {
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x = ~x, y = ~y) %>%
            add_lines(color = ~time, colors = c("lightsteelblue", "royalblue4"), showlegend = FALSE) %>%
            layout(
              title = paste("<b>", main[j], "</b>"), font = list(size = 12), yaxis = list(
                title = ylab[j], ticktext = yticklabels[[j]], tickvals = yticklocs[[j]],
                range = list(yticklocs[[j]][1], yticklocs[[j]][length(yticklocs[[j]])])
              ),
              xaxis = list(title = xlab[j], ticktext = xticklabels[[j]], tickvals = xticklocs[[j]], range = list(xticklocs[[j]][1], xticklocs[[j]][length(xticklocs[[j]])]))
            )
        } else if ((!is.na(xticklabels[[j]][1]) && !is.na(xticklocs[[j]][1])) && (is.na(yticklabels[[j]][1]) || is.na(yticklocs[[j]][1]))) {
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x = ~x, y = ~y) %>%
            add_lines(color = ~time, colors = c("lightsteelblue", "royalblue4"), showlegend = FALSE) %>%
            layout(
              title = paste("<b>", main[j], "</b>"), font = list(size = 12), yaxis = list(title = ylab[j]),
              xaxis = list(title = xlab[j], ticktext = xticklabels[[j]], tickvals = xticklocs[[j]], range = list(xticklocs[[j]][1], xticklocs[[j]][length(xticklocs[[j]])]))
            )
        } else if ((is.na(xticklabels[[j]][1]) || is.na(xticklocs[[j]][1])) && (is.na(yticklabels[[j]][1]) == FALSE && is.na(yticklocs[[j]][1]) == FALSE)) {
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x = ~x, y = ~y) %>%
            add_lines(color = ~time, colors = c("lightsteelblue", "royalblue4"), showlegend = FALSE) %>%
            layout(
              title = paste("<b>", main[j], "</b>"), font = list(size = 12), yaxis = list(
                title = ylab[j], ticktext = yticklabels[[j]], tickvals = yticklocs[[j]],
                range = list(yticklocs[[j]][1], yticklocs[[j]][length(yticklocs[[j]])])
              ),
              xaxis = list(title = xlab[j])
            )
        } else {
          Pl[[j]] <- y %>%
            group_by(time) %>%
            plot_ly(x = ~x, y = ~y) %>%
            add_lines(color = ~time, colors = c("lightsteelblue", "royalblue4"), showlegend = FALSE) %>%
            layout(
              title = paste("<b>", main[j], "</b>"), font = list(size = 12), yaxis = list(title = ylab[j]),
              xaxis = list(title = xlab[j])
            )
        }
      } else if (types[j] == "heatmap") {
        u <- seq(x$argval[[vars[j]]])
        if (is.na(xlab[j])) xlab[j] <- "x"
        if (is.na(tlab[j])) tlab[j] <- "t"
        if (is.na(zlab[j])) zlab[j] <- "z"
        if (is.na(main[j])) main[j] <- paste("Variable", vars[j])
        z0 <- x$B_mat[[vars[j]]] %*% x$coefs[[vars[j]]]
        if (is.na(xticklabels[[j]][1]) || is.na(xticklocs[[j]][1])) {
          Pl[[j]] <- plot_ly(
            z = z0, x = time, y = u, type = "heatmap", colorscale = list(c(0, "#FFFFFAFF"), c(1, "#FF0000FF")),
            showscale = FALSE, hovertemplate = paste0(zlab[j], ":", " %{z}", "\n", tlab[j], ":", " %{x}", "\n", xlab[j], ":", " %{y}")
          ) %>%
            layout(
              title = paste("<b>", main[j], "</b>"), font = list(size = 12), yaxis = list(title = xlab[j]),
              xaxis = list(title = tlab[j])
            )
        } else {
          Pl[[j]] <- plot_ly(
            z = z0, x = time, y = u, type = "heatmap", colorscale = list(c(0, "#FFFFFAFF"), c(1, "#FF0000FF")),
            showscale = FALSE, hovertemplate = paste0(zlab[j], ":", " %{z}", "\n", tlab[j], ":", " %{x}", "\n", xlab[j], ":", " %{y}")
          ) %>%
            layout(title = paste("<b>", main[j], "</b>"), font = list(size = 12), yaxis = list(title = xlab[j], ticktext = xticklabels[[j]], tickvals = xticklocs[[j]]), xaxis = list(title = tlab[j]))
        }
      } else if (types[j] == "3Dsurface") {
        u <- x$argval[[vars[j]]]
        z0 <- x$B_mat[[vars[j]]] %*% x$coefs[[vars[j]]]
        axx <- axy <- axz <- list(
          gridcolor = "rgb(180, 180, 180)",
          zerolinecolor = "rgb(255,255,255)"
        )
        axx$title <- ifelse(is.na(tlab[j]), "t", tlab[j])
        axy$title <- ifelse(is.na(xlab[j]), "x", xlab[j])
        axz$title <- ifelse(is.na(zlab[j]), paste("Variable", vars[j]), zlab[j])
        if (is.na(xticklabels[[j]][1]) == FALSE || is.na(xticklocs[[j]][1]) == FALSE) {
          axy$ticktext <- xticklabels[[j]]
          axy$tickvals <- xticklocs[[j]]
        }

        Pl[[j]] <- plot_ly(
          z = z0, x = time, y = u, colorscale = list(c(0, "#FFFFFAFF"), c(1, "#FF0000FF")),
          hovertemplate = paste0(axz$title[j], ":", " %{z}", "\n", axx$title[j], ":", " %{x}", "\n", axy$title[j], ":", " %{y}")
        ) %>%
          layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>%
          add_surface(showscale = FALSE)
      } else if (types[j] == "3Dline") {
        D0 <- tibble::as_tibble(data.frame(z = c(x$B_mat[[vars[j]]] %*% x$coefs[[vars[j]]])))
        D0$time <- rep(time, each = length(x$argval[[vars[j]]]))
        D0$x <- rep(1:length(x$argval[[vars[j]]]), ncol(x$coefs[[vars[j]]]))
        axx <- axy <- axz <- list(
          gridcolor = "rgb(180, 180, 180)",
          zerolinecolor = "rgb(255,255,255)"
        )
        axx$title <- ifelse(is.na(tlab[j]), "t", tlab[j])
        axy$title <- ifelse(is.na(xlab[j]), "x", xlab[j])
        axz$title <- ifelse(is.na(zlab[j]), paste("Variable", vars[j]), zlab[j])
        if (is.na(xticklabels[[j]][1]) == FALSE || is.na(xticklocs[[j]][1]) == FALSE) {
          axy$ticktext <- xticklabels[[j]]
          axy$tickvals <- xticklocs[[j]]
        }

        Pl[[j]] <- D0 %>%
          group_by(time) %>%
          plot_ly(
            x = ~time, z = ~z, y = ~x, type = "scatter3d", mode = "lines", color = ~z,
            line = list(width = 4), colors = c("#FFFFFAFF", "#FF0000FF"), hovertemplate = paste0(axz$title[j], ":", " %{z}", "\n", axx$title[j], ":", " %{x}", "\n", axy$title[j], ":", " %{y}")
          ) %>%
          layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz)) %>%
          hide_colorbar()
      }
    } else {
      if (is.na(types[j]) == FALSE && types[j] != "heatmap") warning("The only plotting option available for variables observed over two-dimensional domain is \"heatmap\". Other plotting options will be added for these types of variables in the future.")
      count_twod <- count_twod + 1
      if (is.na(ylab[j])) ylab[j] <- "y"
      if (is.na(xlab[j])) xlab[j] <- "x"
      if (is.na(zlab[j])) zlab[j] <- "z"
      if (is.na(main[j])) main[j] <- paste("Variable", vars[j])
      if (is.na(xticklabels[[j]][[1]])) xticklabels[[j]] <- waiver()
      if (is.na(xticklocs[[j]][[1]])) xticklocs[[j]] <- waiver()
      if (is.na(yticklabels[[j]][[1]])) yticklabels[[j]] <- waiver()
      if (is.na(yticklocs[[j]][[1]])) yticklocs[[j]] <- waiver()
      temp_df <- data.frame(matrix(ncol = 3, nrow = ncol(x$coefs[[vars[j]]]) * nrow(x$argval[[vars[j]]])))
      colnames(temp_df) <- c(zlab[j], xlab[j], ylab[j])
      temp_df[, 1] <- c(x$B_mat[[vars[j]]] %*% x$coefs[[vars[j]]])
      temp_df[, 2] <- rep((x$argval[[vars[j]]][, 2]), ncol(x$coefs[[vars[j]]]))
      temp_df[, 3] <- rev(rep((x$argval[[vars[j]]][, 1]), ncol(x$coefs[[vars[j]]])))
      y <- tibble::as_tibble(temp_df)
      # The next line is to help alleviate a bug in ggplot2
      if (time[1] %in% as.character(1:ncol(x$coefs[[vars[j]]]))) time <- as.numeric(time)
      y$Time <- rep(time, each = nrow(x$argval[[vars[j]]]))
      direction <- 1
      if (isTRUE(reverse_color_palette)) direction <- -1
      Pl[[j]] <- ggplotly(ggplot(y, aes_string(xlab[j], ylab[j], fill = zlab[j], frame = "Time")) +
        geom_tile() +
        theme(
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")
        ) +
        scale_fill_distiller(palette = color_palette, direction = direction) +
        labs(fill = zlab[j]) +
        ggtitle(main[j]) +
        scale_y_continuous(
          name = ylab[j], breaks = yticklocs[[j]],
          labels = yticklabels[[j]]
        ) +
        scale_x_continuous(
          name = xlab[j], breaks = xticklocs[[j]],
          labels = xticklabels[[j]]
        ) +
        theme(text = element_text(size = 10))
        +
        theme(
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 15)
        )
        +
        theme(
          axis.text.y = element_text(angle = 90, size = 10),
          axis.title.y = element_text(angle = 90, size = 15)
        )
        +
        theme(plot.title = element_text(size = 30, face = "bold"))
        +
        theme(
          legend.position = "right", legend.key.height = unit(0.67, "cm"), legend.title = element_text(size = 15),
          legend.text = element_text(size = 10)
        ))
    }
  }
  if (p == 1 || "heatmap" %in% types || "3Dsurface" %in% types || "3Dline" %in% types || count_twod >= 1 || subplot == FALSE || length(vars) == 1) {
    for (i in 1:p) print(Pl[[i]])
  } else {
    print(subplot(Pl, titleX = TRUE, titleY = TRUE) %>% layout(title = ""))
  }
  cat("Done.\n")
}

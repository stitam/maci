#' Plot frequencies of character states over time
#' 
#' This is a convenience function for plotting frequencies (counts or ratios)
#' of character states over time. 
#' @param df data.frame; a data frame
#' @param var character; the name of the character variable to plot
#' @param by character; the name of the character variable to group by
#' @param timevar character; the name of the variable which stores the year of
#' isolation. This variable must be numeric.
#' @param minyear numeric; samples isolated before this year will be removed
#' from the data set.
#' @param window numeric; number of preceding years to aggregate.
#' @param type character; defines transformation for the y axis. Can be either
#' \code{"count"} or \code{"ratio"}.
#' @param geom character; defines the plot type. Can be \code{"geom_line"},
#' \code{"geom_col"} or \code{"geom_area"}.
#' @param show_count logical; show sample counts for each year.
#' @param count_size integer; text size when showing counts for each year.
#' @param colors data.frame; a data frame for colors to use for color fill.
#' If \code{NULL}, the colors will be chosen automatically. If specified, the
#' data frame must have two columns, the first variable must have the same name
#' as the \code{var} parameter, the other must be \code{"color"}.  
#' @examples
#' \dontrun{
#' plot_ts(df = aci, var = "continent")
#' }
#' @import ggplot2
#' @importFrom qualpalr qualpal
#' @export
plot_ts <- function(df,
                    var,
                    by = NULL,
                    timevar = "collection_year",
                    minyear = 2000,
                    window = 2,
                    type = "count",
                    geom = "geom_line",
                    show_count = TRUE,
                    count_size = 8,
                    colors = NULL) {
  if ("data.frame" %in% class(df) == FALSE) {
    stop ("Argument 'df' must be a data.frame.")
  }
  if (var %in% names(df) == FALSE) {
    stop(paste0("Variable '", var, "' not found."))
  }
  if (!is.null(by)) {
    if (by %in% names(df) == FALSE) {
      stop(paste0("Variable '", by, "' not found."))
    }
    if (length(unique(df[[var]])) > 2) {
      stop(paste0("If 'by' is not NULL, ",
                  "then 'var' must contain a vector of class 'logical'."))
    }
  }
  if (timevar %in% names(df) == FALSE) {
    stop(paste0("Variable '", timevar, "' not found."))
  }
  if (!is.numeric(df[[timevar]])) {
    stop(paste0("Variable '", timevar, "' must be numeric."))
  }
  if (!is.numeric(minyear)) {
    stop(paste("Variable '", minyear, "' must be numeric."))
  }
  if (!is.numeric(window)) {
    stop(paste("Variable '", window, "' must be numeric."))
  }
  ### stop conditions because not all parameter combinations are compatioble
  if (type == "count" & geom != "geom_col") {
    stop("If 'type' is 'count' then 'geom' must be 'geom_col'.")
  }
  
  ###
  type <- match.arg(type, choices = c("count", "ratio"))
  geom <- match.arg(geom, choices = c("geom_line", "geom_col", "geom_area"))
  if (!is.null(colors) && !is.data.frame(colors)) {
    stop(paste("Variable '", colors, "' must be a data frame."))
  }
  df[[timevar]] <- as.numeric(as.character(df[[timevar]]))
  df <- df[which(df[[timevar]] >= minyear), ]
  if (is.null(by)) {
    df <- as.data.frame(table(df[[timevar]], df[[var]]))
    names(df) <- c("year", var, "count")
    df$year <- as.numeric(as.character(df$year))
    df$yearcount <- NA
    for (i in 1:nrow(df)) {
      df$yearcount[i] <- sum(df$count[which(df$year == df$year[i])])
    }
    df$mean_count <- NA
    df$mean_ratio <- NA
    for (i in 1:nrow(df)) {
      index <- which(
        df[[var]] == df[[var]][i] & 
          df$year <= df$year[i] &
          df$year > df$year[i]-window)
      df$mean_count[i] <- mean(df$count[index])
      index_year <- which(df$year <= df$year[i] & df$year > df$year[i]-window)
      df$mean_ratio[i] <- round(
        sum(df$count[index])/sum(df$count[index_year]), 3)
    }
    if (is.null(colors)) {
      cls <- qualpalr::qualpal(
        length(unique(df[[var]])), colorspace = "pretty")$hex
    } else {
      cls <- colors$color
      names(cls) <- colors[[var]]
      if (length(cls) != length(unique(df[[var]]))) {
        stop(paste0("Variable '", var, "' requires ",
                    length(unique(df[[var]])), " separate colors."))
      }
    }
    df$year <- switch(
      geom,
      geom_line = df$year,
      geom_col = as.character(df$year),
      geom_area = df$year
    )
    ylab <- paste0(
      switch(type, count = "Number", ratio = "Ratio"),
      " of isolates"
    )
    if (window > 1) {
      ylab <- paste0(ylab, " (moving average, past ", window, " years)")
    }
    g <- ggplot(
      df,
      aes(year, switch(type, count = mean_count, ratio = mean_ratio)))+
      switch(
        geom,
        geom_line = geom_line(aes(col = get(var)), alpha = 0.9),
        geom_col = geom_col(aes(fill = get(var)), alpha = 0.9),
        geom_area = geom_area(aes(fill = get(var)), alpha = 0.9)
      )+
      theme(panel.background = element_rect(fill = "white"))+
      theme(axis.line = element_line(colour = "grey"))+
      switch(
        geom,
        geom_line = theme(axis.text.x = element_text(angle = 0)),
        geom_col = theme(axis.text.x = element_text(angle = 90)),
        geom_area = theme(axis.text.x = element_text(angle = 0))
      )+
      scale_fill_manual(values = cls)+
      xlab("Year of isolation")+
      ylab(ylab)+
      guides(fill=guide_legend(title=var, ncol = 1))+
      # no extra space around plot
      scale_y_continuous(expand = c(0,0))+
      switch(
        geom,
        geom_line = scale_x_continuous(expand = c(0,0)),
        geom_col = scale_x_discrete(expand = c(0,0)),
        geom_area = scale_x_continuous(expand = c(0,0))
      )
    if (type == "ratio" & show_count == TRUE) {
      g <- g + annotate("text",
                        label = round(df$yearcount, 0),
                        x = df$year,
                        y = 0.95,
                        size = count_size,
                        hjust = 1,
                        angle = 90)
    }
    g
  } else {
    df <- as.data.frame(table(df[[timevar]], df[[var]], df[[by]]))
    names(df) <- c("year", var, by, "count")
    df$year <- as.numeric(as.character(df$year))
    df$yearcount <- NA
    for (i in 1:nrow(df)) {
      df$yearcount[i] <- sum(
        df$count[which(df$year == df$year[i] & df[[by]] == df[[by]][i])])
    }
    df$mean_count <- NA
    df$mean_ratio <- NA
    for (i in 1:nrow(df)) {
      index <- which(
        df[[var]] == df[[var]][i] & 
          df$year <= df$year[i] &
          df$year > df$year[i]-window &
          df[[by]] == df[[by]][i])
      df$mean_count[i] <- mean(df$count[index])
      index_year <- which(df$year <= df$year[i] &
                            df$year > df$year[i]-window &
                            df[[by]] == df[[by]][i])
      df$mean_ratio[i] <- round(
        sum(df$count[index])/sum(df$count[index_year]), 3)
    }
    if (is.null(colors)) {
      if (length(unique(df[[by]])) > 1) {
        cls <- qualpalr::qualpal(
          length(unique(df[[by]])), colorspace = "pretty")$hex
      } else {
        cls = "black"
      }
    } else {
      cls <- colors
      if (length(cls) != length(unique(df[[by]]))) {
        stop(paste0("Variable '", by, "' requires ",
                    length(unique(df[[by]])), " separate colors."))
      }
    }
    df[[var]] <- as.logical(df[[var]])
    df <- df[which(df[[var]] == TRUE),]
    df$mean_ratio <- ifelse(df$yearcount == 0, 0, df$mean_ratio)
    ylab <- paste0(
      switch(type, count = "Number", ratio = "Ratio"),
      " of isolates"
    )
    if (window > 1) {
      ylab <- paste0(ylab, " (moving average, past ", window, " years)")
    }
    g <- ggplot(
      df,
      aes(year, switch(type, count = mean_count, ratio = mean_ratio)))+
      geom_line(aes(col = get(by)), alpha = 0.9)+
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 90))+
      scale_fill_manual(values = cls)+
      xlab("Year of isolation")+
      ylab(ylab)+
      guides(col=guide_legend(title=var))+
      # no extra space around plot
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0))
    g
  }
}

# TODO
# Function is too complex. There should be two functions, one where var is
# "character" and can have multiple states e.g. ST1, ST2, ST3, etc. and another
# where var is "logical" and can only have two states.

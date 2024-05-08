plot_ts_area <- function(
    df,
    var,
    timevar,
    minyear,
    window,
    show_count = FALSE,
    count_size = 8,
    colors = NULL,
    colorspace = "pretty"
    ) {
  if (!is.numeric(minyear)) {
    stop(paste("Variable '", minyear, "' must be numeric."))
  }
  df <- get_ts_df(df, var = var, timevar = timevar, window = window)
  df <- df[which(df[[timevar]] >= minyear),]
  cls_df <- get_colors(df, var = var, colors = colors, colorspace = colorspace)
  cls_vec <- cls_df$color
  names(cls_vec) <- cls_df[[var]]
  ylab <- "Ratio of isolates"
  if (window > 1) {
    ylab <- paste0(ylab, "\n (moving average, past ", window, " years)")
  }
  g <- ggplot(df, aes(get(timevar), ratio_window)) +
    geom_area(aes(fill = get(var)), alpha = 0.9)+
    scale_fill_manual(values = cls_vec)+
    xlab("Year of isolation")+
    ylab(ylab)+
    guides(fill=guide_legend(title=var))+
    # no extra space around plot
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) + 
    theme(
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(colour = "grey"),
      axis.text.x = element_text(angle = 0)
    )
  if (show_count == TRUE) {
    g <- g + 
      annotate("text",
              label = round(df$yearcount, 0),
              x = df[[timevar]],
              y = 0.95,
              size = count_size,
              hjust = 1,
              angle = 90)
  }
  return(g)
}
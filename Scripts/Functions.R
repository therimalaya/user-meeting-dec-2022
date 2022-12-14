label_vec <- function(vec, label_template = "{..2} (n={format(..1, big.mark = ',')})") {
  imap_chr(vec, ~glue::glue(label_template))
}
tbl2vec <- function(tbl, .dimnames = dimnames(tbl), .dim = dim(tbl), ...) {
  if (is.vector(tbl)) return(tbl)
  vec_names <- expand.grid(.dimnames) %>% interaction(...) %>% as.character()
  vec <- as.vector(tbl)
  names(vec) <- vec_names
  return(vec)
}

get_tstage_label <- function(Tstage) {
  if (missing(Tstage)) {
    message("Argument required:\t", appendLF = FALSE)
    return(NULL)
  }
  new_label <- c(
    "[0,0.5]" = paste0("0.0", "\U2013", "0.5 mm"),
    "(0.5,0.8]" = paste0(">0.5", "\U2013", "0.8 mm"),
    "(0.8,1]" = paste0(">0.8", "\U2013", "1 mm"),
    "(0.8,1)" = paste0(">0.8", "\U2013", "<1 mm"),
    "T1" = paste0("T1 (\U2264", "1.0 mm)"),
    "T1a" = "T1a",
    "T1b" = "T1b",
    "T2" = paste0("T2 (>1.0", "\U2013", "2.0 mm)"),
    "T3" = paste0("T3 (>2.0", "\U2013", "4.0 mm)"),
    "T4" = paste0("T4 (>4.0 mm)"),
    "Unspecified" = "Unspecified"
  )
  get_new_label <- function(x, map = new_label) {
    if (is.factor(x)) {
      y <- levels(x)
      y[y %in% names(map)] <- map[y[y %in% names(map)]]
      levels(x) <- y
    } else {
      x[x %in% names(map)] <- map[x[x %in% names(map)]]
    }
    return(x)
  }
  
  if (is.factor(Tstage)) {
    out <-  factor(
      get_new_label(Tstage),
      levels = get_new_label(levels(Tstage)),
      ordered = is.ordered(Tstage)
    )
  } else {
    out <- get_new_label(Tstage)
  }
  return(out)
}

plot_missing <- function(data, by = NULL, filter = "Sex != 'Overall'") {
  sex_lbl <- attr(data, "sex_label")
  by_lbl <- attr(data, "by_label")
  
  data <- data[eval(parse(text = filter))]
  
  Plot <- data[!(variable == "Ulceration" & DiagYear < 2000)] %>%
    ggplot(aes(DiagYear, Prop, color = variable)) +
    geom_line(na.rm = TRUE) +
    geom_point(shape = 21, fill = "whitesmoke", na.rm = TRUE) +
    geom_vline(xintercept = 2008, linetype = "dashed", color = "grey50") +
    facet_grid(reformulate("Sex", by), labeller = labeller(.rows = by_lbl, .cols = sex_lbl)) +
    theme(legend.position = "bottom") +
    labs(
      x = "Year of diagnosis",
      y = "Percentage of missing cases",
      color = "Missing"
    ) +
    scale_y_continuous(
      labels = function(x) round(x * 100), 
      limits = c(0, 0.85),
      breaks = scales::breaks_width(0.1),
      minor_breaks = scales::breaks_width(0.05)
    ) +
    scale_x_continuous(
      breaks = scales::breaks_width(5),
      minor_breaks = scales::breaks_width(1)
    ) +
    ggthemes::scale_color_stata()
  
  return(Plot)
}

plot_adj_rate <- function(
    data,
    x_var = "DiagYear",
    y_var = "adj.rate",
    type = "spline", 
    group_var = "Tstage", 
    size_var = NULL, 
    extra_label = TRUE,
    row_var = NULL, 
    col_var = NULL, 
    se = TRUE, 
    lrow = 1, 
    log = FALSE, 
    subset = NULL, ... 
) {
  row_var <- if (is.null(row_var)) "." else row_var
  col_var <- if (is.null(col_var)) "." else col_var
  data[, Size := if (!is.null(size_var)) data[[size_var]] else 1]
  group <- do.call(interaction, data[, ..group_var])
  data[, Group := stringr::str_replace(group, "-", "\U2013")]
  
  if (length(group_var) == 1 & "Tstage" %in% group_var) {
    data[, Group := forcats::fct_relabel(Group, get_tstage_label)]
  }
  if ("Sex" %in% names(data)) data[Sex == "Sum", Sex := "Overall"]
  
  if (length(group_var) == 1) {
    legend_title <- if (group_var == "Tstage") "T category" else group_var
  } else {
    legend_title <- paste0(group_var, collapse = ":")
  }
  legend_guide <- guide_legend(
    title = legend_title,
    title.position = "left",
    nrow = lrow
  )
  plot_caption <- "Note: Transparent: complete cases, Opaque: after multiple imputation"
  if (extra_label) {
    row_lbl <- col_lbl <- NULL
    if (row_var != ".") {
      .subset <- ifelse("Case" %in% names(data), "Sex != 'Overall'", "") 
      if (data[, uniqueN(Imp)] == 1) {
        if (data[, unique(Imp)] == "Complete") {
          row_lbl <- data[, xtabs(
            N ~ get(row_var), 
            data = .SD,
            subset = eval(parse(text = .subset))
          )]
        }
      } else {
        row_lbl <- data[
          Imp != "Complete", 
          xtabs(
            N / imp ~ get(row_var), 
            data = .SD, 
            subset = eval(parse(text = .subset))
          )
        ]  
      }
      row_lbl <- row_lbl %>% 
        round() %>% 
        tbl2vec() %>% 
        label_vec()
    }
    if (col_var != ".") {
      .subset <- ifelse("Case" %in% names(data), "Case == 'All Cases'", "")
      
      if (data[, uniqueN(Imp)] == 1) {
        if (data[, unique(Imp)] == "Complete") {
          col_lbl <- data[, xtabs(
            N / imp ~ get(col_var), 
            data = .SD,
            subset = eval(parse(text = .subset))
          )]
        }
      } else {
        col_lbl <- data[
          Imp != "Complete", 
          xtabs(N / imp ~ get(col_var), data = .SD,
                subset = eval(parse(text = .subset)))
        ]
      }
      col_lbl <- col_lbl %>% 
        round() %>% 
        tbl2vec() %>% 
        label_vec()
    }
    facet_label <- labeller(.rows = row_lbl, .cols = col_lbl)
    plot_caption <- glue::glue(
      plot_caption, "\n",
      "Note: n is the number of cases averaged",
      "over {data[Imp == 'Pooled', unique(imp)]} imputed datasets",
      sep = " "
    )
  }
  
  if (!is.null(subset)) {
    data <- data[eval(parse(text = subset))]
  }
  
  plot_frame <- expression({
    Plot <- data %>% 
      ggplot(aes(get(x_var), get(y_var), group = Group)) +
      ggh4x::facet_grid2(reformulate(col_var, row_var), labeller = facet_label, ...) +
      ggthemes::scale_color_stata() +
      ggthemes::scale_fill_stata() +
      guides(color = legend_guide, fill = legend_guide) +
      scale_y_continuous(
        trans = if (log) "log" else "identity",
        breaks = scales::breaks_extended(5),
        minor_breaks = scales::breaks_width(1)
      ) +
      scale_x_continuous(
        limits = c(1983, 2020),
        breaks = scales::breaks_width(5),
        guide = guide_axis(n.dodge = 2),
        labels = \(x) ifelse(as.numeric(x) > 2020, "", x)
      ) +
      scale_size_continuous(range = c(0.2, 2), guide = "none") +
      theme(
        plot.caption = element_text(hjust = 0, size = rel(1), lineheight = rel(1.15)),
        legend.position = "bottom",
        axis.text.x = element_text(vjust = 0.5),
        plot.caption.position = "plot"
      ) +
      labs(
        x = "Year of diagnosis",
        y = "Age-adjusted incidence rate\nper 100,000 person-year"
        # caption = plot_caption
      )
    if (!log) Plot <- Plot + expand_limits(y = 0)
  })
  points_complete_sharp <- expression({
    Plot <- Plot +
      geom_point(
        data = ~subset(.x, Imp %in% c(0, "Complete") & !is.na(Group)),
        aes(color = Group, size = Size),
        fill = "whitesmoke",
        na.rm = TRUE,
        shape = 21
      )
  })
  points_complete <- expression({
    Plot <- Plot +
      geom_point(
        data = ~subset(.x, Imp %in% c(0, "Complete") & !is.na(Group)),
        aes(color = Group, size = Size),
        alpha = 0.25, 
        fill = "whitesmoke",
        na.rm = TRUE
      )
  })
  points_pooled <- expression({
    Plot <- Plot +
      geom_point(
        data = ~subset(.x, !Imp %in% c(0, "Complete")),
        aes(color = Group, size = Size),
        shape = 21,
        fill = "whitesmoke"
      )
  })
  splines <- expression({
    Plot <- Plot +
      geom_line(
        data = ~subset(.x, Imp %in% c(0, "Complete")),
        aes(y = adj.rate_spl, color = Group, group = Group),
        alpha = 0.25,
        na.rm = TRUE
      ) +
      geom_line(
        data = ~subset(.x, !Imp %in% c(0, "Complete")),
        aes(y = adj.rate_spl, color = Group, group = Group),
        na.rm = TRUE
      )
  })
  segmented <- expression({
    Plot <- Plot +
      geom_line(
        data = ~subset(.x, Imp %in% c(0, "Complete")),
        aes(y = fit, color = Group, group = Group),
        alpha = 0.25,
        na.rm = TRUE
      ) +
      geom_line(
        data = ~subset(.x, !Imp %in% c(0, "Complete")),
        aes(y = fit, color = Group, group = Group)
      )
  })
  adj_lines <- expression({
    Plot <- Plot +
      geom_line(
        data = ~subset(.x, !Imp %in% c(0, "Complete")),
        aes(color = Group)
      )
  })
  adj_lines_complete <- expression({
    Plot <- Plot +
      geom_line(
        data = ~subset(.x, Imp %in% c(0, "Complete")),
        aes(color = Group)
      )
  })
  lines <- switch(
    type,
    spline = splines,
    segmented = segmented,
    adj_lines
  )
  ribbons <- expression({
    Plot <- Plot +
      geom_ribbon(
        aes(ymax = uci_spl, ymin = lci_spl, fill = Group, group = Group),
        data = ~subset(.x, !Imp %in% c(0, "Complete")),
        alpha = 0.15, 
        color = NA
      )
  })
  
  eval(plot_frame)
  if (type == "segmented") {
    eval(points_complete)
    eval(points_pooled)
    eval(lines)
  } else if (type == "spline") {
    if (se) eval(ribbons)
    eval(points_complete)
    eval(points_pooled)
    eval(lines)
  } else {
    if (data[, uniqueN(Imp)] == 1) {
      if (data[, unique(Imp)] %in% c(0, "Complete")) {
        eval(adj_lines_complete)
        eval(points_complete_sharp)
      }
    } else {
      eval(lines)
      eval(points_complete)
      eval(points_pooled)
    }
  }
  return(Plot)
}

plot_thist <- function(plot_data, type = c("freq", "density"), ...) {
  library(ggplot2)
  type <- match.arg(type)
  
  group_var <- setdiff(names(plot_data), c("Cuts", "Thickness"))
  cuts <- plot_data[, levels(Cuts)] %>% 
    str_remove_all("[\\(\\)\\[\\]]") %>% 
    str_split(",") %>% 
    unlist() %>% 
    unique() %>% 
    as.numeric()
  
  plt <- ggplot(plot_data, aes(Thickness)) +
    geom_histogram(
      aes(fill = Cuts, y = if (type == "freq") ..count.. else ..density..), 
      binwidth = 0.1, 
      color = "black", 
      size = 0.1
    ) +
    geom_vline(xintercept = first(cuts, -1) + 0.05, size = 0.25, lty = 2) +
    scale_y_continuous(
      expand = expansion(c(0, 0.05))
    ) +
    ggthemes::scale_color_stata() +
    ggthemes::scale_fill_stata() +
    labs(
      x = "Tumour thickness", 
      y = "Number of cases",
      fill = "Cut-points",
      color = "Cut-points"
    ) +
    ggthemes::theme_few() +
    theme(
      legend.position = c(1, 1),
      legend.justification = c(1.05, 1.05),
      legend.direction = "horizontal",
      panel.grid = element_line(color = "#f0f0f0")
    )
  if (length(group_var)) {
    plt <- plt + facet_grid(rows = vars(get(group_var)), ...)
  }
  return(plt)
}

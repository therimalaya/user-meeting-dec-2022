## ----------------------------------------------------------------
source("Scripts/Functions.R")
library(data.table)
library(ggplot2)
library(stringr)
library(purrr)
library(ggh4x)
library(gt)

ggplot2::theme_set(ggthemes::theme_few())
theme_update(panel.grid = element_line(color = "#f0f0f0"))

## -- Local functions ----------------------------------------
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

plot_ns <- function(ns_fit) {
  ns_fit <- ns_fit[!is.na(tstage)]
  groups <- attr(ns_fit, "groups")
  x_scale <- function(plot) {
    plot + scale_x_continuous(
      breaks = scales::breaks_width(2 * 365.241),
      labels = function(x) x/365.241,
      limits = c(0, 30 * 365.241)
    )
  }
  color_scale <- function(plot) {
    plot + 
      ggthemes::scale_color_stata(
        labels = get_tstage_label,
        na.value = "grey"
      ) + ggthemes::scale_fill_stata(
        labels = get_tstage_label,
        na.value = "grey"
      )
  }
  thm <- function(plot) {
    plot + ggthemes::theme_few() +
      theme(
        panel.grid = element_line(color = "#f0f0f0"),
        legend.position = "bottom"
      )
  }
  lbl <- function(plot) {
    plot + labs(
      x = "Survival time (years)",
      y = "Relative survival",
      color = "T category",
      fill = "T category",
      linetype = "Period",
      title = "Relative survival percentage",
      subtitle = paste("By", knitr::combine_words(stringr::str_to_title(groups)))
    )
  }
  anot <- function(plot) {
    plot + annotate(
      geom = "text",
      x = c(1, 5, 10, 15, 25) * 365.241,
      y = Inf,
      label = paste(c(1, 5, 10, 15, 25), "-year\nfollow up"),
      angle = 90,
      hjust = 1,
      vjust = 0
    )
  }
  geoms <- function(plot) {
    plot + geom_step(
      aes(color = tstage),
      na.rm = TRUE, 
    ) + geom_ribbon(
      aes(fill = tstage),
      color = NA, alpha = 0.15
    ) + geom_vline(
      xintercept = c(1, 5, 10, 15, 25) * 365.241, 
      linetype = "dashed", color = "#afafaf"
    )
  }
  
  if ("year_cat" %in% groups) {
    plot <- ggplot(ns_fit, aes(
      time, surv,
      ymin = lower,
      ymax = upper,
      group = interaction(year_cat, tstage),
      linetype = year_cat
    ))
  } else {
    plot <- ggplot(ns_fit, aes(
      time, surv,
      ymin = lower,
      ymax = upper,
      group = tstage
    ))
  }
  plot <- plot %>% geoms %>% thm %>% x_scale %>% color_scale %>% lbl # %>% anot
  
  if ("sex" %in% groups) {
    plot <- plot + facet_grid(
      cols = vars(sex), 
      labeller = labeller(
        sex = stringr::str_to_title
      )
    )
  }
  
  plot + expand_limits(y = c(0, 1.25))
}



## ----------------------------------------------------------------
sgmt_plt_data <- readRDS("Data/adj-rate-sgmt-by-sex-tstage.Rds")[, size := uci - lci]


## ----------------------------------------------------------------
adj_rate_plot <- sgmt_plt_data[Imp == "Complete"] %>%
  plot_adj_rate(
    type = "adj-rate",
    group_var = "Tstage",
    col_var = "Sex",
    size_var = "size"
  ) + coord_cartesian(ylim = c(0, 25)) +
  labs(title = "Age-adjusted incidence rate using only complete cases")

ggsave(
  plot = adj_rate_plot,
  filename = here::here("Images", "adj-rate-plot.svg"),
  width = 12,
  height = 6,
  scale = 0.85
)


## ----------------------------------------------------------------
imp_adj_rate_plot <- sgmt_plt_data %>%
  plot_adj_rate(
    type = "adj-rate",
    group_var = "Tstage",
    col_var = "Sex",
    size_var = "size"
  ) + coord_cartesian(ylim = c(0, 25)) +
  labs(title = "Age-adjusted incidence rate using imputed data")

ggsave(
  plot = imp_adj_rate_plot,
  filename = here::here("Images", "imp-adj-rate-plot.svg"),
  width = 12,
  height = 6,
  scale = 0.85
)

## ----------------------------------------------------------------
sgmt_adj_rate_plot <- sgmt_plt_data %>%
  plot_adj_rate(
    type = "segmented",
    group_var = "Tstage",
    col_var = "Sex",
    size_var = "size"
  ) + coord_cartesian(ylim = c(0, 25)) +
  labs(title = "Age-adjusted incidence rate using imputed data with segmented regression")

ggsave(
  plot = sgmt_adj_rate_plot,
  filename = here::here("Images", "sgmt-adj-rate-plot.svg"),
  width = 12,
  height = 6,
  scale = 0.85
)

## ----------------------------------------------------------------
spl_plt_data <- readRDS("Data/adj-rate-spl-by-sex-tstage.Rds")[, size := uci - lci]

## ----------------------------------------------------------------
spl_adj_rate_plot <- spl_plt_data %>%
  plot_adj_rate(
    type = "spline",
    group_var = "Tstage",
    col_var = "Sex",
    size_var = "size"
  ) + coord_cartesian(ylim = c(0, 25)) +
  labs(title = "Age-adjusted incidence rate using imputed data with splines")

ggsave(
  plot = spl_adj_rate_plot,
  filename = here::here("Images", "spl-adj-rate-plot.svg"),
  width = 12,
  height = 6,
  scale = 0.85
)

## ----------------------------------------------------------------
tbls <- readRDS("Data/aapc-by-sex-tstage.Rds")
color_row <- function(
    tbl, 
    .cols = gt::everything(), 
    .rows = gt::matches("T1"), 
    .color = "royalblue",
    .color_stub = .color,
    ...) {
  tbl %>% 
    gt::tab_style(
      locations = gt::cells_body(columns = .cols, rows = .rows),
      style = cell_text(color = .color, ...)
    ) %>% 
    gt::tab_style(
      locations = gt::cells_stub(rows = .rows),
      style = cell_text(color = .color_stub, ...)
    )
}

apc_tables <- imap(tbls, function(.x, .y) {
  out <- .x %>% 
    gt::fmt_missing(gt::everything(), missing_text = "-") %>% 
    color_row(weight = "bold")
  if (.y == "Trend 2") {
    out <- out %>% 
      color_row(.cols = gt::matches("Trend 2"),
                .color = "firebrick", .color_stub = "royalblue")
  }
  return(out)
})

iwalk(apc_tables, function(.tbl, .name) {
  gt::gtsave(.tbl, here::here(
    "Images", 
    glue::glue("apc-table-{tolower(.name)}.html")
  ))
})

saveRDS(apc_tables, file = here::here("Data/aapc-table-by-sex-tstage.Rds"))

## ----------------------------------------------------------------
site_data <- readRDS("Data/adj-rate-spl-by-sex-tstage-site.Rds")
type_data <- readRDS("Data/adj-rate-spl-by-sex-tstage-type.Rds")


## ----------------------------------------------------------------
site_data[AnatomicSite != "Other"] %>% 
  setnames("DiagYear", "Year", skip_absent = TRUE) %>% 
  split(by = "Sex") %>% 
  iwalk(function(.data, .name) {
    out <- .data %>% 
        plot_adj_rate(
          col_var = "AnatomicSite", 
          row_var = "Sex", 
          x_var = "Year", 
          scales = "free_y", 
          independent = "y"
        ) +
        ggh4x::facetted_pos_scales(
          y = list(
            AnatomicSite == "Trunk" ~ scale_y_continuous(
              breaks = scales::breaks_width(5),
              limits = c(0, 15)
            ),
            AnatomicSite == "Head and neck" ~ scale_y_continuous(
              breaks = scales::breaks_width(1),
              limits = c(0, 3)
            ),
            AnatomicSite == "Upper limbs" ~ scale_y_continuous(
              breaks = scales::breaks_width(1),
              limits = c(0, 4)
            ),
            AnatomicSite == "Lower limbs" ~ scale_y_continuous(
              breaks = scales::breaks_width(2),
              limits = c(0, 8)
            )
          )
        ) +
      labs(
        title = glue::glue("Age-adjusted incidence rate in {tolower(.name)}"),
        subtitle = "By T category and anatomic site"
      ) +
      theme(
        plot.subtitle = element_text(margin = margin(0, 0, 15, 0))
      )
    ggsave(
      plot = out,
      filename = here::here(
        "Images", 
        glue::glue("adj-rate-plot-by-site-{tolower(.name)}.svg")
      ),
      width = 13.5,
      height = 5.5,
      scale = 0.9
    )
  })


## ----------------------------------------------------------------
type_data[MelanomaType != "Other"] %>% 
  setnames("DiagYear", "Year", skip_absent = TRUE) %>% 
  split(by = "Sex") %>% 
  iwalk(function(.data, .name) {
    out <- .data %>% 
      plot_adj_rate(
        col_var = "MelanomaType", 
        row_var = "Sex", 
        x_var = "Year", 
        scales = "free_y", 
        independent = "y"
      ) +
      ggh4x::facetted_pos_scales(
        y = list(
          MelanomaType == "Superficial spreading" ~ scale_y_continuous(
            breaks = scales::breaks_width(5),
            limits = c(0, 25)
          ),
          MelanomaType == "Nodular" ~ scale_y_continuous(
            breaks = scales::breaks_width(1),
            limits = c(0, 4)
          ),
          MelanomaType == "Lentigo maligna" ~ scale_y_continuous(
            breaks = scales::breaks_width(0.5),
            limits = c(0, 1.5)
          )
        )
      ) +
      labs(
        title = glue::glue("Age-adjusted incidence rate in {tolower(.name)}"),
        subtitle = "By T category and melanoma subtype"
      ) +
      theme(
        plot.subtitle = element_text(margin = margin(0, 0, 10, 0))
      )
    ggsave(
      plot = out,
      filename = here::here(
        "Images", 
        glue::glue("adj-rate-plot-by-type-{tolower(.name)}.svg")
      ),
      width = 13.5,
      height = 5.5,
      scale = 0.9
    )
  })

## -- Missing Plot ----------------------------------------
missing_data <- readRDS("Data/missing-data.Rds")
missing_plot <- plot_missing(missing_data, filter = "Sex == 'Overall'")
ggsave(
  plot = missing_plot,
  filename = here::here(
    "Images", 
    glue::glue("missing-plot.svg")
  ),
  width = 6,
  height = 6,
  scale = 0.75
)

## -- Population Pyramid ----------------------------------------
data(stdpop18, package = "popEpi")
pop_data <- copy(stdpop18)
pop_data[, nordic := NULL]
pop_data[, id := .I]

pop_plot <- melt.data.table(
  pop_data,
  id.vars = c("agegroup", "id")
)[variable == "world", value := value * -1][order(id)] %>% 
  ggplot(aes(value, reorder(agegroup, id), fill = variable)) +
  geom_col(width = 1, color = "#f0f0f0", size = 0.2) +
  theme_minimal() +
  ggthemes::scale_fill_economist(
    labels = stringr::str_to_title
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_line(color = "red")
  ) +
  scale_x_continuous(
    labels = abs,
    expand = expansion()
  ) +
  scale_y_discrete(
    expand = expansion()
  ) +
  labs(
    x = "Number of person",
    y = NULL,
    fill = "Population"
  ) +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    label = "World",
    hjust = -0.5,
    vjust = 1.5
  ) +
  annotate(
    geom = "text",
    x = Inf,
    y = Inf,
    label = "Europe",
    hjust = 1.5,
    vjust = 1.5
  )

ggsave(
  pop_plot, 
  filename = "Images/pop-plot.svg", 
  width = 7, 
  height = 7,
  scale = 0.5
)


## -- Plot the historam ----------------------------------------
plot_data_thist <- readRDS("Data/thist.Rds") %>% 
  .[, DiagYear10 := factor(DiagYear10, ordered = FALSE)]
density_df <- plot_data_thist[, with(
  density.default(Thickness, bw = 0.035), 
  data.table(Thickness = x, Density = y)
), by = .(DiagYear10)] %>% 
  rbind(
    plot_data_thist[, with(
      density.default(Thickness, bw = 0.035), 
      data.table(Thickness = x, Density = y)
    )][, DiagYear10 := "Overall (1983–2019)"]
  ) %>% 
  .[, Cuts := cut.default(
    Thickness - 0.05, 
    breaks = c(0, 0.8, 1, 2, 4, 8), 
    include.lowest = TRUE, 
    right = TRUE
  )] %>% .[!is.na(Cuts)]

plts <- plot_data_thist %>% 
  split(by = "DiagYear10") %>% 
  append(list(Overall = plot_data_thist)) %>% 
  imap(function(dta, name) {
    dta[, DiagYear10 := name]
    if (name == "Overall") {
      dta[, DiagYear10 := "Overall (1983–2019)"]
    }
    plt <- ggplot(dta, aes(Thickness)) +
      geom_histogram(
        aes(y = after_stat(count / 30)), 
        binwidth = 0.1,
        fill = "#f0f0f0",
        color = "#e0e0e0",
        linewidth = 0.5,
        na.rm = TRUE
      ) +
      geom_histogram(
        aes(fill = Cuts, y = after_stat(density)),
        binwidth = 0.1,
        alpha = 0.15,
        na.rm = TRUE
      ) +
      geom_density(
        aes(color = Cuts), 
        linewidth = 0.3, 
        linetype = "dashed", 
        na.rm = TRUE
      ) +
      geom_line(
        data = density_df[grepl(name, DiagYear10)],
        aes(y = Density, x = Thickness, color = Cuts),
        na.rm = TRUE
      ) +
      facet_grid(cols = vars(DiagYear10)) +
      scale_x_continuous(
        breaks = scales::breaks_width(0.5),
        minor_breaks = scales::breaks_width(0.1),
        expand = expansion(),
        limits = c(0, 8)
      ) +
      scale_y_continuous(
        name = "Density",
        breaks = c(0, 0.5, 1, 3, 5, 10, 20, 30, 50, 100, 200),
        labels = c(0, 0.5, 1, 3, 5, 10, 20, 30, 50, 100, 200),
        minor_breaks = scales::breaks_log(50)(1:220, n = 100),
        expand = expansion(add = 0.1),
        trans = "log1p",
        position = "left",
        limits = c(0, 200),
        sec.axis = sec_axis(
          name = "Count",
          trans = function(x) x * 30,
          breaks = c(0, 0.5, 1, 3, 5, 10, 20, 30, 50, 100, 200) * 30
        )
      ) +
      ggthemes::theme_few() +
      theme(
        legend.position = "bottom",
        panel.grid = element_line(color = "#f0f0f0")
      ) +
      labs(
        x = "Tumour thickness (mm)",
        title = "Distribution of tumour thickness",
        subtitle = "Background histogram represents the count"
      )
    
    return(plt)
  })

iwalk(plts, function(plt, name) {
  ggsave(
    plt, 
    filename = glue::glue("Images/thist-plot-{name}.svg"), 
    width = 12,
    height = 6,
    scale = 0.85
  )
})

## -- Survival plots ----------------------------------------
ns_fits <- readRDS("Data/ns-fits.Rds")
ns_plots <- map(ns_fits, plot_ns)
iwalk(ns_plots, function(plot, name) {
  ggsave(
    plot, 
    filename = here::here(
      glue::glue("Images/ns-{name}.svg")
    ),
    width = 13.5,
    height = 7.5,
    scale = 0.7
  )
})

## -- Survival trend ----------------------------------------
fit_per_year_df <- readRDS("Data/rs-trend.Rds")
rs_trend_plot <- fit_per_year_df[time %in% (c(5, 10, 15) * 365.241)] %>% 
  copy() %>% 
  .[, diag_year := as.integer(diag_year)] %>% 
  .[, time := time / 365.241] %>%
  .[!(time == 5 & diag_year > 2014)] %>% 
  .[!(time == 10 & diag_year > 2009)] %>% 
  .[!(time == 15 & diag_year > 2004)] %>% 
  ggplot(aes(diag_year, surv, color = tstage, group = interaction(tstage, time))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = tstage), color = NA, alpha = 0.15, na.rm = TRUE) +
  geom_line(na.rm = TRUE) +
  facet_grid(rows = vars(sex), cols = vars(time),
             labeller = labeller(
               time = function(x) paste0("Follow up: ", x),
               sex = stringr::str_to_title
             )) +
  scale_x_continuous(breaks = scales::breaks_width(5), minor_breaks = scales::breaks_width(1)) +
  scale_y_continuous(breaks = scales::breaks_width(0.2), minor_breaks = scales::breaks_width(0.1),
                     labels = scales::label_percent()) +
  ggthemes::scale_color_stata(
    na.value = "grey",
    labels = function(x) fifelse(is.na(x), "Unspecified", get_tstage_label(x))
  ) +
  ggthemes::scale_fill_stata(
    na.value = "grey",
    labels = function(x) fifelse(is.na(x), "Unspecified", get_tstage_label(x))
  ) +
  expand_limits(y = c(0, 1)) +
  ggthemes::theme_few() +
  theme(panel.grid = element_line(color = "#f0f0f0"),
        legend.position = "bottom") +
  labs(
    title = "Relative survival by tumour thickness over time",
    x = "Year of diagnosis",
    y = "Relative survival",
    color = "T category",
    fill = "T category"
  )

ggsave(
  rs_trend_plot, 
  filename = here::here(
    "Images/rs-trend.svg"
  ),
  width = 13.5,
  height = 7.5,
  scale = 0.7
)


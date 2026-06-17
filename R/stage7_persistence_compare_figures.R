# Figure construction for 7.3_persist_cmp.Rmd.
# Requires the 7.3 comparison outputs and creates the combined comparison figure object.

library(ggplot2)

label_method <- function(x) {
  out <- as.character(x)
  out[out == "pipe"] <- "Persistence-based"
  out[out == "rank"] <- rank_label
  out
}

method_cols <- setNames(
  c("#2E78AF", "#B7561D"),
  c("Persistence-based", rank_label)
)

persist_breaks_shared <- seq(0, 1, by = 0.25)
persist_labels_shared <- scales::label_number(accuracy = 0.01)

stage_breaks_shared <- function(x) pretty(x, n = 5)

xlab_plot <- "Cells removed (%)"

theme_persist <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(colour = "grey20"),
      axis.text = element_text(colour = "grey20"),
      axis.line = element_line(colour = "grey20", linewidth = 0.45),
      axis.ticks = element_line(colour = "grey20", linewidth = 0.45),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

# ============================================================
# 04_plot_mean - panel (a) for combined figure
# ============================================================

if (!exists("stage_sum")) stage_sum <- fread(out_stage)

plot_dt_mean <- copy(stage_sum)[
  ,
  method_lab := label_method(method)
][order(method_lab, stage_order)]

x_breaks_mean <- stage_breaks_shared(plot_dt_mean$x)

p_mean <- ggplot(
  plot_dt_mean,
  aes(x = x, y = mean_persist, color = method_lab)
) +
  geom_line(linewidth = 0.95, lineend = "round") +
  geom_point(size = 1.8) +
  scale_color_manual(values = method_cols, name = NULL) +
  scale_x_continuous(
    breaks = x_breaks_mean,
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = persist_breaks_shared,
    labels = persist_labels_shared,
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = xlab_plot,
    y = "Mean persistence"
  ) +
  theme_persist(base_size = 11) +
  theme(
    legend.position = "none"
  )

p_mean_leg <- p_mean +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.margin = margin(t = 2)
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(linewidth = 1.15, shape = 16, size = 3)
    )
  )

# ============================================================
# 04b_focal_species - single focal species used across panels
# ============================================================

focal_name <- if (!is.null(params$focus_species) && length(params$focus_species) >= 1L) {
  as.character(params$focus_species[[1]])
} else {
  "Fossa fossana" # "Hapalemur alaotrensis", "Pterocles personatus", "Galidia elegans", "Fossa fossana"
}
focal_id <- sid(focal_name)

# ============================================================
# 05_plot_focus_panel - panel (b) for combined figure
# focal species persistence across the removal sequence
# ============================================================

if (!exists("sp_long")) sp_long <- fread(out_sp)

dd_focus <- sp_long[species == focal_id]
if (!nrow(dd_focus)) {
  stop(paste0("Species not found in sp_long: ", focal_name), call. = FALSE)
}

plot_dt_focus <- copy(dd_focus)[
  ,
  method_lab := label_method(method)
][order(method_lab, stage_order)]

x_breaks_focus <- stage_breaks_shared(plot_dt_focus$x)

p_focus <- ggplot(
  plot_dt_focus,
  aes(x = x, y = sp_persist, color = method_lab)
) +
  geom_line(linewidth = 0.95, lineend = "round") +
  geom_point(size = 1.8) +
  scale_color_manual(values = method_cols, name = NULL) +
  scale_x_continuous(
    breaks = x_breaks_focus,
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = persist_breaks_shared,
    labels = persist_labels_shared,
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = xlab_plot,
    y = "Species persistence"
  ) +
  theme_persist(base_size = 11) +
    theme(
      legend.position = "none",
      axis.ticks.x = element_line(linewidth = 0.45, colour = "grey20"),
      axis.ticks.length.x = grid::unit(2.2, "pt"),
      plot.caption = element_text(
        hjust = 0.5,
        colour = "grey20",
        margin = margin(t = 4)
      ),
      plot.margin = margin(0, 0, 0, 0)
    )

# ============================================================
# 06_top10_auc_signed - top 10 species by signed area between pipe vs rank
# Uses trapezoidal integration of pipe - rank across x
# Positive values mean pipe > rank on net
# Negative values mean rank > pipe on net
# ============================================================

if (!exists("sp_long")) sp_long <- fread(out_sp)

# ---- species metadata: taxon + body mass + diet/density/dispersal ----
species_header <- names(fread(species_csv, nrows = 0))

need_meta_cols <- c(
  "scientificName",
  "className",
  "BodyMass.Value",
  "Diet",
  "density",
  "dispersal_dist"
)

missing_meta_cols <- setdiff(need_meta_cols, species_header)
if (length(missing_meta_cols)) {
  stop(
    "species_table.csv is missing column(s): ",
    paste(missing_meta_cols, collapse = ", "),
    call. = FALSE
  )
}

sp_meta <- fread(species_csv, select = need_meta_cols)

sp_meta[, scientificName := trimws(as.character(scientificName))]
sp_meta[, species := sid(scientificName)]

sp_meta[, taxon := fifelse(
  tolower(trimws(className)) == "mammalia", "Mammal",
  fifelse(tolower(trimws(className)) == "aves", "Bird", as.character(className))
)]

sp_meta[, body_mass_g := suppressWarnings(as.numeric(BodyMass.Value))]
sp_meta[, density := suppressWarnings(as.numeric(density))]
sp_meta[, dispersal_dist_km := suppressWarnings(as.numeric(dispersal_dist))]
sp_meta[, diet := as.character(Diet)]

sp_meta <- unique(
  sp_meta[, .(
    species,
    taxon,
    body_mass_g,
    diet,
    density,
    dispersal_dist_km
  )]
)

# Keep only the columns needed and reshape so pipe/rank are side by side
auc_dt <- dcast(
  sp_long[, .(species, scientificName, stage, stage_order, x, method, sp_persist)],
  species + scientificName + stage + stage_order + x ~ method,
  value.var = "sp_persist"
)

# Require both methods to be present
auc_dt <- auc_dt[is.finite(pipe) & is.finite(rank)]
setorder(auc_dt, species, stage_order)

# Signed gap between trajectories:
# positive = persistence-based method has higher persistence
# negative = area-ranked method has higher persistence
auc_dt[, signed_gap_pipe_minus_rank := pipe - rank]

# Trapezoidal signed AUC for each species
top10_auc_signed <- auc_dt[, {
  if (.N < 2L) {
    .(auc_signed_pipe_minus_rank = NA_real_)
  } else {
    dx <- diff(x)
    avg_signed_gap <- (
      signed_gap_pipe_minus_rank[-.N] +
        signed_gap_pipe_minus_rank[-1L]
    ) / 2

    .(
      auc_signed_pipe_minus_rank = sum(
        dx * avg_signed_gap,
        na.rm = TRUE
      )
    )
  }
}, by = .(species, scientificName)]

# Add taxon and body mass
top10_auc_signed <- merge(
  top10_auc_signed,
  sp_meta,
  by = "species",
  all.x = TRUE,
  sort = FALSE
)

top10_auc_signed <- top10_auc_signed[
  is.finite(auc_signed_pipe_minus_rank)
][order(-auc_signed_pipe_minus_rank)]

top10_auc_signed[, rank := seq_len(.N)]

# ============================================================
# 07_plot_species_total_area - panel (c) for combined figure
# Simplified stacked bars:
# - show only every other stage
# - color segments by current PU persistence probability
# - use DIFFERENT fill scales for the two methods
# - within each bar, larger / more persistent PUs are stacked lower
# - draw all PU segments with the same border symbology
# - add a light outline around each full bar
# ============================================================

if (!exists("pu_long")) pu_long <- fread(out_pu)

dd <- pu_long[
  species == focal_id & method %in% c("pipe", "rank"),
  .(method, stage, stage_order, pu_id, pu_area_km2, P_pu)
]

if (!nrow(dd)) {
  stop(paste0("Species not found in pu_long: ", focal_name), call. = FALSE)
}

# ---- current PU area + persistence by stage ----
dd_sum <- dd[
  ,
  .(
    pu_area_km2 = sum(pu_area_km2, na.rm = TRUE),
    current_P_pu = max(P_pu, na.rm = TRUE)
  ),
  by = .(method, stage, stage_order, pu_id)
]

# only every other stage
plot_stages <- sort(unique(dd_sum$stage[dd_sum$stage %% 2L == 0L]))
if (!length(plot_stages)) plot_stages <- sort(unique(dd_sum$stage))

dd_sum <- dd_sum[stage %in% plot_stages]

# use all PUs observed for this focal species across the sequence
all_pu_ids <- sort(unique(dd$pu_id))

stage_key <- unique(dd_sum[, .(method, stage, stage_order)])
setorder(stage_key, method, stage_order)

grid_dt <- CJ(
  method = unique(dd_sum$method),
  stage  = plot_stages,
  pu_id  = all_pu_ids,
  unique = TRUE
)

grid_dt <- merge(
  grid_dt,
  stage_key,
  by = c("method", "stage"),
  all.x = TRUE,
  sort = FALSE
)

dd_plot <- merge(
  grid_dt,
  dd_sum,
  by = c("method", "stage", "stage_order", "pu_id"),
  all.x = TRUE,
  sort = FALSE
)

dd_plot[is.na(pu_area_km2), pu_area_km2 := 0]
dd_plot[is.na(current_P_pu), current_P_pu := 0]

# bottom of stack = larger and/or more persistent PUs
dd_plot[
  ,
  stack_score := {
    area_max <- max(pu_area_km2, na.rm = TRUE)
    area_std <- if (is.finite(area_max) && area_max > 0) pu_area_km2 / area_max else 0
    0.5 * area_std + 0.5 * current_P_pu
  },
  by = .(method, stage, stage_order)
]

setorder(dd_plot, method, stage_order, -stack_score, -pu_area_km2, -current_P_pu, pu_id)

dd_plot[
  ,
  `:=`(
    ymin = c(0, cumsum(pu_area_km2)[-.N]),
    ymax = cumsum(pu_area_km2)
  ),
  by = .(method, stage, stage_order)
]

stage_show <- dd_plot[order(stage_order), unique(stage)]
dd_plot[, stage_index := match(stage, stage_show)]
dd_plot[, x_plot := stage_index + fifelse(method == "pipe", -0.22, 0.22)]

bar_half_width <- 0.19
dd_plot[, `:=`(
  xmin = x_plot - bar_half_width,
  xmax = x_plot + bar_half_width
)]

bar_outline <- dd_plot[
  pu_area_km2 > 0,
  .(
    xmin = min(xmin),
    xmax = max(xmax),
    ymin = 0,
    ymax = max(ymax)
  ),
  by = .(method, stage, stage_order)
]

y_max <- dd_plot[, max(ymax, na.rm = TRUE)]

# method-specific persistence palettes
pipe_fill_cols <- c(
  "#F2F4F7",
  "#DCE7F0",
  "#BCD4E6",
  "#8FB8D8",
  "#5F99C6",
  "#2E78AF"
)

rank_fill_cols <- c(
  "#FCF3ED",
  "#F8DDCF",
  "#F2C0A6",
  "#EA9A72",
  "#D9723F",
  "#B7561D"
)

pipe_pal <- scales::col_numeric(
  palette = pipe_fill_cols,
  domain = c(0, 1),
  na.color = "#F2F4F7"
)

rank_pal <- scales::col_numeric(
  palette = rank_fill_cols,
  domain = c(0, 1),
  na.color = "#FCF3ED"
)

dd_plot[
  ,
  fill_col := ifelse(
    method == "pipe",
    pipe_pal(current_P_pu),
    rank_pal(current_P_pu)
  )
]

p_area <- ggplot() +
  geom_rect(
    data = dd_plot[pu_area_km2 > 0],
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill_col
    ),
    color = "white",
    linewidth = 0.22,
    linejoin = "mitre"
  ) +
  geom_rect(
    data = bar_outline,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = NA,
    color = "grey78",
    linewidth = 0.35,
    linejoin = "mitre"
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = seq_along(stage_show),
    labels = as.character(stage_show),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(
    breaks = pretty(c(0, y_max), n = 4),
    labels = scales::label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.03))
  ) +
  coord_cartesian(ylim = c(0, y_max * 1.05), clip = "off") +
  labs(
    x = "Stage",
    y = expression("Total PU area (km"^2*")")
  ) +
  theme_persist(base_size = 11) +
  theme(
    legend.position = "none"
  )

# ---- legend donor: compact stacked colour bars with shared label below ----
legend_x <- seq(0, 1, length.out = 1024)
pipe_raster <- as.raster(matrix(pipe_pal(legend_x), nrow = 1))
rank_raster <- as.raster(matrix(rank_pal(legend_x), nrow = 1))

p_area_leg <- ggplot() +
  annotation_raster(
    raster = pipe_raster,
    xmin = 0, xmax = 1,
    ymin = 0.61, ymax = 0.73
  ) +
  annotation_raster(
    raster = rank_raster,
    xmin = 0, xmax = 1,
    ymin = 0.45, ymax = 0.57
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = persist_breaks_shared,
    labels = persist_labels_shared,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.38, 0.77),
    breaks = NULL,
    labels = NULL,
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = NULL,
    caption = "PU persistence probability"
  ) +
  theme_classic(base_size = 10.1) +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour = "grey20"),
    axis.line.x = element_line(linewidth = 0.45, colour = "grey20"),
    axis.ticks.x = element_line(linewidth = 0.45, colour = "grey20"),
    axis.ticks.length.x = grid::unit(2.2, "pt"),
    plot.caption = element_text(
      hjust = 0.5,
      colour = "grey20",
      margin = margin(t = 4)
    ),
    plot.margin = margin(0, 0, 0, 0)
  )

# ============================================================
# 08_plot_species_persist_vs_total_area - panel (d) for combined figure
# Extend each method trajectory to the origin (0, 0)
# ============================================================

if (!exists("pu_long")) pu_long <- fread(out_pu)
if (!exists("sp_long")) sp_long <- fread(out_sp)

area_dt <- pu_long[
  species == focal_id & method %in% c("pipe", "rank"),
  .(total_area_km2 = sum(pu_area_km2, na.rm = TRUE)),
  by = .(method, stage, stage_order)
]

persist_dt <- sp_long[
  species == focal_id & method %in% c("pipe", "rank"),
  .(method, stage, stage_order, sp_persist)
]

plot_dt_area <- merge(
  area_dt,
  persist_dt,
  by = c("method", "stage", "stage_order"),
  all = FALSE,
  sort = FALSE
)

if (!nrow(plot_dt_area)) {
  stop(paste0("Species not found in pu_long/sp_long: ", focal_name), call. = FALSE)
}

plot_dt_area[, method_lab := label_method(method)]

# add an origin row for each method so the line reaches (0, 0)
origin_dt <- unique(plot_dt_area[, .(method)])
origin_dt[, `:=`(
  stage = NA_integer_,
  stage_order = Inf,
  total_area_km2 = 0,
  sp_persist = 0,
  method_lab = label_method(method)
)]

plot_dt_area_path <- rbindlist(
  list(plot_dt_area, origin_dt),
  use.names = TRUE,
  fill = TRUE
)

setorder(plot_dt_area_path, method, total_area_km2, stage_order)
setorder(plot_dt_area, method, total_area_km2, stage_order)

x_max <- max(plot_dt_area$total_area_km2, na.rm = TRUE)

p_area_persist <- ggplot() +
  geom_path(
    data = plot_dt_area_path,
    aes(x = total_area_km2, y = sp_persist, color = method_lab),
    linewidth = 0.95,
    lineend = "round"
  ) +
  geom_point(
    data = plot_dt_area,
    aes(x = total_area_km2, y = sp_persist, color = method_lab),
    size = 1.8
  ) +
  scale_color_manual(values = method_cols, name = NULL) +
  scale_x_continuous(
    limits = c(0, x_max * 1.02),
    breaks = pretty(c(0, x_max), n = 4),
    labels = scales::label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = persist_breaks_shared,
    labels = persist_labels_shared,
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = expression("Retained PU area (km"^2*")"),
    y = "Species persistence"
  ) +
  theme_persist(base_size = 11) +
  theme(
    legend.position = "none"
  )

# ============================================================
# 09_plot_combined_figure - one 4-panel figure
# Legend layout:
#   left  = stacked colour-bar legend
#   right = method legend
# ============================================================

stopifnot(
  exists("p_mean"),
  exists("p_mean_leg"),
  exists("p_focus"),
  exists("p_area"),
  exists("p_area_persist"),
  exists("p_area_leg")
)

legend_methods <- cowplot::get_legend(p_mean_leg)

p_a <- p_mean +
  theme(
    plot.margin = margin(10, 8, 8, 14),
    axis.title.x = element_text(margin = margin(t = 6))
  )

p_b <- p_focus +
  theme(
    plot.margin = margin(10, 5, 8, 14),
    axis.title.x = element_text(margin = margin(t = 6))
  )

p_c <- p_area +
  theme(
    plot.margin = margin(10, 8, 5, 14)
  )

p_d <- p_area_persist +
  theme(
    plot.margin = margin(10, 5, 5, 14)
  )

# -------------------------
# small header strips above each panel
# keeps the plot panel itself the same size
# -------------------------

n_included_species <- uniqueN(par$species)

hdr_all   <- sprintf("All species (n = %d)", n_included_species)
hdr_focal <- paste0(focal_name)

panel_header <- function(txt, size = 10) {
  cowplot::ggdraw() +
    cowplot::draw_label(
      txt,
      x = 0.6, y = 0.2,
      hjust = 0.5, vjust = 0.5,
      size = size,
      color = "grey20",
      fontface = "plain"
    )
}

add_header_strip <- function(panel, txt, rel_header = 0.09) {
  cowplot::plot_grid(
    panel_header(txt),
    panel,
    ncol = 1,
    rel_heights = c(rel_header, 1)
  )
}

tag_panel <- function(p, tag, draw_x = 0.055, draw_y = 0.02, draw_w = 0.945, draw_h = 0.945) {
  cowplot::ggdraw() +
    cowplot::draw_plot(p, x = draw_x, y = draw_y, width = draw_w, height = draw_h) +
    cowplot::draw_plot_label(
      tag,
      x = 0.01, y = 0.99,
      hjust = 0, vjust = 1,
      fontface = "bold",
      size = 12.5
    )
}

p_a_tag <- tag_panel(p_a, "(a)")
p_b_tag <- tag_panel(p_b, "(b)")
p_c_tag <- tag_panel(p_c, "(c)")
p_d_tag <- tag_panel(p_d, "(d)")

p_a_box <- add_header_strip(p_a_tag, hdr_all)
p_b_box <- add_header_strip(p_b_tag, hdr_focal)
p_c_box <- add_header_strip(p_c_tag, hdr_focal)
p_d_box <- add_header_strip(p_d_tag, hdr_focal)

panel_grid <- cowplot::plot_grid(
  p_a_box, p_b_box, p_c_box, p_d_box,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

# -------------------------
# knobs for legend placement
# -------------------------

# move the ENTIRE legend block right/left
legend_shift_right <- 0.13   # increase -> moves both legends right
legend_total_width <- 0.86   # total width used by both legends together

# relative widths of the left and right legend slots
left_leg_rel_width  <- 0.58  # stacked colour-bar legend slot
right_leg_rel_width <- 0.42  # method legend slot

# LEFT legend (stacked colour-bar legend) placement inside its slot
pu_leg_inner_x      <- 0.11  # increase -> moves left legend right within left slot
pu_leg_inner_width  <- 0.52  # increase -> makes left legend wider

# RIGHT legend (method legend) placement inside its slot
method_leg_inner_x     <- 0.07  # increase -> moves right legend right within right slot
method_leg_inner_width <- 0.82  # increase -> makes right legend wider

# -------------------------
# build legend row
# -------------------------

p_area_leg_left <- cowplot::ggdraw() +
  cowplot::draw_plot(
    p_area_leg,
    x = pu_leg_inner_x,
    y = 0.02,
    width = pu_leg_inner_width,
    height = 0.96
  )

legend_methods_right <- cowplot::ggdraw() +
  cowplot::draw_plot(
    legend_methods,
    x = method_leg_inner_x,
    y = 0.00,
    width = method_leg_inner_width,
    height = 1.00
  )

legend_row_raw <- cowplot::plot_grid(
  p_area_leg_left,
  legend_methods_right,
  nrow = 1,
  rel_widths = c(left_leg_rel_width, right_leg_rel_width)
)

legend_row <- cowplot::ggdraw() +
  cowplot::draw_plot(
    legend_row_raw,
    x = legend_shift_right,
    y = 0,
    width = legend_total_width,
    height = 1
  )

fig_cmp_4panel <- cowplot::plot_grid(
  panel_grid,
  legend_row,
  ncol = 1,
  rel_heights = c(1, 0.14)
)

fig_cmp_4panel

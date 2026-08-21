# Shiny app: Trial 1a x Trial 1b spatial overlap
# Tab 1: SSA overlay map by species (dropdown)
# Tab 2: An. arabiensis zoomed to Mozambique with locator inset

library(shiny)
library(terra)
library(sf)

data_dir <- "data"

paths <- list(
  ara  = file.path(data_dir, "An_arabiensis.tif"),
  col  = file.path(data_dir, "An_coluzzii.tif"),
  mou  = file.path(data_dir, "An_moucheti.tif"),
  gam  = file.path(data_dir, "An_gambiae.tif"),
  ste  = file.path(data_dir, "An_stephensi.tif"),
  pfpr = file.path(data_dir, "PfPR_mean.tif"),
  itn  = file.path(data_dir, "ITN_use_rate.tif"),
  pop  = file.path(data_dir, "POP_MEAN_2000_2020_5km.tif"),
  gadm = file.path(data_dir, "gadm_ssa.gpkg")
)

to01 <- function(r) {
  mx <- suppressWarnings(as.numeric(global(r, "max", na.rm = TRUE)))
  r1 <- if (is.finite(mx) && mx > 1) r / 100 else r
  clamp(r1, 0, 1, values = TRUE)
}
safe_div    <- function(num, den, eps = 1e-9) num / (den + eps)
fast_align  <- function(r, template, method = "bilinear") {
  if (!identical(crs(r), crs(template))) project(r, template, method = method)
  else resample(r, template, method = method)
}

cat("Loading data...\n")
ssa_iso <- c(
  "AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM","COG","CIV","COD","GNQ","ERI",
  "SWZ","ETH","GAB","GMB","GHA","GIN","GNB","KEN","LSO","LBR","MDG","MWI","MLI","MRT","MUS","MOZ",
  "NAM","NER","NGA","RWA","STP","SEN","SYC","SLE","ZAF","SSD","TGO","UGA","TZA","ZMB","ZWE",
  "SDN","SOM","DJI"
)
africa_iso <- unique(c(ssa_iso, "MAR","TUN","LBY","EGY","DZA"))

a0      <- st_read(paths$gadm, layer = "ADM_0", quiet = TRUE)
a1      <- st_read(paths$gadm, layer = "ADM_1", quiet = TRUE)
ADMIN0  <- vect(a0[a0$GID_0 %in% ssa_iso, ])
ADMIN1  <- vect(a1[a1$GID_0 %in% ssa_iso, ])
AFRICA0 <- vect(a0[a0$GID_0 %in% africa_iso, ])
MOZ0    <- vect(a0[a0$GID_0 == "MOZ", ])
MOZ1    <- vect(a1[a1$GID_0 == "MOZ", ])

NEIGH_ISO <- c("ZWE","ZMB","MWI","TZA","ZAF","SWZ")
NEIGH     <- vect(a0[a0$GID_0 %in% NEIGH_ISO, ])

PFPR_MEAN <- rast(paths$pfpr)
ADMIN0_P  <- if (!identical(crs(ADMIN0), crs(PFPR_MEAN))) project(ADMIN0, crs(PFPR_MEAN)) else ADMIN0
TEMPLATE  <- crop(PFPR_MEAN, ADMIN0_P, snap = "out")
ADMIN0_A  <- if (!identical(crs(ADMIN0),  crs(TEMPLATE))) project(ADMIN0,  crs(TEMPLATE)) else ADMIN0
ADMIN1_A  <- if (!identical(crs(ADMIN1),  crs(TEMPLATE))) project(ADMIN1,  crs(TEMPLATE)) else ADMIN1
AFRICA0_A <- if (!identical(crs(AFRICA0), crs(TEMPLATE))) project(AFRICA0, crs(TEMPLATE)) else AFRICA0
MOZ0_A    <- if (!identical(crs(MOZ0),    crs(TEMPLATE))) project(MOZ0,    crs(TEMPLATE)) else MOZ0
MOZ1_A    <- if (!identical(crs(MOZ1),    crs(TEMPLATE))) project(MOZ1,    crs(TEMPLATE)) else MOZ1
NEIGH_A   <- if (!identical(crs(NEIGH),   crs(TEMPLATE))) project(NEIGH,   crs(TEMPLATE)) else NEIGH

PFPR_ALIGNED <- to01(fast_align(PFPR_MEAN, TEMPLATE))
ITN_ALIGNED  <- to01(fast_align(rast(paths$itn), TEMPLATE))

species_rasters <- list(
  arabiensis = rast(paths$ara), coluzzii = rast(paths$col),
  moucheti   = rast(paths$mou), gambiae  = rast(paths$gam)
)
SPECIES_ALIGNED <- lapply(species_rasters, function(r) fast_align(to01(r), TEMPLATE))

pres_floor <- 0.10

# UC1a: low PfPR = Score 1 | UC1b: high PfPR = Score 1
uc1a <- list(dom_t1=0.70, dom_t2=0.40, itn_t1=0.60, itn_t2=0.80, pfpr_t1=0.15, pfpr_t2=0.40)
uc1b <- list(dom_t1=0.70, dom_t2=0.40, itn_t1=0.60, itn_t2=0.80, pfpr_t1=0.40, pfpr_t2=0.15)

compute_tier <- function(dom_v, valid_v, itn_v, pfpr_v, p, pfpr_high = FALSE) {
  dt2 <- if (p$dom_t1 <= p$dom_t2) 0.25 else p$dom_t2
  it2 <- if (p$itn_t1 >= p$itn_t2) 0.90 else p$itn_t2
  tv <- rep(NA_real_, length(dom_v))
  tv[!is.na(valid_v) & !is.na(dom_v) & dom_v >= p$dom_t1] <- 1
  tv[!is.na(valid_v) & !is.na(dom_v) & dom_v >= dt2 & dom_v < p$dom_t1] <- 2
  tv[!is.na(valid_v) & !is.na(dom_v) & dom_v <  dt2] <- 3
  ti <- rep(NA_real_, length(itn_v))
  ti[!is.na(valid_v) & !is.na(itn_v) & itn_v <  p$itn_t1] <- 1
  ti[!is.na(valid_v) & !is.na(itn_v) & itn_v >= p$itn_t1 & itn_v < it2] <- 2
  ti[!is.na(valid_v) & !is.na(itn_v) & itn_v >= it2] <- 3
  tp <- rep(NA_real_, length(pfpr_v))
  if (!pfpr_high) {
    tp[!is.na(valid_v) & !is.na(pfpr_v) & pfpr_v <  p$pfpr_t1] <- 1
    tp[!is.na(valid_v) & !is.na(pfpr_v) & pfpr_v >= p$pfpr_t1 & pfpr_v < p$pfpr_t2] <- 2
    tp[!is.na(valid_v) & !is.na(pfpr_v) & pfpr_v >= p$pfpr_t2] <- 3
  } else {
    tp[!is.na(valid_v) & !is.na(pfpr_v) & pfpr_v >= p$pfpr_t1] <- 1
    tp[!is.na(valid_v) & !is.na(pfpr_v) & pfpr_v >= p$pfpr_t2 & pfpr_v < p$pfpr_t1] <- 2
    tp[!is.na(valid_v) & !is.na(pfpr_v) & pfpr_v <  p$pfpr_t2] <- 3
  }
  sc <- tv + ti + tp
  ov <- rep(NA_real_, length(sc))
  ov[!is.na(sc) & sc >= 3 & sc <= 4] <- 1
  ov[!is.na(sc) & sc >= 5 & sc <= 7] <- 2
  ov[!is.na(sc) & sc >= 8 & sc <= 9] <- 3
  ov
}

zone_cols <- c(
  "1" = "#1B9E77",
  "2" = "#CCEBC5",
  "3" = "#FDDDAC",
  "4" = "#D95F02",
  "5" = "#A6B1C2"
)
zone_labels <- c(
  "1" = "Tier 1 under both trials (dual-purpose)",
  "2" = "Tier 1 under entomological trial only",
  "3" = "Tier 1 under epidemiological trial only",
  "4" = "Tier 2 under both trials",
  "5" = "Excluded"
)

# Pre-compute zone rasters for all species
cat("Pre-computing zone rasters...\n")
ZONE_RASTERS <- lapply(names(SPECIES_ALIGNED), function(sp) {
  base_p  <- SPECIES_ALIGNED[[sp]]
  base_p0 <- base_p; base_p0[is.na(base_p0)] <- 0
  others  <- rast(lapply(setdiff(names(SPECIES_ALIGNED), sp),
                         function(nm) { r <- SPECIES_ALIGNED[[nm]]; r[is.na(r)] <- 0; r }))
  dom     <- safe_div(base_p0, base_p0 + app(others, sum))
  valid_m <- ifel(!is.na(PFPR_ALIGNED) & !is.na(ITN_ALIGNED) &
                    !is.na(base_p) & base_p >= pres_floor, 1, NA)
  valid_v <- values(valid_m)
  itn_v   <- values(ITN_ALIGNED)
  pfpr_v  <- values(PFPR_ALIGNED)
  dom_v   <- values(dom)
  
  t1a_v <- compute_tier(dom_v, valid_v, itn_v, pfpr_v, uc1a, pfpr_high = FALSE)
  t1b_v <- compute_tier(dom_v, valid_v, itn_v, pfpr_v, uc1b, pfpr_high = TRUE)
  
  zone_v <- rep(NA_real_, length(t1a_v))
  zone_v[!is.na(t1a_v) & !is.na(t1b_v) & t1a_v == 1 & t1b_v == 1] <- 1
  zone_v[!is.na(t1a_v) & !is.na(t1b_v) & t1a_v == 1 & t1b_v != 1] <- 2
  zone_v[!is.na(t1a_v) & !is.na(t1b_v) & t1a_v != 1 & t1b_v == 1] <- 3
  zone_v[!is.na(t1a_v) & !is.na(t1b_v) & t1a_v == 2 & t1b_v == 2] <- 4
  zone_v[!is.na(t1a_v) & !is.na(t1b_v) & (t1a_v == 3 | t1b_v == 3)] <- 5
  
  zr <- rast(TEMPLATE); values(zr) <- zone_v
  mask(zr, ADMIN0_A)
})
names(ZONE_RASTERS) <- names(SPECIES_ALIGNED)

# Africa and SSA extents
af      <- ext(AFRICA0_A)
xlim_all <- c(xmin(af), xmax(af))
ylim_all <- c(ymin(af), ymax(af))

# Mozambique extents
moz_xlim  <- c(xmin(MOZ0_A) - 0.5, xmax(MOZ0_A) + 0.5)
moz_ylim  <- c(ymin(MOZ0_A) - 0.5, ymax(MOZ0_A) + 0.5)
e_moz     <- ext(moz_xlim[1], moz_xlim[2], moz_ylim[1], moz_ylim[2])
zr_moz    <- crop(ZONE_RASTERS[["arabiensis"]], e_moz)
v1_moz    <- crop(ADMIN1_A, e_moz)
NEIGH_MOZ <- crop(NEIGH_A, ext(moz_xlim[1]-1, moz_xlim[2]+1, moz_ylim[1]-1, moz_ylim[2]+1))
vals_moz  <- sort(unique(na.omit(values(zr_moz))))
vals_moz  <- vals_moz[vals_moz %in% 1:5]

cat("Ready.\n")

# ── UI ────────────────────────────────────────────────────────
ui <- fluidPage(
  
  tags$head(tags$style(HTML("
    body { font-family: 'Georgia', serif; background: #fafaf8; }
    .well { background: #f0ede8; border: none; border-radius: 8px; }
    h3 { color: #2c5f2e; font-weight: bold; border-bottom: 2px solid #2c5f2e; padding-bottom: 6px; }
    h4 { color: #4a4a4a; margin-top: 14px; }
    .zone-legend { margin-top: 12px; }
    .zone-chip { display: inline-block; width: 14px; height: 14px;
                 border-radius: 2px; margin-right: 6px; vertical-align: middle; }
  "))),
  
  titlePanel(
    div(
      h3("Demonstration trials — spatial overlap"),
      p(style = "color:#666; font-size:0.9em; margin-top:-8px;",
        "Trial 1a \u00d7 Trial 1b eligibility overlap across ", em("Anopheles"), " species")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 
                 h4("Species"),
                 selectInput("species", NULL,
                             choices = c(
                               "An. arabiensis"   = "arabiensis",
                               "An. coluzzii"     = "coluzzii",
                               "An. moucheti"     = "moucheti",
                               "An. gambiae s.s." = "gambiae"
                             ),
                             selected = "arabiensis"
                 ),
                 
                 hr(),
                 div(class = "zone-legend",
                     h4("Legend"),
                     tags$div(
                       tags$span(class="zone-chip", style="background:#1B9E77;"),
                       "Tier 1 under both trials (dual-purpose)", tags$br(),
                       tags$span(class="zone-chip", style="background:#CCEBC5;"),
                       "Tier 1 under entomological trial only", tags$br(),
                       tags$span(class="zone-chip", style="background:#FDDDAC;"),
                       "Tier 1 under epidemiological trial only", tags$br(),
                       tags$span(class="zone-chip", style="background:#D95F02;"),
                       "Tier 2 under both trials", tags$br(),
                       tags$span(class="zone-chip", style="background:#A6B1C2;"),
                       "Excluded"
                     )
                 ),
                 
                 hr(),
                 helpText("Note: Tab 2 (Mozambique zoom) is fixed to An. arabiensis — the only species with meaningful dual-purpose sites.")
    ),
    
    mainPanel(width = 9,
              tabsetPanel(
                
                tabPanel("SSA Overview",
                         br(),
                         plotOutput("map_ssa", height = 620)
                ),
                
                tabPanel("Mozambique zoom \u2014 An. arabiensis",
                         br(),
                         plotOutput("map_moz", height = 620)
                )
              )
    )
  )
)

# ── Server ────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # Helper: get species-level extent
  get_extent <- function(zone_r) {
    all_idx <- which(!is.na(values(zone_r)))
    if (length(all_idx) > 0) {
      xy  <- xyFromCell(zone_r, unique(all_idx))
      xr  <- range(xy[,1]); yr <- range(xy[,2])
      pad <- max(diff(xr), diff(yr)) * 0.05
      xlim <- c(max(xr[1]-pad, xlim_all[1]), min(xr[2]+pad, xlim_all[2]))
      ylim <- c(max(yr[1]-pad, ylim_all[1]), min(yr[2]+pad, ylim_all[2]))
    } else {
      xlim <- xlim_all; ylim <- ylim_all
    }
    list(xlim=xlim, ylim=ylim)
  }
  
  # Plot helper
  plot_zone <- function(zone_r, xlim, ylim, show_moz_box = FALSE) {
    e_view  <- ext(xlim[1], xlim[2], ylim[1], ylim[2])
    af0_v   <- crop(AFRICA0_A, e_view)
    v0_v    <- crop(ADMIN0_A,  e_view)
    v1_v    <- crop(ADMIN1_A,  e_view)
    zr_v    <- crop(zone_r,    e_view)
    vals_v  <- sort(unique(na.omit(values(zr_v))))
    vals_v  <- vals_v[vals_v %in% 1:5]
    
    par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "white")
    plot(af0_v, col = "#e8e8e8", border = "#cccccc", lwd = 0.3,
         axes = FALSE, main = "", xlim = xlim, ylim = ylim)
    plot(v0_v, col = "#f5f5f5", border = NA, add = TRUE)
    
    if (length(vals_v) > 0) {
      plot(zr_v, add = TRUE,
           col    = zone_cols[as.character(vals_v)],
           breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5),
           legend = FALSE, axes = FALSE)
    }
    plot(v1_v,  add = TRUE, lwd = 0.3, col = NA, border = "#cccccc")
    plot(v0_v,  add = TRUE, lwd = 0.5, col = NA, border = "#888888")
    plot(af0_v, add = TRUE, lwd = 0.4, col = NA, border = "#aaaaaa")
    
    if (show_moz_box) {
      moz_ext <- ext(MOZ0_A)
      rect(xmin(moz_ext), ymin(moz_ext), xmax(moz_ext), ymax(moz_ext),
           border = "black", lwd = 1.8)
    }
  }
  
  # ── Tab 1: SSA overview ──────────────────────────────────────
  output$map_ssa <- renderPlot({
    sp     <- input$species
    zone_r <- ZONE_RASTERS[[sp]]
    ex     <- get_extent(zone_r)
    
    plot_zone(zone_r, ex$xlim, ex$ylim,
              show_moz_box = (sp == "arabiensis"))
  })
  
  # ── Tab 2: Mozambique zoom ────────────────────────────────────
  output$map_moz <- renderPlot({
    par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "white")
    
    # Neighbours first
    plot(NEIGH_MOZ, col = "#e8e8e8", border = "#cccccc", lwd = 0.3,
         axes = FALSE, main = "",
         xlim = moz_xlim, ylim = moz_ylim)
    
    # Mozambique background
    plot(MOZ0_A, col = "#f5f5f5", border = NA, add = TRUE)
    
    # Zone raster
    if (length(vals_moz) > 0) {
      plot(zr_moz, add = TRUE,
           col    = zone_cols[as.character(vals_moz)],
           breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5),
           legend = FALSE, axes = FALSE)
    }
    
    # Province borders and country border
    plot(v1_moz, add = TRUE, lwd = 0.7, col = NA, border = "#888888")
    plot(MOZ0_A, add = TRUE, lwd = 2.0, col = NA, border = "black")
    
    # Province labels
    if ("NAME_1" %in% names(MOZ1_A)) {
      cents  <- centroids(MOZ1_A)
      coords <- crds(cents)
      text(coords[,1], coords[,2],
           labels = MOZ1_A$NAME_1,
           cex = 0.75, col = "grey20", font = 1)
    }
    
    # Locator inset — Africa silhouette with Mozambique in red
    fig_b  <- par("fig")
    ix1 <- fig_b[1] + 0.60 * (fig_b[2] - fig_b[1])
    ix2 <- fig_b[2] - 0.01
    iy1 <- fig_b[4] - 0.38 * (fig_b[4] - fig_b[3])
    iy2 <- fig_b[4] - 0.01
    
    par(fig = c(ix1, ix2, iy1, iy2), mar = c(0,0,0,0), new = TRUE, bg = "white")
    plot(AFRICA0_A,
         col    = "#d0d0d0",
         border = "#aaaaaa",
         lwd    = 0.3,
         axes   = FALSE,
         xlim   = c(xmin(af), xmax(af)),
         ylim   = c(ymin(af), ymax(af)))
    plot(MOZ0_A, add = TRUE, col = "#d7191c", border = "#d7191c", lwd = 0.5)
    box(col = "grey50", lwd = 0.8)
  })
}

shinyApp(ui, server)
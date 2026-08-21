# Shiny app for Use Case 4: Prevention of An. stephensi re-emergence
# Urban filter (GHS >= 21) applied first; stephensi occurrence and PfPR scored 1-3

library(shiny)
library(terra)
library(sf)

data_dir <- "data"

paths <- list(
  ste  = file.path(data_dir, "An_stephensi.tif"),
  pfpr = file.path(data_dir, "PfPR_mean.tif"),
  itn  = file.path(data_dir, "ITN_use_rate.tif"),
  inc  = file.path(data_dir, "Pf_Incidence_mean_2000.tif"),
  pop  = file.path(data_dir, "POP_MEAN_2000_2020_5km.tif"),
  ghs  = file.path("Data", "GHS_SMOD_E2025_GLOBE_R2023A_54009_1000_V2_0.tif"),
  gadm = file.path(data_dir, "gadm_ssa.gpkg")
)

to01 <- function(r) {
  mx <- suppressWarnings(as.numeric(global(r, "max", na.rm = TRUE)))
  r1 <- if (is.finite(mx) && mx > 1) r / 100 else r
  clamp(r1, 0, 1, values = TRUE)
}

fast_align <- function(r, template, method = "bilinear") {
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
africa_iso <- unique(c(ssa_iso, "MAR","TUN","LBY","EGY","DZA","SDN","SOM","DJI"))

a0       <- st_read(paths$gadm, layer = "ADM_0", quiet = TRUE)
a1       <- st_read(paths$gadm, layer = "ADM_1", quiet = TRUE)
ADMIN0   <- vect(a0[a0$GID_0 %in% ssa_iso,  ])

# Named vector: country name -> ISO code for zoom dropdown
afro_sf       <- a0[a0$GID_0 %in% ssa_iso, ]
name_col      <- if ("COUNTRY" %in% names(afro_sf)) "COUNTRY" else "NAME_0"
COUNTRY_CHOICES <- c("None (full Africa)" = "",
                     setNames(sort(afro_sf$GID_0),
                              afro_sf[[name_col]][order(afro_sf$GID_0)]))
ADMIN1   <- vect(a1[a1$GID_0 %in% ssa_iso,  ])
AFRICA0  <- vect(a0[a0$GID_0 %in% africa_iso, ])
AFRICA1  <- vect(a1[a1$GID_0 %in% africa_iso, ])

PFPR_MEAN <- rast(paths$pfpr)
ADMIN0_P  <- if (!identical(crs(ADMIN0), crs(PFPR_MEAN))) project(ADMIN0, crs(PFPR_MEAN)) else ADMIN0
TEMPLATE  <- crop(PFPR_MEAN, ADMIN0_P, snap = "out")
ADMIN0_A  <- if (!identical(crs(ADMIN0), crs(TEMPLATE))) project(ADMIN0, crs(TEMPLATE)) else ADMIN0
ADMIN1_A  <- if (!identical(crs(ADMIN1), crs(TEMPLATE))) project(ADMIN1, crs(TEMPLATE)) else ADMIN1
AFRICA0_A <- if (!identical(crs(AFRICA0), crs(TEMPLATE))) project(AFRICA0, crs(TEMPLATE)) else AFRICA0
AFRICA1_A <- if (!identical(crs(AFRICA1), crs(TEMPLATE))) project(AFRICA1, crs(TEMPLATE)) else AFRICA1

PFPR_ALIGNED <- to01(fast_align(PFPR_MEAN, TEMPLATE))
ITN_ALIGNED  <- to01(fast_align(rast(paths$itn), TEMPLATE))
INC_ALIGNED  <- fast_align(rast(paths$inc), TEMPLATE)
POP_ALIGNED  <- fast_align(rast(paths$pop), TEMPLATE)
BURDEN_GLOBAL <- POP_ALIGNED * INC_ALIGNED
CELL_AREA_KM2 <- cellSize(TEMPLATE, unit = "km")

# Stephensi
STE_RAW     <- to01(fast_align(rast(paths$ste), TEMPLATE))
STE_ALIGNED <- mask(STE_RAW, ADMIN0_A)

# GHS-SMOD — aligned to template
GHS_RAW     <- rast(paths$ghs)
GHS_ALIGNED <- fast_align(GHS_RAW, TEMPLATE, method = "mode")
GHS_ALIGNED <- mask(GHS_ALIGNED, ADMIN0_A)

# Stephensi urban presence denominator
# Total An. stephensi predicted area (presence >= ste_floor, regardless of GHS)
# Used as denominator for % of stephensi range that is eligible
STE_TOTAL_MASK <- ifel(!is.na(STE_ALIGNED) & STE_ALIGNED >= 0.10, 1, NA)
STE_TOTAL_MASK <- mask(STE_TOTAL_MASK, ADMIN0_A)
STE_URBAN_AREA <- as.numeric(global(mask(CELL_AREA_KM2, STE_TOTAL_MASK),
                                    "sum", na.rm = TRUE)[1,1])

cat("Data loaded.\n")

ste_floor  <- 0.10   # minimum stephensi presence
ghs_min    <- 21     # urban/peri-urban minimum (suburban)

# Stephensi occurrence tiers (ceiling) — fixed
ste_t1 <- 0.70;  ste_t2 <- 0.40

ui <- fluidPage(
  
  tags$head(tags$style(HTML("
    body { font-family: 'Georgia', serif; background: #fafaf8; color: #2c2c2c; }
    .well { background: #f0ede8; border: none; border-radius: 8px; }
    h3 { color: #6e2f8a; font-weight: bold; border-bottom: 2px solid #6e2f8a; padding-bottom: 6px; }
    h4 { color: #4a4a4a; margin-top: 18px; }
    .method-box { background: #f5eef8; border-left: 4px solid #6e2f8a;
                  padding: 10px 14px; border-radius: 4px; font-size: 0.88em;
                  margin-bottom: 12px; }
    .filter-box { background: #fdfefe; border-left: 4px solid #117a65;
                  padding: 10px 14px; border-radius: 4px; font-size: 0.88em;
                  margin-bottom: 8px; }
    .tier-legend { display: flex; gap: 12px; margin: 8px 0; }
    .tier-chip { padding: 4px 10px; border-radius: 12px; font-size: 0.82em;
                 font-weight: bold; color: white; }
  "))),
  
  titlePanel(
    div(
      h3("Use case 4 — Prevention of re-emergence"),
      p(style = "color:#666; font-size:0.9em; margin-top:-8px;", "Spatial eligibility for preventive Wolbachia deployment targeting ", em("An. stephensi"))
    )
  ),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 
                 div(class = "method-box",
                     strong("How this works:"), br(),
                     "Urban/peri-urban pixels (GHS-SMOD \u2265 21) are filtered first. ",
                     "An. stephensi occurrence probability and PfPR are each scored 1\u20133. ",
                     "Tier 1 (score 2\u20133) = highest priority, Tier 2 (score 4\u20135), Tier 3 (score 6)."
                 ),
                 
                 hr(),
                 h4("An. stephensi occurrence"),
                 helpText("Higher occurrence probability = better. Score 1 = high predicted presence."),
                 sliderInput("ste_t1", "Score 1 minimum occurrence (%)",
                             min = 30, max = 90, value = 70, step = 5),
                 sliderInput("ste_t2", "Score 2 minimum occurrence (%)",
                             min = 10, max = 80, value = 40, step = 5),
                 
                 hr(),
                 h4("Malaria transmission (PfPR)"),
                 helpText("Moderate PfPR = best (15\u201340%). Score 1 = PfPR in this sweet spot."),
                 sliderInput("pfpr_lo", "PfPR lower bound for Score 1 (%)",
                             min = 5, max = 30, value = 15, step = 1),
                 sliderInput("pfpr_hi", "PfPR upper bound for Score 1 (%)",
                             min = 20, max = 60, value = 40, step = 1),
                 
                 hr(),
                 h4("Display"),
                 checkboxInput("auto_zoom", "Auto-zoom to eligible area", TRUE),
                 selectInput("zoom_country", "Zoom to country",
                             choices  = COUNTRY_CHOICES,
                             selected = "",
                             width    = "100%"),
                 
                 div(class = "tier-legend",
                     div(class = "tier-chip", style = "background:#1a9641;", "Tier 1"),
                     div(class = "tier-chip", style = "background:#f4d03f; color:#333;", "Tier 2"),
                     div(class = "tier-chip", style = "background:#d7191c;", "Tier 3")
                 )
    ),
    
    mainPanel(width = 9,
              tabsetPanel(
                tabPanel("Map",
                         plotOutput("map_uc4", height = 650)
                ),
                tabPanel("Impact & Coverage",
                         br(),
                         h4("Coverage"),
                         tableOutput("cov_table"),
                         hr(),
                         h4("Population at risk"),
                         tableOutput("pop_table")
                )
              )
    )
  )
)

server <- function(input, output, session) {
  
  tier_uc4 <- reactive({
    ste_k1   <- input$ste_t1 / 100
    ste_k2   <- input$ste_t2 / 100
    pfpr_lo  <- input$pfpr_lo / 100
    pfpr_hi  <- input$pfpr_hi / 100
    
    urban  <- !is.na(GHS_ALIGNED) & GHS_ALIGNED >= ghs_min
    ste_ok <- !is.na(STE_ALIGNED) & STE_ALIGNED >= ste_floor
    
    ts <- rast(TEMPLATE); values(ts) <- NA
    ts[urban & ste_ok & STE_ALIGNED >= ste_k1] <- 1
    ts[urban & ste_ok & STE_ALIGNED >= ste_k2 & STE_ALIGNED < ste_k1] <- 2
    ts[urban & ste_ok & STE_ALIGNED <  ste_k2] <- 3
    
    tp <- rast(TEMPLATE); values(tp) <- NA
    tp[urban & ste_ok & PFPR_ALIGNED >= pfpr_lo & PFPR_ALIGNED <= pfpr_hi] <- 1
    tp[urban & ste_ok & ((PFPR_ALIGNED >= pfpr_hi * 0.5 & PFPR_ALIGNED < pfpr_lo) |
                           (PFPR_ALIGNED > pfpr_hi & PFPR_ALIGNED <= pfpr_hi * 1.5))] <- 2
    tp[urban & ste_ok & !is.na(PFPR_ALIGNED)] <- 3
    tp[urban & ste_ok & PFPR_ALIGNED >= pfpr_lo & PFPR_ALIGNED <= pfpr_hi] <- 1
    tp[urban & ste_ok & PFPR_ALIGNED >= (pfpr_lo * 0.5) & PFPR_ALIGNED < pfpr_lo] <- 2
    tp[urban & ste_ok & PFPR_ALIGNED > pfpr_hi & PFPR_ALIGNED <= (pfpr_hi * 1.5)] <- 2
    tp[urban & ste_ok & (PFPR_ALIGNED < (pfpr_lo * 0.5) | PFPR_ALIGNED > (pfpr_hi * 1.5))] <- 3
    
    sc  <- ts + tp
    out <- rast(TEMPLATE); values(out) <- NA
    out[!is.na(sc) & sc <= 3] <- 1
    out[!is.na(sc) & sc >= 4 & sc <= 5] <- 2
    out[!is.na(sc) & sc == 6] <- 3
    out <- mask(out, ADMIN0_A)
    out
  })
  
  plot_tier <- function(tier_r, title_txt) {
    vals <- sort(unique(na.omit(as.vector(values(tier_r)))))
    vals <- vals[vals %in% 1:3]
    if (length(vals) == 0) { plot.new(); title("No eligible cells"); return(invisible()) }
    
    tier_cols   <- c("1" = "#1a9641", "2" = "#f4d03f", "3" = "#d7191c")
    tier_breaks <- c(0.5, 1.5, 2.5, 3.5)
    
    af   <- ext(AFRICA0_A)
    xlim <- c(xmin(af), xmax(af)); ylim <- c(ymin(af), ymax(af))
    
    zoom_iso <- input$zoom_country
    if (nchar(zoom_iso) > 0) {
      ctry_v <- vect(afro_sf[afro_sf$GID_0 == zoom_iso, ])
      ctry_v <- if (!identical(crs(ctry_v), crs(TEMPLATE))) project(ctry_v, crs(TEMPLATE)) else ctry_v
      ce    <- ext(ctry_v)
      dx    <- diff(c(xmin(ce), xmax(ce))); dy <- diff(c(ymin(ce), ymax(ce)))
      pad   <- max(dx, dy) * 0.25
      xlim  <- c(xmin(ce) - pad, xmax(ce) + pad)
      ylim  <- c(ymin(ce) - pad, ymax(ce) + pad)
    } else if (isTRUE(input$auto_zoom)) {
      idx <- which(!is.na(values(tier_r)))
      if (length(idx) > 0) {
        xy <- xyFromCell(tier_r, idx)
        xr <- range(xy[,1]); yr <- range(xy[,2])
        dx <- diff(xr); dy <- diff(yr)
        if (!is.finite(dx)||dx==0) dx <- res(tier_r)[1]*2
        if (!is.finite(dy)||dy==0) dy <- res(tier_r)[2]*2
        xlim <- c(max(xr[1]-dx*0.1, xmin(af)), min(xr[2]+dx*0.1, xmax(af)))
        ylim <- c(max(yr[1]-dy*0.1, ymin(af)), min(yr[2]+dy*0.1, ymax(af)))
      }
    }
    e_view    <- ext(xlim[1],xlim[2],ylim[1],ylim[2])
    af0_view  <- crop(AFRICA0_A, e_view)
    v0_view   <- crop(ADMIN0_A, e_view)
    v1_view   <- crop(ADMIN1_A, e_view)
    tier_view <- crop(tier_r, e_view)
    
    par(mar = c(1,1,2.5,1))
    plot(af0_view, col="#e4e4e4", border=NA, axes=FALSE, xlim=xlim, ylim=ylim,
         main=title_txt, cex.main=1.0)
    plot(v0_view, add=TRUE, col="#f0f0f0", border=NA)
    
    if (length(vals) > 0) {
      plot(tier_view, add=TRUE,
           col    = c("#1a9641", "#f4d03f", "#d7191c"),
           breaks = c(0.5, 1.5, 2.5, 3.5),
           legend = FALSE, axes = FALSE)
    }
    
    plot(v1_view,  add=TRUE, lwd=0.3, col=NA, border="#cccccc")
    plot(v0_view,  add=TRUE, lwd=0.6, col=NA, border="#888888")
    plot(af0_view, add=TRUE, lwd=0.4, col=NA, border="#aaaaaa")
    
    legend("bottomleft",
           legend = c("Tier 1", "Tier 2", "Tier 3"),
           fill   = c("#1a9641", "#f4d03f", "#d7191c"),
           border = "grey40",
           bty    = "o", bg = "white", box.col = "grey80",
           cex    = 0.85, xpd = TRUE, title.font = 2)
  }
  
  output$map_uc4 <- renderPlot({ plot_tier(tier_uc4(), "") })
  
  output$cov_table <- renderTable({
    tr  <- tier_uc4()
    m1  <- ifel(tr == 1, 1, NA)
    m12 <- ifel(tr <= 2 & !is.na(tr), 1, NA)
    area <- function(m) {
      a <- as.numeric(global(mask(CELL_AREA_KM2, m), "sum", na.rm=TRUE)[1,1])
      if (!is.finite(a)) 0 else a
    }
    pct_ste <- function(m) {
      a <- area(m)
      paste0(formatC(100 * a / STE_URBAN_AREA, format="f", digits=2), "% of stephensi range")
    }
    data.frame(
      Tier = c("Tier 1", "Tier 1 + 2 (cumulative)"),
      `Area (km²)` = format(round(c(area(m1), area(m12))), big.mark=","),
      `% of stephensi range` = c(pct_ste(m1), pct_ste(m12)),
      check.names = FALSE
    )
  })
  
  output$pop_table <- renderTable({
    tr  <- tier_uc4()
    m1  <- ifel(tr == 1, 1, NA)
    m12 <- ifel(tr <= 2 & !is.na(tr), 1, NA)
    pop <- function(m) {
      p <- as.numeric(global(mask(POP_ALIGNED, m), "sum", na.rm=TRUE)[1,1])
      format(round(if (!is.finite(p)) 0 else p), big.mark=",")
    }
    data.frame(
      Tier = c("Tier 1", "Tier 1 + 2 (cumulative)"),
      `Population at risk` = c(pop(m1), pop(m12)),
      check.names = FALSE
    )
  })
}

shinyApp(ui, server)
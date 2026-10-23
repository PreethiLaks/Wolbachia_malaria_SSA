# Shiny app for Use Case 4: Prevention of An. stephensi re-emergence
# Urban filter (GHS >= 21) applied first; stephensi occurrence and PfPR scored 1-3

library(shiny)
library(terra)
library(sf)

# Data folder must be in the same directory as this script
data_dir <- "data"

to01 <- function(r) {
  mx <- suppressWarnings(as.numeric(global(r, "max", na.rm = TRUE)))
  r1 <- if (is.finite(mx) && mx > 1) r / 100 else r
  clamp(r1, 0, 1, values = TRUE)
}

fast_align <- function(r, template, method = "bilinear") {
  if (!same.crs(r, template)) project(r, template, method = method)
  else resample(r, template, method = method)
}

cat("Loading pre-processed data...\n")

# Load lightweight pre-processed files
bounds          <- readRDS(file.path(data_dir, "boundaries.rds"))
ADMIN0_A        <- vect(bounds$ADMIN0_sf)
ADMIN1_A        <- vect(bounds$ADMIN1_sf)
AFRICA0_A       <- vect(bounds$AFRICA0_sf)
AFRICA1_A       <- vect(bounds$AFRICA1_sf)
COUNTRY_CHOICES <- bounds$COUNTRY_CHOICES
afro_sf         <- bounds$afro_sf

PFPR_ALIGNED  <- rast(file.path(data_dir, "pfpr_aligned.tif"))
POP_ALIGNED   <- rast(file.path(data_dir, "pop_aligned.tif"))
CELL_AREA_KM2 <- rast(file.path(data_dir, "cell_area_km2.tif"))
TEMPLATE      <- PFPR_ALIGNED

# Ensure all boundary vectors share the raster CRS
if (!same.crs(ADMIN0_A,  TEMPLATE)) ADMIN0_A  <- project(ADMIN0_A,  crs(TEMPLATE))
if (!same.crs(ADMIN1_A,  TEMPLATE)) ADMIN1_A  <- project(ADMIN1_A,  crs(TEMPLATE))
if (!same.crs(AFRICA0_A, TEMPLATE)) AFRICA0_A <- project(AFRICA0_A, crs(TEMPLATE))
if (!same.crs(AFRICA1_A, TEMPLATE)) AFRICA1_A <- project(AFRICA1_A, crs(TEMPLATE))

# UC4 specific — stephensi and GHS
STE_ALIGNED   <- rast(file.path(data_dir, "stephensi_aligned.tif"))
GHS_ALIGNED   <- rast(file.path(data_dir, "ghs_aligned.tif"))

ghs_min    <- 21
ste_floor  <- 0.10

# Stephensi range denominator: full range where stephensi is present and population > 0
# (not restricted to urban — the urban filter is for the intervention map, not the denominator)
ste_range  <- ifel(!is.na(STE_ALIGNED) & STE_ALIGNED >= ste_floor &
                     !is.na(POP_ALIGNED) & POP_ALIGNED > 0, 1, NA)
ste_range  <- mask(ste_range, ADMIN0_A)
STE_URBAN_AREA <- as.numeric(global(mask(CELL_AREA_KM2, ste_range),
                                    "sum", na.rm=TRUE)[1,1])

cat("Data loaded.\n")

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
    .tier-summary { font-size: 0.82em; color: #444; margin-top: 4px; margin-bottom: 10px; }
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
                 selectInput("ste_preset", "Threshold preset:",
                             choices  = c("Strict"             = "strict",
                                          "Moderate (default)" = "moderate",
                                          "Relaxed"            = "relaxed",
                                          "Custom"             = "custom"),
                             selected = "moderate"),
                 sliderInput("ste_t1", "Score 1 minimum occurrence (%)",
                             min = 30, max = 90, value = 70, step = 5),
                 sliderInput("ste_t2", "Score 2 minimum occurrence (%)",
                             min = 10, max = 80, value = 40, step = 5),
                 div(class = "tier-summary",
                     tags$span(style="color:#1a9641; font-weight:bold;", "Score 1 (best): "),
                     textOutput("ste_score1_label", inline = TRUE), " | ",
                     tags$span(style="color:#b8860b; font-weight:bold;", "Score 2: "),
                     textOutput("ste_score2_label", inline = TRUE), " | ",
                     tags$span(style="color:#d7191c; font-weight:bold;", "Score 3: "),
                     textOutput("ste_score3_label", inline = TRUE)
                 ),
                 
                 hr(),
                 h4("Malaria transmission (PfPR)"),
                 helpText("Moderate PfPR = best. Score 1 = sweet spot; Score 2 = 5% to lower bound; Score 3 = <5%. Cells above upper bound excluded (high PfPR likely driven by native vectors)."),
                 selectInput("pfpr_preset", "Threshold preset:",
                             choices  = c("Strict"             = "strict",
                                          "Moderate (default)" = "moderate",
                                          "Relaxed"            = "relaxed",
                                          "Custom"             = "custom"),
                             selected = "moderate"),
                 sliderInput("pfpr_lo", "PfPR lower bound for Score 1 (%)",
                             min = 5, max = 30, value = 15, step = 1),
                 sliderInput("pfpr_hi", "PfPR upper bound for Score 1 (%)",
                             min = 20, max = 60, value = 40, step = 1),
                 div(class = "tier-summary",
                     tags$span(style="color:#1a9641; font-weight:bold;", "Score 1 (sweet spot): "),
                     textOutput("pfpr_score1_label", inline = TRUE), " | ",
                     tags$span(style="color:#b8860b; font-weight:bold;", "Score 2 (buffer): "),
                     textOutput("pfpr_score2_label", inline = TRUE), " | ",
                     tags$span(style="color:#d7191c; font-weight:bold;", "Score 3 / excluded: "),
                     textOutput("pfpr_score3_label", inline = TRUE)
                 ),
                 
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
  
  # Preset values for stephensi and PfPR
  ste_presets  <- list(strict=c(80,55), moderate=c(70,40), relaxed=c(60,30))
  pfpr_presets <- list(strict=c(20,35), moderate=c(15,40), relaxed=c(10,50))
  
  # Update stephensi sliders when preset changes
  observeEvent(input$ste_preset, {
    if (input$ste_preset != "custom") {
      vals <- ste_presets[[input$ste_preset]]
      updateSliderInput(session, "ste_t1", value = vals[1])
      updateSliderInput(session, "ste_t2", value = vals[2])
    }
  })
  
  # Update PfPR sliders when preset changes
  observeEvent(input$pfpr_preset, {
    if (input$pfpr_preset != "custom") {
      vals <- pfpr_presets[[input$pfpr_preset]]
      updateSliderInput(session, "pfpr_lo", value = vals[1])
      updateSliderInput(session, "pfpr_hi", value = vals[2])
    }
  })
  
  # Tier summary labels — stephensi
  output$ste_score1_label <- renderText({ paste0("≥", input$ste_t1, "%") })
  output$ste_score2_label <- renderText({ paste0(input$ste_t2, "–", input$ste_t1, "%") })
  output$ste_score3_label <- renderText({ paste0("<", input$ste_t2, "%") })
  
  # Tier summary labels — PfPR
  output$pfpr_score1_label <- renderText({ paste0(input$pfpr_lo, "–", input$pfpr_hi, "%") })
  output$pfpr_score2_label <- renderText({ paste0("5–", input$pfpr_lo, "%") })
  output$pfpr_score3_label <- renderText({ paste0("<5% or >", input$pfpr_hi, "% (excluded)") })
  
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
    # Only cells with PfPR <= pfpr_hi are eligible (>pfpr_hi excluded — likely driven by native vectors)
    pfpr_ok <- !is.na(PFPR_ALIGNED) & PFPR_ALIGNED <= pfpr_hi
    # Score 3: PfPR < 5% (too low — limited Wolbachia impact)
    tp[urban & ste_ok & pfpr_ok & PFPR_ALIGNED < 0.05] <- 3
    # Score 2: PfPR 5% to pfpr_lo (suboptimal but eligible)
    tp[urban & ste_ok & pfpr_ok & PFPR_ALIGNED >= 0.05 & PFPR_ALIGNED < pfpr_lo] <- 2
    # Score 1: PfPR pfpr_lo to pfpr_hi (optimal sweet spot)
    tp[urban & ste_ok & pfpr_ok & PFPR_ALIGNED >= pfpr_lo & PFPR_ALIGNED <= pfpr_hi] <- 1
    
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
      ctry_v <- if (!same.crs(ctry_v, TEMPLATE)) project(ctry_v, crs(TEMPLATE)) else ctry_v
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
    tr   <- tier_uc4()
    n_afro <- nrow(ADMIN0_A)
    m1   <- ifel(tr == 1, 1, NA)
    m12  <- ifel(tr <= 2 & !is.na(tr), 1, NA)
    m123 <- ifel(!is.na(tr), 1, NA)
    
    area_raw <- function(m) {
      a <- as.numeric(global(mask(CELL_AREA_KM2, m), "sum", na.rm=TRUE)[1,1])
      if (!is.finite(a)) 0 else a
    }
    area_fmt <- function(m) {
      a <- area_raw(m)
      paste0(format(round(a), big.mark=","),
             " (", formatC(100 * a / STE_URBAN_AREA, format="f", digits=1), "% of An. stephensi range)")
    }
    cnt_fmt <- function(m) {
      cid <- rasterize(ADMIN0_A, m, field = seq_len(nrow(ADMIN0_A)), touches = FALSE)
      hit <- !is.na(cid) & !is.na(m)
      n   <- length(unique(na.omit(values(ifel(hit, cid, NA)))))
      paste0(n, " / ", n_afro, " (", round(100 * n / n_afro, 1), "%)")
    }
    data.frame(
      Tier      = c("Tier 1", "Tier 1 + 2 (cumulative)", "Tier 1 + 2 + 3 (cumulative)"),
      Countries = c(cnt_fmt(m1), cnt_fmt(m12), cnt_fmt(m123)),
      `Area km² (% of An. stephensi range)` = c(area_fmt(m1), area_fmt(m12), area_fmt(m123)),
      check.names = FALSE
    )
  })
  
  output$pop_table <- renderTable({
    tr   <- tier_uc4()
    m1   <- ifel(tr == 1, 1, NA)
    m12  <- ifel(tr <= 2 & !is.na(tr), 1, NA)
    m123 <- ifel(!is.na(tr), 1, NA)
    pop_total <- as.numeric(global(mask(POP_ALIGNED,
                                        ifel(!is.na(POP_ALIGNED) & POP_ALIGNED > 0, 1, NA)),
                                   "sum", na.rm=TRUE)[1,1])
    pop_fmt <- function(m) {
      p <- as.numeric(global(mask(POP_ALIGNED, m), "sum", na.rm=TRUE)[1,1])
      if (!is.finite(p)) p <- 0
      paste0(format(round(p / 1e6, 2), nsmall = 2), "M",
             " (", formatC(100 * p / pop_total, format="f", digits=1), "% SSA)")
    }
    data.frame(
      Tier = c("Tier 1", "Tier 1 + 2 (cumulative)", "Tier 1 + 2 + 3 (cumulative)"),
      `Urban population at risk` = c(pop_fmt(m1), pop_fmt(m12), pop_fmt(m123)),
      check.names = FALSE
    )
  })
}

shinyApp(ui, server)
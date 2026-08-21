# Shiny app for Use Case 2: Routine implementation
# Additive scoring: each parameter scored 1-3, sum determines tier (T1=3-4, T2=5-7, T3=8-9)
# High ITN and high PfPR = Score 1 (opposite to UC1a)

library(shiny)
library(terra)
library(sf)
library(plotly)

data_dir <- "data"

paths <- list(
  ara  = file.path(data_dir, "An_arabiensis.tif"),
  col  = file.path(data_dir, "An_coluzzii.tif"),
  mou  = file.path(data_dir, "An_moucheti.tif"),
  gam  = file.path(data_dir, "An_gambiae.tif"),
  ste  = file.path(data_dir, "An_stephensi.tif"),
  pfpr = file.path(data_dir, "PfPR_mean.tif"),
  itn  = file.path(data_dir, "ITN_use_rate.tif"),
  inc  = file.path(data_dir, "Pf_Incidence_mean_2000.tif"),
  pop  = file.path(data_dir, "POP_MEAN_2000_2020_5km.tif"),
  gadm = file.path(data_dir, "gadm_ssa.gpkg")
)

# Helper functions
to01 <- function(r) {
  mx <- suppressWarnings(as.numeric(global(r, "max", na.rm = TRUE)))
  r1 <- if (is.finite(mx) && mx > 1) r / 100 else r
  clamp(r1, 0, 1, values = TRUE)
}

safe_div <- function(num, den, eps = 1e-9) num / (den + eps)

fast_align <- function(r, template, method = "bilinear") {
  if (!identical(crs(r), crs(template))) project(r, template, method = method)
  else resample(r, template, method = method)
}

cat("Loading data...\n")

# SSA country list for analysis boundary
ssa_iso <- c(
  "AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM","COG","CIV","COD","GNQ","ERI",
  "SWZ","ETH","GAB","GMB","GHA","GIN","GNB","KEN","LSO","LBR","MDG","MWI","MLI","MRT","MUS","MOZ",
  "NAM","NER","NGA","RWA","STP","SEN","SYC","SLE","ZAF","SSD","TGO","UGA","TZA","ZMB","ZWE",
  "SDN","SOM","DJI"
)

# Full Africa for map display context
africa_iso <- unique(c(ssa_iso, "MAR","TUN","LBY","EGY","DZA"))

a0       <- st_read(paths$gadm, layer = "ADM_0", quiet = TRUE)
a1       <- st_read(paths$gadm, layer = "ADM_1", quiet = TRUE)

# Analysis boundary — sub-Saharan Africa only
ADMIN0   <- vect(a0[a0$GID_0 %in% ssa_iso, ])
ADMIN1   <- vect(a1[a1$GID_0 %in% ssa_iso, ])

# Named vector: country name -> ISO code for zoom dropdown
afro_sf       <- a0[a0$GID_0 %in% ssa_iso, ]
# GADM uses COUNTRY or NAME_0 column for country names
name_col      <- if ("COUNTRY" %in% names(afro_sf)) "COUNTRY" else "NAME_0"
COUNTRY_CHOICES <- c("None (full Africa)" = "",
                     setNames(sort(afro_sf$GID_0),
                              afro_sf[[name_col]][order(afro_sf$GID_0)]))

# Display boundary — all Africa (for map background context)
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

species_rasters <- list(
  arabiensis = rast(paths$ara),
  coluzzii   = rast(paths$col),
  moucheti   = rast(paths$mou),
  gambiae    = rast(paths$gam),
  stephensi  = rast(paths$ste)
)
SPECIES_ALIGNED <- lapply(species_rasters, function(r) fast_align(to01(r), TEMPLATE))

cat("Data loaded.\n")

ui <- fluidPage(
  
  tags$head(tags$style(HTML("
    body { font-family: 'Georgia', serif; background: #fafaf8; color: #2c2c2c; }
    .well { background: #f0ede8; border: none; border-radius: 8px; }
    h3 { color: #2c5f2e; font-weight: bold; border-bottom: 2px solid #2c5f2e; padding-bottom: 6px; }
    h4 { color: #4a4a4a; margin-top: 18px; }
    .method-box { background: #e8f4e8; border-left: 4px solid #2c5f2e;
                  padding: 10px 14px; border-radius: 4px; font-size: 0.88em;
                  margin-bottom: 12px; }
    .tier-legend { display: flex; gap: 12px; margin: 8px 0; }
    .tier-chip { padding: 4px 10px; border-radius: 12px; font-size: 0.82em;
                 font-weight: bold; color: white; }
    .tier-summary { font-size: 0.82em; color: #444; margin-top: 4px; margin-bottom: 10px; }
  "))),
  
  titlePanel(
    div(
      h3("Use case 2 — Routine implementation"),
      p(style = "color:#666; font-size:0.9em; margin-top:-8px;",
        "Spatial eligibility for Wolbachia routine implementation | ",
        em("An. arabiensis, An. coluzzii, An. moucheti, An. gambiae s.s."))
    )
  ),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 
                 div(class = "method-box",
                     strong("How this works:"), br(),
                     "Each parameter is scored 1 (optimal) to 3 (suboptimal). ",
                     "Scores are summed to assign each pixel a tier: ",
                     "Tier 1 (score 3–4) = highest priority, Tier 2 (score 5–7), Tier 3 (score 8–9)."
                 ),
                 
                 h4("Species"),
                 selectInput("species", NULL,
                             choices  = c("An. arabiensis" = "arabiensis",
                                          "An. coluzzii"   = "coluzzii",
                                          "An. moucheti"   = "moucheti",
                                          "An. gambiae s.s." = "gambiae"),
                             selected = "arabiensis"),
                 
                 hr(),
                 
                 # Dominance
                 h4("Vector dominance"),
                 helpText("Higher dominance = better for implementation. Tier 1 requires the highest dominance."),
                 selectInput("dom_preset", "Threshold preset:",
                             choices  = c("Strict (default)" = "strict",
                                          "Moderate"         = "moderate",
                                          "Relaxed"          = "relaxed",
                                          "Custom"           = "custom"),
                             selected = "strict"),
                 sliderInput("dom_t1", "Tier 1 minimum dominance (%)",
                             min = 10, max = 95, value = 70, step = 5),
                 sliderInput("dom_t2", "Tier 2 minimum dominance (%)",
                             min = 5,  max = 90, value = 40, step = 5),
                 div(class = "tier-summary",
                     tags$span(style="color:#1a9641; font-weight:bold;", "Tier 1: "),
                     textOutput("dom_tier1_label", inline = TRUE), " | ",
                     tags$span(style="color:#b8860b; font-weight:bold;", "Tier 2: "),
                     textOutput("dom_tier2_label", inline = TRUE), " | ",
                     tags$span(style="color:#d7191c; font-weight:bold;", "Tier 3: "),
                     textOutput("dom_tier3_label", inline = TRUE)
                 ),
                 
                 hr(),
                 
                 # ITN
                 h4("ITN use rate"),
                 helpText("Higher ITN = better for implementation. Tier 1 requires the highest ITN coverage."),
                 selectInput("itn_preset", "Threshold preset:",
                             choices  = c("Strict (default)" = "strict",
                                          "Moderate"         = "moderate",
                                          "Relaxed"          = "relaxed",
                                          "Custom"           = "custom"),
                             selected = "strict"),
                 sliderInput("itn_t1", "Tier 1 maximum ITN use rate (%)",
                             min = 20, max = 90, value = 60, step = 5),
                 sliderInput("itn_t2", "Tier 2 maximum ITN use rate (%)",
                             min = 25, max = 95, value = 80, step = 5),
                 div(class = "tier-summary",
                     tags$span(style="color:#1a9641; font-weight:bold;", "Tier 1: "),
                     textOutput("itn_tier1_label", inline = TRUE), " | ",
                     tags$span(style="color:#b8860b; font-weight:bold;", "Tier 2: "),
                     textOutput("itn_tier2_label", inline = TRUE), " | ",
                     tags$span(style="color:#d7191c; font-weight:bold;", "Tier 3: "),
                     textOutput("itn_tier3_label", inline = TRUE)
                 ),
                 
                 hr(),
                 
                 # PfPR
                 h4("Malaria transmission (PfPR)"),
                 helpText("Higher PfPR = better for implementation. Tier 1 requires the highest transmission."),
                 selectInput("pfpr_preset", "Threshold preset:",
                             choices  = c("Strict (default)" = "strict",
                                          "Moderate"         = "moderate",
                                          "Relaxed"          = "relaxed",
                                          "Custom"           = "custom"),
                             selected = "strict"),
                 sliderInput("pfpr_t1", "Tier 1 maximum PfPR (%)",
                             min = 5, max = 40, value = 15, step = 1),
                 sliderInput("pfpr_t2", "Tier 2 maximum PfPR (%)",
                             min = 10, max = 60, value = 40, step = 1),
                 div(class = "tier-summary",
                     tags$span(style="color:#1a9641; font-weight:bold;", "Tier 1: "),
                     textOutput("pfpr_tier1_label", inline = TRUE), " | ",
                     tags$span(style="color:#b8860b; font-weight:bold;", "Tier 2: "),
                     textOutput("pfpr_tier2_label", inline = TRUE), " | ",
                     tags$span(style="color:#d7191c; font-weight:bold;", "Tier 3: "),
                     textOutput("pfpr_tier3_label", inline = TRUE)
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
                         plotOutput("map_uc2", height = 650)
                ),
                tabPanel("Impact & Coverage",
                         br(),
                         h4("Coverage"),
                         tableOutput("cov_table"),
                         hr(),
                         h4("Burden"),
                         tableOutput("burden_table")
                ),
                tabPanel("Sensitivity",
                         br(),
                         div(style = "max-width: 300px;",
                             selectInput("sens_param", "Vary which parameter?",
                                         choices = c(
                                           "ITN use rate (Score 1 boundary)"    = "itn",
                                           "Dominance (Score 1 boundary)"       = "dom",
                                           "PfPR (Score 1 boundary)"            = "pfpr"
                                         ),
                                         selected = "itn"
                             )
                         ),
                         helpText("Other parameters are fixed at their current values from the Map tab.
                    The dashed line marks the current baseline."),
                         plotly::plotlyOutput("sens_plot", height = 500)
                )
              )
    )
  )
)

server <- function(input, output, session) {
  
  # Preset values for each parameter
  dom_presets  <- list(strict=c(70,40), moderate=c(50,25), relaxed=c(30,15))
  itn_presets  <- list(strict=c(80,60), moderate=c(70,50), relaxed=c(60,40))
  pfpr_presets <- list(strict=c(40,15), moderate=c(30,10), relaxed=c(20,5))
  
  # Update dominance sliders when preset changes
  observeEvent(input$dom_preset, {
    if (input$dom_preset != "custom") {
      vals <- dom_presets[[input$dom_preset]]
      updateSliderInput(session, "dom_t1", value = vals[1])
      updateSliderInput(session, "dom_t2", value = vals[2])
    }
  })
  
  # Update ITN sliders when preset changes
  observeEvent(input$itn_preset, {
    if (input$itn_preset != "custom") {
      vals <- itn_presets[[input$itn_preset]]
      updateSliderInput(session, "itn_t1", value = vals[1])
      updateSliderInput(session, "itn_t2", value = vals[2])
    }
  })
  
  # Update PfPR sliders when preset changes
  observeEvent(input$pfpr_preset, {
    if (input$pfpr_preset != "custom") {
      vals <- pfpr_presets[[input$pfpr_preset]]
      updateSliderInput(session, "pfpr_t1", value = vals[1])
      updateSliderInput(session, "pfpr_t2", value = vals[2])
    }
  })
  
  # Tier summary labels — dominance
  output$dom_tier1_label <- renderText({ paste0("\u2265", input$dom_t1, "%") })
  output$dom_tier2_label <- renderText({ paste0(input$dom_t2, "\u2013", input$dom_t1, "%") })
  output$dom_tier3_label <- renderText({ paste0("<", input$dom_t2, "%") })
  
  # Tier summary labels — ITN
  output$itn_tier1_label <- renderText({ paste0("\u2265", input$itn_t1, "%") })
  output$itn_tier2_label <- renderText({ paste0(input$itn_t2, "\u2013", input$itn_t1, "%") })
  output$itn_tier3_label <- renderText({ paste0("<", input$itn_t2, "%") })
  
  # Tier summary labels — PfPR
  output$pfpr_tier1_label <- renderText({ paste0("\u2265", input$pfpr_t1, "%") })
  output$pfpr_tier2_label <- renderText({ paste0(input$pfpr_t2, "\u2013", input$pfpr_t1, "%") })
  output$pfpr_tier3_label <- renderText({ paste0("<", input$pfpr_t2, "%") })
  
  dom_layers <- reactive({
    sp      <- input$species
    base_p  <- SPECIES_ALIGNED[[sp]]; base_p[is.na(base_p)] <- 0
    others  <- rast(lapply(setdiff(names(SPECIES_ALIGNED), sp),
                           function(nm) { r <- SPECIES_ALIGNED[[nm]]; r[is.na(r)] <- 0; r }))
    dom     <- safe_div(base_p, base_p + app(others, sum))
    sp_raw  <- SPECIES_ALIGNED[[sp]]
    data_mask <- ifel(!is.na(PFPR_ALIGNED) & !is.na(ITN_ALIGNED) &
                        !is.na(sp_raw) & sp_raw >= 0.10, 1, NA)
    list(dom = dom, sp_raw = sp_raw, data_mask = data_mask)
  }) |> bindCache(input$species)
  
  tier_weighted <- reactive({
    dl <- dom_layers()
    dm <- dl$dom; data_mask <- dl$data_mask
    
    dom_t1  <- as.numeric(input$dom_t1) / 100
    dom_t2  <- as.numeric(input$dom_t2) / 100
    itn_k1  <- input$itn_t1 / 100
    itn_k2  <- input$itn_t2 / 100
    pfpr_k1 <- input$pfpr_t1 / 100
    pfpr_k2 <- input$pfpr_t2 / 100
    
    # T1: dom >= dom_t1 | T2: dom_t2 to dom_t1 | T3: < dom_t2
    # Dynamic T2: if T1 drops to <= dom_t2 (0.40), shift T2 down to 0.25
    dom_t2_eff <- if (dom_t1 <= dom_t2) 0.25 else dom_t2
    valid <- !is.na(data_mask)
    dm_v    <- values(dm); valid_v <- values(data_mask)
    tv_t    <- rep(NA_real_, length(dm_v))
    tv_t[!is.na(valid_v) & !is.na(dm_v) & dm_v >= dom_t1] <- 1
    tv_t[!is.na(valid_v) & !is.na(dm_v) & dm_v >= dom_t2_eff & dm_v < dom_t1] <- 2
    tv_t[!is.na(valid_v) & !is.na(dm_v) & dm_v <  dom_t2_eff] <- 3
    tv <- rast(TEMPLATE); values(tv) <- tv_t
    
    # T1: < itn_k1 | T2: itn_k1 to itn_k2 | T3: >= itn_k2
    ti <- rast(TEMPLATE); values(ti) <- NA
    ti[!is.na(data_mask) & ITN_ALIGNED >= itn_k1] <- 1
    ti[!is.na(data_mask) & ITN_ALIGNED >= itn_k2 & ITN_ALIGNED < itn_k1] <- 2
    ti[!is.na(data_mask) & ITN_ALIGNED <  itn_k2] <- 3
    
    # T1: < pfpr_k1 | T2: pfpr_k1 to pfpr_k2 | T3: >= pfpr_k2
    tp <- rast(TEMPLATE); values(tp) <- NA
    tp[!is.na(data_mask) & PFPR_ALIGNED >= pfpr_k1] <- 1
    tp[!is.na(data_mask) & PFPR_ALIGNED >= pfpr_k2 & PFPR_ALIGNED < pfpr_k1] <- 2
    tp[!is.na(data_mask) & PFPR_ALIGNED <  pfpr_k2] <- 3
    
    # Each criterion scores 1 (best) to 3 (worst)
    # Total score range: 3 (all T1) to 9 (all T3)
    # Score 3-4 → Tier 1 | Score 5-7 → Tier 2 | Score 8-9 → Tier 3
    score <- tv + ti + tp   # sum of dominance + ITN + PfPR tier scores
    
    out <- rast(TEMPLATE); values(out) <- NA
    out[valid & score >= 3 & score <= 4] <- 1
    out[valid & score >= 5 & score <= 7] <- 2
    out[valid & score >= 8 & score <= 9] <- 3
    out <- mask(out, ADMIN0_A)
    score_m <- mask(score, ADMIN0_A)
    list(tier = out, score = score_m)
  })
  
  # score_r: continuous score 3-9 | tier_r: discrete 1/2/3 for tables/legend
  plot_tier <- function(tier_r, title_txt) {
    vals <- sort(unique(na.omit(as.vector(values(tier_r)))))
    vals <- vals[vals %in% 1:3]
    if (length(vals) == 0) { plot.new(); title("No eligible cells"); return(invisible()) }
    
    tier_cols   <- c("1" = "#1a9641", "2" = "#f4d03f", "3" = "#d7191c")
    tier_breaks <- c(0.5, 1.5, 2.5, 3.5)
    
    af   <- ext(AFRICA0_A)
    xlim <- c(xmin(af), xmax(af)); ylim <- c(ymin(af), ymax(af))
    
    # Country zoom takes priority over auto-zoom
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
    e_view   <- ext(xlim[1],xlim[2],ylim[1],ylim[2])
    af0_view <- crop(AFRICA0_A, e_view)
    v0_view  <- crop(ADMIN0_A, e_view)
    v1_view  <- crop(ADMIN1_A, e_view)
    tier_view <- crop(tier_r, e_view)
    
    par(mar = c(1,1,2.5,1))
    plot(af0_view, col="#e4e4e4", border=NA, axes=FALSE, xlim=xlim, ylim=ylim,
         main=title_txt, cex.main=1.0)
    plot(v0_view, add=TRUE, col="#f0f0f0", border=NA)
    
    if (length(vals) > 0) {
      # Use fixed 3-colour palette with full breaks so colours always match tiers
      plot(tier_view, add=TRUE,
           col    = c("#1a9641", "#f4d03f", "#d7191c"),
           breaks = c(0.5, 1.5, 2.5, 3.5),
           legend = FALSE, axes = FALSE)
    }
    
    plot(v1_view, add=TRUE, lwd=0.3, col=NA, border="#cccccc")
    plot(v0_view, add=TRUE, lwd=0.6, col=NA, border="#888888")
    plot(af0_view, add=TRUE, lwd=0.4, col=NA, border="#aaaaaa")
    
    legend("bottomleft",
           legend = c("Tier 1", "Tier 2", "Tier 3"),
           fill   = c("#1a9641", "#f4d03f", "#d7191c"),
           border = "grey40",
           bty    = "o", bg = "white", box.col = "grey80",
           cex    = 0.85, xpd = TRUE, title.font = 2)
  }
  
  output$map_uc2 <- renderPlot({ res <- tier_weighted(); plot_tier(res$tier, "") })
  
  output$cov_table <- renderTable({
    res  <- tier_weighted(); tr <- res$tier
    n_afro <- nrow(ADMIN0_A)   # total sub-Saharan Africa countries = denominator
    
    area <- function(m) {
      a <- as.numeric(global(mask(CELL_AREA_KM2, m), "sum", na.rm=TRUE)[1,1])
      format(round(if (is.finite(a)) a else 0), big.mark=",")
    }
    cnt <- function(m) {
      cid <- rasterize(ADMIN0_A, m, field=seq_len(nrow(ADMIN0_A)), touches=TRUE)
      hit <- !is.na(cid) & !is.na(m)
      length(unique(na.omit(values(ifel(hit, cid, NA)))))
    }
    cnt_fmt <- function(m) {
      n <- cnt(m)
      paste0(n, " / ", n_afro, " (", round(100 * n / n_afro, 1), "%)")
    }
    
    m1   <- ifel(tr == 1, 1, NA)
    m2   <- ifel(tr == 2, 1, NA)
    m3   <- ifel(tr == 3, 1, NA)
    m12  <- ifel(tr <= 2 & !is.na(tr), 1, NA)
    m123 <- ifel(!is.na(tr), 1, NA)
    
    data.frame(
      Tier        = c("Tier 1",
                      "Tier 1 + 2 (cumulative)",
                      "Tier 1 + 2 + 3 (cumulative)"),
      Countries   = c(cnt_fmt(m1), cnt_fmt(m12), cnt_fmt(m123)),
      "Area (km2)" = c(area(m1), area(m12), area(m123)),
      check.names = FALSE
    )
  })
  
  output$burden_table <- renderTable({
    res <- tier_weighted(); tr <- res$tier
    bur_valid <- mask(BURDEN_GLOBAL,
                      ifel(!is.na(INC_ALIGNED) & !is.na(POP_ALIGNED) &
                             POP_ALIGNED>0 & INC_ALIGNED>1e-12, 1, NA))
    afro_bur  <- as.numeric(global(bur_valid, "sum", na.rm=TRUE)[1,1])
    bur <- function(m) {
      v <- as.numeric(global(mask(bur_valid,m),"sum",na.rm=TRUE)[1,1])
      if (!is.finite(v)) v <- 0
      paste0(format(round(v), big.mark=","),
             " (", formatC(100*v/afro_bur, format="f", digits=2), "% sub-Saharan Africa)")
    }
    m1   <- ifel(tr == 1, 1, NA)
    m12  <- ifel(tr <= 2 & !is.na(tr), 1, NA)
    m123 <- ifel(!is.na(tr), 1, NA)
    data.frame(
      Tier = c("Tier 1",
               "Tier 1 + 2 (cumulative)",
               "Tier 1 + 2 + 3 (cumulative)"),
      `Mean annual burden` = c(bur(m1), bur(m12), bur(m123)),
      check.names = FALSE
    )
  })
  
  # Sensitivity — compute on demand when tab is opened
  sens_data <- reactive({
    param     <- input$sens_param
    sp_names  <- c("arabiensis","coluzzii","moucheti","gambiae")
    sp_labels <- c("An. arabiensis","An. coluzzii","An. moucheti","An. gambiae s.s.")
    
    def_dom_t1  <- 0.50; def_dom_t2  <- 0.25
    def_itn_t1  <- 0.80; def_itn_t2  <- 0.60  # UC2: high ITN = Score 1
    def_pfpr_t1 <- 0.40; def_pfpr_t2 <- 0.15  # UC2: high PfPR = Score 1
    pres_floor  <- 0.10
    
    thresholds <- switch(param,
                         itn  = seq(0.50, 0.95, by = 0.05),
                         dom  = seq(0.20, 0.80, by = 0.05),
                         pfpr = seq(0.20, 0.70, by = 0.05)
    )
    
    # Burden valid mask
    bur_mask  <- ifel(!is.na(INC_ALIGNED) & !is.na(POP_ALIGNED) &
                        POP_ALIGNED > 0 & INC_ALIGNED > 1e-12, 1, NA)
    bur_valid <- mask(BURDEN_GLOBAL, bur_mask)
    
    compute_t1_burden <- function(tv_v, ti_v, tp_v) {
      sc     <- tv_v + ti_v + tp_v
      out    <- rep(NA_real_, length(sc))
      out[!is.na(sc) & sc >= 3 & sc <= 4] <- 1
      tier_r <- rast(TEMPLATE); values(tier_r) <- out
      tier_r <- mask(tier_r, ADMIN0_A)
      m1     <- ifel(tier_r == 1, 1, NA)
      as.numeric(global(mask(bur_valid, m1), "sum", na.rm = TRUE)[1,1])
    }
    
    itn_v  <- values(ITN_ALIGNED)
    pfpr_v <- values(PFPR_ALIGNED)
    
    rows <- list()
    withProgress(message = "Computing sensitivity curves...", value = 0, {
      total <- length(sp_names) * length(thresholds)
      done  <- 0
      for (s in seq_along(sp_names)) {
        sp      <- sp_names[s]
        base_sp <- SPECIES_ALIGNED[[sp]]
        base_p0 <- base_sp; base_p0[is.na(base_p0)] <- 0
        others  <- rast(lapply(setdiff(names(SPECIES_ALIGNED), sp),
                               function(nm) { r <- SPECIES_ALIGNED[[nm]]; r[is.na(r)] <- 0; r }))
        dom_sp   <- safe_div(base_p0, base_p0 + app(others, sum))
        dm_sp    <- values(dom_sp)
        valid_sp <- values(ifel(!is.na(PFPR_ALIGNED) & !is.na(ITN_ALIGNED) &
                                  !is.na(base_sp) & base_sp >= pres_floor, 1, NA))
        
        # UC2: high dominance = Score 1, high ITN = Score 1, high PfPR = Score 1
        tv_def <- ifelse(!is.na(valid_sp) & !is.na(dm_sp) & dm_sp >= def_dom_t1, 1,
                         ifelse(!is.na(valid_sp) & !is.na(dm_sp) & dm_sp >= def_dom_t2, 2,
                                ifelse(!is.na(valid_sp) & !is.na(dm_sp), 3, NA)))
        ti_def <- ifelse(!is.na(valid_sp) & !is.na(itn_v) & itn_v >= def_itn_t1, 1,
                         ifelse(!is.na(valid_sp) & !is.na(itn_v) & itn_v >= def_itn_t2, 2,
                                ifelse(!is.na(valid_sp) & !is.na(itn_v), 3, NA)))
        tp_def <- ifelse(!is.na(valid_sp) & !is.na(pfpr_v) & pfpr_v >= def_pfpr_t1, 1,
                         ifelse(!is.na(valid_sp) & !is.na(pfpr_v) & pfpr_v >= def_pfpr_t2, 2,
                                ifelse(!is.na(valid_sp) & !is.na(pfpr_v), 3, NA)))
        
        for (th in thresholds) {
          if (param == "itn") {
            itn_t2_eff <- min(th - 0.10, def_itn_t2)
            ti_v2 <- ifelse(!is.na(valid_sp) & !is.na(itn_v) & itn_v >= th,          1,
                            ifelse(!is.na(valid_sp) & !is.na(itn_v) & itn_v >= itn_t2_eff,  2,
                                   ifelse(!is.na(valid_sp) & !is.na(itn_v), 3, NA)))
            burden <- compute_t1_burden(tv_def, ti_v2, tp_def)
          } else if (param == "dom") {
            dom_t2_eff <- if (th <= def_dom_t2) 0.15 else def_dom_t2
            tv_v2 <- ifelse(!is.na(valid_sp) & !is.na(dm_sp) & dm_sp >= th,         1,
                            ifelse(!is.na(valid_sp) & !is.na(dm_sp) & dm_sp >= dom_t2_eff, 2,
                                   ifelse(!is.na(valid_sp) & !is.na(dm_sp), 3, NA)))
            burden <- compute_t1_burden(tv_v2, ti_def, tp_def)
          } else {
            pfpr_t2_eff <- min(th - 0.10, def_pfpr_t2)
            tp_v2 <- ifelse(!is.na(valid_sp) & !is.na(pfpr_v) & pfpr_v >= th,           1,
                            ifelse(!is.na(valid_sp) & !is.na(pfpr_v) & pfpr_v >= pfpr_t2_eff,  2,
                                   ifelse(!is.na(valid_sp) & !is.na(pfpr_v), 3, NA)))
            burden <- compute_t1_burden(tv_def, ti_def, tp_v2)
          }
          rows[[length(rows)+1]] <- data.frame(
            species = sp_labels[s], threshold = round(th * 100), burden = burden
          )
          done <- done + 1
          setProgress(done / total)
        }
      }
    })
    do.call(rbind, rows)
  })
  
  output$sens_plot <- plotly::renderPlotly({
    df    <- sens_data()
    param <- input$sens_param
    
    baseline_th <- switch(param,
                          itn  = input$itn_t1,
                          dom  = input$dom_t1,
                          pfpr = input$pfpr_t1
    )
    
    x_label <- switch(param,
                      itn  = "Minimum ITN use rate for Score 1 (%)",
                      dom  = "Minimum dominance for Score 1 (%)",
                      pfpr = "Minimum PfPR for Score 1 (%)"
    )
    
    sp_cols <- c(
      "An. arabiensis"   = "#e41a1c",
      "An. coluzzii"     = "#377eb8",
      "An. moucheti"     = "#4daf4a",
      "An. gambiae s.s." = "#984ea3"
    )
    
    # Compute absolute additional burden from baseline for each species
    df_plot <- do.call(rbind, lapply(unique(df$species), function(sp) {
      d        <- df[df$species == sp, ]
      baseline <- d$burden[which.min(abs(d$threshold - baseline_th))]
      d$additional <- round(d$burden - baseline)
      d
    }))
    
    fig <- plotly::plot_ly()
    for (sp in unique(df_plot$species)) {
      d <- df_plot[df_plot$species == sp, ]
      fig <- plotly::add_trace(fig,
                               data = d, x = ~threshold, y = ~additional,
                               type = "scatter", mode = "lines+markers",
                               name = sp,
                               line   = list(color = sp_cols[sp], width = 2),
                               marker = list(color = sp_cols[sp], size = 7),
                               text   = ~paste0(ifelse(additional >= 0, "+", ""),
                                                format(round(additional/1e6, 2), nsmall=2), "M cases/yr"),
                               hovertemplate = paste0(
                                 "<b>", sp, "</b><br>",
                                 x_label, ": %{x}%<br>",
                                 "Additional burden: %{text}<extra></extra>"
                               )
      )
    }
    
    fig <- plotly::layout(fig,
                          xaxis = list(title = x_label),
                          yaxis = list(
                            title = "Additional annual cases captured relative to baseline",
                            tickformat = ".2s"
                          ),
                          shapes = list(
                            list(type = "line", x0 = baseline_th, x1 = baseline_th,
                                 y0 = 0, y1 = 1, yref = "paper",
                                 line = list(color = "grey50", dash = "dash", width = 1.5)),
                            list(type = "line", x0 = min(df$threshold), x1 = max(df$threshold),
                                 y0 = 0, y1 = 0,
                                 line = list(color = "black", width = 1.2))
                          ),
                          annotations = list(list(
                            x = baseline_th, y = 1, yref = "paper",
                            text = "current", showarrow = FALSE,
                            font = list(color = "grey50", size = 11),
                            xanchor = "left", yanchor = "bottom"
                          )),
                          legend        = list(orientation = "h", y = -0.2),
                          hovermode     = "x unified",
                          plot_bgcolor  = "white",
                          paper_bgcolor = "white"
    )
    fig
  })
  
}

shinyApp(ui, server)
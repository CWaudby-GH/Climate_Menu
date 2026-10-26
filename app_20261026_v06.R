###############################################################################
# CLIMATE MENU - SHINY INTERFACE
# Version: 20261026_v06
#
# PURPOSE
# -------
# Thin user-facing interface over the validated Climate Menu backend.
#
# Stage 1:
#   user selections -> capability registry -> validated request
#
# Stage 2:
#   validated request -> targets backend -> results
###############################################################################

suppressPackageStartupMessages({
  library(shiny)
})

source(
  "config/Climate_Menu_Capability_Registry_20261026_v06.R",
  local = FALSE
)

source(
  "config/Config_CM_20261026_v06.R",
  local = FALSE
)


###############################################################################
# DISPLAY HELPERS
###############################################################################

source_labels <- c(
  "AUS-18 | ~20 km | Bias-adjusted" =
    "AUS-18|bias_adjusted",

  "AUS-18 | ~20 km | Raw NARCliM" =
    "AUS-18|raw",

  "SEAus-04 | ~4 km | Raw NARCliM" =
    "SEAus-04|raw",

  "SEAus-04 | ~4 km | Bias-adjusted" =
    "SEAus-04|bias_adjusted"
)

scenario_labels <- c(
  "Historical" = "historical",
  "SSP1-2.6" = "ssp126",
  "SSP2-4.5" = "ssp245",
  "SSP3-7.0" = "ssp370",
  "Compare all SSPs" = "all_ssps"
)

comparison_labels <- stats::setNames(
  climate_menu_comparison_registry$comparison_mode,
  climate_menu_comparison_registry$display_name
)

model_labels <- stats::setNames(
  climate_menu_model_registry$model_mode,
  climate_menu_model_registry$display_name
)

spatial_labels <- stats::setNames(
  climate_menu_spatial_registry$spatial_mode,
  climate_menu_spatial_registry$display_name
)


###############################################################################
# DOMAIN / GEOGRAPHY AVAILABILITY
#
# AUS-18:
#   Australia-wide Climate Menu pathway.
#
# SEAus-04:
#   Current named-GCCSA support is limited to the five capital-city
#   study areas already used in the validated SEAus-04 workflow.
#
# Future extension:
#   LGA, bbox and uploaded geometries should be validated by intersection
#   with an authoritative SEAus-04 domain footprint.
###############################################################################

climate_menu_seaus04_named_geographies <- c(
  "Greater Sydney",
  "Greater Melbourne",
  "Greater Brisbane",
  "Greater Adelaide",
  "Australian Capital Territory",
  "ACT",
  "Canberra"
)


normalise_geography_name <- function(x) {

  x <- trimws(x)
  x <- tolower(x)

  x <- gsub(
    "[^a-z0-9]+",
    " ",
    x
  )

  trimws(x)
}


climate_menu_geography_support <- function(
  domain,
  roi_mode,
  geography_type = NULL,
  geography_names = NULL
) {

  if (identical(domain, "AUS-18")) {

    return(
      list(
        available = TRUE,
        reason = NULL
      )
    )
  }


  if (!identical(domain, "SEAus-04")) {

    return(
      list(
        available = FALSE,
        reason = paste(
          "Unknown Climate Menu domain:",
          domain
        )
      )
    )
  }


  if (
    roi_mode %in%
      c(
        "named_geography",
        "multiple_geographies"
      )
  ) {

    if (!identical(geography_type, "GCCSA")) {

      return(
        list(
          available = FALSE,
          reason = paste(
            "SEAus-04 support for this boundary type",
            "has not yet been enabled in the Shiny layer."
          )
        )
      )
    }


    requested <- normalise_geography_name(
      geography_names
    )

    supported <- normalise_geography_name(
      climate_menu_seaus04_named_geographies
    )


    if (
      length(requested) > 0L &&
      all(requested %in% supported)
    ) {

      return(
        list(
          available = TRUE,
          reason = NULL
        )
      )
    }


    return(
      list(
        available = FALSE,
        reason = paste(
          "SEAus-04 (~4 km) does not support the",
          "selected Climate Menu study area.",
          "Use AUS-18 (~20 km) for this request."
        )
      )
    )
  }


  if (identical(roi_mode, "bbox")) {

    return(
      list(
        available = FALSE,
        reason = paste(
          "SEAus-04 bounding-box requests will be enabled",
          "after domain-footprint intersection is connected.",
          "Use AUS-18 for the current request."
        )
      )
    )
  }


  list(
    available = FALSE,
    reason = paste(
      "The selected study-area definition has not yet",
      "been validated for SEAus-04."
    )
  )
}


climate_menu_source_status <- function(
  source_value,
  capability,
  roi_mode,
  geography_type = NULL,
  geography_names = NULL
) {

  parts <- strsplit(
    source_value,
    "\\|"
  )[[1]]

  domain <- parts[[1]]
  product <- parts[[2]]


  source_row <-
    climate_menu_source_registry[
      climate_menu_source_registry$domain == domain &
      climate_menu_source_registry$product == product,
      ,
      drop = FALSE
    ]


  if (nrow(source_row) != 1L) {

    return(
      list(
        available = FALSE,
        reason = "This source is not registered in Climate Menu."
      )
    )
  }


  # Product may be known scientifically but not currently executable.
  if (!isTRUE(source_row$enabled[[1]])) {

    if (
      identical(domain, "SEAus-04") &&
      identical(product, "bias_adjusted")
    ) {

      return(
        list(
          available = FALSE,
          reason = paste(
            "SEAus-04 bias-adjusted data are known to exist,",
            "but the current remote Climate Menu catalogue path",
            "has not yet been resolved."
          )
        )
      )
    }


    return(
      list(
        available = FALSE,
        reason = "This Climate Menu source is not currently enabled."
      )
    )
  }


  geo <- climate_menu_geography_support(
    domain = domain,
    roi_mode = roi_mode,
    geography_type = geography_type,
    geography_names = geography_names
  )


  if (!isTRUE(geo$available)) {

    return(
      list(
        available = FALSE,
        reason = geo$reason
      )
    )
  }


  # Bias-adjusted NARCliM currently supports the registered adjusted
  # temperature/precipitation variables. Analyses whose preferred product
  # is raw should not silently request an unsupported BA variable.
  if (
    identical(product, "bias_adjusted") &&
    !identical(
      capability$preferred_product[[1]],
      "bias_adjusted"
    )
  ) {

    return(
      list(
        available = FALSE,
        reason = paste(
          "Bias-adjusted NARCliM is not currently registered",
          "for this climate exposure.",
          "Use Raw NARCliM."
        )
      )
    )
  }


  list(
    available = TRUE,
    reason = NULL
  )
}


###############################################################################
# UI
###############################################################################

ui <- fluidPage(

  tags$head(

    tags$style(
      HTML("
        body {
          background-color: #f7f8fa;
        }

        .cm-header {
          background: white;
          border-bottom: 1px solid #dfe3e8;
          padding: 18px 25px;
          margin-bottom: 20px;
        }

        .cm-title {
          margin: 0;
          font-size: 30px;
          font-weight: 600;
        }

        .cm-subtitle {
          color: #5f6b76;
          margin-top: 5px;
          margin-bottom: 0;
        }

        .cm-card {
          background: white;
          border: 1px solid #dfe3e8;
          border-radius: 7px;
          padding: 18px;
          margin-bottom: 16px;
        }

        .cm-step {
          font-size: 17px;
          font-weight: 600;
          margin-bottom: 14px;
        }

        .cm-help {
          color: #66717c;
          font-size: 13px;
          margin-top: -5px;
          margin-bottom: 12px;
        }

        .cm-status-ok {
          background: #eef7ef;
          border-left: 4px solid #4c8b55;
          padding: 12px;
          margin-top: 10px;
        }

        .cm-status-info {
          background: #eef4fa;
          border-left: 4px solid #567da8;
          padding: 12px;
          margin-top: 10px;
        }

        pre {
          background: #f5f6f7;
          border: 1px solid #e0e3e6;
        }
      ")
    )
  ),


  div(
    class = "cm-header",

    h1(
      class = "cm-title",
      "Climate Menu"
    ),

    p(
      class = "cm-subtitle",
      paste(
        "Climate data selection and analysis for",
        "climate and health research"
      )
    )
  ),


  fluidRow(

    column(

      width = 4,


      div(
        class = "cm-card",

        div(
          class = "cm-step",
          "1. Research question"
        ),

        selectInput(
          "exposure_family",
          "Climate exposure",
          choices = climate_menu_exposure_families()
        ),

        selectInput(
          "metric",
          "Analysis",
          choices = NULL
        ),

        uiOutput(
          "metric_description"
        ),

        conditionalPanel(
          condition = "output.needs_threshold",

          numericInput(
            "threshold",
            "Threshold",
            value = 35,
            step = 1
          )
        ),

        conditionalPanel(
          condition = "output.needs_percentile",

          sliderInput(
            "percentile",
            "Percentile",
            min = 50,
            max = 99,
            value = 95,
            step = 1
          )
        ),

        conditionalPanel(
          condition = "output.needs_duration",

          numericInput(
            "min_duration",
            "Minimum consecutive days",
            value = 3,
            min = 2,
            max = 30,
            step = 1
          )
        )
      ),


      div(
        class = "cm-card",

        div(
          class = "cm-step",
          "2. Study area"
        ),

        selectInput(
          "spatial_mode",
          "Study-area definition",
          choices = spatial_labels,
          selected = "named_geography"
        ),

        conditionalPanel(
          condition =
            "input.spatial_mode == 'named_geography' || input.spatial_mode == 'multiple_geographies'",

          selectInput(
            "geography_type",
            "Boundary type",
            choices = c(
              "Greater Capital City Statistical Area" = "GCCSA",
              "Local Government Area" = "LGA"
            ),
            selected = "GCCSA"
          ),

          textInput(
            "geography_names",
            "Area name(s)",
            value = "Greater Perth"
          ),

          helpText(
            "For multiple areas, separate names with semicolons."
          ),

          conditionalPanel(
            condition =
              "input.spatial_mode == 'multiple_geographies'",

            checkboxInput(
              "combine_geographies",
              "Combine selected areas into one study region",
              value = FALSE
            )
          )
        ),

        conditionalPanel(
          condition = "input.spatial_mode == 'bbox'",

          fluidRow(

            column(
              6,
              numericInput(
                "lon_min",
                "West",
                value = 115.5
              )
            ),

            column(
              6,
              numericInput(
                "lon_max",
                "East",
                value = 116.2
              )
            )
          ),

          fluidRow(

            column(
              6,
              numericInput(
                "lat_min",
                "South",
                value = -32.5
              )
            ),

            column(
              6,
              numericInput(
                "lat_max",
                "North",
                value = -31.5
              )
            )
          )
        )
      )
    ),


    column(

      width = 4,


      div(
        class = "cm-card",

        div(
          class = "cm-step",
          "3. Climate data"
        ),

        selectInput(
          "source",
          "NARCliM product",
          choices = source_labels,
          selected = "AUS-18|bias_adjusted"
        ),

        uiOutput(
          "source_availability_note"
        ),

        uiOutput(
          "source_note"
        )
      ),


      div(
        class = "cm-card",

        div(
          class = "cm-step",
          "4. Time periods"
        ),

        selectInput(
          "comparison",
          "Analysis period",
          choices = comparison_labels,
          selected = "historical_future_change"
        ),

        conditionalPanel(
          condition =
            "input.comparison != 'future'",

          fluidRow(

            column(
              6,
              numericInput(
                "hist_start",
                "Historical start",
                value = 1990,
                min = 1951,
                max = 2014
              )
            ),

            column(
              6,
              numericInput(
                "hist_end",
                "Historical end",
                value = 2009,
                min = 1951,
                max = 2014
              )
            )
          )
        ),

        conditionalPanel(
          condition =
            "input.comparison != 'historical'",

          selectInput(
            "scenario",
            "Future emissions pathway",
            choices = scenario_labels[
              names(scenario_labels) != "Historical"
            ],
            selected = "ssp245"
          ),

          fluidRow(

            column(
              6,
              numericInput(
                "future_start",
                "Future start",
                value = 2040,
                min = 2015,
                max = 2100
              )
            ),

            column(
              6,
              numericInput(
                "future_end",
                "Future end",
                value = 2059,
                min = 2015,
                max = 2100
              )
            )
          )
        )
      ),


      div(
        class = "cm-card",

        div(
          class = "cm-step",
          "5. Climate projections"
        ),

        selectInput(
          "model_mode",
          "Model selection",
          choices = model_labels,
          selected = "ensemble_all"
        ),

        conditionalPanel(
          condition = "input.model_mode == 'individual'",

          selectInput(
            "gcm",
            "Global climate model",
            choices = c(
              "ACCESS-ESM1-5",
              "EC-Earth3-Veg",
              "MPI-ESM1-2-HR",
              "NorESM2-MM",
              "UKESM1-0-LL"
            )
          ),

          selectInput(
            "rcm",
            "Regional climate model",
            choices = c(
              "WRF412R3" = "NARCliM2-0-WRF412R3",
              "WRF412R5" = "NARCliM2-0-WRF412R5"
            )
          )
        )
      )
    ),


    column(

      width = 4,


      div(
        class = "cm-card",

        div(
          class = "cm-step",
          "6. Review and run"
        ),

        uiOutput(
          "validation_status"
        ),

        br(),

        actionButton(
          "build_request",
          "Build analysis request",
          class = "btn-primary",
          width = "100%"
        )
      ),


      div(
        class = "cm-card",

        div(
          class = "cm-step",
          "Analysis request"
        ),

        verbatimTextOutput(
          "request_preview"
        )
      ),


      div(
        class = "cm-card",

        div(
          class = "cm-step",
          "Results"
        ),

        div(
          class = "cm-status-info",

          strong(
            "Backend connection is the next step."
          ),

          br(),

          span(
            paste(
              "This interface currently validates and constructs",
              "the Climate Menu request without running the",
              "NARCliM analysis."
            )
          )
        )
      )
    )
  )
)


###############################################################################
# SERVER
###############################################################################

server <- function(input, output, session) {


  ###########################################################################
  # METRIC CHOICES
  ###########################################################################

  observe({

    req(input$exposure_family)

    choices <- climate_menu_metrics(
      input$exposure_family
    )

    updateSelectInput(
      session,
      "metric",
      choices = choices,
      selected = unname(choices[[1]])
    )
  })


  selected_capability <- reactive({

    req(input$metric)

    climate_menu_capability(
      input$metric
    )
  })


  output$metric_description <- renderUI({

    x <- selected_capability()

    div(
      class = "cm-help",
      x$description[[1]]
    )
  })


  output$needs_threshold <- reactive({

    isTRUE(
      selected_capability()$threshold_required[[1]]
    )
  })

  outputOptions(
    output,
    "needs_threshold",
    suspendWhenHidden = FALSE
  )


  output$needs_percentile <- reactive({

    isTRUE(
      selected_capability()$percentile_required[[1]]
    )
  })

  outputOptions(
    output,
    "needs_percentile",
    suspendWhenHidden = FALSE
  )


  output$needs_duration <- reactive({

    isTRUE(
      selected_capability()$event_duration_required[[1]]
    )
  })

  outputOptions(
    output,
    "needs_duration",
    suspendWhenHidden = FALSE
  )


  ###########################################################################
  # SOURCE AVAILABILITY
  ###########################################################################

  current_geography_names <- reactive({

    if (
      !input$spatial_mode %in%
        c(
          "named_geography",
          "multiple_geographies"
        )
    ) {

      return(NULL)
    }


    x <- trimws(
      unlist(
        strsplit(
          input$geography_names,
          ";",
          fixed = TRUE
        )
      )
    )

    x[nzchar(x)]
  })


  source_status_table <- reactive({

    cap <- selected_capability()

    values <- unname(source_labels)

    rows <- lapply(
      values,
      function(value) {

        status <- climate_menu_source_status(
          source_value = value,
          capability = cap,
          roi_mode = input$spatial_mode,
          geography_type =
            if (
              input$spatial_mode %in%
                c(
                  "named_geography",
                  "multiple_geographies"
                )
            )
              input$geography_type
            else
              NULL,
          geography_names =
            current_geography_names()
        )

        data.frame(
          value = value,
          available = isTRUE(status$available),
          reason =
            if (is.null(status$reason))
              ""
            else
              status$reason,
          stringsAsFactors = FALSE
        )
      }
    )

    do.call(
      rbind,
      rows
    )
  })


  observe({

    status <- source_status_table()

    base_labels <- names(source_labels)
    values <- unname(source_labels)

    display_labels <- base_labels


    for (i in seq_along(values)) {

      if (!status$available[status$value == values[[i]]]) {

        display_labels[[i]] <- paste0(
          base_labels[[i]],
          " - unavailable"
        )
      }
    }


    choices <- stats::setNames(
      values,
      display_labels
    )


    selected <- input$source


    if (
      is.null(selected) ||
      !nzchar(selected) ||
      !selected %in% values ||
      !isTRUE(
        status$available[
          match(
            selected,
            status$value
          )
        ]
      )
    ) {

      cap <- selected_capability()

      preferred_value <- paste(
        "AUS-18",
        cap$preferred_product[[1]],
        sep = "|"
      )


      if (
        preferred_value %in% status$value &&
        isTRUE(
          status$available[
            match(
              preferred_value,
              status$value
            )
          ]
        )
      ) {

        selected <- preferred_value

      } else {

        valid_values <- status$value[
          status$available
        ]

        if (length(valid_values)) {

          selected <- valid_values[[1]]

        } else {

          selected <- values[[1]]
        }
      }
    }


    updateSelectInput(
      session,
      "source",
      choices = choices,
      selected = selected
    )
  })


  output$source_availability_note <- renderUI({

    status <- source_status_table()

    unavailable <- status[
      !status$available,
      ,
      drop = FALSE
    ]

    if (!nrow(unavailable)) {
      return(NULL)
    }


    labels <- names(source_labels)[
      match(
        unavailable$value,
        unname(source_labels)
      )
    ]


    div(
      class = "cm-help",

      strong(
        "Unavailable for the current request:"
      ),

      tags$ul(
        lapply(
          seq_len(nrow(unavailable)),
          function(i) {

            tags$li(
              paste0(
                labels[[i]],
                ": ",
                unavailable$reason[[i]]
              )
            )
          }
        )
      )
    )
  })


  ###########################################################################
  # SOURCE INTERPRETATION
  ###########################################################################

  selected_source <- reactive({

    req(input$source)

    parts <- strsplit(
      input$source,
      "\\|"
    )[[1]]

    list(
      domain = parts[[1]],
      product = parts[[2]]
    )
  })


  output$source_note <- renderUI({

    src <- selected_source()
    cap <- selected_capability()

    preferred <- identical(
      src$product,
      cap$preferred_product[[1]]
    )

    if (preferred) {

      div(
        class = "cm-status-ok",
        strong("Selected data product."),
        br(),
        "This product provides the source variable used for the selected analysis."
      )

    } else {

      div(
        class = "cm-status-info",
        strong("Alternative product."),
        br(),
        paste(
          "Climate Menu will use the corresponding raw",
          "NARCliM variable where supported."
        )
      )
    }
  })


  ###########################################################################
  # REQUEST CONSTRUCTION
  ###########################################################################

  build_request_object <- reactive({

    cap <- selected_capability()
    src <- selected_source()

    variable <- if (
      identical(
        src$product,
        "bias_adjusted"
      )
    ) {

      cap$preferred_variable[[1]]

    } else {

      cap$raw_fallback_variable[[1]]
    }


    geography_names <- NULL
    roi_bbox <- NULL

    if (
      input$spatial_mode %in%
        c(
          "named_geography",
          "multiple_geographies"
        )
    ) {

      geography_names <- trimws(
        unlist(
          strsplit(
            input$geography_names,
            ";",
            fixed = TRUE
          )
        )
      )

      geography_names <- geography_names[
        nzchar(geography_names)
      ]
    }


    if (identical(input$spatial_mode, "bbox")) {

      roi_bbox <- c(
        xmin = input$lon_min,
        xmax = input$lon_max,
        ymin = input$lat_min,
        ymax = input$lat_max
      )
    }


    list(

      climate_menu_version =
        "20261026_v06",

      capability_id =
        cap$capability_id[[1]],

      exposure_family =
        cap$exposure_family[[1]],

      metric_name =
        cap$display_name[[1]],

      analysis_type =
        cap$analysis_type[[1]],

      domain =
        src$domain,

      product =
        src$product,

      variable =
        variable,

      comparison =
        input$comparison,

      historical = list(
        start_year =
          if (input$comparison != "future")
            input$hist_start
          else
            NULL,

        end_year =
          if (input$comparison != "future")
            input$hist_end
          else
            NULL
      ),

      future = list(
        scenario =
          if (input$comparison != "historical") {

            if (identical(input$scenario, "all_ssps")) {

              c(
                "ssp126",
                "ssp245",
                "ssp370"
              )

            } else {

              input$scenario
            }

          } else {

            NULL
          },

        start_year =
          if (input$comparison != "historical")
            input$future_start
          else
            NULL,

        end_year =
          if (input$comparison != "historical")
            input$future_end
          else
            NULL
      ),

      model_mode =
        input$model_mode,

      gcm =
        if (input$model_mode == "individual")
          input$gcm
        else
          NULL,

      rcm =
        if (input$model_mode == "individual")
          input$rcm
        else
          NULL,

      analysis_parameters = list(

        threshold =
          if (isTRUE(cap$threshold_required[[1]]))
            input$threshold
          else
            NULL,

        percentile =
          if (isTRUE(cap$percentile_required[[1]]))
            input$percentile / 100
          else
            NULL,

        min_duration =
          if (isTRUE(cap$event_duration_required[[1]]))
            input$min_duration
          else
            NULL
      ),

      roi = list(

        mode =
          input$spatial_mode,

        geography_type =
          if (input$spatial_mode != "bbox")
            input$geography_type
          else
            NULL,

        geography_names =
          geography_names,

        combine_geographies =
          if (
            input$spatial_mode ==
              "multiple_geographies"
          )
            isTRUE(input$combine_geographies)
          else
            TRUE,

        bbox =
          roi_bbox
      )
    )
  })


  ###########################################################################
  # VALIDATION
  ###########################################################################

  validation <- reactive({

    x <- build_request_object()

    problems <- character()


    if (
      !x$domain %in%
        climate_menu_source_registry$domain[
          climate_menu_source_registry$enabled
        ]
    ) {

      problems <- c(
        problems,
        "Selected NARCliM domain is not enabled."
      )
    }


    valid_source <- with(
      climate_menu_source_registry,
      enabled &
        domain == x$domain &
        product == x$product
    )

    if (!any(valid_source)) {

      problems <- c(
        problems,
        paste(
          "The selected domain/product combination",
          "is not currently enabled."
        )
      )
    }


    selected_status <- climate_menu_source_status(
      source_value = input$source,
      capability = selected_capability(),
      roi_mode = x$roi$mode,
      geography_type = x$roi$geography_type,
      geography_names = x$roi$geography_names
    )


    if (!isTRUE(selected_status$available)) {

      problems <- c(
        problems,
        selected_status$reason
      )
    }


    if (
      !is.null(x$historical$start_year) &&
      !is.null(x$historical$end_year) &&
      x$historical$start_year >
        x$historical$end_year
    ) {

      problems <- c(
        problems,
        "Historical start year is after the end year."
      )
    }


    if (
      !is.null(x$future$start_year) &&
      !is.null(x$future$end_year) &&
      x$future$start_year >
        x$future$end_year
    ) {

      problems <- c(
        problems,
        "Future start year is after the end year."
      )
    }


    if (
      x$roi$mode %in%
        c(
          "named_geography",
          "multiple_geographies"
        ) &&
      !length(x$roi$geography_names)
    ) {

      problems <- c(
        problems,
        "At least one study-area name is required."
      )
    }


    if (
      identical(
        x$roi$mode,
        "multiple_geographies"
      ) &&
      length(x$roi$geography_names) < 2L
    ) {

      problems <- c(
        problems,
        paste(
          "Multiple named areas requires at least",
          "two semicolon-separated area names."
        )
      )
    }


    if (
      identical(x$roi$mode, "bbox")
    ) {

      b <- x$roi$bbox

      if (
        b[["xmin"]] >= b[["xmax"]] ||
        b[["ymin"]] >= b[["ymax"]]
      ) {

        problems <- c(
          problems,
          "Bounding-box coordinates are invalid."
        )
      }
    }


    list(
      valid = !length(problems),
      problems = problems
    )
  })


  output$validation_status <- renderUI({

    v <- validation()

    if (v$valid) {

      div(
        class = "cm-help",
        "Review the selections and build the analysis request."
      )

    } else {

      div(
        class = "alert alert-danger",

        strong("Please correct:"),

        tags$ul(
          lapply(
            v$problems,
            tags$li
          )
        )
      )
    }
  })


  ###########################################################################
  # BUILD BUTTON
  ###########################################################################

  final_request <- eventReactive(
    input$build_request,
    {

      v <- validation()

      validate(
        need(
          v$valid,
          paste(
            v$problems,
            collapse = "\n"
          )
        )
      )

      build_request_object()
    }
  )


  output$request_preview <- renderPrint({

    req(final_request())

    str(
      final_request(),
      max.level = 4
    )
  })
}


###############################################################################
# RUN
###############################################################################

shinyApp(
  ui = ui,
  server = server
)
# Minimal Shiny app for exploring treatment effects.
# The app expects a het_ensemble object to be provided in the environment
# as `heteffects_object`. view_effects() will load the app.

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required to run this app.")
}

ui <- shiny::fluidPage(
  shiny::titlePanel("HetEffects: Treatment Effects Explorer"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::selectInput(
        inputId = "plot_type",
        label = "Plot type",
        choices = c("weights", "cate", "diagnostics")
      )
    ),
    shiny::mainPanel(
      shiny::plotOutput("main_plot")
    )
  )
)

server <- function(input, output, session) {
  output$main_plot <- shiny::renderPlot({
    obj <- get("heteffects_object", envir = .GlobalEnv)
    plot(obj, type = input$plot_type)
  })
}

shiny::shinyApp(ui, server)

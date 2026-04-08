library(shiny)
library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Define the UI
ui <- fluidPage(
  titlePanel("Ebola Outbreak Simulator (SEIRD Model)"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Transmission Parameters"),
      sliderInput("beta", "Transmission Rate (beta):", 
                  min = 0, max = 1, value = 0.288, step = 0.01),
      sliderInput("alpha", "Incubation Rate (alpha):", 
                  min = 0.01, max = 0.5, value = 0.047, step = 0.001),
      sliderInput("gamma", "Recovery/Removal Rate (gamma):", 
                  min = 0.01, max = 0.5, value = 0.0714, step = 0.001),
      
      hr(),
      h4("Disease Severity"),
      sliderInput("p", "Case Fatality Rate (CFR):", 
                  min = 0, max = 1, value = 0.613, step = 0.01),
      
      hr(),
      h4("Simulation Settings"),
      sliderInput("days", "Days to Simulate:", 
                  min = 10, max = 1000, value = 400, step = 5)
    ),
    
    mainPanel(
      plotOutput("seirPlot", height = "500px"),
      br(),
      wellPanel(
        h4("Summary Statistics"),
        tableOutput("statTable")
      )
    )
  )
)

# 2. Define the Server Logic
server <- function(input, output) {
  
  # The SEIRD Model Function
  baseSEIRmodel <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- -beta * S * I
      dE <-  beta * S * I - alpha * E
      dI <-  alpha * E - gamma * I
      
      # Using p as the CFR (proportion who die)
      dR <-  gamma * I * (1 - p)
      dD <-  gamma * I * p
      
      # Order must be S, E, I, R, D
      return(list(c(dS, dE, dI, dR, dD)))
    })
  }
  
  # Reactive Calculation
  model_data <- reactive({
    parameters <- c(beta = input$beta, 
                    gamma = input$gamma, 
                    alpha = input$alpha, 
                    p = input$p)
    
    initial_values <- c(S = 0.9999, E = 0, I = 0.0001, R = 0, D = 0)
    times <- seq(0, input$days, by = 1)
    
    out <- ode(y = initial_values, times = times, func = baseSEIRmodel, parms = parameters)
    as.data.frame(out)
  })
  
  # Render the Plot
  output$seirPlot <- renderPlot({
    df <- model_data()
    
    df_long <- df %>%
      pivot_longer(cols = c(S, E, I, R, D), names_to = "Compartment", values_to = "Proportion")
    
    ggplot(df_long, aes(x = time, y = Proportion, color = Compartment)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("S"="#9b59b6", "E"="#f1c40f", "I"="#2ecc71", "R"="#3498db", "D"="#e74c3c")) +
      labs(title = "Epidemic Projection", x = "Days", y = "Proportion of Population") +
      theme_minimal(base_size = 16) +
      theme(legend.position = "bottom")
  })
  
  # Render Summary Table
  output$statTable <- renderTable({
    df <- model_data()
    final_row <- tail(df, 1)
    
    data.frame(
      Metric = c("Basic Reproduction Number (R0)", "Peak Infected Proportion", "Total Survival (R)", "Total Deceased (D)"),
      Value = c(
        input$beta / input$gamma,
        max(df$I),
        final_row$R,
        final_row$D
      )
    )
  }, digits = 4)
}

# 3. Launch
shinyApp(ui = ui, server = server)
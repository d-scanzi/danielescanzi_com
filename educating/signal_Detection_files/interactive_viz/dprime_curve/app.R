#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

fa  <- seq(0, 1, by=0.01)

# Functions
generate_roc_line <- function(d) {
    
    
    # Compute hit rate
    hit <- pnorm(d + qnorm(fa))
    # Create dataframe containing all relevant info
    roc_data <- data.frame(
        dprime  = rep(d, length(fa)),
        hit     = hit,
        fa      = fa
    )
    return(roc_data)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Signal Detection Theory"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("d",
                        "d'",
                        min   = 0,
                        max   = 5,
                        value = 1.5,
                        step  = 0.1),
            sliderInput('criterion',
                        'Criterion',
                        min   = -5,
                        max   = 5,
                        value = 0, 
                        step=0.1),
            sliderInput("noise",
                        "Noise Spread",
                        min   = 1,
                        max   = 10,
                        value = 5,
                        step  = 0.1),
            sliderInput("signal",
                        "Signal Spread",
                        min   = 1,
                        max   = 10,
                        value = 5,
                        step  = 0.1),
            radioButtons('distribution',
                         'Distribution',
                         choices = list('Noise',
                                        'Noise + Signal'))
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("sdtPlot"),
            plotOutput("rocPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    # d' Is the distance between the means of the two distributions
    meanDistr <- reactive({(input$d/2) * 10})
    
    
    
    # Create gaussians
    gaussians <- reactive({
        tibble::tibble(
            grain  = seq(from = -60, to = 60, by = .01),
            noise  = dnorm(grain, mean = -1*meanDistr(), sd = input$noise*2),
            signal = dnorm(grain, mean = meanDistr(), sd = input$signal*2)
        ) %>% 
            tidyr::pivot_longer(cols =  !grain,
                                names_to = 'distribution',
                                values_to = 'value')
    })
    
    
    # Values for plotting
    gaussPeak   <- reactive({max(gaussians()$value)})
    
    selected_distr <- reactive({
        switch(input$distribution,
               'Noise' = 'noise',
               'Noise + Signal' = 'signal')
    })
    left_label <- reactive({
        switch(input$distribution,
               'Noise' = 'Correct Reject',
               'Noise + Signal' = 'Miss')
    })
    
    right_label <- reactive({
        switch(input$distribution,
               'Noise' = 'False Alarm',
               'Noise + Signal' = 'Hit')
    })
    
    
    output$sdtPlot <- renderPlot({
        ggplot(data = gaussians(), aes(x = grain, y = value, fill = distribution)) +
            geom_line(show.legend = FALSE) +
            geom_vline(xintercept = input$criterion*4) +
            geom_ribbon(data = filter(gaussians(), grain < input$criterion*4, distribution == selected_distr()), aes(x=grain, ymax=value, fill = left_label()), ymin=0, alpha = .3) +
            geom_ribbon(data = filter(gaussians(), grain > input$criterion*4, distribution == selected_distr()), aes(x=grain, ymax=value, fill = right_label()), ymin=0, alpha = .3) +
            scale_fill_manual(breaks = c(left_label(), right_label()), values = c( "#ead1dc",'#fff180', 'white', 'white')) +
            annotate('text', x = -1*meanDistr() - 8, y = gaussPeak(), label = 'Noise') +
            annotate('text', x = meanDistr() + 13, y = gaussPeak(), label = 'Signal + Noise') +
            scale_x_continuous(breaks = NULL) + 
            scale_y_continuous(breaks = NULL, limits = c(0, 0.05)) +
            labs(
                x = '',
                y = ''
            ) +
            theme_minimal() +
            theme(legend.text = element_text(size=15),
                  legend.title = element_text(size=18))
    })
    
    # Compute ROC
    triangle_vertex_low <- data.frame(
        x = c(0, 0, 0.5),
        y = c(0.01, 1, 0.51)
    )
    
    triangle_vertex_high <- data.frame(
        x = c(0.5, 0, 1),
        y = c(0.51, 1, 1)
    )
    
    output$rocPlot <- renderPlot({
        selected_d <- input$d
        current_roc <- generate_roc_line(selected_d)
        
        current_hit <- pnorm(input$d / 2 - input$criterion)
        current_fa <- pnorm(-input$d / 2 - input$criterion)
        
        
        ggplot(data=current_roc, aes(x=fa, y=hit)) +
            geom_polygon(data = triangle_vertex_low, 
                         aes(x=x, y=y), 
                         inherit.aes = FALSE,
                         fill = "#fc7b54",
                         alpha = 0.25) +
            geom_polygon(data = triangle_vertex_high, 
                         aes(x=x, y=y), 
                         inherit.aes = FALSE,
                         fill = "#008080",
                         alpha = 0.25) + 
            annotate(geom = "text", x=0.15, y=0.25, label="c>0") +
            annotate(geom = "text", x=0.65, y=0.75, label="c<0") +
            geom_line(color="purple", linewidth=1.5) +
            annotate(geom = "segment", x=0, y=0, xend=1, yend=1) +
            annotate(geom = "point", x = current_fa, y=current_hit, size=3) +
            labs(
                x = "FA RATE",
                y = "HIT RATE",
                title = "d'"
                
            ) +
            theme_minimal() +
            coord_fixed(xlim = c(0,1), 
                        ylim = c(0,1), 
                        expand = TRUE)
        
    })
    
    

    
   
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)

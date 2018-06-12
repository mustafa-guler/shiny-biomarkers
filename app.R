library(shiny)
library(GO.db)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(survival)
library(randomForestSRC)
library(ggRandomForests)
library(gridExtra)

source("loadData.R")
source("goExtraction.R")
source("processData.R")

ui <- fluidPage(
    titlePanel("Cancer Survival BioMarkers"),
    hr(),
    sidebarPanel(
        tags$head(tags$style("#vimp, #variables, #kmplot {height:90vh !important;}")),
        selectInput(
            inputId = "dataset",
            label = "Dataset",
            choices = dir("./data/")
        ),
        textInput(
            inputId = "word",
            label = "GO Filter",
            value = "epigenetic"
        ),
        actionButton(
            inputId = "runButton",
            label = "Run Analysis"
        ),
        br(),
        br(),
        p("Specify a dataset and a GO filter. Genes that have GO annotations containing
          specified GO filter will be considered. If the filter is left blank no gene
          filtering will be done.")
    ),
    mainPanel(
        tabsetPanel(
            type = "tabs",
            tabPanel("Clinical", fluidRow(
                plotOutput("genders"), hr(),
                plotOutput("ages"), hr(),
                plotOutput("vital_status"), hr(),
                plotOutput("race"), hr(),
                plotOutput("ethnicity")
            )),
            tabPanel("Model Variable Importance",
                plotOutput("vimp")
            ),
            tabPanel("Relative Mortality Rates",
                plotOutput("variables")
            ),
            tabPanel("Training and Test Errors",
                htmlOutput("errors")
            ), 
            tabPanel("Kaplan-Meier Plot",
                plotOutput("kmplot")
            ),
            tabPanel("Gene Info",
                tableOutput("gene_info")
            )
        )
    )
)

server <- function(input, output) {
    dataset <- ""
    word <- NA
    
    runClinical <- eventReactive(input$runButton, {
        loadClinical(file.path("data", input$dataset))
    })
    
    runTCGA <- eventReactive(input$runButton, {
        if(dataset != input$dataset) {
            withProgress(
                message = "Loading Data",
                tcga <- loadTCGA(file.path("data", input$dataset))
            )
            dataset <<- input$dataset
        }
        tcga
    })
    
    runGO <- eventReactive(input$runButton, {
        tcga <- runTCGA()
        if(dataset != input$dataset || is.na(word) || word != input$word) {
            withProgress(
                message = "Extracting GO Annotations",
                gores <- goIDs(file.path("data", input$dataset), tcga$gene_id, input$word)
            )
            word <<- input$word
        }
        gores
    })
    
    tr <- integer(0)
    
    runData <- eventReactive(input$runButton, {
        tcga <- runTCGA()
        gores <- runGO()
        clinical <- runClinical()
        filt <- filterData(tcga, gores)
        filt <- addSurvivalData(filt, clinical)
        filt
    })
    
    runModel <- eventReactive(input$runButton, {
        filt <- runData()
        tr <<- trainingData(filt, .7)
        rfsrc(Surv(times, is_dead) ~ ., filt[tr, ], ntree = 1000)
    })
    
    runVIMP <- eventReactive(input$runButton, {
        model <- runModel()
        calcVIMP(model)
    })
    
    output$genders <- renderPlot({
        clinical <- runClinical()
        plotPie(clinical, "gender", "Gender")
    })
    
    output$ages <- renderPlot({
        clinical <- runClinical()
        ggplot(clinical, aes(clinical$age_at_diagnosis)) + geom_histogram() + 
            ggtitle("Ages at Diagnosis") + 
            xlab("Age of Patient at Diagnosis (Years)") #+ 
    })
    
    output$vital_status <- renderPlot({
        clinical <- runClinical()
        plotPie(clinical, "vital_status", "Vital Status")
    })
    
    output$ethnicity <- renderPlot({
        clinical <- runClinical()
        plotPie(clinical, "ethnicity", "Ethnicity")
    })
    
    output$race <- renderPlot({
        clinical <- runClinical()
        plotPie(clinical, "race", "Race")
    })

    output$vimp <- renderPlot({
        plotVIMP(runModel())
    })
    
    output$variables <- renderPlot({
        model <- runModel()
        plotVariables(model, runVIMP())
    })
    
    output$errors <- renderUI({
        model <- runModel()
        filt <- runData()
        trainout <- capture.output(model)
        testout <- capture.output(test(filt, model, tr))
        HTML(
            paste(
                paste(
                    paste("<h2> Training Model </h2>"),
                    paste(trainout, collapse = "<br/>")
                ),
                paste(
                    paste("<h2> Test Model </h2>"),
                    paste(testout, collapse = "<br/>")
                )
            )
            
        )
    })
    
    output$kmplot <- renderPlot({
        dat <- runData()
        vimps <- runVIMP()
        plotKM(dat, vimps)
    })
    
    output$gene_info <- renderTable({
        gores <- runGO()
        # gores %>%
        #     group_by(gores$ensembl_gene_id)
        gores
    })
}

shinyApp(ui = ui, server = server)
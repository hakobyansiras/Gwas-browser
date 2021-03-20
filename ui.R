library(shiny)
library(googleVis)
library(DT)
library(data.table)

shinyUI(
  
  fluidPage(theme = "cerulean.css",
    
    titlePanel (span("GWAS catalog browser", style = "color:blue"), windowTitle = "GWAS"),
    sidebarLayout(
      sidebarPanel(width = 3,
        singleton(tags$head(HTML(
          '
    <script type="text/javascript">
    $(document).ready(function() {
    $("#data_file").attr("disabled", "true").attr("onclick", "return false;");
    
    Shiny.addCustomMessageHandler("download", function(message) {
    $("#data_file").attr("disabled", "true").attr("onclick", "return false;").html(
    "<i class=\\"fa fa-download\\"></i> Download");
    });
    Shiny.addCustomMessageHandler("download_ready", function(message) {
    $("#data_file").removeAttr("disabled").removeAttr("onclick").html(
    "<i class=\\"fa fa-download\\"></i> Download search results");
    });
    })
    </script>
    '
        ))),
        tags$head(tags$script('$(document).keyup(function(event){
            if(event.keyCode == 13){
                      $("#go").click();
            }
        });
        ')),
        
        textInput("a1", label = h5("Please insert searching item"), placeholder = "Gene, SNP, disease"
        ),
        actionButton("go",label =  "Search",icon = icon("search", lib = "glyphicon")),
        
        downloadButton("data_file")
      ),
      
      mainPanel(width = 9,
        div(tabsetPanel(
          tabPanel(title='Initial information',textOutput("text"),
                   textOutput("text2"),tags$head(tags$style("#text2{color: #FBB117;font-size: 20px;}")),
                   tags$head(tags$style("#text{color: #728C00;font-size: 20px;}")),
                   dataTableOutput("table2"),
                   htmlOutput("view"),dataTableOutput("table1")),
          tabPanel(title="Genes and SNPs",textOutput("text1"),
                   div(dataTableOutput("table3"), style = "width: 700px"),br(),
                   div(dataTableOutput("table4"), style = "width: 700px")),
          tabPanel(title="Ontology annotations",textOutput("text3"),dataTableOutput("table5")),
          tabPanel(title="Gwas table",dataTableOutput("table")),
          id="tabs"))
      )
    )
  ))
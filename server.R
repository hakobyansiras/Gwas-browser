library(shiny)
library(DT)
library(gsubfn)
library(data.table)
library(googleVis)
# library(stringr)

gwas <- read.delim("data/gwas_go.tsv",stringsAsFactors = FALSE,encoding = 'UTF-8')
gwas <- data.table(gwas)
gwas$CNV <- NULL

search <- function(y,z) {
  words <- unlist( strsplit( z, "\\ " ) )
  ldt <- apply(y,MARGIN=1, function(x) {
    var <- TRUE
    for (word in words) { 
      var <- grepl(word, x, perl=TRUE,ignore.case = TRUE);
      if ( sum(var) == 0 ) {
        return(FALSE)
      }
    }
    return(TRUE)
  }
  )
  gtac <- y[ldt,]
  return(gtac)
}

sep <- function(x) {
  z <- 0
  x <- gsub(',', '', x)
  matches <- regmatches(x, gregexpr("[[:digit:]]+", x))
  z <- sum(as.numeric(unlist(matches)))
  return(z)
}

sepcikl <- function(x) {
  z <- numeric()
  for (i in 1:(length(x))) {z[i] <- sep(x[i]) }
  return(z)
}

disease <- function(x) {
  dat <- 0
  if(nrow(x) != 0) {
    x <- x[!duplicated(x$DISEASE.TRAIT), ]
    arj <- unlist(x[,9,with = FALSE], recursive = TRUE, use.names = FALSE)
    x$DISEASE.TRAIT <-  apply(x, MARGIN = 1, 
                              function(x)  { names(x) = names(x); paste("<a href=\"", x["LINK"],"\" target=\"_blank\">", x["DISEASE.TRAIT"], "</a>", sep = "") } )
    common_dis <- unlist(x[,8,with = FALSE], recursive = TRUE, use.names = FALSE)
    Participants_sum <- sepcikl(arj)
    dat <- data.table(Disease.Trait = common_dis,Participants_sum)
  }
  return(dat)
}

chr <- function(x) {
  z <- numeric()
  z1 <- character() 
  for(i in 1:23) {z[i] <- nrow(x[CHR_ID == i,]);z1[i] <- i}
  z1[23] <- "X"
  z1[24] <- "Y"
  z[24] <- 0
  dat <- data.table(Chromosome = z1," " = z)
  maxind <- which.max(dat$` `)
  mac <- z1[maxind]
  if(mac == "X") {chrdat <- x[CHR_ID=="23",]}
  else {chrdat <- x[CHR_ID==mac,]}
  snps <- unlist(chrdat[,22,with = FALSE], recursive = TRUE, use.names = FALSE)
  snps <- gsub("x", ",", snps)
  snps <- unlist( strsplit( snps , "\\," ) )
  snps <- snps[! snps %in% "NR"]
  snps <- unique(snps)
  num_snps <- c("(Num_snps = ",length(snps),")")
  if(sum(dat$` `) == 0) {mac <- "Not found";num_snps <- NULL}
  mostchr <- c("Most affected chormosome = ",mac,num_snps)
  return(list(dat,mostchr))
}

numbers <- function(x) {
  study <- unlist(x[,7,with = FALSE], recursive = TRUE, use.names = FALSE)
  study <- unique(study)
  num_study <- length(study)
  snps <- unlist(x[,22,with = FALSE], recursive = TRUE, use.names = FALSE)
  snps <- gsub("x", ",", snps)
  snps <- unlist( strsplit( snps , "\\," ) )
  snps <- snps[! snps %in% "NR"]
  snps <- snps[!is.na(snps)]
  snps <- gsub("(^[[:space:]]+|[[:space:]]+$)", "", snps)
  snps <- unique(snps)
  num_snps <- length(snps)
  genes <- unlist(x[,14,with = FALSE], recursive = TRUE, use.names = FALSE)
  genes <- unlist( strsplit( genes , "\\," ) )
  genes <- genes[! genes %in% "NR"]
  genes <- genes[! genes %in% "intergenic"]
  genes <- genes[!is.na(genes)]
  genes <- gsub("(^[[:space:]]+|[[:space:]]+$)", "", genes)
  genes <- unique(genes)
  num_genes <- length(genes)
  disease <- unlist(x[,8,with = FALSE], recursive = TRUE, use.names = FALSE)
  disease <- unique(disease)
  num_disease <- length(disease)
  return(c("Number of snps = ",num_snps,", studies = ",num_study,", genes = ",num_genes,", disease.traits = ",num_disease))
}

study <- function(x) {
  x <- x[!duplicated(x$STUDY), ]
  x$STUDY <-  apply(x, MARGIN = 1, 
                    function(x)  { names(x) = names(x); paste("<a href=\"", x["LINK"],"\" target=\"_blank\">", x["STUDY"], "</a>", sep = "") } )
  pze <- x[,7,with = FALSE]
  colnames(pze) <- "Studies"
  return(pze)
}

occurrence_gene <- function(x) {
  dat <- NULL
  if(nrow(x) != 0) {
    genes <- unlist(x[,14,with = FALSE], recursive = TRUE, use.names = FALSE)
    genes <- unlist( strsplit( genes , "\\," ) )
    genes <- genes[! genes %in% "NR"]
    genes <- genes[! genes %in% "intergenic"]
    genes <- genes[!is.na(genes)]
    genes <- gsub("(^[[:space:]]+|[[:space:]]+$)", "", genes)
    if(length(genes) != 0) {
      dat <- table(genes)
      dat <- data.table(dat)
      colnames(dat) <- c("Genes","Number_of_occurrences")
      dat$Genes <-  apply(dat, MARGIN = 1, 
                          function(x)  { names(x) = names(dat); paste("<a href=http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",x["Genes"]," target=\"_blank\">", x["Genes"], "</a>", sep = "") } )
      dat <- dat[order(dat$Number_of_occurrences,decreasing = TRUE),] 
    }
  }
  return(dat)
}

occurrence_snp <- function(x) {
  dat <- NULL
  if(nrow(x) != 0) {
    snps <- unlist(x[,22,with = FALSE], recursive = TRUE, use.names = FALSE)
    snps <- gsub("x", ",", snps)
    snps <- unlist( strsplit( snps , "\\," ) )
    snps <- snps[! snps %in% "NR"]
    snps <- snps[!is.na(snps)]
    snps <- gsub("(^[[:space:]]+|[[:space:]]+$)", "", snps)
    if(length(snps) != 0) {
      dat <- table(snps)
      dat <- data.table(dat)
      colnames(dat) <- c("Snps","Number_of_occurrences")
      dat$Snps <-  apply(dat, MARGIN = 1, 
                         function(x)  
                         { names(x) = names(dat); 
                         paste("<a href=http://asia.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=7:87498388-87499388;v=", x["Snps"],";vdb=variation;vf=115741667"," target=\"_blank\">", x["Snps"], "</a>", sep = "") } )
      dat <- dat[order(dat$Number_of_occurrences,decreasing = TRUE),]
    }
  }
  return(dat)
}

ontology <- function(x) {
  dat <- NULL
  if(nrow(x) != 0) {
    x <- x[!duplicated(x$DISEASE.TRAIT), ]
    ont <- apply(x, MARGIN = 1, function(x)  {
      names(x) = names(x);
      paste("<a href=\"", unlist( strsplit( x["MAPPED_TRAIT_URI"] , "\\," ) ),"\" target=\"_blank\">", 
            unlist( strsplit( x["MAPPED_TRAIT"] , "\\," ) ), "</a>", sep = "", collapse = ",")} )
    dat <- data.table(DISEASE = x$DISEASE.TRAIT,MAPPED_TRAIT = ont)
  }
  return(dat)
}

gwastable <- function (x) {
  gwastab <- NULL
  if(nrow(x) != 0) {
    gwastab <- x
    gwastab$PUBMEDID <- apply(gwastab, MARGIN = 1, function(x)  { names(x) = names(gwastab); paste("<a href=\"", x["LINK"],"\" target=\"_blank\">", x["PUBMEDID"], "</a>", sep = "") } )
    gwastab$FIRST.AUTHOR <- apply(gwastab, MARGIN = 1, function(x)  { names(x) = names(gwastab); paste("<a href=\"", x["LINK"],"\" target=\"_blank\">", x["FIRST.AUTHOR"], "</a>", sep = "") } )
    gwastab$MAPPED_TRAIT <- apply(gwastab, MARGIN = 1, function(x)  {
      names(x) = names(gwastab);
      paste("<a href=\"", unlist( strsplit( x["MAPPED_TRAIT_URI"] , "\\," ) ),"\" target=\"_blank\">", 
            unlist( strsplit( x["MAPPED_TRAIT"] , "\\," ) ), "</a>", sep = "", collapse = ",")} )
    gwastab$MAPPED_TRAIT_URI <- NULL
    gwastab$LINK <- NULL
  }
  return(gwastab)
}

shinyServer(function(input, output, session) {
  
  v <- reactiveValues(data = NULL,study = NULL,plot = NULL,num = NULL,most_chr = NULL,
                      inf = NULL,logicdown = FALSE,num_occurrence_g = NULL,num_occurrence_s = NULL,ontology_tab = NULL,gwas_tab = NULL)
  
  observeEvent(input$go, {
    reply <- search(gwas,input$a1)
    if(nrow(reply) != 0) {
      v$data <- disease(reply);
      v$study <- study(reply);
      v$plot <- chr(reply)[[1]];
      v$most_chr <- chr(reply)[[2]];
      v$num <- numbers(reply);
      v$num_occurrence_g <- occurrence_gene(reply);
      v$num_occurrence_s <- occurrence_snp(reply);
      v$ontology_tab <- ontology(reply);
      v$gwas_tab <- gwastable(reply);
      v$logicdown <- TRUE;
      v$inf <- NULL
    } 
    else {
      v$data <- NULL;
      v$study <- NULL;
      v$plot <- NULL;
      v$num <- NULL;
      v$most_chr <- NULL;
      v$num_occurrence_g <- NULL;
      v$num_occurrence_s <- NULL;
      v$ontology_tab <- NULL;
      v$gwas_tab <- NULL;
      v$logicdown <- FALSE;
      v$inf <- "Not found"
    }
  })
  
  observe({
    if (v$logicdown) {
      session$sendCustomMessage("download_ready",TRUE)
    }
    else {
      session$sendCustomMessage("download",TRUE)
    }
  })
  
  myout <- reactive({
    dat <- search(gwas,input$a1)
    return(dat)
  })
  
  output$data_file <- downloadHandler(
    filename = function() { 
      paste('Gwas', '.csv', sep='')
    },
    content = function(file) {
      write.csv(myout(), file)
    }
  )
  
  output$text <- renderText ({
    paste(v$num,v$inf)
  })
  
  output$text2 <- renderText ({
    paste(v$most_chr)
  })
  
  output$text1 <- renderText ({
    paste(v$inf)
  })
  
  output$text3 <- renderText({
    paste(v$inf)
  })
  
  output$table5 <- DT::renderDataTable({
    DT::datatable(
      selection = 'none',
      v$ontology_tab,escape=FALSE,
      filter = 'top',
      class = "cell-border",
      options = list(dom = 'rtp',pageLength = 10,
                     searchHighlight = TRUE,
                     initComplete = JS("function(settings, json) {",
                                       "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                       "}"),
                     scrollX = TRUE
      ),
      style = 'bootstrap'
    )
  })
  
  output$table3 <- DT::renderDataTable({
    DT::datatable(
      selection = 'none',
      v$num_occurrence_g,escape=FALSE,
      filter = 'top',
      class = "cell-border",
      options = list(dom = 'rtp',pageLength = 10,
                     searchHighlight = TRUE,
                     initComplete = JS("function(settings, json) {",
                                       "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                       "}")
      ),
      style = 'bootstrap'
    )
  })
  
  output$table4 <- DT::renderDataTable({
    DT::datatable(
      selection = 'none',
      v$num_occurrence_s,escape=FALSE,
      filter = 'top',
      class = "cell-border",
      options = list(dom = 'rtp',pageLength = 10,
                     searchHighlight = TRUE,
                     initComplete = JS("function(settings, json) {",
                                       "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                       "}")
      ),
      style = 'bootstrap'
    )
  })
  
  output$table2 <- DT::renderDataTable({
    DT::datatable(
      selection = 'none',
      v$study,escape=FALSE,
      class = "cell-border",
      options = list(dom = 'frtp',pageLength = 10,
                     searchHighlight = TRUE,
                     initComplete = JS("function(settings, json) {",
                                       "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                       "}"),
                     scrollX = TRUE
      ),
      style = 'bootstrap'
    )
  })
  
  output$view <- renderGvis({
    validate(need(v$plot, ''))
    gvisColumnChart(v$plot,
                    options=list(backgroundColor = '#fcfcfc',
                                 title = "Affected chromosomes",
                                 titleTextStyle="{color:'blue',
                                 fontSize:18}",
                                 isStacked = "true",
                                 vAxes="[{title:'Number',
                                 textPosition: 'out'}]",
                                 hAxes="[{title:'Chromosome',
                                 textPosition: 'out',
                                 textStyle: {
                                 fontSize: 10}}]",
                                 width=850, height=450)
    )
  })
  
  output$table1 <- DT::renderDataTable({
    DT::datatable(
      selection = 'none',
      v$data,escape=FALSE,
      class = "cell-border",
      options = list(dom = 'frtp',pageLength = 10,
                     searchHighlight = TRUE,
                     initComplete = JS("function(settings, json) {",
                                       "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                       "}"),
                     scrollX = TRUE
      ),
      style = 'bootstrap'
    )
  })
  
  output$table <- DT::renderDataTable({
    DT::datatable(
      selection = 'none',
      v$gwas_tab,escape=FALSE,
      filter = 'bottom',
      class = "cell-border",
      extensions = 'FixedColumns',
      options = list(dom = 'lrtip',
                     searchHighlight = TRUE,
                     initComplete = JS("function(settings, json) {",
                                       "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                       "}"),
                     scrollY = 600,scrollX = TRUE
      ),
      style = 'bootstrap'
    )
  })
  
  })


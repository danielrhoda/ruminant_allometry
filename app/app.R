
# Interactive Shiny App

#
library(shiny)
library(geomorph)
library(Morpho)
library(adephylo)
library(pracma)
library(shinydashboard)
library(phytools)

#load("data/ruminant_app_data.Rdata", envir = .GlobalEnv)
#readRDS("data/ruminant_app_data.rds")

source("app_functions.R")

ordinations <- readRDS(file = "ordinations.rds")
R.df <- readRDS(file = "Rdf.rds")
R2.df <- readRDS("Rdf2.rds")
family.colors <- readRDS("familycolors.rds")
ruminant.taxonomy <- readRDS("ruminanttaxonomy.rds")
full.inventory1 <- readRDS("inventory.rds")
pts.pch <- readRDS("pch.rds")
pts.col <- readRDS("col.rds")
col.pal <- readRDS("colpal.rds")
mR <- readRDS("mR.rds")
R.size <- R.df$csize
R <- R.df$coords





 ui <- dashboardPage(
   dashboardHeader(title = "Ruminant Morphospaces", titleWidth = 'auto'),
   dashboardSidebar(
     selectInput(inputId = 'select',
                 label = "Ordination",
                 choices = c("PCA", "pPCA", "PACA", "PCA (allometry free)", "pPCA (allometry free)", "PACA (allometry free)"),
                 selected = "PCA"),
     fluidRow(
       column(6,
         numericInput(inputId = "xaxis",
                  label = "X Axis PC",
                  value =  1, min = 1, max = 68)),
       column(6,
         numericInput(inputId = 'yaxis',
                  label = "Y Axis PC",
                  value = 2, min = 1, max = 68)))

   ),
   dashboardBody(
   fluidRow(
     column(6,
       plotOutput("plot1",
              width = 900,
              height = 900,
              click = "plot_click")),
     column(6,
       plotOutput("TPS",
                  width = 900,
                  height = 900))
   ),

  fluidRow(
    column(3,
                verbatimTextOutput("info"),
                   tags$head(tags$style(HTML("
                            #info {
                              font-size: 19px;
                            }
                            ")))
    ),
    column(3,
#           div(style = "height:5px"),
           plotOutput("size")),
    column(3,
#           div(style = "height:5px"),
           plotOutput("face")),
    column(3,
#           div(style = "height:5px"),
           plotOutput("diet")),
#     column(2,
# #           div(style = "height:5px"),
#            plotOutput("scree"))
     )
   )
 )



# THINGS TO DO:
# HIGHLIGHT THE BRANCHES IN PHYLOMORPHOSPACE LEADING UP TO THE HILIGHTED TAXON
# PUT A PHYLOGENY TO THE RIGHT, WITH SIZE/FACELENGTH DOTMAPS
# HIGHLIGHT THE BRANCHES LEADING UP TO THE TAXON ON THE PHYLOGENY
# ADD TPS SPLINES TO THE RIGHT
# ADD QUICK (12 FRAME, OR 1 FRAME PER BRANCH) ANIMATIONS OF ANCESTOR -> DESCENDANT TRAJECTORIES THROUGH MORPHOSPACE


 server <- function(input, output){


  ordination.selection <- reactive({ordinations[[input$select]]})
  PC12 <- reactive({ordination.selection()$x[,c(input$xaxis, input$yaxis)]})

  XY <- reactiveValues(xy=data.frame(x=c(0.032,0),y=c(-0.039,0)))
  observeEvent(input$plot_click, {
        XY$xy[1,] <- c(input$plot_click$x, input$plot_click$y)
      })


  taxon.dist <- reactive({apply(PC12(), 1, function(x)euc.dist(x,XY$xy[1,]))})
  taxon.id <- reactive({as.character(names(sort(taxon.dist()))[1])})

  taxon.pc <- reactive({PC12()[which(rownames(PC12())==taxon.id()),]})



  descendant.node <- reactive({which(R.df$tree$tip.label==taxon.id())})
  MRCA <- reactive({nodepath(phy=R.df$tree,from = 131, to = descendant.node())})

  anc.x <- reactive({cbind(fastAnc(tree = R.df$tree, x = PC12()[,1]), fastAnc(tree = R.df$tree, x = PC12()[,2]))})
  anc <- reactive({anc.x()})
  ind <- reactive({match(R.df$tree$tip.label, rownames(PC12()))})
  pts <- reactive({PC12()[ind(), ]})
  z1 <- reactive({rbind(pts(),anc())[MRCA(),]})



  family <- reactive({full.inventory1[taxon.id(), "family"]})
  subfamily <- reactive({full.inventory1[taxon.id(), "subfamily"]})
  tribe <- reactive({full.inventory1[taxon.id(), "tribe"]})
  diet <- reactive({A <- full.inventory1[taxon.id(), "diet_category"]
                      if(A=="BR"){A<-"Browsing"}
                      if(A=="IM"){A<-"Intermediate"}
                      if(A=="GR"){A<-"Grazing"}else{A<-"No data"}})
  size <- reactive({full.inventory1[taxon.id(), "max_BM"]})
  hornlength <- reactive({full.inventory1[taxon.id(), "HL_male_cm"]})
  commonname <- reactive({full.inventory1[taxon.id(), "common_name"]})






   output$plot1 <- renderPlot({

     XLIM <- range(PC12()[,1])*1.3
     YLIM <- range(PC12()[,2])*1.3
     if(diff(XLIM) > diff(YLIM)){YLIM <- XLIM}else{XLIM <- YLIM}


          par(fig=c(0,0.9,0.793,1),new=F, bg = "#ECF0F5")
          plot(NA,xlim= XLIM, ylim=c(0,0.68),xlab=NA,ylab=NA,
               frame = F, axes = F)
          W <- table(ruminant.taxonomy[,1])
          uaf <- c("Cervidae" ,"Bovidae", "Tragulidae" , "Moschidae")
          for(i in 1:4){
          famid <- which(ruminant.taxonomy[,1] == uaf[i])
          dens <- density(PC12()[famid,1])
          dens$y <- dens$y/max(dens$y)
          dens$y <- dens$y * (W[uaf[i]]/sum(W))
          polygon(dens,col = addTrans(family.colors[i],255*0.1), border = family.colors[i], lwd =2)
          arrows(x0=taxon.pc()[1],
                 y0=0,
                 x1=taxon.pc()[1],
                 y1=5,
                 col = "red", lwd = 2, code = 1)
          }


          par(fig=c(0.805,1,0,0.9),new=T)
          plot(NA,ylim=YLIM,xlim=c(0,0.68),xlab=NA,ylab=NA,
               frame = F, axes = F)
          #axis(3, labels = FALSE)
          for(i in 1:4){
          famid <- which(ruminant.taxonomy[,1] == uaf[i])
          dens <- density(PC12()[famid,2])
          dens$y <- dens$y/max(dens$y)
          dens$y <- dens$y * (W[uaf[i]]/sum(W))
          dens2 <- dens
          dens$x <- dens2$y
          dens$y <- dens2$x
          polygon(dens,col = addTrans(family.colors[i],255*0.1), border = family.colors[i], lwd = 2)
          arrows(y0=taxon.pc()[2],
                 x0=0,
                 y1=taxon.pc()[2],
                 x1=5,
                 col = "red", lwd = 2, code = 1)

          }



   cervid.scores <- PC12()[which(ruminant.taxonomy[,1] == "Cervidae"),]
   cervid.chull <- chull(cervid.scores)
   cervid.chull <- c(cervid.chull,cervid.chull[1])

   bovid.scores <- PC12()[which(ruminant.taxonomy[,1] == "Bovidae"),]
   bovid.chull <- chull(bovid.scores)
   bovid.chull <- c(bovid.chull,bovid.chull[1])

   moschid.scores <- PC12()[which(ruminant.taxonomy[,1] == "Moschidae"),]
   moschid.chull <- chull(moschid.scores)
   moschid.chull <- c(moschid.chull,moschid.chull[1])

   tragulid.scores <- PC12()[which(ruminant.taxonomy[,1] == "Tragulidae"),]
   tragulid.chull <- chull(tragulid.scores)
   tragulid.chull <- c(tragulid.chull,tragulid.chull[1])


  if(input$select == "PCA" | input$select == "PCA (allometry free)"){pclab <- "PC"}
  if(input$select == "pPCA" | input$select == "pPCA (allometry free)"){pclab <- "pPC"}
  if(input$select == "PACA" | input$select == "PACA (allometry free)"){pclab <- "PAC"}


            par(fig = c(0,0.9,0,0.9),
                   mar=c(5,5,2,2),
                new = T)
             plot(NA,
                  xlim=XLIM,
                  ylim=YLIM,
                  asp=1,
                  cex.lab=1.5,
                 xlab = paste0(pclab,input$xaxis,": ",
                               toString(round(ordination.selection()$d[input$xaxis]/sum(ordination.selection()$d),digits=3)*100),"% of shape variation"),
                 ylab = paste0(pclab,input$yaxis,": ",
                               toString(round(ordination.selection()$d[input$yaxis]/sum(ordination.selection()$d),digits=3)*100),"% of shape variation"))
          grid(lty=1,col = "lightgray")

            lines(cervid.scores[cervid.chull,1:2],col=family.colors[1],lwd=3)
            lines(bovid.scores[bovid.chull,1:2],col=family.colors[2],lwd=3)
            lines(moschid.scores[moschid.chull,1:2],col=family.colors[4],lwd=3)
            lines(tragulid.scores[tragulid.chull,1:2],col=family.colors[3],lwd=3)

            addtree2(PC12(),tree = R.df$tree)
            points(PC12()[,1],PC12()[,2],pch=pts.pch,bg = pts.col,lwd = 1, cex = rescale.numeric(R.size,c(1,4)))
            points(z1(), type = "l", col = "black", lwd = 4, lty = 1)
            points(rbind(ordination.selection()$anc.x[1,c(input$xaxis, input$yaxis)],
                         ordination.selection()$anc.x[1,c(input$xaxis, input$yaxis)]), pch = 23, bg = 'goldenrod1', cex = 4, lwd = 1)

               points(XY$xy[1,],
                    cex = 2, pch = 8, col = "red", lwd = 1
                    )
              points(taxon.pc()[1], taxon.pc()[2],
                    cex = 4, pch = 19, col = "red",
                    )


            legend("topleft", c("Tragulidae","Cervidae","Moschidae","Bovidae","Common Ancestor", "Clicked Point", "Highlighted Taxon"),
                   pt.lwd = 1,
                   pt.cex = 2.375,
                   col = c(rep("black", 5), "red", "red"),
                   pt.bg = c(family.colors[c(3,1,4,2)],"goldenrod1", "red", "red"),
                   pch = c(24,21,25,22,23, 8, 19),
                   title.adj = -1,
                   cex = 1.5,
                   bg = "transparent")
            box()


   })

   output$TPS <- renderPlot({
     par(bg = "#ECF0F5")
     tt <- gsub("_"," ",taxon.id()) ; ST <- gsub("_"," ",commonname())
     ruminant.TPS.app(r = R.df$coords[,,taxon.id()],
                      title =  tt,
                      subtitle =  ST,
                      meanshape = mR)
   })



   output$info <- renderText({
      if(input$select == "PCA" | input$select == "PCA (allometry free)"){pclab <- "PC"}
      if(input$select == "pPCA" | input$select == "pPCA (allometry free)"){pclab <- "pPC"}
      if(input$select == "PACA" | input$select == "PACA (allometry free)"){pclab <- "PAC"}

      paste0(
        "Species scientific name: ", gsub("_"," ",taxon.id()), "\n",
        "Species common name: ", gsub("_"," ",commonname()), "\n",
        "Family: ", family(), "\n",
        "Subfamily: ", subfamily(), "\n",
        "Tribe: ", tribe(), "\n",
        "Max body size: ", round(size(),2), "kg\n",
        "Diet: ", diet(), "\n",
        "Max horn length: ", hornlength(),"cm"
      )

      }
   )


   output$size <- renderPlot({
    par(bg = "#ECF0F5", mar = c(5,.1,0.1,.1), fig = c(0,1,0.41,1))
    dens <- density(R.df$csize)
    plot(NA, xlim = c(min(R.df$csize)*0.9, max(R.df$csize)*1.1),
         ylim = c(0,max(dens$y)*1.1),
         axes = F,
         frame = T,
         ylab = NA,
         xlab = "log(centroid size)",
         cex.lab = 2,
         col.lab = family.colors[5],
         font.main = 1)
    axis(1)
     polygon(dens, col = addTrans(family.colors[1],255*0.1), border = family.colors[1], lwd = 2)
     box()
     arrows(x0 = R.df$csize[taxon.id()], x1 = R.df$csize[taxon.id()],
            y0 = 0, y1 = max(dens$y)*0.5, code = 1)
   })

  output$face <- renderPlot({
    par(bg = "#ECF0F5", mar = c(5,.1,0.1,.1), fig = c(0,1,0.41,1))
    dens <- density(R.df$facelength)
    plot(NA, xlim = c(min(R.df$facelength)*0.9, max(R.df$facelength)*1.1),
         ylim = c(0,max(dens$y)*1.1),
         axes = F,
         frame = T,
         ylab = NA,
         xlab = "face length / skull length",
         cex.lab = 2,
         col.lab = family.colors[5],
         font.main = 1)
    axis(1)
     polygon(dens, col = addTrans(family.colors[1],255*0.1), border = family.colors[1], lwd = 2)
     box()
     arrows(x0 = R.df$facelength[taxon.id()], x1 = R.df$facelength[taxon.id()],
            y0 = 0, y1 = max(dens$y)*0.5, code = 1)
   })

   output$diet <- renderPlot({
    par(bg = "#ECF0F5", mar = c(5,.1,.1,.1), fig = c(0,1,0.41,1))
    dens <- density(R2.df$percentgrass)
    plot(NA, xlim = c(min(R2.df$percentgrass)*0.9, max(R2.df$percentgrass)*1.1),
         ylim = c(0,max(dens$y)*1.1),
         axes = F,
         frame = T,
         ylab = NA,
         xlab = "% Grass in Diet", cex.lab = 2,
         col.lab = family.colors[5],
         font.main = 1)
    axis(1)
     polygon(dens, col = addTrans(family.colors[1],255*0.1), border = family.colors[1], lwd = 2)
     box()
     grass1 <- setNames(R2.df$percentgrass, names(R2.df$csize))
     grass <- grass1[taxon.id()]
     arrows(x0 = grass, x1 = grass,
            y0 = 0, y1 = max(dens$y)*0.5, code = 1)
   })

  # output$scree <- renderPlot({
  #   par(bg = "#ECF0F5", mar = c(4,.1,4,.1), fig = c(0,1,0.5,1))
  #
  #   pca.d <- vapply(ordination.selection()$d, FUN.VALUE = 1, FUN = function(x) x/sum(R.pca$d))
  #   plot(pca.d[1:20], type = "h", lwd = 4, col = "gray",
  #        xlab="Eigenvalue", ylab = NA,
  #        axes = F,
  #        frame = T, cex.lab = 2,
  #        main = "Chosen PC Axis",
  #        cex.main = 3,
  #        col.main = family.colors[5],
  #        font.main = 1)
  #   axis(1)
  #
  #
  #  })

 }


 shinyApp(ui, server)







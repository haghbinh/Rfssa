# create the shiny application user interface
ui.fssa <- fluidPage(tags$head(tags$style(HTML("body { max-width: 1250px !important; }"))),
                     titlePanel("FSSA Illustration"),
                     sidebarLayout(
                       sidebarPanel(width=3, tags$head(tags$style(type="text/css", ".well { max-width: 300px; }")),
                                    radioButtons("bs.fr","Choose Basis:", choices = c("B-spline", "Fourier"), selected = "B-spline", inline=TRUE),
                                    uiOutput("xdeg", width="250px"), uiOutput("xdf", width="250px"),
                                    tags$hr(style="border-color: red;", width="150px"),
                                    column(6,uiOutput("g")), column(6,uiOutput("sg")), column(6,uiOutput("d")), column(6,uiOutput("dmd")),
                                    sliderInput("mssaL", HTML("Win.L. (MSSA):"), min = 1, max = 50, value = 50, step = 1, width="210px"),
                                    sliderInput("fssaL", HTML("Win.L. (FSSA):"), min = 1, max = 50, value = 20, step = 1, width="210px"),
                                    column(6,uiOutput("run.fpca")), column(6,uiOutput("run.ssa"))),
                       mainPanel(width=9, tags$style(type="text/css", ".shiny-output-error { visibility: hidden; }",".shiny-output-error:before { visibility: hidden; },"),#".nav-tabs {font-size: 10px}"),
                                 tabsetPanel(id="Panel", type = "tabs",
                                             tabPanel(title="Data", value="Data",
                                                      column(12,uiOutput("ts.selected", align = "center"), style="color:red;"),
                                                      fluidRow(column(4,radioButtons('f.choice', 'Choose from:', c("Server" = "server", "Upload" = "upload", "Simulate" = "sim"), selected= "sim", inline = TRUE, width="250px")),
                                                               column(4,uiOutput("s.choice", width="250px")), column(2,uiOutput('noise.t', width="125px")), column(2,uiOutput('noise.p', width="125px")),
                                                               column(4,uiOutput("file")), column(4,uiOutput("sep"),uiOutput("header"))),
                                                      column(4,uiOutput("model")), column(4,uiOutput("t.len")), column(2,uiOutput("a.f")), column(2,uiOutput("n.sd")),
                                                      column(8,plotOutput("data.plot", height = 600, width = 600)), column(4,tableOutput('data'))),
                                             tabPanel("Basis Functions",
                                                      column(8,plotOutput("basis.desc", height = 600, width = 600)), column(4,uiOutput("basis.n", width="300px"))),
                                             tabPanel("Data Description (SSA Summary)",
                                                      column(4,uiOutput("desc", width="250px")), column(4,uiOutput("as.choice", width="400px"), uiOutput("run.fda.gcv", width="200px"), uiOutput("rec.type", width="300px")), column(2,uiOutput("freq")), column(2,uiOutput("sts.choice")),
                                                      fluidRow(column(8,conditionalPanel(condition="output.flag_plot", plotOutput("res.plot", height = 600, width = 600)),
                                                                      conditionalPanel(condition="output.flag_plotly", plotlyOutput("res.ly", height = 600, width = 600))),
                                                               column(4,uiOutput("s.plot"), fluidRow(column(8,uiOutput("b.indx")), column(4,uiOutput("s.CI"))), column(12, uiOutput("comp.obs"), verbatimTextOutput("RMSEs"))))),
                                             tabPanel("Manual", includeMarkdown(system.file("shiny/rmd", "report.Rmd", package = "Rfssa")))))
                     )
)

#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom fda pca.fd eval.fd
#' @importFrom graphics abline axis par plot points polygon
#' @importFrom stats fft integrate rnorm sd ts.plot
#' @importFrom utils head read.table
#' @import Rssa
# Define server logic required to draw a histogram

server.fssa <- function(input, output, clientData, session) {

  load(system.file("shiny/data", "servshiny.rda", package = "Rfssa"));
  iTs <- reactiveVal(list()); iTrs <- reactiveVal(list()); itmp <- reactiveVal(0)
  df <- 100; vf <- 20; T <- 100;
  output$flag_plotly <- reactive(input$desc%in%c("fssa.reconst","ssa.reconst") && input$rec.type%in%c("heatmap","line","3Dline","3Dsurface"));
  output$flag_plot <- reactive(!(input$desc%in%c("fssa.reconst","ssa.reconst") && input$rec.type%in%c("heatmap","line","3Dline","3Dsurface")))
  outputOptions(output, "flag_plotly", suspendWhenHidden = FALSE);
  outputOptions(output, "flag_plot", suspendWhenHidden = FALSE)


  rfar <- function(N,norm,psi,Eps,basis){
    # Create Corresponding matrix of an integral (kernel) operator
    # Get an kernel function k(s,t) corresponds to an integral operator and
    # a basis system incuding d basis functions, return corresponding d*d matrix of the operator with respect
    # to basis system.
    OpsMat <- function(kernel,basis){
      u <- seq(0,1,by=0.01)
      n <- length(u)
      K_mat <- outer(u,u,FUN = kernel)
      K_t <- smooth.basis(u,K_mat,basis)$fd # the kernel function convert
      # to fd object w.r.t the first argument. K_t is an objcect of n fd.
      A <- inprod(K_t,basis) #An n*d matrix.
      K <- smooth.basis(u,A,basis)$fd # d fd object.
      B <- inprod(K,basis) # An d*d matrix.
      return(B) # return to OpsMat.
    }
    Psi_mat0 <- OpsMat(psi,basis)
    Gram <- inprod(basis,basis)
    Psi_mat <- solve(Gram)%*%Psi_mat0
    E <- Eps$coefs
    X <- E;
    for(i in 2:N) X[,i] <- Psi_mat%*%X[,i-1]+E[,i]
    X_fd <- fd(X,basis)
    return(X_fd) # return to rfar function.
  }
  gamma0 <- function(norm){
    f <- function(x){
      g <- function(y) psi0(x,y)^2
      return(integrate(g,0,1)$value) #return into f.
    }
    f <- Vectorize(f)
    A <- integrate(f,0,1)$value
    return(norm/A) #return into gamma.
  }
  psi0 <- function(x,y) 2-(2*x-1)^2-(2*y-1)^2

  fpca_proj <- function(i,U){
    harm <- U$harmonics[i]
    scores <- U$scores[,i]
    m <- nrow(harm$coefs); n <- length(scores)
    coef <- matrix(NA,nrow=m,ncol=n)
    for(i0 in 1:m) for(j in 1:n) coef[i0,j] <- harm$coefs[i0]*scores[j]
    pc <- fd(coef,harm$basis)
    return(pc)
  }
  fpca_rec <- function(d1,d2,U){
    s <- fpca_proj(d1,U);
    if (d2>d1) for(i0 in (d1+1):d2) s <- s + fpca_proj(i0,U)
    return(s)
  }

  Tr <- function(tau, t){
    ifelse(input$model%in%c("f1","f12"),1,0) * exp(tau^2) * cos(2 * pi * t * input$a.f) +
      ifelse(input$model%in%c("f12","f2"),-1,0) * cos(2*pi*tau/0.5) * sin(2 * pi * t * input$a.f)
  }

  simulate <- function() {
    if (is.null(input$a.f)) return()
    tau <- seq(0, 1, length = T); t <- 1:input$t.len;
    Trs <- outer(tau, t, FUN = Tr)
    set.seed(T*input$t.len*input$a.f*input$n.sd); #set.seed(10)
    noise <- Z <- matrix(rnorm(input$t.len*T, 0, input$n.sd), nrow=T);
    if (input$noise.t=="ar1") {
      if (input$noise.p) {
        Z[1, ] <- 0; A <- diag(1,T)
        if (T>2) for (i in 1:(T-2))  diag(A[-(1:i),]) <- input$noise.p^i;
        if (T>1) A[T,1] <- input$noise.p^(T-1)
        noise <- A%*%Z
      }
    }
    if (input$noise.t=="swn") {
      Z <- matrix(rnorm(input$t.len * input$xdf, 0, input$n.sd), ncol=input$t.len);
      basis.Z <- fda::create.bspline.basis(c(0, 1), input$xdf); tau <- seq(0, 1, length = T)
      basis.noise <- fda::fd(Z, basis.Z); noise <- eval.fd(tau, basis.noise)
    }
    if (input$noise.t=="far1") {
      k0 <- gamma0(input$noise.p)
      psi <- function(x,y) k0*psi0(x,y)
      Z[1, ] <- 0; noise <- apply(Z, 2, cumsum)
      if (input$bs.fr=="B-spline") basis.Z <- fda::create.bspline.basis(c(0, 1), input$xdf)
      else basis.Z <- fda::create.fourier.basis(c(0, 1), input$xdf); tau <- seq(0, 1, length = T)
      Eps <- smooth.basis(tau, noise, basis.Z)$fd
      basis.noise <- rfar(input$t.len,input$noise.p,psi,Eps,basis.Z); noise <- eval.fd(tau, basis.noise)
    }
    #Trs <- scale(Trs,scale=F);
    return(list(Trs=Trs, noise=noise))
  }

  observeEvent(input$dmd, {updateTabsetPanel(session, "Panel", selected = "Data"); updateCheckboxGroupInput(session, "s.plot", selected="")})

  output$xdeg <- renderUI({
    if (input$bs.fr=="Fourier") return()
    sliderInput("xdeg", HTML("Degree of B-spline Basis:"), min = 0, max = 5, value = 3, step = 1, width="250px")
  })

  output$xdf <- renderUI({
    sliderInput("xdf", paste("Deg. of freedom of", input$bs.fr, "Basis:"), min = ifelse(input$bs.fr=="B-spline", input$xdeg+1, 3), max = df, value = vf, step = ifelse(input$bs.fr=="B-spline", 1, 2), width="250px")
  })

  output$g <- renderUI({
    textInput("g", "Groups", value = "1:2")
  })

  output$sg <- renderUI({
    m <- length(eval(parse(text=paste0("list(",input$g,")"))))
    sliderInput("sg", "Select G:", min = 1, max = m, value = c(1,m), step = 1)
  })

  output$dmd  <- renderUI({
    checkboxGroupInput('dmd', 'Functions', choices=c('Demean'='dmd','Dbl Range'='dbl'))#,'Heat plot'='hp'))
  })

  output$d <- renderUI({
    sliderInput("d", "d", min = 1, max = min(input$fssaL,input$mssaL), value = c(1,2), step = 1)
  })

  output$run.ssa <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    actionButton('run.ssa', paste('run (F)SSA'))
  })

  output$run.fpca <- renderUI({return()
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    actionButton('run.fpca', paste('run (D)FPCA'))
  })

  run_ssa <- eventReactive(input$run.ssa, {
    withProgress(message = 'SSA.FSSA: Running', value = 0, {
      if(input$bs.fr=="B-spline") bas.fssa <- fda::create.bspline.basis(c(0, 1), nbasis=input$xdf, norder=input$xdeg+1)
      else bas.fssa <- fda::create.fourier.basis(c(0, 1), nbasis=input$xdf);
      tau <- seq(0, 1, length = nrow(iTs()));
      Uf <- fssa(fts(smooth.basis(tau, iTs(), bas.fssa)$fd), input$fssaL);
      Us <- ssa(t(iTs()), input$mssaL, kind = "mssa");
      return(list(Uf=Uf, Us=Us, tau=tau, bas.fssa=bas.fssa))
    })
  })

  run_fpca <- function(){#eventReactive(input$run.fpca, {
    withProgress(message = '(D)FPCA: Running', value = 0, {
      if(input$bs.fr=="B-spline") bas.fssa <- fda::create.bspline.basis(c(0, 1), nbasis=input$xdf, norder=input$xdeg+1)
      else bas.fssa <- fda::create.fourier.basis(c(0, 1), nbasis=input$xdf);
      tau <- seq(0, 1, length = nrow(iTs())); Y <- smooth.basis(tau, iTs(), bas.fssa)$fd
      f.pca <- pca.fd(Y, nharm = min(input$d[2],input$xdf), centerfns = FALSE); fpca.rec <- fpca_rec(min(input$d[1],input$xdf),min(input$d[2],input$xdf), f.pca); fpca.re <- eval.fd(tau, fpca.rec)
      #d.fpca <- fts.dpca(Y, Ndpc = input$d[2]); dfpca.re <- eval.fd(tau,d.fpca$Xhat)
      return(list(fpca=fpca.re, dfpca=fpca.re))
    })
  }#)

  output$s.choice <- renderUI({
    if (input$f.choice!="server") return();
    s.choices <- 1:length(Xs); names(s.choices) <- names(Xs);
    selectInput("s.choice","Select a file from server: ", choices = s.choices, width="250px");
  })

  output$noise.t <- renderUI({
    if (input$f.choice!="sim") return();
    selectInput("noise.t","Type of noice: ", choices = c("AR(1)"="ar1", "FAR(1)"="far1", "Smooth WN"="swn"), width="125px");
  })

  output$noise.p <- renderUI({
    if (input$f.choice!="sim") return(); if (!is.null(input$noise.t)) if (input$noise.t=="swn") return();
    sliderInput("noise.p", "AR Parameter:", min = 0, max = 1, value = 0, step = 0.01, width="125px")
  })

  output$model <- renderUI({
    if (input$f.choice!="sim") return();  choices <- c("f(tau)" = "f1", "g(tau)" = "f2", "f(tau) + g(tau)" = "f12")
    radioButtons('model', 'Model:', choices=choices, selected="f1", inline = TRUE, width="250px")
  })

  output$t.len <- renderUI({
    if (input$f.choice!="sim") return();
    sliderInput("t.len","Length of TS", min = 1, max = 200, value = 50, width="250px")
  })

  output$n.sd <- renderUI({
    if (input$f.choice!="sim") return();
    sliderInput("n.sd","Noise SD:", min=0, max = 1, value = 0.05, width="125px");
  })

  output$a.f <- renderUI({
    if (input$f.choice!="sim") return();
    sliderInput("a.f",HTML("&omega;, Ang. Freq:"), min=0, max = 0.5, value = 0.1, step = 0.01, width="125px");
  })

  output$file <- renderUI({
    if (input$f.choice!="upload") return();
    fileInput('file', 'Choose CSV File', accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
  })
  output$sep <- renderUI({if (input$f.choice!="upload") return(); radioButtons('sep', 'Separator', c(","=',', ":"=':', ";"=';', Tab='\t'), ',', inline = TRUE)})
  output$header <- renderUI({if (input$f.choice!="upload") return(); checkboxInput('header', 'Header', TRUE)})

  output$ts.selected = renderText({
    if (input$f.choice=="upload" && is.null(input$file)) return("<b>Select a 'csv' file that contain the time series in its columns</b>")
    if (input$f.choice=="upload") {Ts <- as.matrix(read.table(input$file$datapath, header=input$header, sep=input$sep))
    } else if (input$f.choice=="server") {if (is.null(input$s.choice)) return(); i <- as.numeric(input$s.choice); Ts <- Xs[[i]];
    } else {simul <- simulate(); Ts <- simul$Trs + simul$noise};
    if (!length(Ts)) return()
    if (is.null(colnames(Ts))) {colnames(Ts) <- paste("fn",1:ncol(Ts))};
    if ("dmd" %in% input$dmd) {Ts <- Ts-mean(Ts); if (input$f.choice=="sim") Trs <- Trs-mean(Trs)} #{Ts <<- scale(Ts,scale=F); if (input$f.choice=="sim") Trs <<- scale(Trs,scale=F)}
    #updateSliderInput(session, "mssaL", max=min(ncol(Ts),trunc(nrow(Ts)/2))); updateSliderInput(session, "fssaL", max=min(ncol(Ts),trunc(nrow(Ts)/2)))
    updateSliderInput(session, "mssaL", max=trunc(ncol(Ts)/2)); updateSliderInput(session, "fssaL", max=min(120,trunc(ncol(Ts)/2)))
    updateSelectInput(session, "desc", selected = "ts"); updateSliderInput(session, "dimn", max=min(10,ncol(Ts)), value=min(2,ncol(Ts)));
    text <- paste("<b>",ncol(Ts),"Time series of length",nrow(Ts),"</b>")
    if (input$f.choice=="sim") iTrs(simul$Trs); iTs(Ts); return(text)
  })

  output$data <- renderTable({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    return(head(as.matrix(iTs()[,1:min(9,ncol(iTs()))]),15))
  })

  output$data.plot = renderPlot({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$f.choice=="server") {i <- as.numeric(input$s.choice); fname <- names(Xs)[i]}
    else if (input$f.choice=="upload") {fname <- input$file$name} else {fname <- "Simulation"};
    ts.plot(iTs(), main=paste("Time Series -", fname), ylab="", ylim=range(iTs()), gpars=list(xaxt="n"), xlab="tau")
    if (input$f.choice=="sim") for (i in 1:ncol(iTrs())) points(iTrs()[,i],type="l",col=2)
  })

  output$basis.n <- renderUI({
    sliderInput("basis.n", "Basis #:", min = 1, max =input$xdf, value = 1, step=1, width="400px")
  })

  output$basis.desc = renderPlot({
    if (is.null(input$basis.n)) return()
    xs <- seq(0,1,length.out=1000);
    if (input$bs.fr=="B-spline") Bx <- fda::bsplineS(xs,breaks=seq(0,1,length.out=input$xdf-input$xdeg+1),norder=input$xdeg+1)
    else Bx <- fda::fourier(xs,nbasis=input$xdf)
    ts.plot(Bx, col=8, main="B-spline Basis", xlab="Grid Points", gpars=list(xaxt="n"))
    points(Bx[,input$basis.n], type="l", lwd=2, col=2)
    axis(1,trunc(summary(1:nrow(Bx))[-4]))
    #b <- plot_ly(x=xs,y=Bx[,1],mode='lines',name = 'B1')
    #for (j in 2:ncol(Bx)) b <- b %>% add_lines(y=Bx[,j], mode = 'lines', name =paste0("B",j))
    #b
  })

  output$desc <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    choices <- list(Summary=c("Functional Time Series" = "ts", "How many basis? (GCV)"="gcv"),
                    FSSA=c("Scree"="fssa.scree", "W.Corr" = "fssa.wcor", "Paired"="fssa.pair", "Singular Vectors"="fssa.singV", "Peiodogram"="fssa.perGr", "Singular Functions"="fssa.singF", "Reconstruction"="fssa.reconst"),
                    SSA=c("Scree"="ssa.scree", "W.Corr" = "ssa.wcor", "Paired"="ssa.pair", "Vectors"="ssa.vec", "Functions" = "ssa.funs", "Reconstruction"="ssa.reconst"));
    selectInput("desc","Select",choices=choices, width="250px")
  })

  output$as.choice <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (is.null(input$desc)) { return() } else if (input$desc!="ts") return()
    radioButtons('as.choice', 'Plot Choices:', c("All" = "all", "Multiple"="mult", "Single" = "single"), selected= "all", inline = TRUE, width="400px")
  })

  output$sts.choice <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (((input$desc=="ts" && input$as.choice!="all")) && !(input$s.plot=="bf" && length(input$s.plot)==1) && length(input$s.plot))
      if (input$as.choice=="single") {
        sliderInput("sts.choice", "Choose function:", min = 1, max = ncol(iTs()), value = ifelse(is.null(input$sts.choice),1,input$sts.choice), step = 1, width="400px")
      } else {
        sliderInput("sts.choice", "Choose clusters:", min = 1, max = input$freq, value = ifelse(is.null(input$sts.choice),0,input$sts.choice), step = 1, width="200px")
      }
  })

  output$rec.type <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (is.null(input$desc)) { return() } else if (!input$desc%in%c("fssa.reconst","fssa.singF","ssa.reconst")) return()
    if (input$desc=="fssa.singF") selectInput("rec.type","Type",choices=c("Heat plot"="lheats","Regular Plot"="lcurves"), width="250px")
    else if (input$desc=="ssa.reconst") selectInput("rec.type","Type",choices=c("Heat Plot"="heatmap","Regular Plot"="line","3D Plot (line)"="3Dline","3D Plot (surface)"="3Dsurface","Old Plot"="1"), width="250px")
    else selectInput("rec.type","Type",choices=c("Heat Plot"="heatmap","Regular Plot"="line","3D Plot (line)"="3Dline","3D Plot (surface)"="3Dsurface","image2D"="3","ribbon3D"="1","ribbon3D-Curtain"="2"), width="250px")
  })

  output$freq <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (((input$desc=="ts" && input$as.choice=="mult")) && !(input$s.plot=="bf" && length(input$s.plot)==1) && length(input$s.plot))
      sliderInput("freq", "Period:", min = 1, max = trunc(ncol(iTs())/2), value = ifelse(is.null(input$freq),1,input$freq), step = 1, width="200px")
  })

  output$s.plot <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (is.null(input$desc)) { return() } else if (input$desc!="ts") return()
    choices <- c("Time Series" = "ts", "True Functions" = "tf", "Basis Func." = "bf", "Smoothing" = "bss", "Functional PCA" = "fpca", "Multivariate SSA" = "ssa", "Functional SSA" = "fssa") #, "Dyn. Func. PCA" = "dfpca"
    if (input$f.choice!="sim") {choices <- choices[-2]};
    checkboxGroupInput('s.plot', 'Plot:', choices=choices, selected=choices[1], width="250px")
  })

  output$b.indx <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || !input$desc%in%c("ts") || !length(intersect(input$s.plot,c("bf","bss")))) return();
    val <- c(1,input$xdf); if(length(input$b.indx)==2) val <- input$b.indx
    sliderInput("b.indx", "Basis Contr.", min = 1, max = input$xdf, value = val, dragRange=TRUE, step = 1, width="200px")
  })

  output$s.CI <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$desc=="ts" && length(intersect(input$s.plot,c("bss"))))
      if (input$as.choice=="single") checkboxInput('s.CI', 'Show CI', ifelse(is.null(input$s.CI), FALSE, input$s.CI), width="200px")
  })

  output$comp.obs  <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (is.null(input$desc)) { return() } else if (input$desc!="ts") return()
    checkboxInput('comp.obs', 'Compare fit. vs obs.', FALSE, width="200px")
  })

  output$run.fda.gcv <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (is.null(input$desc)) { return() } else if (input$desc!="gcv") return()
    actionButton('run.fda.gcv', paste('update GCV'))
  })

  fda.gcv <- eventReactive(input$run.fda.gcv, {
    withProgress(message = 'FDA.GCV: # of Basis', value = 0, { GCV <- NULL; df <- min(df,nrow(iTs()))
    if(input$bs.fr=="B-spline") nbasis <- (input$xdeg+1):(df-1) else nbasis <- seq(3,df,by=2);
    #nbasis <- nbasis[1:min(max(nbasis),nrow(iTs())-(input$xdeg+1))]
    for (l in 1:length(nbasis)) {
      if(input$bs.fr=="B-spline") bas.fssa <- fda::create.bspline.basis(c(0, 1), nbasis=nbasis[l], norder=input$xdeg+1)
      else bas.fssa <- fda::create.fourier.basis(c(0, 1), nbasis=nbasis[l]);
      GCV[l] <- sum(smooth.basis(seq(0,1,length.out = nrow(iTs())),iTs(),bas.fssa)$gcv);
      incProgress(1/length(nbasis), detail = l);
    }; itmp(1); return(list(GCV=GCV, nbasis=nbasis))
    })
  })

  output$res.plot <- renderPlot({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$desc%in%c("fssa.reconst","ssa.reconst") && !is.null(input$rec.type)) if (input$rec.type%in%c("heatmap","line","3Dline","3Dsurface")) return()
    if (input$f.choice=="server") {fname <- names(Xs)[as.numeric(input$s.choice)]} else if (input$f.choice=="upload") {fname <- input$file$name} else {fname <- "Simulation"}
    indx <- as.numeric(input$sts.choice); Ts <- iTs(); name.Ts <- names(iTs())
    if (length(intersect(input$s.plot,c("bf","bss")))) {
      if (is.null(input$b.indx)) return()
      b.indx <- input$b.indx[1]:input$b.indx[2]
    }
    if (length(intersect(input$s.plot,c("bss")))) {
      if(input$bs.fr=="B-spline") B <- fda::bsplineS(seq(0,1,length.out = nrow(Ts)), breaks=seq(0,1,length.out=input$xdf-input$xdeg+1), norder=input$xdeg+1)
      else B <- fda::fourier(seq(0,1,length.out = nrow(Ts)), nbasis=input$xdf)
      if ("bss" %in% input$s.plot) {
        cB <- solve(t(B)%*%B)%*%t(B);
        if (input$desc=="ts") {S <- B%*%cB; vB <- ifelse(nrow(Ts)>ncol(B), sum((Ts-S%*%Ts)^2)/((nrow(Ts)-ncol(B))*ncol(Ts))*diag(S), 0)}
      }
    }
    if ("bf" %in% input$s.plot || input$desc=="basis") {
      if(input$bs.fr=="B-spline") Bs <- fda::bsplineS(seq(0,1,length.out=1000), breaks=seq(0,1,length.out=input$xdf-input$xdeg+1), norder=input$xdeg+1)
      else Bs <- fda::fourier(seq(0,1,length.out=1000), nbasis=input$xdf); Bs <- Bs*sd(Ts)
    }
    #clust <- hclust(dist(t(Ts)),method="ward.D"); clcol <- cutree(clust,k=input$dimn);
    if (substr(input$desc,1,4)=="fssa" || sum(c("fssa","ssa")%in%input$s.plot) || substr(input$desc,1,3)=="ssa") {
      sr <- run_ssa(); input.g <- eval(parse(text=paste0("list(",input$g,")")))
      Q <- Qs <- matrix(0,nrow=nrow(Ts),ncol=ncol(Ts));
      isolate(sr$Qs <- reconstruct(sr$Us, groups=input.g))
      isolate(sr$Qf <- freconstruct(sr$Uf, input.g)); Qf <- sr$Qf[[1]][[1]]; Qf$coefs[,] <- 0
      for (i in input$sg[1]:input$sg[2]) {Qs <- Qs+t(sr$Qs[[i]]); Qf$coefs <- Qf$coefs+sr$Qf[[i]][[1]]$coefs}
      Q <- eval.fd(seq(0, 1, length = nrow(Ts)),Qf)
    }
    if ("fpca" %in% input$s.plot || "dfpca" %in% input$s.plot) {  f.pca <- run_fpca(); fpca <- f.pca$fpca; dfpca <- f.pca$dfpca;  }
    if (substr(input$desc,1,4)=="fssa") {
      if (input$desc=="fssa.scree") plot(sr$Uf,d=input$d[2],type="values")
      else if (input$desc=="fssa.wcor") plot(sr$Uf,d=input$d[2],type="wcor")
      else if (input$desc=="fssa.pair") plot(sr$Uf,d=input$d[2],type="paired")
      else if (input$desc=="fssa.singV") plot(sr$Uf, d=input$d[2], type="vectors")
      else if (input$desc=="fssa.perGr") plot(sr$Uf,d=input$d[2],type="periodogram")
      else if (input$desc=="fssa.singF") {plot(sr$Uf, d=input$d[2], type=ifelse(is.null(input$rec.type),"lheats",input$rec.type))}
      else if (input$desc=="fssa.reconst" && input$rec.type%in%c("1","2","3")) ftsplot(seq(0, 1, length = nrow(Ts)), 1:ncol(Ts), Qf, space = 0.1, type=as.numeric(input$rec.type), ylab = "tau", xlab = "t", main = "Q")
    } else if (substr(input$desc,1,3)=="ssa") {
      if (input$desc=="ssa.scree") plot(sr$Us,type="values",numvalues=input$d[2])
      else if (input$desc=="ssa.wcor") plot(sr$Us,type="wcor",groups=input$d[1]:input$d[2])
      else if (input$desc=="ssa.pair") plot(sr$Us,type="paired",idx=input$d[1]:input$d[2])
      else if (input$desc=="ssa.vec") plot(sr$Us,type="vectors",idx=input$d[1]:input$d[2])
      else if (input$desc=="ssa.funs") plot(sr$Us,type="series",groups=input$d[1]:input$d[2])
      else if (input$desc=="ssa.reconst") ts.plot(Qs,main="Reconstructed", ylim=range(Ts))
    } else if (input$desc=="gcv") {
      res <- fda.gcv(); ind.m <- which(res$GCV==min(res$GCV));
      plot(res$nbasis, res$GCV, type="b", xlab="n.basis", log="y", ylab="GCV", main=paste("Gen. Cross Validation -", fname), cex.lab=1.5, pch=20);
      if (itmp()) {updateSliderInput(session,'xdf',value=res$nbasis[ind.m]); itmp(0)}
      abline(v=input$xdf, col=1, lty=2);
    } else {
      m.lab <- paste(fname,"-",colnames(Ts)[indx])
      if (input$as.choice=="all") {indx <- 1:ncol(Ts); m.lab <- fname} else if (input$as.choice=="mult") {
        indt <- 1:ncol(Ts)%%input$freq; indt[indt==0] <- input$freq; indx <- (1:ncol(Ts))[which(indt==indx)]}
      if ("ts" %in% input$s.plot) {clcol <- rep(1,ncol(Ts))} else {clcol <- rep(0,ncol(Ts))}
      if ("dbl"%in%input$dmd) rng <- range(Ts,-Ts) else rng <- range(Ts)
      ts.plot(Ts[,indx], col=clcol[indx], main=paste("Time Series -",m.lab,"-",name.Ts), ylab="", ylim=rng, gpars=list(xaxt="n"));
      if ("ssa" %in% input$s.plot) {for (i in indx) points(Qs[,i],type="l",col=5)}
      if ("fssa" %in% input$s.plot) {for (i in indx) points(Q[,i],type="l",col=6)}
      if ("fpca" %in% input$s.plot) {for (i in indx) points(fpca[,i],type="l",col=2, lty=2)}
      if ("dfpca" %in% input$s.plot) {for (i in indx) points(dfpca[,i],type="l",col=3, lty=2)}
      if ("bss" %in% input$s.plot) { f.est <- matrix(0, nrow=nrow(Ts), ncol=ncol(Ts))
      for (i in indx) {
        f.est[,i] <- as.matrix(B[,b.indx])%*%(cB%*%Ts[,i])[b.indx];
        if (input$as.choice=="single") if (input$s.CI)
        {polygon(c(1:nrow(Ts),nrow(Ts):1), c(f.est[,i]-2*sqrt(vB),rev(f.est[,i])+2*rev(sqrt(vB))), border=4, lwd=1,col=5);
          if ("ts" %in% input$s.plot) points(Ts[,i],type="l",col=1)}
        points(f.est[,i], col=4, type="l")
      }
      }
      if ("tf" %in% input$s.plot) {for (i in indx) points(iTrs()[,i],type="l",col=2)}
      if ("bf" %in% input$s.plot) {for (i in b.indx) points(seq(1,nrow(Ts),length.out = 1000), Bs[,i],col=8,type="l",lty=2)}
      axis(1,trunc(summary(1:nrow(Ts))[-4]));
      if (input$f.choice=="sim" || input$comp.obs) {
        if (input$comp.obs) Ys <- Ts else Ys <- iTrs()
        if (input$comp.obs) RMSEs <- NULL else RMSEs <- paste(" RMSE.obs =", round(sqrt(mean((Ts[,indx]-Ys[,indx])^2)),4), "\n") #/length(indx)
        if ("bss" %in% input$s.plot) RMSEs <- paste(RMSEs, "RMSE.bs.smooth =", round(sqrt(mean((f.est[,indx]-Ys[,indx])^2)),4), "\n")
        if ("fpca" %in% input$s.plot) RMSEs <- paste(RMSEs, "RMSE.fpca =", round(sqrt(mean((fpca[,indx]-Ys[,indx])^2)),4), "\n")
        if ("ssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "RMSE.ssa =", round(sqrt(mean((Qs[,indx]-Ys[,indx])^2)),4), "\n")
        if ("fssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "RMSE.fssa =", round(sqrt(mean((Q[,indx]-Ys[,indx])^2)),4), "\n")
        if (input$comp.obs) RMSEs <- paste(RMSEs, "\n") else RMSEs <- paste(RMSEs, "\n Ang.obs =", round(max(acos(diag(t(scale(Ts[,indx],center=F))%*%scale(Ys[,indx],center=F))/(nrow(Ts)-1))*180/pi),4), "\n")
        if ("bss" %in% input$s.plot) RMSEs <- paste(RMSEs, "Ang.bs.smooth =", round(max(acos(diag(t(scale(f.est[,indx],center=F))%*%scale(Ys[,indx],center=F))/(nrow(Ts)-1))*180/pi),4), "\n")
        if ("fpca" %in% input$s.plot) RMSEs <- paste(RMSEs, "Ang.fpca =", round(max(acos(diag(t(scale(fpca[,indx],center=F))%*%scale(Ys[,indx],center=F))/(nrow(Ts)-1))*180/pi),4), "\n")
        if ("ssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "Ang.ssa =", round(max(acos(diag(t(scale(Qs[,indx],center=F))%*%scale(Ys[,indx],center=F))/(nrow(Ts)-1))*180/pi),4), "\n")
        if ("fssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "Ang.fssa =", round(max(acos(diag(t(scale(Q[,indx],center=F))%*%scale(Ys[,indx],center=F))/(nrow(Ts)-1))*180/pi),4), "\n")
        output$RMSEs <- renderText({ RMSEs })
      } else output$RMSEs <- renderText({ "Real Data" });
    }
  })

  output$res.ly <- renderPlotly({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$desc%in%c("fssa.reconst","ssa.reconst")) if (!input$rec.type%in%c("heatmap","line","3Dline","3Dsurface")) return()
    sr <- run_ssa(); input.g <- eval(parse(text=paste0("list(",input$g,")")));
    if (input$desc=="fssa.reconst") {
      isolate(sr$Qf <- freconstruct(sr$Uf, input.g)); Qf <- 0;
      for (i in input$sg[1]:input$sg[2]) {  Qf <- Qf + sr$Qf[[i]]  }
      plot(Qf, type=input$rec.type)
    } else {
      isolate(sr$Qs <- reconstruct(sr$Us, groups=input.g)); Qs <- matrix(0,nrow=nrow(iTs()),ncol=ncol(iTs()))
      for (i in input$sg[1]:input$sg[2]) {Qs <- Qs + t(sr$Qs[[i]])}
      plot(fts(smooth.basis(sr$tau,Qs,sr$bas.fssa)$fd),type=input$rec.type)
    }
  })

}

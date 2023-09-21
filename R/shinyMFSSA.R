# create the shiny application user interface
ui.mfssa <- fluidPage(
  tags$head(tags$style(HTML("body { max-width: 1250px !important; }"))),
  titlePanel("MFSSA Illustration"),
  sidebarLayout(
    sidebarPanel(
      width = 3, tags$head(tags$style(type = "text/css", ".well { max-width: 300px; }")),
      tags$div(title = "Pick your functional basis", radioButtons("bs.fr", "Choose Basis:", choices = c("B-spline", "Fourier"), selected = "B-spline", inline = TRUE)),
      tags$div(title = "Pick the degree of polynomial for the B-spline", uiOutput("xdeg", width = "250px")),
      tags$div(title = "Pick the number of basis", uiOutput("xdf", width = "250px")),
      tags$hr(style = "border-color: red;", width = "150px"),
      column(6, tags$div(title = "Grouping used for the SSA algorithms", textInput("g", "Groups", value = "1:2"))),
      column(6, tags$div(title = "Pick the groups fo reconstruction", uiOutput("sg"))),
      column(6, tags$div(title = "The dimensions used in FPCA, SSA and FSSA", uiOutput("d"))),
      column(6, tags$div(title = "See Manual", checkboxGroupInput("dmd.uf", "Functions", choices = c("Demean" = "dmd", "Dbl Range" = "dbl", "Univ. FSSA" = "uf")))),
      tags$div(title = "Window length parameter used in SSA and FSSA", sliderInput("ssaL", HTML("Win.L. (SSA):"), min = 1, max = 50, value = 20, step = 1, width = "210px")),
      column(6, uiOutput("run.ssa"))
    ),
    mainPanel(
      width = 9, tags$style(type = "text/css", ".shiny-output-error { visibility: hidden; }", ".shiny-output-error:before { visibility: hidden; },"), # ".nav-tabs {font-size: 10px}"),
      tabsetPanel(
        id = "Panel", type = "tabs",
        tabPanel(
          title = "Input Data", value = "Data",
          column(12, uiOutput("ts.selected", align = "center"), style = "color:red;"),
          fluidRow(
            column(4, radioButtons("f.choice", "Choose from:", c("Server" = "server", "Upload" = "upload", "Simulate" = "sim"), selected = "sim", inline = TRUE, width = "250px")),
            column(4, uiOutput("s.choice", width = "250px"), column(6, uiOutput("a1.f")), column(6, uiOutput("a1.l"))), column(2, uiOutput("noise.t", width = "125px")), column(2, uiOutput("noise.p", width = "125px")),
            column(4, uiOutput("file")), column(4, uiOutput("sep"), uiOutput("header"))
          ),
          column(4, uiOutput("model")), column(4, uiOutput("t.len")), column(2, uiOutput("a2.f")), column(2, uiOutput("n.sd")),
          column(12, plotOutput("data.plot", height = 500, width = 900))
        ),
        tabPanel(
          "Basis Functions",
          column(8, plotOutput("basis.desc", height = 600, width = 600)), column(4, uiOutput("basis.n", width = "300px"))
        ),
        tabPanel(
          "Data Analysis",
          column(4, uiOutput("desc", width = "250px")), column(4, uiOutput("as.choice", width = "400px"), uiOutput("run.fda.gcv", width = "200px"), uiOutput("rec.type", width = "300px")), column(2, uiOutput("freq")), column(2, uiOutput("sts.choice")),
          fluidRow(
            column(
              8, conditionalPanel(condition = "output.flag_plot", plotOutput("res.plot", height = 600, width = 600)),
              conditionalPanel(condition = "output.flag_plotly", plotlyOutput("res.ly", height = 600, width = 600))
            ),
            column(4, uiOutput("var.which"), uiOutput("s.plot"), fluidRow(column(8, uiOutput("b.indx")), column(4, uiOutput("s.CI"))), column(12, uiOutput("comp.obs"), verbatimTextOutput("RMSEs")))
          )
        ),
        tabPanel(
          "Forecasting",
          column(3, checkboxGroupInput("fcast.method", "Forecasting Method:", choices = c("Recurrent" = "recurrent", "Vector" = "vector"), selected = "recurrent", width = "250px")),
          column(4, uiOutput("fcast.horizon")), column(2, uiOutput("run.fcast")), column(3, uiOutput("fcast.type")),
          fluidRow(column(8, plotlyOutput("fcast.ly", height = 600, width = 600)), column(4, uiOutput("fcast.select"), uiOutput("fcast.var")))
        ),
        tabPanel("Manual", includeMarkdown(system.file("shiny/rmd", "report.Rmd", package = "Rfssa")))
      )
    )
  )
)

# Define server logic required to run mfssa
server.mfssa <- function(input, output, clientData, session) {
  iTs <- reactiveVal(list())
  iTrs <- reactiveVal(list())
  iXs <- reactiveVal(list())
  itmp <- reactiveVal(0)
  previous_s.plot <- reactiveVal(0)
  df <- 100
  vf <- 20
  T <- 100
  output$flag_plotly <- reactive(input$desc %in% c("mfssa.reconst", "ssa.reconst") && input$rec.type %in% c("heatmap", "line", "3Dline", "3Dsurface"))
  output$flag_plot <- reactive(!(input$desc %in% c("mfssa.reconst", "ssa.reconst") && input$rec.type %in% c("heatmap", "line", "3Dline", "3Dsurface")))
  outputOptions(output, "flag_plotly", suspendWhenHidden = FALSE)
  outputOptions(output, "flag_plot", suspendWhenHidden = FALSE)
  hideTab(inputId = "Panel", target = "Forecasting")
  updateTabsetPanel(session, "Panel", selected = "Manual")

  rfar <- function(N, norm, psi, Eps, basis) {
    # Create Corresponding matrix of an integral (kernel) operator
    # Get an kernel function k(s,t) corresponds to an integral operator and
    # a basis system incuding d basis functions, return corresponding d*d matrix of the operator with respect
    # to basis system.
    OpsMat <- function(kernel, basis) {
      u <- seq(0, 1, by = 0.01)
      n <- length(u)
      K_mat <- outer(u, u, FUN = kernel)
      K_t <- smooth.basis(u, K_mat, basis)$fd # the kernel function convert
      # to fd object w.r.t the first argument. K_t is an objcect of n fd.
      A <- inprod(K_t, basis) # An n*d matrix.
      K <- smooth.basis(u, A, basis)$fd # d fd object.
      B <- inprod(K, basis) # An d*d matrix.
      return(B) # return to OpsMat.
    }
    Psi_mat0 <- OpsMat(psi, basis)
    Gram <- inprod(basis, basis)
    Psi_mat <- solve(Gram) %*% Psi_mat0
    E <- Eps$coefs
    X <- E
    for (i in 2:N) X[, i] <- Psi_mat %*% X[, i - 1] + E[, i]
    X_fd <- fd(X, basis)
    return(X_fd) # return to rfar function.
  }
  gamma0 <- function(norm) {
    f <- function(x) {
      g <- function(y) psi0(x, y)^2
      return(integrate(g, 0, 1)$value) # return into f.
    }
    f <- Vectorize(f)
    A <- integrate(f, 0, 1)$value
    return(norm / A) # return into gamma.
  }
  psi0 <- function(x, y) 2 - (2 * x - 1)^2 - (2 * y - 1)^2

  fpca_proj <- function(i, U) {
    harm <- U$harmonics[i]
    scores <- U$scores[, i]
    m <- nrow(harm$coefs)
    n <- length(scores)
    coef <- matrix(NA, nrow = m, ncol = n)
    for (i0 in 1:m) for (j in 1:n) coef[i0, j] <- harm$coefs[i0] * scores[j]
    pc <- fd(coef, harm$basis)
    return(pc)
  }
  fpca_rec <- function(d1, d2, U) {
    s <- fpca_proj(d1, U)
    if (d2 > d1) for (i0 in (d1 + 1):d2) s <- s + fpca_proj(i0, U)
    return(s)
  }

  Tr1 <- function(tau, t) {
    tr1 <- ifelse("f1" %in% input$model, 1, 0) * (cos(2 * pi * input$a1.f * (t + input$a1.l)) * exp(tau^2) - sin(2 * pi * input$a1.f * t) * cos(4 * pi * tau))
    if ("f2" %in% input$model) tr1 <- tr1 + (cos(2 * pi * input$a2.f * t) * exp(1 - tau^2) + sin(2 * pi * input$a2.f * t) * sin(pi * tau))
    return(tr1)
  }
  Tr2 <- function(tau, t) {
    ifelse("f1" %in% input$model, 1, 0) * (sin(2 * pi * input$a1.f * t) * exp(tau^2) + cos(2 * pi * input$a1.f * (t - input$a1.l)) * cos(4 * pi * tau))
  }

  simulate <- function() {
    if (is.null(input$a1.f) || ("f2" %in% input$model && is.null(input$a2.f))) {
      return()
    }
    tau <- seq(0, 1, length = T)
    t <- 1:input$t.len
    noise <- list()
    Trs <- list(outer(tau, t, FUN = Tr1), outer(tau, t, FUN = Tr2))
    set.seed(T * input$t.len * input$n.sd)
    for (i in 1:length(Trs)) {
      noise[[i]] <- Z <- matrix(rnorm(input$t.len * T, 0, input$n.sd), nrow = T)
      if (input$noise.t == "ar1") {
        if (input$noise.p) {
          Z[1, ] <- 0
          A <- diag(1, T)
          if (T > 2) for (j in 1:(T - 2)) diag(A[-(1:j), ]) <- input$noise.p^j
          if (T > 1) A[T, 1] <- input$noise.p^(T - 1)
          noise[[i]] <- A %*% Z
        }
      }
      if (input$noise.t == "swn") {
        Z <- matrix(rnorm(input$t.len * input$xdf, 0, input$n.sd), ncol = input$t.len)
        basis.Z <- fda::create.bspline.basis(c(0, 1), input$xdf)
        tau <- seq(0, 1, length = T)
        basis.noise <- fda::fd(Z, basis.Z)
        noise[[i]] <- eval.fd(tau, basis.noise)
      }
      if (input$noise.t == "far1") {
        k0 <- gamma0(input$noise.p)
        psi <- function(x, y) k0 * psi0(x, y)
        Z[1, ] <- 0
        noise[[i]] <- apply(Z, 2, cumsum)
        if (input$bs.fr == "B-spline") {
          basis.Z <- fda::create.bspline.basis(c(0, 1), input$xdf)
        } else {
          basis.Z <- fda::create.fourier.basis(c(0, 1), input$xdf)
        }
        tau <- seq(0, 1, length = T)
        Eps <- smooth.basis(tau, noise[[i]], basis.Z)$fd
        basis.noise <- rfar(input$t.len, input$noise.p, psi, Eps, basis.Z)
        noise[[i]] <- eval.fd(tau, basis.noise)
      }
    }
    # Trs <- scale(Trs,scale=F);
    return(list(Trs = Trs, noise = noise))
  }

  observeEvent(input$dmd.uf, {
    updateTabsetPanel(session, "Panel", selected = "Data")
    updateCheckboxGroupInput(session, "s.plot", selected = "")
  })
  observeEvent(input$s.plot, {
    if (previous_s.plot() == 0 && (sum(c("mfssa", "fssa", "ssa") %in% input$s.plot))) previous_s.plot(1)
    if (previous_s.plot() == 1 && !(sum(c("mfssa", "fssa", "ssa") %in% input$s.plot))) previous_s.plot(0)
  })
  observeEvent(input$run.ssa, {
    showTab(inputId = "Panel", target = "Forecasting")
  })

  output$xdeg <- renderUI({
    if (input$bs.fr == "Fourier") {
      return()
    }
    sliderInput("xdeg", HTML("Degree of B-spline Basis:"), min = 0, max = 5, value = 3, step = 1, width = "250px")
  })

  output$xdf <- renderUI({
    sliderInput("xdf", paste("Deg. of freedom of", input$bs.fr, "Basis:"), min = ifelse(input$bs.fr == "B-spline", input$xdeg + 1, 3), max = df, value = vf, step = ifelse(input$bs.fr == "B-spline", 1, 2), width = "250px")
  })

  output$sg <- renderUI({
    m <- length(eval(parse(text = paste0("list(", input$g, ")"))))
    sliderInput("sg", "Select G:", min = 1, max = m, value = c(1, m), step = 1)
  })

  output$d <- renderUI({
    sliderInput("d", "d", min = 1, max = input$ssaL, value = c(1, 2), step = 1)
  })

  output$run.ssa <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    actionButton("run.ssa", paste("run M(F)SSA"))
  })

  run_ssa <- eventReactive(
    {
      input$run.ssa
      previous_s.plot()
    },
    {
      withProgress(message = "MSSA.MFSSA: Running", value = 0, {
        if (input$bs.fr == "B-spline") {
          bas.fssa <- fda::create.bspline.basis(c(0, 1), nbasis = input$xdf, norder = input$xdeg + 1)
        } else {
          bas.fssa <- fda::create.fourier.basis(c(0, 1), nbasis = input$xdf)
        }
        tau <- seq(0, 1, length = nrow(iTs()[[1]]))
        uUf <- list()
        if ("uf" %in% input$dmd.uf) {
          for (i in 1:length(iTs())) {
            uUf[[i]] <- fssa(funts(X = iTs()[[i]], basisobj = bas.fssa), input$ssaL)
          }
        }
        mUf <- fssa(funts(X = iTs(), basisobj = rep(list(bas.fssa), length(iTs()))), input$ssaL)
        Ys <- NULL
        for (i in 1:length(iTs())) {
          Ys <- cbind(Ys, t(iTs()[[i]]))
        }
        Us <- ssa(Ys, input$ssaL, kind = "mssa")
        return(list(Us = Us, mUf = mUf, uUf = uUf, tau = tau, bas.fssa = bas.fssa))
      })
    }
  )

  run_fpca <- function() {
    withProgress(message = "(D)FPCA: Running", value = 0, {
      if (input$bs.fr == "B-spline") {
        bas.fssa <- fda::create.bspline.basis(c(0, 1), nbasis = input$xdf, norder = input$xdeg + 1)
      } else {
        bas.fssa <- fda::create.fourier.basis(c(0, 1), nbasis = input$xdf)
      }
      tau <- seq(0, 1, length = nrow(iTs()[[1]]))
      Y <- smooth.basis(tau, iTs()[[ifelse(is.null(input$var.which), 1, as.numeric(input$var.which))]], bas.fssa)$fd
      f.pca <- pca.fd(Y, nharm = min(input$d[2], input$xdf), centerfns = FALSE)
      fpca.rec <- fpca_rec(min(input$d[1], input$xdf), min(input$d[2], input$xdf), f.pca)
      fpca.re <- eval.fd(tau, fpca.rec)
      # d.fpca <- fts.dpca(Y, Ndpc = input$d[2]); dfpca.re <- eval.fd(tau,d.fpca$Xhat)
      return(list(fpca = fpca.re, dfpca = fpca.re))
    })
  }

  output$s.choice <- renderUI({
    if (input$f.choice != "server") {
      return()
    }
    if (!length(iXs())) {
      loadCallcenterData()
      loadJambiData()
      Callcenter <- get("Callcenter", envir = .GlobalEnv)
      Jambi <- get("Jambi", envir = .GlobalEnv)

      Xs <- list()
      Xs[[1]] <- matrix(sqrt(Callcenter$calls), nrow = 240)
      Xs[[2]] <- Xs[[3]] <- matrix(NA, nrow = 128, ncol = dim(Jambi$NDVI)[3])
      for (i in 1:dim(Jambi$NDVI)[3]) {
        Xs[[2]][, i] <- density(Jambi$NDVI[, , i], from = 0, to = 1, n = 128)$y
        Xs[[3]][, i] <- density(Jambi$EVI[, , i], from = 0, to = 1, n = 128)$y
      }
      colnames(Xs[[2]]) <- colnames(Xs[[3]]) <- Jambi$Date
      names(Xs) <- c("Callcenter", "NDVI", "EVI")
      Xs[[4]] <- Xs[2:3]
      names(Xs) <- c(names(Xs[1:3]), "xDI")
      iXs(Xs)
    }
    s.choices <- 1:length(iXs())
    names(s.choices) <- names(iXs())
    selectInput("s.choice", "Select a file from server: ", choices = s.choices, width = "250px")
  })

  output$noise.t <- renderUI({
    if (input$f.choice != "sim") {
      return()
    }
    selectInput("noise.t", "Type of noice: ", choices = c("AR(1)" = "ar1", "FAR(1)" = "far1", "Smooth WN" = "swn"), width = "125px")
  })

  output$noise.p <- renderUI({
    if (input$f.choice != "sim") {
      return()
    }
    if (!is.null(input$noise.t)) {
      if (input$noise.t == "swn") {
        return()
      }
    }
    sliderInput("noise.p", "AR Parameter:", min = 0, max = 1, value = 0, step = 0.01, width = "125px")
  })

  output$model <- renderUI({
    if (input$f.choice != "sim") {
      return()
    }
    choices <- c("A.1t(\u03C91,\u2113) F1 +" = "f1", "A.2t(\u03C92) F2" = "f2")
    checkboxGroupInput("model", "Model:", choices = choices, selected = "f1", inline = TRUE, width = "250px")
  })

  output$t.len <- renderUI({
    if (input$f.choice != "sim") {
      return()
    }
    sliderInput("t.len", "Length of TS", min = 1, max = 200, value = 50, width = "250px")
  })

  output$n.sd <- renderUI({
    if (input$f.choice != "sim") {
      return()
    }
    sliderInput("n.sd", "Noise SD:", min = 0, max = 1, value = 0.05, width = "125px")
  })

  output$a1.f <- renderUI({
    if (input$f.choice != "sim" || !("f1" %in% input$model)) {
      return()
    }
    sliderInput("a1.f", HTML("&omega;1, Ang. Freq:"), min = 0, max = 0.5, value = 0.1, step = 0.01, width = "125px")
  })

  output$a1.l <- renderUI({
    if (input$f.choice != "sim" || !("f1" %in% input$model)) {
      return()
    }
    sliderInput("a1.l", HTML("Par. &ell; in A<sub>1t</sub>:"), min = -10, max = 10, value = 0, step = 1, width = "125px")
  })

  output$a2.f <- renderUI({
    if (input$f.choice != "sim" || !("f2" %in% input$model)) {
      return()
    }
    sliderInput("a2.f", HTML("&omega;2, Ang. Freq:"), min = 0, max = 0.5, value = 0.25, step = 0.01, width = "125px")
  })

  output$file <- renderUI({
    if (input$f.choice != "upload") {
      return()
    }
    fileInput("file", "Choose CSV File", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
  })
  output$sep <- renderUI({
    if (input$f.choice != "upload") {
      return()
    }
    radioButtons("sep", "Separator", c("," = ",", ":" = ":", ";" = ";", Tab = "\t"), ",", inline = TRUE)
  })
  output$header <- renderUI({
    if (input$f.choice != "upload") {
      return()
    }
    checkboxInput("header", "Header", TRUE)
  })

  output$ts.selected <- renderText({
    if (input$f.choice == "upload" && is.null(input$file)) {
      return("<b>Select a 'csv' file that contain the time series in its columns</b>")
    }
    if (input$f.choice == "upload") {
      Ts <- as.matrix(read.table(input$file$datapath, header = input$header, sep = input$sep))
      if (!input$header && !is.numeric(Ts)) {
        headers <- as.factor(Ts[1, ])
        TS <- list()
        for (l in levels(headers)) TS[[l]] <- matrix(as.numeric(Ts[-1, headers == l]), nrow = nrow(Ts) - 1)
        Ts <- TS
      } else {
        Ts <- list(Ts)
      }
    } else if (input$f.choice == "server") {
      if (is.null(input$s.choice)) {
        return()
      }
      i <- as.numeric(input$s.choice)
      Ts <- iXs()[[i]]
      if (is.matrix(Ts)) Ts <- list(Ts)
    } else {
      simul <- simulate()
      Ts <- Map("+", simul$Trs, simul$noise)
    }
    if (!length(Ts)) {
      return()
    }
    if (is.null(colnames(Ts[[1]]))) {
      for (i in 1:length(Ts)) colnames(Ts[[i]]) <- paste("fn", 1:ncol(Ts[[i]]))
    }
    if ("dmd" %in% input$dmd.uf) {
      for (i in 1:length(Ts)) {
        Ts[[i]] <- Ts[[i]] - mean(Ts[[i]])
        if (input$f.choice == "sim") simul$Trs[[i]] <- simul$Trs[[i]] - mean(simul$Trs[[i]])
      }
    }
    if (is.null(names(Ts)) || sum(is.na(names(Ts)))) names(Ts) <- paste("Variable", 1:length(Ts))
    updateSelectInput(session, "desc", selected = "ts")
    updateSliderInput(session, "dimn", max = min(10, ncol(Ts[[1]])), value = min(2, ncol(Ts[[1]])))
    updateSliderInput(session, "ssaL", max = min(120, trunc(ncol(Ts[[1]]) / 2)))
    text <- paste("<b>", ncol(Ts[[1]]), ifelse(length(Ts) == 1, "Univariate", "Multivariate"), "Time series of length", nrow(Ts[[1]]), "</b>")
    if (input$f.choice == "sim") iTrs(simul$Trs)
    iTs(Ts)
    return(text)
  })

  output$data.plot <- renderPlot({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model)) || !length(iTs())) {
      return()
    }
    if (input$f.choice == "server") {
      i <- as.numeric(input$s.choice)
      fname <- names(iXs())[i]
    } else if (input$f.choice == "upload") {
      fname <- input$file$name
    } else {
      fname <- "Simulation"
    }
    par(mfrow = c(1, length(iTs())), mgp = c(1.5, .5, 0), mar = c(3, 3, 2.5, 1.75))
    for (i in 1:length(iTs())) {
      ts.plot(iTs()[[i]], main = paste("Time Series -", fname, "-", names(iTs()[i])), ylab = "", ylim = range(iTs()[[i]]), gpars = list(xaxt = "n"), xlab = "tau")
      if (input$f.choice == "sim") for (j in 1:ncol(iTrs()[[i]])) points(iTrs()[[i]][, j], type = "l", col = 2)
    }
  })

  output$basis.n <- renderUI({
    sliderInput("basis.n", "Basis #:", min = 1, max = input$xdf, value = 1, step = 1, width = "400px")
  })

  output$basis.desc <- renderPlot({
    if (is.null(input$basis.n)) {
      return()
    }
    xs <- seq(0, 1, length.out = 1000)
    if (input$bs.fr == "B-spline") {
      Bx <- fda::bsplineS(xs, breaks = seq(0, 1, length.out = input$xdf - input$xdeg + 1), norder = input$xdeg + 1)
    } else {
      Bx <- fda::fourier(xs, nbasis = input$xdf)
    }
    ts.plot(Bx, col = 8, main = "B-spline Basis", xlab = "Grid Points", gpars = list(xaxt = "n"))
    points(Bx[, input$basis.n], type = "l", lwd = 2, col = 2)
    axis(1, trunc(summary(1:nrow(Bx))[-4]))
  })

  output$desc <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    choices <- list(
      Summary = c("Functional Time Series" = "ts", "How many basis? (GCV)" = "gcv"),
      MFSSA = c("Scree" = "mfssa.scree", "W.Correlation" = "mfssa.wcor", "Paired" = "mfssa.pair", "Singular Vectors" = "mfssa.singV", "Periodogram" = "mfssa.perGr", "Singular Functions" = "mfssa.singF", "Reconstruction" = "mfssa.reconst"),
      MSSA = c("Scree" = "ssa.scree", "W.Correlation" = "ssa.wcor", "Paired" = "ssa.pair", "Singular Vectors" = "ssa.vec", "Functions" = "ssa.funs", "Reconstruction" = "ssa.reconst")
    )
    selectInput("desc", "Select Plot Type", choices = choices, width = "250px")
  })

  output$as.choice <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (is.null(input$desc)) {
      return()
    } else if (input$desc != "ts") {
      return()
    }
    radioButtons("as.choice", "Plot Choices:", c("All" = "all", "Multiple" = "mult", "Single" = "single"), selected = "all", inline = TRUE, width = "400px")
  })

  output$sts.choice <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (((input$desc == "ts" && input$as.choice != "all")) && !("bf" %in% input$s.plot && length(input$s.plot) == 1) && length(input$s.plot)) {
      if (input$as.choice == "single") {
        sliderInput("sts.choice", "Choose function:", min = 1, max = ncol(iTs()[[1]]), value = ifelse(is.null(input$sts.choice), 1, input$sts.choice), step = 1, width = "400px")
      } else {
        sliderInput("sts.choice", "Choose clusters:", min = 1, max = input$freq, value = ifelse(is.null(input$sts.choice), 0, input$sts.choice), step = 1, width = "200px")
      }
    }
  })

  output$rec.type <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (is.null(input$desc)) {
      return()
    } else if (!input$desc %in% c("mfssa.reconst", "mfssa.singF", "ssa.reconst")) {
      return()
    }
    if (input$desc == "mfssa.singF") {
      selectInput("rec.type", "Type", choices = c("Heat plot" = "lheats", "Regular Plot" = "lcurves"), width = "250px")
    } else if (input$desc == "ssa.reconst") {
      selectInput("rec.type", "Type", choices = c("Heat Plot" = "heatmap", "Regular Plot" = "line", "3D Plot (line)" = "3Dline", "3D Plot (surface)" = "3Dsurface"), width = "250px")
    } else {
      selectInput("rec.type", "Type", choices = c("Heat Plot" = "heatmap", "Regular Plot" = "line", "3D Plot (line)" = "3Dline", "3D Plot (surface)" = "3Dsurface"), width = "250px")
    }
  })

  output$freq <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (((input$desc == "ts" && input$as.choice == "mult")) && !("bf" %in% input$s.plot && length(input$s.plot) == 1) && length(input$s.plot)) {
      sliderInput("freq", "Period:", min = 1, max = trunc(ncol(iTs()[[1]]) / 2), value = ifelse(is.null(input$freq), 1, input$freq), step = 1, width = "200px")
    }
  })

  output$var.which <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model)) || is.null(input$desc)) {
      return()
    }
    if (!input$desc %in% c("ts", "gcv", "ssa.reconst", "mfssa.reconst", "mfssa.singF") || length(iTs()) == 1) {
      return()
    }
    choic <- as.character(1:length(iTs()))
    names(choic) <- names(iTs())
    if (!is.null(input$rec.type) && input$desc == "mfssa.reconst" && length(choic) != 1) if (input$rec.type %in% c("line")) choic <- c("All Variables" = "all", choic)
    selectInput("var.which", NULL, choices = choic, selected = "1", width = "125px")
  })

  output$s.plot <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (is.null(input$desc)) {
      return()
    } else if (input$desc != "ts") {
      return()
    }
    choices <- c("Time Series (Raw Data)" = "ts", "True Functions" = "tf", "Multivariate FSSA" = "mfssa", "Multivariate SSA" = "ssa", "Functional PCA" = "fpca", "Basis Function" = "bf", "Smoothing" = "bss") # , "Dyn. Func. PCA" = "dfpca"
    if ("uf" %in% input$dmd.uf) {
      choices <- append(choices, "fssa", 3)
      names(choices)[4] <- "Functional SSA"
    }
    if (input$f.choice != "sim") {
      choices <- choices[-2]
    }
    checkboxGroupInput("s.plot", "Plot:", choices = choices, selected = choices[1], width = "250px")
  })

  output$b.indx <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || !input$desc %in% c("ts") || !length(intersect(input$s.plot, c("bf", "bss")))) {
      return()
    }
    val <- c(1, input$xdf)
    if (length(input$b.indx) == 2) val <- input$b.indx
    sliderInput("b.indx", "Basis Contr.", min = 1, max = input$xdf, value = val, dragRange = TRUE, step = 1, width = "200px")
  })

  output$s.CI <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (input$desc == "ts" && length(intersect(input$s.plot, c("bss")))) {
      if (input$as.choice == "single") checkboxInput("s.CI", "Show CI", ifelse(is.null(input$s.CI), FALSE, input$s.CI), width = "200px")
    }
  })

  output$comp.obs <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (is.null(input$desc)) {
      return()
    } else if (input$desc != "ts") {
      return()
    }
    checkboxInput("comp.obs", "Compare fit. vs obs.", FALSE, width = "200px")
  })

  output$run.fda.gcv <- renderUI({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (is.null(input$desc)) {
      return()
    } else if (input$desc != "gcv") {
      return()
    }
    actionButton("run.fda.gcv", paste("update GCV"))
  })

  fda.gcv <- eventReactive(input$run.fda.gcv, {
    withProgress(message = "FDA.GCV: # of Basis", value = 0, {
      GCV <- NULL
      df <- min(df, nrow(iTs()[[1]]))
      if (input$bs.fr == "B-spline") nbasis <- (input$xdeg + 1):(df - 1) else nbasis <- seq(3, df, by = 2)
      for (l in 1:length(nbasis)) {
        if (input$bs.fr == "B-spline") {
          bas.fssa <- fda::create.bspline.basis(c(0, 1), nbasis = nbasis[l], norder = input$xdeg + 1)
        } else {
          bas.fssa <- fda::create.fourier.basis(c(0, 1), nbasis = nbasis[l])
        }
        GCV[l] <- sum(smooth.basis(seq(0, 1, length.out = nrow(iTs()[[1]])), iTs()[[ifelse(is.null(input$var.which), 1, min(length(iTs()), as.numeric(input$var.which)))]], bas.fssa)$gcv)
        incProgress(1 / length(nbasis), detail = l)
      }
      itmp(1)
      return(list(GCV = GCV, nbasis = nbasis))
    })
  })

  output$res.plot <- renderPlot({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (input$desc %in% c("mfssa.reconst", "ssa.reconst") && !is.null(input$rec.type)) {
      if (input$rec.type %in% c("heatmap", "line", "3Dline", "3Dsurface")) {
        return()
      }
    }
    if (!is.null(input$var.which) && length(iTs()) != 1) var.which <- as.numeric(input$var.which) else var.which <- 1
    if (input$f.choice == "server") {
      fname <- names(iXs())[as.numeric(input$s.choice)]
    } else if (input$f.choice == "upload") {
      fname <- input$file$name
    } else {
      fname <- "Simulation"
    }
    indx <- as.numeric(input$sts.choice)
    Ts <- iTs()[[var.which]]
    name.Ts <- names(iTs())[var.which]
    if (length(intersect(input$s.plot, c("bf", "bss")))) {
      if (is.null(input$b.indx)) {
        return()
      }
      b.indx <- input$b.indx[1]:input$b.indx[2]
    }
    if (length(intersect(input$s.plot, c("bss")))) {
      if (input$bs.fr == "B-spline") {
        B <- fda::bsplineS(seq(0, 1, length.out = nrow(Ts)), breaks = seq(0, 1, length.out = input$xdf - input$xdeg + 1), norder = input$xdeg + 1)
      } else {
        B <- fda::fourier(seq(0, 1, length.out = nrow(Ts)), nbasis = input$xdf)
      }
      if ("bss" %in% input$s.plot) {
        cB <- solve(t(B) %*% B) %*% t(B)
        if (input$desc == "ts") {
          S <- B %*% cB
          vB <- ifelse(nrow(Ts) > ncol(B), sum((Ts - S %*% Ts)^2) / ((nrow(Ts) - ncol(B)) * ncol(Ts)) * diag(S), 0)
        }
      }
    }
    if ("bf" %in% input$s.plot || input$desc == "basis") {
      if (input$bs.fr == "B-spline") {
        Bs <- fda::bsplineS(seq(0, 1, length.out = 1000), breaks = seq(0, 1, length.out = input$xdf - input$xdeg + 1), norder = input$xdeg + 1)
      } else {
        Bs <- fda::fourier(seq(0, 1, length.out = 1000), nbasis = input$xdf)
      }
      Bs <- Bs * sd(Ts)
    }
    if (substr(input$desc, 1, 5) == "mfssa" || sum(c("mfssa", "fssa", "ssa") %in% input$s.plot) || substr(input$desc, 1, 3) == "ssa") {
      sr <- run_ssa()
      input.g <- eval(parse(text = paste0("list(", input$g, ")")))
      mQ <- uQ <- Qs <- matrix(0, nrow = nrow(Ts), ncol = ncol(Ts))
      isolate(sr$Qs <- reconstruct(sr$Us, groups = input.g))
      if ("uf" %in% input$dmd.uf) {
        isolate(sr$uQf <- freconstruct(sr$uUf[[var.which]], input.g))
        uQf <- sr$uQf[[1]]
        uQf$coefs[[1]][, ] <- 0
      }
      isolate(sr$mQf <- freconstruct(sr$mUf, input.g))
      mQf <- sr$mQf[[1]][, var.which]
      mQf$coefs[[1]][, ] <- 0
      for (i in input$sg[1]:input$sg[2]) {
        Qs <- Qs + t(sr$Qs[[i]][, ((var.which - 1) * nrow(Ts) + 1):(var.which * nrow(Ts))])
        if ("uf" %in% input$dmd.uf) {
          uQf$coefs[[1]] <- uQf$coefs[[1]] + sr$uQf[[i]]$coefs[[1]]
        }
        mQf$coefs[[1]] <- mQf$coefs[[1]] + sr$mQf[[i]]$coefs[[var.which]]
      }
      if ("uf" %in% input$dmd.uf) uQ <- uQf$B_mat[[1]] %*% uQf$coefs[[1]]
      mQ <- mQf$B_mat[[1]] %*% mQf$coefs[[1]]
    }
    if ("fpca" %in% input$s.plot || "dfpca" %in% input$s.plot) {
      f.pca <- run_fpca()
      fpca <- f.pca$fpca
      dfpca <- f.pca$dfpca
    }
    if (substr(input$desc, 1, 5) == "mfssa") {
      if (input$desc == "mfssa.scree") {
        plot(sr$mUf, type = "values", d = input$d[2])
      } else if (input$desc == "mfssa.wcor") {
        mfwcor <- fwcor(sr$mUf, groups = input$d[1]:input$d[2])
        wplot(W = mfwcor, cuts = 10)
      } # plot(sr$mUf,type="wcor",groups=input$d[1]:input$d[2])
      else if (input$desc == "mfssa.pair") {
        plot(sr$mUf, type = "paired", idx = input$d[1]:input$d[2])
      } else if (input$desc == "mfssa.singV") {
        plot(sr$mUf, type = "vectors", idx = input$d[1]:input$d[2])
      } else if (input$desc == "mfssa.perGr") {
        plot(sr$mUf, type = "periodogram", idx = input$d[1]:input$d[2])
      } else if (input$desc == "mfssa.singF") {
        plot(sr$mUf, type = ifelse(is.null(input$rec.type), "lheats", input$rec.type), var = var.which, idx = input$d[1]:input$d[2])
      }
    } else if (substr(input$desc, 1, 3) == "ssa") {
      if (input$desc == "ssa.scree") {
        plot(sr$Us, type = "values", numvalues = input$d[2])
      } else if (input$desc == "ssa.wcor") {
        plot(sr$Us, type = "wcor", groups = input$d[1]:input$d[2])
      } else if (input$desc == "ssa.pair") {
        plot(sr$Us, type = "paired", idx = input$d[1]:input$d[2])
      } else if (input$desc == "ssa.vec") {
        plot(sr$Us, type = "vectors", idx = input$d[1]:input$d[2])
      } else if (input$desc == "ssa.funs") plot(sr$Us, type = "series", groups = input$d[1]:input$d[2])
    } else if (input$desc == "gcv") {
      res <- fda.gcv()
      ind.m <- which(res$GCV == min(res$GCV))
      plot(res$nbasis, res$GCV, type = "b", xlab = "n.basis", log = "y", ylab = "GCV", main = paste("Gen. Cross Validation -", fname), cex.lab = 1.5, pch = 20)
      if (itmp()) {
        updateSliderInput(session, "xdf", value = res$nbasis[ind.m])
        itmp(0)
      }
      abline(v = input$xdf, col = 1, lty = 2)
    } else {
      m.lab <- paste(fname, "-", colnames(Ts)[indx])
      if (input$as.choice == "all") {
        indx <- 1:ncol(Ts)
        m.lab <- fname
      } else if (input$as.choice == "mult") {
        indt <- 1:ncol(Ts) %% input$freq
        indt[indt == 0] <- input$freq
        indx <- (1:ncol(Ts))[which(indt == indx)]
      }
      if ("ts" %in% input$s.plot) {
        clcol <- rep(1, ncol(Ts))
      } else {
        clcol <- rep(0, ncol(Ts))
      }
      if ("dbl" %in% input$dmd.uf) rng <- range(Ts, -Ts) else rng <- range(Ts)
      ts.plot(Ts[, indx], col = clcol[indx], main = paste("Time Series -", m.lab, "-", name.Ts), ylab = "", ylim = rng, gpars = list(xaxt = "n"))
      if ("ssa" %in% input$s.plot) {
        for (i in indx) points(Qs[, i], type = "l", col = 5)
      }
      if ("fssa" %in% input$s.plot) {
        for (i in indx) points(uQ[, i], type = "l", col = 7)
      }
      if ("mfssa" %in% input$s.plot) {
        for (i in indx) points(mQ[, i], type = "l", col = 6)
      }
      if ("fpca" %in% input$s.plot) {
        for (i in indx) points(fpca[, i], type = "l", col = 2, lty = 2)
      }
      if ("dfpca" %in% input$s.plot) {
        for (i in indx) points(dfpca[, i], type = "l", col = 3, lty = 2)
      }
      if ("bss" %in% input$s.plot) {
        f.est <- matrix(0, nrow = nrow(Ts), ncol = ncol(Ts))
        for (i in indx) {
          f.est[, i] <- as.matrix(B[, b.indx]) %*% (cB %*% Ts[, i])[b.indx]
          if (input$as.choice == "single") {
            if (input$s.CI) {
              polygon(c(1:nrow(Ts), nrow(Ts):1), c(f.est[, i] - 2 * sqrt(vB), rev(f.est[, i]) + 2 * rev(sqrt(vB))), border = 4, lwd = 1, col = 5)
              if ("ts" %in% input$s.plot) points(Ts[, i], type = "l", col = 1)
            }
          }
          points(f.est[, i], col = 4, type = "l")
        }
      }
      if ("tf" %in% input$s.plot) {
        for (i in indx) points(iTrs()[[var.which]][, i], type = "l", col = 2)
      }
      if ("bf" %in% input$s.plot) {
        for (i in b.indx) points(seq(1, nrow(Ts), length.out = 1000), Bs[, i], col = 8, type = "l", lty = 2)
      }
      axis(1, trunc(summary(1:nrow(Ts))[-4]))
      if (input$f.choice == "sim" || input$comp.obs) {
        if (input$comp.obs) Ys <- Ts else Ys <- iTrs()[[var.which]]
        if (input$comp.obs) RMSEs <- NULL else RMSEs <- paste(" RMSE.obs =", round(sqrt(mean((Ts[, indx] - Ys[, indx])^2)), 4), "\n")
        if ("bss" %in% input$s.plot) RMSEs <- paste(RMSEs, "RMSE.bs.smooth =", round(sqrt(mean((f.est[, indx] - Ys[, indx])^2)), 4), "\n")
        if ("fpca" %in% input$s.plot) RMSEs <- paste(RMSEs, "RMSE.fpca =", round(sqrt(mean((fpca[, indx] - Ys[, indx])^2)), 4), "\n")
        if ("ssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "RMSE.ssa =", round(sqrt(mean((Qs[, indx] - Ys[, indx])^2)), 4), "\n")
        if ("fssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "RMSE.fssa =", round(sqrt(mean((uQ[, indx] - Ys[, indx])^2)), 4), "\n")
        if ("mfssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "RMSE.mfssa =", round(sqrt(mean((mQ[, indx] - Ys[, indx])^2)), 4), "\n")
        if (input$comp.obs) RMSEs <- paste(RMSEs, "\n") else RMSEs <- paste(RMSEs, "\n Ang.obs =", round(max(acos(diag(t(scale(Ts[, indx], center = F)) %*% scale(Ys[, indx], center = F)) / (nrow(Ts) - 1)) * 180 / pi), 4), "\n")
        if ("bss" %in% input$s.plot) RMSEs <- paste(RMSEs, "Ang.bs.smooth =", round(max(acos(diag(t(scale(f.est[, indx], center = F)) %*% scale(Ys[, indx], center = F)) / (nrow(Ts) - 1)) * 180 / pi), 4), "\n")
        if ("fpca" %in% input$s.plot) RMSEs <- paste(RMSEs, "Ang.fpca =", round(max(acos(diag(t(scale(fpca[, indx], center = F)) %*% scale(Ys[, indx], center = F)) / (nrow(Ts) - 1)) * 180 / pi), 4), "\n")
        if ("ssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "Ang.ssa =", round(max(acos(diag(t(scale(Qs[, indx], center = F)) %*% scale(Ys[, indx], center = F)) / (nrow(Ts) - 1)) * 180 / pi), 4), "\n")
        if ("fssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "Ang.fssa =", round(max(acos(diag(t(scale(uQ[, indx], center = F)) %*% scale(Ys[, indx], center = F)) / (nrow(Ts) - 1)) * 180 / pi), 4), "\n")
        if ("mfssa" %in% input$s.plot) RMSEs <- paste(RMSEs, "Ang.mfssa =", round(max(acos(diag(t(scale(mQ[, indx], center = F)) %*% scale(Ys[, indx], center = F)) / (nrow(Ts) - 1)) * 180 / pi), 4), "\n")
        output$RMSEs <- renderText({
          RMSEs
        })
      } else {
        output$RMSEs <- renderText({
          "Real Data"
        })
      }
    }
  })

  output$res.ly <- renderPlotly({
    if ((input$f.choice == "upload" && is.null(input$file)) || (input$f.choice == "sim" && !length(input$model))) {
      return()
    }
    if (input$desc %in% c("mfssa.reconst", "ssa.reconst")) {
      if (!input$rec.type %in% c("heatmap", "line", "3Dline", "3Dsurface")) {
        return()
      }
    }
    sr <- run_ssa()
    input.g <- eval(parse(text = paste0("list(", input$g, ")")))
    if (input$desc == "mfssa.reconst") {
      isolate(sr$mQf <- freconstruct(sr$mUf, input.g))
      mQf <- sr$mQf[[1]]
      for (i in 1:length(mQf$dimSupp)) mQf$coefs[[i]][, ] <- 0
      if (input$var.which == "all" || is.null(input$var.which)) {
        var.which <- NULL
        types <- rep(input$rec.type, length(mQf$dimSupp))
        vars <- 1:length(mQf$dimSupp)
      } else {
        var.which <- as.numeric(input$var.which)
        types <- input$rec.type
        vars <- var.which
      }
      for (i in input$sg[1]:input$sg[2]) {
        mQf <- mQf + sr$mQf[[i]]
      }
      myplot <- plotly_funts(mQf[, vars], types = types, main = paste("Variable", vars), zlab = paste("Variable", vars))
    } else {
      if (is.null(input$var.which)) var.which <- 1 else var.which <- as.numeric(input$var.which)
      isolate(sr$Qs <- reconstruct(sr$Us, groups = input.g))
      Qs <- matrix(0, nrow = nrow(iTs()[[1]]), ncol = ncol(iTs()[[1]]))
      for (i in input$sg[1]:input$sg[2]) {
        Qs <- Qs + t(sr$Qs[[i]][, ((var.which - 1) * nrow(iTs()[[1]]) + 1):(var.which * nrow(iTs()[[1]]))])
      }
      myplot <- plotly_funts(funts(X = Qs, basisobj = sr$bas.fssa), types = input$rec.type, main = paste("Variable", var.which), zlab = paste("Variable", var.which))
    }
    print(myplot[[1]])
  })

  output$fcast.horizon <- renderUI({
    if (is.null(nrow(iTs()[[1]]))) {
      return()
    }
    sliderInput("fcast.horizon", HTML("Forecasting Horizon:"), min = 0, max = ncol(iTs()[[1]]), value = min(50, ncol(iTs()[[1]]) / 2), step = 1, width = "250px")
  })

  output$run.fcast <- renderUI({
    if (is.null(nrow(iTs()[[1]]))) {
      return()
    }
    actionButton("run.fcast", paste("run Forecast"))
  })

  output$fcast.type <- renderUI({
    if (is.null(input$run.fcast)) {
      return()
    }
    selectInput("fcast.type", "Type", choices = c("Heat Plot" = "heatmap", "Regular Plot" = "line", "3D Plot (line)" = "3Dline", "3D Plot (surface)" = "3Dsurface"), width = "250px")
  })

  run_fcast <- eventReactive(input$run.fcast, {
    withProgress(message = "MFSSA Forecast: Running", value = 0, {
      sr <- run_ssa()
      input.g <- eval(parse(text = paste0("list(", input$g, ")")))
      fc <- list()
      if ("recurrent" %in% input$fcast.method) fc$rec <- fforecast(U = sr$mUf, groups = input.g, h = input$fcast.horizon, method = "recurrent")
      if ("vector" %in% input$fcast.method) fc$vector <- fforecast(U = sr$mUf, groups = input.g, h = input$fcast.horizon, method = "vector")
      return(fc)
    })
  })

  output$fcast.var <- renderUI({
    if (is.null(input$run.fcast) || input$run.fcast == 0 || length(iTs()) == 1) {
      return()
    }
    choic <- as.character(1:length(iTs()))
    names(choic) <- names(iTs())
    if (!is.null(input$fcast.type) && length(choic) != 1) if (input$fcast.type %in% c("line")) choic <- c("All Variables" = "all", choic)
    selectInput("fcast.var", NULL, choices = choic, width = "125px")
  })

  output$fcast.select <- renderUI({
    if (is.null(input$run.fcast) || input$run.fcast == 0) {
      return()
    }
    if (length(input$fcast.method) > 1) selectInput("fcast.select", "Select Output", choices = c("Recurrent Forecasting" = "1", "Vector Forecasting" = "2"), width = "250px")
  })

  output$fcast.ly <- renderPlotly({
    if (is.null(input$run.fcast)) {
      return()
    }
    fc <- run_fcast()
    if (length(fc) == 1) i <- 1 else i <- as.numeric(input$fcast.select)
    mQf <- fc[[i]][[1]]
    mQf$coefs[[1]][, ] <- 0
    if (input$fcast.var == "all" || is.null(input$fcast.var)) {
      types <- rep(input$fcast.type, length(fc[[1]][[1]]$coefs))
      vars <- 1:length(fc[[1]][[1]]$coefs)
    } else {
      types <- input$fcast.type
      vars <- as.numeric(input$fcast.var)
    }
    for (j in input$sg[1]:input$sg[2]) {
      mQf <- mQf + fc[[i]][[j]]
    }
    myplot <- plotly_funts(mQf[, vars], types = types, main = paste("Variable", vars), zlab = paste("Variable", vars))
    print(myplot[[1]])
  })
}

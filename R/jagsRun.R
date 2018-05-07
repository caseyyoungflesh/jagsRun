#' Run JAGS
#'
#' Run JAGS in parallel and output output of interest. Number of cores used equals number of chains specified. Be sure that your machine has an adequate number of (virtual) cores available to run the model. Function creates a directory with \code{jagsID} name, saves .rds file with model output, and produces output summary in text file format.
#'
#' @param jagsData List containing data to feed to JAGS
#' @param jagsModel JAGS model file
#' @param jagsInits Initial values for JAGS model. Should be a list of lists (number of embedded lists should equal the number of chains being run in the model). NOTE: each chain should specify a different starting value for a particular parameter and/or use a different seed/RNG to avoid identical chains.
#' @param params Character string or vector of character strings specifying which parameters to track
#' @param jagsID OPTIONAL. Character string with name of jags run (e.g., 'Run_1')
#' @param jagsDsc OPTIONAL. Character string with description of the jags run (e.g., 'First model run')
#' @param db_hash OPTIONAL. Character string with description of data version which will be printed in the output file. Could be latest git commit hash.
#' @param n_chain Numeric specifying number of chains to be run
#' @param n_adapt Numeric specifying how many iterations to use for adaptation
#' @param n_burn Numeric specifying how any iterations to use for burn-in
#' @param n_draw Numeric specifying how many iterations to use for draw (iterations to be kept beyond adaptation and burn-in)
#' @param n_thin Numeric specifying thinning rate
#' @param TYPE Character string specifying which type of model to run. Options are 'STANDARD' (normal model run), COMPILE' (model is comiled only with adapt = 2 [no samples drawn] to check that the likelihood can be calculated [priors and inits are appropriate]), and 'DEBUG'(model run with one chain with 10 samples for adaptation, 10 samples for burn-in, 10 samples kept)
#' @param EXTRA Logical used to specify whether extra iterations should be run if convergence is not met. If \code{TRUE}, up to \code{n_max} iterations are run until convergence is reached (specified by \code{Rhat_max})
#' @param RANDOM Logical specifying whether to use script to generate random inits. If \code{TRUE}, \code{jagsInits} should be a function that generates random initial values.
#' @param Rhat_max Numeric specifying the maximum Rhat value allowed when \code{EXTRA = TRUE}
#' @param n_rburn Numeric specifying how many samples to use for burn in if \code{EXTRA = TRUE} and convergence (defined by \code{Rhat_max}) has not been reached.
#' @param n_max Numeric specifying the maximum number of samples to be drawn when \code{EXTRA = TRUE}
#' @param params_extra Character string or vector of character strings specifying which parameters to evaluate convergence for when \code{EXTRA = TRUE}. Must be a subset of \code{params}.
#' @param params_report Character string or vector of character strings specifying which parameters to report. Must be a subset of \code{params}.
#' @param ppc Character string or vector of character strings specifying the name of elements used for the posteriod predictive check (PPC). If specified, the summary information for these elements will be output in the report.
#' @param obj_out Logical specifying whether MCMC.list object should be output
#' @param save_data Logical specifying whether input data to function should be saved as a .rds object
#' @section Notes: jagsData should be formatted as such: XXXXX. jagsInits should be formatted as such: XXXXX. Jags params should be formatted as such: XXXXX.
#'
#' @export




jagsRun <- function (jagsData,
                     jagsModel,
                     jagsInits,
                     params,
                     jagsID,
                     jagsDsc,
                     db_hash,
                     n_chain = 3,
                     n_adapt = 5000,
                     n_burn,
                     n_draw,
                     n_thin = 1,
                     TYPE = 'STANDARD',
                     EXTRA = FALSE,
                     RANDOM = FALSE,
                     Rhat_max = 1.05,
                     n_rburn = 0,
                     n_max,
                     params_extra = params,
                     params_report = params,
                     ppc = NULL,
                     obj_out = FALSE,
                     save_data = FALSE)
{
  if (TYPE == 'COMPILE')
  {
    if (RANDOM == TRUE) start <- jagsInits(jagsData) else start <- jagsInits[[1]]

    invisible(rjags::load.module('glm'))
    jm = rjags::jags.model(data = jagsData,
                           file = jagsModel,
                           inits = start,
                           n.chains = 1,
                           n.adapt = 2)

    if(jagsModel %in% list.files())
    {
      invisible(file.remove(jagsModel))
    }

    print('Successful!')
  }

  if (TYPE == 'DEBUG')
  {
    if (RANDOM == TRUE) start <- jagsInits(jagsData) else start <- jagsInits[[1]]

    ptm <- proc.time()
    invisible(rjags::load.module('glm'))
    jm = rjags::jags.model(data = jagsData,
                           file = jagsModel,
                           inits = start,
                           n.chains = 1,
                           n.adapt = 10)

    stats::update(jm, n.iter = 10)

    out = rjags::coda.samples(jm,
                              n.iter = 10,
                              variable.names = params,
                              thin = 1)

    tt <- (proc.time() - ptm)[3] / 60

    if(jagsModel %in% list.files())
    {
      invisible(file.remove(jagsModel))
    }
  }

  if (TYPE == 'STANDARD')
  {
    CONVERGE <- FALSE
    cl <- parallel::makeCluster(n_chain)
    pid <- NA

    for (i in 1:n_chain)
    {
      pidNum <- utils::capture.output(cl[[i]])
      start <- regexpr("pid", pidNum)[[1]]
      end <- nchar(pidNum)
      pid[i] <- substr(pidNum, (start + 4), end)
    }

    parallel::clusterExport(cl, c('pid', 'jagsData', 'n_adapt', 'n_burn', 'n_draw', 'n_thin', 'n_rburn', 'params', 'jagsInits', 'jagsModel', 'RANDOM'), envir = environment())

    ptm <- proc.time()
    out.1 <- parallel::clusterEvalQ(cl,{
      require(rjags)
      if (RANDOM == TRUE) start <- jagsInits(jagsData) else {
        processNum <- which(pid == Sys.getpid())
        start <- jagsInits[[processNum]]
      }

      invisible(rjags::load.module('glm'))
      jm = rjags::jags.model(data = jagsData,
                             file = jagsModel,
                             inits = start,
                             n.chains = 1,
                             n.adapt = n_adapt)

      stats::update(jm, n.iter = n_burn)

      samples = rjags::coda.samples(jm,
                                    n.iter = n_draw,
                                    variable.names = params,
                                    thin = n_thin)
      return(samples)
    })

    tt <- (proc.time() - ptm)[3] / 60
    i <- 1
    a <- vector("list", n_chain)
    while(i <= n_chain)
    {
      a[[i]] <- out.1[[i]][[1]]
      i <- i + 1
    }

    out <- coda::mcmc.list(a)
    n_extra <- 0
    n_total <- n_adapt + n_burn + n_draw * n_chain

    if (max(MCMCvis::MCMCsummary(out, params = params_report, Rhat = TRUE)[,6]) <= Rhat_max) CONVERGE <- TRUE

    if (EXTRA == TRUE)
    {

      while(max(MCMCvis::MCMCsummary(out, params = params_extra, Rhat = TRUE)[,6]) > Rhat_max & n_total < n_max) {
        out.2 <- parallel::clusterEvalQ(cl, {
          stats::update(jm, n.iter = n_rburn)
          samples = rjags::coda.samples(jm,
                                        n.iter = n_draw,
                                        variable.names = params,
                                        thin = n_thin)
          return(samples)
        })

        i <- 1
        a <- vector("list", n_chain)
        while(i <= n_chain)
        {
          a[[i]] <- out.2[[i]][[1]]
          i <- i + 1
        }

        out <- coda::mcmc.list(a)
        tt <- (proc.time() - ptm)[3] / 60
        n_extra <- n_extra + n_rburn + n_draw
        n_total <- n_adapt + n_burn + n_draw * n_chain
        if (max(MCMCvis::MCMCsummary(out, params = params_extra, Rhat = TRUE)[,6]) <= Rhat_max) CONVERGE <- TRUE
      }
    }

    parallel::stopCluster(cl)
    closeAllConnections()
    s_out <- MCMCvis::MCMCsummary(out, params = params_report, n.eff = TRUE, digits = 4)
    options(max.print = 50000)

    if (missing(jagsID))
    {
      jagsID <- 'jagsRun_output'
    }
    dir.create(jagsID)

    #move .jags file into dir
    if(jagsModel %in% list.files())
    {
      invisible(file.rename(from = paste0(jagsModel), to = paste0(jagsID, '/', jagsModel)))
    }

    sink(paste0(jagsID, '/results.txt'))
    cat(paste0('jagsID: ', jagsID, ' \n'))
    if (!missing(jagsDsc))
    {
      cat(paste0('jagsDsc: ', jagsDsc, ' \n'))
    } else {
      cat(paste0('jagsDsc: NONE GIVEN', ' \n'))
    }
    if (!missing(db_hash))
    {
      cat(paste0('db_hash: ', db_hash, ' \n'))
    } else {
      cat(paste0('db_hash: NONE GIVEN', ' \n'))
    }
    cat(paste0('Random Inits: ', RANDOM, ' \n'))
    cat(paste0("Inits object: ", as.character(deparse(substitute(jagsInits))), ' \n'))
    cat(paste0('Total minutes: ', round(tt, digits = 2), ' \n'))
    cat(paste0('Total iterations: ', n_total, ' \n'))
    cat(paste0('n_chain: ', n_chain, ' \n'))
    cat(paste0('n_adapt: ', n_adapt, ' \n'))
    cat(paste0('n_burn: ', n_burn, ' \n'))
    cat(paste0('n_draw: ', n_draw, ' \n'))
    cat(paste0('n_thin: ', n_thin, ' \n'))
    cat(paste0('Total samples kept: ', n_chain * (n_draw / n_thin), ' \n'))
    cat(paste0('Extended burnin: ', EXTRA, ' \n'))

    if (EXTRA == TRUE) {
      cat(paste0('Rhat_max: ', Rhat_max, ' \n'))
      cat(paste0('n_extra: ', n_extra, ' \n'))
    }

    cat(paste0('convergence: ', CONVERGE, ' \n'))

    if (is.null(ppc) == FALSE) {
      cat(paste0('ppc: ', MCMCvis::MCMCsummary(out, params = ppc, n.eff = TRUE, digits = 4)[, 1], '\n'))
    }

    cat(' \n')
    print(s_out)
    sink()

    saveRDS(out, paste0(jagsID, '/', jagsID, '.rds'))

    if (save_data == TRUE)
    {
      saveRDS(jagsData, paste0(jagsID, '/jagsData.rds'))
    }

    if(obj_out == TRUE)
    {
      return(out)
    }
  }
}

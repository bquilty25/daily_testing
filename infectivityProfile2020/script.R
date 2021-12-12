#' script.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
library(readxl)
library(tidyverse)
library(parallel)
library(reshape2)
library(ggplot2)
library(cowplot)


#' Parallel lapply function
ParLapply <- function(X, FUN, ..., PARALLEL = TRUE, SEED = NULL) {
  if (PARALLEL) {
    if (is.null(SEED)) SEED <- 111
    num.cores <- min(c(length(X), max(detectCores() - 1, 1)))
    cl <- makeCluster(num.cores, type = "FORK")
    clusterSetRNGStream(cl, SEED)
    out <- parLapply(cl, X, FUN, ...)
    stopCluster(cl)
    return(out)
  } else {
    out <- lapply(X, FUN, ...)
    return(out)
  }
}


# Load and format serial interval data ----
serialData <- data.frame(read_xlsx("data/HeEtAlnatmed2020-fig1cData.xlsx"))
#' Report dates as day of year
refDate <- as.Date("2020-01-01")
serialData$x.lb <- as.numeric(as.Date(serialData$x.lb) - refDate)
serialData$x.ub <- as.numeric(as.Date(serialData$x.ub) - refDate)
serialData$y <- as.numeric(as.Date(serialData$y) - refDate)
names(serialData) <- c("ID","indexOnsetLower","indexOnsetUpper","secondaryOnset")
serialData$minSerialInterval <- serialData$secondaryOnset - serialData$indexOnsetUpper
serialData$maxSerialInterval <- serialData$secondaryOnset - serialData$indexOnsetLower


# Define integration parameters ----
stepSize <- 0.01
tMin <- -50
tMax <- max(serialData$secondaryOnset)


# Define the lognormal distribution for the incubation period ----
getIncubationDistribution <- function(t) {
  #defaultIncPars <- read_csv("data/HeEtAlnatMed2020-incubation.csv", col_types = "dd")
  dlnorm(x = t, meanlog = 1.434065, sdlog = 0.6612)
}


# Define the shifted gamma distribution for the infectivity profile ----
getInfectivityProfile <- function(t, params) {
  ifelse(t + params[["shift"]] <= 1e-8, # Avoid  divergence for shape -> 0
         0,
         dgamma(x = t + params[["shift"]], shape = params[["shape"]], rate = params[["rate"]])
  )
}


# Calculate the probability of observing each serial interval under a range of parameters ----
getProbabilities <- function(all.params, data, tMin, stepSize) {
  #' Function to calculate probabilities of ecah data point for a given set of parameters
  calculateProbabilities <- function(params, tI, serialIntervalIndices, infectivityIndices, incubationDistribution, incubationIndices, stepSize) {
    #' Evaluate the infectivity profile using params
    infectivityProfile <- getInfectivityProfile(tI, params)
    #' Compute the probability of observing each serial interval
    probS <- sapply(seq_along(infectivityIndices), function(sIndex) {
      sum(
        infectivityProfile[seq(infectivityIndices[[sIndex]]["first"],infectivityIndices[[sIndex]]["last"])] *
          incubationDistribution[seq(incubationIndices[[sIndex]]["first"],incubationIndices[[sIndex]]["last"])]
      ) * stepSize
    })
    #' For each data row we integrate the pdf of all available serial intervals
    probRow <- sapply(serialIntervalIndices, function(this.indices) {
      sum(probS[seq(this.indices["first"],this.indices["last"])]) * stepSize
    })
    return(probRow)
  }

  #' Determine all the serial intervals that we need for each row of the data
  serialIntervalsPerRow <- apply(data, 1, function(x) round(x[["secondaryOnset"]] - seq(x[["indexOnsetUpper"]] + 0.5, x[["indexOnsetLower"]] - 0.5, -stepSize), digits = 6))
  #' Keep only the unique ones
  serialIntervals <- sort(unique(unlist(serialIntervalsPerRow)))
  #' Match back the unique serial intervals to the data
  serialIntervalIndices <- lapply(serialIntervalsPerRow, function(this.intervals) {
    indices <- which(serialIntervals %in% this.intervals)
    return(c(first = indices[1], last = indices[length(indices)]))
  })
  #' Define possible infection times (this is the integration variable)
  tI <- seq(tMin, max(serialIntervals), stepSize)

  #' Indices for the infectivity profile that we need to use for each serial interval
  infectivityIndices <- lapply(serialIntervals, function(s) {
    indices <- which(tI <= s)
    return(c(first = indices[1], last = indices[length(indices)]))
  })

  #' Distribution of incubation periods from t = 0
  incubationDistribution <- getIncubationDistribution(tI - tMin)
  #' Indices for the incubation distribution needed for each unique serial interval (reversed for easy convolution with infectivity profile)
  incubationIndices <- lapply(infectivityIndices, function(this.infectivityIndices) {
    c(first = this.infectivityIndices[["last"]] - this.infectivityIndices[["first"]] + 1, last = 1)
  })

  #' Loop over parameter sets and compute probabilities of each serial interval
  probabilities <- ParLapply(seq_len(nrow(all.params)), function(i) {
    calculateProbabilities(
      params = c(shape = all.params[i,"shape"], rate = all.params[i,"rate"], shift = all.params[i,"shift"]),
      tI = tI,
      serialIntervalIndices = serialIntervalIndices,
      infectivityIndices = infectivityIndices,
      incubationDistribution = incubationDistribution,
      incubationIndices = incubationIndices,
      stepSize = stepSize
    )
  })
  return(probabilities)
}


# Grid evaluation functions ----
computeGrid <- function(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz) {
  x <- seq(xmin, xmax, length.out = nx + 1)
  y <- seq(ymin, ymax, length.out = ny + 1)
  z <- seq(zmin, zmax, length.out = nz + 1)
  cells <- lapply(seq_len(nx), function(i) {
    tmp <- lapply(seq_len(ny), function(j) {
      tmptmp <- lapply(seq_len(nz), function(k) {
        data.frame(xmin = x[[i]], xmax = x[[i + 1]],
                   ymin = y[[j]], ymax = y[[j + 1]],
                   zmin = z[[k]], zmax = z[[k + 1]],
                   x = round((x[[i]] + x[[i + 1]])/2, digits = 6),
                   y = round((y[[j]] + y[[j + 1]])/2, digits = 6),
                   z = round((z[[k]] + z[[k + 1]])/2, digits = 6)
        )
      })
      do.call(rbind, tmptmp)
    })
    do.call(rbind, tmp)
  })
  cells <- do.call(rbind, cells)
  return(cells)
}


addLikelihoods <- function(cells) {
  probabilities <- getProbabilities(all.params = data.frame(shape = cells$x, rate = cells$y, shift = cells$z), data = serialData, tMin = tMin, stepSize = stepSize)
  cells$llh.all <- sapply(probabilities, function(p) sum(log(p)))
  cells$llh.missing <- sapply(probabilities, function(p) sum(log(p[p > 0])))
  return(cells)
}


refineGrid <- function(cells, nx, ny, nz, par, threshold) {
  indices <- which(cells[,par] > threshold)
  print(paste0("computing grid for ", length(indices), " new cells"))
  new.cells <- lapply(indices, function(id) {
    out <- computeGrid(cells[id,"xmin"], cells[id,"xmax"], nx, cells[id,"ymin"], cells[id,"ymax"], ny, cells[id,"zmin"], cells[id,"zmax"], nz)
    return(out)
  })
  new.cells <- do.call(rbind, new.cells)
  print(paste0("calculating ", nrow(new.cells), " likelihoods"))
  new.cells <- addLikelihoods(new.cells)
  rbind(cells[-indices,], new.cells)
}


# Evaluate ----
#' We start from a coarse grid and refine
# df <- cbind(computeGrid(0,100, 1, 0,10, 1, 0,40, 1), llh.all = -Inf, llh.missing = Inf)
# df1 <- refineGrid(df, 20, 20, 20, "llh.missing", -Inf)
# df2 <- refineGrid(df1, 2, 2, 2, "llh.missing", -225)
# df3 <- refineGrid(df2, 2, 2, 2, "llh.missing", -225)
# df4 <- refineGrid(df3, 2, 2, 2, "llh.missing", -225)
# save(df,df1,df2,df3,df4, file = "data/gridLikelihoods.RData")


# Analysis ----
load("data/gridLikelihoods.RData")
## Plot likelihood maps ----
plot.df <- df1
plot.df$shift <- factor(plot.df$z, levels = sort(unique(plot.df$z)))

ggplot(plot.df, aes(x = x, y = y, fill = llh.missing, width = xmax - xmin, height = ymax - ymin)) +
  geom_tile() +#colour = "black") +
  geom_point(data = plot.df[plot.df$llh.missing == max(plot.df$llh.missing),], colour = "red") +
  facet_wrap(~shift, labeller = "label_both") +
  scale_x_continuous(limits = c(0,NA)) +
  scale_y_continuous(limits = c(0,NA)) +
  scale_fill_viridis_c() +
  labs(x = "shape", y = "rate") +
  theme_minimal()

ggplot(plot.df, aes(x = x, y = y, fill = llh.all, width = xmax - xmin, height = ymax - ymin)) +
  geom_tile() +#colour = "black") +
  geom_point(data = plot.df[plot.df$llh.all == max(plot.df$llh.all),], colour = "red") +
  facet_wrap(~shift, labeller = "label_both") +
  scale_x_continuous(limits = c(0,NA)) +
  scale_y_continuous(limits = c(0,NA)) +
  scale_fill_viridis_c() +
  labs(x = "shape", y = "rate") +
  theme_minimal()

## Format the likelihod data ----
dff <- df4
likelihoods <- data.frame(shape = dff$x, rate = dff$y, shift = dff$z, llh.missing = dff$llh.missing, llh.all = dff$llh.all)
#' Which parameters have the maximum likelihood
mleIndices <- list(
  missing = which(likelihoods$llh.missing == max(likelihoods$llh.missing)),
  all = which(likelihoods$llh.all == max(likelihoods$llh.all))
)
#' Which parameter combinations have high likleihood
threshold <- qchisq(p = 0.95, df = 3)/2
indices <- list(
  missing = which(likelihoods$llh.missing > likelihoods[mleIndices[["missing"]], "llh.missing"] - threshold),
  all = which(likelihoods$llh.all > likelihoods[mleIndices[["all"]], "llh.all"] - threshold)
)

## Maximum likelihood parameters ----
likelihoods[mleIndices$missing,]
likelihoods[mleIndices$all,]

## Plot parameter distributions ----
#' Plot only the high likelihood parameters
param.df <- sapply(c("missing","all"), function(type) {
  this.df <- likelihoods[indices[[type]], c("shape","rate","shift")]
  if (type == "missing") this.df$llh <- likelihoods[indices[[type]], "llh.missing"]
  else this.df$llh <- likelihoods[indices[[type]], "llh.all"]
  this.df$type <- factor(type, levels = c("missing", "all"))
  return(this.df)
}, simplify = F)
param.df <- do.call(rbind, param.df)
param.df <- melt(param.df, id.vars = c("type","llh"), variable.name = "parameter")
#' Dummy data for consistent plot limits across facets
dummy <- sapply(c("missing","all"), function(type) {
  data.frame(
    type = factor(type, levels = c("missing", "all")),
    llh = ifelse(type == "missing", max(likelihoods$llh.missing), max(likelihoods$llh.all)),
    parameter = factor(rep(c("shape","rate","shift"), each = 2), levels = c("shape","rate","shift")),
    value = c(0, max(likelihoods$shape), 0, max(likelihoods$rate), 0, max(likelihoods$shift))
  )
}, simplify = F, USE.NAMES = T)
#' Plot
plot_grid(
  plotlist = lapply(c("missing","all"), function(type) {
    ggplot(param.df[param.df$type == type, ], aes(x = llh, y = value, colour = llh)) +
      geom_point() +
      geom_point(data = param.df[param.df$llh == max(param.df[param.df$type == type,"llh"]),], colour = "red") +
      geom_blank(data = dummy[[type]]) +
      facet_grid(parameter ~ type, scales = "free") +
      scale_colour_viridis_c() +
      theme_minimal() + theme(legend.position = "bottom")
  }), nrow = 1, align = "hv", axis = "tb")


## Calculate the corresponding infectivity profiles (pdf and CDF) ----
times <- seq(-10, 25, 0.1)
#' First the maximum likelihood functions
mleGamma <- sapply(c("missing","all"), function(type) {
  pars <- likelihoods[mleIndices[[type]], c("shape","rate","shift")]
  data.frame(
    t = times,
    pdf = dgamma(x = times + pars$shift, shape = pars$shape, rate = pars$rate),
    CDF = pgamma(q = times + pars$shift, shape = pars$shape, rate = pars$rate),
    type = factor(type, levels = c("missing", "all"))
  )
}, simplify = F)
mleGamma <- do.call(rbind, mleGamma)
#' Now get the confidence interval from the high-likelihood functions
ribbonGamma <- sapply(c("missing","all"), function(type) {
  ribbon <- ParLapply(times, function(t) {
    pdf <- sapply(indices[[type]], function(i) dgamma(x = t + likelihoods[[i,"shift"]], shape = likelihoods[[i,"shape"]], rate = likelihoods[[i,"rate"]]))
    CDF <- sapply(indices[[type]], function(i) pgamma(q = t + likelihoods[[i,"shift"]], shape = likelihoods[[i,"shape"]], rate = likelihoods[[i,"rate"]]))
    data.frame(
      t = t,
      pdf = mean(pdf),
      pdf.lower = min(pdf),
      pdf.upper = max(pdf),
      CDF = mean(CDF),
      CDF.lower = min(CDF),
      CDF.upper = max(CDF),
      type = factor(type, levels = c("missing", "all"))
    )
  })
  do.call(rbind, ribbon)
}, simplify = F)
ribbonGamma <- do.call(rbind, ribbonGamma)
#' Finally the maximum likelihood estimate from He et al.
heGamma <- data.frame(
  t = times,
  pdf = dgamma(x = times + 2.3066912325330207, shape = 2.1157789506000704, rate = 0.6898582881923863),
  CDF = pgamma(q = times + 2.3066912325330207, shape = 2.1157789506000704, rate = 0.6898582881923863)
)
#' Save these datasets
save(mleGamma, ribbonGamma, heGamma, file = "data/gammaDistributions.RData")
load("data/gammaDistributions.RData")

## Plot infectivity profiles
pdfPlot <- ggplot(ribbonGamma, aes(x = t, y = pdf, colour = type, fill = type)) +
  geom_vline(xintercept = 0, colour = "darkgrey", size = 1) +
  geom_ribbon(aes(ymin = pdf.lower, ymax = pdf.upper, alpha = type), colour = "transparent") +
  scale_alpha_manual(values = c(missing = 0.2, all = 0.5), guide = "none") +
  geom_line(data = mleGamma, aes(linetype = type, size = type)) +
  scale_linetype_manual(values = c(missing = "dashed", all = "solid"), guide = "none") +
  scale_size_manual(values = c(missing = 0.5, all = 1.5), guide = "none") +
  #geom_line(data = heGamma, aes(x = t, y = pdf), colour = "black", size = 0.2, inherit.aes = F) +
  scale_colour_manual(values = c(missing = "#1b7cff", all = "#ff9e1b"), labels = c(missing = "He et al.", all = "corrected"), name = NULL, aesthetics = c("colour","fill")) +
  coord_cartesian(xlim = c(-10,10), ylim = c(0,0.301), expand = F) +
  labs(x = "days after symptom onset", y = "infectivity") +
  #ggtitle("A: probability density") +
  theme_bw() + theme(legend.position = c(0.8,0.8), plot.background = element_rect(colour = "white", fill = "white"), panel.grid = element_blank())

cdfPlot <- ggplot(ribbonGamma, aes(x = t, y = CDF, colour = type, fill = type)) +
  geom_vline(xintercept = 0, colour = "darkgrey", size = 1) +
  geom_ribbon(aes(ymin = CDF.lower, ymax = CDF.upper, alpha = type), colour = "transparent") +
  scale_alpha_manual(values = c(missing = 0.2, all = 0.5), guide = "none") +
  geom_line(data = mleGamma, aes(linetype = type, size = type)) +
  scale_linetype_manual(values = c(missing = "dashed", all = "solid"), guide = "none") +
  scale_size_manual(values = c(missing = 0.5, all = 1.5), guide = "none") +
  #geom_line(data = heGamma, aes(x = t, y = pdf), colour = "black", size = 0.2, inherit.aes = F) +
  scale_colour_manual(values = c(missing = "#1b7cff", all = "#ff9e1b"), labels = c(missing = "He et al.", all = "corrected"), name = NULL, aesthetics = c("colour","fill")) +
  scale_y_continuous(limits = c(0,1), labels = scales::percent) +
  coord_cartesian(xlim = c(-10,10), expand = F) +
  labs(x = "days after symptom onset", y = "cumulative infectivity") +
  #ggtitle("A: probability density") +
  theme_bw() + theme(legend.position = "none", plot.background = element_rect(colour = "white", fill = "white"), panel.grid = element_blank())

plot_grid(pdfPlot, cdfPlot, nrow = 1, align = "hv", axis = "tb", labels = "AUTO")
ggsave("infectivityProfiles.pdf", width = 170, height = 60, units = "mm")
ggsave("infectivityProfiles.eps", device = cairo_ps, fallback_resolution = 600, width = 170, height = 60, units = "mm")
ggsave("infectivityProfiles.jpg", width = 170, height = 60, units = "mm", dpi = 600)

## Tables of results for contact tracing ----
#' Fraction of presymptomatic cases tracked at each time point before symptom onset
times <- c(-1, -2, -3, -4, -5)

f <- sapply(c("missing","all"), function(type) {
  vals <- ParLapply(times, function(t) {
    CDF <- sapply(indices[[type]], function(i) pgamma(q = t + likelihoods[[i,"shift"]], shape = likelihoods[[i,"shape"]], rate = likelihoods[[i,"rate"]]))
    CDF0 <- sapply(indices[[type]], function(i) pgamma(q = 0 + likelihoods[[i,"shift"]], shape = likelihoods[[i,"shape"]], rate = likelihoods[[i,"rate"]]))

    ii <- mleIndices[[type]]
    CDF.mle <- pgamma(q = t + likelihoods[[ii,"shift"]], shape = likelihoods[[ii,"shape"]], rate = likelihoods[[ii,"rate"]])
    CDF0.mle <- pgamma(q = 0 + likelihoods[[ii,"shift"]], shape = likelihoods[[ii,"shape"]], rate = likelihoods[[ii,"rate"]])

    data.frame(
      t = t,
      mle = (CDF0.mle - CDF.mle) / CDF0.mle,
      mean = mean((CDF0 - CDF) / CDF0),
      lower = min((CDF0 - CDF) / CDF0),
      upper = max((CDF0 - CDF) / CDF0),
      presymp.mle = CDF0.mle,
      presymp.mean = mean(CDF0),
      presymp.lower = min(CDF0),
      presymp.lower = max(CDF0),
      type = factor(type, levels = c("missing", "all"))
    )
  })
  do.call(rbind, vals)
}, simplify = F)
f <- do.call(rbind, f)


# Reconstruct the serial interval ----
#' Sample the infection times from the MLE shifted gamma distibrutions
set.seed(42)
numSamples <- 1e6
infectionTimes <- list(
  missing = rgamma(n = numSamples, shape = likelihoods[mleIndices[["missing"]],"shape"], rate = likelihoods[mleIndices[["missing"]],"rate"]) - likelihoods[mleIndices[["missing"]],"shift"],
  all = rgamma(n = numSamples, shape = likelihoods[mleIndices[["all"]],"shape"], rate = likelihoods[mleIndices[["all"]],"rate"]) - likelihoods[mleIndices[["all"]],"shift"],
  he = rgamma(n = numSamples, shape = 2.1157789506000704, rate = 0.6898582881923863) - 2.3066912325330207
)
#' Sample the incubation times from the log-normal
incubationTimes <- rlnorm(n = numSamples, meanlog = 1.434065, sdlog = 0.6612)
#' Add infection and incubation times to generate serial intervals
serialIntervals <- rbind(
  data.frame(s = infectionTimes$missing + incubationTimes, type = factor("missing", levels = c("missing","all","he"))),
  data.frame(s = infectionTimes$all + incubationTimes, type = factor("all", levels = c("missing","all","he"))),
  data.frame(s = infectionTimes$he + incubationTimes, type = factor("he", levels = c("missing","all","he")))
)
#' Calculate all possible serial intervals from the data
histStep <- 0.5
histBreaks <- seq(-10 + histStep/2, 25, histStep)
intervals <- unlist(sapply(seq_len(nrow(serialData)), function(i) seq(serialData[[i,"minSerialInterval"]] - 0.5, serialData[[i,"maxSerialInterval"]] + 0.5, histStep)))
#' Plot the serial interval data
hist.plot <- ggplot() +
  geom_vline(xintercept = 0, colour = "darkgrey", size = 1) +
  geom_histogram(aes(x = intervals, y = ..density..), breaks = histBreaks) +
  labs(x = "all serial intervals (days)", y = "probability") +
  theme_minimal()
#' Overlay the reconstructed serial interval
hist.plot +
  geom_density(data = serialIntervals[serialIntervals$type %in% c("all"),], aes(x = s, colour = type, size = type, linetype = type), fill = "transparent") +
  geom_density(data = serialIntervals[serialIntervals$type %in% c("missing"),], aes(x = s, colour = type, size = type, linetype = type), fill = "transparent") +
  #scale_alpha_manual(values = c(missing = 0, all = 0), guide = "none") +
  scale_linetype_manual(values = c(missing = "dashed", all = "solid"), labels = c(missing = "He et al.", all = "corrected"), name = NULL) +
  scale_size_manual(values = c(missing = 0.75, all = 1.5), labels = c(missing = "He et al.", all = "corrected"), name = NULL) +
  scale_colour_manual(values = c(missing = "#1b7cff", all = "#ff9e1b", he = "black"), aesthetics = c("colour"), labels = c(missing = "He et al.", all = "corrected"), name = NULL) +
  #scale_x_continuous(limits = c(-10,25)) +
  coord_cartesian(xlim = c(-10,25), ylim = c(0,0.151), expand = F) +
  theme_bw() +
  theme(legend.position = c(0.8,0.8), plot.background = element_rect(fill = "white", colour = "white"), panel.grid = element_blank()) +
  labs(x = "serial interval (days)", y = "probability")

ggsave("serialDistribution.pdf", width = 102, height = 60, units = "mm")
ggsave("serialDistribution.eps", device = cairo_ps, fallback_resolution = 600, width = 102, height = 60, units = "mm")
ggsave("serialDistribution.jpg", width = 102, height = 60, units = "mm", dpi = 600)

## Plot the incubation period ----
incubations <- seq(0,20,0.1)
incubationDist <- data.frame(
  t = incubations,
  pdf = dlnorm(x = incubations, meanlog = 1.434065, sdlog = 0.6612),
  CDF = plnorm(q = incubations, meanlog = 1.434065, sdlog = 0.6612)
)

mean.inc <- exp(1.434065 + (0.6612^2)/2)

ggplot(incubationDist, aes(x = t, y = pdf)) +
  geom_vline(xintercept = mean.inc, linetype = "dashed", colour = "darkgrey", size = 1) +
  annotate("text", label = paste0("mean = ", round(mean.inc,1), " days"),
           x = mean.inc + 0.5, y = 0.15,
           hjust = 0, vjust = 0) +
  geom_line(size = 1.5) +
  labs(x = "incubation period (days)", y = "probability") +
  coord_cartesian(ylim = c(0,0.201), expand = F) +
  theme_bw() + theme(plot.background = element_rect(fill = "white", colour = "white"), panel.grid = element_blank())


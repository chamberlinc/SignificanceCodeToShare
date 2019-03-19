my_PlotIMFs <- function (sig, time.span = NULL, imf.list = NULL, original.signal = TRUE, 
                         residue = TRUE, fit.line = FALSE, lwd = 1, cex = 1, ...) 
{
  opts <- list(...)
  if (!"xlab" %in% names(opts)) {
    opts$xlab <- "Time (s)"
  }
  if (!"ylab" %in% names(opts)) {
    opts$ylab <- ""
  }
  if (is.null(time.span)) {
    time.span = c(min(sig$tt), max(sig$tt))
  }
  if (is.null(imf.list)) {
    imf.list = 1:sig$nimf
  }
  if ("averaged.imfs" %in% names(sig)) {
    sig$imf = sig$averaged.imfs
  }
  if ("averaged.residue" %in% names(sig)) {
    sig$residue = sig$averaged.residue
  }
  if (!"main" %in% names(opts)) {
    opts$main <- "Signal and IMFs"
  }
  time.ind = which(sig$tt >= time.span[1] & sig$tt <= time.span[2])
  tt = sig$tt[time.ind]
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = opts$xlab, 
       ylab = opts$ylab, cex.lab = cex, main = opts$main)
  yax.labs = c()
  snum = length(imf.list) + residue + original.signal
  sp = 1/snum
  if (original.signal) {
    snum = snum + 1
    os = sig$original.signal[time.ind] - mean(sig$original.signal[time.ind])
    scale = max(abs(os))
  }
  else {
    scale = max(abs(sig$imf[time.ind, imf.list]))
  }
  if (residue) {
    snum = snum + 1
    res = sig$residue[time.ind] - mean(sig$residue[time.ind])
    res = res * (sp/(2 * scale))
    yax.labs = append(yax.labs, "Residue")
  }
  trace.pos = sp/2
  imfs = sig$imf * (sp/(scale * 2))
  ts = (tt - min(tt)) * (1/(time.span[2] - time.span[1]))
  if (residue) {
    lines(ts, res + trace.pos, lwd = lwd)
    trace.pos = trace.pos + sp
  }
  for (k in rev(imf.list)) {
    lines(ts, imfs[time.ind, k] + trace.pos, lwd = lwd)
    trace.pos = trace.pos + sp
    yax.labs = append(yax.labs, paste("IMF", k))
  }
  if (original.signal) {
    lines(ts, os * (sp/(2 * scale)) + trace.pos, lwd = lwd)
    yax.labs = append(yax.labs, "Signal")
    if (fit.line) {
      if (length(imf.list) > 1) {
        fline = rowSums(imfs[time.ind, imf.list])
      }
      else {
        fline = imfs[time.ind, imf.list]
      }
      if (residue) {
        fline = fline + res
      }
      lines(ts, fline + trace.pos, lwd = lwd, col = "red")
    }
  }
  xax.labs = pretty(seq(min(tt), max(tt), length.out = 11))
  axis(1, pos = 0, at = seq(0, 1, length.out = length(xax.labs)), 
       labels = xax.labs, cex.axis = cex)
  axis(2, pos = 0, at = seq(sp/2, 1, by = sp), labels = yax.labs, 
       cex.axis = cex)
  segments(c(0, 0, 1, 0), c(0, 1, 1, 0), c(0, 1, 1, 1), c(1, 
                                                          1, 0, 0), lwd = lwd)
}

investigate.IMF <- function(IMFs, named.result, orig.timestamps) {
  library(xts)
  library(dygraphs)
  library(foreach)
  
  buzz <- foreach(i = 1:length(IMFs), .combine = "cbind") %do% {
    imf <- IMFs[i]
    xts(named.result$imf[,imf], order.by = orig.timestamps)
  } 
  dygraph(buzz) %>% dyRangeSelector()
  
}

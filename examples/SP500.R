library(bspcov)
library(dplyr)

# fix random seed
set.seed(1234)

# load and preprocess SP500 data
data("SP500")
sectors <- c("Energy", "Financials", "Materials", "Real Estate", "Utilities", "Information Technology")
SP500_idioerr <- proc_SP500(SP500, sectors)
Uhat <- SP500_idioerr$Uhat
sectornames <- SP500_idioerr$sectornames
PPPres <- thresPPP(Uhat, eps = 0, thres = list(value = 0.0020, fun = "hard"), nsample = 100)

# apply thresPPP
postmean <- estimate(PPPres)

# visualization
diag(postmean) <- 0 # hide color
makeheatmap <- function(mat, colgroups, rowgroups = NULL, keepdiag = T) {
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  if (!keepdiag) diag(mat) <- 0
  if (is.null(rowgroups)) rowgroups <- colgroups

  colnames(mat) <- paste0("C", 1:ncol(mat))
  rownames(mat) <- paste0("R", 1:nrow(mat))

  # Data frame with column annotations.
  annotation_col <- data.frame(sector = as.factor(colgroups))
  rownames(annotation_col) <- colnames(mat)
  annotation_row <- data.frame(sector = as.factor(rowgroups))
  rownames(annotation_row) <- rownames(mat)

  paletteLength <- 50
  myColor <- colorRampPalette(c("red", "black", "green"))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  absmaxval <- max(abs(mat))
  myBreaks <- c(
    seq(-absmaxval, 0, length.out = ceiling(paletteLength / 2) + 1),
    seq(absmaxval / paletteLength, absmaxval, length.out = floor(paletteLength / 2))
  )


  mat_colors <- list(sector = brewer.pal(length(unique(colgroups)), "Set1"))
  names(mat_colors[[1]]) <- unique(colgroups)
  pheatmap(
    mat = mat,
    color = myColor,
    breaks = myBreaks,
    show_colnames = FALSE,
    show_rownames = FALSE,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    drop_levels = TRUE,
    main = " ",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    # fontsize          = 14,
    # border_color      = NA,
    annotation_colors = mat_colors
  )
}

g <- makeheatmap(postmean, sectornames, keepdiag = FALSE)
ggplot2::ggsave("figs/thresPPPheatmap.png", g, width = 12, height = 9)

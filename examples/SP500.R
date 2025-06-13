library(bspcov)
library(dplyr)

# fix random seed
set.seed(123)

# load data
data("SP500_idioerr")

# apply thresPPP
PPPres <- thresPPP(SP500_idioerr$Uhat, eps = 0, thres = list(value = 0.0020, fun = "hard"), nsample = 100)
postmean <- estimate(PPPres)

# visualization
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

g <- makeheatmap(cov2cor(postmean), SP500_idioerr$sectornames, keepdiag = FALSE)
ggsave("figs/thresPPPheatmap.png", g, width = 12, height = 9)

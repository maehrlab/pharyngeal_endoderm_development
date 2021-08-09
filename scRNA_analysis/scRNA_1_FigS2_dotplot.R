# Authors: Jack Huey

library(Seurat)
library(ggplot2)
library(cowplot)

atlas_obj = readRDS("whole_atlas_filtered.Robj")
atlas_obj = UpdateSeuratObject(atlas_obj)

plot.genes = c(
"Epcam",
"Krt8",
"Cldn3",
"Cldn4",
"Krt5",
"Krt15",
"Trp63",
"Six1",
"Pax1",
"Psmb11",
"Myl9",
"Hoxb1",
"Hhex",
"Pax8",
"Foxe1",
"Foxa2",
"Sox2",
    
"Plet1",
"Fgf8",
"Edn1",
"Ripply3",
"Sox9",
"Irf6",
"Vim",
"Irx1",
"Irx2",
"Acta2",
"Tagln"
)

PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

# Modifed from Seurat `DotPlot`
object = atlas_obj
features = plot.genes
assay <- DefaultAssay(object = object)
cols = c("lightgrey", "blue")
col.min = -2.5
col.max = 2.5
dot.min = 0
dot.scale = 6
idents = NULL
group.by = NULL
split.by = NULL
cluster.idents = FALSE
scale = FALSE
scale.by = 'radius'
scale.min = NA
scale.max = NA


split.colors = FALSE
scale.func <- switch(EXPR = scale.by, 'size' = scale_size, 'radius' = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
feature.groups <- NULL
cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

data.features <- FetchData(object = object, vars = features, cells = cells)
data.features$id <- Idents(object = object)[cells, drop = TRUE]
data.features$id <- factor(x = data.features$id)
id.levels <- levels(x = data.features$id)
data.features$id <- as.vector(x = data.features$id)
data.plot <- lapply(
X = unique(x = data.features$id),
FUN = function(ident) {
  data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
  avg.exp <- apply(
    X = data.use,
    MARGIN = 2,
    FUN = function(x) {
      return(mean(x = expm1(x = x)))
    }
  )
  pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
  return(list(avg.exp = avg.exp, pct.exp = pct.exp))
}
)
names(x = data.plot) <- unique(x = data.features$id)
data.plot <- lapply(
X = names(x = data.plot),
FUN = function(x) {
  data.use <- as.data.frame(x = data.plot[[x]])
  data.use$features.plot <- rownames(x = data.use)
  data.use$id <- x
  return(data.use)
}
)
data.plot <- do.call(what = 'rbind', args = data.plot)
#data.plot$id <- factor(x = data.plot$id, levels = id.levels)
data.plot$id <-factor(data.plot$id, c(7, 20, 26, 2, 9, 4, 21, 19, 12, 10, 8, 3, 1, 22, 17, 27, 24, 15, 13, 14, 23, 11, 16, 6, 0, 5, 18, 25))
avg.exp.scaled <- sapply(
X = unique(x = data.plot$features.plot),
FUN = function(x) {
  data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
  if (scale) {
    data.use <- scale(x = data.use)
    data.use <- MinMax(data = data.use, min = col.min, max = col.max)
  } else {
    data.use <- log1p(x = data.use)
  }
  return(data.use)
}
)
avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
data.plot$avg.exp.scaled <- avg.exp.scaled
data.plot$features.plot <- factor(
x = data.plot$features.plot,
levels = features
)
data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
data.plot$pct.exp <- data.plot$pct.exp * 100
color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
plot <- ggplot(data = data.plot, mapping = aes_string(y = 'features.plot', x = 'id')) +
geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
guides(size = guide_legend(title = 'Percent Expressed')) +
labs(
  y = 'Features',
  x = 'Identity'
) +
theme_cowplot()
plot <- plot + scale_color_gradientn(colors = c("gray88", "gray87", "darkorchid4", "mediumvioletred", "violetred3", "tomato2", "darkorange2", "goldenrod1", "lightgoldenrod1", "khaki1"))plot <- plot + guides(color = guide_colorbar(title = 'Log normalized expression'))
plot <- plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
plot


pdf("supp-fig-2-dotplot.pdf", width=10, height=10)
plot
dev.off()

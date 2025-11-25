rmDoublets = function(data.filt) {
  # before doing doublet detection we need to run scaling, variable gene selection and PCA, as well as UMAP for visualization. 
  # These steps will be explored in more detail in coming exercises.
  DefaultAssay(data.filt) = "RNA"
  data.filt = NormalizeData(data.filt)
  data.filt = ScaleData(data.filt)
  data.filt = FindVariableFeatures(data.filt)
  data.filt = RunPCA(data.filt, verbose = F, npcs = 100)

  # Then we run doubletFinder, selecting first 10 PCs and a `pK` value of 0.9. To optimize the parameters, you can run the `paramSweep` function in the package.
  # set seed
  # set seed
  set.seed(8)

  # The expected doublet rate as suggested by 10X
  # The percentage of droplets containing two or more cells is influenced by the quantity of cells loaded into the sequencing machine.
  # ref: https://uofuhealth.utah.edu/huntsman/shared-resources/gcb/htg/single-cell/genomics-10x
  cells.nr = ncol(data.filt)
  nExp.pois = if (cells.nr < 750) {
      round(cells.nr * 0.04)
    } else if (cells.nr < 1500) {
      round(cells.nr * 0.008)
    } else if (cells.nr < 2500) {
      round(cells.nr * 0.016)
    } else if (cells.nr < 3500) {
      round(cells.nr * 0.023)
    } else if (cells.nr < 4500) {
      round(cells.nr * 0.031)
    } else if (cells.nr < 5500) {
      round(cells.nr * 0.039)
    } else if (cells.nr < 6500) {
      round(cells.nr * 0.046)
    } else if (cells.nr < 7500) {
      round(cells.nr * 0.054)
    } else if (cells.nr < 8500) {
      round(cells.nr * 0.061)
    } else if (cells.nr < 9500) {
      round(cells.nr * 0.069)
    } else {
      round(cells.nr * 0.076)
    }
  # nExp.pois

  # pK identification (no ground-truth)
  # pK: neighborhood size, i.e. the number of neighbors (in percentage) to consider when calculating pANN
  sweep.list = invisible(paramSweep(data.filt, PCs = 1:10, num.cores = 1))
  sweep.stats = summarizeSweep(sweep.list)
  bcmvn = find.pK(sweep.stats)

  # Optimal pK are selected by maximizing BCmvn (https://www.sciencedirect.com/science/article/pii/S2405471219300730)
  bcmvn.max = bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk = bcmvn.max$pK
  optimal.pk = as.numeric(levels(optimal.pk))[optimal.pk]
  # optimal.pk

  # run DoubletFinder
  # pN: the proportion of artificial doublet cells to generate (default: 0.25)
  data.filt = doubletFinder(seu = data.filt,
                          PCs = 1:10,
                          pK = optimal.pk,
                          nExp = nExp.pois)
  metadata = data.filt@meta.data
  colnames(metadata)[ncol(metadata)] = "doublet_finder"
  data.filt@meta.data = metadata
  DefaultAssay(data.filt) = "SCT"

  return(data.filt)
}
 

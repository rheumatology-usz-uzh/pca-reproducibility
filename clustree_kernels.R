#!/usr/bin/env Rscript
# clustree_kernels.R
#
# At which clustering resolution do OpenBLAS kernels and thread counts start to change the clusters
# of the Seurat PBMC3k tutorial? The script runs the tutorial once per configuration
# (OpenBLAS kernel x number of threads), clusters the cells at a fine grid of resolutions,
# and compares the clustering trees (clustree) across configurations.
#
# Run from the project directory (so that renv activates), inside the container; see clustree_kernels.sh.
#   Rscript clustree_kernels.R run --node=amd --cores=256   # all configurations on this node, in parallel
#   Rscript clustree_kernels.R compare                      # all nodes found in the output directory
#
# Options (defaults in brackets):
#   --data=DIR          10x directory of the tutorial data          [filtered_gene_bc_matrices/hg19]
#   --out=DIR           output directory                            [clustree_results]
#   --node=LABEL        label for this node                         [amd or intel, from the CPU model]
#   --cores=N           cores to use                                [SLURM_CPUS_PER_TASK, else all]
#   --kernels=LIST      OpenBLAS kernels; native = OpenBLAS chooses [native,SapphireRapids,Cooperlake,SkylakeX,Haswell,Zen,Sandybridge]
#   --threads=LIST      OpenBLAS thread counts                      [1,2,4,8,16]
#   --dims=LIST         numbers of principal components             [10,20,30,50,100]
#   --res=MIN:MAX:STEP  resolution grid                             [0.01:2:0.01]
#   --seed=N            seed for PCA, clustering, and UMAP          [23]
#   --umap=TRUE|FALSE   also store UMAP coordinates                 [TRUE]
#   --reference=ID      reference run for the trees                 [first native run with the fewest threads]
#
# Runs whose .rds result already exists are skipped, so rerunning after a failure only redoes what is missing.
#
# Output of compare, in the output directory:
#   fingerprints.csv   md5 letters per step (normalized data, variable features, PCA, graphs, clusters)
#   runs.csv           kernel actually used, variable features vs the reference, groups of identical clusterings
#   partitions.csv     per number of PCs and resolution: distinct clusterings, adjusted Rand index, cells changed
#   onset.csv          first resolution at which any run differs, for each number of PCs
#   onset.pdf          distinct clusterings and cells changed against resolution
#   clustree_NNNpcs.pdf       clustering tree of the reference, nodes coloured by the share of other runs
#                             that contain exactly the same cluster (grey = reproduced everywhere)
#   clustree_NNNpcs_zoom.pdf  reference and most different run side by side around the first difference
# The worker mode (one configuration) is started by the run mode; it is not meant to be called directly.
# The compare mode needs the clustree package: renv::install("clustree").

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0 && !startsWith(args[1], "--")) args[1] else "run"
opt <- list(data = "filtered_gene_bc_matrices/hg19", out = "clustree_results", node = NA, cores = NA,
            kernels = "native,SapphireRapids,Cooperlake,SkylakeX,Haswell,Zen,Sandybridge",
            threads = "1,2,4,8,16", dims = "10,20,30,50,100", res = "0.01:2:0.01",
            seed = "23", umap = "TRUE", reference = NA,
            id = NA, kernel = NA)                     # set by the run mode for each worker
for (a in args[startsWith(args, "--")]) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (!kv[1] %in% names(opt)) stop("unknown option --", kv[1])
  opt[[kv[1]]] <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"
}
as_list <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])
kernels <- as_list(opt$kernels)
threads <- as.integer(as_list(opt$threads))
dims <- as.integer(as_list(opt$dims))
rg <- as.numeric(strsplit(opt$res, ":", fixed = TRUE)[[1]])
res_grid <- round(seq(rg[1], rg[2], by = rg[3]), 6)
seed <- as.integer(opt$seed)
npcs_for <- function(d) max(50L, d)        # Seurat's default of 50 PCs, more only when needed

cpu_model <- function() {
  x <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)[1]
  if (is.na(x)) NA_character_ else trimws(sub("^[^:]*:", "", x))
}
node <- if (is.na(opt$node)) ifelse(grepl("AMD", cpu_model()), "amd", "intel") else opt$node


# ---------------------------------------------------------------------------------------------------
# run: start one worker process per configuration. OpenBLAS reads its variables when it is loaded,
# so each configuration needs its own R process with the variables set before it starts.
# ---------------------------------------------------------------------------------------------------
run_driver <- function() {
  cores <- if (!is.na(opt$cores)) as.integer(opt$cores) else
    as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores()))
  stopifnot(dir.exists(opt$data))
  out_dir <- file.path(opt$out, node)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.unsetenv(c("OPENBLAS_CORETYPE", "OPENBLAS_NUM_THREADS"))   # workers must not inherit these
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
  cfg <- expand.grid(kernel = kernels, threads = threads, stringsAsFactors = FALSE)
  cfg$id <- sprintf("%s_%s_t%02d", node, cfg$kernel, cfg$threads)
  pass <- c(sprintf("--data=%s", normalizePath(opt$data)), sprintf("--out=%s", opt$out),
            sprintf("--node=%s", node), sprintf("--dims=%s", opt$dims), sprintf("--res=%s", opt$res),
            sprintf("--seed=%s", opt$seed), sprintf("--umap=%s", opt$umap))

  run_one <- function(i) {
    env <- c("OPENBLAS_VERBOSE=2", sprintf("OPENBLAS_NUM_THREADS=%d", cfg$threads[i]),
             "OMP_NUM_THREADS=1", "RCPP_PARALLEL_NUM_THREADS=1",   # vary only the BLAS threads
             if (cfg$kernel[i] != "native") sprintf("OPENBLAS_CORETYPE=%s", cfg$kernel[i]))
    log <- file.path(out_dir, paste0(cfg$id[i], ".log"))
    t0 <- Sys.time()
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c(shQuote(script), "worker", pass, sprintf("--id=%s", cfg$id[i]),
                        sprintf("--kernel=%s", cfg$kernel[i]), sprintf("--threads=%d", cfg$threads[i])),
                      env = env, stdout = log, stderr = log)
    c(status = status, minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  }

  todo <- which(!file.exists(file.path(out_dir, paste0(cfg$id, ".rds"))))   # existing results are kept
  cat(sprintf("node %s (%s): %d configurations, %d to run, on %d cores\n",
              node, cpu_model(), nrow(cfg), length(todo), cores))
  res <- vector("list", nrow(cfg))
  if (sum(cfg$threads[todo]) <= cores) {    # everything fits at once
    res[todo] <- parallel::mclapply(todo, run_one, mc.cores = max(1L, length(todo)), mc.preschedule = FALSE)
  } else {                                  # one wave per thread count, as many runs as fit
    for (t in sort(unique(cfg$threads[todo]), decreasing = TRUE)) {
      idx <- todo[cfg$threads[todo] == t]
      res[idx] <- parallel::mclapply(idx, run_one, mc.cores = max(1L, min(length(idx), cores %/% t)),
                                     mc.preschedule = FALSE)
    }
  }
  cfg$exit <- vapply(res, function(r) if (is.numeric(r)) r[["status"]] else NA_real_, 0)
  cfg$minutes <- vapply(res, function(r) if (is.numeric(r)) r[["minutes"]] else NA_real_, 0)
  cfg$core <- vapply(file.path(out_dir, paste0(cfg$id, ".log")), core_from_log, "")
  cfg$result <- file.exists(file.path(out_dir, paste0(cfg$id, ".rds")))
  print(cfg[, c("id", "core", "exit", "minutes", "result")], row.names = FALSE)
  if (!all(cfg$result)) cat("Some runs failed; see their .log files in", out_dir, "\n")
  cat("Next: run the other node type, then  Rscript clustree_kernels.R compare\n")
}

core_from_log <- function(log) {            # the kernel OpenBLAS actually used (OPENBLAS_VERBOSE=2)
  if (!file.exists(log)) return(NA_character_)
  x <- grep("^Core:", readLines(log, warn = FALSE), value = TRUE)
  if (length(x) == 0) NA_character_ else trimws(sub("^Core:", "", x[1]))
}


# ---------------------------------------------------------------------------------------------------
# worker: the Seurat PBMC3k tutorial for one configuration, clustered at every resolution of the grid
# ---------------------------------------------------------------------------------------------------
run_worker <- function() {
  suppressPackageStartupMessages(library(Seurat))
  set.seed(seed)
  md5 <- function(x) {                      # fingerprint: equal md5 = bit-identical object
    f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, compress = FALSE); unname(tools::md5sum(f))
  }
  normalize_signs <- function(object, reduction) {   # as in seurat.qmd: largest loading positive
    L <- Loadings(object, reduction = reduction); E <- Embeddings(object, reduction = reduction)
    for (i in seq_len(ncol(L))) if (L[which.max(abs(L[, i])), i] < 0) { L[, i] <- -L[, i]; E[, i] <- -E[, i] }
    object[[reduction]]@feature.loadings <- L
    object[[reduction]]@cell.embeddings <- E
    object
  }
  fp <- c(); timing <- c(); t0 <- proc.time()[["elapsed"]]
  tick <- function(step) timing[step] <<- round(proc.time()[["elapsed"]] - t0, 1)

  # the tutorial, step by step (https://satijalab.org/seurat/articles/pbmc3k_tutorial)
  pbmc <- CreateSeuratObject(counts = Read10X(data.dir = opt$data), project = "pbmc3k",
                             min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  fp["normalized"] <- md5(LayerData(pbmc, layer = "data"))
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  fp["variable_features"] <- md5(VariableFeatures(pbmc))
  pbmc <- ScaleData(pbmc, features = rownames(pbmc), verbose = FALSE)
  fp["scaled"] <- md5(LayerData(pbmc, layer = "scale.data"))
  tick("preprocessing")

  for (np in sort(unique(vapply(dims, npcs_for, 0L)))) {        # IRLBA, Seurat's default
    red <- paste0("pca", np)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = np, seed.use = seed,
                   reduction.name = red, reduction.key = paste0("PC", np, "_"), verbose = FALSE)
    pbmc <- normalize_signs(pbmc, red)
    fp[red] <- md5(Embeddings(pbmc, red))
  }
  tick("pca")

  clusters <- list(); umap <- list()
  for (d in dims) {
    red <- paste0("pca", npcs_for(d)); g <- paste0("snn", d)
    pbmc <- FindNeighbors(pbmc, reduction = red, dims = 1:d, graph.name = c(paste0("nn", d), g),
                          verbose = FALSE)
    fp[g] <- md5(pbmc[[g]])
    pbmc <- FindClusters(pbmc, graph.name = g, resolution = res_grid, random.seed = seed, verbose = FALSE)
    cols <- paste0(g, "_res.", res_grid)
    cl <- pbmc[[]][, cols, drop = FALSE]
    cl[] <- lapply(cl, function(x) as.integer(as.character(x)))
    clusters[[as.character(d)]] <- cl
    fp[paste0("clusters", d)] <- md5(cl)
    pbmc@meta.data <- pbmc@meta.data[, setdiff(colnames(pbmc@meta.data), cols), drop = FALSE]
    if (as.logical(opt$umap)) {
      pbmc <- RunUMAP(pbmc, reduction = red, dims = 1:d, seed.use = seed, reduction.name = paste0("umap", d),
                      reduction.key = paste0("UMAP", d, "_"), verbose = FALSE)
      umap[[as.character(d)]] <- Embeddings(pbmc, paste0("umap", d))
    }
    tick(paste0("dims", d))
  }

  saveRDS(list(id = opt$id, node = node, kernel = opt$kernel, threads = as.integer(opt$threads),
               cpu = cpu_model(), blas = La_library(), seed = seed, dims = dims, res = res_grid,
               cells = colnames(pbmc), variable_features = VariableFeatures(pbmc),
               fingerprints = fp, clusters = clusters, umap = umap,
               timing = timing, session = utils::sessionInfo()),
          file.path(opt$out, node, paste0(opt$id, ".rds")))
  cat("done:", opt$id, "\n")
}


# ---------------------------------------------------------------------------------------------------
# compare: where do the clusterings start to differ? Tables and clustering trees.
# ---------------------------------------------------------------------------------------------------
run_compare <- function() {
  for (p in c("clustree", "ggraph", "ggplot2", "patchwork"))
    if (!requireNamespace(p, quietly = TRUE)) stop("package ", p, " is needed: renv::install(\"", p, "\")")
  suppressPackageStartupMessages({ library(ggplot2); library(ggraph) })
  files <- list.files(opt$out, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(files) < 2) stop("need at least two results in ", opt$out)
  runs <- lapply(files, readRDS)
  info <- data.frame(id = vapply(runs, `[[`, "", "id"), node = vapply(runs, `[[`, "", "node"),
                     kernel = vapply(runs, `[[`, "", "kernel"), threads = vapply(runs, `[[`, 0L, "threads"),
                     core = vapply(sub("\\.rds$", ".log", files), core_from_log, ""),
                     cpu = vapply(runs, `[[`, "", "cpu"), stringsAsFactors = FALSE)
  stopifnot(all(vapply(runs, function(r) identical(r$cells, runs[[1]]$cells), TRUE)),
            all(vapply(runs, function(r) identical(r$res, runs[[1]]$res), TRUE)))
  dims <- Reduce(intersect, lapply(runs, `[[`, "dims"))
  res <- runs[[1]]$res

  # reference first: the requested one, or the first native run with the fewest threads
  ref <- if (!is.na(opt$reference)) match(opt$reference, info$id) else {
    nat <- which(info$kernel == "native"); if (length(nat) == 0) nat <- seq_along(runs)
    nat[order(info$node[nat], info$threads[nat])][1]
  }
  if (is.na(ref)) stop("reference ", opt$reference, " not found")
  ord <- unique(c(ref, order(info$node, match(info$kernel, kernels), info$threads)))
  runs <- runs[ord]; info <- info[ord, ]; rownames(info) <- NULL
  vf_ref <- runs[[1]]$variable_features                      # variable features: same set? same order?
  info$hvf_not_in_reference <- vapply(runs, function(r)
    if (is.null(r$variable_features)) NA_integer_ else length(setdiff(r$variable_features, vf_ref)), 0L)
  info$hvf_same_order <- vapply(runs, function(r) identical(r$variable_features, vf_ref), TRUE)
  out <- opt$out
  cat(sprintf("%d runs, reference %s; %d numbers of PCs; %d resolutions from %g to %g\n",
              length(runs), info$id[1], length(dims), length(res), min(res), max(res)))

  # 1. fingerprints: at which step do the runs start to differ? (equal letters = bit-identical)
  stages <- unique(unlist(lapply(runs, function(r) names(r$fingerprints))))
  fpm <- vapply(runs, function(r) unname(r$fingerprints[stages]), character(length(stages)))
  codes <- c(LETTERS, letters)
  fpm <- matrix(fpm, nrow = length(stages))                  # stages x runs, also for one stage
  fpl <- matrix(apply(fpm, 1, function(x) codes[match(x, unique(x))]), nrow = length(runs),
                dimnames = list(NULL, stages))               # runs x stages
  fp_tab <- data.frame(info[, c("id", "core")], fpl, check.names = FALSE)
  write.csv(fp_tab, file.path(out, "fingerprints.csv"), row.names = FALSE)
  cat("\nFingerprints (equal letters = bit-identical):\n"); print(fp_tab, row.names = FALSE)

  # 2. partitions at every number of PCs and resolution
  canon <- function(x) match(x, unique(x))                  # labels in order of first appearance
  ari <- function(a, b) {
    tab <- table(a, b); n2 <- length(a) * (length(a) - 1) / 2
    s <- sum(choose(tab, 2)); sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
    e <- sa * sb / n2; m <- (sa + sb) / 2
    if (m == e) 1 else (s - e) / (m - e)
  }
  changed <- function(a, b) {                                # cells outside a one-to-one matching of clusters
    tab <- unclass(table(a, b)); matched <- 0
    while (length(tab) > 0 && max(tab) > 0) {
      k <- which(tab == max(tab), arr.ind = TRUE)[1, ]
      matched <- matched + tab[k[1], k[2]]; tab <- tab[-k[1], -k[2], drop = FALSE]
    }
    length(a) - matched
  }
  labels <- function(d, r) lapply(runs, function(x) x$clusters[[as.character(d)]][[paste0("snn", d, "_res.", r)]])
  part <- do.call(rbind, lapply(dims, function(d) do.call(rbind, lapply(res, function(r) {
    L <- labels(d, r)
    key <- vapply(L, function(x) paste(canon(x), collapse = ","), "")
    differ <- key[-1] != key[1]
    data.frame(dims = d, resolution = r, partitions = length(unique(key)), runs_differing = sum(differ),
               clusters_reference = length(unique(L[[1]])),
               clusters_range = paste(range(lengths(lapply(L, unique))), collapse = "-"),
               min_ari = if (any(differ)) min(vapply(L[-1][differ], function(x) ari(L[[1]], x), 0)) else 1,
               max_cells_changed = if (any(differ)) max(vapply(L[-1][differ], function(x) changed(L[[1]], x), 0)) else 0)
  }))))
  write.csv(part, file.path(out, "partitions.csv"), row.names = FALSE)

  onset <- do.call(rbind, lapply(dims, function(d) {
    p <- part[part$dims == d, ]; j <- which(p$partitions > 1)
    first <- if (length(j)) p$resolution[j[1]] else NA
    who <- if (length(j)) {
      L <- labels(d, first); key <- vapply(L, function(x) paste(canon(x), collapse = ","), "")
      paste(info$id[-1][key[-1] != key[1]], collapse = " ")
    } else ""
    data.frame(dims = d, first_resolution_differing = first,
               clusters_at_first = if (length(j)) p$clusters_reference[j[1]] else NA,
               resolutions_differing = length(j), resolutions = nrow(p), runs_differing_at_first = who)
  }))
  write.csv(onset, file.path(out, "onset.csv"), row.names = FALSE)
  cat("\nFirst resolution at which the runs differ:\n"); print(onset[, 1:5], row.names = FALSE)

  # groups of runs whose clusterings are identical at every resolution, for each number of PCs
  groups <- vapply(dims, function(d) {
    k <- vapply(runs, function(x) unname(x$fingerprints[paste0("clusters", d)]), "")
    c(LETTERS, letters)[match(k, unique(k))]
  }, character(length(runs)))
  groups <- matrix(groups, nrow = length(runs), dimnames = list(NULL, paste0("PCs_", dims)))
  runs_tab <- data.frame(info, groups, check.names = FALSE)
  write.csv(runs_tab, file.path(out, "runs.csv"), row.names = FALSE)
  cat("\nRuns with equal letters have identical clusterings at every resolution:\n")
  print(runs_tab[, c("id", "core", "hvf_not_in_reference", "hvf_same_order", colnames(groups))], row.names = FALSE)

  # 3. plot: number of distinct clusterings and cells that change cluster, by resolution
  pp <- part; pp$pcs <- factor(paste(pp$dims, "PCs"), levels = paste(dims, "PCs"))
  g1 <- ggplot(pp, aes(resolution, partitions)) + geom_step() + facet_wrap(~pcs, nrow = 1) +
    labs(x = NULL, y = "Distinct clusterings") + theme_bw()
  g2 <- ggplot(pp, aes(resolution, max_cells_changed)) + geom_step() + facet_wrap(~pcs, nrow = 1) +
    labs(x = "Resolution", y = "Cells that change cluster\n(largest over runs)") + theme_bw()
  ggsave(file.path(out, "onset.pdf"), patchwork::wrap_plots(g1, g2, ncol = 1),
         width = 3 + 2.5 * length(dims), height = 6)

  # 4. clustering trees, nodes coloured by the share of other runs that contain exactly the same cluster
  reproduced <- function(i, d, r) {                          # for each cluster of run i at resolution r
    L <- labels(d, r); a <- L[[i]]; others <- L[-i]
    sizes <- table(a)
    hits <- vapply(others, function(b) {
      tab <- table(a, b); sb <- table(b)
      vapply(rownames(tab), function(cl) any(tab[cl, ] == sizes[[cl]] & tab[cl, ] == sb[colnames(tab)]), TRUE)
    }, logical(length(sizes)))
    data.frame(resolution = r, cluster = names(sizes), reproduced = rowMeans(matrix(hits, nrow = length(sizes))))
  }
  tree <- function(i, d, rs, title) {
    prefix <- paste0("snn", d, "_res.")
    df <- runs[[i]]$clusters[[as.character(d)]][, paste0(prefix, rs), drop = FALSE]
    lay <- clustree::clustree(df, prefix = prefix, return = "layout")
    rep <- do.call(rbind, lapply(rs, function(r) reproduced(i, d, r)))
    key <- paste(as.numeric(as.character(lay[[prefix]])), lay$cluster)
    lay$reproduced <- rep$reproduced[match(key, paste(rep$resolution, rep$cluster))]
    yb <- tapply(lay$y, as.numeric(as.character(lay[[prefix]])), mean)
    ggraph(lay) +
      geom_edge_link(aes(colour = count, alpha = in_prop), arrow = arrow(length = unit(4, "points")),
                     end_cap = circle(9, "points"), start_cap = circle(9, "points")) +
      scale_edge_colour_gradientn(colours = viridis::viridis(256), limits = c(0, n_cells), name = "Cells") +
      scale_edge_alpha(limits = c(0, 1), name = "In-proportion") +
      geom_node_point(aes(colour = reproduced, size = size)) +
      scale_colour_gradient(low = "#C0392B", high = "grey80", limits = c(0, 1),
                            name = "Share of other\nruns with the\nsame cluster") +
      scale_size(range = c(4, 12), limits = c(0, n_cells), name = "Cells in cluster") +
      geom_node_text(aes(label = cluster), size = 3) +
      scale_y_continuous(breaks = as.numeric(yb), labels = names(yb)) +
      labs(title = title, y = "Resolution", x = NULL) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank(), axis.text.x = element_blank())
  }
  n_cells <- length(runs[[1]]$cells)                        # common scales, so legends merge
  display <- res[abs(res * 10 - round(res * 10)) < 1e-8]      # multiples of 0.1 for the overview
  if (length(display) < 2) display <- res[unique(round(seq(1, length(res), length.out = min(20, length(res)))))]
  for (d in dims) {
    o <- onset[onset$dims == d, ]
    sub <- if (is.na(o$first_resolution_differing)) "no differences on the grid" else
      sprintf("first difference at resolution %g (%d of %d runs differ)", o$first_resolution_differing,
              part$runs_differing[part$dims == d & part$resolution == o$first_resolution_differing], length(runs) - 1)
    g <- tree(1, d, display, sprintf("%d PCs, reference %s: %s", d, info$id[1], sub))
    ggsave(file.path(out, sprintf("clustree_%03dpcs.pdf", d)), g, width = 12, height = 1.5 + 0.55 * length(display))
    if (!is.na(o$first_resolution_differing)) {               # zoom: five grid steps either side of the onset
      j <- match(o$first_resolution_differing, res); rs <- res[max(1, j - 5):min(length(res), j + 5)]
      L <- labels(d, res[j]); worst <- 1 + which.min(vapply(L[-1], function(x) ari(L[[1]], x), 0))
      gz <- patchwork::wrap_plots(tree(1, d, rs, info$id[1]), tree(worst, d, rs, info$id[worst]), nrow = 1,
                                  guides = "collect") +
        patchwork::plot_annotation(title = sprintf("%d PCs around the first difference (resolution %g)", d, res[j]))
      ggsave(file.path(out, sprintf("clustree_%03dpcs_zoom.pdf", d)), gz, width = 18, height = 1.5 + 0.55 * length(rs))
    }
  }
  cat("\nWritten to", out, ": fingerprints.csv, runs.csv, partitions.csv, onset.csv, onset.pdf, clustree_*.pdf\n")
}


switch(mode,
       run = run_driver(),
       worker = run_worker(),
       compare = run_compare(),
       stop("mode must be run, worker, or compare"))

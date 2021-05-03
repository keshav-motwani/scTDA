get_data = function(cluster = "all", n_genes = 2000) {

  sce = ExperimentHub::ExperimentHub()[["EH2259"]]
  sce = sce[, sce$multiplets == "singlet" & !is.na(sce$cell)]
  SingleCellExperiment::logcounts(sce) = normalize_gene(SingleCellExperiment::counts(sce))
  stim = sce$stim
  patient = sce$ind

  counts = t(SingleCellExperiment::logcounts(sce))
  hvg = Seurat::FindVariableFeatures(t(counts))
  indices = rownames(hvg[order(hvg$vst.variance.standardized, decreasing = TRUE)[1:n_genes], ])
  counts = counts[, indices]

  data = mapply(i = expand.grid(unique(stim), unique(patient))[, 1],
                j = expand.grid(unique(stim), unique(patient))[, 2],
                function(i, j) {
                  result = counts[stim == i & patient == j, ]
                  attr(result, "stim") = i
                  attr(result, "patient") = j
                  attr(result, "cell_type") = sce$cell[stim == i & patient == j]
                  result
                })

  pcs = lapply(data, function(x) {
    y = scale(x, scale = FALSE)
    y %*% irlba::irlba(y, nv = 50)$v
  })

  stim = sapply(data, function(x) attr(x, "stim"))
  patient = sapply(data, function(x) attr(x, "patient"))

  if (cluster != "all") {
    pcs = mapply(x = pcs, y = data, function(x, y) x[attr(y, "cell_type") == cluster, ], SIMPLIFY = FALSE)
    data = lapply(data, function(x) x[attr(x, "cell_type") == cluster, ])
  }

  attr(data, "n_genes") = n_genes
  attr(data, "cluster") = cluster
  attr(pcs, "n_genes") = n_genes
  attr(pcs, "cluster") = cluster

  result = list(data = data,
                pcs = pcs,
                stim = stim,
                patient = patient)

  return(result)

}

normalize_gene = function(data) {

  size_factors = Matrix::colSums(data)
  data = log1p(data %*% Matrix::Diagonal(x = 10000/size_factors))

  return(data)

}


compute_distance_matrices = function(data, distance = c("euclidian", "snn")) {

  if (ncol(data[[1]]) < attr(data, "n_genes")) use_pcs = TRUE else use_pcs = FALSE

  result = lapply(data, function(x) as.matrix(Rfast::Dist(as.matrix(x))))

  attr(result, "cluster") = attr(data, "cluster")
  attr(result, "n_genes") = attr(data, "n_genes")
  attr(result, "use_pcs") = use_pcs
  attr(result, "distance") = distance

  return(result)

}

compute_persistence_diagrams = function(distance_matrices) {

  dir.create("results")
  file = paste0("results/persistence_diagrams_", attr(distance_matrices, "cluster"), "_", attr(distance_matrices, "n_genes"), "_", ifelse(attr(distance_matrices, "use_pcs"), "pca", "gene"), "_", attr(distance_matrices, "distance"), ".rds")
  print(file)

  if (file.exists(file)) {

    persistence_diagrams = readRDS(file)

  } else {

    max_filtration = quantile(do.call(c, distance_matrices), 0.05)
    print(paste0("Max filtration value: ", max_filtration))

    persistence_diagrams = parallel::mclapply(distance_matrices, function(x) {
      pd = TDA::ripsDiag(x, maxdimension = 1, maxscale = max_filtration, dist = "arbitrary")[["diagram"]]
      class(pd) = "matrix"
      print("Done")
      pd
    }, mc.cores = 16)

    attr(persistence_diagrams, "max_filtration") = max_filtration

    saveRDS(persistence_diagrams, file)

  }

  return(persistence_diagrams)

}

compute_persistence_landscapes = function(persistence_diagrams) {

  n_steps = 200
  min_t = 0
  max_t = attr(persistence_diagrams, "max_filtration")
  t_vals = seq(min_t, max_t, (max_t-min_t)/n_steps)

  persistence_landscapes = lapply(persistence_diagrams, function(x) {
    result = matrix(0, nrow = 100, ncol = length(t_vals))
    tda_result = t(TDA::landscape(x, dimension = 1, KK = 1:100, tseq = t_vals))
    result[1:nrow(tda_result), 1:ncol(tda_result)] = tda_result
    result
  })
  attr(persistence_landscapes, "t_vals") = t_vals

  persistence_landscape_matrix = landscape_matrix_from_list(persistence_landscapes)

  return(list(persistence_landscapes = persistence_landscapes,
              persistence_landscape_matrix = persistence_landscape_matrix,
              t_vals = t_vals))

}

landscape_matrix_from_list = function(persistence_landscapes){

  n = length(persistence_landscapes)
  m = ncol(persistence_landscapes[[1]])

  max_depth = integer(n)
  for (i in 1:n) {
    max_depth[i] = nrow(persistence_landscapes[[i]])
  }
  K = max(max_depth)

  persistence_landscape_matrix = Matrix::Matrix(0, nrow = n, ncol = m * K, sparse = TRUE)

  for (i in 1:n) {
    for (j in 1:max_depth[i]) {
      persistence_landscape_matrix[i,(1+(j-1)*m):(j*m)] = persistence_landscapes[[i]][j,]
    }
  }

  return(persistence_landscape_matrix)

}

average_landscapes = function(persistence_landscapes) {

  Reduce("+", persistence_landscapes) / length(persistence_landscapes)

}

plot_landscapes = function(persistence_landscapes, t_vals, nrow = 2) {

  if (!is.list(persistence_landscapes)) {
    persistence_landscapes = list(" " = persistence_landscapes)
  }

  if (!is.list(t_vals)) {
    t_vals = replicate(length(persistence_landscapes), t_vals, simplify = FALSE)
    names(t_vals) = names(persistence_landscapes)
  }

  data_list = list()

  for (name in names(persistence_landscapes)) {

    persistence_landscape = persistence_landscapes[[name]]

    persistence = as.numeric(t(persistence_landscape))
    time = rep(t_vals[[name]], nrow(persistence_landscape))
    k = rep(1:nrow(persistence_landscape), each = ncol(persistence_landscape))
    plot_data = data.frame(persistence = persistence,
                           time = time,
                           k = k,
                           name = name)

    data_list = c(data_list, list(plot_data))

  }

  plot_data = do.call(rbind, data_list)

  plot_data$name = factor(plot_data$name, levels = names(persistence_landscapes))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = persistence, group = k, color = k)) +
    ggplot2::geom_line() +
    ggplot2::theme_classic() +
    ggplot2::scale_color_viridis_c(option = "A") +
    ggplot2::theme(strip.background = ggplot2::element_blank(), strip.placement = "outside") +
    ggplot2::facet_wrap(~name, nrow = nrow, dir = "v")

}

prepare_plotting_data_persistence_diagrams = function(persistence_diagrams, degree) {

  data_list = list()

  for (name in names(persistence_diagrams)) {

    data = as.data.frame(persistence_diagrams[[name]])
    colnames(data) = c("degree", "birth", "death")
    data$name = name
    data$y = 1:nrow(data)
    data$degree = as.factor(data$degree)

    data_list = c(data_list, list(data))

  }

  data = do.call(rbind, data_list)

  data = data[data$degree %in% degree, ]

  data$name = factor(data$name, levels = names(persistence_diagrams))

  return(data)

}

plot_barcodes = function(persistence_diagrams, degree) {

  data = prepare_plotting_data_persistence_diagrams(persistence_diagrams, degree)

  ggplot2::ggplot(data = data) +
    ggplot2::xlab("time") +
    ggplot2::geom_segment(ggplot2::aes_string(x = "birth", y = "y", xend = "death", yend = "y"), size = 0.5, alpha = 1) +
    theme_classic() +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   strip.background = element_blank(), strip.placement = "outside") +
    labs(y = NULL) +
    facet_wrap(name~., nrow = 4, scales = "free", dir = "v")

}

plot_persistence_diagrams = function(persistence_diagrams) {

  data = prepare_plotting_data_persistence_diagrams(persistence_diagrams, c(0, 1))

  ggplot(data, aes(x = birth, y = death, color = degree)) + geom_point(alpha = 0.1, size = 0.5) + theme_classic() + labs(x = "birth", y = "death") +
    ggplot2::theme(strip.background = element_blank(), strip.placement = "outside") +
    facet_wrap(name~., nrow = 4, scales = "fixed", dir = "v") +
    scale_color_manual(values = c("red", "black"))

}

paired_sample_permutation_test = function(feature_matrix, pairing, group) {

  n_pairs = length(unique(pairing))
  n = length(pairing)

  feature_matrix = feature_matrix[order(pairing), ]

  tmp = as.matrix(expand.grid(replicate(n_pairs, c(-1, 1), simplify = FALSE)))
  permutations = matrix(nrow = nrow(tmp), ncol = n)
  permutations[, seq(1, n - 1, 2)] = tmp
  permutations[, -seq(1, n - 1, 2)] = permutations[, seq(1, n - 1, 2)] * -1

  null_dist = apply(permutations, 1, function(x) sum((colMeans(feature_matrix[x == 1, ]) - colMeans(feature_matrix[x == -1, ])) ^ 2))
  obs = sum((colMeans(feature_matrix[group == unique(group)[2], ]) - colMeans(feature_matrix[group == unique(group)[1], ])) ^ 2)

  p_value = mean(null_dist > obs)

  data = data.frame(statistic = null_dist, type = ifelse(null_dist == obs, "obs", "perm"), test = paste0("Paired (p-value = ", p_value, ")"))

  return(list(p_value = p_value,
              null_dist = null_dist,
              obs = obs,
              data = data))

}

two_sample_permutation_test = function(feature_matrix, groups, n_repeats = 2^7) {

  n1 = sum(groups == unique(groups)[1])

  null_dist = sapply(1:n_repeats, function(x) {
    order = sample(1:length(groups), length(groups))
    sum((colMeans(feature_matrix[order[1:n1], ]) - colMeans(feature_matrix[order[(n1 + 1):length(groups)], ])) ^ 2)
  })
  obs = sum((colMeans(feature_matrix[groups == unique(groups)[2], ]) - colMeans(feature_matrix[groups == unique(groups)[1], ])) ^ 2)

  p_value = mean(null_dist > obs)

  data = data.frame(statistic = c(obs, null_dist), type = ifelse(c(obs, null_dist) == obs, "obs", "perm"), test = paste0("Two sample (p-value = ", p_value, ")"))

  return(list(p_value = p_value,
              null_dist = null_dist,
              obs = obs,
              data = data))

}

plot_permutation_tests = function(permutation_tests, ncol = 2) {

  data_list = list()

  for (name in names(permutation_tests)) {

    data = permutation_tests[[name]]
    data$name = name

    data_list = c(data_list, list(data))

  }

  data = do.call(rbind, data_list)

  ggplot(data, aes(x = statistic)) +
    geom_histogram(bins = 10, color = "black", fill = "white") +
    theme_classic() +
    geom_rug(aes(color = type), alpha = 0.5) +
    scale_color_manual(values = c("red", "black")) +
    geom_vline(data = data[data$type == "obs", ], aes(xintercept = statistic), color = "red") +
    theme(legend.position = "none") +
    facet_wrap(~name + test, ncol = ncol, scales = "free") +
    theme(strip.background = element_blank(), strip.placement = "outside")

}

plot_persistence_landscape_pca = function(landscape_matrices, pairing, group) {

  data_list = list()

  for (name in names(landscape_matrices)) {

    mat = landscape_matrices[[name]]
    mat = scale(mat, scale = FALSE)
    pc = mat %*% svd(mat)$v[, 1:2]
    data = data.frame(pc)
    colnames(data) = paste0("PC", 1:2)
    data$patient = as.factor(pairing)
    data$treatment = group
    data$name = name

    data_list = c(data_list, list(data))

  }

  data = do.call(rbind, data_list)

  ggplot(data, aes(x = PC1, y = PC2, color = treatment, group = patient, shape = treatment)) +
    geom_point(size = 3) +
    geom_line(color = "gray") +
    theme_classic() +
    scale_color_manual(values = c("black", "red")) +
    facet_wrap(~name, nrow = 2, scales = "free") +
    theme(strip.background = element_blank(), strip.placement = "outside")

}

run_analysis = function(clusters) {

  difference_landscapes = list()
  landscape_matrices = list()
  permutation_test_data = list()
  t_vals_list = list()

  for (cluster in clusters) {

    data = get_data(cluster = cluster)
    stim = data$stim
    patient = data$patient
    names = paste(stim, patient)

    distance_matrices = compute_distance_matrices(data$pcs, distance = "euclidian")
    persistence_diagrams = compute_persistence_diagrams(distance_matrices)
    names(persistence_diagrams) = names

    persistence_landscapes = compute_persistence_landscapes(persistence_diagrams)
    names(persistence_landscapes$persistence_landscapes) = paste(stim, patient)

    difference = average_landscapes(persistence_landscapes$persistence_landscapes[stim == "stim"]) - average_landscapes(persistence_landscapes$persistence_landscapes[stim == "ctrl"])
    t_vals = persistence_landscapes$t_vals

    plot_persistence_diagrams(persistence_diagrams)
    ggsave(file.path("results", paste0(cluster, "_persistence_diagrams.pdf")), height = 6, width = 6)
    plot_barcodes(persistence_diagrams, 1)
    ggsave(file.path("results", paste0(cluster, "_barcodes.pdf")), height = 5, width = 5)

    plot_landscapes(persistence_landscapes$persistence_landscapes, persistence_landscapes$t_vals, nrow = 4)
    ggsave(file.path("results", paste0(cluster, "_persistence_landscapes.pdf")), height = 6, width = 9)
    plot_landscapes(difference, persistence_landscapes$t_vals)
    ggsave(file.path("results", paste0(cluster, "_persistence_landscape_difference.pdf")), height = 4, width = 6)
    plot_persistence_landscape_pca(list(" " = persistence_landscapes$persistence_landscape_matrix), patient, stim)
    ggsave(file.path("results", paste0(cluster, "_persistence_landscape_pca.pdf")), height = 4, width = 5)

    paired = paired_sample_permutation_test(persistence_landscapes$persistence_landscape_matrix, patient, stim)
    two_sample = two_sample_permutation_test(persistence_landscapes$persistence_landscape_matrix, stim)
    data = rbind(paired$data, two_sample$data)
    plot_permutation_tests(list(" " = data))
    ggsave(file.path("results", paste0(cluster, "_permutation_test.pdf")), height = 4, width = 8)

    difference_landscapes = c(difference_landscapes, list(difference))
    landscape_matrices = c(landscape_matrices, list(persistence_landscapes$persistence_landscape_matrix))
    permutation_test_data = c(permutation_test_data, list(data))
    t_vals_list = c(t_vals_list, list(t_vals))

  }

  names(difference_landscapes) = clusters
  names(permutation_test_data) = clusters
  names(landscape_matrices) = clusters
  names(t_vals_list) = clusters

  plot_permutation_tests(permutation_test_data) +
    facet_wrap(~name + test, nrow = 2, dir = "v", scales = "free")
  ggsave("results/combined_permutation_test.pdf", height = 6, width = 14)

  plot_landscapes(difference_landscapes, t_vals_list)
  ggsave("results/combined_persistence_landscape_difference.pdf", height = 6, width = 9)

  plot_persistence_landscape_pca(landscape_matrices, patient, stim)
  ggsave("results/combined_persistence_landscape_pca.pdf", height = 6, width = 9)

}

library(ggplot2)
run_analysis("all")
run_analysis(c("B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells"))

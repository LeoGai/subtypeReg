if (getRversion() >= "2.15.1") utils::globalVariables(c("State","Cluster","ID","OptimalState"))

.valid_k_range <- function(k_range, n_id) {
  k_range <- unique(as.integer(k_range))
  k_range <- k_range[!is.na(k_range) & k_range >= 2]
  if (length(k_range) == 0L) k_range <- 2L
  k_max <- max(2L, n_id - 1L)
  k_range[k_range <= k_max]
}

#' Subtype-aware Registration of Longitudinal Data
#'
#' @param data data.frame with columns ID, Feature, Time, Value
#' @param alpha in (0,1], quantile for selecting representative points per cluster
#' @param k_range integer vector, candidate numbers of clusters, e.g., 2:8
#' @param timepos_option numeric shifts, e.g., c(4,3,2,1)
#' @param tau silhouette threshold; if initial > tau, do fewer iterations
#' @param tmin,tmax numeric, time window
#' @param knots numeric, internal knots for splines::bs
#' @param degree integer, spline degree, default 3
#' @param lambda numeric, glmnet ridge lambda, default 1e-7
#' @param verbose logical, print silhouettes and centers if TRUE
#'
#' @return data.frame with registered Time within [tmin, tmax]
#' @import dplyr
#' @importFrom cluster pam silhouette
#' @importFrom stats setNames
#' @export
SubtypeAware_Registration <- function(
  data,
  alpha = 0.95,
  k_range = 2:8,
  timepos_option = c(4,3,2,1),
  tau = 0.45,
  tmin,
  tmax,
  knots,
  degree = 3,
  lambda = 1e-7,
  verbose = FALSE
) {
  stopifnot(all(c("ID","Time","Value") %in% names(data)))
  ids <- unique(data$ID)
  n_id <- length(ids)
  k_range <- .valid_k_range(k_range, n_id)

  pick_best_state <- function(center, df_all_states_one_id, timepos_option, num_cols) {
    cols <- num_cols
    lv <- c("Original", paste("Shift", timepos_option))
    df <- df_all_states_one_id
    df$State <- factor(df$State, levels = lv, ordered = TRUE)
    df <- df[order(df$State), , drop = FALSE]

    m  <- as.matrix(df[, cols, drop = FALSE])
    d2 <- rowSums((m - matrix(center, nrow = nrow(m), ncol = length(center), byrow = TRUE))^2)
    d2[!is.finite(d2)] <- NA_real_
    if (all(is.na(d2))) return(c(state = "Original", mind = NA_real_))
    idx <- which.min(d2)
    c(state = as.character(df$State[idx]), mind = d2[idx])
  }

  get_coef_one <- function(df_id, shift = NULL) {
    tmp <- df_id
    if (!is.null(shift)) tmp$Time <- tmp$Time + shift
    tmp <- tmp[tmp$Time >= tmin & tmp$Time <= tmax, , drop = FALSE]
    calculate_bspline_coefficients(tmp, knots = knots, degree = degree, lambda = lambda)
  }

  cf_template <- calculate_bspline_coefficients(
    data[data$ID == ids[1], , drop = FALSE],
    knots = knots, degree = degree, lambda = lambda
  )
  coef_names <- names(cf_template)
  num_cols <- coef_names

  rows <- list()
  for (id in ids) {
    df_id <- data[data$ID == id, , drop = FALSE]
    cf0 <- get_coef_one(df_id, shift = NULL)
    rows[[length(rows) + 1L]] <- c(ID = id, State = "Original", stats::setNames(cf0, coef_names))
    for (s in timepos_option) {
      cfs <- get_coef_one(df_id, shift = s)
      rows[[length(rows) + 1L]] <- c(ID = id, State = paste("Shift", s), stats::setNames(cfs, coef_names))
    }
  }
  all_coef <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  all_coef$ID <- if (is.numeric(all_coef$ID)) as.numeric(all_coef$ID) else all_coef$ID
  all_coef[num_cols] <- lapply(all_coef[num_cols], function(x) as.numeric(x))

  base_coef <- dplyr::filter(all_coef, State == "Original")
  data_clustering <- base_coef[, num_cols, drop = FALSE]

  sil_vec <- numeric(length(k_range))
  for (i in seq_along(k_range)) {
    k <- k_range[i]
    km <- cluster::pam(data_clustering, k = k)
    sil <- cluster::silhouette(km$clustering, stats::dist(data_clustering))
    sil_vec[i] <- mean(sil[, "sil_width"])
  }
  if (verbose) print(sil_vec)

  optimal_k <- k_range[which.max(sil_vec)]
  km0 <- cluster::pam(data_clustering, k = optimal_k)
  base_coef$Cluster <- km0$clustering

  selected_list <- list()
  for (k in seq_len(optimal_k)) {
    cd <- dplyr::filter(base_coef, Cluster == k)
    dists <- sqrt(rowSums((cd[, num_cols, drop = FALSE] - km0$medoids[k, ])^2))
    thr <- stats::quantile(dists, alpha, na.rm = TRUE)
    selected_list[[k]] <- cd[dists <= thr, , drop = FALSE]
  }

  optimal_states <- data.frame(ID = base_coef$ID[0], Cluster = integer(0),
                               OptimalState = character(0), MinDistance = numeric(0))
  for (k in seq_len(optimal_k)) {
    repeat {
      any_update <- FALSE
      sel <- selected_list[[k]]
      if (!nrow(sel)) { selected_list[[k]] <- sel; break }
      center <- colMeans(sel[, num_cols, drop = FALSE], na.rm = TRUE)
      if (any(!is.finite(center))) { selected_list[[k]] <- sel; break }
      for (i in seq_len(nrow(sel))) {
        pid <- sel$ID[i]
        all_states_pid <- dplyr::filter(all_coef, ID == pid)
        pick <- pick_best_state(center, all_states_pid, timepos_option, num_cols)
        if (sel$State[i] != pick["state"]) { sel$State[i] <- pick["state"]; any_update <- TRUE }
        optimal_states <- rbind(optimal_states,
                                data.frame(ID = pid, Cluster = k,
                                           OptimalState = as.character(pick["state"]),
                                           MinDistance = as.numeric(pick["mind"])))
      }
      selected_list[[k]] <- sel
      if (!isTRUE(any_update)) break
    }
  }
  optimal_states <- dplyr::distinct(optimal_states, ID, .keep_all = TRUE)

  apply_best_back <- function(df_base, df_all, opt_map) {
    merged <- dplyr::inner_join(df_all, opt_map[, c("ID","OptimalState")], by = "ID")
    keep <- dplyr::filter(merged, State == OptimalState)
    res <- df_base
    row_pos <- match(keep$ID, res$ID)
    res[row_pos, c(num_cols, "State")] <- keep[, c(num_cols, "State")]
    res
  }
  data_clustering_adj <- apply_best_back(base_coef, all_coef, optimal_states)

  iter_count     <- 0
  previous_k     <- optimal_k
  s_old          <- sil_vec
  max_iterations <- if (max(sil_vec, na.rm = TRUE) > tau) 1L else 100L

  repeat {
    iter_count <- iter_count + 1
    sil_vec2 <- numeric(length(k_range))
    for (i in seq_along(k_range)) {
      k <- k_range[i]
      km <- cluster::pam(data_clustering_adj[, num_cols, drop = FALSE], k = k)
      sil <- cluster::silhouette(km$clustering, stats::dist(data_clustering_adj[, num_cols, drop = FALSE]))
      sil_vec2[i] <- mean(sil[, "sil_width"])
    }
    if (verbose) print(sil_vec2)

    k2 <- k_range[which.max(sil_vec2)]
    top2_new <- sort(sil_vec2, decreasing = TRUE)[1:2]
    top2_old <- sort(s_old,     decreasing = TRUE)[1:2]
    if (k2 == previous_k && all(top2_new <= top2_old)) break
    if (iter_count >= max_iterations) break
    if (iter_count > 50) break

    optimal_k  <- k2
    previous_k <- k2
    s_old      <- sil_vec2

    km2 <- cluster::pam(data_clustering_adj[, num_cols, drop = FALSE], k = optimal_k)
    data_clustering_adj$Cluster <- km2$clustering

    selected_list <- list()
    for (k in seq_len(optimal_k)) {
      cd <- dplyr::filter(data_clustering_adj, Cluster == k)
      dists <- sqrt(rowSums((cd[, num_cols, drop = FALSE] - km2$medoids[k, ])^2))
      thr <- stats::quantile(dists, alpha, na.rm = TRUE)
      selected_list[[k]] <- cd[dists <= thr, , drop = FALSE]
    }

    optimal_states <- optimal_states[0, ]
    for (k in seq_len(optimal_k)) {
      repeat {
        any_update <- FALSE
        sel <- selected_list[[k]]
        if (!nrow(sel)) { selected_list[[k]] <- sel; break }
        center <- colMeans(sel[, num_cols, drop = FALSE], na.rm = TRUE)
        if (any(!is.finite(center))) { selected_list[[k]] <- sel; break }
        for (i in seq_len(nrow(sel))) {
          pid <- sel$ID[i]
          all_states_pid <- dplyr::filter(all_coef, ID == pid)
          pick <- pick_best_state(center, all_states_pid, timepos_option, num_cols)
          if (sel$State[i] != pick["state"]) { sel$State[i] <- pick["state"]; any_update <- TRUE }
          optimal_states <- rbind(optimal_states,
                                  data.frame(ID = pid, Cluster = k,
                                             OptimalState = as.character(pick["state"]),
                                             MinDistance = as.numeric(pick["mind"])))
        }
        selected_list[[k]] <- sel
        if (!isTRUE(any_update)) break
      }
    }
    optimal_states <- dplyr::distinct(optimal_states, ID, .keep_all = TRUE)
    data_clustering_adj <- apply_best_back(data_clustering_adj, all_coef, optimal_states)
  }

  if (nrow(optimal_states) > 0) {
    centers_mean <- data_clustering_adj %>%
      dplyr::filter(ID %in% optimal_states$ID) %>%
      dplyr::group_by(Cluster) %>%
      dplyr::summarise(dplyr::across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  } else {
    centers_mean <- data_clustering_adj %>%
      dplyr::group_by(Cluster) %>%
      dplyr::summarise(dplyr::across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  }

  if (isTRUE(verbose)) {
    cat("[SubtypeAware_Registration] final k =", optimal_k, "\n")
    cat("[SubtypeAware_Registration] Centers (mean):\n")
    print(round(centers_mean[, c("Cluster", num_cols), drop = FALSE], 4))
  }

  remain_ids <- setdiff(data_clustering_adj$ID, optimal_states$ID)
  if (length(remain_ids)) {
    for (pid in remain_ids) {
      all_states_pid <- dplyr::filter(all_coef, ID == pid)[, c("State", num_cols), drop = FALSE]
      ord_num <- suppressWarnings(as.numeric(gsub("^Shift\\s*", "", all_states_pid$State)))
      ord     <- ifelse(all_states_pid$State == "Original", -Inf, ord_num)
      all_states_pid <- all_states_pid[order(ord), , drop = FALSE]

      best <- NULL
      for (cidx in seq_len(nrow(centers_mean))) {
        center_vec <- as.numeric(centers_mean[cidx, num_cols, drop = TRUE])
        if (any(!is.finite(center_vec))) next

        M  <- as.matrix(all_states_pid[, num_cols, drop = FALSE])
        d2 <- rowSums((M - matrix(center_vec, nrow = nrow(M), ncol = length(center_vec), byrow = TRUE))^2)

        if (all(!is.finite(d2))) next
        j <- which.min(replace(d2, !is.finite(d2), Inf))

        rec <- data.frame(
          ID = pid,
          Cluster = centers_mean$Cluster[cidx],
          OptimalState = as.character(all_states_pid$State[j]),
          MinDistance = as.numeric(d2[j])
        )
        if (is.null(best) || rec$MinDistance < best$MinDistance) best <- rec
      }
      if (is.null(best)) {
        best <- data.frame(ID = pid, Cluster = NA_integer_,
                           OptimalState = "Original", MinDistance = NA_real_)
      }
      optimal_states <- rbind(optimal_states, best)
    }
  }

  parse_shift <- function(state) {
    if (identical(state, "Original")) return(0)
    as.numeric(trimws(sub("^Shift\\s*", "", state)))
  }
  shift_map <- setNames(vapply(optimal_states$OptimalState, parse_shift, numeric(1)), optimal_states$ID)

  out <- data
  for (id in names(shift_map)) {
    sel <- out$ID == id
    out$Time[sel] <- out$Time[sel] + shift_map[[id]]
  }
  out <- out[out$Time >= tmin & out$Time <= tmax, , drop = FALSE]
  rownames(out) <- NULL

  out
}

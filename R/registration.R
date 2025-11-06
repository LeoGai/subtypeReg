#' Subtype-aware temporal registration via spline-coefficient clustering
#'
#' @param data data.frame with columns ID, Feature, Time, Value
#' @param alpha in (0,1], quantile for selecting representative points per cluster
#' @param k_range integer vector, candidate numbers of clusters, e.g., 2:8
#' @param timepos_option numeric shifts (can be negative), e.g., c(-1,0,1) or c(4,3,2,1)
#' @param tau silhouette threshold; if initial > tau, do fewer iterations
#' @param tmin,tmax numeric, time window
#' @param knots numeric, internal knots for splines::bs
#' @param degree integer, spline degree, default 3
#' @param lambda numeric, glmnet ridge lambda, default 1e-7
#'
#' @return data.frame with registered Time within [tmin, tmax]
#' @import dplyr
#' @importFrom cluster pam silhouette
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
  lambda = 1e-7
) {
  stopifnot(all(c("ID","Time","Value") %in% names(data)))
  ids <- unique(data$ID)
  n_id <- length(ids)
  k_range <- .valid_k_range(k_range, n_id)

  get_coef_one <- function(df_id, shift = NULL) {
    tmp <- df_id
    if (!is.null(shift)) {
      tmp$Time <- tmp$Time + shift
      tmp <- tmp[tmp$Time >= tmin & tmp$Time <= tmax, , drop = FALSE]
    }
    calculate_bspline_coefficients(tmp, knots = knots, degree = degree, lambda = lambda)
  }

  cf_template <- calculate_bspline_coefficients(
    data[data$ID == ids[1], , drop = FALSE],
    knots = knots, degree = degree, lambda = lambda
  )
  coef_names <- names(cf_template)

  rows <- list()
  for (id in ids) {
    df_id <- data[data$ID == id, , drop = FALSE]
    cf0 <- get_coef_one(df_id, shift = NULL)
    rows[[length(rows) + 1L]] <- c(ID = id, State = "Original", setNames(cf0, coef_names))
    for (s in timepos_option) {
      cfs <- get_coef_one(df_id, shift = s)
      rows[[length(rows) + 1L]] <- c(ID = id, State = paste0("Shift ", s), setNames(cfs, coef_names))
    }
  }
  all_coef <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  all_coef$ID <- if (is.numeric(all_coef$ID)) as.numeric(all_coef$ID) else all_coef$ID
  num_cols <- coef_names
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

  pick_best_state <- function(center, df_all_states_one_id) {
    m <- as.matrix(df_all_states_one_id[, num_cols, drop = FALSE])
    d2 <- rowSums((m - matrix(center, nrow = nrow(m), ncol = length(center), byrow = TRUE))^2)
    idx <- which.min(d2)
    c(state = df_all_states_one_id$State[idx], mind = d2[idx])
  }

  optimal_states <- data.frame(ID = base_coef$ID[0], Cluster = integer(0), OptimalState = character(0), MinDistance = numeric(0))
  for (k in seq_len(optimal_k)) {
    repeat {
      any_update <- FALSE
      sel <- selected_list[[k]]
      center <- colMeans(sel[, num_cols, drop = FALSE], na.rm = TRUE)
      for (i in seq_len(nrow(sel))) {
        pid <- sel$ID[i]
        all_states_pid <- dplyr::filter(all_coef, ID == pid)
        pick <- pick_best_state(center, all_states_pid)
        if (sel$State[i] != pick["state"]) {
          sel$State[i] <- pick["state"]; any_update <- TRUE
        }
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

  iter_count <- 0
  if (max(sil_vec, na.rm = TRUE) <= tau) {
    repeat {
      iter_count <- iter_count + 1
      sil_vec2 <- numeric(length(k_range))
      for (i in seq_along(k_range)) {
        k <- k_range[i]
        km <- cluster::pam(data_clustering_adj[, num_cols, drop = FALSE], k = k)
        sil <- cluster::silhouette(km$clustering, stats::dist(data_clustering_adj[, num_cols, drop = FALSE]))
        sil_vec2[i] <- mean(sil[, "sil_width"])
      }
      k2 <- k_range[which.max(sil_vec2)]
      if (k2 == optimal_k || iter_count > 50) break
      optimal_k <- k2
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
          center <- colMeans(sel[, num_cols, drop = FALSE], na.rm = TRUE)
          for (i in seq_len(nrow(sel))) {
            pid <- sel$ID[i]
            all_states_pid <- dplyr::filter(all_coef, ID == pid)
            pick <- pick_best_state(center, all_states_pid)
            if (sel$State[i] != pick["state"]) {
              sel$State[i] <- pick["state"]; any_update <- TRUE
            }
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
  }

  km_final <- cluster::pam(data_clustering_adj[, num_cols, drop = FALSE], k = optimal_k)
  data_clustering_adj$Cluster <- km_final$clustering
  centers <- as.data.frame(km_final$medoids)
  centers$Cluster <- seq_len(nrow(centers))

  remain_ids <- setdiff(data_clustering_adj$ID, optimal_states$ID)
  if (length(remain_ids)) {
    for (pid in remain_ids) {
      all_states_pid <- dplyr::filter(all_coef, ID == pid)
      best <- None <- NULL
      for (cidx in seq_len(nrow(centers))) {
        center <- as.numeric(centers[cidx, colnames(km_final$medoids), drop = TRUE])
        pick <- pick_best_state(center, all_states_pid)
        rec <- data.frame(ID = pid, Cluster = cidx,
                          OptimalState = as.character(pick["state"]),
                          MinDistance = as.numeric(pick["mind"]))
        if (is.null(best) || rec$MinDistance < best$MinDistance) best <- rec
      }
      optimal_states <- rbind(optimal_states, best)
    }
  }

  parse_shift <- function(state) {
    if (identical(state, "Original")) return(0)
    as.numeric(sub("^Shift\s+", "", state))
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

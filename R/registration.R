if (getRversion() >= "2.15.1") utils::globalVariables(c("State","Cluster","ID","OptimalState"))

.valid_k_range <- function(k_range, n_id) {
  k_range <- unique(as.integer(k_range))
  k_range <- k_range[!is.na(k_range) & k_range >= 2]
  if (length(k_range) == 0L) k_range <- 2L
  k_max <- max(2L, n_id - 1L)
  k_range[k_range <= k_max]
}

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
  timepos_option = c(4, 3, 2, 1),
  tau = 0.45,
  tmin,
  tmax,
  knots,
  degree = 3,
  lambda = 1e-7,
  verbose = FALSE
) {
  stopifnot(all(c("ID","Time","Value") %in% names(data)))

  # use character IDs internally to avoid join mismatches
  ids <- unique(as.character(data$ID))

  get_coef_one <- function(df_id, shift = NULL) {
    tmp <- df_id
    if (!is.null(shift)) tmp$Time <- tmp$Time + shift
    tmp <- tmp[tmp$Time >= tmin & tmp$Time <= tmax, , drop = FALSE]
    calculate_bspline_coefficients(
      tmp, knots = knots, degree = degree, lambda = lambda, standardize = TRUE
    )
  }

  # build all states: Original, then Shift <timepos_option>
  cf_template <- get_coef_one(data[as.character(data$ID) == ids[1], , drop = FALSE], shift = NULL)
  coef_names <- names(cf_template)
  num_cols <- coef_names

  rows <- list()
  for (id in ids) {
    df_id <- data[as.character(data$ID) == id, , drop = FALSE]
    cf0 <- get_coef_one(df_id, shift = NULL)
    rows[[length(rows) + 1L]] <- c(ID = id, State = "Original",
                                   stats::setNames(cf0, coef_names))
    for (s in timepos_option) {
      cfs <- get_coef_one(df_id, shift = s)
      rows[[length(rows) + 1L]] <- c(ID = id, State = paste("Shift", s),
                                     stats::setNames(cfs, coef_names))
    }
  }
  all_coef <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  all_coef$ID <- as.character(all_coef$ID)
  all_coef[num_cols] <- lapply(all_coef[num_cols], function(x) as.numeric(x))

  base_coef <- all_coef[all_coef$State == "Original", , drop = FALSE]
  data_clustering <- base_coef[, num_cols, drop = FALSE]

  # choose k via silhouette
  kr <- unique(as.integer(k_range)); kr <- kr[!is.na(kr) & kr >= 2]
  if (!length(kr)) kr <- 2L
  sil_vec <- numeric(length(kr))
  for (i in seq_along(kr)) {
    k <- kr[i]
    km <- cluster::pam(data_clustering, k = k)
    sil <- cluster::silhouette(km$clustering, stats::dist(data_clustering))
    sil_vec[i] <- mean(sil[, "sil_width"])
  }
  if (verbose) print(sil_vec)

  optimal_k <- kr[which.max(sil_vec)]
  km0 <- cluster::pam(data_clustering, k = optimal_k)
  base_coef$Cluster <- km0$clustering

  # representatives by alpha-quantile, na.rm = TRUE
  selected_list <- vector("list", optimal_k)
  for (k in seq_len(optimal_k)) {
    cd <- base_coef[base_coef$Cluster == k, , drop = FALSE]
    d <- sqrt(rowSums((cd[, num_cols, drop = FALSE] - km0$medoids[k, ])^2))
    thr <- stats::quantile(d, alpha, na.rm = TRUE)
    selected_list[[k]] <- cd[d <= thr, , drop = FALSE]
  }

  # tie-break order
  pick_best_state <- function(center, df_states) {
    lv <- c("Original", paste("Shift", timepos_option))
    df <- df_states
    df$State <- factor(df$State, levels = lv, ordered = TRUE)
    df <- df[order(df$State), , drop = FALSE]
    M <- as.matrix(df[, num_cols, drop = FALSE])
    d2 <- rowSums((M - matrix(center, nrow = nrow(M), ncol = length(center), byrow = TRUE))^2)
    idx <- which.min(d2)
    c(state = as.character(df$State[idx]), mind = d2[idx])
  }

  # assign best states for representatives
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
        pid <- as.character(sel$ID[i])
        all_states_pid <- all_coef[all_coef$ID == pid, , drop = FALSE]
        pick <- pick_best_state(center, all_states_pid)
        if (sel$State[i] != pick["state"]) { sel$State[i] <- pick["state"]; any_update <- TRUE }
        optimal_states <- rbind(optimal_states,
                                data.frame(ID = pid, Cluster = k,
                                           OptimalState = as.character(pick["state"]),
                                           MinDistance = as.numeric(pick["mind"])))
      }
      selected_list[[k]] <- sel
      if (!is.TRUE(any_update)) break
    }
  }
  optimal_states <- dplyr::distinct(optimal_states, ID, .keep_all = TRUE)
  optimal_states$ID <- as.character(optimal_states$ID)

  # apply chosen states back to base_coef
  merged <- dplyr::inner_join(all_coef, optimal_states[, c("ID","OptimalState")], by = "ID")
  keep <- merged[merged$State == merged$OptimalState, , drop = FALSE]
  data_clustering_adj <- base_coef
  row_pos <- match(keep$ID, data_clustering_adj$ID)
  data_clustering_adj[row_pos, c(num_cols, "State")] <- keep[, c(num_cols, "State")]

  # legacy stopping rule loop (no final PAM)
  iter_count <- 0L
  previous_k <- optimal_k
  s_old <- sil_vec
  max_iterations <- if (max(sil_vec, na.rm = TRUE) > tau) 1L else 100L

  repeat {
    iter_count <- iter_count + 1L
    sil_vec2 <- numeric(length(kr))
    for (i in seq_along(kr)) {
      k <- kr[i]
      km <- cluster::pam(data_clustering_adj[, num_cols, drop = FALSE], k = k)
      sil <- cluster::silhouette(km$clustering, stats::dist(data_clustering_adj[, num_cols, drop = FALSE]))
      sil_vec2[i] <- mean(sil[, "sil_width"])
    }
    if (verbose) print(sil_vec2)
    k2 <- kr[which.max(sil_vec2)]

    if (k2 == previous_k &&
        all(sort(sil_vec2, decreasing = TRUE)[1:2] <= sort(s_old, decreasing = TRUE)[1:2])) {
      break
    }
    if (max_iterations == 1L) break
    if (iter_count > 50L) break

    previous_k <- k2
    s_old <- sil_vec2

    km2 <- cluster::pam(data_clustering_adj[, num_cols, drop = FALSE], k = k2)
    data_clustering_adj$Cluster <- km2$clustering

    selected_list <- vector("list", k2)
    for (k in seq_len(k2)) {
      cd <- data_clustering_adj[data_clustering_adj$Cluster == k, , drop = FALSE]
      d <- sqrt(rowSums((cd[, num_cols, drop = FALSE] - km2$medoids[k, ])^2))
      thr <- stats::quantile(d, alpha, na.rm = TRUE)
      selected_list[[k]] <- cd[d <= thr, , drop = FALSE]
    }

    optimal_k <- k2
    optimal_states <- optimal_states[0, ]
    for (k in seq_len(optimal_k)) {
      repeat {
        any_update <- FALSE
        sel <- selected_list[[k]]
        if (!nrow(sel)) { selected_list[[k]] <- sel; break }
        center <- colMeans(sel[, num_cols, drop = FALSE], na.rm = TRUE)
        if (any(!is.finite(center))) { selected_list[[k]] <- sel; break }
        for (i in seq_len(nrow(sel))) {
          pid <- as.character(sel$ID[i])
          all_states_pid <- all_coef[all_coef$ID == pid, , drop = FALSE]
          pick <- pick_best_state(center, all_states_pid)
          if (sel$State[i] != pick["state"]) { sel$State[i] <- pick["state"]; any_update <- TRUE }
          optimal_states <- rbind(optimal_states,
                                  data.frame(ID = pid, Cluster = k,
                                             OptimalState = as.character(pick["state"]),
                                             MinDistance = as.numeric(pick["mind"])))
        }
        selected_list[[k]] <- sel
        if (!is.TRUE(any_update)) break
      }
    }
    optimal_states <- dplyr::distinct(optimal_states, ID, .keep_all = TRUE)
    optimal_states$ID <- as.character(optimal_states$ID)

    merged <- dplyr::inner_join(all_coef, optimal_states[, c("ID","OptimalState")], by = "ID")
    keep <- merged[merged$State == merged$OptimalState, , drop = FALSE]
    row_pos <- match(keep$ID, data_clustering_adj$ID)
    data_clustering_adj[row_pos, c(num_cols, "State")] <- keep[, c(num_cols, "State")]
  }

  # centers from representatives (fallback to all points if reps empty)
  if (nrow(optimal_states) > 0) {
    centers_mean <- data_clustering_adj %>%
      dplyr::filter(ID %in% optimal_states$ID) %>%
      dplyr::group_by(Cluster) %>%
      dplyr::summarise(dplyr::across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)),
                       .groups = "drop")
  } else {
    centers_mean <- data_clustering_adj %>%
      dplyr::group_by(Cluster) %>%
      dplyr::summarise(dplyr::across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)),
                       .groups = "drop")
  }

  if (verbose) {
    message("[SubtypeAware_Registration] final k = ", optimal_k)
    message("[SubtypeAware_Registration] Centers (mean):")
    print(centers_mean)
  }

  # classify remaining IDs vs centers (robust to NA/Inf)
  remain_ids <- setdiff(data_clustering_adj$ID, optimal_states$ID)
  if (length(remain_ids)) {
    robust_dist2 <- function(state_row, center_vec) {
      a <- as.numeric(state_row); b <- as.numeric(center_vec)
      mask <- is.finite(a) & is.finite(b)
      if (!any(mask)) return(NA_real_)
      sum((a[mask] - b[mask])^2)
    }
    for (pid in remain_ids) {
      all_states_pid <- all_coef[all_coef$ID == pid, c("State", num_cols), drop = FALSE]
      lv <- c("Original", paste("Shift", timepos_option))
      all_states_pid$State <- factor(all_states_pid$State, levels = lv, ordered = TRUE)
      all_states_pid <- all_states_pid[order(all_states_pid$State), , drop = FALSE]

      best <- NULL
      for (cidx in seq_len(nrow(centers_mean))) {
        center_vec <- as.numeric(centers_mean[cidx, num_cols, drop = TRUE])
        d2 <- apply(all_states_pid[, num_cols, drop = FALSE], 1, robust_dist2, center_vec = center_vec)
        if (all(!is.finite(d2))) next
        j <- which.min(replace(d2, !is.finite(d2), Inf))
        rec <- data.frame(
          ID = pid,
          Cluster = centers_mean$Cluster[cidx],
          OptimalState = as.character(all_states_pid$State[j]),
          MinDistance = as.numeric(d2[j])
        )
        if (is.null(best)) {
          if (is.finite(rec$MinDistance)) best <- rec
        } else {
          if (is.finite(rec$MinDistance) &&
              (!is.finite(best$MinDistance) || rec$MinDistance < best$MinDistance)) {
            best <- rec
          }
        }
      }
      if (is.null(best)) {
        best <- data.frame(ID = pid, Cluster = NA_integer_,
                           OptimalState = as.character(all_states_pid$State[1]),
                           MinDistance = NA_real_)
      }
      optimal_states <- rbind(optimal_states, best)
    }
  }

  # apply shifts to original data using character ID comparison
  parse_shift <- function(state) {
    if (identical(state, "Original")) return(0)
    as.numeric(trimws(sub("^Shift\\s*", "", state)))
  }
  shift_map <- stats::setNames(vapply(optimal_states$OptimalState, parse_shift, numeric(1)),
                               optimal_states$ID)

  out <- data
  id_chr <- as.character(out$ID)
  for (id in names(shift_map)) {
    sel <- id_chr == id
    out$Time[sel] <- out$Time[sel] + shift_map[[id]]
  }
  out <- out[out$Time >= tmin & out$Time <= tmax, , drop = FALSE]
  rownames(out) <- NULL
  out
}











# Mesh-to-polygon projection for cortical atlas creation ----
#
# Projects 3D inflated mesh triangles directly to 2D sf polygons,
# replacing the screenshot-based pipeline (steps 2-5).

#' @noRd
cross3 <- function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

#' Compute orthonormal view basis from camera position
#'
#' Replicates the Three.js OrthographicCamera look-at with up = (0,0,1).
#' @param camera_pos Numeric length-3 camera position vector.
#' @return Named list with `right` and `up` unit vectors (each length 3).
#' @noRd
compute_view_basis <- function(camera_pos) {
  forward <- -camera_pos / sqrt(sum(camera_pos^2))
  up <- c(0, 0, 1)

  right <- cross3(forward, up)
  right_len <- sqrt(sum(right^2))

  if (right_len < 1e-10) {
    up <- c(0, 1, 0)
    right <- cross3(forward, up)
    right_len <- sqrt(sum(right^2))
  }
  right <- right / right_len

  up_corrected <- cross3(right, forward)
  list(right = right, up = up_corrected)
}

#' Project 3D vertices to 2D using orthographic projection
#'
#' @param verts_3d N x 3 matrix of vertex positions.
#' @param basis View basis from `compute_view_basis()`.
#' @return N x 2 matrix of (x, y) screen coordinates.
#' @noRd
project_vertices_2d <- function(verts_3d, basis) {
  cbind(
    verts_3d %*% basis$right,
    verts_3d %*% basis$up
  )
}

#' Determine which mesh faces are front-facing for a given camera
#'
#' @param verts_3d N x 3 matrix of vertex positions.
#' @param faces F x 3 matrix of 1-indexed face vertex indices.
#' @param camera_pos Length-3 camera position vector.
#' @return Logical vector of length F.
#' @noRd
cull_backfaces <- function(verts_3d, faces, camera_pos) {
  v0 <- verts_3d[faces[, 1], , drop = FALSE]
  v1 <- verts_3d[faces[, 2], , drop = FALSE]
  v2 <- verts_3d[faces[, 3], , drop = FALSE]

  e1 <- v1 - v0
  e2 <- v2 - v0
  fn <- cbind(
    e1[, 2] * e2[, 3] - e1[, 3] * e2[, 2],
    e1[, 3] * e2[, 1] - e1[, 1] * e2[, 3],
    e1[, 1] * e2[, 2] - e1[, 2] * e2[, 1]
  )

  centers <- (v0 + v1 + v2) / 3
  to_cam <- matrix(camera_pos, nrow = nrow(fn), ncol = 3, byrow = TRUE) -
    centers

  rowSums(fn * to_cam) > 0
}

#' Build per-vertex label vector from vertices_df
#'
#' @param vertices_df Data frame with `label` and `vertices` (list-column of
#'   0-indexed integer vectors).
#' @param n_vertices Total number of vertices in the mesh.
#' @param hemi_short "lh" or "rh" — only labels starting with this prefix are
#'   included.
#' @return Character vector of length `n_vertices` (NA for unlabelled).
#' @noRd
build_vertex_label_vector <- function(vertices_df, n_vertices, hemi_short) {
  vertex_labels <- rep(NA_character_, n_vertices)
  prefix <- paste0(hemi_short, "_")

  for (i in seq_len(nrow(vertices_df))) {
    lbl <- vertices_df$label[i]
    if (!startsWith(lbl, prefix)) next
    idx <- vertices_df$vertices[[i]] + 1L
    idx <- idx[idx >= 1L & idx <= n_vertices]
    vertex_labels[idx] <- lbl
  }
  vertex_labels
}

#' Assign a face label based on its three vertex labels
#'
#' Uniform faces (all 3 vertices same label) get that label directly.
#' Boundary faces are assigned to the **smallest** neighboring region
#' so that small regions gain area and become more visible in plots.
#'
#' @param vertex_labels Character vector from `build_vertex_label_vector()`.
#' @param faces F x 3 matrix of 1-indexed vertex indices.
#' @return Character vector of length F (NA only if all vertices are NA).
#' @noRd
assign_face_labels <- function(vertex_labels, faces) {
  l1 <- vertex_labels[faces[, 1]]
  l2 <- vertex_labels[faces[, 2]]
  l3 <- vertex_labels[faces[, 3]]

  result <- rep(NA_character_, nrow(faces))

  all_same <- !is.na(l1) & !is.na(l2) & !is.na(l3) & l1 == l2 & l2 == l3
  result[all_same] <- l1[all_same]

  region_sizes <- table(vertex_labels[!is.na(vertex_labels)])

  boundary <- !all_same
  idx <- which(boundary)
  for (i in idx) {
    labs <- c(l1[i], l2[i], l3[i])
    labs <- unique(labs[!is.na(labs)])
    if (length(labs) == 0) next
    if (length(labs) == 1) {
      result[i] <- labs
      next
    }
    sizes <- region_sizes[labs]
    result[i] <- names(which.min(sizes))
  }

  result
}

#' Build sf polygon for one region from projected triangles
#'
#' @param verts_2d N x 2 projected vertex coordinates.
#' @param faces F x 3 matrix of 1-indexed vertex indices.
#' @param face_mask Logical vector selecting which faces to include.
#' @return An `sfc_POLYGON`/`sfc_MULTIPOLYGON`, or NULL if no faces.
#' @noRd
#' @importFrom sf st_polygon st_sfc st_union st_make_valid
build_region_polygon <- function(verts_2d, faces, face_mask) {
  sel <- faces[face_mask, , drop = FALSE]
  if (nrow(sel) == 0) return(NULL)

  triangles <- vector("list", nrow(sel))
  for (i in seq_len(nrow(sel))) {
    coords <- verts_2d[sel[i, ], , drop = FALSE]
    coords <- rbind(coords, coords[1, , drop = FALSE])
    triangles[[i]] <- sf::st_polygon(list(coords))
  }

  geom <- sf::st_make_valid(sf::st_union(sf::st_sfc(triangles)))
  geom
}

#' Project one hemisphere + one view to sf polygons
#'
#' @param mesh Brain mesh list with `vertices` and `faces` data frames.
#' @param vertex_labels Character vector from `build_vertex_label_vector()`.
#' @param camera_pos Length-3 camera position.
#' @param hemi_short "lh" or "rh".
#' @param view View name ("lateral", "medial", etc.).
#' @return sf data.frame with columns: filenm, hemi_short, hemi, view,
#'   label, geometry.
#' @noRd
#' @importFrom sf st_sf st_sfc
project_mesh_view <- function(mesh, vertex_labels, camera_pos,
                              hemi_short, view) {
  verts_3d <- as.matrix(mesh$vertices)
  faces_1idx <- as.matrix(mesh$faces) + 1L

  basis <- compute_view_basis(camera_pos)
  verts_2d <- project_vertices_2d(verts_3d, basis)

  visible <- cull_backfaces(verts_3d, faces_1idx, camera_pos)
  face_labels <- assign_face_labels(vertex_labels, faces_1idx)

  region_labels <- unique(face_labels[!is.na(face_labels) & visible])
  hemi_long <- if (hemi_short == "lh") "left" else "right"

  results <- vector("list", length(region_labels))
  for (j in seq_along(region_labels)) {
    lbl <- region_labels[j]
    mask <- visible & !is.na(face_labels) & face_labels == lbl
    poly <- build_region_polygon(verts_2d, faces_1idx, mask)
    if (!is.null(poly)) {
      results[[j]] <- data.frame(
        filenm = paste0(hemi_short, "_", view, "_", lbl),
        hemi_short = hemi_short,
        hemi = hemi_long,
        view = view,
        label = lbl,
        stringsAsFactors = FALSE
      )
      results[[j]]$geometry <- poly
    }
  }

  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0) return(NULL)

  combined <- do.call(rbind, results)
  sf::st_as_sf(combined)
}


#' Project mesh to 2D polygons for all hemispheres and views
#'
#' Replaces the screenshot pipeline (steps 2-5) with direct geometric
#' projection of mesh triangles.
#'
#' @param components Atlas components list with `vertices_df`.
#' @param hemisphere Character vector of hemisphere codes ("lh", "rh").
#' @param views Character vector of view names.
#' @param smooth_refinements Number of Chaikin corner-cutting refinements
#'   to apply before simplification. 0 = no smoothing.
#' @param verbose Logical.
#' @return sf data.frame with columns: filenm, hemi_short, hemi, view,
#'   label, geometry.
#' @noRd
#' @importFrom sf st_simplify st_make_valid
project_mesh_to_polygons <- function(components, hemisphere, views,
                                     tolerance = 0,
                                     smooth_refinements = 2,
                                     verbose = FALSE) {
  all_results <- list()

  for (hemi in hemisphere) {
    mesh <- ggseg.formats::get_brain_mesh(hemi, "inflated")
    n_verts <- nrow(mesh$vertices)
    vertex_labels <- build_vertex_label_vector(
      components$vertices_df, n_verts, hemi
    )

    for (view in views) {
      key <- paste(hemi, view, sep = "_")
      cam <- camera_presets[[key]]
      if (is.null(cam)) next

      if (verbose) {
        cli::cli_alert_info("Projecting {.val {hemi}} {.val {view}}")
      }

      result <- project_mesh_view(mesh, vertex_labels, cam, hemi, view)
      if (!is.null(result)) {
        all_results <- c(all_results, list(result))
      }
    }
  }

  if (length(all_results) == 0) {
    cli::cli_abort("No polygons generated from mesh projection")
  }

  sf_data <- do.call(rbind, all_results)

  if (smooth_refinements > 0) {
    rlang::check_installed("smoothr", reason = "for polygon smoothing")
    sf_data <- smoothr::smooth(sf_data, method = "chaikin",
                               refinements = smooth_refinements)
    sf_data <- sf::st_make_valid(sf_data)
  }

  if (tolerance > 0) {
    sf_data <- sf::st_simplify(sf_data, preserveTopology = TRUE,
                               dTolerance = tolerance)
    sf_data <- sf::st_make_valid(sf_data)
  }

  sf_data
}


#' Build cortical sf data from mesh projection
#'
#' Build cortical sf data using mesh projection.
#'
#' @param components Atlas components.
#' @param hemisphere Hemisphere codes.
#' @param views View names.
#' @param tolerance Simplification tolerance.
#' @param smooth_refinements Chaikin corner-cutting refinements.
#' @param verbose Logical.
#' @return sf data.frame with label, view, geometry columns.
#' @noRd
#' @importFrom dplyr group_by mutate ungroup select
#' @importFrom sf st_combine st_as_sf
cortical_build_sf_projected <- function(components, hemisphere, views,
                                        tolerance = 0,
                                        smooth_refinements = 2,
                                        verbose = FALSE) {
  projected <- project_mesh_to_polygons(
    components, hemisphere, views,
    tolerance = tolerance,
    smooth_refinements = smooth_refinements,
    verbose = verbose
  )

  projected |>
    layout_cortical_views() |>
    dplyr::group_by(view, label) |>
    dplyr::mutate(geometry = sf::st_combine(geometry)) |>
    dplyr::ungroup() |>
    dplyr::select(label, view, geometry) |>
    sf::st_as_sf()
}

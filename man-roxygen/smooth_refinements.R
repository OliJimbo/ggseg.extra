#' @param smooth_refinements Number of Chaikin corner-cutting refinements
#'   to apply to 2D polygons. Higher values produce smoother region
#'   boundaries (typical range: 0--3). 0 disables smoothing.
#'   If not specified, uses `options("ggseg.extra.smooth_refinements")` or the
#'   `GGSEG_EXTRA_SMOOTH_REFINEMENTS` environment variable. Default is 2.

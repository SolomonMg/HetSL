#' @title Plotting Methods for HetEffects
#' @description Visualization helpers for ensemble weights and treatment effects.
#' @name plotting
#' @keywords internal
NULL

#' Plot het_ensemble
#'
#' @param x A fitted \code{het_ensemble} object.
#' @param type Type of plot: "weights", "cate", "mcate", or "diagnostics".
#' @param moderator Moderator variable for MCATE plots.
#' @param ... Additional arguments.
#'
#' @export
plot.het_ensemble <- function(x,
                              type = c("weights", "cate", "mcate", "diagnostics"),
                              moderator = NULL,
                              ...) {
  type <- match.arg(type)

  if (type == "weights") {
    .plot_weights(x, ...)
  } else if (type == "cate") {
    .plot_cate(x, ...)
  } else if (type == "mcate") {
    if (is.null(moderator)) {
      stop("'moderator' is required for MCATE plots.", call. = FALSE)
    }
    .plot_mcate(x, moderator, ...)
  } else {
    .plot_diagnostics(x, ...)
  }
}

#' Plot ensemble weights
#' @keywords internal
.plot_weights <- function(object, ...) {
  weights <- sort(object$weights, decreasing = TRUE)
  df <- data.frame(learner = names(weights), weight = as.numeric(weights))

  if (.is_installed("ggplot2")) {
    ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = stats::reorder(rlang::.data$learner, rlang::.data$weight),
        y = rlang::.data$weight
      )
    ) +
      ggplot2::geom_col(fill = "#2a6f97") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Weight", title = "Ensemble Weights") +
      ggplot2::theme_minimal()
  } else {
    graphics::barplot(
      height = df$weight,
      names.arg = df$learner,
      las = 2,
      horiz = TRUE,
      col = "steelblue",
      xlab = "Weight",
      main = "Ensemble Weights"
    )
  }
}

#' Plot CATE distribution
#' @keywords internal
.plot_cate <- function(object, ...) {
  cate_values <- cate(object)

  if (.is_installed("ggplot2")) {
    df <- data.frame(cate = cate_values)
    ggplot2::ggplot(df, ggplot2::aes(x = rlang::.data$cate)) +
      ggplot2::geom_histogram(bins = 30, fill = "#2a6f97", color = "white") +
      ggplot2::labs(x = "CATE", y = "Count", title = "CATE Distribution") +
      ggplot2::theme_minimal()
  } else {
    graphics::hist(
      cate_values,
      breaks = 30,
      col = "steelblue",
      main = "CATE Distribution",
      xlab = "CATE"
    )
  }
}

#' Plot MCATE by moderator
#' @keywords internal
.plot_mcate <- function(object, moderator, ...) {
  mc <- mcate(object, moderator = moderator, ...)

  if (.is_installed("ggplot2")) {
    ggplot2::ggplot(
      mc,
      ggplot2::aes(x = rlang::.data$moderator_value, y = rlang::.data$estimate)
    ) +
      ggplot2::geom_line(color = "#2a6f97") +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = rlang::.data$ci_lower,
          ymax = rlang::.data$ci_upper
        ),
        alpha = 0.2,
        fill = "#2a6f97"
      ) +
      ggplot2::labs(
        x = moderator,
        y = "MCATE",
        title = sprintf("MCATE by %s", moderator)
      ) +
      ggplot2::theme_minimal()
  } else {
    plot(
      mc$moderator_value, mc$estimate,
      type = "l",
      col = "steelblue",
      xlab = moderator,
      ylab = "MCATE",
      main = sprintf("MCATE by %s", moderator)
    )
    if (!all(is.na(mc$ci_lower))) {
      graphics::polygon(
        c(mc$moderator_value, rev(mc$moderator_value)),
        c(mc$ci_lower, rev(mc$ci_upper)),
        col = grDevices::rgb(0.27, 0.49, 0.65, 0.2),
        border = NA
      )
      graphics::lines(mc$moderator_value, mc$estimate, col = "steelblue")
    }
  }
}

#' Diagnostics plot
#' @keywords internal
.plot_diagnostics <- function(object, ...) {
  fitted <- fitted.het_ensemble(object, type = "ensemble")
  y <- object$data$y

  if (.is_installed("ggplot2")) {
    df <- data.frame(fitted = fitted, y = y)
    ggplot2::ggplot(df, ggplot2::aes(x = fitted, y = y)) +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#2a6f97") +
      ggplot2::labs(x = "Fitted", y = "Observed", title = "Observed vs Fitted") +
      ggplot2::theme_minimal()
  } else {
    plot(
      fitted, y,
      xlab = "Fitted",
      ylab = "Observed",
      main = "Observed vs Fitted",
      pch = 19,
      col = grDevices::rgb(0.2, 0.4, 0.6, 0.4)
    )
    graphics::abline(0, 1, col = "gray50", lty = 2)
  }
}

#' Launch Shiny viewer for effects
#'
#' @param object A fitted \code{het_ensemble} object.
#' @param ... Additional arguments passed to \code{shiny::runApp}.
#'
#' @export
view_effects <- function(object, ...) {
  if (!inherits(object, "het_ensemble")) {
    stop("'object' must be a het_ensemble object.", call. = FALSE)
  }
  if (!.is_installed("shiny")) {
    stop("Package 'shiny' is required for view_effects(). Please install it.",
         call. = FALSE)
  }

  app_dir <- system.file("shiny", package = "HetEffects")
  if (app_dir == "") {
    stop("Shiny app not found in package installation.", call. = FALSE)
  }

  # Provide object to the app via global environment for simplicity.
  assign("heteffects_object", object, envir = .GlobalEnv)
  shiny::runApp(app_dir, ...)
}

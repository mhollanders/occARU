#' Set up CmdStan for use with occARU
#'
#' Installs CmdStan and checks the C++ toolchain required to compile Stan
#' models. This only needs to be run once. If CmdStan is already installed
#' and working, this function does nothing.
#'
#' @param ... Additional arguments passed to [cmdstanr::install_cmdstan()].
#' @export
setup_occARU <- function(...) {
  # pass immediately if a sufficient version is already installed
  current <- tryCatch(cmdstanr::cmdstan_version(), error = \(e) NULL)
  if (
    !is.null(current) && package_version(current) >= package_version("2.36.0")
  ) {
    cli::cli_inform(c(
      "v" = "CmdStan {current} already installed. Nothing to do."
    ))
  }

  # check toolchain before attempting install
  cli::cli_inform(c("i" = "Checking C++ toolchain..."))
  check <- tryCatch(
    cmdstanr::check_cmdstan_toolchain(fix = TRUE, quiet = TRUE),
    error = \(e) e
  )
  if (inherits(check, "error")) {
    cli::cli_abort(c(
      "C++ toolchain check failed.",
      "i" = "See {.url https://mc-stan.org/cmdstanr/articles/cmdstanr.html}
            for platform-specific setup instructions.",
      "x" = "{check$message}"
    ))
  }
  cli::cli_inform(c("v" = "Toolchain OK."))

  if (is.null(current)) {
    cli::cli_inform(c("i" = "CmdStan not found. Installing..."))
  } else {
    cli::cli_inform(c(
      "i" = "CmdStan {current} is too old (>= 2.36.0 required).
            Installing latest..."
    ))
  }
  cmdstanr::install_cmdstan(...)
  cli::cli_inform(c(
    "v" = "CmdStan {cmdstanr::cmdstan_version()} installed successfully."
  ))
}

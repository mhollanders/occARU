# Set up CmdStan for use with occARU

Installs CmdStan and checks the C++ toolchain required to compile Stan
models. This only needs to be run once. If CmdStan is already installed
and working, this function does nothing.

## Usage

``` r
setup_occARU(...)
```

## Arguments

- ...:

  Additional arguments passed to
  [`cmdstanr::install_cmdstan()`](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html).

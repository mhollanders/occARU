# The occARU model

occARU implements a Bayesian multispecies occupancy model with count
observation model for autonomous recording unit (ARU) data such as
camera traps and passive acoustic monitoring (PAM). Sites are defined as
individual ARUs, each producing time-stamped images or recordings with
species labels and individual counts. These counts are then binned into
surveys spanning an arbitrary number of days, with detections thinned to
appropriate intervals. The model is implemented in Stan ([Carpenter et
al. 2017](#ref-carpenter2017)) via
[cmdstanr](https://mc-stan.org/cmdstanr/) ([Gabry et al.
2025](#ref-gabry2025)), using Hamiltonian Monte Carlo (HMC) to sample
from the posterior. When only one species is provided, a single-species
veresion of the model is used automatically. The following sections
describe each model component in turn; Stan programs and associated
functions are available on
[GitHub](https://github.com/mhollanders/occARU/blob/main/inst/stan/occARU.stan).

## Occupancy

To accommodate structural zeros arising from species absences at the
site level, we model the occupancy of species s \in 1, \dots, S at site
i \in 1, \dots, I as a latent Bernoulli variable, where z\_{is} = 1
indicates occupancy:

\begin{aligned} z\_{is} &\sim \mathrm{Bernoulli} \left( \psi\_{is}
\right) \\ \psi\_{is} &= \operatorname{logit}^{-1} \left(
{\alpha\_\psi}\_{\[s\]} + {\boldsymbol{\beta\_\psi}}\_{\[s\]}
{\boldsymbol{X_1}}\_{\[i\]} \right), \end{aligned} \tag{1}

where \boldsymbol{\alpha\_\psi} are species-level intercepts,
\boldsymbol{\beta\_\psi} are species-level coefficients and
\boldsymbol{X_1} represents a design matrix. occARU accommodates three
types of predictors: continuous, categorical, and ordinal (ordered
categorical). Categorical coefficients are modeled as zero-sum Gaussian
variables, and ordinal coefficients are parameterised via a simplex
decomposition of the category differences ([Bürkner and Charpentier
2020](#ref-burkner2020)), where the coefficient corresponds to the
effect of the highest category.

Species-specific coefficients are modeled as multivariate normal,
facilitating inference on correlations in predictor responses across
species. Spatial effects on occupancy are not included by design because
they are weakly identified in the absence of repeated measures of
occupancy at the site level ([Doser et al. 2022](#ref-doser2022)).
Future versions of occARU will include multi-season occupancy models,
where spatial random effects are more appropriate.

Stan does not accommodate discrete parameters (z\_{is}) and they are
marginalised out of the likelihood; however, posterior inference on
underlying occupancy states is performed through the forward-backward
sampling algorithm by setting the argument `latent = TRUE` in
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

## Detection

The species counts during surveys j \in 1, \dots, J produce a
multi-species time series \boldsymbol{y}\_{ij} = (y\_{ij1}, \dots,
y\_{ijS}). Rather than reducing counts to binary detection/non-detection
as in traditional occupancy models ([MacKenzie et al.
2002](#ref-mackenzie2002)), occARU models the counts directly,
capitalising on the magnitude of detections. Conditional on z\_{is} = 1,
detections are modeled as Poisson by default:

\begin{aligned} y\_{ijs} \mid z\_{is} = 1 &\sim \mathrm{Poisson} \left(
\mu\_{ijs} \right) \\ \mu\_{ijs} &= \exp \Bigl( {\alpha\_\mu}\_{\[s\]} +
{\boldsymbol{\beta\_\mu}}\_{\[s\]} {\boldsymbol{X_2}}\_{\[i\]} +
\boldsymbol{\gamma}\_{\[s\]} {\boldsymbol{X_3}}\_{\[ij\]} \\ & \quad
\quad \quad \\ + \iota\_{is}\left(\boldsymbol{I} -
{\boldsymbol{P\_{X_2}}}\_{\[i\]} \right) +
\kappa\_{js}\left(\boldsymbol{I} - {\boldsymbol{P\_{X_3}}}\_{\[j\]}
\right) \\ & \quad \quad \quad \\ + \log \Delta\_{ij}\Bigr).
\end{aligned} \tag{2}

Here, \Delta\_{ij} are the proportion of each survey that the site was
actively sampling, taking on values between 0 and 1, used as offsets for
the detection rates \boldsymbol{\mu}\_{ij}. The species-specific
intercepts \boldsymbol{\alpha\_\mu} are modeled jointly with occupancy
intercepts as bivariate normal to facilitate estimating the
occupancy/detection correlation across species:

\begin{pmatrix} {\alpha\_\psi}\_{\[s\]} \\ {\alpha\_\mu}\_{\[s\]}
\end{pmatrix} \sim \mathrm{MVNormal} \left( \begin{pmatrix}
\operatorname{logit} \bar{\psi} \\ \log \bar{\mu} \end{pmatrix},
\begin{pmatrix} \tau\_\psi^2 & \rho \tau\_\psi \tau\_\mu \\ \rho
\tau\_\psi \tau\_\mu & \tau\_\mu^2 \end{pmatrix} \right). \tag{3}

As with occupancy ([Equation 1](#eq-z)), \boldsymbol{\beta\_\mu} are
species-level coefficients and \boldsymbol{X_2} represents a site-level
design matrix. Additionally, \boldsymbol{\gamma} are species-level
coefficients for the site-by-survey design matrix \boldsymbol{X_3}, used
for predictors that vary spatio-temporally. Again,
\boldsymbol{\beta\_\mu} and \boldsymbol{\gamma} are modeled as
multivariate normal to explore interspecific correlations in detection
rate responses to site-level and site-by-survey-level predictors. Both
\boldsymbol{X_2} and \boldsymbol{X_3} accommodate continuous,
categorical, and ordinal predictors.

### Site and survey random effects

To handle repeated measures at the level of sites and surveys, site
(\boldsymbol{\iota}\_i) and survey (\boldsymbol{\kappa}\_j) random
effects are included by default. These are modeled as matrix normal
distributions, where the row covariances are Gaussian process (GP)
kernels capturing spatial and temporal structure, and the column
covariances are interspecific covariance matrices capturing correlations
across species. This formulation accommodates effects that may differ in
magnitude across species but are nevertheless correlated:

\begin{aligned} \boldsymbol{\iota}\_i &\sim \mathrm{MatrixNormal} \left(
\bar{\iota}\_i, \boldsymbol{K\_{\iota}}, \Sigma\_{\iota} \right) \\
\boldsymbol{\kappa}\_j &\sim \mathrm{MatrixNormal} \left(
\bar{\kappa}\_j, \boldsymbol{K\_{\kappa}}, \Sigma\_{\kappa} \right).
\end{aligned}

Here, \bar{\boldsymbol{\iota}} and \bar{\boldsymbol{\kappa}} are mean
site and survey effects, themselves modeled as GPs;
\boldsymbol{K\_{\iota}} and \boldsymbol{K\_{\kappa}} are I \times I and
J \times J GP kernels; and \Sigma\_\iota and \Sigma\_\kappa are S \times
S interspecific covariance matrices, respectively. occARU uses
exponentiated quadratic kernels:

K\_{n, n'} = \tau^2 \exp \left(- \frac{\lVert \boldsymbol{x}\_n -
\boldsymbol{x}\_{n'} \rVert^2}{2\ell^2} \right), \tag{4}

where \lVert \boldsymbol{x}\_n - \boldsymbol{x}\_{n'} \rVert is the
Euclidean distance between inputs (site UTMs in kms or survey index),
\tau is the marginal scale of the GP, and \ell is the length scale
controlling the rate of decay. occARU facilitates fitting
species-specific GP kernels; however, this requires S additional
Cholesky decompositions for both site and survey effects. By default,
the same kernel is used for both the mean and species-specific GP
realisations, with realisations differing only in magnitude through
species-specific scales, requiring only one Cholesky decomposition per
GP.

occARU has the option to include a periodic kernel in the temporal GP
kernel \boldsymbol{K\_\kappa} (by setting `periodic_gp = TRUE` in
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md)):

K\_{n, n'} = \tau^2 \exp \left(- \frac{2 \sin^2 \left( \pi \lvert x_n -
x\_{n'} \rvert / p \right)}{\ell^2} \right), \tag{5}

where p is the period in survey units (e.g., 52 for an annual cycle with
weekly surveys) and \tau and \ell have the same interpretation as above.
Because GP kernels are closed under addition, a single Cholesky
decomposition is required. The variance partition of the two GP kernels
is handled via a simplex \boldsymbol{K\_\phi} =
\left({K\_\phi}\_{\[1\]}, {K\_\phi}\_{\[2\]}\right), \sum\_{g=1}^2
{K\_\phi}\_{\[g\]} = 1.

Site and survey random effects can be fully disabled or fitted without
GPs, collapsing to multivariate normals, with the arguments
`spatial = c("gp", "mvn", "none")` and
`temporal = c("gp", "mvn", "none")` in
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).
Identifiability of the matrix normal scales is ensured by fixing the
marginal scale of the GP kernels to 1.

### Orthogonal projection of random effects

When predictors are included alongside spatial or temporal random
effects, the fixed and random effects can be confounded because the
random effects absorb part of the variance that the fixed effects are
trying to explain ([Hodges and Reich 2010](#ref-hodges2010)). One remedy
is to construct the random effects to be orthogonal to the column space
of the design matrix. Taking the spatial random effects as an example,
instead of sampling \boldsymbol{\iota}, occARU samples
\boldsymbol{\iota} \left(\boldsymbol{I} - \boldsymbol{P\_{X_2}} \right),
where \boldsymbol{P\_{X_2}} := \boldsymbol{X_2^\prime}
\left(\boldsymbol{X_2} \boldsymbol{X_2^\prime} \right)^{-1}
\boldsymbol{X_2} and \boldsymbol{I} - \boldsymbol{P\_{X_2}} is the
orthogonal complement of the column space of the site-level design
matrix \boldsymbol{X_2}. Using so-called restricted random effects
essentially constrain random effects to the residual of the design
matrix and thus prioritise inference on the fixed effects.

occARU uses orthogonal projection for site effects and optionally survey
effects, by setting the argument `project_kappa = TRUE` (default) in
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md),
which is also more efficient for HMC. Because \boldsymbol{X_3} varies by
both site and survey, the site-averaged matrix is used to produce the
orthogonal projection. This assumes that survey predictors will be
correlated across sites; if these predictors vary considerably between
sites, the confounding will be limited and orthogonal projection is
unnecessary.

“Unconditional” fixed effect coefficients can still be recovered as,
e.g., \boldsymbol{\beta\_{\mu}} + \boldsymbol{\iota}
\left(\boldsymbol{X_2} \boldsymbol{X_2^\prime} \right)^{-1}
\boldsymbol{X_2} ([Hanks et al. 2015](#ref-hanks2015)). This facilitates
comparison of coefficients with and without the orthogonal projection of
the random effects, without refitting the model. occARU fit objects
contain both the restricted (e.g., `mu_beta`) and unrestricted (e.g.,
`mu_beta2`) coefficients, and both can be plotted with
[`plot_coefficients()`](https://mhollanders.github.io/occARU/reference/plot_coefficients.md)
using the `restricted` argument. Currently, “unrestricted” coefficients
are only available for continuous and ordinal predictors, but
categorical predictors are included in the projection.

### Overdispersion and residual species interactions

Two extensions beyond the basic Poisson observation model are available
via the `overdispersion = c("none", "nb", "olre")` argument in
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md):

1.  Negative binomial (`overdispersion = "nb"`): a species-specific
    overdispersion parameter \phi_s is estimated, so that y\_{ijs} \mid
    z\_{is} = 1 \sim \mathrm{NegBin}(\mu\_{ijs}, \phi_s).
2.  Observation-level random effects (OLREs)
    (`overdispersion = "olre"`): multivariate normal residuals
    \boldsymbol{\epsilon}\_{ij} = (\epsilon\_{ij1}, \dots,
    \epsilon\_{ijS}) are added to the linear predictor, yielding a
    Poisson-lognormal model. The OLREs are modeled hierarchically
    \boldsymbol{\epsilon}\_{ij} \sim \mathrm{MVN} \left(
    \bar\epsilon\_{ij}, \Sigma\_\epsilon \right), where
    \boldsymbol{\bar\epsilon} are mean residuals shared across species
    at each site-survey combination, capturing detection anomalies
    common to all species, and the covariance matrix \Sigma\_\epsilon
    captures residual interspecific correlations at the observation
    level—the site-by-survey-specific detection counts—which may reflect
    direct species interactions such as predator–prey dynamics.

## Variance decomposition

The model uses global-local shrinkage priors for variance components of
both occupancy and detection submodels ([Aguilar and Bürkner
2023](#ref-aguilar2023)). This accommodates a wide range of model
complexity (i.e., various number of covariates with different random
effects) in a principled unified prior framework. occARU places
half-Student-*t* priors on the variance of the linear predictor of
occupancy log odds (W\_\psi) and log detection rates (W\_\mu). These
variances are then simplex-decomposed between the different model
components (e.g., species-level intercepts, hierarchical coefficients,
random site and survey effects, and potentially Poisson OLREs). The
occARU model automatically determines the number of variance partitions
for each submodel (V\_\psi and V\_\mu) and produces the prior scales
using the V-length simplexes \boldsymbol{\phi\_\psi} and
\boldsymbol{\phi\_\mu} as, e.g., \boldsymbol{\tau\_\mu} := \left(\tau_1,
\dots, \tau\_{V\_\mu} \right) = \sqrt{W\_\mu \cdot
\boldsymbol{\phi\_\mu}}.

Variance decomposition can be done with either Dirichlet or zero-sum
logistic-normal distributions for \boldsymbol{\phi}, set via the
`variance_decomposition = c("dirichlet", "logistic-normal")` argument in
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).
In either case, the hyperparameters \theta\_\psi and \theta\_\mu control
the concentration or sparsity of the simplex: in Dirichlet models,
\theta is the inverse concentration, and in logistic-normal models, it
is the scale. For both underlying distributions, larger values of \theta
correspond to a more sparse simplex (i.e., where a small number of
components account for most of the variance), with smaller values
corresponding to more concentrated simplexes (i.e., where different
model components more evenly account for the variance).

Veriance decomposition is also used for species-specific scales in site
and survey effects (and potentially OLREs). The scales vector
\boldsymbol{\phi\_\mu} contains one partition for the species-level
deviations from the mean function for each, which are then simplex
decomposed to produce the species-specific scales. Parameter named in
the fit objects are `iota_t`, `kappa_t`, and `epsilon_t`, respectively.

## Identifiability and sum-to-zero constraints

All random effects (species-level intercepts, fixed effect coefficients,
site effects, survey effects, and OLRE residuals) are fit with
sum-to-zero constraints for parameter identifiability ([Ogle and Barber
2020](#ref-ogle2020)). Without these constraints, the overall intercepts
(\bar{\psi} and \bar{\mu}), mean effects (i.e.,
\bar{\boldsymbol{\beta}}, \bar{\boldsymbol{\gamma}},
\boldsymbol{\bar{\iota}}, \boldsymbol{\bar{\kappa}}, and
\boldsymbol{\bar{\epsilon}}) and the species-level deviations are weakly
identified. The constraints are implemented directly in the Stan
parameterisation using `sum_to_zero_vector` and `sum_to_zero_matrix`
types, which induce the constraint with a favourable posterior geometry
leading to increased computational efficiency. Note that technically,
the resulting distributions are no longer (multivariate/matrix) normal,
but rather the zero-sum variants.

## Assessing model fit

occARU produces site-by-species level log likelihoods (`log_lik`) as
generated quantities; these values can be used with Pareto-smoothed
importance sampling leave-one-out cross-validation (PSIS-LOO-CV,
[Vehtari et al. 2017](#ref-vehtari2017)) as implemented in the
[loo](https://mc-stan.org/loo/) package ([Vehtari et al.
2025](#ref-vehtari2025)). PSIS-LOO-CV is known to be problematic when
there are observation-level parameters, i.e. \iota\_{is} and, when OLREs
are included, \epsilon\_{ijs}. As a remedy, Monte Carlo integration of
these random effects is used to produce a second log likelihood object
(`log_lik2`), where the user can set the number of Monte Carlo draws to
perform per HMC iteration (argument `loo_draws` in
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md)).

Posterior predictions can be produced for both y\_{ijs} and Q\_{is} :=
\sum\_{j=1}^J y\_{ijs}, by setting the argument
`ppc = c("Q", "y", "both", "none")` in
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).
This facilitates posterior predictive checking, i.e. via the
[bayesplot](https://mc-stan.org/bayesplot/) package ([Gabry et al.
2025](#ref-gabry2025)). Additionally, the log prior density of the prior
hyperparameters is stored (`lprior`), facilitating prior predictive
checking via the [priorsense](https://n-kall.github.io/priorsense/)
package ([Kallioinen et al. 2023](#ref-kallioinen2023)).

## Priors

Default priors for hyperparameters are chosen to be weakly informative
on the natural scale of the outcome (species counts). Prior families of
parameters are fixed, but hyperparameters can be configured. They can be
inspected and modified with
[`set_priors()`](https://mhollanders.github.io/occARU/reference/set_priors.md),
which produces an `occARU_priors` object to be passed to
[`fit_model()`](https://mhollanders.github.io/occARU/reference/fit_model.md).

| Parameter                               | Interpretation                                        | Default Prior                             |
|-----------------------------------------|-------------------------------------------------------|-------------------------------------------|
| \bar{\psi}                              | Mean occupancy probability                            | \mathrm{Beta} \left( 1, 1 \right)         |
| \bar{\mu}                               | Mean detection rate                                   | \mathrm{Gamma} \left( 1, 1 \right)        |
| W\_{\psi}                               | Variance of occupancy log odds linear predictor       | \mathit{t}^+ \left( 3, 0, 1 \right)       |
| W\_{\mu}                                | Variance of log detection rate linear predictor       | \mathit{t}^+ \left( 3, 0, 2.5 \right)     |
| \theta\_{\psi}, \theta\_{\mu}           | Variance partition simplex sparsity                   | \mathrm{Gamma} \left( 1, 1 \right)        |
| \ell\_\iota, \boldsymbol{\ell\_\kappa}  | GP length scales                                      | \mathrm{InvGamma} \left( 1, 1 \right)     |
| \phi_s                                  | Negative binomial overdispersion                      | \mathrm{InvGamma} \left( 0.4, 0.3 \right) |
| Cholesky factor of correlation matrices | Intercepts and coefficient/random effect correlations | \mathrm{LKJ} \left( 1 \right)             |

## References

Aguilar, Javier Enrique, and Paul-Christian Bürkner. 2023. “Intuitive
Joint Priors for Bayesian Linear Multilevel Models: The R2D2M2 Prior.”
*Electronic Journal of Statistics* 17 (1): 1711–67.
<https://doi.org/10.1214/23-EJS2136>.

Bürkner, Paul-Christian, and Emmanuel Charpentier. 2020. “Modelling
Monotonic Effects of Ordinal Predictors in Bayesian Regression Models.”
*British Journal of Mathematical & Statistical Psychology* 73 (3):
420–51. <https://doi.org/10.1111/bmsp.12195>.

Carpenter, Bob, Andrew Gelman, Matthew D. Hoffman, et al. 2017. “Stan: A
Probabilistic Programming Language.” *Journal of Statistical Software*
76: 1. <https://doi.org/10.18637/jss.v076.i01>.

Doser, Jeffrey W., Andrew O. Finley, Marc Kéry, and Elise F. Zipkin.
2022. “spOccupancy: An R Package for Single-Species, Multi-Species, and
Integrated Spatial Occupancy Models.” *Methods in Ecology and Evolution*
13 (8): 1670–78. <https://doi.org/10.1111/2041-210X.13897>.

Gabry, Jonah, Rok Češnovar, Andrew Johnson, and Steve Bronder. 2025.
*cmdstanr: R Interface to ’CmdStan’*. Manual.

Hanks, Ephraim M., Erin M. Schliep, Mevin B. Hooten, and Jennifer A.
Hoeting. 2015. “Restricted Spatial Regression in Practice:
Geostatistical Models, Confounding, and Robustness Under Model
Misspecification.” *Environmetrics* 26 (4): 243–54.
<https://doi.org/10.1002/env.2331>.

Hodges, James S., and Brian J. Reich. 2010. “Adding Spatially-Correlated
Errors Can Mess up the Fixed Effect You Love.” *The American
Statistician* 64 (4): 325–34. <https://doi.org/10.1198/tast.2010.10052>.

Kallioinen, Noa, Topi Paananen, Paul-Christian Bürkner, and Aki Vehtari.
2023. “Detecting and Diagnosing Prior and Likelihood Sensitivity with
Power-Scaling.” *Statistics and Computing* 34 (57).
<https://doi.org/10.1007/s11222-023-10366-5>.

MacKenzie, Darryl I., James D. Nichols, Gideon B. Lachman, Sam Droege,
J. Andrew Royle, and Catherine A. Langtimm. 2002. “Estimating Site
Occupancy Rates When Detection Probabilities Are Less Than One.”
*Ecology* 83 (8): 2248–55.
<https://doi.org/10.1890/0012-9658(2002)083%5B2248:ESORWD%5D2.0.CO;2>.

Ogle, Kiona, and Jarrett J. Barber. 2020. “Ensuring Identifiability in
Hierarchical Mixed Effects Bayesian Models.” *Ecological Applications*
30 (7): e02159. <https://doi.org/10.1002/eap.2159>.

Vehtari, Aki, Jonah Gabry, Måns Magnusson, et al. 2025. *loo: Efficient
Leave-One-Out Cross-Validation and WAIC for Bayesian Models*.

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017. “Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and WAIC.”
*Statistics and Computing* 27 (5): 1413–32.
<https://doi.org/10.1007/s11222-016-9696-4>.

---
editor_options: 
  markdown: 
    wrap: sentence
---

# Relational Mortality Models: Rates and Counts

This repository provides tools for fitting relational mortality models using a flexible, two-step Generalized Additive Model (GAM) approach.

The primary utility of these scripts is to convert abridged, sparse, or noisy empirical mortality data into a smooth, continuous, single-year mortality schedule.
It achieves this by anchoring the empirical data to a known standard demographic schedule and using a penalized spline to smooth the localized deviations.

This repository offers two distinct implementations: a **Count-Based framework (Recommended)** and a **Rate-Based framework** (if you don't have access to death counts and exposures).

------------------------------------------------------------------------

### The Shared Framework: Solving the Identifiability Problem

Both models (for counts and rates) face a fundamental identifiability challenge: the standard mortality schedule and the smoothing spline are highly collinear over age $x$.

If we attempt to estimate the global linear trend against the standard schedule and the flexible spline simultaneously, the optimizer faces a nearly flat likelihood surface.
I most cases, the model will correctly use the standard curve to explain most of the variation and the spline-based terms to model the residual variation.
However, in some cases the spline will "steal" variance from the standard curve, effectively ignoring the structural demographic information.

Here I enforce a strict hierarchy using a **two-step sequential fitting strategy**: 1.
**Baseline Model:** We force the standard demographic curve to explain as much variance as possible.
2.
**Flexible Residuals:** We fit the GAM spline strictly to the residuals (or offsets) of the baseline model, smoothing out only the localized bumps and dips.

This acts as a rigid, practical constraint: the model *must* use the standard demographic schedule first, and only deploys the spline to mop up the remaining localized deviations.

------------------------------------------------------------------------

### Approach 1: The Count-Based Model (Negative Binomial)

**File:** `mortality_models_counts.R`

This is the preferred approach because it models the actual data-generating process: discrete death counts ($D_x$) occurring over a continuous volume of exposure ($E_x$), expressed in person-years.

Using summary rates ($m_x = D_x / E_x$) discards critical information about the variance of the estimate (because all the information on the population size is lost).
By modeling counts directly, we naturally handle age intervals with exactly zero deaths, and we correctly weight observations based on their population size.

$$D_x \sim \text{NB}(\mu_x, \theta)$$ $$\log(\mu_x) = \log(E_x) + \alpha + \beta \log(m_x^{\text{std}}) + f(x)$$

**How the Two-Step Model Works Here:** We fit the baseline model to the counts using $\log(E_x)$ as a fixed offset.
For the second step, we extract the baseline linear predictor and add it to the exposure, passing the sum `log(E_x) + baseline_rate` as the new offset to the GAM.
This guarantees the spline models pure structural deviations.

------------------------------------------------------------------------

### Approach 2: The Rate-Based Model (OLS)

**File:** `mortality_models_rates.R`

This approach follows the traditional demographic practice of treating the log-transformed mortality rate as the fundamental unit of observation.

$$\log(m_x) = \alpha + \beta \log(m_x^{\text{std}}) + f(x) + \epsilon$$

**How the Two-Step Model Works Here:** The function first fits an Ordinary Least Squares (OLS) model to extract the global trend.
The penalized spline is then fit directly to the continuous residuals of that baseline model.

**Limitations:**

-   **Non-Zero Rates:** Because it relies on $\log(m_x)$, the input data cannot contain rates of exactly zero.

-   **Ignores Population Size:** The model has no notion of the reliability of the rates and treats rates derived from a population of 100 and a population of 1,000,000 as having the same information content.
    It is best used for large populations without structural zeros.

------------------------------------------------------------------------

### Core Functions

Both `.R` files contain the same main functions:

-   `singularize()`: The primary user-facing function. Takes empirical data and standard schedules, handles missing bounds, fits the sequence, and outputs predicted continuous single-year rates.
-   `fit_two_step_model()`: The internal engine handling the sequential fits and offsets.
-   `predict_two_step_model()`: Handles safe prediction over granular age grids.
-   `fit_kannisto()`: Fits a logistic extrapolation to older ages (default 80-95) to ensure biologically plausible asymptotic closure at the oldest ages.

#### Handling Data Artifacts

Both models include a `skip_ages` argument within `singularize()`.
Users can pass a vector `c(min, max)` or a list of vectors `list(c(0, 5), c(20, 25))` to flag localized reporting artifacts or age-heaping.
The models treat these bands as missing data and interpolate over them using the standard curve and the spline.

------------------------------------------------------------------------

### Testing and Examples

This repository includes two testing scripts that reproduce standard life table outputs:

1.  **`test_count_models.R`**: Demonstrates the Negative Binomial model.
    -   *Test 1* replicates life tables from the Human Mortality Database, converting abridged inputs to single-year continuous schedules.
    -   *Test 2* stress-tests the model using problematic data. It explicitly introduces a $0$ death count at age 10 to demonstrate how the count-based likelihood processes structural zeros, and utilizes the `skip_ages` argument to bridge over selected age ranges.
2.  **`test_rate_models.R`**: Demonstrates the OLS rate-based model using the same datasets (excluding the zero-count test, which violates the log-transform assumption), providing a direct methodological comparison. Both scripts utilize the `tinyplot` package to visualize the resulting mortality curves and survival functions.

### Utility Functions: Classical Life Table Construction

In addition to the relational mortality models, this repository includes standard demographic utility functions to construct life tables from abridged or single-year data.

**File:** `lt_functions.R`

This file provides two functions for computing standard life table columns ($q_x, l_x, d_x, L_x, T_x, e_x$):

-   `lt_abridged(nMx)`: Constructs a life table for classical abridged age groups (0, 1-4, 5-9, 10-14, ..., 90+).

-   `lt_single(Mx)`: Constructs a life table for complete, single-year age groups (0, 1, 2, ..., 100+).

#### The Mathematical Translation and $a_x$ Assumptions

The fundamental step in building a life table is converting the observed mortality *rate* ($m_x$) into a probability of death ($q_x$).
This conversion relies on $a_x$, which represents the average fraction of the interval lived by those who die within that interval.

These functions use standard demographic approximations for $a_x$:

-   **Constant Mid-Point Assumption:** For most closed intervals, deaths are assumed to be distributed uniformly, meaning $a_x = 0.5$.

-   **Infant/Childhood Mortality:** Because infant deaths are heavily skewed toward the first few days of life, $a_x$ must be adjusted.

    -   In `lt_single()`, $a_0$ is hardcoded to $0.14$.

    -   In `lt_abridged()`, $a_0$ is estimated using a standard linear approximation: $0.07 + 1.7 \cdot m_0$, and $a_{1-4}$ is set to $0.4$.

-   **Open-Ended Interval:** For the final age group, the functions use the reciprocal of the mortality rate ($1/m_x$) to close the table, which mathematically assumes an exponentially distributed time of death for the oldest old.

*Note: These* $a_x$ approximations are sex-neutral and generalized.
If you are working with specific regional data where infant mortality distribution is highly distinct, you may need to manually adjust the $a_x$ vectors inside the functions.

#### Testing and Validation

**File:** `test_lt.R`

This script provides validation for the life table functions by comparing their output against official tables from the Human Mortality Database (HMD).

**Single-Year Test:** Computes $e_x$ using `lt_single()` and plots it against HMD single-year life expectancy.
Because the $a_x = 0.5$ assumption cleanly matches single-year data, the results are essentially identical.

**Abridged Test:** Computes $e_x$ using `lt_abridged()` and plots it against HMD abridged life expectancy.
You will notice slight, expected deviations here.
This is because the HMD uses highly specific, country- and sex-dependent formulas for $a_0$ and $a_{1-4}$, whereas this utility function relies on a generalized constant approximation.

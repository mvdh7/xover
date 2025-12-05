from collections import namedtuple

import jax
import numpy as onp
import statsmodels.api as sm
from jax import numpy as np
from jax.lax import scan
from scipy.stats import t as t_dist


jax.config.update("jax_enable_x64", True)

# namedtuples for function results
FurthestFirstSteps = namedtuple(
    "FurthestFirstSteps",
    ("index", "adjustment", "uncertainty"),
)
FurthestFirstResult = namedtuple(
    "FurthestFirstResult",
    (
        "adjustments",
        "uncertainties",
        "steps",
        "xovers_adjusted",
        "weights",
        "allowed",
        "dof",
        "t_crit",
        "niter",
    ),
)
TrendsWLS = namedtuple("TrendsWLS", ("intercept", "slope", "slope_se"))


def rms(values, **kwargs):
    """Calculate the root-mean-square of a set of values.

    Additional `kwargs` are passed on to `np.mean`.
    """
    return np.sqrt(np.mean(values**2, **kwargs))


def rms_weighted(values, weights, **kwargs):
    """Calculate the weighted root-mean-square of a set of values.

    Additional `kwargs` are passed on to `np.sum`.
    """
    return np.sqrt(
        np.sum(weights * values**2, **kwargs) / np.sum(weights, **kwargs)
    )


def adjust_xovers(xovers, adjustments):
    """Apply a set of adjustments to a crossover matrix.

    Parameters
    ----------
    xovers : np.array
        A square `(n, n)` matrix of the crossover differences, where `n` is the
        total number of cruises being compared.  `xovers[i, j]` contains the
        crossover difference for cruise `i` - cruise `j`, or zero if there is
        no crossover for that pair.
    adjustments : np.array
        A size `(n,)` array of the current set of adjustments to be applied.

    Returns
    -------
    np.array
        A square `(n, n)` matrix of the adjusted crossover differences.
    """
    return xovers + np.vstack(adjustments) - adjustments


def offsets_weighted(xovers, weights):
    """Calculate weighted mean offsets for a set of crossovers.

    Parameters
    ----------
    xovers : np.array
        A square `(n, n)` matrix of the crossover differences, where `n` is the
        total number of cruises being compared.  `xovers[i, j]` contains the
        crossover difference for cruise `i` - cruise `j`, or zero if there is
        no crossover for that pair.
    weights : np.array
        A square `(n, n)` matrix of the crossover weights.  `weights[i, j]`
        contains the weight for the corresponding crossover, and should be set
        to zero where there are no crossovers and along the main diagonal.

    Returns
    -------
    np.array
        A size `(n,)` array of the mean weighted offsets.
    """
    return np.nansum(xovers * weights, axis=1) / np.sum(weights, axis=1)


def dof_kish(weights):
    """Kish's effective degrees of freedom for weighted data.

    Parameters
    ----------
    weights : np.array
        A square `(n, n)` matrix of the crossover weights.  `weights[i, j]`
        contains the weight for the corresponding crossover, and should be set
        to zero where there are no crossovers and along the main diagonal.

    Returns
    -------
    np.array
        A size `(n,)` vector of the effective degrees of freedom for each row
        of `weights`.
    """
    return np.sum(weights, axis=1) ** 2 / np.sum(weights**2, axis=1)


def std_weighted(xovers, weights, offsets, dof):
    """Calculate weighted standard deviations for a set of crossover offsets,
    for use with `offsets_weighted`.

    Parameters
    ----------
    xovers : np.array
        A square `(n, n)` matrix of the crossover differences, where `n` is the
        total number of cruises being compared.  `xovers[i, j]` contains the
        crossover difference for cruise `i` - cruise `j`, or zero if there is
        no crossover for that pair.
    weights : np.array
        A square `(n, n)` matrix of the crossover weights.  `weights[i, j]`
        contains the weight for the corresponding crossover, and should be set
        to zero where there are no crossovers and along the main diagonal.
    offsets : np.array
        A size `(n,)` vector of crossover offsets, calculated from `xovers`
        and `weights` with `offsets_weighted`.
    dof : np.array
        A size `(n,)` vector of the degrees of freedom for each cruise, for
        example calculated with `dof_kish`.

    Returns
    -------
    np.array
        A size `(n,)` array of the weighted standard deviations of the offsets.
    """
    residuals = xovers - np.vstack(offsets)
    variance = (
        np.sum(weights * residuals**2, axis=1)
        / (np.sum(weights, axis=1))
        * (dof / (dof - 1))
    )
    return np.sqrt(variance)


def sem_weighted(xovers, weights, offsets, dof):
    """Calculate the standard errors in weighted means.

    Parameters
    ----------
    xovers : np.array
        A square `(n, n)` matrix of the crossover differences, where `n` is the
        total number of cruises being compared.  `xovers[i, j]` contains the
        crossover difference for cruise `i` - cruise `j`, or zero if there is
        no crossover for that pair.
    weights : np.array
        A square `(n, n)` matrix of the crossover weights.  `weights[i, j]`
        contains the weight for the corresponding crossover, and should be set
        to zero where there are no crossovers and along the main diagonal.
    offsets : np.array
        A size `(n,)` vector of crossover offsets, calculated from `xovers`
        and `weights` with `offsets_weighted`.
    dof : np.array
        A size `(n,)` vector of the degrees of freedom for each cruise, for
        example calculated with `dof_kish`.

    Returns
    -------
    np.array
        A size `(n,)` array of the standard errors in the weighted means.
    """
    return std_weighted(xovers, weights, offsets, dof) / np.sqrt(dof)


def t_critical(dof, window=0.68):
    """Critical t values corresponding to the provided degrees of freedom and
    confidence interval window.

    Parameters
    ----------
    dof : np.array
        A size `(n,)` vector of the degrees of freedom for each cruise, for
        example calculated with `dof_kish`.
    window : float, optional
        Confidence interval to define the critical t, by default 0.68.

    Returns
    -------
    np.array
        A size `(n,)` array of the critical t value for each `dof` and for the
        specified confidence `window`.
    """
    alpha = 1 - window  # window = 0.68 for 1-sigma uncertainty
    return t_dist.ppf(1 - alpha / 2, df=dof)


def sem_weighted_t(xovers, weights, offsets, dof, t_crit):
    """Standard error in the mean accounting for number of samples by using
    the t-distribution.

    Parameters
    ----------
    xovers : np.array
        A square `(n, n)` matrix of the crossover differences, where `n` is the
        total number of cruises being compared.  `xovers[i, j]` contains the
        crossover difference for cruise `i` - cruise `j`, or zero if there is
        no crossover for that pair.
    weights : np.array
        A square `(n, n)` matrix of the crossover weights.  `weights[i, j]`
        contains the weight for the corresponding crossover, and should be set
        to zero where there are no crossovers and along the main diagonal.
    offsets : np.array
        A size `(n,)` vector of crossover offsets, calculated from `xovers`
        and `weights` with `offsets_weighted`.
    dof : np.array
        A size `(n,)` vector of the degrees of freedom for each cruise, for
        example calculated with `dof_kish`.
    t_crit : np.array
        A size `(n,)` array of the critical t value for each `dof` and for the
        specified confidence `window`.

    Returns
    -------
    np.array
        A size `(n,)` array of the standard errors in the weighted means.
    """
    return t_crit * sem_weighted(xovers, weights, offsets, dof)


def offset_uncertainties(xovers, weights, offsets, dof, t_crit):
    """Final uncertainties in the offsets, including applying a default value
    for cruises with only one crossover.

    Parameters
    ----------
    xovers : np.array
        A square `(n, n)` matrix of the crossover differences, where `n` is the
        total number of cruises being compared.  `xovers[i, j]` contains the
        crossover difference for cruise `i` - cruise `j`, or zero if there is
        no crossover for that pair.
    weights : np.array
        A square `(n, n)` matrix of the crossover weights.  `weights[i, j]`
        contains the weight for the corresponding crossover, and should be set
        to zero where there are no crossovers and along the main diagonal.
    offsets : np.array
        A size `(n,)` vector of crossover offsets, calculated from `xovers`
        and `weights` with `offsets_weighted`.
    dof : np.array
        A size `(n,)` vector of the degrees of freedom for each cruise, for
        example calculated with `dof_kish`.
    t_crit : np.array
        A size `(n,)` array of the critical t value for each `dof` and for the
        specified confidence `window`.

    Returns
    -------
    np.array
        A size `(n,)` array of the final uncertainties in the offsets.
    """
    offsets_u = sem_weighted_t(xovers, weights, offsets, dof, t_crit)
    offsets_u = np.where(
        np.isnan(offsets_u) | np.isinf(offsets_u),
        2 * np.std(offsets, ddof=1),  # this is the single-crossover value
        offsets_u,
    )
    return offsets_u


def get_trends_wls(xovers, weights, dates, min_xovers=3):
    """Calculate weighted least squares trends in crossovers.

    A minimum number of crossovers must be available in order to calculate each
    trend, otherwise zero is returned.

    Parameters
    ----------
    xovers : np.array
        A square `(n, n)` matrix of the crossover differences, where `n` is the
        total number of cruises being compared.  `xovers[i, j]` contains the
        crossover difference for cruise `i` - cruise `j`, or zero if there is
        no crossover for that pair.
    weights : np.array
        A square `(n, n)` matrix of the crossover weights.  `weights[i, j]`
        contains the weight for the corresponding crossover, and should be set
        to zero where there are no crossovers and along the main diagonal.
    dates : np.array
        A square `(n, n)` matrix of the difference in time between crossover
        pairs, in decimal years.  `dates[i, j]` contains the date difference
        for cruise `i` - cruise `j`, or zero if there is no crossover for that
        pair.
    min_xovers : int, optional
        The minimum number of crossovers required to calculate each trend, by
        default 3.  The number of crossovers is determined from where the
        `weights` are non-zero.

    Returns
    -------
    TrendsWLS
        A namedtuple with the following fields:
            intercept : np.array
                A size `(n,)` array of the y-axis intercept value for each
                trend, in the unit of `xovers`.
            slope : np.array
                A size `(n,)` array of the trends, in `xovers` unit per year.
                Value is zero where no trend could be computed due to
                insufficient crossovers.
            slope_se : np.array
                A size `(n,)` array of the standard error in the trends, in
                `xovers` unit per year.  Value is `np.nan` where no trend could
                be computed due to insufficient crossovers.
    """
    offsets = offsets_weighted(xovers, weights)
    intercept = onp.zeros(offsets.shape)
    slope = onp.zeros(offsets.shape)
    slope_se = onp.zeros(offsets.shape)
    for i in range(offsets.size):
        L = weights[i] > 0
        if L.sum() >= min_xovers:
            wls = sm.WLS(
                onp.array(xovers[i][L]),
                sm.add_constant(dates[i][L]),
                weights=weights[i][L],
            ).fit()
            intercept[i], slope[i] = wls.params
            slope_se[i] = wls.bse[1]
        else:
            intercept[i] = offsets[i]
            slope[i] = 0
            slope_se[i] = np.nan
    return TrendsWLS(intercept, slope, slope_se)


def furthest_first(
    xovers,
    niter=None,
    allowed=None,
    dof=None,
    t_crit=None,
    weights=None,
):
    """Furthest-first iterative adjustment calculator.

    Parameters
    ----------
    xovers : np.array
        A square `(n, n)` matrix of the crossover differences, where `n` is the
        total number of cruises being compared.  `xovers[i, j]` contains the
        crossover difference for cruise `i` - cruise `j`, or zero if there is
        no crossover for that pair.
    niter : int, optional
        Number of iterations of the algorithm to run, by default `n`.
    allowed : np.array or None, optional
        A size `(n,)` Boolean array of which cruises are allowed to receive
        adjustments.  If not provided, all cruises are allowed.
    dof : np.array or None, optional
        (Effective) degrees of freedom for each row of `xovers`. If not
        provided, it is calculated from the `weights` and Kish's formula.
    t_crit : np.array or None, optional
        Critical value for the t-distribution. If not provided, it is
        calculated from `dof` with a 1-sigma uncertainty.
    weights : np.array, optional
        A square `(n, n)` matrix of the crossover weights.  `weights[i, j]`
        contains the weight for the corresponding crossover, and should be set
        to zero where there are no crossovers and along the main diagonal.
        If not provided, all crossovers are weighted equally, noting that every
        non-NaN element of the crossover matrix is considered to represent a
        crossover, including zero values.

    Returns
    -------
    FurthestFirstResult
        A namedtuple with the following fields:
            adjustments : np.array
                A `(n,)` array of the suggested adjustment for each cruise.
            uncertainties : np.array
                A `(n,)` array of the 1-sigma uncertainty in the suggested
                adjustment for each cruise.
            steps : FurthestFirstSteps
                A namedtuple with the following fields:
                    index : np.array
                        A `(niter,)` array of the index of the cruise adjusted
                        at each step of the iteration.
                    adjustment : np.array
                        A `(niter,)` array of the adjustment applied at each
                        step.
                    uncertainty : np.array
                        A `(niter,)` array of the uncertainty in the adjustment
                        at each step.
            xovers_adjusted : np.array
                A square `(n, n)` matrix of the crossover differences after
                the final set of suggested adjustments has been applied.
            weights : np.array
                The input or internally generated `weights`.
            allowed : np.array
                The input or internally generated `allowed`.
            dof : np.array
                The input or internally generated `dof`.
            t_crit : np.array
                The input or internally generated `t_crit`.
            niter : int
                The input or internally generated `niter`.
    """

    def scan_ff(adjustments, x):
        xovers_adj = adjust_xovers(xovers, adjustments)
        offsets = offsets_weighted(xovers_adj, weights)
        # Calculate offset uncertainties, assigning a value for cruises that
        # have only one offset
        offsets_u = offset_uncertainties(
            xovers_adj, weights, offsets, dof, t_crit
        )
        # Determine which cruise is to be adjusted
        offsets_norm = np.abs(offsets / offsets_u)
        a = np.argmax(allowed * offsets_norm)
        # Make the adjustment
        adjustments = adjustments.at[a].subtract(offsets[a])
        return adjustments, [a, -offsets[a], offsets_u[a]]

    # Assign defaults if not provided
    if niter is None:
        niter = xovers.shape[0]
    if allowed is None:
        allowed = np.full(xovers.shape[0], True)
    if weights is None:
        weights = ~np.isnan(xovers)
    if dof is None:
        dof = dof_kish(weights)
    if t_crit is None:
        t_crit = t_critical(dof)
    # Run the inversion
    adjustments = np.zeros(xovers.shape[0])
    adjustments, steps = scan(scan_ff, adjustments, length=niter)
    xovers_adjusted = adjust_xovers(xovers, adjustments)
    # Package up outputs
    steps = FurthestFirstSteps(
        onp.array(steps[0].astype(int)),
        onp.array(steps[1]),
        onp.array(steps[2]),
    )
    ffr = FurthestFirstResult(
        onp.array(adjustments),
        None,
        steps,
        onp.array(xovers_adjusted),
        weights,
        allowed,
        dof,
        t_crit,
        niter,
    )
    # Compute adjustment uncertainties
    uncertainties = uncertainties_at_step(ffr)
    return FurthestFirstResult(
        onp.array(adjustments),
        uncertainties,
        steps,
        onp.array(xovers_adjusted),
        weights,
        allowed,
        dof,
        t_crit,
        niter,
    )


def adjustments_at_step(ffr, step=None):
    """Calculate the set of adjustments at a specific iteration step.

    Parameters
    ----------
    ffr : FurthestFirstResult
        Output from `furthest_first`.
    step : int, optional
        Which iteration step to calculate adjustments at, by default `None`,
        in which case the final step is calculated.

    Returns
    -------
    np.array
        The set of adjustments at the requested iteration `step`.
    """
    if step is None:
        step = ffr.niter

    def scan_adj(adjustments, x):
        ix = np.array(ffr.steps.index)[x]
        adjustments = adjustments.at[ix].add(np.array(ffr.steps.adjustment)[x])
        return adjustments, x

    adjustments = np.zeros(ffr.adjustments.size)
    adjustments = scan(scan_adj, adjustments, np.arange(step))[0]
    return onp.array(adjustments)


def uncertainties_at_step(ffr, step=None):
    """Calculate the adjustment uncertainties at a specific iteration step.

    Parameters
    ----------
    ffr : FurthestFirstResult
        Output from `furthest_first`.
    step : int, optional
        Which iteration step to calculate uncertainties at, by default `None`,
        in which case the final step is calculated.

    Returns
    -------
    np.array
        The set of adjustment uncertainties at the requested iteration `step`.
    """
    if step is None:
        step = ffr.niter

    def scan_u(uncertainties, x):
        ix = np.array(ffr.steps.index)[x]
        adjustments_u = uncertainties.at[ix].set(
            np.hypot(uncertainties[ix], np.array(ffr.steps.uncertainty)[x])
        )
        return adjustments_u, x

    uncertainties = np.zeros(ffr.adjustments.size)
    uncertainties = scan(scan_u, uncertainties, np.arange(step))[0]
    # Add final step of adjustment uncertainty
    adjustments = adjustments_at_step(ffr, step=step)
    xovers_adj = adjust_xovers(ffr.xovers_adjusted, -ffr.adjustments)
    xovers_adj = adjust_xovers(xovers_adj, adjustments)
    # ^ remove final adjustments and add adjustments at the step of interest
    offsets = offsets_weighted(xovers_adj, ffr.weights)
    offsets_u = offsets_u = offset_uncertainties(
        xovers_adj, ffr.weights, offsets, ffr.dof, ffr.t_crit
    )
    uncertainties = np.hypot(uncertainties, offsets_u)
    return onp.array(uncertainties)

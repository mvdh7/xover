# %%
import numpy as np

import xover.inversion as xinv


# Universal settings for np.isclose
isclose = dict(rtol=1e-5, atol=1e-8)  # rtol=1e-5, atol=1e-8

# Import test data
xovers = np.loadtxt("tests/data/xovers.txt")
weights = np.loadtxt("tests/data/weights.txt")
dates = np.loadtxt("tests/data/dates.txt")

# Preliminary calculations
offsets = xinv.offsets_weighted(xovers, weights)
dof = xinv.dof_kish(weights)
t_crit = xinv.t_critical(dof, window=0.68)  # window=0.68
swt = xinv.sem_weighted_t(xovers, weights, offsets, dof, t_crit)
unc = xinv.offset_uncertainties(xovers, weights, offsets, dof, t_crit)
ff = xinv.furthest_first(xovers)
ff_weights = xinv.furthest_first(xovers, weights=weights)
trends = xinv.get_trends_wls(xovers, weights, dates)


def test_offsets_weighted():
    assert offsets.dtype == np.dtype("float64")
    assert offsets.shape == (xovers.shape[0],)


def test_dof_kish():
    assert dof.dtype == np.dtype("float64")
    assert dof.shape == (xovers.shape[0],)


def test_t_crit():
    assert t_crit.dtype == np.dtype("float64")
    assert t_crit.shape == (xovers.shape[0],)


def test_sem_weighted_t():
    assert swt.dtype == np.dtype("float64")
    assert swt.shape == (xovers.shape[0],)


def test_offset_uncertainties():
    assert unc.dtype == np.dtype("float64")
    assert unc.shape == (xovers.shape[0],)
    assert not np.isinf(unc).any()
    assert not np.isnan(unc).any()


def test_trends():
    assert isinstance(trends, xinv.TrendsWLS)
    assert (
        trends.intercept.shape
        == trends.slope.shape
        == trends.slope_se.shape
        == (xovers.shape[0],)
    )
    assert np.all(trends.slope[np.isnan(trends.slope_se)] == 0)


def test_furthest_first():
    assert isinstance(ff, xinv.FurthestFirstResult)
    assert isinstance(ff.steps, xinv.FurthestFirstSteps)
    assert ff.adjustments.shape == ff.uncertainties.shape == (xovers.shape[0],)
    assert ff.xovers_adjusted.shape == xovers.shape
    assert (
        ff.steps.index.shape
        == ff.steps.adjustment.shape
        == ff.steps.uncertainty.shape
    )


def test_furthest_first_niter():
    niter = 100
    ff_niter = xinv.furthest_first(xovers, niter=niter)
    # If niter not provided, it should be the height/width of xovers
    assert ff.niter == xovers.shape[0]
    # If niter provided, it should be used
    assert ff_niter.niter == niter
    # Order of adjustment steps should be the same, regardless of niter
    assert (ff.steps.index[:niter] == ff_niter.steps.index).all()
    assert np.allclose(
        ff.steps.adjustment[:niter], ff_niter.steps.adjustment, **isclose
    )
    assert np.allclose(
        ff.steps.uncertainty[:niter], ff_niter.steps.uncertainty, **isclose
    )
    # Number of adjustment steps should match niter
    assert ff.steps.index.shape == (ff.niter,)
    assert ff.steps.adjustment.shape == (ff.niter,)
    assert ff.steps.uncertainty.shape == (ff.niter,)
    assert ff_niter.steps.index.shape == (niter,)
    assert ff_niter.steps.adjustment.shape == (niter,)
    assert ff_niter.steps.uncertainty.shape == (niter,)


def test_furthest_first_weights():
    # If weights not provided, dof should be constant everywhere
    assert np.allclose(ff.dof, ff.dof[0], **isclose)
    # If weights provided, dof should be different, as should the adjustments
    assert np.allclose(ff_weights.dof, dof, **isclose)
    assert not np.allclose(ff.adjustments, ff_weights.adjustments, **isclose)
    assert not np.allclose(
        ff.uncertainties, ff_weights.uncertainties, **isclose
    )


def test_consistent_values():
    # offsets_weighted
    assert np.isclose(offsets[0], -2.38454116, **isclose)
    assert np.isclose(offsets[-1], -0.94992911, **isclose)
    assert np.isclose(offsets.mean(), -0.10216455, **isclose)
    # dof_kish
    assert np.isclose(dof[0], 10.73941413, **isclose)
    assert np.isclose(dof[-1], 9.46370074, **isclose)
    assert np.isclose(dof.mean(), 22.03026303, **isclose)
    # t_crit
    assert np.isclose(t_crit[0], 1.0426798577949656, **isclose)
    assert np.isclose(t_crit[-1], 1.0495219377342695, **isclose)
    assert np.isclose(t_crit.mean(), 1.082709593765716, **isclose)
    # sem_weighted_t
    assert np.isclose(swt[0], 0.481113, **isclose)
    assert np.isclose(swt[-1], 0.69189412, **isclose)
    assert np.isclose(
        np.nanmean(swt[~np.isinf(swt)]), 0.839602451616156, **isclose
    )
    assert np.isinf(swt).sum() == 5
    assert np.isnan(swt).sum() == 16
    # offset_uncertainties
    assert np.isclose(unc[0], 0.481113, **isclose)
    assert np.isclose(unc[-1], 0.69189412, **isclose)
    assert np.isclose(unc.mean(), 1.0070331, **isclose)
    assert np.allclose(unc[np.isnan(swt)], 5.24860945, **isclose)
    assert np.allclose(unc[np.isinf(swt)], 5.24860945, **isclose)
    # furthest_first
    assert np.isclose(ff.adjustments[0], 0.06967276524982273, **isclose)
    assert np.isclose(ff.adjustments[-1], 0.02718525013925123, **isclose)
    assert np.isclose(ff.adjustments.mean(), 0.000661363590579137, **isclose)
    assert np.isclose(ff.uncertainties[0], 0.02902267668787203, **isclose)
    assert np.isclose(ff.uncertainties[-1], 0.026583029081349067, **isclose)
    assert np.isclose(ff.uncertainties.mean(), 0.04340880182066582, **isclose)
    assert np.isclose(ff_weights.adjustments[0], 2.2778438277752033, **isclose)
    assert np.isclose(
        ff_weights.adjustments[-1], 1.5379851124822732, **isclose
    )
    assert np.isclose(
        ff_weights.adjustments.mean(), 0.2321543340009231, **isclose
    )
    assert np.isclose(
        ff_weights.uncertainties[0], 0.6028799028896545, **isclose
    )
    assert np.isclose(
        ff_weights.uncertainties[-1], 0.5450323067626754, **isclose
    )
    assert np.isclose(
        ff_weights.uncertainties.mean(), 0.6905165969831086, **isclose
    )
    # trends_wls
    assert np.isclose(trends.slope[0], 0.03493679817827933, **isclose)
    assert np.isclose(trends.slope[-1], 0.15093132912539575, **isclose)
    assert np.isclose(trends.slope.mean(), 0.07313471170123866, **isclose)
    assert np.isnan(trends.slope_se).sum() == 40


def test_at_step():
    niter = 100
    ff_niter = xinv.furthest_first(xovers, niter=niter)
    adj_niter = xinv.adjustments_at_step(ff, step=niter)
    unc_niter = xinv.uncertainties_at_step(ff, step=niter)
    assert np.allclose(adj_niter, ff_niter.adjustments, **isclose)
    assert np.allclose(unc_niter, ff_niter.uncertainties, **isclose)


# test_offsets_weighted()
# test_dof_kish()
# test_t_crit()
# test_sem_weighted_t()
# test_offset_uncertainties()
# test_trends()
# test_furthest_first()
# test_furthest_first_niter()
# test_furthest_first_weights()
# test_consistent_values()
# test_at_step()

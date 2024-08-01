import numpy as np
from scipy.integrate import simpson
from scipy.special import spherical_jn
from cOpt.object.jlzeros import JLZEROS

from scipy.interpolate import CubicSpline
from scipy.linalg import rq

def _inner_prod(chi1, chi2, r):
    '''
    Inner product of two numerical radial functions.

    Parameters
    ----------
        chi1 : array of float
            First radial function.
        chi2 : array of float
            Second radial function.
        r : array of float
            Radial grid.

    Returns
    -------
        float
            Inner product of chi1 and chi2.

    '''
    return simpson(r**2 * chi1 * chi2, x=r)


def _rad_norm(chi, r):
    '''
    Norm of a radial function.

    Parameters
    ----------
        chi : array of float
            Radial function.
        r : array of float
            Radial grid.

    Returns
    -------
        float
            Norm of chi.

    '''
    return np.sqrt(_inner_prod(chi, chi, r))


def _smooth(r, rcut, sigma):
    '''
    Smoothing function used in the generation of numerical radial functions.

    Parameters
    ----------
        r : array of float
            Radial grid.
        rcut : int or float
            Cutoff radius.
        sigma : float
            Smoothing parameter.

    Returns
    -------
        g : array of float
            Smoothing function on the radial grid.
    
    References
    ----------
        Chen, M., Guo, G. C., & He, L. (2010).
        Systematically improvable optimized atomic basis sets for ab initio calculations.
        Journal of Physics: Condensed Matter, 22(44), 445501.
    
    '''
    g = 1. - np.exp(-0.5*((r-rcut)/sigma)**2) if sigma != 0 else np.ones_like(r)
    g[r >= rcut] = 0.0
    return g


def jl_raw(l, q, r, rcut=None, deriv=0):
    '''
    Truncated spherical Bessel functions and derivatives.

    The q-th rcut-truncated l-th order spherical Bessel function is defined as

                -
                |   spherical_jn(l, JLZEROS[l][q] * r / rcut)   r <= rcut
        f(r) =  |
                |   0                                           r > rcut
                -

    where JLZEROS[l][q] is the q-th positive zero of the l-th order spherical
    Besesl function.

    Parameters
    ----------
        l : int
            Order of the spherical Bessel function.
        q : int
            Wavenumber index. Also the number of nodes.
        r : array of float
            Grid where the function is evaluated.
        rcut : int or float, optional
            Cutoff radius. If not given, the last element of "r" is used.
        deriv : int
            Order of the derivative. 0 for the function itself.

    Returns
    -------
        array of float

    Notes
    -----
    Functions of the same l & rcut but different q are orthogonal in terms of
    _inner_prod, but they are not normalized.

    '''
    if rcut is None:
        rcut = r[-1] if hasattr(r, '__len__') else r

    k = JLZEROS[l][q] / rcut

    def _recur(l, m):
        if m == 0:
            if hasattr(r, '__len__'):
                tmp = spherical_jn(l, k * r)
                tmp[r > rcut] = 0.0
                return tmp
            else:
                return 0.0 if r > rcut else spherical_jn(l, k * r)
        else:
            if l == 0:
                return - k * _recur(1, m-1)
            else:
                return ( l * _recur(l-1, m-1) - (l+1) * _recur(l+1, m-1) ) \
                        * k / (2*l+1)

    return _recur(l, deriv)


def jl_raw_norm(l, q, rcut):
    '''
    Norm of a truncated spherical Bessel function.

    Note
    ----
    The integral

        \int_0^{rcut} [ r * spherical_jn(l, JLZEROS[l][q] * r / rcut) ]**2 dr

    has a simple analytical form. See, e.g.,

    Arfken, Weber and Harris, Mathematical Methods for Physicists, 7th ed., p. 704.

    '''
    return ( rcut**1.5 * np.abs(spherical_jn(l+1, JLZEROS[l][q])) ) / np.sqrt(2)


def jl_reduce(l, n, rcut, from_raw=True):
    '''
    This function returns a transformation matrix from truncated spherical Bessel
    functions (raw or normalized) to smoothly reduced spherical Bessel functions.

    Consider a set of n truncated spherical Bessel functions {f} with cutoff radius
    rcut, this function returns an n-by-(n-1) matrix T such that the transformed basis

                [e1, e2, ..., e_{n-1}] = [f1, f2, ..., fn] @ T

    are orthonormal and have vanishing first and second derivatives at rcut.

    Parameters
    ----------
        l : int
            Order of the spherical Bessel function.
        n : int
            Number of initial truncated spherical Bessel functions. Note that the
            size of the transformed basis is n-1.
        rcut : int or float
            Cutoff radius.
        from_raw : bool
            If False, the truncated functions are assumed to be normalized.

    Notes
    -----
    For a normalized truncated spherical Bessel function, the ratio between its
    first and second derivative at rcut is always a constant (-rcut/2). Therefore,
    in order to have vanishing first and second derivatives, it is sufficient to
    work on the first derivative only.

    '''
    if n == 1:
        return np.zeros((1,0))

    inv_raw_norm = np.array([1.0 / jl_raw_norm(l, q, rcut) for q in range(n)])

    # first derivative of the normalized truncated spherical Bessel function at rcut
    D = np.array([[jl_raw(l, q, rcut, deriv=1) * inv_raw_norm[q] for q in range(n)]])

    # null space of D
    C = np.linalg.svd(D, full_matrices=True)[2].T[:,1:]

    # Instead of a "canonicalization" in terms of the kinetic energy,
    # we choose to maintain the consistency of results w.r.t. different
    # numbers of spherical Bessel functions, i.e., the result from N
    # spherical Bessel functions (which is an N-by-(N-1) matrix) should
    # be identical to the upper-left N-by-(N-1) block of the result from
    # any M (M>N) spherical Bessel functions.
    T = inv_raw_norm.reshape(-1,1) * rq(C)[0] if from_raw else rq(C)[0]

    # make sure the largest-magnitude element in each column is positive
    idx = np.argmax(np.abs(T), axis=0)
    return T * np.sign(T[idx, range(n-1)])


def coeff_reduced2raw(coeff, rcut):
    '''
    Converts the coefficients in the smoothly reduced spherical Bessel basis
    to those w.r.t the (raw) truncated spherical Bessel functions.

    Parameters
    ----------
        coeff: list of list of list of float
            A nested list of basis coefficients organized as coeff[l][zeta][p]
            where l, zeta and p label the angular momentum, zeta number and basis
            index respectively.
        rcut : int or float
            Cutoff radius.

    Returns
    -------
        coeff_raw : list of list of list of float
            A list of coefficients in the (raw) truncated spherical Bessel basis.

    Notes
    -----
    Items of the same l must consist of the same number of basis functions.

    '''
    coeff_basis = [ np.array(coeff_l).T for coeff_l in coeff ]

    # note that (array([]) @ array([])).tolist() would give 0 where we want []
    return [(jl_reduce(l, coeff_l.shape[0] + 1, rcut, True) @ coeff_l).T.tolist()
            if coeff_l.size > 0 else []
            for l, coeff_l in enumerate(coeff_basis)]


def coeff_raw2normalized(coeff, rcut):
    '''
    Converts the coefficients in the raw truncated spherical Bessel basis
    to those w.r.t the normalized truncated spherical Bessel functions.

    Parameters
    ----------
        coeff: list of list of list of float
            A nested list of basis coefficients organized as coeff[l][zeta][p]
            where l, zeta and p label the angular momentum, zeta number and basis
            index respectively.
        rcut : int or float
            Cutoff radius.

    Returns
    -------
        coeff_normalized: list of list of list of float
            A list of coefficients in the normalized truncated spherical Bessel basis.

    '''
    return [[[coeff_lzq * jl_raw_norm(l, q, rcut)
              for q, coeff_lzq in enumerate(coeff_lz)]
             for coeff_lz in coeff_l]
            for l, coeff_l in enumerate(coeff)]


def coeff_normalized2raw(coeff, rcut):
    '''
    Converts the coefficients in the normalized truncated spherical Bessel basis
    to those w.r.t the raw truncated spherical Bessel functions.

    Parameters
    ----------
        coeff: list of list of list of float
            A nested list of basis coefficients organized as coeff[l][zeta][p]
            where l, zeta and p label the angular momentum, zeta number and basis
            index respectively.
        rcut : int or float
            Cutoff radius.

    Returns
    -------
        coeff_raw : list of list of list of float
            A list of coefficients in the (raw) truncated spherical Bessel basis.

    '''
    return [[[coeff_lzq / jl_raw_norm(l, q, rcut)
              for q, coeff_lzq in enumerate(coeff_lz)]
             for coeff_lz in coeff_l]
            for l, coeff_l in enumerate(coeff)]


def build_raw(coeff, rcut, r, sigma=0.0, orth=False, normalize=False):
    '''
    Builds a set of numerical radial functions by linear combinations of
    truncated spherical Bessel functions.

    Parameters
    ----------
        coeff : list of list of list of float
            A nested list of spherical Bessel coefficients organized as coeff[l][zeta][p]
            where l, zeta and p label the angular momentum, zeta number and basis index
            respectively.
        rcut : int or float
            Cutoff radius.
        r : array of float
            Grid where the radial functions are evaluated.
        sigma : float
            Smoothing parameter.
        orth : bool
            Whether to Gram-Schmidt orthogonalize the radial functions. If True,
            the resulting radial functions may not be consistent with the given
            spherical Bessel coefficients.
        normalize : bool
            Whether to normalize the radial functions. If True, the resulting
            radial functions may not be consistent with the given coefficients.

    Returns
    -------
        chi : list of list of array of float
            A nested list of numerical radial functions organized as chi[l][zeta][ir].

    '''

    g = _smooth(r, rcut, sigma)
    chi = [[None for _ in coeff_l] for coeff_l in coeff]

    for l, coeff_l in enumerate(coeff):
        for zeta, coeff_lz in enumerate(coeff_l):
            chi[l][zeta] = sum(coeff_lzq * spherical_jn(l, JLZEROS[l][q]*r/rcut) \
                    for q, coeff_lzq in enumerate(coeff_lz))

            chi[l][zeta] *= g # smooth & truncate

            if orth: # Gram-Schmidt
                chi[l][zeta] -= sum(chi[l][y] * _inner_prod(chi[l][y], chi[l][zeta], r) \
                        for y in range(zeta))

            if normalize:
                chi[l][zeta] /= _rad_norm(chi[l][zeta], r)

    return chi


def build_reduced(coeff, rcut, r, orthonormal=False):
    '''
    Builds a set of numerical radial functions by linear combinations of
    orthonormal end-smoothed mixed spherical Bessel basis.

    Parameters
    ----------
        coeff : list of list of list of float
            A nested list of spherical Bessel coefficients organized as coeff[l][zeta][p]
            where l, zeta and p label the angular momentum, zeta number and basis index
            respectively.
        rcut : int or float
            Cutoff radius.
        r : array of float
            Grid where the radial functions are evaluated.
        orthonormal : bool
            Whether to orthonormalize the radial functions. If True, the resulting radial
            functions may not be consistent with the given coefficients.
    
    Returns
    -------
        chi : list of list of array of float
            A nested list of numerical radial functions organized as chi[l][zeta][ir].

    Notes
    -----
    Items of the same l must consist of the same number of basis functions.

    '''
    if orthonormal:
        coeff = [np.linalg.qr(np.array(coeff_l).T)[0].T.tolist() if coeff_l else [] for coeff_l in coeff]

    return build_raw(coeff_reduced2raw(coeff, rcut), rcut, r, 0.0, False, False)
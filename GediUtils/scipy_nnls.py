#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright © 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright © 2003-2013 SciPy Developers.
All rights reserved.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from scipy.optimize import _nnls
from numpy import asarray_chkfinite, zeros, double
__all__ = ['nnls']

def nnls(A, b, maxiter):
    """
    Solve ``argmin_x || Ax - b ||_2`` for ``x>=0``. This is a wrapper
    for a FORTRAN non-negative least squares solver.

    Parameters
    ----------
    A : ndarray
        Matrix ``A`` as shown above.
    b : ndarray
        Right-hand side vector.
    maxiter: int, optional
        Maximum number of iterations, optional.
        Default is ``3 * A.shape[1]``.

    Returns
    -------
    x : ndarray
        Solution vector.
    rnorm : float
        The residual, ``|| Ax-b ||_2``.

    Notes
    -----
    The FORTRAN code was published in the book below. The algorithm
    is an active set method. It solves the KKT (Karush-Kuhn-Tucker)
    conditions for the non-negative least squares problem.

    References
    ----------
    Lawson C., Hanson R.J., (1987) Solving Least Squares Problems, SIAM

    """


    A, b = map(asarray_chkfinite, (A, b))


    if len(A.shape) != 2:
        raise ValueError("expected matrix")
    if len(b.shape) != 1:
        raise ValueError("expected vector")


    m, n = A.shape


    if m != b.shape[0]:
        raise ValueError("incompatible dimensions")


    maxiter = -1 if maxiter is None else int(maxiter)


    w = zeros((n,), dtype=double)
    zz = zeros((m,), dtype=double)
    index = zeros((n,), dtype=int)


    x, rnorm, mode = _nnls.nnls(A, m, n, b, w, zz, index, maxiter)
    if mode != 1:
        raise RuntimeError("too many iterations")


    return x, rnorm
"""Scalar weak-guidance step-index LP modes; SI units, infinite cladding."""
import numpy as np
from scipy.special import jv, jvp, kv, kvp
from scipy.optimize import brentq


def solve(radius=5e-6, wavelength=1.55e-6, na=.26, n_core=1.45):
    if min(radius, wavelength, na, n_core) <= 0 or na >= n_core:
        raise ValueError('Require positive parameters and NA < n_core.')
    n_clad = np.sqrt(n_core**2-na**2)
    v = 2*np.pi*radius*na/wavelength
    if not .05 <= v <= 30:
        raise ValueError('This educational scan is limited to 0.05 <= V <= 30.')
    modes = []
    for ell in range(int(np.ceil(v))+1):
        # Cross-multiplied boundary condition avoids poles at zeros of J_l.
        def residual(u):
            w = np.sqrt(v*v-u*u)
            return u*jvp(ell, u)*kv(ell, w)-w*kvp(ell, w)*jv(ell, u)
        grid = np.linspace(1e-7, v*(1-1e-10), 6000)
        values = residual(grid)
        roots = []
        for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if fa*fb < 0:
                root = brentq(residual, a, b, xtol=1e-13)
                if not roots or abs(root-roots[-1]) > 1e-6:
                    roots.append(root)
        for m, u in enumerate(roots, 1):
            neff = np.sqrt(n_core*n_core-(u*wavelength/(2*np.pi*radius))**2)
            modes.append(dict(mode=f'LP{ell}{m}', l=ell, m=m, u=float(u),
                              neff=float(neff), spatial_degeneracy=1 if ell == 0 else 2))
    return dict(V=float(v), n_core=n_core, n_clad=float(n_clad),
                modes=sorted(modes, key=lambda m: -m['neff']))


if __name__ == '__main__':
    import json
    print(json.dumps(solve(), indent=2))

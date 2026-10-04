"""Finite-grid von Karman phase screens and PSFs. No subharmonic correction."""
import numpy as np


def phase_screen(n=128, length=2., r0=.1, outer_scale=27., seed=42):
    if n < 4 or min(length, r0, outer_scale) <= 0:
        raise ValueError('Positive physical sizes and n >= 4 required.')
    rng = np.random.default_rng(seed)
    f = np.fft.fftfreq(n, d=length/n)
    fx, fy = np.meshgrid(f, f)
    psd = .023*r0**(-5/3)*(fx*fx+fy*fy+outer_scale**-2)**(-11/6)
    psd[0, 0] = 0  # remove piston
    spectrum = (rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*np.sqrt(psd)/length
    phase = np.fft.ifft2(spectrum).real*n*n
    return phase-phase.mean()


def pupil(n=128, obstruction=0.):
    if n < 4 or not 0 <= obstruction < 1:
        raise ValueError('Invalid pupil dimensions or obstruction.')
    y, x = np.indices((n,n))-(n-1)/2
    r = np.hypot(x,y)/(n/4)
    return ((r <= 1)&(r >= obstruction)).astype(float)


def psf(aperture, phase):
    aperture, phase = np.asarray(aperture), np.asarray(phase)
    if aperture.shape != phase.shape or not aperture.any():
        raise ValueError('Matching shapes and a nonempty aperture required.')
    intensity = np.abs(np.fft.fftshift(np.fft.fft2(aperture*np.exp(1j*phase))))**2
    return intensity/intensity.sum()


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    p = pupil(); phase = phase_screen()
    fig, axes = plt.subplots(1, 3, figsize=(12,4))
    for ax, data, title in zip(axes, [phase, psf(p, np.zeros_like(p)), psf(p,phase)], ['Phase (rad)', 'Ideal PSF', 'Turbulent PSF']):
        ax.imshow(data); ax.set_title(title); ax.axis('off')
    fig.tight_layout(); fig.savefig('adaptive_optics.png', dpi=150)

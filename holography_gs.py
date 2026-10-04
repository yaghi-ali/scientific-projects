"""Gerchberg-Saxton phase-only hologram, reconstructed from the nanotech project."""
import numpy as np


def forward(field):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(field), norm='ortho'))


def inverse(field):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(field), norm='ortho'))


def target_pattern(n=64):
    y, x = np.indices((n, n))
    return (((x-n/2)**2+(y-n/2)**2 < (n/4)**2)
            & ((x-n/2)**2+(y-n/2)**2 > (n/7)**2)).astype(float)


def optimize(target, iterations=100, seed=42):
    target = np.asarray(target, dtype=float)
    if target.ndim != 2 or not np.isfinite(target).all() or np.any(target < 0) or not target.any():
        raise ValueError('Expected a nonnegative, nonzero finite 2D intensity target.')
    if iterations < 1:
        raise ValueError('iterations must be positive')
    # Match total target energy to a unit-amplitude illumination field.
    amplitude = np.sqrt(target * target.size / target.sum())
    rng = np.random.default_rng(seed)
    field = np.exp(1j*rng.uniform(-np.pi, np.pi, target.shape))
    errors = []
    for _ in range(iterations):
        far = forward(field)
        errors.append(float(np.linalg.norm(np.abs(far)-amplitude)/np.linalg.norm(amplitude)))
        field = np.exp(1j*np.angle(inverse(amplitude*np.exp(1j*np.angle(far)))))
    return np.angle(field), np.abs(forward(field))**2, np.array(errors)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    target = target_pattern()
    phase, reconstruction, errors = optimize(target)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, data, name in zip(axes, [target, phase, reconstruction], ['Target intensity', 'Phase (rad)', 'Reconstruction']):
        ax.imshow(data); ax.set_title(name); ax.axis('off')
    fig.tight_layout(); fig.savefig('holography.png', dpi=150)
    print(f'Amplitude error: {errors[0]:.4f} -> {errors[-1]:.4f}')

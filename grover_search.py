"""Exact NumPy state-vector model of Grover search; single marked state."""
import numpy as np


def search(n_qubits=3, target=5, iterations=2):
    if not isinstance(n_qubits, int) or not 1 <= n_qubits <= 16:
        raise ValueError('Use 1 to 16 qubits.')
    n = 2**n_qubits
    if not 0 <= target < n or not isinstance(iterations, int) or iterations < 0:
        raise ValueError('Invalid target or iteration count.')
    state = np.ones(n)/np.sqrt(n)
    for _ in range(iterations):
        state[target] *= -1  # oracle
        state = 2*state.mean()-state  # diffusion
    return state**2


if __name__ == '__main__':
    for n in (2, 3):
        for r in (1, 2, 3):
            print(f'{n} qubits, {r} iterations: P(target)={search(n, 1, r)[1]:.6f}')

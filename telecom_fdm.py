"""Three-channel OOK/FDM teaching example. No laboratory dataset required."""
import numpy as np


def encode(text):
    return np.unpackbits(np.frombuffer(text.encode('utf-8'), dtype=np.uint8))


def synthesize(texts=('Optics  ', 'Python  ', 'Signals '), samples_per_bit=400):
    bits = [encode(t) for t in texts]
    if len({len(b) for b in bits}) != 1 or not bits[0].size:
        raise ValueError('All channels must contain the same nonzero byte count.')
    # Integer carrier periods per bit give orthogonal coherent integration.
    if samples_per_bit < 100:
        raise ValueError('Use at least 100 samples per bit.')
    carriers = np.array([19, 25, 31]) / samples_per_bit
    n = np.arange(len(bits[0]) * samples_per_bit)
    signal = sum(np.repeat(b, samples_per_bit) * np.cos(2*np.pi*f*n)
                 for b, f in zip(bits, carriers))
    return signal, carriers, np.array(bits)


def amplitudes(signal, carriers, samples_per_bit):
    signal = np.asarray(signal, dtype=float)
    if signal.ndim != 1 or not np.isfinite(signal).all():
        raise ValueError('Expected a finite 1D signal.')
    if samples_per_bit < 1 or len(signal) % samples_per_bit:
        raise ValueError('Signal must contain an integer number of bits.')
    blocks = signal.reshape(-1, samples_per_bit)
    n = np.arange(samples_per_bit)
    return np.array([2*np.abs(blocks @ np.exp(-2j*np.pi*f*n))/samples_per_bit
                     for f in carriers])


def decode(bits):
    bits = np.asarray(bits, dtype=np.uint8)
    if bits.ndim != 1 or len(bits) % 8 or not np.isin(bits, [0, 1]).all():
        raise ValueError('Expected complete binary bytes.')
    return np.packbits(bits).tobytes().decode('utf-8')


def experiment(seed=42):
    signal, carriers, reference = synthesize()
    rng = np.random.default_rng(seed)
    results = []
    for snr in (3, 6, 10):
        power = np.var(signal)
        noise = rng.normal(0, np.sqrt(power / 10**(snr/10)), signal.shape)
        bits = (amplitudes(signal+noise, carriers, 400) > .5).astype(np.uint8)
        results.append(dict(snr_db=snr, measured_snr_db=float(10*np.log10(power/np.mean(noise**2))),
                            errors=int(np.sum(bits != reference)), bits=int(bits.size),
                            ber=float(np.mean(bits != reference))))
    return results


if __name__ == '__main__':
    import json
    x, f, b = synthesize()
    print(''.join(decode(channel) for channel in (amplitudes(x, f, 400) > .5)))
    print(json.dumps(experiment(), indent=2))

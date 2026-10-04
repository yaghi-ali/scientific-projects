import unittest
import numpy as np
from telecom_fdm import synthesize, amplitudes, decode, experiment
from fiber_modes import solve
from holography_gs import optimize, target_pattern
from adaptive_optics import phase_screen, pupil, psf
from grover_search import search
from tracking_counter import LineCounter


class ScientificChecks(unittest.TestCase):
    def test_fdm_roundtrip_and_noise(self):
        x, f, bits = synthesize()
        recovered = (amplitudes(x, f, 400) > .5).astype(np.uint8)
        np.testing.assert_array_equal(bits, recovered)
        self.assertEqual(''.join(decode(b) for b in recovered), 'Optics  Python  Signals ')
        self.assertTrue(all(abs(r['snr_db']-r['measured_snr_db']) < .2 for r in experiment()))

    def test_lp_modes_and_single_mode_limit(self):
        result = solve()
        self.assertEqual({m['mode'] for m in result['modes']}, {'LP01','LP11','LP21','LP02','LP31'})
        self.assertTrue(all(result['n_clad'] < m['neff'] < result['n_core'] for m in result['modes']))
        single = solve(radius=2e-6, na=.15)
        self.assertEqual([m['mode'] for m in single['modes']], ['LP01'])

    def test_hologram_energy_and_improvement(self):
        phase, image, errors = optimize(target_pattern(), 60)
        self.assertAlmostEqual(image.sum(), phase.size, places=7)
        self.assertLess(errors[-1], errors[0])
        self.assertLessEqual(np.abs(phase).max(), np.pi)

    def test_phase_scaling_and_psf(self):
        a = phase_screen(n=32, r0=.1); b = phase_screen(n=32, r0=.2)
        np.testing.assert_allclose(b, a*2**(-5/6), atol=1e-10)
        self.assertAlmostEqual(a.mean(), 0, places=10)
        p = pupil(32)
        self.assertAlmostEqual(psf(p,a).sum(), 1)
        self.assertLess(psf(p,a).max(), psf(p,np.zeros_like(p)).max())

    def test_grover_exact_formula(self):
        for n in (2,3):
            for target in range(2**n):
                for r in range(5):
                    p = search(n,target,r)
                    self.assertAlmostEqual(p.sum(),1)
                    self.assertAlmostEqual(p[target], np.sin((2*r+1)*np.arcsin(1/np.sqrt(2**n)))**2)

    def test_counter_hysteresis_reverse_and_expiry(self):
        c = LineCounter(max_gap=3)
        events = []
        for frame,y in enumerate([.2,.49,.51,.49,.6,.51,.49,.4]):
            events.extend(c.update(frame,{1:y}))
        self.assertEqual([e['direction'] for e in events], ['IN','OUT'])
        self.assertEqual(c.update(20,{1:.8}), [])
        with self.assertRaises(ValueError): c.update(19,{})

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError): solve(na=2)
        with self.assertRaises(ValueError): optimize(np.zeros((4,4)))
        with self.assertRaises(ValueError): search(target=100)
        with self.assertRaises(ValueError): amplitudes(np.arange(5), [.1], 4)


if __name__ == '__main__':
    unittest.main()

"""Polarisation-multiplexed time-bin topology, built from generic parts.

Gobby et al. route Alice's two arms onto orthogonal polarisations with a
beam combiner, and Bob's beam splitter sends them into *opposite* arms:

    "Polarising beam combiner and splitter are used so that photons from
     Alice's short arm are directed into Bob's long arm (S-L) and vice
     versa (L-S)."                                    -- Gobby 2004 [1]

Only S-L and L-S exist.  A balanced 50:50 pair would instead produce four
paths, of which S-S and L-L are non-interfering satellites at t = 0 and
t = 2*delay carrying half the launched energy -- energy the detection gate
then throws away.  That mu/2 loss is an artefact of modelling the wrong
interferometer, and these tests pin the difference.

Nothing here is Gobby-specific machinery: the topology is assembled from
`AsymmetricMZI` (with an unbalanced `split_ratio` and its arms taken out
via `recombine=False`), `optics.pbc`, `optics.pbs` and `optics.voa`.

References
----------
[1] Gobby, Yuan & Shields, Appl. Phys. Lett. 84(19), 3762-3764, 2004.
[2] Zeilinger, Am. J. Phys. 49(9), 882-883, 1981 -- coupler conventions.
"""
import numpy as np
import pytest

from src.channel.interferometer import AsymmetricMZI
from src.channel.optics import pbs, pbc, voa

DELAY = 5.8e-9
DT = 2.9e-10          # 20 samples per delay
N = 200
DS = int(DELAY / DT)  # adjacent time bins are DS samples apart
HALF = DS // 2 - 1    # window strictly inside half the bin spacing
CENTRE = 50
SIGMA = 2.5           # HALF/SIGMA = 3.6 sigma -> >99.9% of a bin's energy

SPLIT_RATIO = 1.6                       # reference:encoded, Gobby [1]
KAPPA = 1.0 / (1.0 + SPLIT_RATIO)       # encoded (short) arm share
BALANCE_DB = 10.0 * np.log10(SPLIT_RATIO)

PHI_A = {'X': {0: 0.0, 1: np.pi}, 'Y': {0: np.pi / 2, 1: 3 * np.pi / 2}}
PHI_B = {'X': 0.0, 'Y': np.pi / 2}


def bin_window(peak):
    """Slice isolating one time bin (never straddles its neighbours)."""
    return slice(max(0, peak - HALF), peak + HALF + 1)


@pytest.fixture
def pulse():
    t = np.arange(N, dtype=float)
    E = np.zeros((N, 2), dtype=complex)
    E[:, 0] = np.exp(-0.5 * ((t - CENTRE) / SIGMA) ** 2)
    return E


def polmux_chain(E, phi_a, phi_b, balance=True, phase_arm='short',
                 bob_sign=+1.0):
    """Alice -> PBC -> PBS -> Bob, from generic components.

    Mirrors the wiring in `bb84_time_bin.simulate_bb84_time_bin`.

    `phase_arm` and `bob_sign` are exposed so the coupling between them
    can be tested; production uses ('short', +1).
    """
    enc = AsymmetricMZI(delay=DELAY, mode='encoder', split_ratio=KAPPA,
                        phase_arm=phase_arm)
    dec = AsymmetricMZI(delay=DELAY, mode='decoder')

    E_enc, E_ref = enc.modulate(E, DT, phase=phi_a, recombine=False)
    field = pbc(E_ref[:, 0], E_enc[:, 0])          # reference -> H, encoded -> V

    h, v = pbs(field)                              # Bob routes by polarisation
    arm_ref = np.column_stack([h, np.zeros_like(h)])
    arm_enc = np.column_stack([v, np.zeros_like(v)])
    if balance:
        arm_ref = voa(arm_ref, BALANCE_DB)
    arm_ref = arm_ref * np.exp(bob_sign * 1j * phi_b)  # Bob modulates the reference
    return dec.modulate((arm_ref, arm_enc), DT)


def balanced_chain(E, phi_a, phi_b):
    enc = AsymmetricMZI(delay=DELAY, mode='encoder')
    dec = AsymmetricMZI(delay=DELAY, mode='decoder')
    return dec.modulate(enc.modulate(E, DT, phase=phi_a), DT, phase=phi_b)


class TestNoSatellites:
    """The structural difference: two paths, not four."""

    def test_single_interference_peak(self, pulse):
        E_c, E_d = polmux_chain(pulse, 0.0, 0.0)
        total = np.abs(E_c[:, 0]) ** 2 + np.abs(E_d[:, 0]) ** 2
        peak = CENTRE + DS
        assert total[bin_window(peak)].sum() / total.sum() > 0.999
        # the slots a balanced pair would populate are empty
        for sat in (CENTRE, CENTRE + 2 * DS):
            assert total[bin_window(sat)].sum() / total.sum() < 1e-3

    def test_balanced_pair_does_produce_satellites(self, pulse):
        """Control: the 50:50 chain splits energy into three bins."""
        E_c, E_d = balanced_chain(pulse, 0.0, 0.0)
        total = (np.abs(E_c) ** 2).sum(axis=1) + (np.abs(E_d) ** 2).sum(axis=1)
        central = total[bin_window(CENTRE + DS)].sum() / total.sum()
        sats = sum(total[bin_window(s)].sum() / total.sum()
                   for s in (CENTRE, CENTRE + 2 * DS))
        assert 0.4 < central < 0.6, f"central share {central}"
        assert 0.4 < sats < 0.6, f"satellite share {sats}"


class TestTransmission:
    """Balancing the arms costs 1 - 2*kappa; that is a prediction.

    Tolerance note: the window spans +/-3.6 sigma, so ~3e-4 of a Gaussian
    bin's tail falls outside it and the measured ratios sit a few times
    1e-8 below the exact values.  rtol=1e-5 is far above that windowing
    residual and far below any physically meaningful deviation.
    """

    def test_gated_transmission_is_two_kappa(self, pulse):
        launched = (np.abs(pulse) ** 2).sum()
        E_c, E_d = polmux_chain(pulse, 0.0, 0.0)
        w = bin_window(CENTRE + DS)
        gated = (np.abs(E_c[w]) ** 2).sum() + (np.abs(E_d[w]) ** 2).sum()
        assert np.isclose(gated / launched, 2 * KAPPA, rtol=1e-5)

    def test_balanced_chain_loses_half(self, pulse):
        """The mu/2 artefact the polarisation routing removes."""
        launched = (np.abs(pulse) ** 2).sum()
        E_c, E_d = balanced_chain(pulse, 0.0, 0.0)
        w = bin_window(CENTRE + DS)
        gated = (np.abs(E_c[w]) ** 2).sum() + (np.abs(E_d[w]) ** 2).sum()
        assert np.isclose(gated / launched, 0.5, rtol=1e-5)

    def test_prediction_matches_gobby_stated_visibilities(self):
        """Geometry 0.769 vs 0.793 implied by the paper's own fringe data.

        Inverting V = S/(S + 2*P_e) against Gobby's stated visibilities
        gives mu_eff/mu = 0.793.  The geometry predicts 2*kappa with no
        free parameter.  Agreement to 3% is the cross-check that this is
        not a fitted correction -- if a future change makes these diverge,
        something has been tuned.
        """
        assert abs(2 * KAPPA - 0.793) / 0.793 < 0.05


class TestBB84Mapping:
    """The protocol mapping. Absent before, and two real bugs slipped past.

    Twenty tests that swept phase continuously or checked energy all
    passed while (a) the coupler used the imaginary convention, making
    fringes go as sin(phi) so the {0, pi} encoding landed on the zeros,
    and (b) Bob's phase sat on the encoded arm, giving phi_A + phi_B and
    inverting the Y basis.  Both produced ~50% QBER end to end.
    """

    @pytest.mark.parametrize("basis", ['X', 'Y'])
    @pytest.mark.parametrize("bit", [0, 1])
    def test_matched_basis_maps_bit_to_port(self, pulse, basis, bit):
        E_c, E_d = polmux_chain(pulse, PHI_A[basis][bit], PHI_B[basis])
        p_c = (np.abs(E_c) ** 2).sum()
        p_d = (np.abs(E_d) ** 2).sum()
        read = 0 if p_c > p_d else 1
        assert read == bit, (
            f"basis {basis} bit {bit}: P_c={p_c:.4f} P_d={p_d:.4f} -> {read}")

    @pytest.mark.parametrize("basis", ['X', 'Y'])
    @pytest.mark.parametrize("bit", [0, 1])
    def test_matched_basis_is_deterministic(self, pulse, basis, bit):
        """Equalised arms give V = 1, so the wrong port is dark."""
        E_c, E_d = polmux_chain(pulse, PHI_A[basis][bit], PHI_B[basis])
        p = sorted([(np.abs(E_c) ** 2).sum(), (np.abs(E_d) ** 2).sum()])
        assert p[0] / p[1] < 1e-9, "wrong port should be extinguished"

    @pytest.mark.parametrize("a_basis,b_basis", [('X', 'Y'), ('Y', 'X')])
    @pytest.mark.parametrize("bit", [0, 1])
    def test_mismatched_basis_carries_no_information(self, pulse, a_basis,
                                                     b_basis, bit):
        E_c, E_d = polmux_chain(pulse, PHI_A[a_basis][bit], PHI_B[b_basis])
        p_c = (np.abs(E_c) ** 2).sum()
        p_d = (np.abs(E_d) ** 2).sum()
        assert abs(p_c - p_d) / (p_c + p_d) < 0.02


class TestPhaseArmCoupling:
    """Which arm carries the modulator is coupled to Bob's sign.

    Gobby encodes on Alice's *short* arm, so phi_A rides the encoded path
    and phi_B the reference path, giving a relative phase of phi_A - phi_B.
    An earlier version put phi_A on the long (reference) arm together with
    Bob's exp(-i*phi_B), giving phi_B - phi_A: identical intensities,
    because cos is even, but the apparatus modelled backwards.

    Correcting it requires *both* the arm move and Bob's sign flip.  These
    tests pin the pair, and the negative control below is what gives the
    equivalence check its teeth.
    """

    ALL = [(b, k, bb) for b in ('X', 'Y') for k in (0, 1) for bb in ('X', 'Y')]

    def _table(self, pulse, phase_arm, bob_sign):
        out = {}
        for a_basis, bit, b_basis in self.ALL:
            E_c, E_d = polmux_chain(pulse, PHI_A[a_basis][bit],
                                    PHI_B[b_basis], phase_arm=phase_arm,
                                    bob_sign=bob_sign)
            out[(a_basis, bit, b_basis)] = ((np.abs(E_c) ** 2).sum(),
                                            (np.abs(E_d) ** 2).sum())
        return out

    def test_both_changes_are_intensity_equivalent(self, pulse):
        """('long', -1) and ('short', +1) must agree to floating point."""
        old = self._table(pulse, 'long', -1.0)
        new = self._table(pulse, 'short', +1.0)
        for key in old:
            for i in (0, 1):
                a, b = old[key][i], new[key][i]
                assert abs(a - b) <= 1e-12 * max(abs(a), abs(b), 1e-300), key

    def test_arm_move_without_sign_flip_inverts_Y_only(self, pulse):
        """Negative control: a check that cannot fail proves nothing.

        Moving the phase without flipping Bob's sign gives phi_A + phi_B.
        X is unaffected (0 and pi are their own negatives mod 2pi); Y
        inverts.  If this ever stops failing, the equivalence test above
        has lost its sensitivity.
        """
        old = self._table(pulse, 'long', -1.0)
        broken = self._table(pulse, 'short', -1.0)

        for key in [k for k in self.ALL if k[0] == 'X' and k[2] == 'X']:
            for i in (0, 1):
                a, b = old[key][i], broken[key][i]
                assert abs(a - b) <= 1e-12 * max(abs(a), abs(b), 1e-300), \
                    f"X basis should be unaffected, {key} moved"

        for key in [('Y', 0, 'Y'), ('Y', 1, 'Y')]:
            a_c, a_d = old[key]
            b_c, b_d = broken[key]
            assert abs(a_c - b_c) > 0.5 * max(a_c, b_c), \
                f"Y basis must invert without the sign flip, {key} did not"


class TestVisibility:

    def test_equalised_arms_give_unit_visibility(self, pulse):
        p = [( np.abs(polmux_chain(pulse, phi, 0.0)[0]) ** 2).sum()
             for phi in np.linspace(0, 2 * np.pi, 201)]
        p = np.array(p)
        V = (p.max() - p.min()) / (p.max() + p.min())
        assert np.isclose(V, 1.0, atol=1e-9)

    def test_unequalised_arms_cap_visibility_at_0973(self, pulse):
        """Skipping the VOA caps V at 2*sqrt(r)/(1+r) -- which the paper's
        stated >0.99 excludes.  That measurement is what forces the
        balancing, so it is worth a test rather than a comment."""
        p = [(np.abs(polmux_chain(pulse, phi, 0.0, balance=False)[0]) ** 2).sum()
             for phi in np.linspace(0, 2 * np.pi, 401)]
        p = np.array(p)
        V = (p.max() - p.min()) / (p.max() + p.min())
        expected = 2 * np.sqrt(SPLIT_RATIO) / (1 + SPLIT_RATIO)
        assert np.isclose(V, expected, atol=1e-3)
        assert V < 0.99


class TestBiasErrorIsNotGobbyShaped:
    """Generality check: a bias-error parameter must be usable by a
    protocol that has nothing to do with Gobby.

    The risk this guards against is overfitting the simulator to one
    paper -- adding an impairment that only makes sense in the chain it
    was measured for.  Imperfect modulator bias is not that: it is
    universal to phase-modulated QKD, and `PhaseModulator` is already
    constructed by five modules of which the Gobby chain is not one.
    Duplinskiy is polarisation-encoded, a different protocol entirely.
    """

    def test_duplinskiy_accepts_a_bias_error(self):
        from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy
        r = simulate_bb84_duplinskiy(2000, fiber_length=10,
                                     bias_offset_v=0.45, seed=11)
        assert r['n_total'] == 2000

    def test_duplinskiy_default_is_inert(self):
        """Zero bias must reproduce the protocol bit-for-bit."""
        from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy
        a = simulate_bb84_duplinskiy(2000, fiber_length=10, seed=11)
        b = simulate_bb84_duplinskiy(2000, fiber_length=10,
                                     bias_offset_v=0.0, seed=11)
        assert a['qber'] == b['qber']
        assert a['n_sifted'] == b['n_sifted']

    def test_bias_error_raises_qber_in_a_non_gobby_protocol(self):
        """A half-V_pi bias is a gross error and must show up as one,
        confirming the parameter is wired to the physics and not merely
        accepted and ignored."""
        from src.protocols.bb84_duplinskiy import simulate_bb84_duplinskiy
        from src.channel.phase_modulator import PhaseModulator
        vpi = PhaseModulator(crystal_cut='X').Vpi
        clean = simulate_bb84_duplinskiy(20000, fiber_length=0, seed=3)
        biased = simulate_bb84_duplinskiy(20000, fiber_length=0,
                                          bias_offset_v=vpi / 2, seed=3)
        assert biased['qber'] > clean['qber']


class TestDriftClockIsNotTheDetectorClock:
    """`run_duration` decouples drift from the pulse budget.

    Drift used to advance as `pulse_idx / repetition_rate`, so asking for
    tighter error bars silently lengthened the simulated *experiment*.  In
    the Gobby sweep the 122 km point needs 1e9 pulses = 500 s at 2 MHz,
    which accumulated 25 deg of drift against the 6 deg of the paper's
    stated two-minute transfer and put the QBER at 13.52% against a stated
    8.9%.  More pulses must mean a better estimate of the same experiment.

    The detector clock is deliberately left alone: dead time and
    afterpulsing are defined against real elapsed time.
    """

    RATE = np.radians(0.05)

    def _sim(self, **kw):
        from src.protocols.bb84_time_bin import simulate_bb84_time_bin
        base = dict(num_bits=40000, fiber_length=0, alpha_dB=0.2,
                    repetition_rate=2e6, mu=0.1, spad_eta=0.045,
                    gate_width=3.5e-9, pulse_width=80e-12, seed=42,
                    phase_drift_rad_s=self.RATE)
        base.update(kw)
        return simulate_bb84_time_bin(**base)

    def test_none_reproduces_the_pulse_budget_clock(self):
        """Default must be bit-identical to the historical behaviour."""
        a = self._sim(run_duration=None)
        b = self._sim()
        assert (a['qber'], a['n_sifted']) == (b['qber'], b['n_sifted'])

    def test_accumulated_drift_is_independent_of_pulse_count(self):
        """The invariant that was missing: a fixed duration must give a
        fixed final drift however many pulses sample it."""
        d = AsymmetricMZI(delay=5.8e-9, mode='decoder',
                          phase_drift_rad_s=self.RATE)
        for n in (1000, 60_000, 3_000_000, 1_000_000_000):
            scale = 120.0 / (n - 1)
            final = d.arm_phase_offset((n - 1) * scale)
            assert np.degrees(final) == pytest.approx(6.0, abs=1e-9)

    def test_expected_qber_does_not_drift_with_pulse_budget(self):
        """Doubling the pulse count at fixed duration must not move the
        expectation -- it should only tighten the error bar."""
        import math
        out = []
        for n in (600_000, 1_200_000):
            r = self._sim(num_bits=n, run_duration=120.0)
            q, k = r['qber'], r['n_sifted']
            out.append((q, math.sqrt(max(q * (1 - q), 1e-12) / k)))
        (q1, s1), (q2, s2) = out
        assert abs(q1 - q2) < 3.0 * math.hypot(s1, s2)

    def test_a_longer_declared_run_gives_more_drift(self):
        """Sanity in the other direction: duration must still matter."""
        short = self._sim(num_bits=200_000, run_duration=1.0)
        long_ = self._sim(num_bits=200_000, run_duration=6000.0)
        assert long_['qber'] > short['qber']


class TestGobbyBiasIsSolvedJointlyWithDrift:
    """The stated 3.3% is an aggregate of BOTH named mechanisms.

    Gobby attribute it to "slight inaccuracies of the phase modulator
    biases, **as well as phase drift during the experiment**".  Taking
    `arccos(1 - 2*e_mod)` gives the bias-ONLY value, 20.93 deg; applying
    drift on top of that counts the drift twice, which put the 0 km floor
    at 4.30% once the run duration was set to the paper's own 120 s.
    """

    def test_bias_plus_drift_reproduces_the_stated_aggregate(self):
        import math
        G = _validate_gobby()
        d0, r, T = G.PHASE_ERROR_RAD, G.PHASE_DRIFT_RAD_S, G.KEY_TRANSFER_S
        D = r * T
        mean = 0.5 - (math.sin(d0 + D) - math.sin(d0)) / (2.0 * D)
        assert mean == pytest.approx(G.E_MOD, abs=1e-9)

    def test_bias_is_below_the_bias_only_reading(self):
        """17.86 deg, not 20.93 -- a ramp centred on the naive value."""
        import math
        G = _validate_gobby()
        assert np.degrees(G.PHASE_ERROR_RAD) == pytest.approx(17.864, abs=0.01)
        assert G.PHASE_ERROR_RAD < math.acos(1.0 - 2.0 * G.E_MOD)

    def test_falls_back_to_the_closed_form_without_drift(self):
        import math
        G = _validate_gobby()
        assert G._bias_for_aggregate(G.E_MOD, 0.0, 120.0) == pytest.approx(
            math.acos(1.0 - 2.0 * G.E_MOD), abs=1e-12)


def _validate_gobby():
    import os
    import sys
    d = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..', 'analysis', 'val_gobby'))
    if d not in sys.path:
        sys.path.insert(0, d)
    import validate_gobby
    return validate_gobby

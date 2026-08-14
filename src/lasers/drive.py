"""CW vs gain-switched operation from one rate-equation model.

The only difference between the two modes is the injection-current waveform.
Everything else -- pulse shape, chirp, phase statistics -- emerges from the
same coupled-wave + carrier equations.
"""
import numpy as np
from dfb_laser import DFBLaser, LaserParams


def cw(level):
    return lambda t: level


def gain_switched(i_bias, i_peak, period, width, t_rise=20e-12):
    """Square-ish current pulses from below to above threshold."""
    def drive(t):
        phase = t % period
        if phase < t_rise:
            return i_bias + (i_peak - i_bias) * phase / t_rise
        if phase < width:
            return i_peak
        if phase < width + t_rise:
            return i_peak - (i_peak - i_bias) * (phase - width) / t_rise
        return i_bias
    return drive


def pulse_phases(res, period, t_skip=2e-9):
    """Optical phase at each pulse peak, and the phase step between pulses."""
    t, E = res.t, res.field_right
    P = np.abs(E) ** 2
    phases, peaks = [], []
    k = int(np.ceil(t_skip / period))
    while (k + 1) * period < t[-1]:
        m = (t >= k * period) & (t < (k + 1) * period)
        if m.sum() > 3:
            i = np.argmax(np.where(m, P, -1))
            if P[i] > 0:
                phases.append(np.angle(E[i]))
                peaks.append(P[i])
        k += 1
    return np.array(phases), np.array(peaks)


if __name__ == "__main__":
    P_BIAS, P_PEAK = 0.040, 0.220
    PERIOD, WIDTH = 500e-12, 150e-12

    p = LaserParams(kappa_i0=50., quarter_wave_shift=True, N_bragg=2.61e18)

    print("=" * 68)
    print("CW MODE  (constant current)")
    las = DFBLaser(p, n_sections=30, seed=3)
    rc = las.simulate(cw(0.120), t_end=12e-9, record_every=4)
    tot = (rc.power_right + rc.power_left) * 1e3
    n = len(tot) // 3
    print("  steady total power   : %.3f mW  (rms ripple %.3f)"
          % (tot[-n:].mean(), tot[-n:].std()))
    print("  carrier density      : %.4f e18 cm^-3" % (rc.carrier_mean[-1] / 1e18))
    nu = rc.instantaneous_frequency()[-n:]
    print("  residual chirp        : %.1f MHz rms" % (nu.std() / 1e6))
    phi = np.unwrap(np.angle(rc.field_right))
    dphi = np.diff(phi[-n:])
    print("  phase step / sample   : mean %+.4f rad, std %.4f rad"
          % (dphi.mean(), dphi.std()))

    print()
    print("=" * 68)
    print("GAIN-SWITCHED MODE  (%.0f -> %.0f mA, %.0f ps pulses at %.2f GHz)"
          % (P_BIAS * 1e3, P_PEAK * 1e3, WIDTH * 1e12, 1e-9 / PERIOD))
    las = DFBLaser(p, n_sections=30, seed=3)
    rg = las.simulate(gain_switched(P_BIAS, P_PEAK, PERIOD, WIDTH),
                      t_end=30e-9, record_every=2)
    tot_g = (rg.power_right + rg.power_left) * 1e3
    print("  peak total power     : %.3f mW" % tot_g.max())
    print("  mean total power     : %.3f mW" % tot_g.mean())
    print("  carrier swing        : %.4f -> %.4f e18 cm^-3"
          % (rg.carrier_mean[len(rg.t) // 2:].min() / 1e18,
             rg.carrier_mean[len(rg.t) // 2:].max() / 1e18))

    nu_g = rg.instantaneous_frequency()
    half = len(nu_g) // 2
    print("  chirp excursion      : %.2f GHz peak-to-peak"
          % ((nu_g[half:].max() - nu_g[half:].min()) / 1e9))

    ph, pk = pulse_phases(rg, PERIOD, t_skip=6e-9)
    print("  pulses analysed      : %d" % len(ph))
    if len(ph) > 4:
        d = np.diff(np.unwrap(ph))
        d = (d + np.pi) % (2 * np.pi) - np.pi
        print("  pulse-to-pulse phase : std %.3f rad  (uniform = %.3f)"
              % (d.std(), 2 * np.pi / np.sqrt(12)))
        r_vec = np.abs(np.mean(np.exp(1j * ph)))
        print("  phase concentration  : |<e^{i phi}>| = %.3f "
              "(0 = fully random, 1 = locked)" % r_vec)

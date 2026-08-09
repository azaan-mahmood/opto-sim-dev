"""Closed-form QBER for the Gobby, Yuan & Shields (2004) link.

QBER(L) = e_mod + (1 - V_fringe(L)) / 2
V_fringe(L) = S(L) / (S(L) + 2 * P_e)
S(L) = mu * T_int * 10**(-alpha * L / 10) * eta_bob

NO PARAMETER IN THIS MODULE IS FITTED
=====================================
This is a rule, not a description of the current values.

No parameter -- in the physical models or in the analytic comparisons --
may be obtained by fitting to data the model exists to reproduce.  The one
exception is a value the source *itself* states as fitted, and then it must
be cited as such.  **Every parameter must carry a citable source.**

This rule was written after being burned by its absence.  An earlier
version of this module set `MU_EFF = 0.0793`, obtained by inverting
V = S/(S + 2*P_e) against Gobby's *measured* fringe visibilities -- fitting
to the very data the model was supposed to predict -- while the docstring
claimed "Every parameter is from the source paper. Nothing is fitted."
That single fitted number produced a phantom discrepancy which cost
substantial effort to chase, and it resolved only when the quantity it was
standing in for was identified: the interferometer transmission.

  mu       = 0.1      photons/clock leaving Alice                       [1]
  r        = 1.6      Alice's reference:encoded intensity ratio         [1]
  T_int    = 2/(1+r)  = 0.7692, interferometer transmission -- DERIVED
                      from r, not fitted.  Alice's 1.6:1 split must be
                      undone before the paths interfere (an unequal pair
                      caps visibility at 2*sqrt(r)/(1+r) = 0.973, which
                      the paper's stated ">99% up to 65 km" excludes).
                      Alice's S-L and Bob's L-S compensate one another --
                      they must, or phase encoding could not survive the
                      PBS -- with the fibre stretcher trimming the
                      residual.  Equalising discards the excess reference
                      power, costing 1 - 2*kappa_A.
  mu_eff   = mu*T_int = 0.0769, the effective signal at the central
                      interference peak.  Every factor is stated; nothing
                      is inverted from a measurement.
  alpha    = 0.2      dB/km @ 1550 nm                                   [1]
  eta_bob  = 0.045    detector QE (0.12) x Bob apparatus (-4.26 dB)     [1]
  P_e      = 8.5e-7   /clock = 3.2e-7 dark + 5.3e-7 stray clock light   [1]
  e_mod    = 0.033    modulation error, Gobby Fig. 3 arrow             [1]

As a cross-check -- *not* as an input -- inverting Gobby's stated
visibilities gives mu_eff/mu = 0.793 against the 0.769 this geometry
predicts: 3% agreement between two independent routes.  If a future change
makes those diverge, something has been tuned.

This is the analytic reference implementation of GOBBY-2 (section 19 of
opto-sim-issues-and-fixes.md).  It exists to be the *prediction*: the
Monte Carlo field chain in analysis/val_gobby/validate_gobby.py must
reproduce it before its comparison against Gobby's measured points means
anything.  See section 19.13(a): a closed form matching four points is
evidence the physics is understood, not evidence the simulator works.

References
----------
[1] Gobby, C., Yuan, Z. L., & Shields, A. J. (2004). Quantum key
    distribution over 122 km of standard telecom fiber. Appl. Phys.
    Lett. 84(19), 3762-3764.
"""
import numpy as np

# --- Gobby's stated parameters -------------------------------------------
# Every value below is either stated in [1] or derived in closed form from
# values stated in [1].  Nothing here may be fitted -- see the module
# docstring.  If you find yourself wanting to adjust one of these to close
# a residual, that is the failure mode this rule exists to prevent.
MU = 0.1           # photons per clock cycle leaving Alice            [1]
SPLIT_RATIO = 1.6  # Alice's reference:encoded intensity ratio        [1]
T_INT = 2.0 / (1.0 + SPLIT_RATIO)   # = 0.7692, DERIVED from SPLIT_RATIO
MU_EFF = MU * T_INT                 # = 0.0769, DERIVED -- never a literal
ALPHA_DB = 0.2     # dB/km @ 1550 nm                                  [1]
ETA_BOB = 0.045    # detector QE (0.12) x Bob apparatus (-4.26 dB)    [1]
P_E = 8.5e-7       # total error probability /clock: dark + stray     [1]
P_E_DARK = 3.2e-7  # dark-count component of P_E (section 19.11)      [1]
E_MOD = 0.033      # modulation error, constant in L (Fig. 3 arrow)   [1]


def signal(dist_km, mu_eff=MU_EFF):
    """Signal click probability per clock cycle, S(L)."""
    return mu_eff * 10.0 ** (-ALPHA_DB * np.asarray(dist_km) / 10.0) * ETA_BOB


def fringe_visibility(dist_km, p_e=P_E, mu_eff=MU_EFF):
    """Erroneous-count fringe visibility, V_fringe(L) = S / (S + 2*P_e)."""
    s = signal(dist_km, mu_eff=mu_eff)
    return s / (s + 2.0 * p_e)


def erroneous_counts(dist_km, p_e=P_E, mu_eff=MU_EFF):
    """Erroneous-count contribution, (1 - V_fringe)/2, in percent."""
    return (1.0 - fringe_visibility(dist_km, p_e=p_e, mu_eff=mu_eff)) / 2.0 * 100.0


def qber(dist_km, p_e=P_E, mu_eff=MU_EFF):
    """Total QBER in percent: e_mod + erroneous counts.

    p_e overrides the background budget (section 19.11's out-of-sample
    test passes p_e = P_E_DARK to remove the stray-light term); mu_eff
    overrides the effective signal (the section 19.3 sensitivity rows).
    """
    return E_MOD * 100.0 + erroneous_counts(dist_km, p_e=p_e, mu_eff=mu_eff)

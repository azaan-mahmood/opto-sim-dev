from .eye import eye_diagram
from .fields import plot_field
from .stokes import compute_stokes_parameters, poincare, cos_sim


def polarimeter(*args, **kwargs):
    """Lazy import - pulls in tkinter only when actually called.

    Enables headless use (CI, containers, servers without a display).
    """
    from .polarimeter import polarimeter as _polarimeter
    return _polarimeter(*args, **kwargs)

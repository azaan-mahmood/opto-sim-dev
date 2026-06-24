import pytest
import numpy as np
import random


def pytest_configure(config):
    """Pin the RNG seed for all tests unless --seed is passed."""
    seed = config.getoption("seed", None)
    if seed is None:
        seed = 42
    random.seed(seed)
    np.random.seed(seed)
    config.seed = seed


def pytest_addoption(parser):
    parser.addoption(
        "--seed", action="store", default=42, type=int,
        help="RNG seed for reproducible tests (default 42)."
    )


@pytest.fixture
def rng_seed(request):
    """Return the seed used for this test run."""
    return request.config.seed

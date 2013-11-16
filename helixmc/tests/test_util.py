from helixmc import util
import numpy as np


def test_params_join():
    data = np.array([[0, 0, 3, 0, 0, 0.2], [0, 0, 1, 0, 0, 0.1]])
    params_new = util.params_join(data)
    expected = np.array([0, 0, 4, 0, 0, 0.3])
    assert np.sum((params_new - expected) ** 2) == 0

from helixmc import util
import numpy as np
from numpy.testing import assert_allclose


def test_params_join():
    data = np.array([[0, 0, 3, 0, 0, 0.2], [0, 0, 1, 0, 0, 0.1]])
    params_new = util.params_join(data)
    expected = np.array([0, 0, 4, 0, 0, 0.3])
    assert_allclose(params_new, expected)


def test_params_coords_data():
    n = 100
    params = np.random.rand(n, 6)
    o, R = util.params2coords(params)
    assert_allclose(params, util.coords2params(o, R))
    dr, frames = util.params2data(params)
    assert_allclose(params, util.data2params(dr, frames))


def test_coords_dr():
    n = 100
    dr = np.random.rand(n, 3)
    coords = util.dr2coords(dr)
    assert_allclose(dr, util.coords2dr(coords))


def test_unitarize():
    n = 100
    theta = np.random.rand(n)
    axis = np.random.rand(n, 3)
    R = util.R_axis(theta, axis)
    assert_allclose(R, util.unitarize(R))


def test_R_axis():
    n = 100
    theta = np.random.rand(n)
    x = np.tile([1, 0, 0], (n, 1))
    y = np.tile([0, 1, 0], (n, 1))
    z = np.tile([0, 0, 1], (n, 1))
    assert_allclose(util.Rx(theta), util.R_axis(theta, x))
    assert_allclose(util.Ry(theta), util.R_axis(theta, y))
    assert_allclose(util.Rz(theta), util.R_axis(theta, z))


def test_read_seq_from_fasta():
    import tempfile
    tf = tempfile.NamedTemporaryFile()
    tf.write('> Here are comments\n')
    tf.write('; Here are comments\n')
    tf.write('# Here are comments\n')
    tf.write('AUCG\n')
    tf.write('GCCA\n')
    tf.flush()
    assert 'AUCGGCCA' == util.read_seq_from_fasta(tf.name)
    tf.close()


def test_writhe_twist():
    n = 100
    dr = np.tile([0, 0, 3.0], (n, 1))
    assert util.writhe_exact(dr) == 0
    assert util.writhe_fuller(dr) == 0

    rb_axis = np.tile([1.0, 0, 0], (n+1, 1))
    assert util.ribbon_twist(dr, rb_axis) == 0

    angle = np.pi * 0.1
    total_angle = angle * n
    for i in xrange(1, n+1):
        rb_axis[i] = util.Rz(angle).dot(rb_axis[i-1])
    assert_allclose(util.ribbon_twist(dr, rb_axis), total_angle)

    dr = np.random.rand(n, 3)
    assert_allclose(util.writhe_exact(dr) % (4 * np.pi),
                    util.writhe_fuller(dr) % (4 * np.pi))


def test_MC_acpt_rej():
    assert util.MC_acpt_rej(10.0, 1.0)
    assert util.MC_acpt_rej(1.0, 1.0)
    assert not util.MC_acpt_rej(-1e10, 1e10)

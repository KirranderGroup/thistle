"""The three builders of 1- and 2-RDMs from determinant expansions (pyxscat.ci): dense (full
active space), bitwise (determinant pairs, Fortran) and sparse (bit-string excitations), against
each other and against PySCF (FCI and selected CI)."""
import numpy as np
import pytest

from pyxscat import ci

pyscf = pytest.importorskip("pyscf")
from pyscf import fci  # noqa: E402
from pyscf.fci import cistring, selected_ci  # noqa: E402

METHODS = ("dense", "bitwise", "sparse")


def _cas(norb, na, nb, density, seed):
    rng = np.random.default_rng(seed)
    sa, sb = cistring.make_strings(range(norb), na), cistring.make_strings(range(norb), nb)
    c = rng.standard_normal((len(sa), len(sb))) * (rng.random((len(sa), len(sb))) < density)
    c.flat[0] = 1.0
    c /= np.linalg.norm(c)
    dets = [(int(a), int(b)) for a in sa for b in sb]
    return c, dets


@pytest.mark.parametrize("norb,na,nb", [(8, 4, 4), (7, 4, 3), (9, 5, 2), (6, 6, 0), (10, 3, 3)])
@pytest.mark.parametrize("method", METHODS)
def test_against_pyscf_fci(norb, na, nb, method):
    c, dets = _cas(norb, na, nb, 0.4, seed=norb + na)
    r1, r2 = fci.direct_spin1.make_rdm12(c, norb, (na, nb))
    keep = np.flatnonzero(c.ravel())                       # only the nonzero determinants are given
    d1, d2 = ci.make_rdm12(norb, na, nb, [dets[i] for i in keep], c.ravel()[keep], method=method)
    assert np.abs(d1 - r1).max() < 1e-13
    assert np.abs(d2 - r2).max() < 1e-13


def _selected(norb, ne, nstr, ndet, seed):
    rng = np.random.default_rng(seed)

    def strings():
        out = set()
        while len(out) < nstr:
            out.add(sum(1 << int(i) for i in rng.choice(norb, ne, replace=False)))
        return np.array(sorted(out), dtype=np.uint64)
    sa, sb = strings(), strings()
    c = np.zeros((nstr, nstr))
    idx = rng.choice(c.size, ndet, replace=False)
    c.flat[idx] = rng.standard_normal(ndet)
    c /= np.linalg.norm(c)
    dets = [(int(sa[i // nstr]), int(sb[i % nstr])) for i in idx]
    return sa, sb, c, idx, dets


@pytest.mark.parametrize("method", ("bitwise", "sparse"))
def test_selected_ci_in_a_large_orbital_space(method):
    """27 orbitals (the full CI space has 8.5e10 determinants), 1500 selected determinants."""
    norb, ne = 27, 6
    sa, sb, c, idx, dets = _selected(norb, ne, 80, 1500, seed=11)
    civ = selected_ci._as_SCIvector(c, (sa.astype(np.int64), sb.astype(np.int64)))
    r1 = selected_ci.make_rdm1(civ, norb, (ne, ne))
    r2 = selected_ci.make_rdm2(civ, norb, (ne, ne))
    d1, d2 = ci.make_rdm12(norb, ne, ne, dets, c.flat[idx], method=method)
    assert np.abs(d1 - r1).max() < 1e-13
    assert np.abs(d2 - r2).max() < 1e-13
    assert ci._choose(norb, ne, ne, len(dets), "auto") != "dense"


def test_64_orbitals_and_sum_rules():
    norb, ne = 64, 3
    _, _, c, idx, dets = _selected(norb, ne, 60, 800, seed=4)
    b1, b2 = ci.make_rdm12(norb, ne, ne, dets, c.flat[idx], method="bitwise")
    s1, s2 = ci.make_rdm12(norb, ne, ne, dets, c.flat[idx], method="sparse")
    assert np.abs(b2 - s2).max() < 1e-13 and np.abs(b1 - s1).max() < 1e-13
    n = 2 * ne
    assert np.trace(b1) == pytest.approx(n, abs=1e-12)
    assert np.einsum("ppqq->", b2) == pytest.approx(n * (n - 1), abs=1e-11)
    assert np.allclose(b2, b2.transpose(2, 3, 0, 1), atol=1e-14)       # dm2[p,q,r,s] = dm2[r,s,p,q]


def test_bitwise_threads_and_repeated_determinants():
    c, dets = _cas(8, 4, 4, 0.5, seed=2)
    v = c.ravel()
    ref = ci.make_rdm12(8, 4, 4, dets, v, method="dense")
    one = ci._rdm12_bitwise(8, dets, v, True, nthreads=1)
    assert np.abs(one[1] - ref[1]).max() < 1e-13
    # the same determinant given twice: coefficients add up (as in the dense builder)
    dets2 = dets + dets[:10]
    v2 = np.concatenate([v, v[:10]])
    for m in METHODS:
        assert np.abs(ci.make_rdm12(8, 4, 4, dets2, v2, method=m)[1]
                      - ci.make_rdm12(8, 4, 4, dets2, v2, method="dense")[1]).max() < 1e-13


def test_sparse_in_several_passes(monkeypatch):
    c, dets = _cas(9, 5, 4, 1.0, seed=5)
    ref = ci.make_rdm12(9, 5, 4, dets, c.ravel(), method="dense")
    monkeypatch.setattr(ci, "_ENTRIES_PER_PASS", 3000)               # about 140 passes
    d1, d2 = ci.make_rdm12(9, 5, 4, dets, c.ravel(), method="sparse")
    assert np.abs(d2 - ref[1]).max() < 1e-13 and np.abs(d1 - ref[0]).max() < 1e-13


def test_dense_in_blocks(monkeypatch):
    c, dets = _cas(9, 5, 4, 1.0, seed=6)
    r1, r2 = fci.direct_spin1.make_rdm12(c, 9, (5, 4))
    monkeypatch.setattr(ci, "_DENSE_BLOCK_BYTES", 8 * 81 * 126 * 7)      # blocks of 7 beta strings
    d1, d2 = ci.make_rdm12(9, 5, 4, dets, c.ravel(), method="dense")
    assert np.abs(d2 - r2).max() < 1e-13 and np.abs(d1 - r1).max() < 1e-13


def test_transition_rdm1_sparse_vs_dense():
    a, dets = _cas(7, 4, 3, 0.6, seed=8)
    b, _ = _cas(7, 4, 3, 0.6, seed=9)
    td = ci.make_transition_rdm1(7, 4, 3, dets, a.ravel(), dets, b.ravel(), method="dense")
    ts = ci.make_transition_rdm1(7, 4, 3, dets, a.ravel(), dets, b.ravel(), method="sparse")
    tp = fci.direct_spin1.trans_rdm1(a, b, 7, (4, 3))
    assert np.abs(td - ts).max() < 1e-14
    assert np.abs(td - tp.T).max() < 1e-13              # the transpose of PySCF's convention


def test_bitwise_beyond_64_orbitals():
    """Multiword bit strings: a CAS(10,8) vector placed on scattered orbitals of a 128-orbital
    space reproduces PySCF's RDMs (atomic accumulation path: 2 GB array)."""
    from pyscf import gto, scf, mcscf, fci
    from pyscf.fci import cistring
    mol = gto.M(atom="N 0 0 0; N 0 0 1.1", basis="cc-pvdz", verbose=0)
    mc = mcscf.CASCI(scf.RHF(mol).run(), 8, 10).run()
    d1, d2 = fci.direct_spin1.make_rdm12(mc.ci, 8, (5, 5))
    sa = cistring.make_strings(range(8), 5)
    pos = [3, 30, 61, 64, 65, 90, 100, 127]
    remap = lambda x: sum(1 << pos[k] for k in range(8) if x >> k & 1)
    dets = [(remap(int(a)), remap(int(b))) for a in sa for b in sa]
    e1, e2 = ci.make_rdm12(128, 5, 5, dets, mc.ci.ravel())
    assert np.abs(e1[np.ix_(pos, pos)] - d1).max() < 1e-12
    assert np.abs(e2[np.ix_(pos, pos, pos, pos)] - d2).max() < 1e-12
    assert abs(np.abs(e2).sum() - np.abs(d2).sum()) < 1e-10

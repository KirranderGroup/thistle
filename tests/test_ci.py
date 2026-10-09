"""Determinant RDMs, CSF expansion and core embedding (pyxscat.ci) against PySCF."""
import numpy as np
import pytest

from pyxscat import ci

fci = pytest.importorskip("pyscf.fci")
from pyscf.fci import cistring  # noqa: E402


@pytest.mark.parametrize("norb,na,nb", [(4, 2, 2), (6, 3, 2), (6, 3, 3), (5, 1, 0), (7, 4, 3)])
def test_rdms_match_pyscf_fci(norb, na, nb):
    rng = np.random.default_rng(norb * 10 + na)
    sa = cistring.make_strings(range(norb), na)
    sb = cistring.make_strings(range(norb), nb)
    c = rng.normal(size=(len(sa), len(sb)))
    c /= np.linalg.norm(c)
    d1, d2 = fci.direct_spin1.make_rdm12(c, norb, (na, nb))
    dets = [(int(a), int(b)) for a in sa for b in sb]
    m1, m2 = ci.make_rdm12(norb, na, nb, dets, c.reshape(-1))
    assert np.allclose(m1, d1, atol=1e-13)
    assert np.allclose(m2, d2, atol=1e-13)


def _s2(norb, na, nb, dets, co):
    c, (sa, ia), (sb, ib) = ci.ci_matrix(norb, na, nb, dets, co)
    pa = cistring.make_strings(range(norb), na)
    pb = cistring.make_strings(range(norb), nb)
    cp = c[np.ix_([ia[int(s)] for s in pa], [ib[int(s)] for s in pb])]
    cp /= np.linalg.norm(cp)
    return fci.spin_op.spin_square0(cp, norb, (na, nb))[0]


@pytest.mark.parametrize("csf", ["ud", "2ud0", "uudd", "udud", "u2udd", "uuddud", "u", "uud", "2u0", "uuu", "uudu"])
def test_csfs_are_normalised_spin_eigenfunctions(csf):
    nopen = csf.count("u") + csf.count("d")
    nel = nopen + 2 * csf.count("2")
    s = (csf.count("u") - csf.count("d")) / 2
    for ms2 in range(nopen % 2, int(2 * s) + 1, 2):
        dets, co = ci.csfs_to_determinants([csf], [1.0], ms2)
        na = (nel + ms2) // 2
        assert np.linalg.norm(co) == pytest.approx(1.0, abs=1e-12)
        assert _s2(len(csf), na, nel - na, dets, co) == pytest.approx(s * (s + 1), abs=1e-12)


def test_csfs_of_one_configuration_are_orthogonal():
    csfs = ["uudd", "udud"]           # the two singlets of four open shells
    vecs = []
    for csf in csfs:
        dets, co = ci.csfs_to_determinants([csf], [1.0], 0)
        vecs.append(dict(zip(dets, co)))
    keys = set(vecs[0]) | set(vecs[1])
    assert abs(sum(vecs[0].get(k, 0) * vecs[1].get(k, 0) for k in keys)) < 1e-12


def test_core_embedding_matches_pyscf_casci():
    from pyscf import gto, scf, mcscf
    mol = gto.M(atom="N 0 0 0; N 0 0 1.12", basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run()
    mc = mcscf.CASCI(mf, 6, 6).run()
    d1a, d2a = fci.direct_spin1.make_rdm12(mc.ci, 6, (3, 3))
    d1, d2 = ci.embed_core(d1a, d2a, mc.ncore, mol.nao)
    c = mc.mo_coeff
    r1, r2 = mcscf.addons.make_rdm12(mc)                      # AO basis
    assert np.allclose(c @ d1 @ c.T, r1, atol=1e-10)
    r2_mo = np.einsum("pqrs,pa,qb,rc,sd->abcd", r2, *(mol.intor("int1e_ovlp") @ c,) * 4, optimize=True)
    assert np.allclose(d2, r2_mo, atol=1e-9)


def test_transition_rdm1_matches_pyscf():
    rng = np.random.default_rng(4)
    norb, na, nb = 5, 3, 2
    sa = cistring.make_strings(range(norb), na)
    sb = cistring.make_strings(range(norb), nb)
    ca = rng.normal(size=(len(sa), len(sb))); ca /= np.linalg.norm(ca)
    cb = rng.normal(size=(len(sa), len(sb))); cb /= np.linalg.norm(cb)
    dets = [(int(a), int(b)) for a in sa for b in sb]
    t = ci.make_transition_rdm1(norb, na, nb, dets, ca.reshape(-1), dets, cb.reshape(-1))
    assert np.allclose(t, fci.direct_spin1.trans_rdm1(ca, cb, norb, (na, nb)).T, atol=1e-13)


def test_interleaved_phase():
    # |1a 1b| is the same in both orders; |1b 2a| (beta below alpha) changes sign
    assert ci.interleaved_phase(0b01, 0b01) == 1.0
    assert ci.interleaved_phase(0b10, 0b01) == -1.0
    assert ci.interleaved_phase(0b01, 0b10) == 1.0

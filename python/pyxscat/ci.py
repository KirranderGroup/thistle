# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Reduced density matrices from configuration-interaction vectors.

Pure numpy and independent of any electronic structure package.  A CI wavefunction is given
as Slater determinants over ``norb`` active orbitals, each as a pair of occupation bit
strings (bit k set = orbital k occupied), with the determinant defined as::

    |alpha string, beta string> = a+_{a1 alpha} a+_{a2 alpha} ... a+_{b1 beta} a+_{b2 beta} ... |0>

(all alpha creation operators first, each set in ascending orbital order; PySCF's convention).

Spin-summed densities follow PySCF:  dm1[p, q] = <E_pq>,  dm2[p, q, r, s] = <p+ r+ s q>
(chemists' notation, trace N(N-1)), with E_pq = sum_sigma a+_{p sigma} a_{q sigma}.
"""
import itertools as it

import numpy as np

_MAX_BYTES = 4e9
_DENSE_BLOCK_BYTES = 5e8                   # work array of one block of the dense builder


# --------------------------------------------------------------------------- determinant RDMs
def _strings(norb, nel):
    out = [sum(1 << k for k in c) for c in it.combinations(range(norb), nel)]
    return np.array(out, dtype=np.int64), {s: i for i, s in enumerate(out)}


def _excitation_maps(norb, strings, index):
    """For each (r, s): arrays (I, J, sign) with  a+_r a_s |I> = sign |J>."""
    maps = {}
    for r in range(norb):
        for s in range(norb):
            src, dst, sg = [], [], []
            for i, st in enumerate(strings):
                if not (st >> s) & 1:
                    continue
                t = st & ~(1 << s)
                if (t >> r) & 1:
                    continue
                t |= 1 << r
                lo, hi = min(r, s), max(r, s)
                between = 0 if r == s else bin(st & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))).count("1")
                src.append(i)
                dst.append(index[t])
                sg.append(-1.0 if between % 2 else 1.0)
            maps[r, s] = (np.array(src, dtype=np.int64), np.array(dst, dtype=np.int64), np.array(sg))
    return maps


def ci_matrix(norb, nalpha, nbeta, dets, coeffs):
    """Arrange determinant coefficients in the alpha x beta string matrix (full active space)."""
    sa, ia = _strings(norb, nalpha)
    sb, ib = _strings(norb, nbeta)
    c = np.zeros((sa.size, sb.size))
    for (a, b), v in zip(dets, coeffs):
        try:
            c[ia[int(a)], ib[int(b)]] += v
        except KeyError:
            raise ValueError("determinant inconsistent with norb/nalpha/nbeta") from None
    return c, (sa, ia), (sb, ib)


def make_rdm12(norb, nalpha, nbeta, dets, coeffs, normalize=True, method="auto"):
    """Spin-summed 1- and 2-particle density matrices of a determinant expansion.

    ``dets`` is a sequence of (alpha_bits, beta_bits) integers, ``coeffs`` the CI coefficients.
    Returns dm1 (norb, norb) and dm2 (norb, norb, norb, norb).

    Three exact builders (they agree to 1e-15)::

        "dense"    dm2[p,q,r,s] = <E_qp Psi|E_rs Psi> - delta_qr dm1[p,s] with the intermediates
                       x_rs = E_rs|Psi> over the complete active space (alpha x beta strings, in
                       blocks); cost ~ norb^4 N_FCI in BLAS.  Best for CASSCF/CASCI vectors.
        "bitwise"  the determinant-pair algorithm of the earlier code (Fortran, OpenMP): pairs
                       screened with popcnt(ieor), Slater-Condon rules for those differing by at
                       most two spin-orbitals; cost ~ N_det^2 (a few instructions per pair),
                       independent of N_FCI.  Best for truncated or selected expansions (MCCI, a
                       print threshold, RAS) in large orbital spaces.  Any number of orbitals
                       (the only builder beyond 64: multiword bit strings).
        "sparse"   the x_rs intermediates only on the determinants reached from the given ones
                       (sparse product, SciPy); cost ~ N_det norb N_el.  Linear in N_det: for
                       very long expansions on few cores.
        "auto"     the cheapest by an estimate of the three costs.
    """
    method = _choose(norb, nalpha, nbeta, len(dets), method)
    if method == "bitwise" or norb > 64:                  # beyond 64 orbitals only the pair loop
        return _rdm12_bitwise(norb, dets, coeffs, normalize)
    if method == "sparse":
        return _rdm12_sparse(norb, nalpha, nbeta, dets, coeffs, normalize)
    return _rdm12_dense(norb, nalpha, nbeta, dets, coeffs, normalize)


def _n_fci(norb, nalpha, nbeta):
    from math import comb
    return comb(norb, nalpha) * comb(norb, nbeta)


def _choose(norb, nalpha, nbeta, ndet, method, nthreads=None):
    if method in ("dense", "sparse", "bitwise"):
        return method
    if method != "auto":
        raise ValueError("method must be 'auto', 'dense', 'bitwise' or 'sparse'")
    if nthreads is None:
        from ._lib import lib
        nthreads = int(lib.pyxscat_max_threads())
    nfci = _n_fci(norb, nalpha, nbeta)
    nexc = nalpha * (norb - nalpha + 1) + nbeta * (norb - nbeta + 1)
    # seconds, measured (x86-64, 2026): ~1 ns per screened pair and thread for the pair loop
    cost = {"sparse": ndet * nexc * 3.4e-7,                                  # generation, sort, product
            "bitwise": 0.5 * ndet * ndet * 1.1e-9 / nthreads + ndet * nexc * 2e-8}
    if nfci * 8 <= _MAX_BYTES / 4:                   # the full CI vector fits comfortably
        cost["dense"] = norb ** 4 * nfci / 2e10 + norb ** 2 * nfci * 2e-8       # GEMM + scatter
    if norb > 64:
        return "bitwise"
    return min(cost, key=cost.get)


def _words(x, nw):
    """A Python-int bit string -> nw little-endian 64-bit words (int64 view)."""
    out = np.zeros(nw, dtype=np.uint64)
    for k in range(nw):
        out[k] = (x >> (64 * k)) & 0xFFFFFFFFFFFFFFFF
    return out.view(np.int64)


def _rdm12_bitwise(norb, dets, coeffs, normalize, nthreads=None):
    import ctypes
    from ._lib import lib
    nw = (norb + 63) // 64
    if nw == 1:
        a, b, c = _determinant_arrays(norb, dets, coeffs)
        a = a.view(np.int64).reshape(-1, 1)
        b = b.view(np.int64).reshape(-1, 1)
    else:                                                   # any number of orbitals: Python ints
        acc = {}
        for (da, db), cv in zip(dets, coeffs):
            acc[int(da), int(db)] = acc.get((int(da), int(db)), 0.0) + float(cv)
        keys = list(acc)
        a = np.array([_words(k[0], nw) for k in keys]).reshape(len(keys), nw)
        b = np.array([_words(k[1], nw) for k in keys]).reshape(len(keys), nw)
        c = np.array([acc[k] for k in keys])
    if normalize:
        c = c / np.linalg.norm(c)
    dm1 = np.zeros((norb, norb))
    dm2 = np.zeros((norb, norb, norb, norb))
    ierr = ctypes.c_int(0)
    lib.pyxscat_rdm12_bitwise(len(c), norb, nw, np.ascontiguousarray(a), np.ascontiguousarray(b),
                              np.ascontiguousarray(c), int(nthreads or 0), dm1, dm2, ctypes.byref(ierr))
    if ierr.value:
        raise RuntimeError(f"pyxscat_rdm12_bitwise failed (ierr={ierr.value})")
    return dm1, dm2


def _rdm12_dense(norb, nalpha, nbeta, dets, coeffs, normalize):
    c, (sa, ia), (sb, ib) = ci_matrix(norb, nalpha, nbeta, dets, coeffs)
    if normalize:
        c /= np.linalg.norm(c)
    na, nb = c.shape
    ma = _excitation_maps(norb, sa, ia)
    mb = _excitation_maps(norb, sb, ib)
    # blocks of beta strings so that x (norb^2 x na x block) stays below ~0.5 GB
    block = int(max(1, min(nb, _DENSE_BLOCK_BYTES // (8 * norb * norb * na))))
    dm1 = np.zeros((norb, norb))
    g = np.zeros((norb * norb, norb * norb))
    for j0 in range(0, nb, block):
        j1 = min(nb, j0 + block)
        cj = c[:, j0:j1]
        x = np.zeros((norb, norb, na, j1 - j0))              # x[r, s] = (E_rs |Psi>) on these columns
        for (r, s), (src, dst, sg) in ma.items():
            if src.size:
                x[r, s][dst, :] += sg[:, None] * cj[src, :]
        for (r, s), (src, dst, sg) in mb.items():
            m = (dst >= j0) & (dst < j1)
            if np.any(m):
                x[r, s][:, dst[m] - j0] += sg[m][None, :] * c[:, src[m]]
        xf = x.reshape(norb * norb, -1)
        dm1 += (xf @ cj.reshape(-1)).reshape(norb, norb)      # <E_pq>
        g += xf @ xf.T                                         # <E_ab Psi | E_rs Psi>
    dm2 = g.reshape(norb, norb, norb, norb).transpose(1, 0, 2, 3).copy()
    for q in range(norb):
        dm2[:, q, q, :] -= dm1
    return dm1, dm2


# --------------------------------------------------------------------------- external-space structure
def make_rdm12_external(norb, nint, dets, coeffs, normalize=True, nthreads=None, max_external=2):
    """RDMs of an expansion in which at most ``max_external`` (2) electrons occupy the
    "external" orbitals ``nint, ..., norb - 1`` (uncontracted MRCI, CISD-like expansions):
    the pair loop is restricted to determinant pairs that can couple, and the block of the
    2RDM with four external indices is formed by matrix products instead.

    Determinants are grouped by their external label (the set of external spin-orbitals).
    Two determinants couple only if the labels differ by at most four spin-orbitals; with
    exactly four the internal parts are identical and only the external-external-external-
    external block of dm2 receives a contribution,

        Gamma[c,a,d,b] = sum_K sum_{sigma tau} T_K^{sigma tau}[c,d] T_K^{sigma tau}[a,b],

    where T_K^{sigma tau}[a,b] is the (signed) coefficient of the determinant with internal
    part K and externals (a sigma, b tau).  So the kernel runs over pairs whose labels differ
    by at most three spin-orbitals (orders of magnitude fewer than all pairs) and skips that
    block, which is then added from the products.  Exact; any number of orbitals.
    """
    import ctypes
    from ._lib import lib
    if max_external != 2:
        raise NotImplementedError("max_external must be 2")
    if not (0 <= nint <= norb):
        raise ValueError("nint must be between 0 and norb")
    nw = (norb + 63) // 64
    acc = {}
    for (da, db), cv in zip(dets, coeffs):
        key = (int(da), int(db))
        acc[key] = acc.get(key, 0.0) + float(cv)
    keys = list(acc)
    c = np.array([acc[k] for k in keys])
    if normalize:
        c = c / np.linalg.norm(c)
    extmask = ((1 << norb) - 1) ^ ((1 << nint) - 1)
    # external label of each determinant: tuple of (orbital, spin), sorted
    labels = []
    for a, b in keys:
        ea, eb = a & extmask, b & extmask
        lab = tuple(sorted([(k, 0) for k in _setbits(ea)] + [(k, 1) for k in _setbits(eb)]))
        if len(lab) > max_external:
            raise ValueError(f"a determinant has {len(lab)} external electrons (> {max_external})")
        labels.append(lab)
    uniq = sorted(set(labels))
    gid = {lab: g for g, lab in enumerate(uniq)}
    order = np.argsort([gid[l] for l in labels], kind="stable")
    keys = [keys[i] for i in order]
    c = c[order]
    labels = [labels[i] for i in order]
    grp = np.array([gid[l] for l in labels], dtype=np.int32)              # 0-based, sorted
    ng = len(uniq)
    gstart = np.searchsorted(grp, np.arange(ng + 1)).astype(np.int32)
    # neighbour groups g' >= g with |L xor L'| <= 3 (4 is the products block, > 4 never couples)
    by_size = {0: [], 1: [], 2: []}
    contains = {}
    for g, lab in enumerate(uniq):
        by_size[len(lab)].append(g)
        for so in lab:
            contains.setdefault(so, []).append(g)
    nbl = []
    nbstart = np.zeros(ng + 1, dtype=np.int32)
    for g, lab in enumerate(uniq):
        if len(lab) == 2:
            cand = set(by_size[0]) | set(by_size[1])
            for so in lab:
                cand.update(contains[so])          # size-2 labels sharing a spin-orbital (and g itself)
        else:                                      # size 0 or 1: everything can couple
            cand = set(range(ng))
        cand = sorted(x for x in cand if x >= g)
        nbstart[g] = len(nbl)
        nbl.extend(cand)
    nbstart[ng] = len(nbl)
    nb = np.array(nbl, dtype=np.int32)
    wa = np.array([_words(k[0], nw) for k in keys]).reshape(len(keys), nw)
    wb = np.array([_words(k[1], nw) for k in keys]).reshape(len(keys), nw)
    dm1 = np.zeros((norb, norb))
    dm2 = np.zeros((norb, norb, norb, norb))
    ierr = ctypes.c_int(0)
    lib.pyxscat_rdm12_bitwise_groups(len(c), norb, nw, np.ascontiguousarray(wa), np.ascontiguousarray(wb),
                                     np.ascontiguousarray(c), nint, ng, gstart + 1, grp + 1, len(nb),
                                     nbstart + 1, nb + 1, int(nthreads or 0), dm1, dm2, ctypes.byref(ierr))
    if ierr.value:
        raise RuntimeError(f"pyxscat_rdm12_bitwise_groups failed (ierr={ierr.value})")
    # ---- the external block from products over the internal parts K
    nv = norb - nint
    if nv > 0:
        internal = {}
        for (a, b), cv, lab in zip(keys, c, labels):
            if len(lab) != 2:
                continue
            ka, kb = a & ~extmask, b & ~extmask
            internal.setdefault((ka, kb), []).append((lab, cv))
        if internal:
            klist = list(internal)
            # T[K, spin case, p, q] over ordered external pairs (p sigma, q tau); T antisymmetric
            t = {case: np.zeros((len(klist), nv, nv)) for case in ((0, 0), (0, 1), (1, 0), (1, 1))}
            for ik, k in enumerate(klist):
                ka, kb = k
                for lab, cv in internal[k]:
                    (p, sp), (q, sq) = lab
                    sign = _sign_two(ka, kb, p, sp, q, sq)       # D = sign a+_{p sp} a+_{q sq} |K>
                    t[sp, sq][ik, p - nint, q - nint] += sign * cv
                    t[sq, sp][ik, q - nint, p - nint] -= sign * cv
            g = dm2[nint:, nint:, nint:, nint:]                # Gamma[c, a, d, b], written in place
            step = max(1, int(5e8 // (8 * nv ** 3)))            # rows c per block: ~0.5 GB temporaries
            for case, m in t.items():
                if not np.any(m):
                    continue
                mm = m.reshape(len(klist), nv * nv)             # [K, (c,d)]
                for c0 in range(0, nv, step):
                    c1 = min(nv, c0 + step)
                    blk = mm[:, c0 * nv:c1 * nv].T @ mm         # [(c,d), (a,b)] for c in the block
                    g[c0:c1] += blk.reshape(c1 - c0, nv, nv, nv).transpose(0, 2, 1, 3)   # (c,d,a,b) -> (c,a,d,b)
    return dm1, dm2


def _setbits(x):
    out = []
    while x:
        low = x & -x
        out.append(low.bit_length() - 1)
        x ^= low
    return out


def _sign_two(ka, kb, p, sp, q, sq):
    """D = sign * a+_{p sp} a+_{q sq} |K> with D and K in canonical (alpha ascending, then beta
    ascending) order; (p, sp) and (q, sq) are not occupied in K."""
    def before(a, b, k, spin):
        if spin == 0:
            n = bin(a & ((1 << k) - 1)).count("1")
        else:
            n = bin(a).count("1") + bin(b & ((1 << k) - 1)).count("1")
        return -1 if n % 2 else 1
    s = before(ka, kb, q, sq)
    if sq == 0:
        ka |= 1 << q
    else:
        kb |= 1 << q
    s *= before(ka, kb, p, sp)
    return s


# --------------------------------------------------------------------------- sparse (bit strings)
_M1, _M2, _M4 = np.uint64(0x5555555555555555), np.uint64(0x3333333333333333), np.uint64(0x0F0F0F0F0F0F0F0F)
_H01 = np.uint64(0x0101010101010101)


def _popcount(x):
    if hasattr(np, "bitwise_count"):                         # NumPy >= 2.0
        return np.bitwise_count(x).astype(np.uint64)
    x = x - ((x >> np.uint64(1)) & _M1)
    x = (x & _M2) + ((x >> np.uint64(2)) & _M2)
    x = (x + (x >> np.uint64(4))) & _M4
    return (x * _H01) >> np.uint64(56)


def _determinant_arrays(norb, dets, coeffs):
    if norb > 64:
        raise ValueError("the sparse RDM builder handles up to 64 active orbitals")
    a = np.array([int(d[0]) for d in dets], dtype=np.uint64)
    b = np.array([int(d[1]) for d in dets], dtype=np.uint64)
    c = np.asarray(coeffs, dtype=float)
    key, inv = np.unique(np.stack([a, b], axis=1), axis=0, return_inverse=True)   # merge repeats
    return key[:, 0].copy(), key[:, 1].copy(), np.bincount(inv.ravel(), weights=c, minlength=len(key))


_GOLDEN = np.uint64(0x9E3779B97F4A7C15)
_ENTRIES_PER_PASS = 2.5e7                  # ~1.5 GB of work arrays per pass of the sparse builder


def _labels(norb, a, b):
    """One sortable label per determinant: a single uint64 when 2 norb <= 64, else a pair."""
    if 2 * norb <= 64:
        return (a << np.uint64(norb)) | b
    return np.ascontiguousarray(np.stack([a, b], axis=1)).view(
        [("a", np.uint64), ("b", np.uint64)]).ravel()


def _part(a, b, npass):
    """Which pass a determinant belongs to (a well-mixed hash of its strings)."""
    if npass == 1:
        return None
    with np.errstate(over="ignore"):
        h = (a * _GOLDEN) ^ (b * np.uint64(0xC2B2AE3D27D4EB4F))
    return (h >> np.uint64(40)) % np.uint64(npass)


def _excitations(norb, a, b, c, npass=1, k=0):
    """Nonzero entries of x_rs = E_rs |Psi> = sum_sigma a+_{r sigma} a_{s sigma} |Psi> whose
    determinant falls in pass k of npass: (row r*norb + s, alpha bits, beta bits, value), in the
    alpha-first convention (a beta excitation crosses the alpha operators twice: no extra sign)."""
    one = np.uint64(1)
    rows, aa, bb, vv = [], [], [], []
    for spin, st in ((0, a), (1, b)):
        occ = [((st >> np.uint64(i)) & one).astype(bool) for i in range(norb)]
        for s in range(norb):
            if not occ[s].any():
                continue
            for r in range(norb):
                if r == s:
                    m = occ[s]
                    new = st[m]
                    sg = np.ones(new.size)
                else:
                    m = occ[s] & ~occ[r]
                    if not m.any():
                        continue
                    old = st[m]
                    new = (old & ~(one << np.uint64(s))) | (one << np.uint64(r))
                    lo, hi = min(r, s), max(r, s)
                    between = np.uint64(((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                    sg = 1.0 - 2.0 * (_popcount(old & between) & one).astype(float)
                na_, nb_ = (new, b[m]) if spin == 0 else (a[m], new)
                val = sg * c[m]
                part = _part(na_, nb_, npass)
                if part is not None:
                    keep = part == np.uint64(k)
                    na_, nb_, val = na_[keep], nb_[keep], val[keep]
                if val.size:
                    rows.append(np.full(val.size, r * norb + s, dtype=np.int64))
                    aa.append(na_)
                    bb.append(nb_)
                    vv.append(val)
    if not rows:
        e = np.zeros(0, dtype=np.uint64)
        return np.zeros(0, dtype=np.int64), e, e, np.zeros(0)
    return np.concatenate(rows), np.concatenate(aa), np.concatenate(bb), np.concatenate(vv)


def _sparse_passes(norb, nalpha, nbeta, dets, coeffs, normalize):
    """Yield (x, labels) per pass: x = E_rs |Psi> as a sparse (norb^2, ncol) matrix restricted
    to one share of the determinants it reaches, labels = its sorted column labels.  The passes
    split the work so that memory stays bounded (~1.5 GB) for any number of determinants."""
    import scipy.sparse as sp
    a, b, c = _determinant_arrays(norb, dets, coeffs)
    if normalize:
        c = c / np.linalg.norm(c)
    nexc = nalpha * (norb - nalpha + 1) + nbeta * (norb - nbeta + 1)
    npass = int(np.ceil(len(c) * nexc / _ENTRIES_PER_PASS))
    keys = _labels(norb, a, b)
    for k in range(npass):
        rows, xa, xb, val = _excitations(norb, a, b, c, npass, k)
        labels, col = np.unique(_labels(norb, xa, xb), return_inverse=True)
        del xa, xb
        x = sp.csr_matrix((val, (rows, col.ravel())), shape=(norb * norb, len(labels)))
        del rows, col, val
        x.sum_duplicates()
        yield x, labels, keys, c


def _project(labels, keys, values):
    """A vector given on determinants ``keys`` as a dense vector over the columns ``labels``."""
    out = np.zeros(len(labels))
    if not len(labels):
        return out
    pos = np.searchsorted(labels, keys)
    ok = pos < len(labels)
    ok[ok] = labels[pos[ok]] == keys[ok]
    out[pos[ok]] = values[ok]
    return out


def _rdm12_sparse(norb, nalpha, nbeta, dets, coeffs, normalize):
    dm1 = np.zeros((norb, norb))
    g = np.zeros((norb * norb, norb * norb))
    for x, labels, keys, c in _sparse_passes(norb, nalpha, nbeta, dets, coeffs, normalize):
        dm1 += (x @ _project(labels, keys, c)).reshape(norb, norb)      # <Psi|E_pq|Psi>
        g += (x @ x.T).toarray()                                         # <E_ab Psi|E_rs Psi>
    dm2 = g.reshape(norb, norb, norb, norb).transpose(1, 0, 2, 3).copy()
    for q in range(norb):
        dm2[:, q, q, :] -= dm1
    return dm1, dm2


def make_transition_rdm1(norb, nalpha, nbeta, dets_a, coeffs_a, dets_b, coeffs_b, normalize=True,
                         method="auto"):
    """Spin-summed transition density matrix t[p, q] = <A|E_pq|B> between two determinant
    expansions in the same orbitals (the transpose of PySCF's ``trans_rdm1`` convention).
    ``method`` as in :func:`make_rdm12`."""
    method = _choose(norb, nalpha, nbeta, max(len(dets_a), len(dets_b)), method)
    if method == "bitwise":                       # (one-particle: the sparse route is linear and fast)
        method = "sparse"
    if method == "sparse":
        a, b, ca = _determinant_arrays(norb, dets_a, coeffs_a)
        if normalize:
            ca = ca / np.linalg.norm(ca)
        keys_a = _labels(norb, a, b)
        t = np.zeros((norb, norb))
        for x, labels, _, _ in _sparse_passes(norb, nalpha, nbeta, dets_b, coeffs_b, normalize):
            t += (x @ _project(labels, keys_a, ca)).reshape(norb, norb)
        return t
    ca, (sa, ia), (sb, ib) = ci_matrix(norb, nalpha, nbeta, dets_a, coeffs_a)
    cb = ci_matrix(norb, nalpha, nbeta, dets_b, coeffs_b)[0]
    if normalize:
        ca /= np.linalg.norm(ca)
        cb /= np.linalg.norm(cb)
    ma = _excitation_maps(norb, sa, ia)
    mb = _excitation_maps(norb, sb, ib)
    t = np.zeros((norb, norb))
    for (r, s), (src, dst, sg) in ma.items():
        if src.size:
            t[r, s] += np.sum(ca[dst, :] * (sg[:, None] * cb[src, :]))
    for (r, s), (src, dst, sg) in mb.items():
        if src.size:
            t[r, s] += np.sum(ca[:, dst] * (sg[None, :] * cb[:, src]))
    return t


# --------------------------------------------------------------------------- core embedding
def factor_active(dm2_act, tol=1e-12):
    """Active 2RDM as a sum of products, Gamma[p,q,r,s] = sum_k w_k V_k[p,q] V_k[r,s].

    Gamma is averaged over its 8 index permutations (the scattering only sees that part) and
    diagonalised as a symmetric (pq) x (rs) matrix; eigenvalues below ``tol`` x the largest are
    dropped.  At most n_act (n_act + 1) / 2 terms.  Returns (w, V) with V of shape (K, na, na)."""
    from .core import symmetrize_dm2
    na = dm2_act.shape[0]
    g = symmetrize_dm2(dm2_act).reshape(na * na, na * na)
    w, v = np.linalg.eigh(0.5 * (g + g.T))
    keep = np.abs(w) > tol * max(np.abs(w).max(), 1e-300)
    return w[keep], v[:, keep].T.reshape(-1, na, na)


def embed_core(dm1_act, dm2_act, ncore, nmo=None, factored=False, tol=1e-12):
    """Add ``ncore`` doubly occupied orbitals in front of an active space; optionally pad with
    empty orbitals up to ``nmo``.

    ``factored=True`` returns the 2RDM as a :class:`pyxscat.core.FactoredDM2` instead of the
    n^4 array: core-core and core-active blocks as products of the core and active 1RDMs
    (3 Coulomb + 3 exchange terms) and the active block from :func:`factor_active`.  Its
    basis transformation to AOs then costs K n_AO^2 memory instead of n_AO^4, and the kernels
    contract it directly (mode 3), as for CCSD."""
    na = dm1_act.shape[0]
    if factored:
        from .core import FactoredDM2
        n = ncore + na if nmo is None else nmo
        a = slice(ncore, ncore + na)
        dc = np.zeros((n, n))
        dc[:ncore, :ncore] = 2.0 * np.eye(ncore)
        da = np.zeros((n, n))
        da[a, a] = dm1_act
        w, v = factor_active(dm2_act, tol)
        va = np.zeros((len(w), n, n))
        va[:, a, a] = v
        ac, bc, ax, bx = [va * w[:, None, None]], [va], [], []
        if ncore:
            # Dc x Dc - 1/2 Dc.Dc  +  Dc x Da + Da x Dc - 1/2 (Dc.Da + Da.Dc): exchange pairs (pr)(qs)
            ac += [dc[None], dc[None], da[None]]
            bc += [dc[None], da[None], dc[None]]
            ax += [-0.5 * dc[None], -0.5 * dc[None], -0.5 * da[None]]
            bx += [dc[None], da[None], dc[None]]
        cat = lambda x: np.concatenate(x) if x else None
        return dc + da, FactoredDM2(cat(ac), cat(bc), cat(ax), cat(bx))
    n = ncore + na if nmo is None else nmo
    a = slice(ncore, ncore + na)
    dc = np.zeros((n, n))
    dc[:ncore, :ncore] = 2.0 * np.eye(ncore)
    da = np.zeros((n, n))
    da[a, a] = dm1_act
    dm1 = dc + da
    dm2 = np.zeros((n, n, n, n))
    dm2 += np.einsum("pq,rs->pqrs", dc, dc) - 0.5 * np.einsum("ps,rq->pqrs", dc, dc)
    dm2 += np.einsum("pq,rs->pqrs", dc, da) + np.einsum("pq,rs->pqrs", da, dc)
    dm2 -= 0.5 * (np.einsum("ps,rq->pqrs", dc, da) + np.einsum("ps,rq->pqrs", da, dc))
    dm2[a, a, a, a] += dm2_act
    return dm1, dm2


# --------------------------------------------------------------------------- CSF -> determinants
def csf_to_determinants(csf, ms2=None):
    """Expand a spin-adapted configuration (GUGA step string such as ``"22ud00"``, as printed by
    Molcas/OpenMolcas) in Slater determinants.

    ``ms2`` is 2 M_S of the determinants (default: the lowest non-negative value).  Returns a
    list of (coefficient, alpha_bits, beta_bits) in the determinant convention of this module.
    The coupling coefficients follow Shavitt's graphical unitary group approach (the convention
    of the Molcas CI printout); the reordering phase into alpha-first order is included.
    """
    csf = csf.strip()
    n = len(csf)
    nopen = csf.count("u") + csf.count("d")
    if ms2 is None:
        ms2 = nopen % 2
    if (nopen - ms2) % 2 or ms2 > nopen:
        raise ValueError(f"M_S = {ms2}/2 impossible for {csf}")
    # Paldus a, b at each level
    a = np.zeros(n, dtype=int)
    b = np.zeros(n, dtype=int)
    aa = bb = 0
    for i, ch in enumerate(csf):
        if ch == "u":
            bb += 1
        elif ch == "d":
            aa += 1
            bb -= 1
        elif ch == "2":
            aa += 1
        elif ch != "0":
            raise ValueError(f"not a CSF step string: {csf!r}")
        a[i], b[i] = aa, bb
    if b.min() < 0:
        raise ValueError(f"invalid spin coupling (S < 0) in {csf!r}")

    n_alpha_open = (nopen + ms2) // 2
    out = []
    for alpha_open in it.combinations(range(nopen), n_alpha_open):
        alpha_open = set(alpha_open)
        coeff = 1.0
        i_s = i_a = i_b = 0
        abits = bbits = 0
        seq = []                                   # spins of the creation operators in orbital order
        for i, ch in enumerate(csf):
            if ch == "0":
                continue
            if ch == "2":
                if b[i] % 2:
                    coeff = -coeff
                i_a += 1
                i_b += 1
                abits |= 1 << i
                bbits |= 1 << i
                seq += [0, 1]
                continue
            is_alpha = i_s in alpha_open
            if ch == "u":
                coeff *= (a[i] + b[i] - (i_b if is_alpha else i_a)) / b[i]
            else:
                coeff *= ((i_b if is_alpha else i_a) - a[i] + 1) / (b[i] + 2)
                if (b[i] % 2 == 0) == is_alpha:
                    coeff = -coeff
            if is_alpha:
                i_a += 1
                abits |= 1 << i
                seq.append(0)
            else:
                i_b += 1
                bbits |= 1 << i
                seq.append(1)
            i_s += 1
        if coeff == 0.0:
            continue
        # phase to bring the orbital-ordered product into alpha-first order
        inv = 0
        nbeta_seen = 0
        for sp in seq:
            if sp == 1:
                nbeta_seen += 1
            else:
                inv += nbeta_seen
        val = np.sign(coeff) * np.sqrt(abs(coeff))
        if inv % 2:
            val = -val
        out.append((val, abits, bbits))
    return out


def interleaved_phase(abits, bbits):
    """Phase that converts a determinant written in orbital order (a+_1alpha a+_1beta
    a+_2alpha ...; Molpro's convention) into the alpha-first order used here."""
    inv = nbeta = 0
    k = 0
    top = max(abits, bbits).bit_length()
    while k < top:
        if (abits >> k) & 1:
            inv += nbeta
        if (bbits >> k) & 1:
            nbeta += 1
        k += 1
    return -1.0 if inv % 2 else 1.0


def csfs_to_determinants(csfs, coeffs, ms2=None):
    """Convert a CSF expansion to a determinant expansion (dets, coeffs)."""
    acc = {}
    for csf, cv in zip(csfs, coeffs):
        for v, ab, bb in csf_to_determinants(csf, ms2):
            acc[ab, bb] = acc.get((ab, bb), 0.0) + v * cv
    dets = list(acc)
    return dets, np.array([acc[d] for d in dets])

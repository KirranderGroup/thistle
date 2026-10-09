! Part of THISTLE (kernel library pyxscat).
! SPDX-License-Identifier: Apache-2.0
! Copyright 2026 Andres Moreno Carrascosa
! Authors: A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander
!-----------------------------------------------------------------------
! PyXSCAT - rotationally averaged elastic scattering by density-weighted pairs
!
!   I(q) = sum_{mu nu la si} A_mu,nu B_la,si < chi_mu chi_nu | K_l(q r12) | chi_la chi_si >
!
! with K_0 = j0 (isotropic) or K_2 = j2 P2 (second Legendre moment about z).  A = B = D gives
! the elastic scattering <|F(q)|^2>; A /= B the coherent mixed (transition) terms (symmetrised).
!
! The pair density of the elastic term is a product, D (x) D, so the density can be folded into
! each primitive-shell pair BEFORE the pairs are combined: every pair becomes one Hermite charge
! distribution W(tuv) at its centre P,
!
!   rho_P(r) = (pi/p)^{3/2} K_ab sum_tuv W(tuv) Lambda_tuv(r - P),
!
! and the elastic term is a sum over pairs of pairs,
!
!   I(q) = sum_{P,Q} pref_P pref_Q exp(-q^2/4p - q^2/4p') sum_m S_m(P - Q) f_m(q, |P - Q|),
!   S_m  = sum_{TUV} U(TUV) sum_{ijk: T+U+V-i-j-k = m} a(T,i) a(U,j) a(V,k) hx^(T-2i) hy^(U-2j) hz^(V-2k),
!   U(TUV) = sum_{tuv + t'u'v' = TUV} W_P(tuv) W_Q(t'u'v') (-1)^(t'+u'+v'),
!   f_m  = (-q^2)^m j_m(qh)/(qh)^m.
!
! Compared with the quartet kernel (which must loop over every Cartesian component of the four
! functions because a general 2RDM does not factorise), the work per pair of pairs is one
! convolution of two Hermite tables and one radial sum.  The radial functions j_m(x)/x^m are
! evaluated stably at every x, so coincident and nearly coincident pair centres need no
! special treatment and lose no digits.
!-----------------------------------------------------------------------
module xs_elastic
    use xs_pairft, only: hermite_pair
    use xs_radial, only: jn_over_xn, hermite_radial_coeffs, upward_threshold
    use iso_fortran_env, only: int64
    implicit none
    integer, parameter, private :: dp = selected_real_kind(15)
    integer, parameter, private :: ikind = selected_int_kind(8)
contains

    subroutine elastic_pairs(nprim, nao, nq, x, y, z, ga, l, m, n, group, contr, coeff, q, &
            da, db, lorder, screen, res, npairs_out, ierr)
        integer(kind=ikind), intent(in) :: nprim, nao, nq, lorder
        real(kind=dp), intent(in) :: x(nprim), y(nprim), z(nprim), ga(nprim), coeff(nprim), q(nq)
        integer(kind=ikind), intent(in) :: l(nprim), m(nprim), n(nprim), group(nprim), contr(nprim)
        real(kind=dp), intent(in) :: da(nao, nao), db(nao, nao), screen
        real(kind=dp), intent(out) :: res(nq), npairs_out(2)
        integer(kind=ikind), intent(out) :: ierr

        integer(kind=ikind) :: ng, g, i, a, b, gi, gj, lmax, lp, ltot, t, u, v, k, ca, cb, np, ip, jp, iq
        integer(kind=ikind) :: ia, ib, ic, mm, la, lb, ntuv
        integer(kind=ikind), allocatable :: gstart(:), gcount(:), lg(:), idx(:, :, :), lpair(:)
        integer(kind=ikind), allocatable :: it(:), iu(:), iv(:)
        real(kind=dp), allocatable :: cmax(:), ex(:, :, :), ey(:, :, :), ez(:, :, :)
        real(kind=dp), allocatable :: pc(:, :), pexp(:), ppref(:), wa(:, :), wb(:, :), wmax(:)
        real(kind=dp), allocatable :: uu(:), s0(:), s3(:), ax(:, :), ay(:, :), az(:, :), arad(:, :)
        real(kind=dp), allocatable :: fr(:), ga2(:), gtab(:, :), qpow(:, :), xq(:), sx(:), cx(:), gm(:, :)
        real(kind=dp) :: rx, rx2
        real(kind=dp), allocatable :: ub(:)
        integer(kind=ikind) :: nlast
        real(kind=dp), allocatable :: xup(:)
        integer(kind=ikind) :: kk
        real(kind=dp) :: pi, p, px, py, pz, kab, pref, cc, dfa, dfb, hx, hy, hz, hh, sgn, bound, w, qq, sym
        real(kind=dp) :: val, cterm, eq2
        integer(kind=int64) :: npp, npp_done

        ierr = 0
        res = 0.0_dp
        npairs_out = 0.0_dp
        if (lorder /= 0 .and. lorder /= 2) then
            ierr = 8; return
        end if
        pi = acos(-1.0_dp)
        ng = maxval(group)
        if (group(1) /= 1) then
            ierr = 5; return
        end if
        do i = 2, nprim
            if (group(i) /= group(i-1) .and. group(i) /= group(i-1) + 1) then
                ierr = 5; return
            end if
        end do
        allocate(gstart(ng), gcount(ng), lg(ng), cmax(ng))
        gcount = 0
        do i = 1, nprim
            g = group(i)
            if (gcount(g) == 0) then
                gstart(g) = i
                lg(g) = l(i) + m(i) + n(i)
            end if
            gcount(g) = gcount(g) + 1
        end do
        do g = 1, ng
            cmax(g) = maxval(abs(coeff(gstart(g):gstart(g)+gcount(g)-1)))
        end do
        lmax = maxval(lg)
        lp = 2*lmax                       ! largest Hermite order of one pair
        ltot = 2*lp + lorder              ! largest order of a pair of pairs
        allocate(idx(0:ltot, 0:ltot, 0:ltot))
        idx = 0
        ntuv = 0
        do k = 0, ltot
            do t = 0, k
                do u = 0, k - t
                    ntuv = ntuv + 1
                    idx(t, u, k - t - u) = ntuv
                end do
            end do
        end do
        allocate(it(ntuv), iu(ntuv), iv(ntuv))
        do k = 0, ltot
            do t = 0, k
                do u = 0, k - t
                    it(idx(t, u, k-t-u)) = t; iu(idx(t, u, k-t-u)) = u; iv(idx(t, u, k-t-u)) = k - t - u
                end do
            end do
        end do

        ! ---- 1. density-weighted Hermite distributions of all significant pairs ----------------
        np = ng * (ng + 1) / 2
        allocate(pc(3, np), pexp(np), ppref(np), lpair(np), wmax(np))
        allocate(wa(ncum(lp), np), wb(ncum(lp), np))
        allocate(ex(0:lp+1, 0:lmax, 0:lmax), ey(0:lp+1, 0:lmax, 0:lmax), ez(0:lp+1, 0:lmax, 0:lmax))
        np = 0
        do gi = 1, ng
            do gj = gi, ng
                a = gstart(gi); b = gstart(gj)
                p = ga(a) + ga(b)
                kab = exp(-ga(a)*ga(b)/p * ((x(a)-x(b))**2 + (y(a)-y(b))**2 + (z(a)-z(b))**2))
                pref = (pi/p)**1.5_dp * kab
                if (pref * cmax(gi) * cmax(gj) < 1.0e-3_dp * screen) cycle
                px = (ga(a)*x(a) + ga(b)*x(b)) / p
                py = (ga(a)*y(a) + ga(b)*y(b)) / p
                pz = (ga(a)*z(a) + ga(b)*z(b)) / p
                call hermite_pair(lg(gi), lg(gj), px - x(a), px - x(b), p, ex(0:lg(gi)+lg(gj)+1, 0:lg(gi), 0:lg(gj)))
                call hermite_pair(lg(gi), lg(gj), py - y(a), py - y(b), p, ey(0:lg(gi)+lg(gj)+1, 0:lg(gi), 0:lg(gj)))
                call hermite_pair(lg(gi), lg(gj), pz - z(a), pz - z(b), p, ez(0:lg(gi)+lg(gj)+1, 0:lg(gi), 0:lg(gj)))
                np = np + 1
                wa(:, np) = 0.0_dp
                wb(:, np) = 0.0_dp
                do a = gstart(gi), gstart(gi) + gcount(gi) - 1
                    do b = gstart(gj), gstart(gj) + gcount(gj) - 1
                        ca = contr(a); cb = contr(b)
                        dfa = da(ca, cb); dfb = db(ca, cb)
                        if (gi /= gj) then
                            dfa = dfa + da(cb, ca); dfb = dfb + db(cb, ca)
                        end if
                        cc = coeff(a) * coeff(b)
                        if (dfa == 0.0_dp .and. dfb == 0.0_dp) cycle
                        do t = 0, l(a) + l(b)
                            do u = 0, m(a) + m(b)
                                do v = 0, n(a) + n(b)
                                    w = cc * ex(t, l(a), l(b)) * ey(u, m(a), m(b)) * ez(v, n(a), n(b))
                                    wa(idx(t, u, v), np) = wa(idx(t, u, v), np) + dfa * w
                                    wb(idx(t, u, v), np) = wb(idx(t, u, v), np) + dfb * w
                                end do
                            end do
                        end do
                    end do
                end do
                wmax(np) = pref * max(maxval(abs(wa(:, np))), maxval(abs(wb(:, np))))
                if (wmax(np) == 0.0_dp) then
                    np = np - 1; cycle
                end if
                pc(:, np) = [px, py, pz]
                pexp(np) = p
                ppref(np) = pref
                lpair(np) = lg(gi) + lg(gj)
            end do
        end do
        deallocate(ex, ey, ez)

        ! ---- 2. pairs of pairs ---------------------------------------------------------------
        ! q tables shared by all threads: pref_P exp(-q^2/4p) and (-q^2)^m
        allocate(gtab(nq, np), qpow(0:ltot, nq))
        do ip = 1, np
            do iq = 1, nq
                gtab(iq, ip) = ppref(ip) * exp(-q(iq)**2 / (4.0_dp * pexp(ip)))
            end do
        end do
        do iq = 1, nq
            qpow(0, iq) = 1.0_dp
            do mm = 1, ltot
                qpow(mm, iq) = -qpow(mm - 1, iq) * q(iq)**2
            end do
        end do
        ! upward recurrence only where it is stable (see xs_radial: upward_threshold)
        allocate(xup(0:ltot))
        do kk = 0, ltot
            xup(kk) = upward_threshold(kk)
        end do
        npp = 0; npp_done = 0
        !$OMP PARALLEL default(shared) private(ip, jp, la, lb, k, ia, ib, ic, mm, hx, hy, hz, &
        !$OMP& hh, sgn, bound, uu, s0, s3, ax, ay, az, arad, fr, ga2, iq, qq, val, sym, cterm, eq2, &
        !$OMP& xq, sx, cx, gm, rx, rx2, kk, ub, nlast) &
        !$OMP& reduction(+:res, npp, npp_done)
        allocate(uu(ntuv), s0(0:ltot), s3(0:ltot), ax(0:ltot, 0:ltot/2), ay(0:ltot, 0:ltot/2), &
                 az(0:ltot, 0:ltot/2), arad(0:ltot, 0:ltot/2), fr(0:ltot), ga2(nq), xq(nq), sx(nq), cx(nq), &
                 gm(0:ltot, nq), ub(0:ltot))
        call hermite_radial_coeffs(ltot, arad)
        !$OMP DO schedule(dynamic, 4)
        do ip = 1, np
            do jp = ip, np
                npp = npp + 1
                sym = 2.0_dp
                if (ip == jp) sym = 1.0_dp
                la = lpair(ip); lb = lpair(jp)
                ! U(TUV) = 1/2 sum [Wa_P Wb_Q + Wb_P Wa_Q](tuv, t'u'v') (-1)^(t'+u'+v')
                k = la + lb
                uu(1:ncum(k)) = 0.0_dp
                do ia = 1, ncum(la)
                    if (wa(ia, ip) == 0.0_dp .and. wb(ia, ip) == 0.0_dp) cycle
                    do ib = 1, ncum(lb)
                        sgn = 1.0_dp
                        if (mod(it(ib) + iu(ib) + iv(ib), 2) == 1) sgn = -1.0_dp
                        ic = idx(it(ia) + it(ib), iu(ia) + iu(ib), iv(ia) + iv(ib))
                        uu(ic) = uu(ic) + 0.5_dp * sgn * (wa(ia, ip) * wb(ib, jp) + wb(ia, ip) * wa(ib, jp))
                    end do
                end do
                ! Screening.  |d^tuv j0(q|h|)| <= q^(t+u+v), so the contribution at q is bounded by
                ! sym pref_P pref_Q exp(-q^2/4p - q^2/4p') sum_m |U_m| q^m (x 2 for j2 P2).  With q in
                ! ascending order, drop the tail of q points where the bound is below screen.
                do mm = 0, k
                    ub(mm) = sum(abs(uu(ncum(mm - 1) + 1:ncum(mm))))
                end do
                nlast = nq
                do while (nlast >= 1)
                    qq = q(nlast)
                    bound = 0.0_dp
                    do mm = k, 0, -1
                        bound = bound * qq + ub(mm)
                    end do
                    bound = bound * sym * gtab(nlast, ip) * gtab(nlast, jp) * real(1 + lorder, dp)
                    if (bound >= screen) exit
                    nlast = nlast - 1
                end do
                if (nlast == 0) cycle
                npp_done = npp_done + 1
                hx = pc(1, ip) - pc(1, jp)
                hy = pc(2, ip) - pc(2, jp)
                hz = pc(3, ip) - pc(3, jp)
                hh = sqrt(hx*hx + hy*hy + hz*hz)
                ! per-direction tables a(T,i) h^(T-2i)
                call powers(ltot, hx, arad, ax)
                call powers(ltot, hy, arad, ay)
                call powers(ltot, hz, arad, az)
                ! S_m for K_0 (and the shifted table for the 3 d^2/dz^2 of Hobson's theorem)
                call radial_sum(k, uu, 0, ax, ay, az, ltot, s0)
                if (lorder == 2) call radial_sum(k, uu, 2, ax, ay, az, ltot, s3)
                ! q dependence: sum_m S_m (-q^2)^m j_m(qh)/(qh)^m
                if (lorder == 0 .and. k == 0) then
                    ! s-type charge distributions (the most common case): S_0 j0(qh)
                    do iq = 1, nlast
                        qq = q(iq) * hh
                        if (qq > 1.0e-3_dp) then
                            val = sin(qq) / qq
                        else
                            val = 1.0_dp - qq*qq / 6.0_dp * (1.0_dp - qq*qq / 20.0_dp)
                        end if
                        res(iq) = res(iq) + sym * gtab(iq, ip) * gtab(iq, jp) * s0(0) * val
                    end do
                    cycle
                end if
                ! general case, vectorised over q: upward recurrence where it is stable (x > K + 1),
                ! the careful scalar routine elsewhere
                kk = k + lorder
                do iq = 1, nlast
                    xq(iq) = q(iq) * hh
                end do
                do iq = 1, nlast
                    sx(iq) = sin(xq(iq))
                    cx(iq) = cos(xq(iq))
                end do
                do iq = 1, nlast
                    if (xq(iq) > xup(kk)) then
                        rx = 1.0_dp / xq(iq)
                        rx2 = rx * rx
                        gm(0, iq) = sx(iq) * rx
                        if (kk >= 1) gm(1, iq) = (gm(0, iq) - cx(iq)) * rx2
                        do mm = 1, kk - 1
                            gm(mm + 1, iq) = (real(2*mm + 1, dp) * gm(mm, iq) - gm(mm - 1, iq)) * rx2
                        end do
                    else
                        call jn_over_xn(xq(iq), kk, gm(0:kk, iq))
                    end if
                end do
                do iq = 1, nlast
                    val = 0.0_dp
                    qq = q(iq)
                    if (lorder == 0) then
                        do mm = 0, k
                            val = val + s0(mm) * qpow(mm, iq) * gm(mm, iq)
                        end do
                    else if (qq > 0.0_dp) then
                        do mm = 0, k + 2
                            if (mm <= k) val = val + qq * qq * s0(mm) * qpow(mm, iq) * gm(mm, iq)
                            val = val + 3.0_dp * s3(mm) * qpow(mm, iq) * gm(mm, iq)
                        end do
                        val = val / (2.0_dp * qq * qq)
                    end if
                    res(iq) = res(iq) + sym * gtab(iq, ip) * gtab(iq, jp) * val
                end do
            end do
        end do
        !$OMP END DO
        deallocate(uu, s0, s3, ax, ay, az, arad, fr, ga2, xq, sx, cx, gm, ub)
        !$OMP END PARALLEL
        deallocate(gtab, qpow, xup)
        npairs_out(1) = real(npp, dp)
        npairs_out(2) = real(npp_done, dp)

    contains

        ! number of (t,u,v) with t+u+v <= kk (they are packed first by order)
        pure integer(kind=ikind) function ncum(kk)
            integer(kind=ikind), intent(in) :: kk
            ncum = (kk + 1) * (kk + 2) * (kk + 3) / 6       ! 0 for kk = -1
        end function ncum

        pure subroutine powers(lt, hc, ar, out)
            integer(kind=ikind), intent(in) :: lt
            real(kind=dp), intent(in) :: hc, ar(0:lt, 0:lt/2)
            real(kind=dp), intent(out) :: out(0:lt, 0:lt/2)
            real(kind=dp) :: hp(0:lt)
            integer(kind=ikind) :: tt, ii
            hp(0) = 1.0_dp
            do tt = 1, lt
                hp(tt) = hp(tt - 1) * hc
            end do
            out = 0.0_dp
            do tt = 0, lt
                do ii = 0, tt / 2
                    out(tt, ii) = ar(tt, ii) * hp(tt - 2*ii)
                end do
            end do
        end subroutine powers

        ! s(m) = sum_TUV U(TUV) sum_ijk ax(T',i) ay(U,j) az(V',k), m = T'+U+V'-i-j-k, with the z order
        ! shifted by dz (0, or 2 for d^2/dz^2)
        ! (the private work arrays of the OpenMP region are passed explicitly: a contained
        ! procedure would otherwise see the host's shared copies)
        subroutine radial_sum(kk, uin, dz, bx, by, bz, lt, s)
            integer(kind=ikind), intent(in) :: kk, lt
            integer, intent(in) :: dz
            real(kind=dp), intent(in) :: uin(:), bx(0:lt, 0:lt/2), by(0:lt, 0:lt/2), bz(0:lt, 0:lt/2)
            real(kind=dp), intent(out) :: s(0:)
            integer(kind=ikind) :: c, tt, uv, vv, i1, j1, k1, ord
            real(kind=dp) :: c1, c2
            s(0:kk + dz) = 0.0_dp
            do c = 1, ncum(kk)
                if (uin(c) == 0.0_dp) cycle
                tt = it(c); uv = iu(c); vv = iv(c) + dz
                do i1 = 0, tt / 2
                    c1 = uin(c) * bx(tt, i1)
                    if (c1 == 0.0_dp) cycle
                    do j1 = 0, uv / 2
                        c2 = c1 * by(uv, j1)
                        if (c2 == 0.0_dp) cycle
                        do k1 = 0, vv / 2
                            ord = tt + uv + vv - i1 - j1 - k1
                            s(ord) = s(ord) + c2 * bz(vv, k1)
                        end do
                    end do
                end do
            end do
        end subroutine radial_sum

    end subroutine elastic_pairs

end module xs_elastic

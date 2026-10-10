! Part of THISTLE (kernel library pyxscat).
! SPDX-License-Identifier: Apache-2.0
! Copyright 2026 Andres Moreno Carrascosa
! Authors: A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander
!-----------------------------------------------------------------------
!> PyXSCAT - electron-nucleus term of electron scattering
!>
!> In the first Born approximation a fast electron is scattered by the full
!> charge distribution, nuclei and electrons:
!>
!>      dsigma/dOmega = 4 / (a0^2 q^4) < | sum_A Z_A exp(iq.R_A) - F(q) |^2 >  + inelastic
!>
!> which, rotationally averaged, needs one quantity beyond X-ray scattering,
!> the electron-nucleus cross term
!>
!>      I_ne(q) = sum_A Z_A sum_mu,nu D_mu,nu int chi_mu chi_nu K_l(q, r - R_A) dr
!>
!> with K_0 = j0(q|h|) (isotropic signal) or K_2 = j2(q|h|) P2(cos theta_h)
!> (second Legendre moment about z, aligned molecules).
!>
!> Each primitive pair is expanded in Hermite Gaussians at the pair centre P:
!>
!>      int Lambda_tuv(r - P) j0(q|r - A|) dr = (pi/p)^{3/2} exp(-q^2/4p) d^{t+u+v}/dh^tuv j0(q|h|),
!>      h = P - A,
!>
!> and the Cartesian derivatives of the radial function j0(q|h|) follow from
!> the McMurchie-Davidson recursion with the auxiliary functions
!>
!>      f_n(h) = (1/h d/dh)^n j0(q h) = (-q^2)^n j_n(qh) / (qh)^n.
!>
!> j2 P2 follows from Hobson's theorem: j2(qh) P2(cos theta) = [3 d^2/dhz^2 + q^2] j0(qh) / (2 q^2).
!> Nothing is singular when a nucleus sits on a pair centre (h = 0).
!-----------------------------------------------------------------------
module xs_nuclear
    use xs_pairft, only: hermite_pair
    use xs_radial, only: jn_over_xn
    implicit none
    integer, parameter, private :: dp = selected_real_kind(15)
    integer, parameter, private :: ikind = selected_int_kind(8)
contains

    !> Hermite derivatives R(k) = d^{t+u+v}/dh_x^t dh_y^u dh_z^v j0(q|h|) for t+u+v <= lmax,
    !> packed by idx(t,u,v).
    pure subroutine j0_derivatives(q, hx, hy, hz, lmax, ldim, idx, ntuv, r)
        real(kind=dp), intent(in) :: q, hx, hy, hz
        integer(kind=ikind), intent(in) :: lmax, ldim, ntuv
        integer(kind=ikind), intent(in) :: idx(0:ldim, 0:ldim, 0:ldim)
        real(kind=dp), intent(out) :: r(ntuv)
        real(kind=dp) :: g(0:lmax), old(ntuv), fn, q2n
        integer(kind=ikind) :: n, t, u, v, k

        call jn_over_xn(q * sqrt(hx*hx + hy*hy + hz*hz), lmax, g)
        old = 0.0_dp
        r = 0.0_dp
        do n = lmax, 0, -1
            q2n = (-q*q)**n
            fn = q2n * g(n)
            r(idx(0, 0, 0)) = fn
            do k = 1, lmax - n
                do t = 0, k
                    do u = 0, k - t
                        v = k - t - u
                        if (t > 0) then
                            r(idx(t, u, v)) = hx * old(idx(t-1, u, v))
                            if (t > 1) r(idx(t, u, v)) = r(idx(t, u, v)) + (t-1) * old(idx(t-2, u, v))
                        else if (u > 0) then
                            r(idx(t, u, v)) = hy * old(idx(t, u-1, v))
                            if (u > 1) r(idx(t, u, v)) = r(idx(t, u, v)) + (u-1) * old(idx(t, u-2, v))
                        else
                            r(idx(t, u, v)) = hz * old(idx(t, u, v-1))
                            if (v > 1) r(idx(t, u, v)) = r(idx(t, u, v)) + (v-1) * old(idx(t, u, v-2))
                        end if
                    end do
                end do
            end do
            if (n > 0) old = r
        end do
    end subroutine j0_derivatives

    !> I_ne(q) = sum_A zat(A) sum_mu,nu dm(mu,nu) int chi_mu chi_nu K_lorder(q, r - R_A) dr
    subroutine electron_nucleus(nprim, nao, nq, x, y, z, ga, l, m, n, group, contr, coeff, q, &
            natom, zat, rat, dm, lorder, screen, res, ierr)
        integer(kind=ikind), intent(in) :: nprim, nao, nq, natom, lorder
        real(kind=dp), intent(in) :: x(nprim), y(nprim), z(nprim), ga(nprim), coeff(nprim), q(nq)
        integer(kind=ikind), intent(in) :: l(nprim), m(nprim), n(nprim), group(nprim), contr(nprim)
        real(kind=dp), intent(in) :: zat(natom), rat(3, natom), dm(nao, nao), screen
        real(kind=dp), intent(out) :: res(nq)
        integer(kind=ikind), intent(out) :: ierr

        integer(kind=ikind) :: ng, g, i, a, b, gi, gj, lmax, lder, iq, t, u, v, ca, cb, iat, ntuv, k
        integer(kind=ikind), allocatable :: gstart(:), gcount(:), lg(:), idx(:, :, :)
        real(kind=dp), allocatable :: cmax(:), ex(:, :, :), ey(:, :, :), ez(:, :, :), rr(:, :), w(:)
        real(kind=dp) :: pi, p, px, py, pz, kab, pref, cc, dmax, dfac, s, gq

        ierr = 0
        res = 0.0_dp
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
        lder = 2*lmax + lorder
        allocate(idx(0:lder, 0:lder, 0:lder))
        idx = 0
        ntuv = 0
        do k = 0, lder
            do t = 0, k
                do u = 0, k - t
                    ntuv = ntuv + 1
                    idx(t, u, k - t - u) = ntuv
                end do
            end do
        end do
        dmax = maxval(abs(dm))

        !$OMP PARALLEL default(shared) private(gi, gj, a, b, p, px, py, pz, kab, pref, ex, ey, ez, &
        !$OMP& rr, w, t, u, v, ca, cb, cc, iq, iat, dfac, s, gq) reduction(+:res)
        allocate(ex(0:2*lmax+1, 0:lmax, 0:lmax), ey(0:2*lmax+1, 0:lmax, 0:lmax), ez(0:2*lmax+1, 0:lmax, 0:lmax))
        allocate(rr(ntuv, nq), w(ntuv))
        !$OMP DO schedule(dynamic)
        do gi = 1, ng
            do gj = gi, ng
                a = gstart(gi); b = gstart(gj)
                p = ga(a) + ga(b)
                kab = exp(-ga(a)*ga(b)/p * ((x(a)-x(b))**2 + (y(a)-y(b))**2 + (z(a)-z(b))**2))
                pref = (pi/p)**1.5_dp * kab
                if (pref * cmax(gi) * cmax(gj) * dmax * (1.0_dp + real(lg(gi) + lg(gj), dp)) < screen) cycle
                px = (ga(a)*x(a) + ga(b)*x(b)) / p
                py = (ga(a)*y(a) + ga(b)*y(b)) / p
                pz = (ga(a)*z(a) + ga(b)*z(b)) / p
                call hermite_pair(lg(gi), lg(gj), px - x(a), px - x(b), p, ex(0:lg(gi)+lg(gj)+1, 0:lg(gi), 0:lg(gj)))
                call hermite_pair(lg(gi), lg(gj), py - y(a), py - y(b), p, ey(0:lg(gi)+lg(gj)+1, 0:lg(gi), 0:lg(gj)))
                call hermite_pair(lg(gi), lg(gj), pz - z(a), pz - z(b), p, ez(0:lg(gi)+lg(gj)+1, 0:lg(gi), 0:lg(gj)))

                ! density-weighted Hermite coefficients of this shell pair: w(tuv)
                w = 0.0_dp
                do a = gstart(gi), gstart(gi) + gcount(gi) - 1
                    do b = gstart(gj), gstart(gj) + gcount(gj) - 1
                        ca = contr(a); cb = contr(b)
                        dfac = dm(ca, cb)
                        if (gi /= gj) dfac = dfac + dm(cb, ca)
                        cc = coeff(a) * coeff(b) * dfac
                        if (cc == 0.0_dp) cycle
                        do t = 0, l(a) + l(b)
                            do u = 0, m(a) + m(b)
                                do v = 0, n(a) + n(b)
                                    w(idx(t, u, v)) = w(idx(t, u, v)) + &
                                            cc * ex(t, l(a), l(b)) * ey(u, m(a), m(b)) * ez(v, n(a), n(b))
                                end do
                            end do
                        end do
                    end do
                end do
                if (maxval(abs(w)) * pref < screen) cycle

                do iat = 1, natom
                    do iq = 1, nq
                        call j0_derivatives(q(iq), px - rat(1, iat), py - rat(2, iat), pz - rat(3, iat), &
                                lg(gi) + lg(gj) + lorder, lder, idx, ntuv, rr(:, iq))
                    end do
                    do iq = 1, nq
                        gq = zat(iat) * pref * exp(-q(iq)**2 / (4.0_dp*p))
                        s = 0.0_dp
                        do t = 0, lg(gi) + lg(gj)
                            do u = 0, lg(gi) + lg(gj) - t
                                do v = 0, lg(gi) + lg(gj) - t - u
                                    if (w(idx(t, u, v)) == 0.0_dp) cycle
                                    if (lorder == 0) then
                                        s = s + w(idx(t, u, v)) * rr(idx(t, u, v), iq)
                                    else if (q(iq) > 0.0_dp) then
                                        s = s + w(idx(t, u, v)) * (3.0_dp * rr(idx(t, u, v + 2), iq) &
                                                + q(iq)**2 * rr(idx(t, u, v), iq)) / (2.0_dp * q(iq)**2)
                                    end if
                                end do
                            end do
                        end do
                        res(iq) = res(iq) + gq * s
                    end do
                end do
            end do
        end do
        !$OMP END DO
        deallocate(ex, ey, ez, rr, w)
        !$OMP END PARALLEL
    end subroutine electron_nucleus

end module xs_nuclear

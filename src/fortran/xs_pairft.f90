! Part of THISTLE (kernel library pyxscat).
! SPDX-License-Identifier: Apache-2.0
! Copyright 2026 Andres Moreno Carrascosa
! Authors: A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander
!-----------------------------------------------------------------------
! PyXSCAT - Fourier transforms of basis-function products at fixed q
!
!   F(q, mu, nu) = int chi_mu(r) chi_nu(r) exp(+i q.r) dr
!
! for an arbitrary list of scattering vectors (no orientational average).
! Everything that depends on the molecular orientation (aligned or
! oriented molecules, detector images) follows from these by contraction
! with density matrices, because the plane wave factorises:
!
!   elastic      F(q)  = sum D_mu,nu F_mu,nu(q)
!   two-electron I2(q) = sum Gamma_{mu nu, la si} F_mu,nu(q) F_la,si(q)^*
!
! Each primitive pair is expanded in Hermite Gaussians at the pair centre P;
! the transform of a Hermite Gaussian is analytic:
!   int Lambda_tuv(r-P) exp(iq.r) = (pi/p)^{3/2} exp(-q^2/4p) exp(iq.P)
!                                   (i qx)^t (i qy)^u (i qz)^v
! so there is no Bessel expansion and no special treatment of coincident
! centres.
!-----------------------------------------------------------------------
module xs_pairft
    use iso_fortran_env, only: int64
    implicit none
    integer, parameter, private :: dp = selected_real_kind(15)
    integer, parameter, private :: ikind = selected_int_kind(8)
contains

    ! Hermite coefficients E(t; i, j) for one primitive pair and one direction (Gaussian
    ! prefactor excluded, E(0;0,0) = 1).
    pure subroutine hermite_pair(la, lb, pa, pb, p, e)
        integer(kind=ikind), intent(in) :: la, lb
        real(kind=dp), intent(in) :: pa, pb, p
        real(kind=dp), intent(out) :: e(0:la+lb+1, 0:la, 0:lb)
        ! work array padded in t (see xs_hermite: no guarded loads in the recurrence)
        real(kind=dp) :: w(-1:la+lb+2, 0:la, 0:lb)
        integer(kind=ikind) :: i, j, t
        real(kind=dp) :: a2
        a2 = 0.5_dp / p
        w = 0.0_dp
        w(0, 0, 0) = 1.0_dp
        do i = 0, la - 1
            do t = 0, i + 1
                w(t, i+1, 0) = a2*w(t-1, i, 0) + pa*w(t, i, 0) + real(t+1, dp)*w(t+1, i, 0)
            end do
        end do
        do j = 0, lb - 1
            do i = 0, la
                do t = 0, i + j + 1
                    w(t, i, j+1) = a2*w(t-1, i, j) + pb*w(t, i, j) + real(t+1, dp)*w(t+1, i, j)
                end do
            end do
        end do
        e = w(0:la+lb+1, :, :)
    end subroutine hermite_pair

    subroutine pair_transform(nprim, nao, nq, x, y, z, ga, l, m, n, group, contr, coeff, qv, &
            screen, f, ierr)
        integer(kind=ikind), intent(in) :: nprim, nao, nq
        real(kind=dp), intent(in) :: x(nprim), y(nprim), z(nprim), ga(nprim), coeff(nprim)
        integer(kind=ikind), intent(in) :: l(nprim), m(nprim), n(nprim), group(nprim), contr(nprim)
        real(kind=dp), intent(in) :: qv(3, nq), screen
        complex(kind=dp), intent(out) :: f(nq, nao, nao)
        integer(kind=ikind), intent(out) :: ierr

        integer(kind=ikind) :: ng, g, i, a, b, gi, gj, lmax, iq, t, u, v, ca, cb
        integer(kind=ikind), allocatable :: gstart(:), gcount(:), lg(:)
        real(kind=dp), allocatable :: cmax(:)
        real(kind=dp) :: pi, p, px, py, pz, kab, pref, cc
        real(kind=dp), allocatable :: ex(:, :, :), ey(:, :, :), ez(:, :, :)
        complex(kind=dp), allocatable :: g0(:), iqx(:, :), iqy(:, :), iqz(:, :), acc(:)
        complex(kind=dp) :: ci

        ierr = 0
        f = (0.0_dp, 0.0_dp)
        pi = acos(-1.0_dp)
        ci = (0.0_dp, 1.0_dp)
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

        ! powers (i q)^t for every q
        allocate(iqx(nq, 0:2*lmax), iqy(nq, 0:2*lmax), iqz(nq, 0:2*lmax))
        iqx(:, 0) = 1.0_dp; iqy(:, 0) = 1.0_dp; iqz(:, 0) = 1.0_dp
        do t = 1, 2*lmax
            iqx(:, t) = iqx(:, t-1) * ci * qv(1, :)
            iqy(:, t) = iqy(:, t-1) * ci * qv(2, :)
            iqz(:, t) = iqz(:, t-1) * ci * qv(3, :)
        end do

        ! Threads split the primitive-shell pairs; several pairs feed the same contracted
        ! element, so the accumulation into f is serialised (it is a small part of the work).
        !$OMP PARALLEL default(shared) private(gi, gj, a, b, p, px, py, pz, kab, pref, ex, ey, ez, &
        !$OMP& g0, acc, t, u, v, ca, cb, cc, iq)
        allocate(ex(0:2*lmax+1, 0:lmax, 0:lmax), ey(0:2*lmax+1, 0:lmax, 0:lmax), ez(0:2*lmax+1, 0:lmax, 0:lmax))
        allocate(g0(nq), acc(nq))
        !$OMP DO schedule(dynamic)
        do gi = 1, ng
            do gj = gi, ng
                a = gstart(gi); b = gstart(gj)
                p = ga(a) + ga(b)
                kab = exp(-ga(a)*ga(b)/p * ((x(a)-x(b))**2 + (y(a)-y(b))**2 + (z(a)-z(b))**2))
                pref = (pi/p)**1.5_dp * kab
                if (pref * cmax(gi) * cmax(gj) * (1.0_dp + real(lg(gi) + lg(gj), dp)) < screen) cycle
                px = (ga(a)*x(a) + ga(b)*x(b)) / p
                py = (ga(a)*y(a) + ga(b)*y(b)) / p
                pz = (ga(a)*z(a) + ga(b)*z(b)) / p
                call hermite_pair(lg(gi), lg(gj), px - x(a), px - x(b), p, ex(0:lg(gi)+lg(gj)+1, 0:lg(gi), 0:lg(gj)))
                call hermite_pair(lg(gi), lg(gj), py - y(a), py - y(b), p, ey(0:lg(gi)+lg(gj)+1, 0:lg(gi), 0:lg(gj)))
                call hermite_pair(lg(gi), lg(gj), pz - z(a), pz - z(b), p, ez(0:lg(gi)+lg(gj)+1, 0:lg(gi), 0:lg(gj)))
                do iq = 1, nq
                    g0(iq) = pref * exp(-(qv(1,iq)**2 + qv(2,iq)**2 + qv(3,iq)**2) / (4.0_dp*p)) &
                            * exp(ci * (qv(1,iq)*px + qv(2,iq)*py + qv(3,iq)*pz))
                end do
                do a = gstart(gi), gstart(gi) + gcount(gi) - 1
                    do b = gstart(gj), gstart(gj) + gcount(gj) - 1
                        cc = coeff(a) * coeff(b)
                        acc = (0.0_dp, 0.0_dp)
                        do t = 0, l(a) + l(b)
                            do u = 0, m(a) + m(b)
                                do v = 0, n(a) + n(b)
                                    acc = acc + (ex(t, l(a), l(b)) * ey(u, m(a), m(b)) * ez(v, n(a), n(b))) &
                                            * iqx(:, t) * iqy(:, u) * iqz(:, v)
                                end do
                            end do
                        end do
                        acc = acc * g0 * cc
                        ca = contr(a); cb = contr(b)
                        !$OMP CRITICAL (pairft_acc)
                        f(:, ca, cb) = f(:, ca, cb) + acc
                        if (gi /= gj) f(:, cb, ca) = f(:, cb, ca) + acc
                        !$OMP END CRITICAL (pairft_acc)
                    end do
                end do
            end do
        end do
        !$OMP END DO
        deallocate(ex, ey, ez, g0, acc)
        !$OMP END PARALLEL
    end subroutine pair_transform

end module xs_pairft

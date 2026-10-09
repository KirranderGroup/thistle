! Part of THISTLE (kernel library pyxscat).
! SPDX-License-Identifier: Apache-2.0
! Copyright 2026 Andres Moreno Carrascosa
! Authors: A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander
!-----------------------------------------------------------------------
! PyXSCAT - radial functions shared by the one- and two-centre j0 kernels
!
!   g_n(x) = j_n(x) / x^n                       (stable for every x >= 0)
!   a(T, i) = T! / (i! (T-2i)! 2^i)               (Cartesian derivatives of a function of h^2/2)
!
! For any radial function f(|h|) with f_n = (1/h d/dh)^n f,
!   d^T/dhx^T d^U/dhy^U d^V/dhz^V f = sum_ijk a(T,i) a(U,j) a(V,k) hx^(T-2i) hy^(U-2j) hz^(V-2k)
!                                     f_(T+U+V-i-j-k),
! and for f = j0(q h):  f_n = (-q^2)^n g_n(q h).
!-----------------------------------------------------------------------
module xs_radial
    implicit none
    integer, parameter, private :: dp = selected_real_kind(15)
    integer, parameter, private :: ikind = selected_int_kind(8)
    integer, parameter, private :: NJJ = 128
contains

    ! g(n) = j_n(x) / x^n, n = 0..nmax, accurate for every x >= 0.
    pure subroutine jn_over_xn(x, nmax, g)
        real(kind=dp), intent(in) :: x
        integer(kind=ikind), intent(in) :: nmax
        real(kind=dp), intent(out) :: g(0:nmax)
        integer(kind=ikind) :: n, k, nstart, nup
        real(kind=dp) :: term, s, x2, dfac, jm, j0, jp, j1, scale
        real(kind=dp) :: jj(0:NJJ), jd(0:NJJ)

        x2 = x * x
        if (x < 1.0_dp) then
            ! power series for the two highest orders, sum_k (-x^2/2)^k / (k! (2n+2k+1)!!),
            ! then the downward recurrence g(n-1) = (2n+1) g(n) - x^2 g(n+1), stable for x < n
            dfac = 1.0_dp
            do n = 0, nmax
                dfac = dfac * real(2*n + 1, dp)          ! (2n+1)!!
            end do
            do n = nmax, nmax + 1
                if (n == nmax + 1) dfac = dfac * real(2*n + 1, dp)
                term = 1.0_dp / dfac
                s = term
                do k = 1, 30
                    term = -term * x2 / (2.0_dp * k * (2*n + 2*k + 1))
                    s = s + term
                    if (abs(term) < 1.0e-17_dp * abs(s)) exit
                end do
                if (n == nmax) then
                    g(nmax) = s
                else
                    jp = s                                 ! g(nmax + 1)
                end if
            end do
            do n = nmax, 1, -1
                if (n == nmax) then
                    g(n - 1) = real(2*n + 1, dp) * g(n) - x2 * jp
                else
                    g(n - 1) = real(2*n + 1, dp) * g(n) - x2 * g(n + 1)
                end if
            end do
            return
        end if

        j0 = sin(x) / x
        if (x > real(nmax, dp) + 1.0_dp) then
            ! upward recurrence: stable while n < x
            g(0) = j0
            if (nmax == 0) return
            j1 = sin(x) / x2 - cos(x) / x
            jm = j0; j0 = j1
            g(1) = j1 / x
            scale = x
            do n = 1, nmax - 1
                jp = real(2*n + 1, dp) / x * j0 - jm
                jm = j0; j0 = jp
                scale = scale * x
                g(n + 1) = jp / scale
            end do
            return
        end if

        ! 1 <= x <= nmax + 1: combine both directions.  Upward (from j0, j1) for the orders
        ! n <= x, where it is stable; Miller's backward recurrence for the orders above x, where
        ! only the downward direction is stable, started high enough above nmax that its error
        ! (which falls by x^2/(4m^2 - 1) <= 1/4 per order) is below 1e-17 and matched to the upward value at
        ! the junction (at nup or nup - 1, whichever is farther from a zero of j_n).
        nup = min(int(x, kind=ikind), nmax)
        jj(0) = j0
        jj(1) = sin(x) / x2 - cos(x) / x
        do n = 1, nup - 1
            jj(n + 1) = real(2*n + 1, dp) / x * jj(n) - jj(n - 1)
        end do
        if (nmax > nup) then
            ! start where the neglected minimal solution is below 1e-17 relative at order nmax:
            ! its weight drops by x^2/(4m^2 - 1) per order above nmax
            nstart = nmax
            s = 1.0_dp
            do while (s > 1.0e-17_dp .and. nstart < NJJ - 2)
                nstart = nstart + 1
                s = s * x2 / real(4*nstart*nstart - 1, dp)
            end do
            nstart = nstart + 2
            jd(nstart + 1) = 0.0_dp
            jd(nstart) = 1.0e-30_dp
            do n = nstart, nup, -1
                jd(n - 1) = real(2*n + 1, dp) / x * jd(n) - jd(n + 1)
                if (abs(jd(n - 1)) > 1.0e200_dp) jd(n - 1:nstart) = jd(n - 1:nstart) * 1.0e-200_dp
            end do
            if (abs(jj(nup)) >= abs(jj(nup - 1))) then
                scale = jj(nup) / jd(nup)
            else
                scale = jj(nup - 1) / jd(nup - 1)
            end if
            do n = nup + 1, nmax
                jj(n) = jd(n) * scale
            end do
        end if
        s = 1.0_dp
        do n = 0, nmax
            g(n) = jj(n) / s
            s = s * x
        end do
    end subroutine jn_over_xn


    ! a(T, i) for T <= lmax, i <= T/2
    pure subroutine hermite_radial_coeffs(lmax, a)
        integer(kind=ikind), intent(in) :: lmax
        real(kind=dp), intent(out) :: a(0:lmax, 0:lmax/2)
        integer(kind=ikind) :: t, i
        a = 0.0_dp
        do t = 0, lmax
            a(t, 0) = 1.0_dp
            do i = 1, t / 2
                ! a(t,i) = a(t,i-1) * (t-2i+2)(t-2i+1) / (2 i)
                a(t, i) = a(t, i - 1) * real((t - 2*i + 2) * (t - 2*i + 1), dp) / real(2 * i, dp)
            end do
        end do
    end subroutine hermite_radial_coeffs

    ! The upward recurrence g(m+1) = [(2m+1) g(m) - g(m-1)] / x^2 is stable only while the
    ! order stays below x (above it the irregular solution y_m takes over, and rounding errors
    ! grow like (2m-1)!! (2m+1)!! / x^(2m+1)).  It is therefore used only when x > kk + 1, for
    ! all orders at once; below that, jn_over_xn propagates downward (series + downward
    ! recurrence for x < 1, Miller's normalised backward recurrence for 1 <= x <= kk + 1).
    pure real(kind=dp) function upward_threshold(kk)
        integer(kind=ikind), intent(in) :: kk
        upward_threshold = real(kk, dp) + 1.0_dp
    end function upward_threshold

    ! out(iq) = sum_{b=0}^{order} coef(b+1) q^(2b) j_b(q h)/(q h)^b
    !         = sum_b coef(b+1) j_b(q h) (q/h)^b,
    ! the radial sum of the quartet kernel, accurate for every q and h (including q h -> 0, where
    ! the closed forms of j_1, j_2, j_3 lose 2, 4, 6 powers of q h to cancellation).
    subroutine radial_series(nq, q, h, order, coef, out)
        integer(kind=ikind), intent(in) :: nq, order
        real(kind=dp), intent(in) :: q(nq), h, coef(order + 1)
        real(kind=dp), intent(out) :: out(nq)
        real(kind=dp) :: g(0:order), x, rx2, s, c, q2, qp, xup
        integer(kind=ikind) :: iq, b
        xup = upward_threshold(order)
        do iq = 1, nq
            x = q(iq) * h
            if (x > xup) then
                s = sin(x); c = cos(x)
                g(0) = s / x
                rx2 = 1.0_dp / (x * x)
                if (order >= 1) g(1) = (g(0) - c) * rx2
                do b = 1, order - 1
                    g(b + 1) = (real(2*b + 1, dp) * g(b) - g(b - 1)) * rx2
                end do
            else
                call jn_over_xn(x, order, g)
            end if
            q2 = q(iq) * q(iq)
            qp = 1.0_dp
            s = 0.0_dp
            do b = 0, order
                s = s + coef(b + 1) * qp * g(b)
                qp = qp * q2
            end do
            out(iq) = s
        end do
    end subroutine radial_series

end module xs_radial

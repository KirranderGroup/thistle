! Part of THISTLE (kernel library pyxscat).
! SPDX-License-Identifier: Apache-2.0
! Copyright 2026 Andres Moreno Carrascosa
! Authors: A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander
!-----------------------------------------------------------------------
! PyXSCAT - Hermite (McMurchie-Davidson) expansion coefficients
!
! For a pair of primitive shells A (centre xa, exponent a) and B (xb, b)
! the product of two 1D Cartesian Gaussians is expanded in Hermite
! Gaussians centred at P = (a xa + b xb)/(a + b):
!
!   (x-xa)^i (x-xb)^j exp(-a(x-xa)^2 - b(x-xb)^2)
!        = K_ab  sum_t E(t; i, j) Lambda_t(x - P)
!
! The Gaussian prefactor K_ab is NOT included here (it is carried by the
! pair prefactor of the kernel), so E(0; 0, 0) = 1.
!
! Two-term recurrences, valid for any angular momentum:
!   E(t; i+1, j) = E(t-1; i, j)/(2p) + (P-xa) E(t; i, j) + (t+1) E(t+1; i, j)
!   E(t; i, j+1) = E(t-1; i, j)/(2p) + (P-xb) E(t; i, j) + (t+1) E(t+1; i, j)
!
! This replaces the explicit tables of the original code (limited to g).
!-----------------------------------------------------------------------
module xs_hermite
    implicit none
    integer, parameter, private :: dp = selected_real_kind(15)
    integer, parameter, private :: ikind = selected_int_kind(8)
contains

    ! dd(t+1, i+1, j+1, gb, ga): coefficients for the shell pair (ga, gb),
    ! i <= lg(ga), j <= lg(gb).  The index order (.., gb, ga) is the one the
    ! kernel expects.
    subroutine hermite_table(ng, maxl, lg, xg, ag, dd)
        integer(kind=ikind), intent(in) :: ng, maxl
        integer(kind=ikind), intent(in) :: lg(ng)
        real(kind=dp), intent(in) :: xg(ng), ag(ng)
        real(kind=dp), intent(out) :: dd(2*maxl+1, maxl+1, maxl+1, ng, ng)
        integer(kind=ikind) :: ga, gb, i, j, t
        real(kind=dp) :: p, pc, pa, pb, a2
        ! The t index is padded on both sides (t = -1 and t = 2 maxl + 2 always hold zero), so
        ! the recurrence needs no guarded loads.  ifx 2026.1 at -O2/-O3 miscompiled the earlier
        ! form (with "if (t > 0) lo = e(t-1, ...)" inside an OpenMP loop over a private
        ! automatic array): the whole quartet kernel came out wrong.  The table is cheap
        ! (ng^2 short recurrences), so it is also no longer threaded.
        real(kind=dp) :: e(-1:2*maxl+2, 0:maxl, 0:maxl)

        do ga = 1, ng
            do gb = 1, ng
                p = ag(ga) + ag(gb)
                pc = (ag(ga)*xg(ga) + ag(gb)*xg(gb)) / p
                pa = pc - xg(ga)
                pb = pc - xg(gb)
                a2 = 0.5_dp / p
                e = 0.0_dp
                e(0, 0, 0) = 1.0_dp
                ! raise i with j = 0
                do i = 0, lg(ga) - 1
                    do t = 0, i + 1
                        e(t, i+1, 0) = a2*e(t-1, i, 0) + pa*e(t, i, 0) + real(t+1, dp)*e(t+1, i, 0)
                    end do
                end do
                ! raise j for every i
                do j = 0, lg(gb) - 1
                    do i = 0, lg(ga)
                        do t = 0, i + j + 1
                            e(t, i, j+1) = a2*e(t-1, i, j) + pb*e(t, i, j) + real(t+1, dp)*e(t+1, i, j)
                        end do
                    end do
                end do
                dd(:, :, :, gb, ga) = 0.0_dp
                do j = 0, lg(gb)
                    do i = 0, lg(ga)
                        do t = 0, i + j
                            dd(t+1, i+1, j+1, gb, ga) = e(t, i, j)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine hermite_table

end module xs_hermite

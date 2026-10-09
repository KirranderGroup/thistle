! Part of THISTLE (kernel library pyxscat).
! SPDX-License-Identifier: Apache-2.0
! Copyright 2026 Andres Moreno Carrascosa
! Authors: A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander
!-----------------------------------------------------------------------
! PyXSCAT - spin-summed 1- and 2-RDMs of a determinant expansion, bitwise
!
! The determinant-pair ("bitwise") algorithm of the earlier code (twordm_jeremy.f90):
! every pair (I, J) is screened with popcnt(ieor(.)) and only pairs differing by at
! most two spin-orbitals contribute, through the Slater-Condon rules.  Cost ~ N_det^2
! for the screening (a few machine instructions per pair, OpenMP over I) plus the
! coupled pairs; independent of the size of the full CI space.
!
! Any number of orbitals: each spin string is nw = ceil(norb / 64) 64-bit words
! (orbital k in word k/64, bit mod(k, 64)).  The RDMs are accumulated in thread-private
! copies when nthreads x 8 norb^4 bytes is small (< 2 GB), otherwise straight into the
! shared arrays with atomic updates (norb = 300: a 6.5 GB dm2, one copy).
!
! Conventions (as pyxscat.ci): determinant = alpha creators (ascending) then beta
! creators (ascending) on the vacuum;
!   dm1[p,q]     = sum_sigma <a+_{p sigma} a_{q sigma}>
!   dm2[p,q,r,s] = sum_{sigma tau} <a+_{p sigma} a+_{r tau} a_{s tau} a_{q sigma}>
! (C order on the Python side: here d1(q,p) and d2(s,r,q,p)).
!-----------------------------------------------------------------------
module xs_rdm
    use iso_fortran_env, only: int64, real64
    implicit none
    private
    public :: rdm12_bitwise, rdm12_bitwise_groups
    integer, parameter :: dp = real64, i8 = int64
    integer, parameter :: maxocc = 1024          ! occupied spin-orbitals of one determinant

    type :: accum
        logical :: atomic = .false.
        integer :: nint = -1                     ! >= 0: skip dm2 elements with all four orbitals >= nint
        real(dp), pointer :: a1(:, :) => null()
        real(dp), pointer :: a2(:, :, :, :) => null()
    end type accum

contains

    ! number of set bits of w(1:nw) below spin-orbital k (0-based)
    pure integer function count_below(w, nw, k)
        integer, intent(in) :: nw, k
        integer(i8), intent(in) :: w(nw)
        integer :: iw, kw
        iw = k / 64 + 1
        kw = mod(k, 64)
        count_below = popcnt(iand(w(iw), maskr(kw, i8)))
        if (iw > 1) count_below = count_below + sum(popcnt(w(1:iw-1)))
    end function count_below

    ! sign (-1)^(number of occupied spin-orbitals before (k, spin)) in alpha-first order
    pure integer function before(wa, wb, nw, k, spin)
        integer, intent(in) :: nw, k, spin
        integer(i8), intent(in) :: wa(nw), wb(nw)
        integer :: n
        if (spin == 0) then
            n = count_below(wa, nw, k)
        else
            n = sum(popcnt(wa)) + count_below(wb, nw, k)
        end if
        before = 1 - 2 * modulo(n, 2)
    end function before

    ! apply a_k (create = .false.) or a+_k (create = .true.) to (wa, wb); returns the sign
    integer function apply_op(wa, wb, nw, k, spin, create)
        integer, intent(in) :: nw, k, spin
        integer(i8), intent(inout) :: wa(nw), wb(nw)
        logical, intent(in) :: create
        integer :: iw, kw
        apply_op = before(wa, wb, nw, k, spin)
        iw = k / 64 + 1
        kw = mod(k, 64)
        if (spin == 0) then
            if (create) then
                wa(iw) = ibset(wa(iw), kw)
            else
                wa(iw) = ibclr(wa(iw), kw)
            end if
        else
            if (create) then
                wb(iw) = ibset(wb(iw), kw)
            else
                wb(iw) = ibclr(wb(iw), kw)
            end if
        end if
    end function apply_op

    ! the set spin-orbitals of the alpha and beta strings, as (orbital, spin) lists
    subroutine bits(wa, wb, nw, orb, spn, n)
        integer, intent(in) :: nw
        integer(i8), intent(in) :: wa(nw), wb(nw)
        integer, intent(out) :: orb(:), spn(:), n
        integer(i8) :: w
        integer :: iw
        n = 0
        do iw = 1, nw
            w = wa(iw)
            do while (w /= 0_i8)
                n = n + 1
                orb(n) = 64 * (iw - 1) + trailz(w)
                spn(n) = 0
                w = iand(w, w - 1_i8)
            end do
        end do
        do iw = 1, nw
            w = wb(iw)
            do while (w /= 0_i8)
                n = n + 1
                orb(n) = 64 * (iw - 1) + trailz(w)
                spn(n) = 1
                w = iand(w, w - 1_i8)
            end do
        end do
    end subroutine bits

    ! first set bit of w(1:nw) (-1 if none)
    pure integer function first_bit(w, nw)
        integer, intent(in) :: nw
        integer(i8), intent(in) :: w(nw)
        integer :: iw
        first_bit = -1
        do iw = 1, nw
            if (w(iw) /= 0_i8) then
                first_bit = 64 * (iw - 1) + trailz(w(iw))
                return
            end if
        end do
    end function first_bit

    subroutine rdm12_bitwise(ndet, norb, nw, wa, wb, c, d1, d2, nthreads, ierr)
        integer, intent(in) :: ndet, norb, nw, nthreads
        integer(i8), intent(in) :: wa(nw, ndet), wb(nw, ndet)
        real(dp), intent(in) :: c(ndet)
        real(dp), intent(out), target :: d1(norb, norb), d2(norb, norb, norb, norb)
        integer, intent(out) :: ierr
        real(dp), allocatable, target :: a1(:, :), a2(:, :, :, :)
        type(accum) :: acc
        integer :: i, j, nd, nthr, p, q, r, s, iw
        real(dp) :: bytes, tmp
        logical :: atomic

        ierr = 0
        if (norb < 1 .or. nw < (norb + 63) / 64) then
            ierr = 1
            return
        end if
        d1 = 0.0_dp
        d2 = 0.0_dp
        bytes = 8.0_dp * real(norb, dp)**4
        atomic = bytes * real(max(1, nthreads), dp) > 2.0e9_dp
        nthr = max(1, nthreads)
        !$omp parallel num_threads(nthr) default(shared) private(a1, a2, acc, i, j, nd, iw)
        acc%atomic = atomic
        acc%nint = -1                      ! (private copies: do not rely on default initialisation)
        if (atomic) then
            acc%a1 => d1
            acc%a2 => d2
        else
            allocate(a1(norb, norb), a2(norb, norb, norb, norb))
            a1 = 0.0_dp
            a2 = 0.0_dp
            acc%a1 => a1
            acc%a2 => a2
        end if
        !$omp do schedule(dynamic, 16)
        do i = 1, ndet
            call diagonal(wa(:, i), wb(:, i), nw, 0.5_dp * c(i) * c(i), acc)
            do j = i + 1, ndet
                nd = 0
                do iw = 1, nw
                    nd = nd + popcnt(ieor(wa(iw, i), wa(iw, j))) + popcnt(ieor(wb(iw, i), wb(iw, j)))
                end do
                if (nd > 4) cycle
                if (nd == 2) then
                    call single(wa(:, i), wb(:, i), wa(:, j), wb(:, j), nw, c(i) * c(j), acc)
                else if (nd == 4) then
                    call double(wa(:, i), wb(:, i), wa(:, j), wb(:, j), nw, c(i) * c(j), acc)
                end if
            end do
        end do
        !$omp end do
        if (.not. atomic) then
            !$omp critical
            d1 = d1 + a1
            d2 = d2 + a2
            !$omp end critical
            deallocate(a1, a2)
        end if
        !$omp end parallel
        ! the pairs (J, I) are the Hermitian conjugates of (I, J): D + D^T with
        ! dm1[p,q] <-> dm1[q,p] and dm2[p,q,r,s] <-> dm2[q,p,s,r] (diagonal entered with 1/2)
        d1 = d1 + transpose(d1)
        !$omp parallel do default(shared) private(p, q, r, s, tmp) schedule(dynamic)
        do p = 1, norb
            do q = 1, norb
                do r = 1, norb
                    do s = 1, norb
                        if (lower(p, q, r, s)) then
                            tmp = d2(s, r, q, p) + d2(r, s, p, q)
                            d2(s, r, q, p) = tmp
                            d2(r, s, p, q) = tmp
                        else if (p == q .and. r == s) then
                            d2(s, r, q, p) = 2.0_dp * d2(s, r, q, p)
                        end if
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    contains
        ! (p,q,r,s) strictly before its partner (q,p,s,r) in lexicographic order
        pure logical function lower(p, q, r, s)
            integer, intent(in) :: p, q, r, s
            if (p /= q) then
                lower = p < q
            else
                lower = r < s
            end if
        end function lower
    end subroutine rdm12_bitwise

    ! The same pair loop restricted to candidate pairs: determinants are sorted into groups
    ! (contiguous, gstart(g):gstart(g+1)-1), and for group g only the groups
    ! nb(nbstart(g):nbstart(g+1)-1), all >= g, are paired with it (g itself first when it
    ! is to be included).  Elements of dm2 with all four orbitals >= nint are not accumulated
    ! (they are formed elsewhere); nint < 0 keeps them.  The (I, J) / (J, I) symmetrisation
    ! is applied as in rdm12_bitwise, so pairs must be listed once (g <= g').
    subroutine rdm12_bitwise_groups(ndet, norb, nw, wa, wb, c, nint, ng, gstart, grp, nnb, nbstart, nb, &
                                    d1, d2, nthreads, ierr)
        integer, intent(in) :: ndet, norb, nw, nint, ng, nnb, nthreads
        integer(i8), intent(in) :: wa(nw, ndet), wb(nw, ndet)
        real(dp), intent(in) :: c(ndet)
        integer, intent(in) :: gstart(ng + 1), grp(ndet), nbstart(ng + 1), nb(nnb)
        real(dp), intent(out), target :: d1(norb, norb), d2(norb, norb, norb, norb)
        integer, intent(out) :: ierr
        real(dp), allocatable, target :: a1(:, :), a2(:, :, :, :)
        type(accum) :: acc
        integer :: i, j, nd, nthr, p, q, r, s, iw, g, k, g2, j0, j1
        real(dp) :: bytes, tmp
        logical :: atomic

        ierr = 0
        if (norb < 1 .or. nw < (norb + 63) / 64 .or. ng < 1) then
            ierr = 1
            return
        end if
        d1 = 0.0_dp
        d2 = 0.0_dp
        bytes = 8.0_dp * real(norb, dp)**4
        atomic = bytes * real(max(1, nthreads), dp) > 2.0e9_dp
        nthr = max(1, nthreads)
        !$omp parallel num_threads(nthr) default(shared) private(a1, a2, acc, i, j, nd, iw, g, k, g2, j0, j1)
        acc%atomic = atomic
        acc%nint = nint
        if (atomic) then
            acc%a1 => d1
            acc%a2 => d2
        else
            allocate(a1(norb, norb), a2(norb, norb, norb, norb))
            a1 = 0.0_dp
            a2 = 0.0_dp
            acc%a1 => a1
            acc%a2 => a2
        end if
        !$omp do schedule(dynamic, 16)
        do i = 1, ndet
            call diagonal(wa(:, i), wb(:, i), nw, 0.5_dp * c(i) * c(i), acc)
            g = grp(i)
            do k = nbstart(g), nbstart(g + 1) - 1
                g2 = nb(k)
                j0 = gstart(g2)
                if (g2 == g) j0 = i + 1
                j1 = gstart(g2 + 1) - 1
                do j = j0, j1
                    nd = 0
                    do iw = 1, nw
                        nd = nd + popcnt(ieor(wa(iw, i), wa(iw, j))) + popcnt(ieor(wb(iw, i), wb(iw, j)))
                    end do
                    if (nd > 4) cycle
                    if (nd == 2) then
                        call single(wa(:, i), wb(:, i), wa(:, j), wb(:, j), nw, c(i) * c(j), acc)
                    else if (nd == 4) then
                        call double(wa(:, i), wb(:, i), wa(:, j), wb(:, j), nw, c(i) * c(j), acc)
                    end if
                end do
            end do
        end do
        !$omp end do
        if (.not. atomic) then
            !$omp critical
            d1 = d1 + a1
            d2 = d2 + a2
            !$omp end critical
            deallocate(a1, a2)
        end if
        !$omp end parallel
        d1 = d1 + transpose(d1)
        !$omp parallel do default(shared) private(p, q, r, s, tmp) schedule(dynamic)
        do p = 1, norb
            do q = 1, norb
                do r = 1, norb
                    do s = 1, norb
                        if (lower2(p, q, r, s)) then
                            tmp = d2(s, r, q, p) + d2(r, s, p, q)
                            d2(s, r, q, p) = tmp
                            d2(r, s, p, q) = tmp
                        else if (p == q .and. r == s) then
                            d2(s, r, q, p) = 2.0_dp * d2(s, r, q, p)
                        end if
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    contains
        pure logical function lower2(p, q, r, s)
            integer, intent(in) :: p, q, r, s
            if (p /= q) then
                lower2 = p < q
            else
                lower2 = r < s
            end if
        end function lower2
    end subroutine rdm12_bitwise_groups

    ! <I| ... |I>: dm1[k,k] += w; dm2[u,u,v,v] += w over ordered pairs u /= v of occupied
    ! spin-orbitals, dm2[u,v,v,u] -= w when u and v have the same spin
    subroutine diagonal(wa, wb, nw, w, acc)
        integer, intent(in) :: nw
        integer(i8), intent(in) :: wa(nw), wb(nw)
        real(dp), intent(in) :: w
        type(accum), intent(inout) :: acc
        integer :: orb(maxocc), spn(maxocc), n, u, v
        call bits(wa, wb, nw, orb, spn, n)
        do u = 1, n
            call add1(acc, orb(u), orb(u), w)
            do v = 1, n
                if (u == v) cycle
                call add2(acc, orb(u), orb(u), orb(v), orb(v), w)
                if (spn(u) == spn(v)) call add2(acc, orb(u), orb(v), orb(v), orb(u), -w)
            end do
        end do
    end subroutine diagonal

    ! I = Phi a+_x a_h J (one spin-orbital differs); w = c_I c_J
    subroutine single(ia, ib, ja, jb, nw, w, acc)
        integer, intent(in) :: nw
        integer(i8), intent(in) :: ia(nw), ib(nw), ja(nw), jb(nw)
        real(dp), intent(in) :: w
        type(accum), intent(inout) :: acc
        integer(i8) :: ta(nw), tb(nw)
        integer :: orb(maxocc), spn(maxocc), n, x, h, sp, phi, k
        real(dp) :: v
        if (any(ia /= ja)) then
            sp = 0
            x = first_bit(iand(ia, not(ja)), nw)
            h = first_bit(iand(ja, not(ia)), nw)
        else
            sp = 1
            x = first_bit(iand(ib, not(jb)), nw)
            h = first_bit(iand(jb, not(ib)), nw)
        end if
        ta = ja
        tb = jb
        phi = apply_op(ta, tb, nw, h, sp, .false.)
        phi = phi * apply_op(ta, tb, nw, x, sp, .true.)
        v = real(phi, dp) * w
        call add1(acc, x, h, v)                                   ! dm1[x, h]
        call bits(iand(ia, ja), iand(ib, jb), nw, orb, spn, n)    ! common spin-orbitals
        do k = 1, n
            call add2(acc, x, h, orb(k), orb(k), v)
            call add2(acc, orb(k), orb(k), x, h, v)
            if (spn(k) == sp) then
                call add2(acc, x, orb(k), orb(k), h, -v)
                call add2(acc, orb(k), h, x, orb(k), -v)
            end if
        end do
    end subroutine single

    ! I = Phi a+_p1 a+_p2 a_h2 a_h1 J (two spin-orbitals differ); w = c_I c_J
    subroutine double(ia, ib, ja, jb, nw, w, acc)
        integer, intent(in) :: nw
        integer(i8), intent(in) :: ia(nw), ib(nw), ja(nw), jb(nw)
        real(dp), intent(in) :: w
        type(accum), intent(inout) :: acc
        integer :: po(4), ps(4), np, ho(4), hs(4), nh, phi
        integer(i8) :: ta(nw), tb(nw)
        real(dp) :: v
        call bits(iand(ia, not(ja)), iand(ib, not(jb)), nw, po, ps, np)     ! particles (in I only)
        call bits(iand(ja, not(ia)), iand(jb, not(ib)), nw, ho, hs, nh)     ! holes (in J only)
        ta = ja
        tb = jb
        phi = apply_op(ta, tb, nw, ho(1), hs(1), .false.)
        phi = phi * apply_op(ta, tb, nw, ho(2), hs(2), .false.)
        phi = phi * apply_op(ta, tb, nw, po(2), ps(2), .true.)
        phi = phi * apply_op(ta, tb, nw, po(1), ps(1), .true.)
        v = real(phi, dp) * w
        ! <I| a+_x a+_y a_z a_w |J>, {x,y} = {p1,p2}, {z,w} = {h1,h2}, spin(x)=spin(w), spin(y)=spin(z)
        if (ps(1) == hs(1) .and. ps(2) == hs(2)) then
            call add2(acc, po(1), ho(1), po(2), ho(2), v)
            call add2(acc, po(2), ho(2), po(1), ho(1), v)
        end if
        if (ps(1) == hs(2) .and. ps(2) == hs(1)) then
            call add2(acc, po(1), ho(2), po(2), ho(1), -v)
            call add2(acc, po(2), ho(1), po(1), ho(2), -v)
        end if
    end subroutine double

    ! dm1[p,q] (0-based, C order) += v
    subroutine add1(acc, p, q, v)
        type(accum), intent(inout) :: acc
        integer, intent(in) :: p, q
        real(dp), intent(in) :: v
        if (acc%atomic) then
            !$omp atomic update
            acc%a1(q + 1, p + 1) = acc%a1(q + 1, p + 1) + v
        else
            acc%a1(q + 1, p + 1) = acc%a1(q + 1, p + 1) + v
        end if
    end subroutine add1

    ! dm2[p,q,r,s] (0-based orbitals, C order) += v
    subroutine add2(acc, p, q, r, s, v)
        type(accum), intent(inout) :: acc
        integer, intent(in) :: p, q, r, s
        real(dp), intent(in) :: v
        if (acc%nint >= 0) then
            if (min(p, q, r, s) >= acc%nint) return
        end if
        if (acc%atomic) then
            !$omp atomic update
            acc%a2(s + 1, r + 1, q + 1, p + 1) = acc%a2(s + 1, r + 1, q + 1, p + 1) + v
        else
            acc%a2(s + 1, r + 1, q + 1, p + 1) = acc%a2(s + 1, r + 1, q + 1, p + 1) + v
        end if
    end subroutine add2

end module xs_rdm

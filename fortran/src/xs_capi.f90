! Part of THISTLE (kernel library pyxscat).
! SPDX-License-Identifier: Apache-2.0
! Copyright 2026 Andres Moreno Carrascosa
! Authors: A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander
!-----------------------------------------------------------------------
!> PyXSCAT - C-compatible entry points (called from Python through ctypes)
!>
!> One bind(C) wrapper per kernel; the ctypes signatures on the Python side are declared in
!> python/pyxscat/_lib.py.  Scalars are passed by value, arrays by reference in Fortran order,
!> optional arrays as NULL pointers; nthreads > 0 sets the OpenMP thread count for the call.
!-----------------------------------------------------------------------
module xs_capi
    use iso_c_binding
    use iso_fortran_env, only: compiler_version, compiler_options
    use xs_kernel_j0
    use xs_pairft
    use xs_nuclear
    use xs_elastic
    use xs_radial, only: radial_series
    use xs_rdm, only: rdm12_bitwise, rdm12_bitwise_groups
    use omp_lib
    implicit none
contains

    !> Rotationally averaged two-electron scattering.
    !>
    !>      mode 1 : dm2 is an (ncontr)^4 array with 8-fold symmetry; da, db unused
    !>      mode 2 : dm2 is ignored (pass NULL); da, db are symmetric (ncontr, ncontr)
    !>
    !> ierr = 0 on success (see xs_kernel_j0 for the codes).
    subroutine pyxscat_two_electron_j0(nprim, ncontr, nq, x, y, z, ga, l, m, n, group, contr, coeff, q, &
            mode, dm2, da, db, wc, wx, nkc, fac, fbc, nkx, fax, fbx, dgemm_ptr, &
            cutoffcentre, cutoffz, cutoffmd, screen_tol, lorder, nthreads, &
            res, stats, ierr) bind(C, name='pyxscat_two_electron_j0')
        integer(c_int), value, intent(in) :: nprim, ncontr, nq, mode, nthreads, lorder, nkc, nkx
        real(c_double), intent(in) :: x(nprim), y(nprim), z(nprim), ga(nprim), coeff(nprim), q(nq)
        integer(c_int), intent(in) :: l(nprim), m(nprim), n(nprim), group(nprim), contr(nprim)
        type(c_ptr), value, intent(in) :: dm2, fac, fbc, fax, fbx
        type(c_funptr), value, intent(in) :: dgemm_ptr
        real(c_double), intent(in) :: da(ncontr, ncontr), db(ncontr, ncontr)
        real(c_double), value, intent(in) :: wc, wx, cutoffcentre, cutoffz, cutoffmd, screen_tol
        real(c_double), intent(out) :: res(nq), stats(4)
        integer(c_int), intent(out) :: ierr
        real(c_double), pointer :: zbig(:, :, :, :), ac(:, :, :), bc(:, :, :), ax(:, :, :), bx(:, :, :)
        real(c_double), target :: none3(1, 1, 1)
        integer :: nt_old, mkc, mkx

        nullify(zbig)
        if (mode == 1 .and. c_associated(dm2)) call c_f_pointer(dm2, zbig, [ncontr, ncontr, ncontr, ncontr])
        ! factored densities (mode 3): arrays (k, i, j) with k fastest; absent ones -> a dummy
        mkc = 0; mkx = 0
        ac => none3; bc => none3; ax => none3; bx => none3
        if (mode == 3) then
            if (nkc > 0 .and. c_associated(fac) .and. c_associated(fbc)) then
                mkc = nkc
                call c_f_pointer(fac, ac, [nkc, ncontr, ncontr])
                call c_f_pointer(fbc, bc, [nkc, ncontr, ncontr])
            end if
            if (nkx > 0 .and. c_associated(fax) .and. c_associated(fbx)) then
                mkx = nkx
                call c_f_pointer(fax, ax, [nkx, ncontr, ncontr])
                call c_f_pointer(fbx, bx, [nkx, ncontr, ncontr])
            end if
        end if

        nt_old = omp_get_max_threads()
        if (nthreads > 0) call omp_set_num_threads(nthreads)
        if (mode == 3) then
            call two_electron_j0(nprim, ncontr, nq, x, y, z, ga, l, m, n, group, contr, coeff, q, &
                    mode, zbig, da, db, wc, wx, mkc, ac, bc, mkx, ax, bx, dgemm_ptr, &
                    cutoffcentre, cutoffz, cutoffmd, screen_tol, lorder, res, stats, ierr)
        else
            call two_electron_j0(nprim, ncontr, nq, x, y, z, ga, l, m, n, group, contr, coeff, q, &
                    mode, zbig, da, db, wc, wx, 0, none3, none3, 0, none3, none3, c_null_funptr, &
                    cutoffcentre, cutoffz, cutoffmd, screen_tol, lorder, res, stats, ierr)
        end if
        if (nthreads > 0) call omp_set_num_threads(nt_old)
    end subroutine pyxscat_two_electron_j0

    !> Hermite coefficients for one Cartesian direction (exposed for unit tests).
    subroutine pyxscat_hermite_table(ng, maxl, lg, xg, ag, dd) bind(C, name='pyxscat_hermite_table')
        integer(c_int), value, intent(in) :: ng, maxl
        integer(c_int), intent(in) :: lg(ng)
        real(c_double), intent(in) :: xg(ng), ag(ng)
        real(c_double), intent(out) :: dd(2*maxl+1, maxl+1, maxl+1, ng, ng)
        call hermite_table(ng, maxl, lg, xg, ag, dd)
    end subroutine pyxscat_hermite_table

    !> Fourier transforms of basis-function products, F(q, mu, nu), at fixed scattering vectors.
    subroutine pyxscat_pair_transform(nprim, nao, nq, x, y, z, ga, l, m, n, group, contr, coeff, qv, &
            screen, nthreads, f, ierr) bind(C, name='pyxscat_pair_transform')
        integer(c_int), value, intent(in) :: nprim, nao, nq, nthreads
        real(c_double), intent(in) :: x(nprim), y(nprim), z(nprim), ga(nprim), coeff(nprim), qv(3, nq)
        integer(c_int), intent(in) :: l(nprim), m(nprim), n(nprim), group(nprim), contr(nprim)
        real(c_double), value, intent(in) :: screen
        complex(c_double_complex), intent(out) :: f(nq, nao, nao)
        integer(c_int), intent(out) :: ierr
        integer :: nt_old
        nt_old = omp_get_max_threads()
        if (nthreads > 0) call omp_set_num_threads(nthreads)
        call pair_transform(nprim, nao, nq, x, y, z, ga, l, m, n, group, contr, coeff, qv, screen, f, ierr)
        if (nthreads > 0) call omp_set_num_threads(nt_old)
    end subroutine pyxscat_pair_transform

    !> Electron-nucleus term of electron scattering (rotational average, l = 0 or 2).
    subroutine pyxscat_electron_nucleus(nprim, nao, nq, x, y, z, ga, l, m, n, group, contr, coeff, q, &
            natom, zat, rat, dm, lorder, screen, nthreads, res, ierr) bind(C, name='pyxscat_electron_nucleus')
        integer(c_int), value, intent(in) :: nprim, nao, nq, natom, lorder, nthreads
        real(c_double), intent(in) :: x(nprim), y(nprim), z(nprim), ga(nprim), coeff(nprim), q(nq)
        integer(c_int), intent(in) :: l(nprim), m(nprim), n(nprim), group(nprim), contr(nprim)
        real(c_double), intent(in) :: zat(natom), rat(3, natom), dm(nao, nao)
        real(c_double), value, intent(in) :: screen
        real(c_double), intent(out) :: res(nq)
        integer(c_int), intent(out) :: ierr
        integer :: nt_old
        nt_old = omp_get_max_threads()
        if (nthreads > 0) call omp_set_num_threads(nthreads)
        call electron_nucleus(nprim, nao, nq, x, y, z, ga, l, m, n, group, contr, coeff, q, &
                natom, zat, rat, dm, lorder, screen, res, ierr)
        if (nthreads > 0) call omp_set_num_threads(nt_old)
    end subroutine pyxscat_electron_nucleus

    !> Elastic (and coherent mixed) scattering by density-weighted pairs (l = 0 or 2).
    subroutine pyxscat_elastic_pairs(nprim, nao, nq, x, y, z, ga, l, m, n, group, contr, coeff, q, &
            da, db, lorder, screen, nthreads, res, stats, ierr) bind(C, name='pyxscat_elastic_pairs')
        integer(c_int), value, intent(in) :: nprim, nao, nq, lorder, nthreads
        real(c_double), intent(in) :: x(nprim), y(nprim), z(nprim), ga(nprim), coeff(nprim), q(nq)
        integer(c_int), intent(in) :: l(nprim), m(nprim), n(nprim), group(nprim), contr(nprim)
        real(c_double), intent(in) :: da(nao, nao), db(nao, nao)
        real(c_double), value, intent(in) :: screen
        real(c_double), intent(out) :: res(nq), stats(2)
        integer(c_int), intent(out) :: ierr
        integer :: nt_old
        nt_old = omp_get_max_threads()
        if (nthreads > 0) call omp_set_num_threads(nthreads)
        call elastic_pairs(nprim, nao, nq, x, y, z, ga, l, m, n, group, contr, coeff, q, &
                da, db, lorder, screen, res, stats, ierr)
        if (nthreads > 0) call omp_set_num_threads(nt_old)
    end subroutine pyxscat_elastic_pairs

    !> sum_b coef(b+1) j_b(q h) (q/h)^b for a list of q (exposed for unit tests).
    subroutine pyxscat_radial_series(nq, q, h, order, coef, res) bind(C, name='pyxscat_radial_series')
        integer(c_int), value, intent(in) :: nq, order
        real(c_double), intent(in) :: q(nq), coef(order + 1)
        real(c_double), value, intent(in) :: h
        real(c_double), intent(out) :: res(nq)
        call radial_series(int(nq, kind=4), q, h, int(order, kind=4), coef, res)
    end subroutine pyxscat_radial_series

    !> Spin-summed 1- and 2-RDMs of a determinant expansion by the bitwise determinant-pair
    !> algorithm (one 64-bit word per spin).  d1, d2 are C-ordered [p,q] and [p,q,r,s].
    subroutine pyxscat_rdm12_bitwise(ndet, norb, nw, wa, wb, c, nthreads, d1, d2, ierr) &
            bind(C, name='pyxscat_rdm12_bitwise')
        integer(c_int), value, intent(in) :: ndet, norb, nw, nthreads
        integer(c_int64_t), intent(in) :: wa(nw, ndet), wb(nw, ndet)
        real(c_double), intent(in) :: c(ndet)
        real(c_double), intent(out) :: d1(norb, norb), d2(norb, norb, norb, norb)
        integer(c_int), intent(out) :: ierr
        integer :: nt, e
        nt = nthreads
        if (nt <= 0) nt = omp_get_max_threads()
        call rdm12_bitwise(int(ndet), int(norb), int(nw), wa, wb, c, d1, d2, nt, e)
        ierr = e
    end subroutine pyxscat_rdm12_bitwise

    !> rdm12_bitwise restricted to candidate pairs: determinants sorted into ng groups with the
    !> neighbour lists nb (see xs_rdm); dm2 elements with all four orbitals >= nint are skipped.
    subroutine pyxscat_rdm12_bitwise_groups(ndet, norb, nw, wa, wb, c, nint, ng, gstart, grp, nnb, nbstart, nb, &
                                            nthreads, d1, d2, ierr) bind(C, name='pyxscat_rdm12_bitwise_groups')
        integer(c_int), value, intent(in) :: ndet, norb, nw, nint, ng, nnb, nthreads
        integer(c_int64_t), intent(in) :: wa(nw, ndet), wb(nw, ndet)
        real(c_double), intent(in) :: c(ndet)
        integer(c_int), intent(in) :: gstart(ng + 1), grp(ndet), nbstart(ng + 1), nb(nnb)
        real(c_double), intent(out) :: d1(norb, norb), d2(norb, norb, norb, norb)
        integer(c_int), intent(out) :: ierr
        integer :: nt, e
        nt = nthreads
        if (nt <= 0) nt = omp_get_max_threads()
        call rdm12_bitwise_groups(int(ndet), int(norb), int(nw), wa, wb, c, int(nint), int(ng), gstart, grp, &
                                  int(nnb), nbstart, nb, d1, d2, nt, e)
        ierr = e
    end subroutine pyxscat_rdm12_bitwise_groups

    !> Compiler and options the library was built with (for bug reports): copied into buf.
    subroutine pyxscat_build_info(buf, n) bind(C, name='pyxscat_build_info')
        integer(c_int), value, intent(in) :: n
        character(kind=c_char), intent(out) :: buf(n)
        character(len=:), allocatable :: s
        integer :: i
        s = compiler_version() // ' | ' // compiler_options()
        buf = c_null_char
        do i = 1, min(len(s), n - 1)
            buf(i) = s(i:i)
        end do
    end subroutine pyxscat_build_info

    !> Number of OpenMP threads the library uses when nthreads <= 0 is passed (omp_get_max_threads).
    integer(c_int) function pyxscat_max_threads() bind(C, name='pyxscat_max_threads')
        pyxscat_max_threads = omp_get_max_threads()
    end function pyxscat_max_threads

end module xs_capi

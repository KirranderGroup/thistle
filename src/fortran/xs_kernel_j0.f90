! Part of THISTLE (kernel library pyxscat).
! SPDX-License-Identifier: Apache-2.0
! Copyright 2026 Andres Moreno Carrascosa
! Authors: A. Moreno Carrascosa, M. Simmermacher, P. Y. Wang, E. Moen Liane, J. Cooper, N. Zotev, T. Northey, A. Kirrander
!-----------------------------------------------------------------------
! PyXSCAT - rotationally averaged (j0) two-electron scattering kernel
!
!   I2(q) = sum_{ijkr} Z(ij,kr) < chi_i chi_j | j0(q r12) | chi_k chi_r >
!
! chi are contracted Cartesian Gaussians and Z is a pair density in that
! basis.  The routine knows nothing about how Z was obtained: it is
! supplied by a "density provider"
!
!   mode 1 : explicit four-index array Z (8-fold symmetric)
!   mode 2 : built on the fly from two symmetric one-particle matrices
!   mode 3 : a sum of products of one-particle matrices (factored densities, e.g. CCSD):
!              Z(a,b,c,d) = sum_k Ac_k(a,b) Bc_k(c,d)  [Coulomb pairing, Ac/Bc symmetric]
!                         + sum_k Ax_k(a,c) Bx_k(b,d)  [exchange pairing, any matrices],
!            symmetrised over the 8 index permutations inside the kernel
!              Z(ij,kr) = wc/2 [A_ij B_kr + B_ij A_kr]
!                       + wx   [A_ik A_jr + A_ir A_jk]
!            wc = 1, wx = 0     -> elastic (A = B = D) or coherent mixed
!            wc = 0, wx = -1/4  -> closed-shell exchange
!            wc = 1, wx = -1/4  -> closed-shell mean field (Hartree-Fock)
!
! Algorithm: canonical loop over contracted-shell quartets with primitive
! quartets inside, shell-pair screening, thread-local buffers.  It is the
! "one loop" kernel of Moreno Carrascosa and co-workers, taking its input
! in memory instead of from files.
!
! Basis description (all arrays over primitive Cartesian functions):
!   x, y, z, ga   centre and exponent
!   l, m, n       Cartesian powers
!   group         primitive shell (same centre, exponent and l+m+n);
!                 numbered 1..ng, members contiguous
!   contr         index of the contracted function the primitive belongs to
!   coeff         contraction coefficient times normalisation
!-----------------------------------------------------------------------
module xs_kernel_j0
    use Bessels_j0
    use p0_cases
    use xs_hermite
    use xs_radial, only: radial_series
    use iso_fortran_env, only: int64
    use iso_c_binding, only: c_char, c_int, c_double, c_funptr, c_associated, c_f_procpointer
    implicit none
    integer, parameter, private :: dp = selected_real_kind(15)
    integer, parameter, private :: ikind = selected_int_kind(8)
    integer(kind=ikind), parameter, private :: bdim = BDIM_MAX   ! from the Bessel module (h functions)

    abstract interface
        ! BLAS dgemm in its Fortran calling convention (everything by reference)
        subroutine dgemm_iface(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(C)
            import :: c_char, c_int, c_double
            character(kind=c_char), intent(in) :: transa, transb
            integer(c_int), intent(in) :: m, n, k, lda, ldb, ldc
            real(c_double), intent(in) :: alpha, beta
            real(c_double), intent(in) :: a(*), b(*)
            real(c_double), intent(inout) :: c(*)
        end subroutine dgemm_iface
    end interface

contains

    subroutine two_electron_j0(nprim, ncontr, nq, x, y, z, ga, l, m, n, group, contr, coeff, qin, &
            mode, zbig, da, db, wc, wx, nkc, ac, bc, nkx, ax, bx, dgemm_ptr, &
            cutoffcentre, cutoffz, cutoffmd, screen_tol, lorder, tsi, stats, ierr)
        integer(kind=ikind), intent(in) :: nprim, ncontr, nq, mode, lorder, nkc, nkx
        type(c_funptr), intent(in) :: dgemm_ptr            ! BLAS dgemm (Fortran ABI) or null: matmul
        real(kind=dp), intent(in) :: x(nprim), y(nprim), z(nprim), ga(nprim), coeff(nprim), qin(nq)
        integer(kind=ikind), intent(in) :: l(nprim), m(nprim), n(nprim), group(nprim), contr(nprim)
        real(kind=dp), pointer, intent(in) :: zbig(:, :, :, :)
        real(kind=dp), intent(in) :: da(ncontr, ncontr), db(ncontr, ncontr)
        real(kind=dp), intent(in) :: ac(nkc, ncontr, ncontr), bc(nkc, ncontr, ncontr)     ! k fastest
        real(kind=dp), intent(in) :: ax(nkx, ncontr, ncontr), bx(nkx, ncontr, ncontr)
        real(kind=dp), intent(in) :: wc, wx, cutoffcentre, cutoffz, cutoffmd, screen_tol
        real(kind=dp), intent(out) :: tsi(nq)
        real(kind=dp), intent(out) :: stats(4)
        integer(kind=ikind), intent(out) :: ierr

        real(kind=dp), allocatable :: q(:)
        real(kind=dp), allocatable :: p0matrix(:, :, :, :)
        real(kind=dp), allocatable, target :: ddx(:, :, :, :, :), ddy(:, :, :, :, :), ddz(:, :, :, :, :)
        real(kind=dp), allocatable :: px(:, :), py(:, :), pz(:, :), e12(:, :, :)
        real(kind=dp), allocatable :: xg(:), yg(:), zg(:), ag(:)
        integer(kind=ikind), allocatable :: lg(:), gstart(:), gcount(:)
        integer(kind=ikind), allocatable :: shell(:), shell_start(:), shell_end(:)
        integer(kind=ikind), allocatable :: tpi(:), tpj(:), order(:)
        real(kind=dp), allocatable :: wsh(:), psum(:), tsum(:), cost(:)
        integer(kind=ikind) :: npair, ip
        ! progress report on stderr (PYXSCAT_PROGRESS=1): cost-weighted share of the shell pairs started
        logical :: progress
        integer :: penv
        integer(int64) :: clk0, clk1, clkrate
        real(kind=dp) :: costall, costdone, mydone, mycost, elapsed
        real(kind=dp), allocatable :: e12max(:, :), zpair(:, :), zs(:, :), e12s(:, :), cmaxg(:)

        ! thread-local work space
        real(kind=dp), allocatable :: zc(:, :, :, :), c1(:), c2(:), c3(:), c4(:), hpre(:, :, :, :), f(:)
        integer(kind=ikind), allocatable :: posI(:), posJ(:), posK(:), posR(:)
        real(kind=dp), pointer :: dx1(:, :, :), dy1(:, :, :), dz1(:, :, :), dx2(:, :, :), dy2(:, :, :), dz2(:, :, :)

        integer(kind=ikind) :: ng, maxl, max_gc, nshell, n4
        integer(kind=ikind) :: i, j, k, r, ii, jj, a, b, c, d, g
        integer(kind=ikind) :: si, sj, sk, sr, sij, skr, srmax, gi0, gj0, gk0, gr0
        integer(kind=ikind) :: spi, spj, spk, spr, ij, kr
        real(kind=dp) :: pi, gap, hx, hy, hz, h, sym_factor, dmax, dmaxx
        real(kind=dp) :: bnd_ij, bnd, zblock, e12ij, maxe12s
        integer(kind=int64) :: nshq_tot, nshq_done, nquart_tot, nquart_done
        logical :: same
        procedure(dgemm_iface), pointer :: dgemm => null()

        ierr = 0
        tsi = 0.0_dp
        stats = 0.0_dp
        pi = acos(-1.0_dp)

        ! ---- checks ----------------------------------------------------------------------
        if (nprim < 1 .or. ncontr < 1 .or. nq < 1) then
            ierr = 1; return
        end if
        if (mode == 1 .and. .not. associated(zbig)) then
            ierr = 2; return
        end if
        if (mode /= 1 .and. mode /= 2 .and. mode /= 3) then
            ierr = 3; return
        end if
        if (mode == 3 .and. nkc + nkx < 1) then
            ierr = 3; return
        end if
        ng = maxval(group)
        maxl = maxval(l + m + n)
        if (lorder /= 0 .and. lorder /= 2) then
            ierr = 8; return
        end if
        if (4*maxl + lorder + 1 > bdim) then
            ierr = 4; return          ! beyond the derivative orders of the Bessel module
        end if
        if (group(1) /= 1) then
            ierr = 5; return
        end if
        do i = 2, nprim
            if (group(i) /= group(i-1) .and. group(i) /= group(i-1) + 1) then
                ierr = 5; return      ! groups must be contiguous and ascending
            end if
        end do
        if (minval(contr) < 1 .or. maxval(contr) > ncontr) then
            ierr = 6; return
        end if

        ! ---- primitive shells ------------------------------------------------------------
        allocate(gstart(ng), gcount(ng), lg(ng), xg(ng), yg(ng), zg(ng), ag(ng))
        gcount = 0
        do i = 1, nprim
            g = group(i)
            if (gcount(g) == 0) then
                gstart(g) = i
                lg(g) = l(i) + m(i) + n(i)
                xg(g) = x(i); yg(g) = y(i); zg(g) = z(i); ag(g) = ga(i)
            else
                if (l(i) + m(i) + n(i) /= lg(g) .or. ga(i) /= ag(g) .or. x(i) /= xg(g) &
                        .or. y(i) /= yg(g) .or. z(i) /= zg(g)) then
                    ierr = 7; return  ! members of a group must share centre, exponent and L
                end if
            end if
            gcount(g) = gcount(g) + 1
        end do
        max_gc = maxval(gcount)

        ! ---- q grid, coincident-centre table, Hermite tables, pair prefactors -------------
        allocate(q(nq))
        q = qin
        n4 = 4*maxl + lorder + 1
        allocate(p0matrix(nq, n4, n4, n4))
        call p0_limits(p0matrix, 4*maxl + lorder, q)

        allocate(ddx(2*maxl+1, maxl+1, maxl+1, ng, ng), ddy(2*maxl+1, maxl+1, maxl+1, ng, ng), &
                 ddz(2*maxl+1, maxl+1, maxl+1, ng, ng))
        call hermite_table(ng, maxl, lg, xg, ag, ddx)
        call hermite_table(ng, maxl, lg, yg, ag, ddy)
        call hermite_table(ng, maxl, lg, zg, ag, ddz)

        allocate(px(ng, ng), py(ng, ng), pz(ng, ng), e12(nq, ng, ng))
        do jj = 1, ng
            do ii = 1, ng
                gap = ag(ii) + ag(jj)
                px(ii, jj) = (ag(ii)*xg(ii) + ag(jj)*xg(jj)) / gap
                py(ii, jj) = (ag(ii)*yg(ii) + ag(jj)*yg(jj)) / gap
                pz(ii, jj) = (ag(ii)*zg(ii) + ag(jj)*zg(jj)) / gap
                e12(:, ii, jj) = (pi/gap)**1.5_dp * exp(-q*q*0.25_dp/gap) &
                        * exp(-ag(ii)*ag(jj)/gap * ((xg(ii) - xg(jj))**2 + (yg(ii) - yg(jj))**2 + (zg(ii) - zg(jj))**2))
            end do
        end do

        ! ---- contracted shells: consecutive primitive shells addressing the same functions -
        allocate(shell(ng), shell_start(ng), shell_end(ng))
        nshell = 0
        do i = 1, ng
            same = .false.
            if (i > 1) then
                if (gcount(i) == gcount(i-1)) then
                    same = all(contr(gstart(i):gstart(i)+gcount(i)-1) == contr(gstart(i-1):gstart(i-1)+gcount(i-1)-1))
                end if
            end if
            if (.not. same) then
                nshell = nshell + 1
                shell_start(nshell) = i
            end if
            shell(i) = nshell
            shell_end(nshell) = i
        end do

        ! ---- screening tables ---------------------------------------------------------------
        !   |contribution(i,j,k,r)| <~ e12max(i,j) e12max(k,r) c_i c_j c_k c_r max|Z block|
        allocate(e12max(ng, ng), cmaxg(ng), zs(nshell, nshell), e12s(nshell, nshell))
        e12max = maxval(e12, dim=1)
        do i = 1, ng
            cmaxg(i) = maxval(abs(coeff(gstart(i):gstart(i)+gcount(i)-1)))
        end do
        if (mode == 1) then
            allocate(zpair(ncontr, ncontr))
            zpair = 0.0_dp
            do jj = 1, ncontr
                do ii = 1, ncontr
                    zpair = max(zpair, abs(zbig(:, :, ii, jj)))
                end do
            end do
        else if (mode == 2) then
            dmax = max(maxval(abs(da)), maxval(abs(db)))
        else
            ! sum_k max|A_k| max|B_k| bounds every element of the pair-density block
            dmax = 0.0_dp
            do k = 1, nkc
                dmax = dmax + maxval(abs(ac(k, :, :)))*maxval(abs(bc(k, :, :)))
            end do
            dmaxx = 0.0_dp
            do k = 1, nkx
                dmaxx = dmaxx + maxval(abs(ax(k, :, :)))*maxval(abs(bx(k, :, :)))
            end do
            if (c_associated(dgemm_ptr)) call c_f_procpointer(dgemm_ptr, dgemm)
        end if
        do sj = 1, nshell
            do si = 1, nshell
                ii = shell_start(si); jj = shell_start(sj)
                if (mode == 1) then
                    zs(si, sj) = 0.0_dp
                    do r = gstart(jj), gstart(jj) + gcount(jj) - 1
                        do k = gstart(ii), gstart(ii) + gcount(ii) - 1
                            zs(si, sj) = max(zs(si, sj), zpair(contr(k), contr(r)))
                        end do
                    end do
                else if (mode == 2) then
                    zs(si, sj) = (abs(wc) + 2.0_dp*abs(wx)) * dmax * dmax      ! rigorous, loose
                else
                    zs(si, sj) = dmax
                end if
                e12s(si, sj) = 0.0_dp
                do j = shell_start(sj), shell_end(sj)
                    do i = shell_start(si), shell_end(si)
                        e12s(si, sj) = max(e12s(si, sj), e12max(i, j)*cmaxg(i)*cmaxg(j))
                    end do
                end do
            end do
        end do
        maxe12s = maxval(e12s)

        nshq_tot = 0; nshq_done = 0; nquart_tot = 0; nquart_done = 0

        ! ---- tasks: one per shell pair (si <= sj), largest first -------------------------------
        ! (one task per si left a few very large tasks: 15x at most on 32 threads for benzene)
        ! cost ~ w_i w_j sum_{(sk,sr) <= (si,sj)} w_k w_r, w = functions x primitives of a shell
        npair = merge(0_ikind, nshell*(nshell + 1)/2, mode == 3 .and. nkc == 0)
        allocate(tpi(npair), tpj(npair), order(npair), cost(npair), wsh(nshell), psum(0:nshell), tsum(0:nshell))
        do si = 1, nshell
            wsh(si) = real((shell_end(si) - shell_start(si) + 1)*gcount(shell_start(si)), dp)
        end do
        psum(0) = 0.0_dp; tsum(0) = 0.0_dp
        do si = 1, nshell
            psum(si) = psum(si - 1) + wsh(si)
        end do
        do si = 1, nshell
            tsum(si) = tsum(si - 1) + wsh(si)*(psum(nshell) - psum(si - 1))
        end do
        ip = 0
        if (npair > 0) then
            do si = 1, nshell
                do sj = si, nshell
                    ip = ip + 1
                    tpi(ip) = si; tpj(ip) = sj
                    cost(ip) = wsh(si)*wsh(sj)*(tsum(si - 1) + wsh(si)*(psum(sj) - psum(si - 1)))
                end do
            end do
            call sort_desc(npair, cost, order)
        end if
        call get_environment_variable("PYXSCAT_PROGRESS", length=penv)
        progress = penv > 0 .and. npair > 0
        costall = 0.0_dp; costdone = 0.0_dp
        if (progress) then
            costall = sum(cost)
            call system_clock(clk0, clkrate)
            write(0, '(a, i0, a)') "  quartet kernel: ", npair, " shell pairs"
            flush(0)
        end if

        ! =========================================================================================
        ! CANONICAL LOOP OVER CONTRACTED-SHELL QUARTETS, PRIMITIVE QUARTETS INSIDE
        ! =========================================================================================
        !$OMP PARALLEL default(shared) private(posI, posJ, posK, posR, spi, spj, spk, spr, zc, &
        !$OMP& f, h, hx, hy, hz, i, j, k, r, a, b, c, d, dx1, dx2, dy1, dy2, dz1, dz2, c1, c2, c3, c4, hpre, ij, kr, &
        !$OMP& sym_factor, si, sj, sk, sr, sij, skr, srmax, gi0, gj0, gk0, gr0, bnd_ij, bnd, zblock, e12ij, ip, &
        !$OMP& mydone, mycost, clk1, elapsed) &
        !$OMP& reduction(+:tsi, nshq_tot, nshq_done, nquart_tot, nquart_done)

        allocate(zc(max_gc, max_gc, max_gc, max_gc))
        allocate(c1(max_gc), c2(max_gc), c3(max_gc), c4(max_gc))
        allocate(posI(max_gc), posJ(max_gc), posK(max_gc), posR(max_gc))
        allocate(hpre(bdim, n4, n4, n4), f(nq))

        !$OMP DO schedule(dynamic)
        do ip = 1, npair
                si = tpi(order(ip)); sj = tpj(order(ip))
                if (progress) then
                    mycost = cost(order(ip))
                    !$OMP ATOMIC CAPTURE
                    costdone = costdone + mycost
                    mydone = costdone
                    !$OMP END ATOMIC
                    if (floor(20.0_dp*mydone/costall) > floor(20.0_dp*(mydone - mycost)/costall)) then
                        call system_clock(clk1)
                        elapsed = real(clk1 - clk0, dp)/real(clkrate, dp)
                        !$OMP CRITICAL (xs_progress)
                        write(0, '(a, i3, a, f9.0, a, f9.0, a)') "  quartet kernel: ", nint(100.0_dp*mydone/costall), &
                            "% started, ", elapsed, " s, about ", &
                            elapsed*(costall - mydone)/max(mydone, 1.0e-300_dp), " s left"
                        flush(0)
                        !$OMP END CRITICAL (xs_progress)
                    end if
                end if
                sij = si*nshell + sj
                bnd_ij = e12s(si, sj)
                if (bnd_ij*zs(si, sj)*maxe12s < screen_tol) cycle
                gi0 = shell_start(si); gj0 = shell_start(sj)
                spi = gcount(gi0); spj = gcount(gj0)
                posI(1:spi) = contr(gstart(gi0):gstart(gi0)+spi-1)
                posJ(1:spj) = contr(gstart(gj0):gstart(gj0)+spj-1)

                do sk = 1, si
                    srmax = min(nshell, sij - sk*nshell)                      ! exactly (sk,sr) <= (si,sj)
                    do sr = sk, srmax
                        skr = sk*nshell + sr
                        nshq_tot = nshq_tot + 1

                        bnd = bnd_ij*e12s(sk, sr)*min(zs(si, sj), zs(sk, sr))
                        if (bnd < screen_tol) cycle

                        gk0 = shell_start(sk); gr0 = shell_start(sr)
                        spk = gcount(gk0); spr = gcount(gr0)
                        posK(1:spk) = contr(gstart(gk0):gstart(gk0)+spk-1)
                        posR(1:spr) = contr(gstart(gr0):gstart(gr0)+spr-1)

                        ! --- the pair-density block, once per shell quartet ---------------------
                        if (mode == 1) then
                            ! explicit loops rather than a vector-subscripted section of a pointer
                            ! array (the simplest form for every compiler)
                            do d = 1, spr
                                do c = 1, spk
                                    do b = 1, spj
                                        do a = 1, spi
                                            zc(a, b, c, d) = zbig(posI(a), posJ(b), posK(c), posR(d))
                                        end do
                                    end do
                                end do
                            end do
                        else if (mode == 2) then
                            do d = 1, spr
                                do c = 1, spk
                                    do b = 1, spj
                                        do a = 1, spi
                                            zc(a, b, c, d) = &
                                                0.5_dp*wc*(da(posI(a), posJ(b))*db(posK(c), posR(d)) &
                                                         + db(posI(a), posJ(b))*da(posK(c), posR(d))) &
                                                + wx*(da(posI(a), posK(c))*da(posJ(b), posR(d)) &
                                                    + da(posI(a), posR(d))*da(posJ(b), posK(c)))
                                        end do
                                    end do
                                end do
                            end do
                        else
                            ! factored density: the Coulomb-paired terms here (their block is 8-fold
                            ! symmetric by itself), the exchange-paired ones in exchange_pass below
                            call coulomb_block(nkc, ncontr, ac, bc, posI, posJ, posK, posR, spi, spj, spk, spr, zc)
                        end if
                        zblock = maxval(abs(zc(1:spi, 1:spj, 1:spk, 1:spr)))
                        if (zblock*bnd_ij*e12s(sk, sr) < screen_tol) cycle
                        nshq_done = nshq_done + 1

                        ! --- primitive quartets of this shell quartet ---------------------------
                        do i = shell_start(si), shell_end(si)
                            c1(1:spi) = coeff(gstart(i):gstart(i)+spi-1)
                            do j = max(i, shell_start(sj)), shell_end(sj)                 ! i <= j
                                ij = i*ng + j
                                c2(1:spj) = coeff(gstart(j):gstart(j)+spj-1)
                                e12ij = e12max(i, j)*cmaxg(i)*cmaxg(j)*zblock
                                do k = shell_start(sk), shell_end(sk)
                                    c3(1:spk) = coeff(gstart(k):gstart(k)+spk-1)
                                    do r = max(k, shell_start(sr)), shell_end(sr)         ! k <= r
                                        kr = k*ng + r
                                        if (sij == skr .and. kr > ij) cycle               ! one representative per pair of pairs
                                        nquart_tot = nquart_tot + 1

                                        if (e12ij*e12max(k, r)*cmaxg(k)*cmaxg(r) < screen_tol) cycle
                                        c4(1:spr) = coeff(gstart(r):gstart(r)+spr-1)

                                        sym_factor = 1.0_dp
                                        if (i /= j) sym_factor = sym_factor*2.0_dp
                                        if (k /= r) sym_factor = sym_factor*2.0_dp
                                        if (ij /= kr) sym_factor = sym_factor*2.0_dp

                                        hx = px(k, r) - px(i, j)
                                        hy = py(k, r) - py(i, j)
                                        hz = pz(k, r) - pz(i, j)
                                        h = sqrt(hx*hx + hy*hy + hz*hz)

                                        dx1 => ddx(:, :, :, j, i)
                                        dy1 => ddy(:, :, :, j, i)
                                        dz1 => ddz(:, :, :, j, i)
                                        dx2 => ddx(:, :, :, r, k)
                                        dy2 => ddy(:, :, :, r, k)
                                        dz2 => ddz(:, :, :, r, k)

                                        if (h < cutoffcentre) then
                                            call quartet_coincident(nq, q, lorder, l, m, n, gstart, gcount, p0matrix, &
                                                    dx1, dy1, dz1, dx2, dy2, dz2, i, j, k, r, zc, c1, c2, c3, c4, &
                                                    cutoffz, cutoffmd, f)
                                        else
                                            call quartet_distant(q, lorder, l, m, n, gstart, gcount, hx, hy, hz, h, &
                                                    dx1, dy1, dz1, dx2, dy2, dz2, i, j, k, r, zc, c1, c2, c3, c4, &
                                                    cutoffz, cutoffmd, hpre, f)
                                        end if
                                        nquart_done = nquart_done + 1

                                        tsi = tsi + sym_factor*f*e12(:, i, j)*e12(:, k, r)
                                    end do
                                end do
                            end do
                        end do

                    end do
                end do
        end do
        !$OMP END DO

        deallocate(zc, c1, c2, c3, c4, posI, posJ, posK, posR, hpre, f)
        !$OMP END PARALLEL

        if (mode == 3 .and. nkx > 0) then
            call exchange_pass(nq, q, nkx, ncontr, ax, bx, dgemm, nshell, shell_start, shell_end, gstart, gcount, contr, &
                    coeff, max_gc, n4, lorder, l, m, n, p0matrix, ddx, ddy, ddz, px, py, pz, e12, e12s, e12max, cmaxg, ng, &
                    dmaxx, cutoffcentre, cutoffz, cutoffmd, screen_tol, tsi, nshq_tot, nshq_done, nquart_tot, nquart_done)
        end if

        stats(1) = real(nshq_tot, dp)
        stats(2) = real(nshq_done, dp)
        stats(3) = real(nquart_tot, dp)
        stats(4) = real(nquart_done, dp)

    end subroutine two_electron_j0

    ! ---------------------------------------------------------------------------------
    ! Exchange-paired factored terms, Z(a,b,c,d) = sum_k Ax_k(I_a, K_c) Bx_k(J_b, L_d), symmetrised
    ! over the 8 index permutations the canonical quartet loop assumes:
    !   1/8 [A(I,K)B(J,L) + A(J,K)B(I,L) + A(I,L)B(J,K) + A(J,L)B(I,K) + the same with A, B transposed].
    ! The shell quartets are the canonical ones (si <= sj, sk <= sr, (si,sj) >= (sk,sr)), visited
    ! with sk outside and a batch of si next, so that every factor block carries either batch
    ! indices (gathered once per batch or per sr) or inner indices (gathered once per sj and
    ! shared by the batch); each image is one matrix product over k (BLAS dgemm when given).
    ! ---------------------------------------------------------------------------------
    subroutine exchange_pass(nq, q, nkx, ncontr, ax, bx, dgemm, nshell, shell_start, shell_end, gstart, gcount, contr, &
            coeff, max_gc, n4, lorder, l, m, n, p0matrix, ddx, ddy, ddz, px, py, pz, e12, e12s, e12max, cmaxg, ng, &
            dmaxx, cutoffcentre, cutoffz, cutoffmd, screen_tol, tsi, nshq_tot, nshq_done, nquart_tot, nquart_done)
        integer(kind=ikind), intent(in) :: nq, nkx, ncontr, nshell, max_gc, n4, lorder, ng
        real(kind=dp), intent(in), allocatable :: q(:)
        real(kind=dp), intent(in) :: ax(nkx, ncontr, ncontr), bx(nkx, ncontr, ncontr)
        procedure(dgemm_iface), pointer, intent(in) :: dgemm
        integer(kind=ikind), intent(in) :: shell_start(:), shell_end(:), gstart(:), gcount(:), contr(:), l(:), m(:), n(:)
        real(kind=dp), intent(in) :: coeff(:), p0matrix(:, :, :, :)
        real(kind=dp), intent(in), target :: ddx(:, :, :, :, :), ddy(:, :, :, :, :), ddz(:, :, :, :, :)
        real(kind=dp), intent(in) :: px(:, :), py(:, :), pz(:, :), e12(:, :, :), e12s(:, :), e12max(:, :), cmaxg(:)
        real(kind=dp), intent(in) :: dmaxx, cutoffcentre, cutoffz, cutoffmd, screen_tol
        real(kind=dp), intent(inout) :: tsi(:)
        integer(kind=int64), intent(inout) :: nshq_tot, nshq_done, nquart_tot, nquart_done

        integer(kind=ikind), parameter :: batch = 12
        integer(kind=ikind) :: ntask, it, sk, s0, s1, si, sj, sr, nci, ci0, spi, spj, spk, spr, mjl, mjk, mik, mil
        integer(kind=ikind) :: i, j, k, r, ij, kr, a, b, c, d, sij, skr, nb, t0, t1, ncj, cj0
        integer(kind=ikind), allocatable :: tsk(:), ts0(:), ts1(:), cstart(:), posI(:), posJ(:), posK(:), posR(:), posB(:)
        integer(kind=ikind), allocatable :: jstart(:), posC(:)
        real(kind=dp), allocatable :: aik(:, :), bik(:, :), atik(:, :), btik(:, :)         ! per batch (I, K)
        real(kind=dp), allocatable :: ail(:, :), bil(:, :), atil(:, :), btil(:, :)         ! per sr    (I, L)
        real(kind=dp), allocatable :: ajk(:, :), bjk(:, :), ajl(:, :), bjl(:, :), atjk(:, :), btjk(:, :), atjl(:, :), btjl(:, :)
        real(kind=dp), allocatable :: y1(:, :), y2(:, :), y3(:, :), y4(:, :), y5(:, :), y6(:, :), y7(:, :), y8(:, :)
        real(kind=dp), allocatable :: zc(:, :, :, :), c1(:), c2(:), c3(:), c4(:), hpre(:, :, :, :), f(:)
        real(kind=dp) :: bnd, zblock, sym_factor, hx, hy, hz, h, e12ij
        real(kind=dp), pointer :: dx1(:, :, :), dy1(:, :, :), dz1(:, :, :), dx2(:, :, :), dy2(:, :, :), dz2(:, :, :)

        ! tasks: (sk, batch of si >= sk)
        ntask = 0
        do sk = 1, nshell
            ntask = ntask + (nshell - sk + batch)/batch
        end do
        allocate(tsk(ntask), ts0(ntask), ts1(ntask))
        it = 0
        do sk = 1, nshell
            do s0 = sk, nshell, batch
                it = it + 1
                tsk(it) = sk; ts0(it) = s0; ts1(it) = min(nshell, s0 + batch - 1)
            end do
        end do

        !$OMP PARALLEL default(shared) private(it, sk, s0, s1, si, sj, sr, nci, ci0, spi, spj, spk, spr, mjl, mjk, mik, mil, &
        !$OMP& i, j, k, r, ij, kr, a, b, c, d, sij, skr, nb, t0, t1, ncj, cj0, cstart, posI, posJ, posK, posR, posB, &
        !$OMP& jstart, posC, aik, bik, atik, btik, ail, bil, atil, btil, ajk, bjk, ajl, bjl, atjk, btjk, atjl, btjl, &
        !$OMP& y1, y2, y3, y4, y5, y6, y7, y8, zc, c1, c2, c3, c4, hpre, f, bnd, zblock, sym_factor, hx, hy, hz, h, e12ij, &
        !$OMP& dx1, dy1, dz1, dx2, dy2, dz2) reduction(+:tsi, nshq_tot, nshq_done, nquart_tot, nquart_done)
        nb = batch*max_gc
        allocate(cstart(batch + 1), jstart(batch + 1), posI(max_gc), posJ(max_gc), posK(max_gc), posR(max_gc), posB(nb), posC(nb))
        allocate(zc(max_gc, max_gc, max_gc, max_gc), c1(max_gc), c2(max_gc), c3(max_gc), c4(max_gc))
        allocate(hpre(bdim, n4, n4, n4), f(nq))
        allocate(aik(nkx, nb*max_gc), bik(nkx, nb*max_gc), atik(nkx, nb*max_gc), btik(nkx, nb*max_gc))
        allocate(ail(nkx, nb*max_gc), bil(nkx, nb*max_gc), atil(nkx, nb*max_gc), btil(nkx, nb*max_gc))
        allocate(ajk(nkx, nb*max_gc), bjk(nkx, nb*max_gc), ajl(nkx, nb*max_gc), bjl(nkx, nb*max_gc), &
                 atjk(nkx, nb*max_gc), btjk(nkx, nb*max_gc), atjl(nkx, nb*max_gc), btjl(nkx, nb*max_gc))
        allocate(y1(nb*max_gc, nb*max_gc), y2(nb*max_gc, nb*max_gc), y3(nb*max_gc, nb*max_gc), y4(nb*max_gc, nb*max_gc), &
                 y5(nb*max_gc, nb*max_gc), y6(nb*max_gc, nb*max_gc), y7(nb*max_gc, nb*max_gc), y8(nb*max_gc, nb*max_gc))

        !$OMP DO schedule(dynamic)
        do it = 1, ntask
            sk = tsk(it); s0 = ts0(it); s1 = ts1(it)
            spk = gcount(shell_start(sk))
            posK(1:spk) = contr(gstart(shell_start(sk)):gstart(shell_start(sk))+spk-1)
            nci = 0                                            ! functions of the si batch, one list
            do si = s0, s1
                cstart(si - s0 + 1) = nci + 1
                spi = gcount(shell_start(si))
                posB(nci+1:nci+spi) = contr(gstart(shell_start(si)):gstart(shell_start(si))+spi-1)
                nci = nci + spi
            end do
            cstart(s1 - s0 + 2) = nci + 1
            if (maxval(e12s(s0:s1, :))*maxval(e12s(sk, :))*dmaxx < screen_tol) cycle
            mik = nci*spk
            do c = 1, spk
                do a = 1, nci
                    aik(:, a + (c-1)*nci) = ax(:, posB(a), posK(c))
                    bik(:, a + (c-1)*nci) = bx(:, posB(a), posK(c))
                    atik(:, a + (c-1)*nci) = ax(:, posK(c), posB(a))
                    btik(:, a + (c-1)*nci) = bx(:, posK(c), posB(a))
                end do
            end do
            do sr = sk, nshell
                spr = gcount(shell_start(sr))
                skr = sk*nshell + sr
                if (maxval(e12s(s0:s1, :))*e12s(sk, sr)*dmaxx < screen_tol) cycle
                posR(1:spr) = contr(gstart(shell_start(sr)):gstart(shell_start(sr))+spr-1)
                mil = nci*spr
                do d = 1, spr
                    do a = 1, nci
                        ail(:, a + (d-1)*nci) = ax(:, posB(a), posR(d))
                        bil(:, a + (d-1)*nci) = bx(:, posB(a), posR(d))
                        atil(:, a + (d-1)*nci) = ax(:, posR(d), posB(a))
                        btil(:, a + (d-1)*nci) = bx(:, posR(d), posB(a))
                    end do
                end do
                do t0 = s0, nshell, batch                      ! batches of sj (sj >= si for some si)
                    t1 = min(nshell, t0 + batch - 1)
                    if (maxval(e12s(s0:s1, t0:t1))*e12s(sk, sr)*dmaxx < screen_tol) cycle
                    ncj = 0
                    do sj = t0, t1
                        jstart(sj - t0 + 1) = ncj + 1
                        spj = gcount(shell_start(sj))
                        posC(ncj+1:ncj+spj) = contr(gstart(shell_start(sj)):gstart(shell_start(sj))+spj-1)
                        ncj = ncj + spj
                    end do
                    jstart(t1 - t0 + 2) = ncj + 1
                    mjl = ncj*spr; mjk = ncj*spk
                    do d = 1, spr
                        do b = 1, ncj
                            ajl(:, b + (d-1)*ncj) = ax(:, posC(b), posR(d))
                            bjl(:, b + (d-1)*ncj) = bx(:, posC(b), posR(d))
                            atjl(:, b + (d-1)*ncj) = ax(:, posR(d), posC(b))
                            btjl(:, b + (d-1)*ncj) = bx(:, posR(d), posC(b))
                        end do
                    end do
                    do c = 1, spk
                        do b = 1, ncj
                            ajk(:, b + (c-1)*ncj) = ax(:, posC(b), posK(c))
                            bjk(:, b + (c-1)*ncj) = bx(:, posC(b), posK(c))
                            atjk(:, b + (c-1)*ncj) = ax(:, posK(c), posC(b))
                            btjk(:, b + (c-1)*ncj) = bx(:, posK(c), posC(b))
                        end do
                    end do
                    ! the eight images for the whole (si batch) x (sj batch), one matrix product each
                    call gemm_tn(dgemm, nkx, mik, mjl, aik, bjl, y1)        ! y1[(I,K),(J,L)] = A(I,K) B(J,L)
                    call gemm_tn(dgemm, nkx, mjk, mil, ajk, bil, y2)        ! y2[(J,K),(I,L)] = A(J,K) B(I,L)
                    call gemm_tn(dgemm, nkx, mil, mjk, ail, bjk, y3)        ! y3[(I,L),(J,K)] = A(I,L) B(J,K)
                    call gemm_tn(dgemm, nkx, mjl, mik, ajl, bik, y4)        ! y4[(J,L),(I,K)] = A(J,L) B(I,K)
                    call gemm_tn(dgemm, nkx, mik, mjl, atik, btjl, y5)      ! transposed factors
                    call gemm_tn(dgemm, nkx, mjk, mil, atjk, btil, y6)
                    call gemm_tn(dgemm, nkx, mil, mjk, atil, btjk, y7)
                    call gemm_tn(dgemm, nkx, mjl, mik, atjl, btik, y8)

                    do sj = t0, t1
                        spj = gcount(shell_start(sj))
                        cj0 = jstart(sj - t0 + 1) - 1
                        posJ(1:spj) = posC(cj0+1:cj0+spj)
                        do si = max(s0, sk), min(s1, sj)               ! canonical: si <= sj, (si,sj) >= (sk,sr)
                            sij = si*nshell + sj
                            if (sij < skr) cycle
                            nshq_tot = nshq_tot + 1
                            bnd = e12s(si, sj)*e12s(sk, sr)
                            if (bnd*dmaxx < screen_tol) cycle
                            spi = gcount(shell_start(si))
                            ci0 = cstart(si - s0 + 1) - 1
                            posI(1:spi) = posB(ci0+1:ci0+spi)
                            do d = 1, spr
                                do c = 1, spk
                                    do b = 1, spj
                                        do a = 1, spi
                                            zc(a, b, c, d) = 0.125_dp*( &
                                                  y1(ci0 + a + (c-1)*nci, cj0 + b + (d-1)*ncj) &
                                                + y5(ci0 + a + (c-1)*nci, cj0 + b + (d-1)*ncj) &
                                                + y2(cj0 + b + (c-1)*ncj, ci0 + a + (d-1)*nci) &
                                                + y6(cj0 + b + (c-1)*ncj, ci0 + a + (d-1)*nci) &
                                                + y3(ci0 + a + (d-1)*nci, cj0 + b + (c-1)*ncj) &
                                                + y7(ci0 + a + (d-1)*nci, cj0 + b + (c-1)*ncj) &
                                                + y4(cj0 + b + (d-1)*ncj, ci0 + a + (c-1)*nci) &
                                                + y8(cj0 + b + (d-1)*ncj, ci0 + a + (c-1)*nci))
                                        end do
                                    end do
                                end do
                            end do
                            zblock = maxval(abs(zc(1:spi, 1:spj, 1:spk, 1:spr)))
                            if (zblock*bnd < screen_tol) cycle
                            nshq_done = nshq_done + 1

                            ! --- canonical primitive quartets, as in the main loop ---
                            do i = shell_start(si), shell_end(si)
                                c1(1:spi) = coeff(gstart(i):gstart(i)+spi-1)
                                do j = max(i, shell_start(sj)), shell_end(sj)
                                    ij = i*ng + j
                                    c2(1:spj) = coeff(gstart(j):gstart(j)+spj-1)
                                    e12ij = e12max(i, j)*cmaxg(i)*cmaxg(j)*zblock
                                    do k = shell_start(sk), shell_end(sk)
                                        c3(1:spk) = coeff(gstart(k):gstart(k)+spk-1)
                                        do r = max(k, shell_start(sr)), shell_end(sr)
                                            kr = k*ng + r
                                            if (sij == skr .and. kr > ij) cycle
                                            nquart_tot = nquart_tot + 1
                                            if (e12ij*e12max(k, r)*cmaxg(k)*cmaxg(r) < screen_tol) cycle
                                            c4(1:spr) = coeff(gstart(r):gstart(r)+spr-1)
                                            sym_factor = 1.0_dp
                                            if (i /= j) sym_factor = sym_factor*2.0_dp
                                            if (k /= r) sym_factor = sym_factor*2.0_dp
                                            if (ij /= kr) sym_factor = sym_factor*2.0_dp
                                            hx = px(k, r) - px(i, j)
                                            hy = py(k, r) - py(i, j)
                                            hz = pz(k, r) - pz(i, j)
                                            h = sqrt(hx*hx + hy*hy + hz*hz)
                                            dx1 => ddx(:, :, :, j, i); dy1 => ddy(:, :, :, j, i); dz1 => ddz(:, :, :, j, i)
                                            dx2 => ddx(:, :, :, r, k); dy2 => ddy(:, :, :, r, k); dz2 => ddz(:, :, :, r, k)
                                            if (h < cutoffcentre) then
                                                call quartet_coincident(nq, q, lorder, l, m, n, gstart, gcount, p0matrix, &
                                                        dx1, dy1, dz1, dx2, dy2, dz2, i, j, k, r, zc, c1, c2, c3, c4, &
                                                        cutoffz, cutoffmd, f)
                                            else
                                                call quartet_distant(q, lorder, l, m, n, gstart, gcount, hx, hy, hz, h, &
                                                        dx1, dy1, dz1, dx2, dy2, dz2, i, j, k, r, zc, c1, c2, c3, c4, &
                                                        cutoffz, cutoffmd, hpre, f)
                                            end if
                                            nquart_done = nquart_done + 1
                                            tsi = tsi + sym_factor*f*e12(:, i, j)*e12(:, k, r)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$OMP END DO
        deallocate(cstart, jstart, posI, posJ, posK, posR, posB, posC, zc, c1, c2, c3, c4, hpre, f, aik, bik, atik, btik, &
                   ail, bil, atil, btil, ajk, bjk, ajl, bjl, atjk, btjk, atjl, btjl, y1, y2, y3, y4, y5, y6, y7, y8)
        !$OMP END PARALLEL
        deallocate(tsk, ts0, ts1)
    end subroutine exchange_pass

    ! z(1:ncol, 1:n) = a(1:nk, 1:ncol)^T b(1:nk, 1:n)
    subroutine sort_desc(n, key, idx)
        ! idx = permutation that orders key descending (iterative merge sort, stable)
        integer(kind=ikind), intent(in) :: n
        real(kind=dp), intent(in) :: key(n)
        integer(kind=ikind), intent(out) :: idx(n)
        integer(kind=ikind), allocatable :: tmp(:)
        integer(kind=ikind) :: width, lo, mid, hi, i, j, k
        allocate(tmp(n))
        do i = 1, n
            idx(i) = i
        end do
        width = 1
        do while (width < n)
            lo = 1
            do while (lo <= n)
                mid = min(lo + width - 1, n); hi = min(lo + 2*width - 1, n)
                i = lo; j = mid + 1; k = lo
                do while (i <= mid .and. j <= hi)
                    if (key(idx(i)) >= key(idx(j))) then
                        tmp(k) = idx(i); i = i + 1
                    else
                        tmp(k) = idx(j); j = j + 1
                    end if
                    k = k + 1
                end do
                do while (i <= mid)
                    tmp(k) = idx(i); i = i + 1; k = k + 1
                end do
                do while (j <= hi)
                    tmp(k) = idx(j); j = j + 1; k = k + 1
                end do
                lo = lo + 2*width
            end do
            idx = tmp
            width = 2*width
        end do
        deallocate(tmp)
    end subroutine sort_desc

    subroutine gemm_tn(dgemm, nk, ncol, n, a, b, z)
        procedure(dgemm_iface), pointer, intent(in) :: dgemm
        integer(kind=ikind), intent(in) :: nk, ncol, n
        real(kind=dp), intent(in) :: a(:, :), b(:, :)
        real(kind=dp), intent(inout) :: z(:, :)
        integer(c_int) :: m_, n_, k_, lda, ldb, ldc
        if (associated(dgemm)) then
            m_ = int(ncol, c_int); n_ = int(n, c_int); k_ = int(nk, c_int)
            lda = int(size(a, 1), c_int); ldb = int(size(b, 1), c_int); ldc = int(size(z, 1), c_int)
            call dgemm('T', 'N', m_, n_, k_, 1.0_dp, a, lda, b, ldb, 0.0_dp, z, ldc)
        else
            z(1:ncol, 1:n) = matmul(transpose(a(1:nk, 1:ncol)), b(1:nk, 1:n))
        end if
    end subroutine gemm_tn



    ! Pair-density block of a shell quartet from the Coulomb-paired factored terms,
    !   sum_k 1/2 [Ac_k(I,J) Bc_k(K,L) + Bc_k(I,J) Ac_k(K,L)]     (Ac, Bc symmetric),
    ! which has the full 8-fold index symmetry of the canonical loop.
    subroutine coulomb_block(nkc, ncontr, ac, bc, posI, posJ, posK, posR, spi, spj, spk, spr, zc)
        integer(kind=ikind), intent(in) :: nkc, ncontr, spi, spj, spk, spr
        real(kind=dp), intent(in) :: ac(nkc, ncontr, ncontr), bc(nkc, ncontr, ncontr)
        integer(kind=ikind), intent(in) :: posI(:), posJ(:), posK(:), posR(:)
        real(kind=dp), intent(inout) :: zc(:, :, :, :)
        real(kind=dp) :: gij_a(nkc, spi*spj), gij_b(nkc, spi*spj), gkl_a(nkc, spk*spr), gkl_b(nkc, spk*spr)
        real(kind=dp) :: y(spi*spj, spk*spr)
        integer(kind=ikind) :: a, b, c, d
        do b = 1, spj
            do a = 1, spi
                gij_a(:, a + (b-1)*spi) = ac(:, posI(a), posJ(b))
                gij_b(:, a + (b-1)*spi) = bc(:, posI(a), posJ(b))
            end do
        end do
        do d = 1, spr
            do c = 1, spk
                gkl_a(:, c + (d-1)*spk) = ac(:, posK(c), posR(d))
                gkl_b(:, c + (d-1)*spk) = bc(:, posK(c), posR(d))
            end do
        end do
        y = matmul(transpose(gij_a), gkl_b) + matmul(transpose(gij_b), gkl_a)
        do d = 1, spr
            do c = 1, spk
                do b = 1, spj
                    do a = 1, spi
                        zc(a, b, c, d) = 0.5_dp*y(a + (b-1)*spi, c + (d-1)*spk)
                    end do
                end do
            end do
        end do
    end subroutine coulomb_block

    ! =========================================================================
    ! Pair centres coincide (h < cutoffcentre): tabulated limits
    ! =========================================================================
    subroutine quartet_coincident(nq, q, lorder, l, m, n, gs, gc, p0mat, dx1, dy1, dz1, dx2, dy2, dz2, gi, gj, gk, gr, &
            zcontr, coeff1, coeff2, coeff3, coeff4, cutoff1, cutoff2, f)
        integer(kind=ikind), intent(in) :: nq, gi, gj, gk, gr, lorder
        real(kind=dp), intent(in) :: q(nq)
        real(kind=dp), intent(in) :: cutoff1, cutoff2
        real(kind=dp), intent(in) :: p0mat(:, :, :, :)
        real(kind=dp), intent(in) :: dx1(:, :, :), dy1(:, :, :), dz1(:, :, :), dx2(:, :, :), dy2(:, :, :), dz2(:, :, :)
        real(kind=dp), intent(in) :: zcontr(:, :, :, :)
        real(kind=dp), intent(in) :: coeff1(:), coeff2(:), coeff3(:), coeff4(:)
        integer(kind=ikind), intent(in) :: gs(:), gc(:), l(:), m(:), n(:)
        real(kind=dp), intent(out) :: f(nq)
        integer(kind=ikind) :: ll, mm, nn, llp, mmp, nnp, i, j, k, r, posi, posj, posk, posr
        real(kind=dp) :: ztot, mdl, mdm, mdn, mdlp, mdmp, mdnp, prod6, prod5, prod4, prod3, prod2, prod1
        real(kind=dp) :: h1, h1_base
        real(kind=dp) :: fzz(nq)

        f = 0.0_dp
        fzz = 0.0_dp
        posi = 1
        do i = gs(gi), gs(gi) + gc(gi) - 1
            posj = 1
            do j = gs(gj), gs(gj) + gc(gj) - 1
                posk = 1
                do k = gs(gk), gs(gk) + gc(gk) - 1
                    posr = 1
                    do r = gs(gr), gs(gr) + gc(gr) - 1
                        ztot = zcontr(posi, posj, posk, posr)*coeff1(posi)*coeff2(posj)*coeff3(posk)*coeff4(posr)
                        if (abs(ztot) < cutoff1) then
                            posr = posr + 1
                            cycle
                        end if
                        do ll = 0, l(i) + l(j)
                            mdl = dx1(ll + 1, l(i) + 1, l(j) + 1)
                            if (abs(mdl) < 1.0e-30_dp) cycle
                            prod1 = mdl*ztot
                            do mm = 0, m(i) + m(j)
                                mdm = dy1(mm + 1, m(i) + 1, m(j) + 1)
                                if (abs(mdm) < 1.0e-30_dp) cycle
                                prod2 = mdm*prod1
                                if (mod(ll + mm, 2) == 0) then
                                    h1_base = 1.0_dp
                                else
                                    h1_base = -1.0_dp
                                end if
                                do nn = 0, n(i) + n(j)
                                    if (mod(nn, 2) == 0) then
                                        h1 = h1_base
                                    else
                                        h1 = -h1_base
                                    end if
                                    mdn = dz1(nn + 1, n(i) + 1, n(j) + 1)
                                    if (abs(mdn) < 1.0e-30_dp) cycle
                                    prod3 = mdn*h1*prod2
                                    do llp = 0, l(k) + l(r)
                                        mdlp = dx2(llp + 1, l(k) + 1, l(r) + 1)
                                        if (abs(mdlp) < 1.0e-30_dp) cycle
                                        prod4 = mdlp*prod3
                                        do mmp = 0, m(k) + m(r)
                                            mdmp = dy2(mmp + 1, m(k) + 1, m(r) + 1)
                                            if (abs(mdmp) < 1.0e-30_dp) cycle
                                            prod5 = mdmp*prod4
                                            do nnp = 0, n(k) + n(r)
                                                mdnp = dz2(nnp + 1, n(k) + 1, n(r) + 1)
                                                prod6 = mdnp*prod5
                                                if (abs(prod6) > cutoff2) then
                                                    f = f + prod6*p0mat(:, ll + llp + 1, mm + mmp + 1, nn + nnp + 1)
                                                    if (lorder == 2) fzz = fzz &
                                                            + prod6*p0mat(:, ll + llp + 1, mm + mmp + 1, nn + nnp + 3)
                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                        posr = posr + 1
                    end do
                    posk = posk + 1
                end do
                posj = posj + 1
            end do
            posi = posi + 1
        end do
        if (lorder == 2) call combine_j2(q, f, fzz)
    end subroutine quartet_coincident

    ! =========================================================================
    ! Distinct pair centres: derivatives of j0(q h) through the Bessel module
    ! =========================================================================
    subroutine quartet_distant(mu, lorder, l, m, n, gs, gc, hx, hy, hz, h, dx1, dy1, dz1, dx2, dy2, dz2, gi, gj, gk, gr, &
            zcontr, coeff1, coeff2, coeff3, coeff4, cutoff1, cutoff2, hpre, f)
        integer(kind=ikind), intent(in) :: gi, gj, gk, gr, lorder
        integer(kind=ikind), intent(in) :: l(:), m(:), n(:), gs(:), gc(:)
        real(kind=dp), intent(in) :: cutoff1, cutoff2, hx, hy, hz, h
        real(kind=dp), intent(in) :: zcontr(:, :, :, :)
        real(kind=dp), intent(in) :: coeff1(:), coeff2(:), coeff3(:), coeff4(:)
        real(kind=dp), intent(in), allocatable :: mu(:)
        real(kind=dp), intent(in) :: dx1(:, :, :), dy1(:, :, :), dz1(:, :, :), dx2(:, :, :), dy2(:, :, :), dz2(:, :, :)
        real(kind=dp), intent(inout) :: hpre(:, :, :, :)
        real(kind=dp), intent(out) :: f(size(mu))

        integer(kind=ikind) :: posi, posj, posk, posr, i, j, k, r
        integer(kind=ikind) :: ll, mm, nn, llp, mmp, nnp, llmax
        real(kind=dp) :: ztot, mdn, mdl, mdm, mdlp, mdmp, mdnp
        real(kind=dp) :: prod1, prod2, prod3, prod4, prod5, prod6
        real(kind=dp) :: a(bdim, bdim), b(bdim, bdim), c(bdim, bdim)
        real(kind=dp) :: h_saved(bdim), h_saved2(bdim), bd(bdim), fzz(size(mu))
        integer(kind=ikind) :: lmx
        real(kind=dp) :: h1, h1_base

        llmax = l(gs(gi)) + m(gs(gi)) + n(gs(gi)) + l(gs(gj)) + m(gs(gj)) + n(gs(gj)) + &
                l(gs(gk)) + m(gs(gk)) + n(gs(gk)) + l(gs(gr)) + m(gs(gr)) + n(gs(gr))

        lmx = llmax + lorder
        call rrdj0(lmx, hx, a)
        call rrdj0(lmx, hy, b)
        call rrdj0(lmx, hz, c)

        do k = 0, lmx
            do j = 0, lmx - k
                do i = 0, lmx - k - j
                    bd = 0.0_dp
                    call BesselDeriv(bd, k, j, i, a, b, c, lmx)
                    hpre(:, k + 1, j + 1, i + 1) = bd
                end do
            end do
        end do

        h_saved = 0.0_dp
        h_saved2 = 0.0_dp
        posi = 1
        do i = gs(gi), gs(gi) + gc(gi) - 1
            posj = 1
            do j = gs(gj), gs(gj) + gc(gj) - 1
                posk = 1
                do k = gs(gk), gs(gk) + gc(gk) - 1
                    posr = 1
                    do r = gs(gr), gs(gr) + gc(gr) - 1
                        ztot = zcontr(posi, posj, posk, posr)*coeff1(posi)*coeff2(posj)*coeff3(posk)*coeff4(posr)
                        if (abs(ztot) < cutoff1) then
                            posr = posr + 1
                            cycle
                        end if
                        do ll = 0, l(i) + l(j)
                            mdl = dx1(ll + 1, l(i) + 1, l(j) + 1)
                            prod1 = mdl*ztot
                            do mm = 0, m(i) + m(j)
                                mdm = dy1(mm + 1, m(i) + 1, m(j) + 1)
                                prod2 = mdm*prod1
                                if (mod(ll + mm, 2) == 0) then
                                    h1_base = 1.0_dp
                                else
                                    h1_base = -1.0_dp
                                end if
                                do nn = 0, n(i) + n(j)
                                    if (mod(nn, 2) == 0) then
                                        h1 = h1_base
                                    else
                                        h1 = -h1_base
                                    end if
                                    mdn = dz1(nn + 1, n(i) + 1, n(j) + 1)
                                    prod3 = mdn*h1*prod2
                                    do llp = 0, l(k) + l(r)
                                        mdlp = dx2(llp + 1, l(k) + 1, l(r) + 1)
                                        prod4 = mdlp*prod3
                                        do mmp = 0, m(k) + m(r)
                                            mdmp = dy2(mmp + 1, m(k) + 1, m(r) + 1)
                                            prod5 = mdmp*prod4
                                            do nnp = 0, n(k) + n(r)
                                                mdnp = dz2(nnp + 1, n(k) + 1, n(r) + 1)
                                                prod6 = mdnp*prod5
                                                if (abs(prod6) > cutoff2) then
                                                    h_saved = h_saved + hpre(:, ll + llp + 1, mm + mmp + 1, nn + nnp + 1)*prod6
                                                    if (lorder == 2) h_saved2 = h_saved2 &
                                                            + hpre(:, ll + llp + 1, mm + mmp + 1, nn + nnp + 3)*prod6
                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                        posr = posr + 1
                    end do
                    posk = posk + 1
                end do
                posj = posj + 1
            end do
            posi = posi + 1
        end do
        ! sum_b h_saved(b+1) j_b(q h) (q/h)^b, evaluated stably for every q h
        call radial_series(size(mu, kind=ikind), mu, h, llmax, h_saved(1:llmax+1), f)
        if (lorder == 2) then
            call radial_series(size(mu, kind=ikind), mu, h, lmx, h_saved2(1:lmx+1), fzz)
            call combine_j2(mu, f, fzz)
        end if
    end subroutine quartet_distant

    ! =========================================================================
    ! j2(q r) P2(cos theta) from j0 by the solid-harmonic (Hobson) theorem:
    !   (3 d^2/dz^2 - Laplacian) j0(q r) = 2 q^2 j2(q r) P2(cos theta),  Laplacian j0 = -q^2 j0
    !   => j2 P2 = (3 d^2/dz^2 j0 + q^2 j0) / (2 q^2)
    ! f holds the j0 integral, fzz the same with two extra z derivatives.
    ! =========================================================================
    subroutine combine_j2(q, f, fzz)
        real(kind=dp), intent(in) :: q(:), fzz(:)
        real(kind=dp), intent(inout) :: f(:)
        where (q > 1.0e-5_dp)
            f = (3.0_dp*fzz + q*q*f) / (2.0_dp*q*q)
        elsewhere
            f = 0.0_dp
        end where
    end subroutine combine_j2

    ! =========================================================================
    ! Limits of the derivatives of j0(q h) at h = 0, any order:
    !   d^L/dx^L d^M/dy^M d^N/dz^N j0(q|h|) at 0 = <(i q n_x)^L (i q n_y)^M (i q n_z)^N>_n
    !   = (-1)^(S/2) q^S (L-1)!! (M-1)!! (N-1)!! / (S+1)!!   (L, M, N even; S = L+M+N), else 0
    ! =========================================================================
    subroutine p0_limits(p0, smax, q)
        real(kind=dp), intent(out) :: p0(:, :, :, :)
        integer(kind=ikind), intent(in) :: smax
        real(kind=dp), intent(in) :: q(:)
        integer(kind=ikind) :: i, j, k, s
        real(kind=dp) :: c
        p0 = 0.0_dp
        do k = 0, smax, 2
            do j = 0, smax, 2
                do i = 0, smax, 2
                    s = i + j + k
                    c = dfact(i - 1) * dfact(j - 1) * dfact(k - 1) / dfact(s + 1)
                    if (mod(s / 2, 2) == 1) c = -c
                    p0(:, i + 1, j + 1, k + 1) = c * q**s
                end do
            end do
        end do
    end subroutine p0_limits

    pure real(kind=dp) function dfact(n)
        integer(kind=ikind), intent(in) :: n
        integer(kind=ikind) :: k
        dfact = 1.0_dp
        k = n
        do while (k > 1)
            dfact = dfact * k
            k = k - 2
        end do
    end function dfact

end module xs_kernel_j0

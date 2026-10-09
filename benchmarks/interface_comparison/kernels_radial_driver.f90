program drv
    use Bessels_j0
    implicit none
    integer, parameter :: dp = selected_real_kind(15)
    integer(kind=selected_int_kind(8)) :: order, b, i
    real(dp), allocatable :: mu(:)
    real(dp) :: hs(20), res(1), xs(12)
    xs = [1d-3, 1d-2, 1d-1, 1d0, 1d1, 5d1, 9d1, 1.2d2, 2d2, 3d2, 5d2, 8d2]
    allocate(mu(1))
    do order = 1, 8
        do b = 0, order
            do i = 1, size(xs)
                hs = 0d0; hs(b + 1) = 1d0
                mu(1) = xs(i)                      ! q = x, H = 1
                call bessels0rr(res, order, mu, 1d0, hs)
                print '(i3,i3,es12.4,es25.16)', order, b, xs(i), res(1)
            end do
        end do
    end do
end program

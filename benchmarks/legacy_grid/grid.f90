program grid
    use omp_lib
    use qgrid3d
    implicit none
    integer, parameter   ::   dp = selected_real_kind(15)
    integer, parameter   ::   ikind = selected_int_kind(8)
    real(kind=dp)                ::   e0, thsmax
    integer                      ::   n1, n2, npts
    real(kind=dp), allocatable   ::   qxyz(:,:)
    real(kind=dp)                ::   qnorm
    integer                      ::   i 

    ! --- Definition of Input Parameters ---
    
    ! Mean-Photon Energy in keV
    e0 = 10.0_dp
    ! Maximum Scattering Angle in Degrees
    thsmax = 60.0_dp
    ! Number of Rings
    n1 = 49
    ! Number of Points on Each Ring
    n2 = 73
    
    ! Total Number of Grid Points
    npts = n1 * (n2 - 1) + 1
    
    ! --- Dynamical Allocation ---
    allocate(qxyz(npts, 3))

    call qcoord3(e0,thsmax,n1,n2,qxyz)

    do i = 1, npts
        qnorm = sqrt(qxyz(i,1)**2 + qxyz(i,2)**2 + qxyz(i,3)**2)
        print '(4f19.12)', qxyz(i,1), qxyz(i,2), qxyz(i,3), qnorm
    end do
    
    open(unit=15, file='qgrid.dat')
    do i = 1, npts
        write(15,'(3f19.12)') qxyz(i,1), qxyz(i,2), qxyz(i,3)
    end do
    close(15)

    deallocate(qxyz)

end program grid

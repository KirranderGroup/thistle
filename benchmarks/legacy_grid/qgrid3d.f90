module qgrid3d
    implicit none
    contains

        subroutine qcoord3(e0,thsmax,n1,n2,qxyz)
       
        integer, parameter   ::   dp = selected_real_kind(15)
        integer, parameter   ::   ikind = selected_int_kind(8)
        real(kind=dp), intent(in)                              ::   e0, thsmax
        integer, intent(in)                                    ::   n1, n2
        integer                                                ::   npts
        real(kind=dp)                                          ::   pi
        real(kind=dp)                                          ::   keV2au, k0
        real(kind=dp)                                          ::   qxymax, dqxy, dthq
        real(kind=dp), dimension(n1+1)                         ::   qxy
        real(kind=dp), dimension(n2-1)                         ::   thq
        real(kind=dp), dimension(n2-1)                         ::   qx, qy
        real(kind=dp), dimension(n1*(n2-1)+1,3), intent(out)   ::   qxyz(:,:)
        integer                                                ::   i, j, idx 

        ! Total Number of Grid Points
        
        npts = n1 * (n2 - 1) + 1
    
        ! --- Definition of Code Parameters ---
    
        pi = 4.0_dp * atan(1.0_dp) 
        keV2au = 0.2681727606190511_dp 
        k0 = e0 * keV2au;
        qxymax = k0 * sin(thsmax / 180_dp * pi)
        dqxy = qxymax / n1
        dthq = 2.0_dp * pi / (n2 - 1);
     
        ! --- Definition of Grid ---
    
        do i = 1, (n1+1)
            qxy(i) = real(i - 1, dp) * dqxy
        end do

        do j = 1, (n2-1)
            thq(j) = real(j - 1, dp) * dthq
        end do
    
        qx = sin(thq)
        qy = cos(thq)

        qxyz(1, :) = [ 0.0_dp, 0.0_dp, 0.0_dp ]

        idx = 2
        do i = 2, (n1+1)
            do j = 1, (n2-1)
                qxyz(idx, 1) = qx(j) * qxy(i)
                qxyz(idx, 2) = qy(j) * qxy(i)
                qxyz(idx, 3) = k0 * (1.0_dp - sqrt(1.0_dp - (qxyz(idx,1)**2 + qxyz(idx,2)**2) / (k0**2)))
                idx = idx + 1
            end do
        end do
    
       end subroutine qcoord3

end module qgrid3d

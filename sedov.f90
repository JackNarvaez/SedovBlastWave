module sedov
    implicit none
    public :: sedov_blast
    private :: funct_r2, funct_rho2, funct_rdr2, funct_rhodrho2

    contains

    subroutine sedov_blast(t, e0, gamm, rho0, ndim)
        real, intent(in) :: t, e0, gamm, rho0
        integer, intent(in) :: ndim

        real, parameter :: alpha = 0.851060
        integer, parameter :: npts = 1000
        real :: rshock, rhoshock
        real :: gamm1
        real :: gamm2
        real :: e
        real :: a_const
        real :: V, V2, V0, dV
        integer :: ndim2
        real, dimension(npts) :: r_ad
        real, dimension(npts) :: rho_ad

        real :: pa, pb, pc, pd, pe
        real :: alpha0, alpha1, alpha2, alpha3, alpha4, alpha5
        integer :: ii


        gamm1 = gamm + 1.
        gamm2 = gamm - 1.
        ndim2 = ndim + 2
        e = alpha*e0
        a_const = (e/rho0)**(1./ndim2)


        !
        ! Define paramethers
        ! 
        pa = ndim2*gamm1/4.0
        pb = gamm1 / gamm2
        pc = 0.5 * ndim2 * gamm
        pd = ndim2*gamm1 / (ndim2*gamm1 - 2.*(2.+ndim*gamm2))
        pe = 1. + 0.5*ndim*gamm2

        !
        ! Define powers
        !
        alpha0 = 2./ndim2
        alpha2 = - gamm2 / (2.*gamm2+ndim)
        alpha1 = ndim2*gamm / (2. + ndim*gamm2) * (2*ndim*(2.-gamm) / (gamm*ndim2**2) - alpha2)
        alpha3 = ndim / (2.*gamm2 + ndim)
        alpha4 = ndim2/(2.-gamm)*alpha1
        alpha5 = - 2./(2.-gamm)

        !
        ! Calculate radius of the shock front
        !
        rshock = funct_r2(t, a_const, ndim2)
        rhoshock = funct_rho2(rho0, gamm1, gamm2)

        V0 = 2.0/(ndim2*gamm)
        V2 = 4.0/(ndim2*gamm1)
        V = V0

        dV = (v2 - V0)/REAL(npts-1)  ! need to check the denominator
        do ii=1,npts
            r_ad(ii)  = funct_rdr2(V, pa, pb, pc, pd, pe, alpha0, alpha1, alpha2)*rshock
            rho_ad(ii) = funct_rhodrho2(V, gamm, pb, pc, pd, pe, alpha3, alpha4, alpha5)*rhoshock
            V = V0 + ii*dV
            print *, r_ad(ii), rho_ad(ii)
        enddo

    end subroutine sedov_blast

    real function funct_r2(t, a_const, ndim2)
        real, intent(in) :: t, a_const
        integer, intent(in) :: ndim2

        funct_r2 = a_const * t**(2./ndim2)
    end function funct_r2

    real function funct_rdr2(V, a, b, c, d, e, alpha0, alpha1, alpha2)
        real, intent(in) :: V, a, b, c, d, e, alpha0, alpha1, alpha2
        real :: x1, x2, x3
        real :: temp

        x1 = a*V 
        x2 = b*(c*V-1.)
        x3 = d*(1.-e*V)

        !print *, x2, x3, x4, alpha3, alpha4, alpha5
        temp = x1**alpha0 * x2**alpha2 * x3**alpha1
        funct_rdr2 = 1.0/temp

    end function funct_rdr2

    real function funct_rho2(rho1, gamm1, gamm2)
        real, intent(in) :: rho1, gamm1, gamm2

        funct_rho2 = (gamm1/gamm2) * rho1
    end function funct_rho2

    real function funct_rhodrho2(V, gamm, b, c, d, e, alpha3, alpha4, alpha5)
        real, intent(in) :: V, gamm, b, c, d, e, alpha3, alpha4, alpha5
        real :: x2, x3, x4

        x2 = b*(c*V-1.)
        x3 = d*(1.-e*V)
        x4 = b*(1.-c/gamm*V) 

        !print *, x2, x3, x4, alpha3, alpha4, alpha5
        funct_rhodrho2 = x2**alpha3 * x3**alpha4 * x4**alpha5

    end function funct_rhodrho2

end module sedov
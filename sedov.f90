!==============================================================================
! Module Sedov
! Calculate the analytical solution for the Sedov Blast wave in the standard
! case. Current vertion is valid only for a gamma value equal to 1.4
! Written by Narvaez J. 2023
!==============================================================================

module sedov
    implicit none
    public :: sedov_blast
    private :: funct_r2, funct_rho2, funct_rdr2, funct_rhodrho2

    contains

    subroutine sedov_blast(t, e0, gamm, rho0, ndim, rmax, prop)
        real, intent(in) :: t, e0, gamm, rho0, rmax
        integer, intent(in) :: ndim, prop

        real, parameter :: alpha = 0.851060 ! only valid for gamma = 1.4
        integer, parameter :: npts = 1000
        integer, parameter :: npts2 = 100
        real :: rshock, yshock, y0
        real :: gamm1
        real :: gamm2
        real :: e
        real :: a_const
        real :: V, V2, V0, dV
        real :: dr2
        real :: r, y
        integer :: ndim2

        real :: pa, pb, pc, pd, pe
        real :: alpha0, alpha1, alpha2, alpha3, alpha4, alpha5
        integer :: ii

        character(len=15) :: outfile ! Output file.

        gamm1 = gamm + 1.
        gamm2 = gamm - 1.
        ndim2 = ndim + 2
        e = alpha*e0
        a_const = (e/rho0)**(1./ndim2)

        !
        ! Define paramethers
        ! 
        pa = ndim2*gamm1/4.
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
        select case(prop)
        case(1) ! density
            yshock = funct_rho2(rho0, gamm1, gamm2)
            y0 = rho0
            outfile = "density.dat"
        case(2) !Pressure
            yshock = funct_p2(e, rshock, gamm1, ndim, ndim2)
            y0 = rho0*gamm2*e
            outfile = "pressure.dat"
        end select

        open(unit=1, file=outfile)

        V0 = 2.0/(ndim2*gamm)
        V2 = 4.0/(ndim2*gamm1)
        V = V0

        dV = (v2 - V0)/REAL(npts-1)

        select case(prop)
        case(1) ! density
            do ii=1,npts
                r = funct_rdr2(V, pa, pb, pc, pd, pe, alpha0, alpha1, alpha2)*rshock
                y = funct_rhodrho2(V, gamm, pb, pc, pd, pe, alpha3, alpha4, alpha5)*yshock
                V = V0 + ii*dV
                write(1,*) r, y
            enddo
        case(2) ! pressure
            do ii=1,npts
                r = funct_rdr2(V, pa, pb, pc, pd, pe, alpha0, alpha1, alpha2)*rshock
                y = funct_pdp2(V, gamm, pa, pb, pc, pd, pe, alpha0, alpha1, alpha4, alpha5, ndim)*yshock
                V = V0 + ii*dV
                write(1,*) r, y
            enddo
        end select
        
        if (rmax > rshock) then
            dr2 = (rmax-rshock)/REAL(npts2-1)
            do ii=1, npts2
                r = rshock + ii*dr2
                y = y0
                write(1,*) r, y0
            enddo
        end if

        close(1)
    
    end subroutine sedov_blast

    real function funct_r2(t, a_const, ndim2)
    ! Calculate radius of the shock front
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

        temp = x1**alpha0 * x2**alpha2 * x3**alpha1
        funct_rdr2 = 1.0/temp

    end function funct_rdr2

    real function funct_rho2(rho1, gamm1, gamm2)
    ! Calculate density at the shock front
        real, intent(in) :: rho1, gamm1, gamm2

        funct_rho2 = (gamm1/gamm2) * rho1
    end function funct_rho2

    real function funct_rhodrho2(V, gamm, b, c, d, e, alpha3, alpha4, alpha5)
        real, intent(in) :: V, gamm, b, c, d, e, alpha3, alpha4, alpha5
        real :: x2, x3, x4

        x2 = b*(c*V-1.)
        x3 = d*(1.-e*V)
        x4 = b*(1.-c/gamm*V) 

        funct_rhodrho2 = x2**alpha3 * x3**alpha4 * x4**alpha5

    end function funct_rhodrho2

    real function funct_p2(e, r2, gamm1, ndim, ndim2)
    ! Calculate pressure at the shock front
        real, intent(in) :: e, r2, gamm1
        integer, intent(in) :: ndim, ndim2

        funct_p2 = 8.*e/(ndim2*gamm1*r2**ndim)
    end function funct_p2

    real function funct_pdp2(V, gamm, a, b, c, d, e, alpha0, alpha1, alpha4, alpha5, ndim)
        real, intent(in) :: V, gamm, a, b, c, d, e, alpha0, alpha1, alpha4, alpha5
        integer, intent(in) :: ndim
        real :: x1, x3, x4

        x1 = a*V
        x3 = d*(1.-e*V)
        x4 = b*(1.-c/gamm*V) 

        funct_pdp2 = x1**(alpha0*ndim) * x3**(alpha4-2.*alpha1) * x4**(1.+alpha5)

    end function funct_pdp2

end module sedov
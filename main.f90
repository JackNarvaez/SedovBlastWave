program main
    use sedov, only:sedov_blast
    implicit none

    real :: t       = 1.0
    real :: e0      = 1.0
    real :: gamm    = 1.4
    real :: rho0    = 1.0
    real :: rmax    = 2.0
    integer, parameter :: ndim = 3

    call sedov_blast(t, e0, gamm, rho0, ndim, rmax)

end program main
!==============================================================================
! The Sedov Blast Wave
! Calculate the analytical solution for the Sedov Blast wave in the standard
! case. Current vertion is valid only for a gamma value equal to 1.4
! Written by Narvaez J. 2023
!==============================================================================

program main
    use sedov, only:sedov_blast
    implicit none

    real :: t
    real :: e0
    real :: rho0
    real :: rmax
    real :: gamm = 1.4
    integer, parameter :: ndim = 3

    print *, '--------------------------------------------------'
    print *, '-----        The 3D Sedov Blast Wave         -----' 
    print *, '--------------------------------------------------'
    print *, 'Enter time t:'
    read  *, t
    print *, 'Enter initial energy e0:'
    read  *, e0
    print *, 'Enter initial density:'
    read  *, rho0
    print *, 'Enter maximum radius:'
    read  *, rmax

    call sedov_blast(t, e0, gamm, rho0, ndim, rmax)

    print *, '--------------------------------------------------'
    print *, '-----                  End                   -----' 
    print *, '--------------------------------------------------'

end program main
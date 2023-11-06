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
    integer :: prop
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
    print 2

    2 format( 'Analytical Solution:' ,/, &
              '  1) Density' ,/, &
              '  2) Pressure', /, &
              'Property:')
    read  *, prop

    call sedov_blast(t, e0, gamm, rho0, ndim, rmax, prop)

    ! Write Gnuplot file

    open(unit=3, file='plot.gp')

    write(3, *) 'set xlabel "r"'
    
    select case(prop)
    case(1)
        write(3, *) 'set ylabel "{/Symbol r}"'
        write(3, *) 'plot [0:', rmax, '] [-0.1:*] "density.dat" using 1:2 with lines title "{/Symbol r}(r)"'
        write(3, *) 'set terminal png'
    case(2)
        write(3, *) 'set ylabel "P"'
        write(3, *) 'plot [0:', rmax, '] [0:*] "pressure.dat" using 1:2 with lines title "P(r)"'
        write(3, *) 'set terminal png'
    end select

    close(3)

    print *, '--------------------------------------------------'
    print *, '-----                  End                   -----' 
    print *, '--------------------------------------------------'

end program main
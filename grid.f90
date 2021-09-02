program enclosure
! this is a program for 2D-simulation of the rayleigh-bernard convection in an incompressible plane Poiseullie flow.
! this is a sequential program in fortran-90.
! this program uses the modified SMAC scheme with Rhie-Chow momentum interpolation to maintain the pressure-velocity coupling on the grid points.
! this program uses the Kuwahara 3rd order upwind scheme for convective terms.
! 
!
implicit none
!
!
integer :: i,j                                !local loop variables
double precision :: dx,dy
double precision :: Pi
double precision, parameter :: a = 0.5d0
double precision, parameter :: delta = 9.00d0
!************************************************************************************
integer, parameter :: Nx=150,Ny=150
!************************************************************************************
integer ,parameter   :: nstep =0   
!************************************************************************************
double precision, allocatable, dimension(:) :: x,y,Xi,Eta
!************************************************************************************
   allocate( x(Nx) )
   allocate( Xi(Nx) )
   allocate( y(Ny) )
   allocate( Eta(Ny) )
!***********************************************************************************
! the computational domain generation takes place here.
!
               Xi(1)=0.0d0
               Eta(1)=0.0d0

               Xi(Nx)=1.0d0
               Eta(Ny)=1.0d0

! calculating the number of mesh points.
               dx =( Xi(Nx)-Xi(1) )/dble(Nx-1) 
               dy =( Eta(Ny)-Eta(1) )/dble(Ny-1)

       do i=1,Nx-1
   Xi(i+1) = Xi(i)+dx
       end do
 
       do j=1,Ny-1
   Eta(j+1) = Eta(j)+dy
       end do

!***********************************************************************************
! now mapping the computational domain in the physical domain using the stretching function.
! using the tan-hyperbolic mapping.
    do i=1,Nx
        x(i)  = a*(1.0d0+( dtanh(0.5d0*delta*(Xi(i)-a))/dtanh(0.5d0*a*delta) ) )
   end do
    do j=1,Ny
        y(j)  = a*(1.0d0+( dtanh(0.5d0*delta*(Eta(j)-a))/dtanh(0.5d0*a*delta) ) )
   end do
! now we save the mesh data in ASCII tecplot format.

       open(unit=2,file='2D_enclosure.tp',status='unknown')
 write(2,03) 'title ="',nstep,'"'
 write(2,*) 'VARIABLES = "x","y" '
 write(2,04) 'ZONE i=',Nx,',j=',Ny,',DATAPACKING="POINT"'
         do i=1,Nx
              write(2,*)
                  do j=1,Ny
                      write(2,87) x(i),y(j)
                  end do
         end do 
 close(2)
03 format(A,I5,A)
04 format(A,I5,A,I5,A)
87 format(1X,F22.10,1X,F22.10)
!***********************************************************************************
    open(unit=3,file='mesh.xy',status='unknown')

    write(3,*) Nx,Ny           

          do i=1,Nx
       write(3,33) i,x(i)
          end do
         
             do j=1,Ny
      write(3,33) j,y(j)
             end do

    close(3)

33 format(1X,I4,1X,F23.15)
!*************finally the clean-up act********************************************
!***********************************************************************************
deallocate(x)
deallocate(y)
deallocate(Xi)
deallocate(Eta)
!***********************************************************************************
end program 

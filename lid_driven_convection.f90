!**********************************************
 program incompressible
!**********************************************
!this program treat the viscous terms implicitly.
!this program uses the correction-pressure based modified SMAC scheme (refer Hasan & Sanghi,JHT,126,963-984) alongwith QUICK-upwind to deal convective-terms.
!this program uses Rhie-Chow momentum interpolation to deal with the grid-scale oscillations due to CPPE.
!this program uses fourth order accuracy in diffusion terms.
!**********************************************
  implicit none
!**********************************************
 logical, parameter :: restart = .false.
 integer, parameter :: maxstep = 15000000
 double precision, parameter :: dt=1.0d-05
 integer, parameter :: maxiter = 10000
 integer, parameter :: snapshot = 10000
 integer, parameter :: fil_cnt = 1
 double precision, parameter :: errtol=1.0d-05
 integer :: Nx,Ny,nstep
 integer :: i,j,jcnt
 double precision, parameter :: Re = 10000.0d0
 double precision, parameter :: Pr = 7.1d0
 double precision, parameter :: Ri = 50.0d0
 double precision, parameter :: kappa = 0.0d0
 double precision :: xcell_Pe,ycell_Pe,Pi
 double precision, allocatable, dimension(:) :: x,y,de,dw,dee,dww,dn,ds,dnn,dss
 double precision, allocatable, dimension(:) :: etaEE_2cen,etaE_2cen,etaW_2cen,etaWW_2cen
 double precision, allocatable, dimension(:) :: etaNN_2cen,etaN_2cen,etaS_2cen,etaSS_2cen
 double precision, allocatable, dimension(:) :: etaEE_1cen,etaE_1cen,etaW_1cen,etaWW_1cen
 double precision, allocatable, dimension(:) :: etaNN_1cen,etaN_1cen,etaS_1cen,etaSS_1cen
 double precision, allocatable, dimension(:) :: etaE_gradP,etaW_gradP,etaN_gradP,etaS_gradP
 double precision, allocatable, dimension(:) :: etaEE_upw_pos,etaE_upw_pos,etaW_upw_pos,etaWW_upw_pos
 double precision, allocatable, dimension(:) :: etaNN_upw_pos,etaN_upw_pos,etaS_upw_pos,etaSS_upw_pos
 double precision, allocatable, dimension(:) :: etaEE_upw_neg,etaE_upw_neg,etaW_upw_neg,etaWW_upw_neg
 double precision, allocatable, dimension(:) :: etaNN_upw_neg,etaN_upw_neg,etaS_upw_neg,etaSS_upw_neg
 double precision, allocatable, dimension(:,:) :: Avp,Ave,Avw,Avn,Avs,Avee,Avww,Avnn,Avss
 double precision, allocatable, dimension(:,:) :: Atp,Ate,Atw,Atn,Ats,Atee,Atww,Atnn,Atss
 double precision :: dx_int,dy_int,x_conv,y_conv,xintPLUS,xintMINUS,yintPLUS,yintMINUS
 double precision :: Ap,Ae,Aw,An,As,Aee,Aww,Ann,Ass
 double precision :: de_bnd,dw_bnd,dee_bnd,dww_bnd,dn_bnd,ds_bnd,dnn_bnd,dss_bnd
 double precision :: etaE_bnd2,etaEE_bnd2
 double precision :: etaW_bnd2,etaWW_bnd2
 double precision :: etaN_bnd2,etaNN_bnd2
 double precision :: etaS_bnd2,etaSS_bnd2
 double precision :: etaE,etaEE
 double precision :: etaW,etaWW
 double precision :: etaN,etaNN
 double precision :: etaS,etaSS
 double precision, allocatable, dimension(:,:) :: App,Ape,Apw,Apn,Aps
 double precision, allocatable, dimension(:,:,:) :: Ut,forc,Q,conv
 double precision, allocatable, dimension(:,:) :: press,pcorr,Qp
 double precision, dimension(2) :: gradP
 double precision :: dx1,dx2,dy1,dy2,Nuss_e,Nuss_w,Nuss_s,dTw_by_dx,dTe_by_dx,dTs_by_dy,tmpp
 double precision :: time,aa,at
 character(len=30) :: filname
 integer :: waste
!**********************************************
      Pi = 4.0d0*datan(1.0d0)
 print*, Pi
!**********************************************
!*****read the grid-data
  open(unit=1,file='mesh.xy',status='unknown')

  read(1,*) Nx,Ny

  print*, Nx,Ny
!**************allocate the flow variables***************

   allocate( x(Nx) )
   allocate( y(Ny) )
   allocate( de(Nx) )
   allocate( dw(Nx) )
   allocate( dee(Nx) )
   allocate( dww(Nx) )
   allocate( dn(Ny) )
   allocate( ds(Ny) )
   allocate( dnn(Ny) )
   allocate( dss(Ny) )
   allocate( Ut(Nx,Ny,5) )
   allocate( press(Nx,Ny) )
   allocate( pcorr(Nx,Ny) )
   allocate( forc(Nx,Ny,2) )
   allocate( Q(Nx,Ny,7) )
   allocate( conv(Nx,Ny,3) )
   allocate( Qp(Nx,Ny) )
   allocate( etaEE_2cen(Nx) )
   allocate( etaE_2cen(Nx) )
   allocate( etaW_2cen(Nx) )
   allocate( etaWW_2cen(Nx) )
   allocate( etaNN_2cen(Ny) )
   allocate( etaN_2cen(Ny) )
   allocate( etaS_2cen(Ny) )
   allocate( etaSS_2cen(Ny) )
   allocate( etaEE_1cen(Nx) )
   allocate( etaE_1cen(Nx) )
   allocate( etaW_1cen(Nx) )
   allocate( etaWW_1cen(Nx) )
   allocate( etaNN_1cen(Ny) )
   allocate( etaN_1cen(Ny) )
   allocate( etaS_1cen(Ny) )
   allocate( etaSS_1cen(Ny) )
   allocate( etaEE_upw_pos(Nx) )
   allocate( etaE_upw_pos(Nx) )
   allocate( etaW_upw_pos(Nx) )
   allocate( etaWW_upw_pos(Nx) )
   allocate( etaNN_upw_pos(Ny) )
   allocate( etaN_upw_pos(Ny) )
   allocate( etaS_upw_pos(Ny) )
   allocate( etaSS_upw_pos(Ny) )
   allocate( etaEE_upw_neg(Nx) )
   allocate( etaE_upw_neg(Nx) )
   allocate( etaW_upw_neg(Nx) )
   allocate( etaWW_upw_neg(Nx) )
   allocate( etaNN_upw_neg(Ny) )
   allocate( etaN_upw_neg(Ny) )
   allocate( etaS_upw_neg(Ny) )
   allocate( etaSS_upw_neg(Ny) )
   allocate( etaE_gradP(Nx) )
   allocate( etaW_gradP(Nx) )
   allocate( etaN_gradP(Ny) )
   allocate( etaS_gradP(Ny) )
   allocate( App(Nx,Ny) )
   allocate( Ape(Nx,Ny) )
   allocate( Apw(Nx,Ny) )
   allocate( Apn(Nx,Ny) )
   allocate( Aps(Nx,Ny) )
   allocate( Avp(Nx,Ny) )
   allocate( Ave(Nx,Ny) )
   allocate( Avw(Nx,Ny) )
   allocate( Avn(Nx,Ny) )
   allocate( Avs(Nx,Ny) )
   allocate( Avee(Nx,Ny) )
   allocate( Avww(Nx,Ny) )
   allocate( Avnn(Nx,Ny) )
   allocate( Avss(Nx,Ny) )

   allocate( Atp(Nx,Ny) )
   allocate( Ate(Nx,Ny) )
   allocate( Atw(Nx,Ny) )
   allocate( Atn(Nx,Ny) )
   allocate( Ats(Nx,Ny) )
   allocate( Atee(Nx,Ny) )
   allocate( Atww(Nx,Ny) )
   allocate( Atnn(Nx,Ny) )
   allocate( Atss(Nx,Ny) )
!********************************************************
     do i=1,Nx
   read(1,11) waste,x(i)
     end do

     do j=1,Ny
   read(1,11) waste,y(j)
     end do

  close(1)
11 format(1X,I4,1X,F23.15)

 print*, "grid data read"
!********************************************************
!*****calculate the step-sizes

       do i=2,Nx-1
   de(i) = x(i+1)-x(i)
   dw(i) = x(i)-x(i-1)
  if(i.ge.3.and.i.le.(Nx-2))then
  dee(i) = x(i+2)-x(i)
  dww(i) = x(i)-x(i-2)
  else
  dee(i) = 0.0d0
  dww(i) = 0.0d0
  end if
       end do

       do j=2,Ny-1
   dn(j) = y(j+1)-y(j)
   ds(j) = y(j)-y(j-1)
  if(j.ge.3.and.j.le.(Ny-2))then
  dnn(j) = y(j+2)-y(j)
  dss(j) = y(j)-y(j-2)
  else
  dnn(j) = 0.0d0
  dss(j) = 0.0d0
  end if
       end do

!********************************************************
!*****calculate the weights for different discretization schemes
!********************************************************
!
!*****for one-sided approximations of the first and second-order derivative at boundaries
!****************************************
         de_bnd = x(2)-x(1)
        dee_bnd = x(3)-x(1)
!****************************************
       etaE = dee_bnd/(de_bnd*(dee_bnd-de_bnd))
      etaEE = (-1.0d0*de_bnd)/(dee_bnd*(dee_bnd-de_bnd))

       etaE_bnd2 = -2.0d0/(de_bnd*(dee_bnd-de_bnd))
      etaEE_bnd2 = 2.0d0/(dee_bnd*(dee_bnd-de_bnd))
!****************************************
!****************************************
         dw_bnd = x(Nx)-x(Nx-1)
        dww_bnd = x(Nx)-x(Nx-2)
!****************************************
       etaW = (-1.0d0*dww_bnd)/(dw_bnd*(dww_bnd-dw_bnd))
      etaWW = dw_bnd/(dww_bnd*(dww_bnd-dw_bnd))

       etaW_bnd2 = -2.0d0/(dw_bnd*(dww_bnd-dw_bnd))
      etaWW_bnd2 = 2.0d0/(dww_bnd*(dww_bnd-dw_bnd))
!****************************************
!****************************************
        dn_bnd = y(2)-y(1)
       dnn_bnd = y(3)-y(1)
!****************************************
       etaN = dnn_bnd/(dn_bnd*(dnn_bnd-dn_bnd))
      etaNN = (-1.0d0*dn_bnd)/(dnn_bnd*(dnn_bnd-dn_bnd))

       etaN_bnd2 = -2.0d0/(dn_bnd*(dnn_bnd-dn_bnd))
      etaNN_bnd2 = 2.0d0/(dnn_bnd*(dnn_bnd-dn_bnd))
!****************************************
!****************************************
        ds_bnd = y(Ny)-y(Ny-1)
       dss_bnd = y(Ny)-y(Ny-2)
!****************************************
       etaS = (-1.0d0*dss_bnd)/(ds_bnd*(dss_bnd-ds_bnd))
      etaSS = ds_bnd/(dss_bnd*(dss_bnd-ds_bnd))

       etaS_bnd2 = -2.0d0/(ds_bnd*(dss_bnd-ds_bnd))
      etaSS_bnd2 = 2.0d0/(dss_bnd*(dss_bnd-ds_bnd))
!****************************************
!****************************************
!
!*****for Second derivative central-approximations
!***close all 4th order discrete-molecules for time-being.
       do i=2,Nx-1
!*****
  if(i.ge.4.and.i.le.(Nx-3))then
   etaWW_2cen(i) = -2.0d0*( (dee(i)*dw(i))-(dee(i)*de(i))+(de(i)*dw(i)) )/( dww(i)*(dww(i)-dw(i))*(dww(i)+dee(i))*(dww(i)+de(i)) )
    etaW_2cen(i) = 2.0d0*( (dee(i)*dww(i))-(dee(i)*de(i))+(de(i)*dww(i)) )/( dw(i)*(dww(i)-dw(i))*(dw(i)+dee(i))*(dw(i)+de(i)) )
    etaE_2cen(i) = 2.0d0*( (dee(i)*dw(i))-(dww(i)*dw(i))+(dee(i)*dww(i)) )/( de(i)*(dee(i)-de(i))*(dww(i)+de(i))*(dw(i)+de(i)) )
   etaEE_2cen(i) = -2.0d0*( (de(i)*dw(i))-(dww(i)*dw(i))+(de(i)*dww(i)) )/( dee(i)*(dee(i)-de(i))*(dw(i)+dee(i))*(dww(i)+dee(i)) )
  else
   etaEE_2cen(i) = 0.0d0
    etaE_2cen(i) = 2.0d0/( de(i)*( de(i)+dw(i) ) )
    etaW_2cen(i) = 2.0d0/( dw(i)*( de(i)+dw(i) ) )
   etaWW_2cen(i) = 0.0d0
  end if
       end do
!**********************
       do j=2,Ny-1
!*****
  if(j.ge.4.and.j.le.(Ny-3))then
   etaSS_2cen(j) = -2.0d0*( (dnn(j)*ds(j))-(dnn(j)*dn(j))+(dn(j)*ds(j)) )/( dss(j)*(dss(j)-ds(j))*(dss(j)+dnn(j))*(dss(j)+dn(j)) )
    etaS_2cen(j) = 2.0d0*( (dnn(j)*dss(j))-(dnn(j)*dn(j))+(dn(j)*dss(j)) )/( ds(j)*(dss(j)-ds(j))*(ds(j)+dnn(j))*(ds(j)+dn(j)) )
    etaN_2cen(j) = 2.0d0*( (dnn(j)*ds(j))-(dss(j)*ds(j))+(dnn(j)*dss(j)) )/( dn(j)*(dnn(j)-dn(j))*(dss(j)+dn(j))*(ds(j)+dn(j)) )
   etaNN_2cen(j) = -2.0d0*( (dn(j)*ds(j))-(dss(j)*ds(j))+(dn(j)*dss(j)) )/( dnn(j)*(dnn(j)-dn(j))*(ds(j)+dnn(j))*(dss(j)+dnn(j)) )
  else
   etaNN_2cen(j) = 0.0d0
    etaN_2cen(j) = 2.0d0/( dn(j)*( dn(j)+ds(j) ) )
    etaS_2cen(j) = 2.0d0/( ds(j)*( dn(j)+ds(j) ) )
   etaSS_2cen(j) = 0.0d0
  end if
       end do
!**********************

!*****for First derivative central-approximations
       do i=2,Nx-1
  if(i.ge.3.and.i.le.(Nx-2))then
   etaWW_1cen(i) = ( dee(i)*dw(i)*de(i) )/( dww(i)*(dww(i)-dw(i))*(dww(i)+dee(i))*(dww(i)+de(i)) )
    etaW_1cen(i) = ( -dee(i)*dww(i)*de(i) )/( dw(i)*(dww(i)-dw(i))*(dw(i)+dee(i))*(dw(i)+de(i)) )
    etaE_1cen(i) = ( dee(i)*dw(i)*dww(i) )/( de(i)*(dee(i)-de(i))*(dww(i)+de(i))*(dw(i)+de(i)) )
   etaEE_1cen(i) = ( -de(i)*dw(i)*dww(i) )/( dee(i)*(dee(i)-de(i))*(dw(i)+dee(i))*(dww(i)+dee(i)) )
  end if
       end do
!**********************
       do j=2,Ny-1
  if(j.ge.3.and.j.le.(Ny-2))then
   etaSS_1cen(j) = ( dnn(j)*ds(j)*dn(j) )/( dss(j)*(dss(j)-ds(j))*(dss(j)+dnn(j))*(dss(j)+dn(j)) )
    etaS_1cen(j) = ( -dnn(j)*dss(j)*dn(j) )/( ds(j)*(dss(j)-ds(j))*(ds(j)+dnn(j))*(ds(j)+dn(j)) )
    etaN_1cen(j) = ( dnn(j)*ds(j)*dss(j) )/( dn(j)*(dnn(j)-dn(j))*(dss(j)+dn(j))*(ds(j)+dn(j)) )
   etaNN_1cen(j) = ( -dn(j)*ds(j)*dss(j) )/( dnn(j)*(dnn(j)-dn(j))*(ds(j)+dnn(j))*(dss(j)+dnn(j)) )
  end if
       end do
!**********************
       do i=2,Nx-1
    etaE_gradP(i) = dw(i)/( de(i)*( de(i)+dw(i) ) )
    etaW_gradP(i) = -de(i)/( dw(i)*( de(i)+dw(i) ) )
       end do
!**********************
       do j=2,Ny-1
    etaN_gradP(j) = ds(j)/( dn(j)*( dn(j)+ds(j) ) )
    etaS_gradP(j) = -dn(j)/( ds(j)*( dn(j)+ds(j) ) )
       end do
!*****for First derivative upwind-approximations
!***for up>0.0
       do i=2,Nx-1
  if(i.ge.3.and.i.le.(Nx-2))then
   etaWW_upw_pos(i) = ( dw(i)*de(i) )/( dww(i)*(dww(i)-dw(i))*(dww(i)+de(i)) )
    etaW_upw_pos(i) = ( -dww(i)*de(i) )/( dw(i)*(dww(i)-dw(i))*(dw(i)+de(i)) )
    etaE_upw_pos(i) = ( dw(i)*dww(i) )/( de(i)*(dww(i)+de(i))*(dw(i)+de(i)) )
   etaEE_upw_pos(i) = 0.0d0
  end if
       end do
!***for vp>0.0
       do j=2,Ny-1
  if(j.ge.3.and.j.le.(Ny-2))then
   etaSS_upw_pos(j) = ( ds(j)*dn(j) )/( dss(j)*(dss(j)-ds(j))*(dss(j)+dn(j)) )
    etaS_upw_pos(j) = ( -dss(j)*dn(j) )/( ds(j)*(dss(j)-ds(j))*(ds(j)+dn(j)) )
    etaN_upw_pos(j) = ( ds(j)*dss(j) )/( dn(j)*(dss(j)+dn(j))*(ds(j)+dn(j)) )
   etaNN_upw_pos(j) = 0.0d0
  end if
       end do
!***for up<0.0
       do i=2,Nx-1
  if(i.ge.3.and.i.le.(Nx-2))then
   etaWW_upw_neg(i) = 0.0d0
    etaW_upw_neg(i) = ( -dee(i)*de(i) )/( dw(i)*(dw(i)+dee(i))*(dw(i)+de(i)) )
    etaE_upw_neg(i) = ( dee(i)*dw(i) )/( de(i)*(dee(i)-de(i))*(dw(i)+de(i)) )
   etaEE_upw_neg(i) = ( -de(i)*dw(i) )/( dee(i)*(dee(i)-de(i))*(dw(i)+dee(i)) )
  end if
       end do
!***for vp<0.0
       do j=2,Ny-1
  if(j.ge.3.and.j.le.(Ny-2))then
   etaSS_upw_neg(j) = 0.0d0
    etaS_upw_neg(j) = ( -dnn(j)*dn(j) )/( ds(j)*(ds(j)+dnn(j))*(ds(j)+dn(j)) )
    etaN_upw_neg(j) = ( dnn(j)*ds(j) )/( dn(j)*(dnn(j)-dn(j))*(ds(j)+dn(j)) )
   etaNN_upw_neg(j) = ( -dn(j)*ds(j) )/( dnn(j)*(dnn(j)-dn(j))*(ds(j)+dnn(j)) )
  end if
       end do
!********************************************************
!*****calculate the Stencil-parameters for different elliptic-equations.
!********************************************************
!*****for the transport-equation of velocity
   aa = dt/Re
   print*, aa
!***use second order central discretization
      do i=2,Nx-1
      do j=2,Ny-1
!******************************
   Avp(i,j) = 1.0d0+( aa*( etaE_2cen(i)+etaW_2cen(i)+etaN_2cen(j)+etaS_2cen(j) ) )
   Ave(i,j) = (-aa)*etaE_2cen(i)
   Avw(i,j) = (-aa)*etaW_2cen(i)
   Avn(i,j) = (-aa)*etaN_2cen(j)
   Avs(i,j) = (-aa)*etaS_2cen(j)
  Avee(i,j) = (-aa)*etaEE_2cen(i)
  Avww(i,j) = (-aa)*etaWW_2cen(i)
  Avnn(i,j) = (-aa)*etaNN_2cen(j)
  Avss(i,j) = (-aa)*etaSS_2cen(j)
!******************************
       end do
       end do
!******************************
   at = dt/(Re*Pr)
   print*, at
!***use second order central discretization
      do i=2,Nx-1
      do j=2,Ny-1
!******************************
   Atp(i,j) = 1.0d0+( at*( etaE_2cen(i)+etaW_2cen(i)+etaN_2cen(j)+etaS_2cen(j) ) )
   Ate(i,j) = (-at)*etaE_2cen(i)
   Atw(i,j) = (-at)*etaW_2cen(i)
   Atn(i,j) = (-at)*etaN_2cen(j)
   Ats(i,j) = (-at)*etaS_2cen(j)
  Atee(i,j) = (-at)*etaEE_2cen(i)
  Atww(i,j) = (-at)*etaWW_2cen(i)
  Atnn(i,j) = (-at)*etaNN_2cen(j)
  Atss(i,j) = (-at)*etaSS_2cen(j)
!******************************
       end do
       end do
!*****for the poisson-equation of correction-pressure
      do i=2,Nx-1
      do j=2,Ny-1
!******************************
     Ape(i,j) = 2.0d0/(de(i)*(de(i)+dw(i)))
     Apw(i,j) = 2.0d0/(dw(i)*(de(i)+dw(i)))
     Apn(i,j) = 2.0d0/(dn(j)*(dn(j)+ds(j)))
     Aps(i,j) = 2.0d0/(ds(j)*(dn(j)+ds(j)))
     App(i,j) = (-2.0d0/(de(i)*dw(i)))+(-2.0d0/(dn(j)*ds(j)))
!******************************
       end do
       end do
!******************************
   if(restart.eq..true.)then   !***!
!******************************
!**************************************************
  print*,'restart is true so you might have faced a power failure.'
!**************************************************
 open(unit=1,file='2D.tp',status='unknown')
  read(1,23) time,nstep
    read(1,*)
    read(1,*)

    do i=1,Ny
    read(1,*)
           do j=1,Nx
 read(1,24) x(i),y(j),Ut(i,j,1),Ut(i,j,2),Ut(i,j,3),Ut(i,j,4),Ut(i,j,5),press(i,j),pcorr(i,j)
           end do
    end do

   close(1)
23 format(8X,F10.6,I8,1X)
24 format(2(1X,F10.7),7(1X,E14.7))
!**************************************************
   else  !***!
!**************************************************
   time=0.0d0
   nstep=0
   call initialize()   !the subroutine for initializing the flow-field is called here.
   call savesnapshot()
!**************************************************
    end if  !***!
!**********************************************************************************
!**********
  do while(nstep.lt.maxstep)

    nstep=nstep+1
    time=time+dt

        call convective()

!**********

        call solve_transport()

!**********

        call rhei_chow()

!**********

        call pressure_poisson()

!**********

        call update_primitives()

!**********

        call boundary_condition()

        call divergence()

!***compute the local Nusselt number at western boundary
    tmpp = 0.0d0

    do j=2,Ny-1
!************************************
   yintPLUS = ( y(j+1)+y(j) )/2.0d0
  yintMINUS = ( y(j)+y(j-1) )/2.0d0

     dy_int = yintPLUS-yintMINUS
!************************************

    dTw_by_dx = etaE*( Ut(2,j,3)-Ut(1,j,3) )+etaEE*( Ut(3,j,3)-Ut(1,j,3) )
    tmpp = tmpp+(dTw_by_dx*dy_int)

    end do

    Nuss_w = -tmpp/( (y(Ny-1)-y(2))*(Ny-2) )

!***compute the local Nusselt number at eastern boundary
    tmpp = 0.0d0

    do j=2,Ny-1
!************************************
   yintPLUS = ( y(j+1)+y(j) )/2.0d0
  yintMINUS = ( y(j)+y(j-1) )/2.0d0

     dy_int = yintPLUS-yintMINUS
!************************************

    dTe_by_dx = etaW*( Ut(Nx,j,3)-Ut(Nx-1,j,3) )+etaWW*( Ut(Nx,j,3)-Ut(Nx-2,j,3) )
    tmpp = tmpp+(dTe_by_dx*dy_int)

    end do

    Nuss_e = -tmpp/( (y(Ny-1)-y(2))*(Ny-2) )

!***compute the local Nusselt number at southern boundary
    tmpp = 0.0d0

    do i=2,Nx-1
!************************************
   xintPLUS = ( x(i+1)+x(i) )/2.0d0
  xintMINUS = ( x(i)+x(i-1) )/2.0d0

     dx_int = xintPLUS-xintMINUS
!************************************

    dTs_by_dy = etaS*( Ut(i,2,3)-Ut(i,1,3) )+etaSS*( Ut(i,3,3)-Ut(i,1,3) )
    tmpp = tmpp+(dTs_by_dy*dx_int)

    end do

    Nuss_s = -tmpp/( (x(Nx)-x(1))*Nx )

!************************************

    if(mod(nstep,fil_cnt)==0)then
!*****write the file for time history*****
    open(unit=2,file='signal.dat',status='unknown',access='append')
    write(2,21) time,Ut(63,63,1),Ut(63,63,2),Ut(63,63,3),Ut(177,177,1),Ut(177,177,2),Ut(177,177,3),Ut(63,177,1),Ut(63,177,2),Ut(63,177,3),Ut(177,63,1),Ut(177,63,2),Ut(177,63,3),Ut(101,101,1),Ut(101,101,2),Ut(101,101,3)
    close(2)

21 format(16(1X,F12.7))

    open(unit=3,file='Nusselt.dat',status='unknown',access='append')
    write(3,33) time,Nuss_s,Nuss_w,Nuss_e
    close(3)

33 format(4(1X,F10.5))
!*****close file write*****
    end if

!**********
    if(mod(nstep,snapshot)==0)then
    call savesnapshot()
    end if

  end do

!********************************************************
!*****finally the clean-up act

   deallocate(x)
   deallocate(y)
   deallocate(de)
   deallocate(dw)
   deallocate(dee)
   deallocate(dww)
   deallocate(dn)
   deallocate(ds)
   deallocate(dnn)
   deallocate(dss)
   deallocate(Ut)
   deallocate(press)
   deallocate(pcorr)
   deallocate(forc)
   deallocate(Q)
   deallocate(conv)
   deallocate(Qp)
   deallocate(etaEE_2cen)
   deallocate(etaE_2cen)
   deallocate(etaW_2cen)
   deallocate(etaWW_2cen)
   deallocate(etaNN_2cen)
   deallocate(etaN_2cen)
   deallocate(etaS_2cen)
   deallocate(etaSS_2cen)
   deallocate(etaEE_1cen)
   deallocate(etaE_1cen)
   deallocate(etaW_1cen)
   deallocate(etaWW_1cen)
   deallocate(etaNN_1cen)
   deallocate(etaN_1cen)
   deallocate(etaS_1cen)
   deallocate(etaSS_1cen)
   deallocate(etaEE_upw_pos)
   deallocate(etaE_upw_pos)
   deallocate(etaW_upw_pos)
   deallocate(etaWW_upw_pos)
   deallocate(etaNN_upw_pos)
   deallocate(etaN_upw_pos)
   deallocate(etaS_upw_pos)
   deallocate(etaSS_upw_pos)
   deallocate(etaEE_upw_neg)
   deallocate(etaE_upw_neg)
   deallocate(etaW_upw_neg)
   deallocate(etaWW_upw_neg)
   deallocate(etaNN_upw_neg)
   deallocate(etaN_upw_neg)
   deallocate(etaS_upw_neg)
   deallocate(etaSS_upw_neg)
   deallocate(etaE_gradP)
   deallocate(etaW_gradP)
   deallocate(etaN_gradP)
   deallocate(etaS_gradP)
   deallocate(App)
   deallocate(Ape)
   deallocate(Apw)
   deallocate(Apn)
   deallocate(Aps)
   deallocate(Avp)
   deallocate(Ave)
   deallocate(Avw)
   deallocate(Avn)
   deallocate(Avs)
   deallocate(Avee)
   deallocate(Avww)
   deallocate(Avnn)
   deallocate(Avss)
   deallocate(Atp)
   deallocate(Ate)
   deallocate(Atw)
   deallocate(Atn)
   deallocate(Ats)
   deallocate(Atee)
   deallocate(Atww)
   deallocate(Atnn)
   deallocate(Atss)

!********************************************************
   contains
!********************************************************
!
!
!
!********************************************************
    subroutine initialize()
!************************************
    implicit none
!************************************
         do i=1,Nx
         do j=1,Ny

    Ut(i,j,1) = 0.0d0
    Ut(i,j,2) = 0.0d0
    Ut(i,j,3) = 0.0d0
    Ut(i,j,4) = 0.0d0
    Ut(i,j,5) = 0.0d0
   press(i,j) = 0.0d0
   pcorr(i,j) = 0.0d0
  forc(i,j,1) = 0.0d0
  forc(i,j,2) = 0.0d0

         end do
         end do
!********************************
!***at the top-layer


!************************************
    end subroutine initialize
!********************************************************
!
!
!
!********************************************************
     subroutine convective()
!************************************
    implicit none
!************************************
!************************************

       do i=2,Nx-1
       do j=2,Ny-1
       do jcnt=1,3

!************************************
   xintPLUS = ( x(i+1)+x(i) )/2.0d0
  xintMINUS = ( x(i)+x(i-1) )/2.0d0

   yintPLUS = ( y(j+1)+y(j) )/2.0d0
  yintMINUS = ( y(j)+y(j-1) )/2.0d0
!************************************
     dx_int = xintPLUS-xintMINUS
     dy_int = yintPLUS-yintMINUS
!************************************
    xcell_Pe = dabs(Re*Ut(i,j,1)*dx_int)
    ycell_Pe = dabs(Re*Ut(i,j,2)*dy_int)
!*****************************
!***along x-direction
!*****************************
   if(i.le.2.or.i.ge.(Nx-1))then       !*****###*****!
 x_conv=Ut(i,j,1)*( ( etaE_gradP(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt)) )+( etaW_gradP(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) ) )
   elseif(xcell_Pe.ge.2.0d0.and.Ut(i,j,1).gt.0.0d0)then
 x_conv = Ut(i,j,1)*( ( etaEE_upw_pos(i)*(Ut(i+2,j,jcnt)-Ut(i,j,jcnt)) )+( etaE_upw_pos(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt)) )+( etaW_upw_pos(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) )+( etaWW_upw_pos(i)*(Ut(i-2,j,jcnt)-Ut(i,j,jcnt)) ) )
   elseif(xcell_Pe.ge.2.0d0.and.Ut(i,j,1).lt.0.0d0)then
   x_conv = Ut(i,j,1)*( ( etaEE_upw_neg(i)*(Ut(i+2,j,jcnt)-Ut(i,j,jcnt)) )+( etaE_upw_neg(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt)) )+( etaW_upw_neg(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) )+( etaWW_upw_neg(i)*(Ut(i-2,j,jcnt)-Ut(i,j,jcnt)) ) )
   else
   x_conv = Ut(i,j,1)*( ( etaEE_1cen(i)*(Ut(i+2,j,jcnt)-Ut(i,j,jcnt)) )+( etaE_1cen(i)*(Ut(i+1,j,jcnt)-Ut(i,j,jcnt)) )+( etaW_1cen(i)*(Ut(i-1,j,jcnt)-Ut(i,j,jcnt)) )+( etaWW_1cen(i)*(Ut(i-2,j,jcnt)-Ut(i,j,jcnt)) ) )
   end if
!*****************************
!***along y-direction
!*****************************
   if(j.le.2.or.j.ge.(Ny-1))then       !*****###*****!
   y_conv = Ut(i,j,2)*( ( etaN_gradP(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt)) )+( etaS_gradP(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) ) )
   elseif(ycell_Pe.ge.2.0d0.and.Ut(i,j,2).gt.0.0d0)then
   y_conv = Ut(i,j,2)*( ( etaNN_upw_pos(j)*(Ut(i,j+2,jcnt)-Ut(i,j,jcnt)) )+( etaN_upw_pos(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt)) )+( etaS_upw_pos(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) )+( etaSS_upw_pos(j)*(Ut(i,j-2,jcnt)-Ut(i,j,jcnt)) ) )
   elseif(ycell_Pe.ge.2.0d0.and.Ut(i,j,2).lt.0.0d0)then
   y_conv = Ut(i,j,2)*( ( etaNN_upw_neg(j)*(Ut(i,j+2,jcnt)-Ut(i,j,jcnt)) )+( etaN_upw_neg(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt)) )+( etaS_upw_neg(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) )+( etaSS_upw_neg(j)*(Ut(i,j-2,jcnt)-Ut(i,j,jcnt)) ) )
   else
   y_conv = Ut(i,j,2)*( ( etaNN_1cen(j)*(Ut(i,j+2,jcnt)-Ut(i,j,jcnt)) )+( etaN_1cen(j)*(Ut(i,j+1,jcnt)-Ut(i,j,jcnt)) )+( etaS_1cen(j)*(Ut(i,j-1,jcnt)-Ut(i,j,jcnt)) )+( etaSS_1cen(j)*(Ut(i,j-2,jcnt)-Ut(i,j,jcnt)) ) )
   end if
!*****************************
!***obtain the final-convective term

   conv(i,j,jcnt) = x_conv+y_conv

!*****************************
       end do
       end do
       end do
!*****************************

       do i=2,Nx-1
       do j=2,Ny-1
!*****************************

    gradP(1) = (etaE_gradP(i)*(press(i+1,j)-press(i,j)))+(etaW_gradP(i)*(press(i-1,j)-press(i,j)))
    gradP(2) = (etaN_gradP(j)*(press(i,j+1)-press(i,j)))+(etaS_gradP(j)*(press(i,j-1)-press(i,j)))

    Q(i,j,4) = Ut(i,j,1)+dt*(-conv(i,j,1))
    Q(i,j,1) = Q(i,j,4)+(dt*-gradP(1))

    Q(i,j,5) = Ut(i,j,2)+dt*( -conv(i,j,2)+(Ri*Ut(i,j,3)) )
    Q(i,j,2) = Q(i,j,5)+(dt*-gradP(2) )

    Q(i,j,3) = Ut(i,j,3)+dt*(-conv(i,j,3))
!*****************************
       end do
       end do
!************************************
    end subroutine convective
!********************************************************
!
!
!
!********************************************************
     subroutine solve_transport()
!************************************
    implicit none
!************************************
 integer :: iter
 double precision :: sumav,sumnor,tmp,res
!************************************
!    iter=0
    sumav=0.0d0
!************************************
  do jcnt=1,5
!************************************
    do iter=1,maxiter
!************************************
    tmp=0.0d0

   do i=2,Nx-1
   do j=2,Ny-1
!*****************************
   if(jcnt.eq.3)then
   Ae = Ate(i,j)
   Aw = Atw(i,j)
   An = Atn(i,j)
   As = Ats(i,j)
   Ap = Atp(i,j)
  Aee = Atee(i,j)
  Aww = Atww(i,j)
  Ann = Atnn(i,j)
  Ass = Atss(i,j)
   else
   Ae = Ave(i,j)
   Aw = Avw(i,j)
   An = Avn(i,j)
   As = Avs(i,j)
   Ap = Avp(i,j)
  Aee = Avee(i,j)
  Aww = Avww(i,j)
  Ann = Avnn(i,j)
  Ass = Avss(i,j)
   end if
!*****************************
!   if(i.ge.4.and.j.ge.4.and.i.le.(Nx-3).and.j.le.(Ny-3))then       !*****###*****!
!*****************************
!   Ut(i,j,jcnt)=(Q(i,j,jcnt)-( (Aee*Ut(i+2,j,jcnt))+(Aww*Ut(i-2,j,jcnt))+(Ann*Ut(i,j+2,jcnt))+(Ass*Ut(i,j-2,jcnt))+(Ae*Ut(i+1,j,jcnt))+(Aw*Ut(i-1,j,jcnt))+(An*Ut(i,j-1,jcnt))+(As*Ut(i,j+1,jcnt)) ))/Ap
!*****************************
!   else                                 !*****###*****!
!*****************************
   Ut(i,j,jcnt)=( Q(i,j,jcnt)-( (Ae*Ut(i+1,j,jcnt))+(Aw*Ut(i-1,j,jcnt))+(An*Ut(i,j+1,jcnt))+(As*Ut(i,j-1,jcnt)) ) )/Ap
!*****************************
!   end if                               !*****###*****!
!*****************************
   end do
   end do

   do i=2,Nx-1
   do j=2,Ny-1
!*****************************
   if(jcnt.eq.3)then
   Ae = Ate(i,j)
   Aw = Atw(i,j)
   An = Atn(i,j)
   As = Ats(i,j)
   Ap = Atp(i,j)
  Aee = Atee(i,j)
  Aww = Atww(i,j)
  Ann = Atnn(i,j)
  Ass = Atss(i,j)
   else
   Ae = Ave(i,j)
   Aw = Avw(i,j)
   An = Avn(i,j)
   As = Avs(i,j)
   Ap = Avp(i,j)
  Aee = Avee(i,j)
  Aww = Avww(i,j)
  Ann = Avnn(i,j)
  Ass = Avss(i,j)
   end if
!*****************************
!*****************************
!   if(i.ge.4.and.j.ge.4.and.i.le.(Nx-3).and.j.le.(Ny-3))then       !*****###*****!
!*****************************
! res=Ut(i,j,jcnt)-(Q(i,j,jcnt)-( (Aee*Ut(i+2,j,jcnt))+(Aww*Ut(i-2,j,jcnt))+(Ann*Ut(i,j+2,jcnt))+(Ass*Ut(i,j-2,jcnt))+(Ae*Ut(i+1,j,jcnt))+(Aw*Ut(i-1,j,jcnt))+(An*Ut(i,j-1,jcnt))+(As*Ut(i,j+1,jcnt)) ))/Ap
!*****************************
!   else                                 !*****###*****!
!*****************************
   res = Ut(i,j,jcnt)-((Q(i,j,jcnt)-( (Ae*Ut(i+1,j,jcnt))+(Aw*Ut(i-1,j,jcnt))+(An*Ut(i,j+1,jcnt))+(As*Ut(i,j-1,jcnt)) ))/Ap)
!*****************************
!   end if                               !*****###*****!
!*****************************
    tmp = tmp+dabs(res)
! print*, iter,sumav
   end do
   end do

  if(tmp.eq.0) goto 111
  if(iter.eq.1) sumnor=tmp
!  if(sumnor.ne.0.0d0) sumav=tmp/sumnor
  sumav=tmp/sumnor

  if(sumav.le.errtol)goto 111
!************************************
   end do
!************************************
      print *,sumav,iter,jcnt
111    continue
  end do
!************************************
    end subroutine solve_transport
!********************************************************
!
!
!
!********************************************************
     subroutine rhei_chow()
!************************************
    implicit none
!************************************
    double precision, dimension(2) :: divVELstar
    double precision :: Ustar_xint_for,Ustar_xint_bck,Vstar_yint_for,Vstar_yint_bck
!************************************
!*****now we have got up,vp,wp and provisional u,v,w********************
!****************************************************************************************
!now applying Rhie-Chow Momentum Interpolation to avoid spurious oscillations in pressure
!****************************************************************************************
!*****boundary conditions for the w/o pressure field velocity (up,vp)*****************
         Ut(1:Nx,1,4) = 0.0d0    ! predicted up is zero at the bottom wall.
         Ut(1:Nx,1,5) = 0.0d0    ! predicted vp is zero at the bottom wall.

        Ut(1:Nx,Ny,4) = 1.0d0    ! predicted up is unity at the top wall.
        Ut(1:Nx,Ny,5) = 0.0d0    ! predicted vp is zero at the top wall.

        Ut(1,2:Ny-1,4)= 0.0d0    ! predicted up is zero at the back wall.
        Ut(1,2:Ny-1,5)= 0.0d0    ! predicted vp is zero at the back wall.

       Ut(Nx,2:Ny-1,4)= 0.0d0    ! predicted up is zero at the front wall.
       Ut(Nx,2:Ny-1,5)= 0.0d0    ! predicted vp is zero at the front wall.
!****************************************************************************************
!***********************
         do i=2,Nx-1
         do j=2,Ny-1


   Ustar_xint_for = ( 0.5d0*(Ut(i+1,j,4)+Ut(i,j,4)) )-( dt*(press(i+1,j)-press(i,j))/de(i) )
   Ustar_xint_bck = ( 0.5d0*(Ut(i,j,4)+Ut(i-1,j,4)) )-( dt*(press(i,j)-press(i-1,j))/dw(i) )


   Vstar_yint_for = ( 0.5d0*(Ut(i,j+1,5)+Ut(i,j,5)) )-( dt*(press(i,j+1)-press(i,j))/dn(j) )
   Vstar_yint_bck = ( 0.5d0*(Ut(i,j,5)+Ut(i,j-1,5)) )-( dt*(press(i,j)-press(i,j-1))/ds(j) )

    divVELstar(1) =  2.0d0*(Ustar_xint_for-Ustar_xint_bck)/(de(i)+dw(i))
    divVELstar(2) =  2.0d0*(Vstar_yint_for-Vstar_yint_bck)/(dn(j)+ds(j))

          Qp(i,j) = ( divVELstar(1)+divVELstar(2) )/dt
!****************************************************************************************
         end do
         end do
!************************************
    end subroutine rhei_chow
!********************************************************
!
!
!
!********************************************************
     subroutine pressure_poisson()
!************************************
    implicit none
!************************************
    double precision :: sumav,tmp,sumnor,res
    integer :: iter
!************************************
!*****initial condition for correction pressure*****************
         pcorr(1:Nx,1) = 0.0d0    ! pcorr is zero at the bottom wall.
         pcorr(1:Nx,1) = 0.0d0    ! pcorr is zero at the bottom wall.

        pcorr(1:Nx,Ny) = 0.0d0    ! pcorr is zero at the top wall.
        pcorr(1:Nx,Ny) = 0.0d0    ! pcorr is zero at the top wall.

        pcorr(1,2:Ny-1)= 0.0d0    ! pcorr is zero at the back wall.
        pcorr(1,2:Ny-1)= 0.0d0    ! pcorr is zero at the back wall.

       pcorr(Nx,2:Ny-1)= 0.0d0    ! pcorr is zero at the front wall.
       pcorr(Nx,2:Ny-1)= 0.0d0    ! pcorr is zero at the front wall.
!************************************
   iter=0
   sumav=0.0d0
!************************************
   do iter = 1,maxiter
    tmp=0.0d0

   do i=2,Nx-1
   do j=2,Ny-1
    res = pcorr(i,j)-( (Qp(i,j)-( (Ape(i,j)*pcorr(i+1,j))+(Apw(i,j)*pcorr(i-1,j))+(Apn(i,j)*pcorr(i,j+1))+(Aps(i,j)*pcorr(i,j-1)) ))/App(i,j) )
    tmp = tmp+dabs(res)
   end do
   end do

  if(tmp.eq.0)then
  goto 191
  end if
  if(iter.eq.1)then
  sumnor=tmp
  end if
!  if(sumnor.ne.0.0d0) sumav=tmp/sumnor
  sumav=tmp/sumnor

   do i=2,Nx-1
   do j=2,Ny-1
   pcorr(i,j)=(Qp(i,j)-( (Ape(i,j)*pcorr(i+1,j))+(Apw(i,j)*pcorr(i-1,j))+(Apn(i,j)*pcorr(i,j+1))+(Aps(i,j)*pcorr(i,j-1)) ))/App(i,j)
   end do
   end do

  if(sumav.lt.errtol)then
  goto 191
  end if

       end do

 print*, 'P-P',iter,sumav
191    continue
!****************************************
           do j=2,Ny-1
           pcorr(1,j) = ( (etaE*pcorr(2,j))+(etaEE*pcorr(3,j)) )/(etaE+etaEE)
           pcorr(Nx,j) = ( (etaW*pcorr(Nx-1,j))+(etaWW*pcorr(Nx-2,j)) )/(etaW+etaWW)
           end do

           do i=1,Nx
           pcorr(i,1) = ( (etaN*pcorr(i,2))+(etaNN*pcorr(i,3)) )/(etaN+etaNN)
           pcorr(i,Ny) = ( (etaS*pcorr(i,Ny-1))+(etaSS*pcorr(i,Ny-2)) )/(etaS+etaSS)
           end do
!************************************
    end subroutine pressure_poisson
!************************************
!
!
!
!********************************************************
     subroutine update_primitives()
!************************************
    implicit none
!************************************
    double precision :: ucorr,vcorr
    double precision, dimension(2) :: gradPcorr
!************************************
       do i=2,Nx-1
       do j=2,Ny-1
!*****************************
      press(i,j) = press(i,j)+pcorr(i,j)

    gradPcorr(1) = (etaE_gradP(i)*(pcorr(i+1,j)-pcorr(i,j)))+(etaW_gradP(i)*(pcorr(i-1,j)-pcorr(i,j)))
    gradPcorr(2) = (etaN_gradP(j)*(pcorr(i,j+1)-pcorr(i,j)))+(etaS_gradP(j)*(pcorr(i,j-1)-pcorr(i,j)))

           ucorr = -dt*gradPcorr(1)
           vcorr = -dt*gradPcorr(2)

       Ut(i,j,1) = Ut(i,j,1)+ucorr
       Ut(i,j,2) = Ut(i,j,2)+vcorr
!*****************************
         end do
         end do
!************************************
    end subroutine update_primitives
!************************************
!
!
!
!************************************
     subroutine boundary_condition()
!************************************
    implicit none
!************************************
    double precision :: del2U_delx2,del2U_dely2
    double precision :: del2V_delx2,del2V_dely2
    double precision :: sumtmp
    double precision :: conv_v_1,conv_v_2
!************************************
!*****boundary conditions for the velocity (u,v)*****************
         Ut(1:Nx,1,1) = 0.0d0    ! u is zero at the bottom wall.
         Ut(1:Nx,1,2) = 0.0d0    ! v is zero at the bottom wall.
         Ut(1:Nx,1,3) = 0.0d0    ! T is zero at the bottom wall.

        Ut(1:Nx,Ny,1) = 1.0d0    ! u is unity at the top wall.
        Ut(1:Nx,Ny,2) = 0.0d0    ! v is zero at the top wall.

        Ut(1,2:Ny-1,1)= 0.0d0                  ! u is zero at the back wall.
        Ut(1,2:Ny-1,2)= 0.0d0                  ! v is zero at the back wall.

        do j=2,Ny-1
        Ut(1,j,3) = 1.0d0 !dsin(kappa*Pi*y(j))        ! T is sinusoidal at the back wall.
        end do

       Ut(Nx,2:Ny-1,1)= 0.0d0                  ! u is zero at the front wall.
       Ut(Nx,2:Ny-1,2)= 0.0d0                  ! v is zero at the front wall.

       do j=2,Ny-1
       Ut(Nx,j,3) = 1.0d0 !dsin(kappa*Pi*y(j))        ! T is sinusoidal at the front wall.
       end do
!****************************************************************************************
!apply adiabatic boundary condition for temperature on the moving fluid-layer.

     do i=1,Nx
     Ut(i,Ny,3) = ( (etaS*Ut(i,Ny-1,3))+(etaSS*Ut(i,Ny-2,3)) )/(etaS+etaSS)
     end do

!***at i=1 and for all j=2,Ny-1

     do j=2,Ny-1

     del2U_delx2 = ( etaE_bnd2*Ut(2,j,1) )+( etaEE_bnd2*Ut(3,j,1) )
     del2U_dely2 = 0.0d0
     sumtmp = (del2U_delx2+del2U_dely2)/Re
     press(1,j) = ( (etaE*press(2,j))+(etaEE*press(3,j))-sumtmp )/(etaE+etaEE)

     end do

!***at i=Nx and for all j=2,Ny-1

     do j=2,Ny-1

     del2U_delx2 = ( etaW_bnd2*Ut(Nx-1,j,1) )+( etaWW_bnd2*Ut(Nx-2,j,1) )
     del2U_dely2 = 0.0d0
     sumtmp = (del2U_delx2+del2U_dely2)/Re
     press(Nx,j) = ( (etaW*press(Nx-1,j))+(etaWW*press(Nx-2,j))-sumtmp )/(etaW+etaWW)

     end do

!***at j=1 and for all i=1,Nx

     do i=1,Nx

     del2V_dely2 = ( etaN_bnd2*Ut(i,2,2) )+( etaNN_bnd2*Ut(i,3,2) )
     del2V_delx2 = 0.0d0
     sumtmp = (del2V_delx2+del2V_dely2)/Re
     press(i,1) = ( (etaN*press(i,2))+(etaNN*press(i,3))-sumtmp )/(etaN+etaNN)

     end do

!***at j=Ny and for all i=1,Nx

     do i=1,Nx

     del2V_dely2 = ( etaS_bnd2*Ut(i,Ny-1,2) )+( etaSS_bnd2*Ut(i,Ny-2,2) )
     del2V_delx2 = 0.0d0
     sumtmp = (del2V_delx2+del2V_dely2)/Re
     press(i,Ny) = ( (etaS*press(i,Ny-1))+(etaSS*press(i,Ny-2))-sumtmp )/(etaS+etaSS)

     end do
!************************************
    end subroutine boundary_condition
!************************************
!
!
!
!************************************
     subroutine divergence()
!************************************
    implicit none
!************************************
    double precision :: div_x,div_y,summ,div
!************************************

       summ = 0.0d0

       do i=2,Nx-1
       do j=2,Ny-1
!*****************************
     div_x = (etaE_gradP(i)*(Ut(i+1,j,1)-Ut(i,j,1)))+(etaW_gradP(i)*(Ut(i-1,j,1)-Ut(i,j,1)))
     div_y = (etaN_gradP(j)*(Ut(i,j+1,2)-Ut(i,j,2)))+(etaS_gradP(j)*(Ut(i,j-1,2)-Ut(i,j,2)))

     summ = summ+div_x+div_y
!*****************************
         end do
         end do

      div = dsqrt( dabs(summ)/dble((Nx-2)*(Ny-2)) )

      print*, nstep,time,div
!************************************
    end subroutine divergence
!************************************
!
!
!
!************************************
     subroutine savesnapshot()
!************************************
    implicit none
!************************************
!*****write the flow-field data in ASCII-Tecplot format directly*****
 write(filname,'(a,i7.7,a)')'2D',nstep/snapshot,'.tp'
 open(unit=7,file=filname,status='unknown')

 write(7,13) 'TITLE ="',time,nstep,'"'
 write(7,*) 'VARIABLES = "x","y","U","V","temp","Up","Vp","press","pcorr"'
 write(7,14) 'ZONE j=',Nx,',i=',Ny,',DATAPACKING ="POINT"'

    do i=1,Nx
    write(7,*)
           do j=1,Ny
!*******************
  write(7,59) x(i),y(j),Ut(i,j,1),Ut(i,j,2),Ut(i,j,3),Ut(i,j,4),Ut(i,j,5),press(i,j),pcorr(i,j)
           end do
    end do

   close(7)

13 format(A,F10.6,I8,A)
14 format(A,I3,A,I3,A)
59 format(2(1X,F10.7),7(1X,E14.7))
!************************************
    end subroutine savesnapshot
!************************************
!
!
!
!********************************************************
  end program incompressible
!********end of program**********************************

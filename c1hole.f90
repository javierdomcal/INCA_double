subroutine c1hole(nrad,sfalpha)
    use wfxinfo
    use geninfo
    implicit none
    !radial quadrature!!!!!!!!
    integer, intent(in) :: nrad
    double precision, intent(in) :: sfalpha
    !Local variables
    double precision, dimension(nrad) :: r, rad, wr
    !angular quadrature!!!!!!!
    double precision, dimension(1202) :: wa, lbx, lby, lbz, phi,theta
    double precision, dimension(1202) :: cosin, sinsin, cost
    integer :: nang
    double precision :: xs
    !grid points!!!!!!!!!!!!!!
    double precision :: px,py,pz
    double precision :: x,y,z
    !functions!!!!!!!!!!!!!!!!
    double precision :: rDM1, beckex, Density, rdm1_alf
    double precision, dimension(40) :: hc1
    double precision :: dens, xchole
    !!!!!!!!!!!!!!!!!!!!!!!!!/
    double precision :: rdmv1, rho2, rho2s, rdmv2 !vector 1rdm and pair dens
    integer :: i,j,k,s
    double precision :: kk, xcint, xcint2, radint, densint, rdm1int, radint2
    double precision :: b,alf
    character*40 :: outputname
    integer :: i1, npoints
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!obtain radial nodes and weights for integrations!!!!!!!!!!!

    call sub_GauLeg(-1.d0,1.d0,rad,wr,nrad)
    do i=1,nrad
         r(i)=(1.d0+rad(i))/(1.d0-rad(i))*sfalpha
         wr(i)=(2.d0*sfalpha/((1.d0-rad(i))**2.d0))*wr(i)
    end do

    !obtain angular nodes (angles) and weights
    lbx=0.d0
    lby=0.d0 
    lbz=0.d0
    nang=1202
    !call LD0006(lbx,lby,lbz,Wa,nAng)  !6 grid points
    !call LD0110(lbx,lby,lbz,wa,nang)
    !call LD0590(lbx,lby,lbz,wa,nAng) !590
    call LD1202(lbx,lby,lbz,wa,nAng)

    do i=1,nang
         theta(i)=acos(lbz(i))
         if (sin(theta(i)).ne.0.d0) then
              xs=lbx(i)/sin(theta(i))
              if (xs.gt.1.d0) xs=1.d0
              if (xs.lt.-1.d0) xs=-1.d0
              phi(i)=acos(xs)
              if (lby(i).lt.0.d0) phi(i)=-phi(i)
         else
              phi(i)=0.d0
         end if   
         wa(i)=wa(i)*4.d0*pi 
    end do

    open(unit=3, file="angular_grid")
    do i=1,nang
           x=cos(phi(i))*sin(theta(i))
           y=sin(phi(i))*sin(theta(i))
           z=cos(theta(i))
           write(3,*) x,y,z
           cosin(i)=x
           sinsin(i)=y
           cost(i)=z
    end do
    close(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!Sanity checks!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Test 1: integrate the density and 1rdm diagonal (must be = nelec)
    radint=0.d0
    radint2=0.d0
    !do i=1,nrad
    !   densint=0.d0
    !   rdm1int=0.d0
    !   do k=1,nang
    !       x=r(i)*cos(phi(k))*sin(theta(k))
    !       y=r(i)*sin(phi(k))*sin(theta(k))
    !       z=r(i)*cos(theta(k))
    !       densint=densint+wa(k)*Density(x,y,z)
    !       rdm1int=rdm1int+wa(k)*rDM1(x,y,z,x,y,z)
    !   end do    
    !   radint=radint+wr(i)*densint*r(i)**2.d0   
    !   radint2=radint2+wr(i)*rdm1int*r(i)**2.d0
   ! end do  
   ! write(*,*) "----------Test 1-------------" 
   ! write(*,*) "integral of density=", radint,"=?", nelec
   ! write(*,*) "integral of diagonal 1rdm=", radint2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!end of test 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     open(unit=3, file='refpoints.inp')
     read(3,*) npoints
     do i1=1,npoints      
     read(3,*) x,y,z, outputname
    !!!!!!Test 2: Compute exchange hole at reference point:!!!!!!!!!!!!!
    !!!!!!!!!!-With Becke model
    !!!!!!!!! -With 1rdm^2/density
      open(unit=4, file=outputname)
      write(*,*) "Reference point"
      write(*,*) x,y,z
      write(*,*) "------------------"
      call becke1(x,y,z,b,alf)
      write(*,*) "b=", b
      write(*,*) "a=", alf
      kk=0.d0 !interelectronic distance
      write(*,*) "Loop over distances"
      do i=1,1000
          kk=kk+0.005d0
          xchole=beckex(b,alf,kk)   !exchange hole model
          rdmv1=0.d0
          !do k=1,nang !angular integral of the squared 1rdm
          !     px=x+kk*cosin(k)  !cos(phi(k))*sin(theta(k))
          !     py=y+kk*sinsin(k) !sin(phi(k))*sin(theta(k))
          !     pz=z+kk*cost(k)   !cos(theta(k))
          !     rdmv1=rdmv1+wa(k)*(rDM1_alf(x,y,z,px,py,pz)*rDM1_alf(px,py,pz,x,y,z))
          !end do
          write(4,*) kk, xchole!, rdmv1/((Dens)*4.d0*pi) !write both models
      end do
      STOP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!end of test 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Test 3: Integrate xchole, becke model and 1rdm^2/dens!!!!!!!!!!!!!!!!
      write(*,*) "End of test 2"
      xcint=0.d0
      xcint2=0.d0
      do i=1,nrad
          rdmv1=0.d0
          rdmv2=0.d0
          do k=1,nang
             px=x+r(i)*cos(phi(k))*sin(theta(k))
             py=y+r(i)*sin(phi(k))*sin(theta(k))
             pz=z+r(i)*cos(theta(k))
             rdmv1=rdmv1+wa(k)*(rDM1(x,y,z,px,py,pz)*rDM1(px,py,pz,x,y,z))
          end do
          xcint2=xcint2+wr(i)*r(i)**2.d0*(rdmv1)
          xchole=beckex(b,alf,r(i))
          xcint=xcint+wr(i)*xchole*(r(i)**2.d0)*4.d0*pi
      end do
      write(*,*) "----------Test 3-------------"
      write(*,*) "Integral of exchange hole="
      write(*,*) "Becke Model=", xcint
      if (uhf) then 
            write(*,*) "1rdm^2/dens=", xcint2/(dens)
      else        
            write(*,*) "1rdm^2/dens", xcint2/(dens/2.d0)
      end if        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    close(4)  
    !see the plot of the exchange hole vs the distance at (x,y,z) point.
    end do
    close(3)
    write(*,*) "-----Test xchole done-----------"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Compute C1 part of Coulomb Hole!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(unit=4, file="c1hole.out")
    kk=0.0001d0   !interelectronic distance
    do s=1,40 !for 20 interelectronic distances (EXAMPLE)
         hc1(s)=0.d0     
         do i=1,nrad       !loop for radial nodes
              rho2s=0.d0   !spherically averaged pair density at r dist and s dist.
              do j=1,nang  !loop for angular nodes (to average r)
                   x=r(i)*cos(phi(j))*sin(theta(j))  !compute x,y,z points (\Vec{r})
                   y=r(i)*sin(phi(j))*sin(theta(j))
                   z=r(i)*cos(theta(j))
                   rdmv1=0.d0 !spherically averaged density at s distance at r point. 
                   do k=1,nang
                        !compute interelectronic separation vector (px,py,pz) (to average s)
                        px=x+kk*cos(phi(k))*sin(theta(k))
                        py=y+kk*sin(phi(k))*sin(theta(k))
                        pz=z+kk*cos(theta(k))                                           
                        !spherical average of s
                        rdmv1=rdmv1+wa(k)*rDM1(x,y,z,px,py,pz)*rDM1(px,py,pz,x,y,z)
                   end do
                   call becke1(x,y,z,b,alf)
                   xchole=beckex(b,alf,kk)*4.d0*pi
                   rho2=rdmv1-xchole*(Density(x,y,z)) !c1 part of Coulomb hole: hc(\Vec{r},s)
                   rho2s=rho2s+wa(j)*rho2   !spherical average of \Vec{r}
              end do
              hc1(s)=hc1(s)+wr(i)*rho2s*r(i)**2.d0 !radial average of r
         end do
         write(4,*) kk, hc1(s)
         kk=kk+0.02d0
    end do 
end subroutine c1hole

subroutine becke1(x,y,z,b,alf) !give a point and returns a and b from eq 17 of Becke-Rousell
  use fractions
  use geninfo
  double precision, intent(in) :: x,y,z
  double precision, intent(out) :: b,alf
  !!!!!local variables!!!!!!!!!!!!!!!!!!
  double precision :: xx ! x in Becke-Rousell paper. 
  double precision :: dens, E_kin, Grad, Lapl, newtr
  double precision, parameter :: trsh=1d-16
  write(*,*) "In Becke 1"
  call dens_ops(x,y,z,Dens,E_kin,Grad,Lapl) !subroutine to compute densities (one spin)
  write(*,*) "Dens, E_kin, Grad, Lapl" !this quantities are for one spin
  write(*,*) Dens, E_kin, Grad, Lapl
  if (abs(Dens).ge.trsh) then
        xx=newtr(Dens,E_kin,Grad,Lapl)
        write(*,*) "xx=", xx
        if (xx.lt.0.d0) then
                write(*,*) "Danger, negative root"
                xx=abs(xx)
        end if
        write(*,*) "xx=", xx
        b=(xx**(3.d0)*exp(-xx)*0.125d0*((pi*Dens)**(-1.d0)))**(hr)
        alf=xx*(b**(-1.d0))
  else
        b=0.d0
        alf=0.d0
  end if
end subroutine becke1        



function beckex(b,alf,s) !Becke-Roussel Exchange Hole, eqn 16
    use fractions
    use geninfo
    implicit none    
    double precision :: beckex
    double precision, intent(in) :: b,alf,s  !r and s are distances
    !!!local variables!!!!!!!!!!!!!!!!!!
    double precision :: p1, p2, p3
        p1=0.0625d0*alf*(pi*b*s)**(-1.d0)
        p2=(alf*abs(b-s)+1.d0)*exp(-alf*abs(b-s))
        p3=(alf*abs(b+s)+1.d0)*exp(-alf*abs(b+s))
        beckex=p1*(p2-p3)
end function beckex 

function newtr(Dens,E_kin,Grad,Lapl) !performs Newton-Raphson to find x
   use fractions  
   use geninfo   
   implicit none
   double precision :: newtr
   double precision, intent(in) :: Dens, E_kin, Grad, Lapl
   !!!local variables!!!!!
   double precision :: dif
   double precision :: x0
   double precision :: qval
   double precision :: f, fp, fpr, fr,k
   double precision :: q
   integer, parameter :: maxit=1000
   double precision, parameter :: trsh=1d-8
   integer :: i
   double precision :: sf !scaling factor (step size)
   qval=q(Dens,E_kin,Grad,Lapl) !compute eq 20b
   write(*,*) "Q=", qval
   qval=qval
   sf=0.5d0

   k=bih*pi**(bih)*Dens**(boh)*qval**(-1.d0) !right hand part of eq 21
   if (qval.gt.0.d0) then
           x0=3.d0 !see plot in mathematica (we want to stard from the right)
           do i=1,maxit
               fpr=fp(x0)
               fr=f(x0,k)
               if (fr/fpr.ge.1.d0) then
                    sf=10*fpr/fr
               end if        
               newtr=x0-sf*(fr*(fpr**(-1.d0)))
               fr=f(newtr,k)
               if (newtr.lt.2.d0) then
                   newtr=2.5d0
               end if        
               write(*,*) "xx=", newtr
               dif=abs(newtr-x0)
               if (dif.lt.trsh) then
                  write(*,*) "Converged after", i, "iterations"
                  exit
               end if
               if (i.eq.maxit) then
                       write(*,*) "Newton-Raphson1 not converged"
                       STOP
               end if        
               x0=newtr
           end do
   else
           x0=-1.d0
           do i=1,maxit
              fpr=fp(x0)
              fr=f(x0,k)
              newtr=x0-sf*fr*(fpr**(-1.d0))
              write(*,*) "xx=", newtr
              if (newtr.gt.2.d0) then
                    newtr=1.5d0
              end if
              dif=abs(newtr-x0)
              if (dif.lt.trsh) then
                  write(*,*) "Converger after", i, "iterations"
                  exit
              end if
              if (i.eq.maxit) then
                   write(*,*) "Newton-Raphson2 not converged"
                   STOP
              end if
              x0=newtr
           end do
   end if         
end function

function f(x,k) !equation 21 of Becke-Roussel
    use fractions 
    use geninfo   
    implicit none  
    double precision, parameter :: trsh=1d-8
    double precision, intent(in) :: x,k 
    double precision :: f1,f 
    f1=x*exp(-bih*x)*(x-2.d0)**(-1.d0)                             !nomes cal calcular q una vegada
    f=f1-k
end function

function fp(x) !derivative of equation 22
    use fractions
    use geninfo !for pi variable
    implicit none    
    double precision :: fp
    double precision, intent(in) :: x
    fp=-bih*exp(-bih*x)*(x**2.d0-2.d0*x+3.d0)*(x-2.d0)**(-2.d0)

end function

function q(Dens,E_kin,Grad,Lapl) !equation 20b
    use fractions
    implicit none    
    double precision, intent(in) :: Dens,E_kin,Grad,Lapl
    double precision :: D
    double precision :: q
    double precision, parameter :: gmm=0.8d0
    q=bs*(Lapl-2.d0*gmm*D(Dens,E_kin,Grad))
end function

function D(Dens,E_kin,Grad) !equation 13b
    implicit none
    double precision,intent(in) :: Dens,E_kin,Grad
    double precision :: D
    double precision, parameter :: trsh=1d-15   
    if (abs(Dens).lt.trsh) then
           D=E_kin
    else       
           D=E_kin-0.25d0*(Grad)*Dens**(-1.d0) !Grad-->scalar product of gradient
    end if       
    write(*,*) "D=", D
end function



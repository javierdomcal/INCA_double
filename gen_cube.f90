!*********************************************************************
!!Creates a cubefile (.cube extension) which can be opened on vmd    !
!!and see the 3D representation of the magnitude we want i.e.        !
!!density, gradient of the density, laplacian, MOs, AOs, etc         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cubefile(a,cubename) 
!creates a cubefile with a magnitude we choose f(a)
use geninfo
use cubeinfo
implicit none
character*40 :: cubename
character*5 :: cubepart
integer, intent(in) :: a !a=1 primitive, a=2 AO, a=3 MO, a=5 Density, a=6 Laplacian
integer :: i, j, aa 
double precision, parameter :: zero=0.d0

aa=0

if ((uhf).and.((a.eq.3).or.(a.eq.5))) then    !if we have an unrestricted wavefunction
       do i=1,3                              !do a loop for alpha, beta and total wf.
                aa=aa+1                       !3 cubefiles will be generated
                if (i.eq.1) cubepart="alpha"  !insert description in the name of the .cube
                if (i.eq.2) cubepart="beta_"  !so that we know what we have represented
                if (i.eq.3) cubepart="all__"
                open(unit=2,file=cubepart//cubename) 
                write(2,*) "CUBE FILE"
                write(2,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
                write(2,*) natoms, center(1), center(2), center(3)
                write(2,*) np(1), step(1), zero, zero
                write(2,*) np(2), zero, step(2), zero
                write(2,*) np(3), zero, zero, step(3)
                do j=1,natoms
                     write(2,*) an(j), chrg(j), cartes(j,1), cartes(j,2), cartes(j,3)
                end do
                call fwrite(a,aa)  !writes the values of the function in the cubepoints
                close(2)
       end do
else  
       open(unit=2,file=cubename)
       write(2,*) "CUBE FILE"
       write(2,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
       write(2,*) natoms, center(1), center(2), center(3)
       write(2,*) np(1), step(1), zero, zero
       write(2,*) np(2), zero, step(2), zero
       write(2,*) np(3), zero, zero, step(3)
       do i=1,natoms
            write(2,*) an(i), chrg(i), cartes(i,1), cartes(i,2), cartes(i,3)
        end do
        call fwrite(a,aa) !writes the value of a function in each point
        close(2)
end if
end subroutine cubefile


subroutine fwrite(a,aa)    !evaluates a given function in the specified points.
use geninfo
use cubeinfo
implicit none
integer, intent(in) :: a, aa
integer :: i,j,k
double precision :: Prim, AO, MoOr, MO_a, MO_b, Density, Dens_a, Dens_b, Lapl !Laplacian
double precision :: Dens, MO_t
double precision, dimension(3) :: xm !minimum point when np odd or even
double precision :: x,y,z
logical :: num

num=.true. !interruptor manual per numerical or analitical !change this to choose in the input
x=0.
y=0.
z=0.
xm=0.d0
do i=1,3
    if (mod(np(i),2).eq.0)then !odd number
         xm(i)=center(i)-(step(i)/2)-((np(i)-2)/2)*step(i)
    else !even number  !Calculates the starting point acording to cubeinfo data
         xm(i)=center(i)-((np(i)-1)/2)*step(i)
    end if  
end do

do i=1,np(1)         !Depending on a calculates the point with a different function
     x=xm(1)+step(1)*(i-1)
     do j=1,np(2) 
          y=xm(2)+step(2)*(j-1)
          do k=1,np(3)
              z=xm(3)+step(3)*(k-1)
              if (a.eq.1) then                                           
                    write(2,40) Prim(x,y,z,npr)
              else if (a.eq.2) then
                    write(2,40) AO(x,y,z,cao)   
              else if (a.eq.3) then
                    if (uhf) then
                         if (aa.eq.1) write(2,40) MO_a(x,y,z,mo)
                         if (aa.eq.2) write(2,40) MO_b(x,y,z,mo)
                         if (aa.eq.3) then
                              MO_t=0.d0         
                              MO_t=MO_a(x,y,z,mo)+MO_b(x,y,z,mo)
                              write(2,40) MO_t
                         end if
                    else
                         write(2,40) MoOr(x,y,z,mo) 
                    end if 
              else if (a.eq.5) then
                    if (uhf) then
                         if (aa.eq.1) write(2,40) Dens_a(x,y,z)
                         if (aa.eq.2) write(2,40) Dens_b(x,y,z)
                         if (aa.eq.3) then 
                              Dens=0.d0
                              Dens=Dens_a(x,y,z)+Dens_b(x,y,z)
                              write(2,40) Dens
                         end if            
                    else
                         write(2,40) Density(x,y,z)    
                    end if 
              else if (a.eq.6) then !laplacian
                    if (num) then
                        write(2,40) Lapl(x,y,z) !numerically
                    else
                    !write(2,40) Laplacian(x,y,z) !analitically (does not work)
                        write(*,*) "Not available"       
                    end if                
              else   
                    write(*,*) "Error, no cubefile will be generated"
              end if       
          end do
     end do     
end do  

40 format(6(E16.6E3))

end subroutine fwrite

!*******************************************************************************
 function Prim(x,y,z,npr) 
   use geninfo
   implicit none
   double precision :: Prim
   double precision, intent(in) :: x, y, z
   integer, intent(in) :: npr
   integer :: i,j,k
   i=1
   j=2
   k=3
   Prim=0.d0
   Prim=(x-(cartes(Ra(npr),i)))**(TMN(npr,i))* &
             (y-(cartes(Ra(npr),j)))**(TMN(npr,j))* &
            (z-(cartes(Ra(npr),k)))**(TMN(npr,k))* &
          exp(-Alpha(npr)*((x-cartes(Ra(npr),i))**2+(y-(cartes(Ra(npr),j)))**2+ & 
           (z-(cartes(Ra(npr),k)))**2))
              
 end function
   
!**************************************************************************
function sd(x,y,z,npr) !second derivative of primitives analitically
use geninfo
implicit none
double precision :: sd, Prim
double precision, intent(in) :: x, y, z
integer, intent(in) :: npr
double precision :: laplx, laply, laplz
double precision :: t,m,n
double precision :: xa,ya,za

t=0.d0 
m=0.d0 
n=0.d0
xa=0.d0
ya=0.d0 
za=0.d0
laplx=0.d0 
laply=0.d0
laplz=0.d0
sd=0.d0

t=TMN(npr,1)
m=TMN(npr,2)
n=TMN(npr,3)
xa=cartes(Ra(npr),1)
ya=cartes(Ra(npr),2)
za=cartes(Ra(npr),3)

laplx=(Prim(x,y,z,npr)/(x**t))*(t*(t-1)*x**(t-2) &
+4*t*(x**(t-1))*Alpha(npr)*(xa-x)+((x**t)*Alpha(npr) &
*(4*Alpha(npr)*(xa-x)-2)))

laply=(Prim(x,y,z,npr)/(y**m))*(m*(m-1)*y**(m-2) &
+4*m*(y**(m-1))*Alpha(npr)*(ya-y)+((y**m)*Alpha(npr) &
*(4*Alpha(npr)*(ya-y)-2)))

laplz=(Prim(x,y,z,npr)/(z**n))*(n*(n-1)*z**(n-2) &
+4*n*(z**(n-1))*Alpha(npr)*(za-z)+((z**n)*Alpha(npr) &
*(4*Alpha(npr)*(za-z)-2)))

sd=laplx+laply+laplz

end function
    
      
!**************************************************************************
   function AO(x,y,z,cao)
   use geninfo
   use loginfo      !cao from 1 to nao, 1->1S, 2->2S, 3->2Px, 4->2Py, ...
   implicit none
   double precision :: AO
   double precision :: Prim
   double precision, intent(in) :: x, y, z
   integer, intent(in) :: cao !choosen ao in the input file
   integer :: smm !sum for all the primitives
   integer :: i
   integer :: sm, kk, k, kp
   smm=0
   AO=0.d0
   sm=0
   k=1
   if (cao.eq.1) then 
     k=1
   else if (cao.eq.2) then
     kp=1
     k=k+npao(1)
   else  
     kp=cao-1    !kp=previous AO 
     do i=1,kp  !set k starting primitive, Flg(k)
      k=k+npao(i)  !sum the primitives until the current AO
     end do
   end if
     
   kk=k+npao(cao)-1 !kk=last primitive of the current AO (cao)
    do i=k,kk        
      AO= AO + Flg(i) * Prim(x,y,z,i)    
    end do
 
  end function AO
   
!*************************************************************************************   
   function MoOr(x,y,z,mo)
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MoOr
   double precision :: Prim
   integer, intent(in) :: mo
   integer :: i    
   MoOr=0.d0
   do i=1,nprim
      MoOr= MoOr+ T(mo,i)*Prim(x,y,z,i)
   end do   
   end function
 !**************************************************************************
  function MO_a(x,y,z,mo) 
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x,y,z
   integer, intent(in) :: mo
   double precision :: MO_a
   double precision :: Prim
   integer :: i
   MO_a=0.d0
   do i=1,nprim
     MO_a= MO_a + T_a(mo,i) * Prim(x,y,z,i)
   end do
   end function
   
   function MO_b(x,y,z,mo)
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x,y,z
   integer, intent(in) :: mo
   double precision :: MO_b
   double precision :: Prim
   integer :: i
   MO_b=0.d0
   do i=1,nprim
     MO_b=MO_b+ T_b(mo,i)*Prim(x,y,z,i)
   end do
   end function   
    
 !**********************************************************************
   function Density(x,y,z)   !from wfx file
   use wfxinfo 
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MoOr
   double precision :: Density
   integer :: i
    Density=0.d0
    if (corr) then
      do i=1,noccmo
        Density= Density + Occ(i) * MoOr(x,y,z,i)**2  
      end do        
    else
     do i=1,noccmo 
      Density=Density+MoOr(x,y,z,i)**2
     end do   
    end if
    end function
!*******************************************************************     
   function Dens_a(x,y,z) 
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MO_a
   double precision :: Dens_a
   integer :: i
     Dens_a=0.d0
     do i=1,nalfae
       Dens_a=Dens_a+MO_a(x,y,z,i)**2
     end do
    end function  
 !****************************************************************
   function Dens_b(x,y,z) 
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MO_b
   double precision :: Dens_b
   integer :: i
     Dens_b=0.d0
     do i=1,nalfae
       Dens_b=Dens_b+MO_b(x,y,z,i)**2
     end do
    end function  
!*****************************************************************
    
function lapl(x,y,z) !computes Laplacian numerically
use wfxinfo
implicit none 
double precision :: lapl                           !this could be a function!!
double precision, intent(in) :: x,y,z
double precision :: spx, spy, spz !laplacian, square partial x, s.p.y, spz
double precision :: density,h

h=0.0000001d0

spx=(Density(x+h,y,z)+Density(x-h,y,z)-2.d0*Density(x,y,z))/(h**2)
spy=(Density(x,y+h,z)+Density(x,y-h,z)-2.d0*Density(x,y,z))/(h**2)
spz=(Density(x,y,z+h)+Density(x,y,z-h)-2.d0*Density(x,y,z))/(h**2)

lapl=spx+spy+spz

end function

!************************************************************************* 

!function Laplacian(x,y,z)  !laplacian analitically using fd and sd functions (first derivative and second derivative of primitives  !no funtziona
!use wfxinfo
!use geninfo
!implicit none
!double precision :: Laplacian, Laplx, Laply, Laplz, dp, d2p
!double precision :: prim
!double precision, intent(in) :: x, y, z  
!integer :: i,j,k

!Laplacian=0.d0                      !consider the UHF case!
!Laplx=0.d0
!Laply=0.d0
!Laplz=0.d0

!do i=1,noccmo    
!          !current molecular orbital value (also valid for UHF)
!   do j=1,nprim
!    do k=1,nprim    
!     Laplx=Laplx+2*(T(i,j)*T(i,k)*d2p(x,y,z,j,1)*Prim(x,y,z,k)+(T(i,j)*dp(x,y,z,j,1))**2)   
!     Laply=Laply+2*(T(i,j)*T(i,k)*d2p(x,y,z,j,2)*Prim(x,y,z,k)+(T(i,j)*dp(x,y,z,j,2))**2)   
!     Laplz=Laplz+2*(T(i,j)*T(i,k)*d2p(x,y,z,j,3)*Prim(x,y,z,k)+(T(i,j)*dp(x,y,z,j,3))**2)            
!     end do                
!    end do   
!end do
!Laplacian=laplx+laply+laplz
!end function

function dp(x,y,z,j,part)   !no funtziona
use geninfo
implicit none 
double precision, intent(in) :: x,y,z
double precision :: dp, prim
integer, intent(in) :: part !1->x, 2->2, 3->z
integer, intent(in) :: j !number of primitive
double precision :: ti, m, n, xa, ya, za

  ti=TMN(j,1)
  m=TMN(j,2)
  n=TMN(j,3)
  xa=cartes(Ra(j),1)
  ya=cartes(Ra(j),2)
  za=cartes(Ra(j),3)

 if (part.eq.1) then !first partial derivative of prim with respect to x
   dp=(Prim(x,y,z,j)*(ti*x**(ti-1))/(x**ti))+Prim(x,y,z,j)*2*Alpha(j)*(xa-x)
 else if (part.eq.2) then !y
   dp=(Prim(x,y,z,j)*(m*y**(m-1))/(x**m))+Prim(x,y,z,j)*2*Alpha(j)*(ya-y)
 else if (part.eq.3) then !z
   dp=(Prim(x,y,z,j)*(n*z**(n-1))/(x**n))+Prim(x,y,z,j)*2*Alpha(j)*(za-z)
 end if   
 
end function

function d2p(x,y,z,j,part)   !no funtziona
use geninfo
implicit none 
double precision, intent(in) :: x,y,z
double precision :: d2p, prim
integer, intent(in) :: part !1->x, 2->2, 3->z
integer, intent(in) :: j !number of primitive
double precision :: ti, m, n, xa, ya, za
  ti=TMN(j,1)
  m=TMN(j,2)
  n=TMN(j,3)
  xa=cartes(Ra(j),1)
  ya=cartes(Ra(j),2)
  za=cartes(Ra(j),3)

 if (part.eq.1) then !first partial derivative of prim with respect to x
   d2p=(Prim(x,y,z,j)/(x**ti))*(ti*(ti-1)*x**(ti-2) &
    +4*ti*(x**(ti-1))*Alpha(j)*(xa-x)+((x**ti)*Alpha(j) &
    *(4*Alpha(j)*(xa-x)-2)))
 else if (part.eq.2) then !y
   d2p=(Prim(x,y,z,j)/(y**m))*(m*(m-1)*y**(m-2) &
    +4*m*(y**(m-1))*Alpha(j)*(ya-y)+((y**m)*Alpha(j) &
    *(4*Alpha(j)*(ya-y)-2)))
 else if (part.eq.3) then !z
   d2p=(Prim(x,y,z,j)/(z**n))*(n*(n-1)*z**(n-2) &
    +4*n*(z**(n-1))*Alpha(j)*(za-z)+((z**n)*Alpha(j) &
    *(4*Alpha(j)*(za-z)-2)))
 end if   

end function
              
!*******************************************************************

subroutine gradient(x,y,z,gradx,grady,gradz) !computes the gradient numerically
                                             !now we won't use that but could be useful
 implicit none
 double precision, intent(in) :: x,y,z
 double precision, intent(out) :: gradx, grady, gradz
 double precision :: density, h

 h=0.0000000000001d0 !change in the derivative

 gradx=(Density(x+h,y,z)-Density(x,y,z))/h
 grady=(Density(x,y+h,z)-Density(x,y,z))/h
 gradz=(Density(x,y,z+h)-Density(x,y,z))/h
  
end subroutine
    
!*****************************************************



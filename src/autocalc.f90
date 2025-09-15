subroutine autocalc()
use geninfo !contains cartes and natoms
use intrainfo !contains nquad and cent
implicit none
integer :: i,j,k,nqmax,sm
double precision, allocatable, dimension(:) :: x0, Alp, N
double precision, parameter :: trsh=1.d-1, zero=1.d-15, trs=5.d-1
double precision :: Density, rDM1, Gauss
double precision :: Val_r,Dens_a,Dens_b, r
double precision :: Pk,bn
double precision :: x1,x2,y1,y2,z1,z2
 nqmax=(natoms*(natoms-1))/2
 allocate(x0(nqmax))
 allocate(Alp(nqmax))
 allocate(N(nqmax))
 write(*,*) "In subroutine Autocalc"  
 write(*,*) "natoms=", natoms, nqmax  
 sm=0
 write(*,*) "start loop"
 do i=1,natoms-1
    do j=i+1,natoms
         x1=cartes(i,1)
         y1=cartes(i,2)
         z1=cartes(i,3)
         x2=cartes(j,1)
         y2=cartes(j,2)
         z2=cartes(j,3)
         sm=sm+1
         x0(sm)=dsqrt((x1-x2)**2.d0+(y1-y2)**2.d0+(z1-z2)**2.d0)
         !Dens_A=Density(x1,y1,z1)
         !Dens_B=Density(x2,y2,z2)
         !Pk=Dens_A*Dens_B-(rDM1(x1,y1,z1,x2,y2,z2)**2.d0)
         bn=an(i)*an(j)
         Pk=Pk/bn 
         Alp(sm)=pi*0.0625d0*((Pk*dexp(1.d0))**2.d0)*(bn**(-2.d0))
         N(sm)=dsqrt(pi)*0.25d0*Alp(sm)**(-1.5d0) !normalization factor if integral is 1
         N(sm)=bn/N(sm)                           !full factor for exponential
    end do
 end do

 open(unit=3,file='kaka')
 r=0.d0
 Val_r=0.d0
 do i=1,500
     r=r+0.03d0
     do j=1,nqmax
        !write(*,*) "x0,Alp,Nj"
        write(*,*) x0(j), Alp(j), N(j)
        Val_r=Val_r+Gauss(r,x0(j),Alp(j),N(j))
        write(*,*) "Gauss=", Gauss(r,x0(j),Alp(j),N(j))
     end do
     write(3,*) r, Val_r
     Val_r=0.d0
 end do      
      
end subroutine autocalc               

  function Gauss(r,x0,Alp,N)
  double precision :: Gauss
  double precision :: r,x0,Alp,N
  double precision :: minim, a1  
    a1=(dsqrt(Alp))**(-1.d0)
    minim=x0-a1
    if (r.ge.minim) then
        Gauss=N*((r-x0+a1)**2.d0)*dexp(-Alp*(r-x0+a1)**2.d0)
    else
        Gauss=0.d0 !neglect the peak from the left
    end if
   end function


!  function rDM1(x1,y1,z1,x2,y2,z2) !1-electron reduced density matrix (HF case)
!   use wfxinfo
!   implicit none
!   double precision :: rDM1
!   double precision :: MoOr
!   double precision :: x1, y1, z1, x2, y2, z2
!   integer :: i   !don't cross alpha and beta, if it's closed shell do /2
!       do i=1,noccmo
!             rDM1=rDM1+MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
!       end do
!       rDM1=rDM1/2.d0
!   end function

















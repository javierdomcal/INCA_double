!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine angular_nodes(gpt,nAng)                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Converts angular node option to number of angular points (nAng)   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 
integer, intent(in) :: gpt
integer, intent(out) :: nAng
 if (gpt.eq.1) nAng=6
 if (gpt.eq.2) nAng=14
 if (gpt.eq.3) nAng=26
 if (gpt.eq.4) nAng=38
 if (gpt.eq.5) nAng=50
 if (gpt.eq.6) nAng=74
 if (gpt.eq.7) nAng=86
 if (gpt.eq.8) nAng=110
 if (gpt.eq.9) nAng=146
 if (gpt.eq.10) nAng=170
 if (gpt.eq.11) nAng=194
 if (gpt.eq.12) nAng=230
 if (gpt.eq.13) nAng=266
 if (gpt.eq.14) nAng=302
 if (gpt.eq.15) nAng=350
 if (gpt.eq.16) nAng=434
 if (gpt.eq.17) nAng=590
 if (gpt.eq.18) nAng=770
 if (gpt.eq.19) nAng=974
 if (gpt.eq.20) nAng=1202
 if (gpt.eq.21) nAng=1454
 if (gpt.eq.22) nAng=1730
 if (gpt.eq.23) nAng=2030
 if (gpt.eq.24) nAng=2354
 if (gpt.eq.25) nAng=2702
 if (gpt.eq.26) nAng=3074
 if (gpt.eq.27) nAng=3470
 if (gpt.eq.28) nAng=3890
 if (gpt.eq.29) nAng=4334
 if (gpt.eq.30) nAng=4802
 if (gpt.eq.31) nAng=5294
 if (gpt.eq.32) nAng=5810
end subroutine angular_nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gridpoints(nradc,nAngc,sfalpha,nquad,cent,a,b,Ps) 
! Compute grid points to perform a multicenter radial integral from 0 to infinity    
!or from a to b (specify in the input file).
! Radial quadrature:Gauss-Legendre ; Angular quadrature: Gauss-Levedeb.
! Using Becke's weights to compute multicenter integrals.
! Neglects grid points with negative z using intracule symmetry, I(r)=I(-r)
! Neglects grid points with 0 weight.
! The total (reduced) grid ponints are stored in a 3*rrgrid matrix
use quadratures !contains output variables rg, rrg, and rgrid that will be used in intracule.f90       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
!global variables
integer, intent(in) :: nquad  !number of quad.
integer, dimension(nquad), intent(in) :: nradc, nAngc !number of radial points, number of anglular points
double precision, dimension(nquad),intent(in) :: sfalpha  !scaling factor for the radial points
double precision, dimension(nquad), intent(in) :: Ps
double precision, dimension(3,nquad) :: cent !center of the quadratures
double precision :: a,b !limits of the integral, 0 if the limits are from 0 to infty
!local variables
!1-levedeb quadrature
integer :: nang
double precision, allocatable, dimension(:,:) :: r_lb
double precision, allocatable, dimension(:) :: Wlb
double precision, allocatable, dimension(:) :: theta, phi
!2-legendre quadrature
integer :: nrad
double precision, allocatable, dimension(:) :: xl_i, wl_i
double precision, allocatable, dimension(:) :: radius !variable change from 0 to infty
!Legendre+Lebedev combination
double precision, allocatable, dimension(:) :: weight, srweight !total weight for each point (before neglecting points)
!3-becke vi
double precision, allocatable, dimension(:) :: w_beck
!-starting grid points 
double precision, allocatable, dimension(:,:) :: fgr, brrg, rrg 


integer :: i, ia, ir, sm, smp, sma, smnn, j, k, smr, smpr, ngrid, i1
integer :: np !number of points
double precision, parameter :: pi=4.d0*atan(1.d0)
double precision, parameter :: trsh=1.d-15, trsh2=1.d-16
double precision :: xs
!double precision :: a,b
!ntrsh=-trsh !set negative threshold

!count maximum number of points
write(*,*) "Ps=,", Ps(:)
np=0
do i=1,nquad
  np=np+nradc(i)*nAngc(i)
end do

maxgrid=np

write(*,*) "Error?"
allocate(weight(np))
allocate(fgr(3,np)) !first grid points
allocate(brrg(3,np))
allocate(smn(nquad)) !counts number of grid points per quadrature
write(*,*) "Error2="
sm=0 !count total grid points (1-->np)
smp=0 !count last grid point of previous centre
smr=0  !count sym reduced grid points
smpr=0  !count sym reduced gp of previous centre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i1=1,nquad       !loop over centres
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nrad=nradc(i1)
   nAng=nAngc(i1)
   write(*,*) "nrad=,nang=", nrad, nang
   !!Obtain radial nodes and weights (Gauss-Legendre)!!!!!!!!!!!!!!!!!!!!!
   allocate(xl_i(nrad))
   allocate(wl_i(nrad))
   call sub_GauLeg(-1.d0,1.d0,xl_i,wl_i,nrad)
   !do i=1,nrad
   !    write(*,*) xl_i(i), wl_i(i)
   !end do      
   !!!Obtain angular nodes and weights (Gauss-Lebedev)!!!!!!!!!!!!!!!!!!!!!!!  
   allocate(r_lb(3,nang))
   allocate(Wlb(nang))
   Wlb=0.d0
   r_lb=0.d0
   if (nAng.eq.6) call LD0006(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !6 grid points
   if (nAng.eq.14) call LD0014(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !14 grid points
   if (nAng.eq.26) call LD0026(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !26 grid points 
   if (nAng.eq.38) call LD0038(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !38 grid points
   if (nAng.eq.50) call LD0050(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !50 grid points
   if (nAng.eq.74) call LD0074(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !74 grid points
   if (nAng.eq.86) call LD0086(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !86 grid points
   if (nAng.eq.110) call LD0110(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !110 grid points
   if (nAng.eq.146) call LD0146(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !146 grid points
   if (nAng.eq.170) call LD0170(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !170 grid points
   if (nAng.eq.194) call LD0194(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !194 grid points 
   if (nAng.eq.230) call LD0230(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !230 grid points
   if (nAng.eq.266) call LD0266(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !266 grid points 
   if (nAng.eq.302) call LD0302(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !302 grid points
   if (nAng.eq.350) call LD0350(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !350 grid points
   if (nAng.eq.434) call LD0434(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !434
   if (nAng.eq.590) call LD0590(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !590
   if (nAng.eq.770) call LD0770(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !770
   if (nAng.eq.974) call LD0974(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !974
   if (nAng.eq.1202) call LD1202(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !1202
   if (nAng.eq.1454) call LD1454(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !1454
   if (nAng.eq.1730) call LD1730(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !1730
   if (nAng.eq.2030) call LD2030(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !2030
   if (nAng.eq.2354) call LD2354(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !2354
   if (nAng.eq.2702) call LD2702(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !2702
   if (nAng.eq.3074) call LD3074(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !3074
   if (nAng.eq.3470) call LD3470(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !3470
   if (nAng.eq.3890) call LD3890(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !3890
   if (nAng.eq.4334) call LD4334(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !4334
   if (nAng.eq.4802) call LD4802(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !4802
   if (nAng.eq.5294) call LD5294(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !5294
   if (nAng.eq.5810) call LD5810(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !5810
   allocate(theta(nang))
   allocate(phi(nang))
   do i=1,nAng   !compute angles     
        theta(i)= acos(r_lb(3,i)) 
        if (sin(theta(i)).ne.0.d0) then
             xs=r_lb(1,i)/sin(theta(i))
             if (xs.gt.1.d0) xs=1.d0
             if (xs.lt.-1.d0) xs=-1.d0
              phi(i)=acos(xs)
              if (r_lb(2,i).lt.0.d0) phi(i)=-phi(i)
        else
              phi(i)=0.d0
        end if  
        Wlb(i)=Wlb(i)*2.d0*pi !store 2pi factor on the weight (see solid angle integral)
   end do
   deallocate(r_lb)
   !!!!!!!!!!!!!Compute grid points of quadrature!!!!!!!!!!!
   nGrid=nAng*nrad
   allocate(radius(nrad))
   
   !compute radial points for the quadrature from 0 to infty or a to b
   if ((abs(a)+abs(b)).ne.0.d0) then
       write(*,*) "radius from a to b"    
       do ir=1,nrad    
          radius(ir)=(b-a)*0.5d0*xl_i(ir)+(a+b)*0.5d0
          wl_i(ir)=(b-a)*0.5d0*wl_i(ir) !radial integration (see pdf) 
       end do   
   else
       do ir=1,nrad    
          radius(ir)=(1.d0+xl_i(ir))/(1.d0-xl_i(ir))*sfalpha(i1) 
          wl_i(ir)=(2.d0*sfalpha(i1)/((1.d0-xl_i(ir))**2.d0))*wl_i(ir)*radius(ir)**2.d0    
        end do   
   end if  
   
   !compute grid points for all the becke centers
   smn(i1)=0
   do ir=1,nrad
        do ia=1,nAng
               sm=sm+1   !count total grid points    
               fgr(1,sm)=radius(ir)*cos(phi(ia))*sin(theta(ia))+cent(1,i1)
               fgr(2,sm)=radius(ir)*sin(phi(ia))*sin(theta(ia))+cent(2,i1)
               fgr(3,sm)=radius(ir)*cos(theta(ia))+cent(3,i1)
               !store total grid points (all quadratures)              
               if (fgr(3,sm).ge.0.d0) then
                  smn(i1)=smn(i1)+1 !sum the number of sym reduced points of each quadrature                  
               end if                  
        end do   
   end do  
   write(*,*) "Total number of gp after center",i1,"=", sm, ngrid
   !neglect points by symmetry!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   write(*,*) nGrid-(smn(i1)), "Points will be neglected in center", i1     
   rGrid=smn(i1)  !reduced grid points of current center
   !store reduced points and their total weight
   sm=smp   !start again from 1st point of center
   smr=smpr !start from 1st reduced point of center
   smnn=0   
   write(*,*) "*********Neglecting points by symmetry*******************"
   write(*,*) "zeta zero center"      
   do j=1,nrad
       do k=1,nAng  
          sm=sm+1                     !sum total grid points
          if (fgr(3,sm).ge.0.d0) then !reduce points by symmetry
              smr=smr+1               !sum sym reduced grid points   
              brrg(:,smr)=fgr(:,sm)   !store sym reduced points
              if (fgr(3,sm).lt.trsh) then    !z is 0, do not multiply by 2
                   weight(smr)=Wlb(k)*Wl_i(j)
              else
                   weight(smr)=2.d0*Wlb(k)*Wl_i(j) !z is positive, use sym (I(z)=I(-z))
              end if
          else if (fgr(3,sm).lt.0.d0) then 
              !sym neglected point    
              smnn=smnn+1
          end if  
       end do             
   end do   
   smpr=smr !store last reduced point of the quadrature
   smp=sm   !store last point of the quadrature
   write(*,*) smnn, "points have been neglected in center number", i1  
   deallocate(radius) 
   write(*,*) "Total number of gp after center", i1, "=", sm
   write(*,*) "Total number of reduced gp after center", i1, "=", smr
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   deallocate(theta)
   deallocate(phi)
   deallocate(Wlb)
   deallocate(xl_i)
   deallocate(wl_i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end do !end loop over quadratures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rgrid=sum(smn) !set total number of grid points
write(*,*) "There are", rgrid, "symmetry reduced grid points"
open(unit=3,file="all_points.dat")
do i=1,np
    write(3,*) fgr(1,i), fgr(2,i), fgr(3,i)
end do
close(3)
open(unit=3,file="sym_red_gp.dat")
do i=1,rgrid
    write(3,*) brrg(1,i), brrg(2,i), brrg(3,i)
end do
close(3)

allocate(rrg(3,rgrid))
allocate(srweight(rgrid))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sm=0
do i=1,nquad                   !store the sym reduced points in a smaller matrix 
 do j=1,smn(i)
   sm=sm+1
   rrg(:,sm)=brrg(:,sm)
   srweight(sm)=weight(sm)
 end do  
end do 

deallocate(brrg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (nquad.gt.1) then         !we have more than one quadrature centre
        write(*,*) "number of points per centre", smn(:)
        allocate(w_beck(rgrid)) !compute becke weights for each point    
        call becke(rrg,smn,rgrid,nquad,cent,w_beck,Ps)
        do i=1,rgrid
            srweight(i)=srweight(i)*w_beck(i) !store becke weight into the total one
        end do
        deallocate(w_beck)
       !reduce points with 0 weight
        sm=0
        sma=0
        do i=1,rgrid
            if ((srweight(i).le.trsh2)) then !not 1st center and little weight
               sm=sm+1                     !sum number of neglected points  
               if (i.gt.smn(1)) then
                  sma=sma+1 !count neglected points in 2nd center 
               end if        
            end if        
        end do
        rrgrid=rgrid-sm !number of sym and weight reduced points
else !only 1 quadrature centre--> No Becke    
         rrgrid=rgrid
end if

write(*,*) sma, "neglected points of zero weight"
write(*,*) sm, "neglected points of zero weight"

allocate(rweight(rrgrid)) !reduced weight
allocate(rrrg(3,rrgrid))  !doubly reduced points
sm=0
write(*,*) "Final number of grid points=", rrgrid
do i=1,rgrid
     if ((srweight(i).gt.trsh2)) then !remove points with low (zero) weight
           sm=sm+1 
           rweight(sm)=srweight(i) 
           rrrg(:,sm)=rrg(:,i)           
     end if
end do
!
open(unit=3,file="sym_weight_red_gp.dat")
do i=1,rrgrid
    write(3,*) rrrg(1,i), rrrg(2,i), rrrg(3,i)
end do
close(3)




deallocate(weight)
deallocate(rrg)
close(4)
end subroutine 


subroutine centercalc()
use geninfo !contains cartes and natoms
use intrainfo !contains nquad and cent
implicit none
integer :: i,j,k, nqmax, sm, sm1, smp
double precision, allocatable, dimension(:,:) :: c1
double precision, allocatable, dimension(:) :: dist
integer, allocatable, dimension(:) :: nelec_c1, nelec_cent
double precision, parameter :: trsh=1.d-1, zero=1.d-15, trs=5.d-1
double precision :: z, k_val
logical :: equal
!integer, allocatable, dimension(:) :: eqc
nqmax=(natoms*(natoms-1))+1
allocate(c1(3,nqmax)); allocate(nelec_c1(nqmax))
c1(:,1)=0.d0 !first center always at zero    
sm=1
do j=1,natoms-1
   do k=j+1,natoms
      if (j.ne.k) then
         if ((an(j).gt.1).and.(an(k).gt.1)) then !only for non-hydrogen atoms      
            sm=sm+1      
            c1(:,sm)=cartes(j,:)-cartes(k,:)
            nelec_c1(sm)=(an(i)-chrg(i))*(an(j)-chrg(j))
            z=dabs(c1(3,sm))
            if (z.lt.trs) c1(3,sm)=0.d0 !to exploit symmetry
            smp=sm  
            sm=sm+1 !count centres  
            c1(3,sm)=-c1(3,smp)
            c1(1,sm)=-c1(1,smp) !I(x,y,z)=I(-x,-y,-z)
            c1(2,sm)=-c1(2,smp)      
            nelec_c1(sm)=nelec_c1(smp)
         end if  
      end if  
   end do
end do
!check if there are equal centers
nqmax=sm
sm=0
do i=2,nqmax-1
   do j=i+1,nqmax
      sm1=sm1+1
      z=dabs(c1(1,i))+dabs(c1(2,i))+dabs(c1(3,i)) !sum x,y,z axis (if 0->equivalent centres exist)
      if (z.ne.0.d0) then                      !this center was not repeated
         equal=all(abs(c1(:,i)-c1(:,j)).lt.trs)
         if (equal) then             !center j is equal to center i
              sm=sm+1
              c1(:,j)=0.d0                     !set center j to zero
         end if   
      end if
   end do
end do
nquad=nqmax-sm
allocate(cent(3,nquad)); allocate(nelec_cent(nquad))
sm=1
sm1=1
do i=2,nqmax
   z=dabs(c1(1,i))+dabs(c1(2,i))+dabs(c1(3,i))
   if (z.ne.0.d0) then  !store all non-equivalent centres in a lower-dimensional array
      sm1=sm1+1
      cent(:,sm1)=c1(:,i)
      nelec_cent(sm1)=nelec_c1(i)
   end if
end do 
allocate(nradc(nquad)); allocate(nangc(nquad)); allocate(sfalpha(nquad)); allocate(Ps(nquad))
allocate(dist(nquad))
do i=1,nquad !same nrad and nang for all centers (you might want to change this in the future)
   nradc(i)=nrad
   nangc(i)=nang
end do     
do i=1,nquad
   dist(i)=sqrt(sum(c1(:,i)**2))
end do
sfalpha(1)=minval(dist,dim=1,mask=(dist.gt.zero))/2.d0
Ps(1)=1.d0
sm=0
k_val=nrad/nquad
do i=2,nquad
   if (dist(i).lt.zero) then
      sfalpha(i)=((dist(i)-dist(i-1))/2.d0)*nradc(i)/K_val
      Ps(i)=nelec_cent(i)
   end if   
end do
write(*,*) "Inteegration parameters for each center:"
write(*,*)  "Center (x,y,z)    nrad    nang    sfalpha    Weight"
do i=1,nquad
   write(*,*) cent(1,i), cent(2,i), cent(3,i), nradc(i), nangc(i), sfalpha(i), Ps(i)
end do

end subroutine centercalc               


subroutine centercalc()
use geninfo !contains cartes and natoms
use intrainfo !contains nquad and cent
implicit none
integer, allocatable, dimension(:) :: eqc
integer :: i,j,k, nqmax, sm, sm1, smp
double precision, allocatable, dimension(:,:) :: c1
double precision, parameter :: trsh=1.d-2, zero=1.d-15
double precision :: diff, z
!integer, allocatable, dimension(:) :: eqc
 nqmax=(natoms*(natoms-1))+1
 allocate(c1(3,nqmax))
 c1(:,1)=0.d0 !first center always at zero    
 sm=1
 do j=1,natoms-1
    do k=j+1,natoms
       if (j.ne.k) then
         if ((an(j).gt.1).and.(an(k).gt.1)) then !only for non-hydrogen atoms      
            sm=sm+1      
            z=cartes(j,3)-cartes(k,3)
            c1(:,sm)=cartes(j,:)-cartes(k,:)
            if (abs(z).gt.zero) then
               smp=sm  
               sm=sm+1 !count centres  
               c1(3,sm)=-c1(3,smp)
               c1(1,sm)=-c1(1,smp) !I(x,y,z)=I(-x,-y,-z)
               c1(2,sm)=-c1(2,smp)
             end if        
         end if  
       end if  
    end do
 end do

 !check if there are equal centers
 allocate(eqc(sm))
 nqmax=sm
 eqc=0
 sm=0
 do i=2,nqmax-1
   do j=i+1,nqmax
      sm1=sm1+1
      diff=abs(sum(c1(:,i)-c1(:,j)))
      if (diff.lt.trsh) then
              sm=sm+1
              eqc(sm)=j !store equal center
      end if
    end do
  end do

  nquad=nqmax-sm
  allocate(cent(3,nquad))
  sm=1
  sm1=0
  do i=1,nqmax
      if (eqc(sm).eq.i) then   
           sm=sm+1
      else
         sm1=sm1+1
         cent(:,sm1)=c1(:,i)
      end if
  end do
  
end subroutine centercalc               




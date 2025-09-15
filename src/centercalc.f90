subroutine centercalc()
use geninfo !contains cartes and natoms
use intrainfo !contains nquad and cent
implicit none
integer :: i,j,k, nqmax, sm, sm1, smp
double precision, allocatable, dimension(:,:) :: c1
double precision, parameter :: trsh=1.d-1, zero=1.d-15, trs=5.d-1
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
            c1(:,sm)=cartes(j,:)-cartes(k,:)
            z=dabs(c1(3,sm))
            if (z.lt.trs) c1(3,sm)=0.d0 !to exploit symmetry
            smp=sm  
            sm=sm+1 !count centres  
            c1(3,sm)=-c1(3,smp)
            c1(1,sm)=-c1(1,smp) !I(x,y,z)=I(-x,-y,-z)
            c1(2,sm)=-c1(2,smp)      
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
         diff=dabs(sum(c1(:,i)-c1(:,j)))
         if (diff.lt.trsh) then             !center j is equal to center i
              sm=sm+1
              c1(:,j)=0.d0                     !set center j to zero
         end if   
      end if
    end do
 end do

 nquad=nqmax-sm
 allocate(cent(3,nquad))
 sm=1
 sm1=1
 do i=2,nqmax
      z=dabs(c1(1,i))+dabs(c1(2,i))+dabs(c1(3,i))
      if (z.ne.0.d0) then  !store all non-equivalent centres in a lower-dimensional array
           sm1=sm1+1
           cent(:,sm1)=c1(:,i)
      end if
 end do
  
end subroutine centercalc               




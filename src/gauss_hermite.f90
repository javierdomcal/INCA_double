subroutine gauherm(Lrtot,n) 
use quadratures !module contains gauss-hermite nodes and weights
!performs gauss-hermite quadrature                          
!we will obtain the nodes and weights depending on the degree of the polinomial (2n-1) --> (n) 
!(n is nx ny or nz in intracule.f90)
 implicit none
 !global variable
 integer, intent(in) :: Lrtot !degree of polynomial
 integer, intent(out) :: n    !number of gauss-hermite nodes
 !local variables
 integer :: factn, i2
 double precision :: Cons
 double precision, parameter :: pi=4.d0*atan(1.d0)
 
 if (Lrtot.gt.0) then
 !compute the number of nodes for exact Hermite quadrature
     if(MOD(Lrtot,2).eq.0) then !even 
           n=(Lrtot/2) +1 !nx->number of nodes
     else 
           n=(Lrtot+1)/2 
     end if
 else
     !only one node
           n=1
 end if  
 
 allocate(rh(n))  !nodes
 allocate(w_r(n)) !weights
 
 !store rh values depending on pol. degree
  
  if (n.eq.1) then !1 node-->degree of the pol is 0 or 1 
   rh(1)=0.0000000000000d0   
  else if (n.eq.2) then !2 nodes-->degree of the pol is 3 or 2                
   rh(2)=0.7071067811865475 
   rh(1)=-0.7071067811865475 
  else if (n.eq.3) then !3 nodes-->degree of pol is 5 or 4
   rh(2)=0.d0
   rh(3)= 1.224744871391589 
   rh(1)=-1.224744871391589  
  else if (n.eq.4) then   !7 or 6
   rh(3)=0.5246476232752903 
   rh(4)=1.650680123885785 
   rh(2)=-0.5246476232752903 
   rh(1)=-1.650680123885785 
  else if (n.eq.5) then    !9 or 8
   rh(1)=-2.020182870456086
   rh(2)=-0.9585724646138185
   rh(3)=0.000000000000000d0
   rh(4)=2.020182870456086
   rh(5)=0.9585724646138185
  else if (n.eq.6) then    !11 or 10
   rh(4)=0.4360774119276165 
   rh(5)=1.335849074013697 
   rh(6)=2.350604973674492   
   rh(3)=-0.4360774119276165 
   rh(2)=-1.335849074013697 
   rh(1)=-2.350604973674492
  else if (n.eq.7) then   !13 or 12
   rh(4)=0.d0
   rh(5)=0.8162878828589647 
   rh(6)=1.673551628767471 
   rh(7)=2.651961356835233
   rh(3)=-0.8162878828589647 
   rh(2)=-1.673551628767471 
   rh(1)=-2.651961356835233
  else if (n.eq.8) then    !15 or 14
   rh(5)=0.3811869902073221 
   rh(6)=1.157193712446780 
   rh(7)=1.981656756695843 
   rh(8)=2.930637420257244
   rh(4)=-0.3811869902073221 
   rh(3)=-1.157193712446780 
   rh(2)=-1.981656756695843 
   rh(1)=-2.930637420257244
  else if (n.eq.9) then     !17 or 16
   rh(6)=0.7235510187528376 
   rh(7)=1.468553289216668 
   rh(8)=2.266580584531843 
   rh(9)=3.190993201781528 
   rh(5)=0.d0
   rh(4)=-0.7235510187528376 
   rh(3)=-1.468553289216668 
   rh(2)=-2.266580584531843 
   rh(1)=-3.190993201781528 
  else if (n.eq.10) then     !19 or 18
    rh(6)=0.3429013272237046 
    rh(7)=1.036610829789514 
    rh(8)=1.756683649299882 
    rh(9)=2.532731674232790 
    rh(10)=3.436159118837738 
    rh(5)=-0.3429013272237046 
    rh(4)=-1.036610829789514 
    rh(3)=-1.756683649299882 
    rh(2)=-2.532731674232790 
    rh(1)=-3.436159118837738 
  else if (n.eq.11) then     !degree of pol is 21 or 20
    rh(6)=0.d0
    rh(7)=0.6568095668820998 
    rh(8)=1.326557084494933 
    rh(9)=2.025948015825755 
    rh(10)=2.783290099781652 
    rh(11)=3.668470846559583 
    rh(5)=-0.6568095668820998 
    rh(4)=-1.326557084494933 
    rh(3)=-2.025948015825755 
    rh(2)=-2.783290099781652 
    rh(1)=-3.668470846559583 
  else 
   write(*,*) "Error, total ang. momenta greater than 3"
  end if
  
  factn=1     !compute weights using the equation
  do i2=1,n   
   factn=factn*i2  
  end do
  Cons=((2.d0)**(dble(n-1)) * dble(factn)* sqrt(pi))/dble(n**2)
   
  do i2=1,n
   w_r(i2)=cons*(1.d0/(Hermite(rh(i2),n))**2.d0)  
  end do  

  contains
    function Hermite(x,n) !
    double precision :: Hermite
    integer, intent(in) :: n !numer of nodes of pol
    double precision, intent(in) :: x !roots of Hermite pol
   
   if (n.eq.1) then !0 deg
     Hermite=1.d0
   else if (n.eq.2) then
     Hermite=2.d0*x 
   else if (n.eq.3) then
     Hermite=4.d0*(x**2.d0) -2.d0
   else if (n.eq.4) then
     Hermite=8.d0*(x**3.d0) -12.d0*x
   else if (n.eq.5) then
     Hermite=16.d0*(x**4.d0)-48.d0*(x**2.d0)+12.d0  
   else if (n.eq.6) then
     Hermite=32.d0*(x**5.d0)-160.d0*(x**3.d0)+120.d0*x  
   else if (n.eq.7) then
     Hermite=64.d0*(x**6.d0)-480.d0*(x**4.d0)+720.d0*(x**2.d0)-120.d0
   else if (n.eq.8) then
     Hermite=128.d0*(x**7.d0)-1344.d0*(x**5.d0)+3360.d0*(x**3.d0)-1680.d0*x    
   else if (n.eq.9) then
     Hermite=256.d0*(x**8.d0)-3584.d0*(x**6.d0)+13440.d0*(x**4.d0)-13440.d0*(x**2.d0)+1680.d0  
   else if (n.eq.10) then
     Hermite=512.d0*(x**9.d0)-9216.d0*(x**7.d0)+48384.d0*(x**5.d0)-80640.d0*(x**3.d0)+30240.d0*x
   else if (n.eq.11) then
   Hermite=1024.d0*(x**10.d0)-23040.d0*(x**8.d0)+161280.d0*(x**6.d0)-403200.d0*(x**4.d0)+302400.d0*(x**2.d0)-30240.d0  
   end if       
   end function

end subroutine gauherm 
   

   
   
   


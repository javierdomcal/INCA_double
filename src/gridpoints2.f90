!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gridpoints2(nblock,tart,step,n_an_per_part)                                                    !
! Compute grid points for with manually given radius (in the input file)!   
use quadratures !contains weights, nradi, radi !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
!global variables
integer, intent(in) :: nblock !number of parts
double precision, dimension(nblock), intent(in) :: step
integer, intent(in), dimension(nblock) :: n_an_per_part !number of steps per part
double precision, intent(in), dimension(2,nblock) :: tart !range of each Block
!local variables
!1-levedeb quadrature
double precision, allocatable, dimension(:,:) :: r_lb
double precision, allocatable, dimension(:) :: Wlb
!double precision, allocatable, dimension(:) :: theta, phi
double precision, allocatable, dimension(:) :: pweight
!grid points
integer, allocatable,dimension(:) :: npb
integer :: nr_total
double precision, allocatable, dimension(:,:) :: gpt
double precision :: z
integer :: i, ia, ir, sm, smnn, j, smp, n_an, smrad
double precision, parameter :: pi=4.d0*datan(1.d0)
double precision, parameter :: trsh=1.d-15, trsh2=1.d-16, tol=dsqrt(epsilon(1.d0))
!double precision :: a,b
double precision :: xs
double precision :: r_start,r_end
write(*,*) nblock
write(*,*) "tol=", tol
allocate(npb(nblock))
nr_total=0
do i=1,nblock
    r_start=tart(1,i)
    r_end=tart(2,i)
    !compute number of points
    npb(i)=int((r_end-r_start)/step(i))+1
    nr_total = nr_total + npb(i)
end do    

nradi = nr_total
allocate(radi(nradi))
allocate(smn(nradi))
radi = 0.d0
smn = 0.d0

sm = 0.d0
rpgrid = 0
do i = 1, nblock
    do j = 0, npb(i)-1
        sm = sm + 1
        radi(sm) = tart(1,i) + step(i)*dble(j)
    end do
    rpgrid = rpgrid + npb(i) * n_an_per_part(i)
end do

allocate(gpt(3,rpgrid))
allocate(pweight(rpgrid))
smrad=0
sm=0 !sum for grid points
smnn=0 !sum symmetry reduced points
do i=1,nblock !loop for each radius fragment
    n_an=n_an_per_part(i)
    allocate(r_lb(3,n_an))
    allocate(Wlb(n_an))
    !!!Obtain angular nodes and weights (Gauss-Lebedev)!!!!!!!!!!!!!!!!!!!!!!!   
    if (n_an.eq.6) call LD0006(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !6 grid points
    if (n_an.eq.14) call LD0014(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !14 grid points
    if (n_an.eq.26) call LD0026(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !26 grid points 
    if (n_an.eq.38) call LD0038(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !38 grid points
    if (n_an.eq.50) call LD0050(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !50 grid points
    if (n_an.eq.74) call LD0074(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !74 grid points
    if (n_an.eq.86) call LD0086(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !86 grid points
    if (n_an.eq.110) call LD0110(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !110 grid points
    if (n_an.eq.146) call LD0146(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !146 grid points
    if (n_an.eq.170) call LD0170(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !170 grid points
    if (n_an.eq.194) call LD0194(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !194 grid points 
    if (n_an.eq.230) call LD0230(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !230 grid points
    if (n_an.eq.266) call LD0266(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !266 grid points 
    if (n_an.eq.302) call LD0302(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !302 grid points
    if (n_an.eq.350) call LD0350(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !350 grid points
    if (n_an.eq.434) call LD0434(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !434
    if (n_an.eq.590) call LD0590(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !590
    if (n_an.eq.770) call LD0770(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !770
    if (n_an.eq.974) call LD0974(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !974
    if (n_an.eq.1202) call LD1202(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !1202
    if (n_an.eq.1454) call LD1454(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !1454
    if (n_an.eq.1730) call LD1730(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !1730
    if (n_an.eq.2030) call LD2030(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !2030
    if (n_an.eq.2354) call LD2354(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !2354
    if (n_an.eq.2702) call LD2702(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !2702
    if (n_an.eq.3074) call LD3074(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !3074
    if (n_an.eq.3470) call LD3470(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !3470
    if (n_an.eq.3890) call LD3890(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !3890
    if (n_an.eq.4334) call LD4334(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !4334
    if (n_an.eq.4802) call LD4802(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !4802
    if (n_an.eq.5294) call LD5294(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !5294
    if (n_an.eq.5810) call LD5810(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !5810 
    !allocate(theta(n_an))
    !allocate(phi(n_an))
    !do j=1,n_an   !compute angles     
    !    theta(j)= dacos(r_lb(3,j))
        !if (1.d0 - dabs(r_lb(3,j)).gt.tol) then
            !xs=r_lb(1,j)/dsin(theta(j))
            !if (xs.gt. 1.d0) xs=1.d0
            !if (xs.lt.-1.d0) xs=-1.d0
    !        phi(j) = datan2(r_lb(2,j), r_lb(1,j))
            !phi(j)=dacos(xs)
            !if (r_lb(2,j).lt. 0.d0) phi(j)=-1.0d0*phi(j)
        !else
        !      phi(j)=0.d0
        !end if
        !Wlb(j)=Wlb(j)*2.d0*pi !store 2pi factor on the weight (see solid angle integral)
    !end do   
     !angles and weights are computed for each block

    !!!!!!!!Compute grid points!!!!!!!!!!!!!!!!!!!!!!!!!!    
    do ir=1,npb(i) !loop for each radius inside the block
        smrad=smrad+1 !sum over radius
        do ia=1,n_an  !sum over angular points
            sm=sm+1 !sum grid points
            z=radi(smrad)*r_lb(3,ia) !compute z
            if (z.ge.-tol) then
                smnn=smnn+1 !sum the number of reduced points of each quadrature
                smn(smrad)=smn(smrad)+1 !sum number of points per radius
                gpt(1,smnn)=radi(smrad)*r_lb(1,ia)!dcos(phi(ia))*dsin(theta(ia))
                gpt(2,smnn)=radi(smrad)*r_lb(2,ia)!dsin(phi(ia))*dsin(theta(ia))
                gpt(3,smnn)=z
                if (abs(z) < tol) then
                    pweight(smnn)=Wlb(ia)
                else
                    pweight(smnn)=2.d0*Wlb(ia)
                end if   
            end if
        end do   
    end do  
    deallocate(Wlb)
    deallocate(r_lb)
end do   !end loop for radius fragment
 
write(*,*) "Original number of grid points", sm
write(*,*) "Number of grid points after I(r)=I(-r)", smnn
write(*,*) "Number of grid points per radi", smn
allocate(rpg(3,smnn))
allocate(w_ang(smnn))
w_ang = 0.0d0
rpg = 0.0d0 
do i=1,smnn  !store the data in reduced size matrix
    rpg(:,i)=gpt(:,i) 
    w_ang(i)=pweight(i)          
end do
rgrid=smnn
end subroutine 

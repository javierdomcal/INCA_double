!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine intracule()  !Computes vector or radial intracule, or radial integration
                        !need .wfx and .dm2p as input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use geninfo !information about the primitive functions
use intrastuff !subroutines and functions to compute the intracule 
use intrainfo !information from the input
use quadratures !nodes and weights of the quadratures, +total grid points
implicit none
!!!!!!!!!!!!!!!local variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: kk1, kk2 !integers we do not want to read from .dm2
integer :: dfact    !double factorial
integer :: ii,iii,sum, ig, ir, summ, np, sm, smm !some integers
double precision :: pot 
double precision :: R_i_k_2, R_j_l_2  !R_ik=(R_i-R_k)^2
double precision :: Aa_ijkl !grid independent part of the intracule
double precision :: n_prim_i, n_prim_j,n_prim_k, n_prim_l !normalization of primitives for DM2
double precision :: n_prim_t !total normalization factor
double precision :: DMval !value of DM2
double precision :: A_ind !eq.18 (grid independent part)
double precision :: screen1, screen2, lim !integral screenings
double precision :: U_x, U_y, U_z !coef. for 2nd integral screening
double precision, allocatable, dimension(:) :: C_x, C_y, C_z !coefficients of V 
double precision :: V_x, V_y, V_z !eq 16
!Gauss-hermite quadrature
integer :: nn !number of gauss hermite nodes (for x y and z)
integer, dimension(3) :: Lrtot  !total angular momentum and number of nodes
!GRID POINTS
integer :: ngrid !number of grid points
double precision, allocatable, dimension(:,:) :: r !grid points
double precision, dimension(3) :: rp  !prima points(X', Y', Z')
!intracules
double precision :: r_integral                         !radial integral
double precision, allocatable, dimension(:) :: r_intra !radial intracule
double precision, allocatable, dimension(:) :: I_vec   !intracule at a point
!check accuracy of calculations
double precision :: trace_DM2prim, trDM2 !normalized and not normalized DM2prim
integer :: npairs      !number of electron pairs
real :: T1, T2, T3, T4, TT1, TT2, TT3, TT4, TT5 !time check
real :: Tread, T1screen, T2screen, Tgrid

 lim=thresh*(dble(nprim)*(dble(nprim)+1.d0)*0.5d0)**(-1.d0) !limit for the 1st integral screening
 
 call cpu_time(T1)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!Obtain grid points!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (radial_integral) then
     call gridpoints(nradc,nAngc,sfalpha,nquad,cent,a,b,Ps) !obtain the grid points
     allocate(r(3,rrgrid))
     r=rrrg
     ngrid=rrgrid
     deallocate(rrrg)
 else if (radial_plot) then
     call gridpoints2(nblock,tart,stp,n_an_per_part)
     allocate(r(3,rgrid))
     r=rpg
     ngrid=rgrid
 else if (cubeintra) then !vectorial plot
     call gridpoints3(center_i,step_i,np_i)    
     allocate(r(3,rgrid))
     r=rg
     ngrid=rgrid
 end if    
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 allocate(I_vec(ngrid))
 I_vec=0.d0
 trace_DM2prim=0.d0 
 trDM2=0.d0
 call cpu_time(T2)

! call intracalc(r,ngrid)
 !subroutine intracalc(Computes intracule function with the given points)
 open(unit=5,file=dm2name, form='unformatted',access='stream') !open binary file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!Start loop over primitive quartets!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!Following Cioslowski and Liu algorithm!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 sum=0     !quartets skipped in the 1st screening
 summ=0    !quartets skipped in the 2nd screeening         
 write(*,*) "starting loop for primitives"
 write(*,*) "Loop over", ngrid, "grid points" 
 rewind(5)    !start reading from the begining of the dm2 file
 smm=0
 trDM2=0.d0  !sum of all the DM2 terms.
 trace_DM2prim=0.d0 !sum of all the normalized DM2 terms.
!Time check
Tread=0.d0
T1screen=0.d0
T2screen=0.d0
Tgrid=0.d0

 do while (.true.)  !loop for primitive quartets.
 
          call cpu_time(TT1)
    
          read(5,end=100) kk1,i,j,k,l,DMval,kk2 !read a line from binary file .dm2
          if (i.eq.0) goto 100 !file is finished   
          
          call cpu_time(TT2) 
          Tread=Tread+(TT2-TT1)

          trDM2=trDM2+DMval

          smm=smm+1
          !normalization of DM2--> Normalize the primitives
          N_prim_i=(2.d0*Alpha(i)/pi)**(0.75d0)&
          *sqrt(((4.d0*Alpha(i))**(dble(TMN(i,1)+TMN(i,2)+TMN(i,3))))*&
          dble(dfact(2*TMN(i,1)-1)*dfact(2*TMN(i,2)-1)*dfact(2*TMN(i,3)-1))**(-1.d0))  
   
          N_prim_j=(2.d0*Alpha(j)/pi)**(0.75d0)&
          *sqrt(((4.d0*Alpha(j))**(dble(TMN(j,1)+TMN(j,2)+TMN(j,3))))*&
          dble(dfact(2*TMN(j,1)-1)*dfact(2*TMN(j,2)-1)*dfact(2*TMN(j,3)-1))**(-1.d0)) 
   
          N_prim_k=(2.d0*Alpha(k)/pi)**(0.75d0)&
          *sqrt(((4.d0*Alpha(k))**(dble(TMN(k,1)+TMN(k,2)+TMN(k,3))))*&
          dble(dfact(2*TMN(k,1)-1)*dfact(2*TMN(k,2)-1)*dfact(2*TMN(k,3)-1))**(-1.d0)) 
         
          N_prim_l=(2.d0*Alpha(l)/pi)**(0.75d0)&
          *sqrt(((4.d0*Alpha(l))**(dble(TMN(l,1)+TMN(l,2)+TMN(l,3))))*&
          dble(dfact(2*TMN(l,1)-1)*dfact(2*TMN(l,2)-1)*dfact(2*TMN(l,3)-1))**(-1.d0)) 
   
          !compute DMval with the normalization of primitives
          n_prim_t=N_prim_i*N_prim_j*N_prim_k*N_prim_l
          DMval=n_prim_t*DMval 
          trace_DM2prim=trace_DM2prim+DMval !sum all the DM2 quartets to check accuracy
           
          !compute the first variables

          a_ik=Alpha(i)+Alpha(k)  
          a_jl=Alpha(j)+Alpha(l)                  !eqn. 10                         
          e_ik=Alpha(i)*Alpha(k)*a_ik**(-1.d0)
          e_jl=Alpha(j)*Alpha(l)*a_jl**(-1.d0)   
                
          R_i_k_2=(Cartes(Ra(i),1)-Cartes(Ra(k),1))**2.d0+& !this is needed to compute A_ijkl (eqn. 18)
                  (Cartes(Ra(i),2)-Cartes(Ra(k),2))**2.d0+&      !R_i_k_2=(R_i-R_k)²
                  (Cartes(Ra(i),3)-Cartes(Ra(k),3))**2.d0    
                
          R_j_l_2=(Cartes(Ra(j),1)-Cartes(Ra(l),1))**2.d0+&
                  (Cartes(Ra(j),2)-Cartes(Ra(l),2))**2.d0+&
                  (Cartes(Ra(j),3)-Cartes(Ra(l),3))**2.d0
                   
         !!!!1st integral screening!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                   
         screen1=abs(DMval)*sqrt(J_ik(i,k,R_i_K_2)*J_jl(j,l,R_j_l_2))       
         
         call cpu_time(TT3)
         T1screen=T1screen-(TT3-TT2)         
         if (screen1.ge.lim) then    
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 
                 !compute the other variables (eqn. 12)
                 a_ijkl=a_ik+a_jl
                 e_ijkl=(a_ik*a_jl)*a_ijkl**(-1.d0)
                 sqe=sqrt(e_ijkl)
                 do ii=1,3                   !X_ik,Y_ik,Z_ik,...
                       R_ik(ii)=(Alpha(i)*Cartes(Ra(i),ii)+Alpha(k)*Cartes(Ra(k),ii))*a_ik**(-1.d0)
                       R_jl(ii)=(Alpha(j)*Cartes(Ra(j),ii)+Alpha(l)*Cartes(Ra(l),ii))*a_jl**(-1.d0)
                       R_ijkl(ii)=(a_ik*R_ik(ii)+a_jl*R_jl(ii))*a_ijkl**(-1.d0)
                 end do                      
                 Alf_ijkl=0.5d0*(a_ik-a_jl)*a_ijkl**(-1.d0)

                 !compute Aa_ijkl (the grid independent part)  !eq.18
                 A_ind=(a_ijkl)**(-1.5d0)* exp(-e_ik*R_i_k_2-e_jl*R_j_l_2)
                 Aa_ijkl=DMval*A_ind 



                 !!!!!!!!!!!!Calculate coeficients of V!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 sm=0
                 U_x=0.d0; U_y=0.d0; U_z=0.d0
                 do ii=1,3
                       Lrtot(ii)=TMN(i,ii)+TMN(j,ii)+TMN(k,ii)+TMN(l,ii) !Lrtot=degree of eqn 15 (I don't use Lmax)
                       call gauherm(Lrtot(ii), nn) !obtain Gauss-Hermite nodes(rh) and weights(w_r) (2L+1)   
                       np=Lrtot(ii)+1              !np is the number of coefficients to represent the polynomial V
                       if (ii.eq.1) call polycoef(C_x, np, nn, ii) 
                       if (ii.eq.2) call polycoef(C_y, np, nn, ii)
                       if (ii.eq.3) call polycoef(C_z, np, nn, ii)
                       deallocate(rh)
                       deallocate(w_r)                            
                       !evaluates Vr at Lxtot+1=np points to create the augmented matrix M
                       !diagonalizes M to obtain the coeficients of V --> C 
                        do iii=1,np
                            if (ii.eq.1) U_x=U_x+C_x(iii)*w_m(iii)  !compute U coefficients
                            if (ii.eq.2) U_y=U_y+C_y(iii)*w_m(iii)  !for second integral screening
                            if (ii.eq.3) U_z=U_z+C_z(iii)*w_m(iii)  !eq.43 of the paper  
                        end do
                 end do
                                  
                 !!!!!!!!!!2nd Integral Screening!!!!!!!!!!!!!!!!!!!!!                      
                 screen2=abs(Aa_ijkl*U_x*U_y*U_z) !eqn.40
                 call CPU_time(TT4)
                 T2screen=T2screen+(TT4-TT3)                
                 if (screen2.ge.lim) then 
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                             
                            !loop over grid points                                
                            do ig=1,nGrid                                             
                                   V_x=0.d0; V_y=0.d0; V_z=0.d0  
                                   do ii=1,3                                       !loop for x y and z  
                                        rp(ii)=sqe*(r(ii,ig)+r_ik(ii)-r_jl(ii))   !compute R'(eq. 19), 
                                                                                   !using total grid points
                                        pot=dble(Lrtot(ii))
                                        !!!!!!!Compute Vx, Vy and Vz!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                        do iii=1,Lrtot(ii)+1   !Ltot+1 is the number of coefficients we have
                                                if (pot.ne.0) then
                                                     if (ii.eq.1) V_x=V_x+C_x(iii)*rp(1)**(pot)
                                                     if (ii.eq.2) V_y=V_y+C_y(iii)*rp(2)**(pot)
                                                     if (ii.eq.3) V_z=V_z+C_z(iii)*rp(3)**(pot)   
                                                     pot=pot-1.d0
                                                else
                                                     if (ii.eq.1) V_x=V_x+C_x(iii) !when pot is 0
                                                     if (ii.eq.2) V_y=V_y+C_y(iii)
                                                     if (ii.eq.3) V_z=V_z+C_z(iii)        
                                                end if 
                                        end do  !end loop over the polynomial coefficients
                                   end do       !end loop over x,y,z
                                   !!!!!!!!!!!!!Calculate intracule at a point!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   I_vec(ig)=I_vec(ig)+Aa_ijkl*exp(-(rp(1)**2.d0+rp(2)**2.d0+rp(3)**2.d0))*V_x*V_y*V_z
                                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  
                            end do       !End loop over grid points
                            call CPU_time(TT5)
                            Tgrid=Tgrid+(TT5-TT4)
                 else                     
                              summ=summ+1  !count quartets that do not pass 2nd screening 
                 end if 
                 deallocate(C_x)
                 deallocate(C_y)   !deallocate polynomial coefficients
                 deallocate(C_z)                                                                
       else               
               sum=sum+1    !count quartets that do not pass 1st screening                                 
       end if                
 end do !end loop over quartets 
 100 continue           !it comes here when .dm2 file is finished
 write(*,*) "Ended loop for primitives"
 !end subroutine intracalc 
 !as output gives a vector with the intracule at the given points
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!Compute the integrals using the quadrature !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!weights and the vectorial intracule        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call cpu_time(T3)
 if (radial_plot) then        !compute radial intracule
      open(unit=3, file=r_plot_name)   
      allocate(r_intra(nradi)) 
      ig=0
      sm=0
      r_intra=0.d0
      ir=0
      do i=1,nradi        
            do k=1,smn(i) !sum reduced angular points per radi  
                sm=sm+1            
                r_intra(i)=r_intra(i)+w_ang(sm)*I_vec(sm)!Perform the angular quadrature 
            end do
            r_intra(i)=radi(i)**(2.d0)*r_intra(i)
            write(3,*) radi(i), r_intra(i)
      end do 
      deallocate(r_intra)
      close(3)
 end if     
 if (radial_integral) then !integral of the intracule
     r_integral=0.d0
     do i=1,rrgrid
         r_integral=r_integral+rweight(i)*I_vec(i)
     end do 
     deallocate(rweight)
 end if
 if (cubeintra) then
         open(unit=3, file=cubeintraname)
         write(3,*) "CUBE FILE"
         write(3,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
         write(3,*) natoms, center_i(1), center_i(2), center_i(3)
         write(3,*) np_i(1), step_i(1), 0.d0, 0.d0
         write(3,*) np_i(2), 0.d0, step_i(2), 0.d0
         write(3,*) np_i(3), 0.d0, 0.d0, step_i(3)
         do i=1,natoms
            write(3,*) an(i), chrg(i), cartes(i,1), cartes(i,2), cartes(i,3)
         end do
         do i=1,ngrid
           write(3,*) I_vec(i)
         end do  
         close(3)
 end if        
! 40 format(6(E16.6E3))
 call cpu_time(T4)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!Print output!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 open(unit=3,file=outname) !open output file 
 write(3,*) "*******Job Finished************"
 write(3,*) "-----------Computational time----------------"
 write(3,*) "Total CPU time", T4-T1 
 write(3,*) "Grid points computing time", T2-T1
 write(3,*) "Primitive quartet loop time", T3-T2
 write(3,*) "CPU time for intracule integrations", T4-T3
 write(3,*) "---Primitive quartet time analysis-----------"
 write(3,*) "Reading .dm2p file-->", Tread
 write(3,*) "1st primitive screening-->", T1screen
 write(3,*) "2nd primitive screening-->", T2screen
 write(3,*) "Grid points loop-->", Tgrid
 write(3,*) "---------------------------------------------"
 write(3,*) "------Grid point + primitive quartet info----"
 write(3,*) "Original grid points", maxgrid
 write(3,*) "Symmetry reduced grid points", rgrid
 write(3,*) "Total number of reduced grid points", rrgrid  
 write(3,*) "---------------------------------------------"
 write(3,*) "--------------Accuracy check-----------------"
 write(3,*) "Sum of all DM2prim terms=", trDM2, trace_DM2prim
 write(3,*) "Thresholds used for DM2prim", trsh1, trsh2
 write(3,*) "Threshold used for screenings(Tau)=", thresh
 write(3,*) "Computed limit for the screenings", lim
 write(3,*) "---------------------------------------------"
 if (radial_integral) then   
    write(3,*) "--------------RADIAL INTEGRAL----------------" 
    write(3,*) "Number of centres", nquad
    do i=1,nquad
         write(3,*) i, ":::", cent(:,i)
    end do
    write(3,*) "Weights for centres", Ps(:)
    write(3,*) "Alpha parameter", sfalpha(:)
    write(3,*) "Gauss-Legendre nodes", nradc(:)
    write(3,*) "Gauss-Lebedev nodes", nangc(:)
    write(3,*) "---------------------------------------------"
    write(3,*) "---------------Final result------------------"
    write(3,*) "Intrac_integration=", r_integral
    npairs=(nelec)*(nelec-1)/2
    write(3,*) "Radial_integral error=", r_integral-dble(npairs)
 else if (radial_plot) then
     write(3,*) "RADIAL PLOT, I(S) vs s"    
     write(3,*) "Number of distances", nradi
     write(3,*) "From", radi(1), "to", radi(nradi)
 else if (cubeintra) then
     write(3,*) "CUBEFILE GENERATED"
 end if
 close(3)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function dfact(a) !computes double factorial
integer :: dfact
integer, intent(in) :: a
dfact=1
if (mod(a,2).eq.0) then !n is even
  do i=1,a/2
      dfact=dfact*(2*i)
  end do
else 
  do i=1,(a+1)/2
      dfact=dfact*(2*i-1)
  end do
end if
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


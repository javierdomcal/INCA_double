!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine intracule(normalize_dm2p)  !Computes vector or radial intracule, or radial integration
                        !need .wfx and .dm2p as input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use geninfo !information about the primitive functions
use intrastuff !subroutines and functions to compute the intracule 
use wfxinfo
use intrainfo !information from the input
use quadratures !nodes and weights of the quadratures, +total grid points
use output_manager  ! Unified output system
implicit none
logical, intent(in) :: normalize_dm2p

! Variables for unified output system
type(output_config) :: config
type(data_container) :: data
integer :: output_status

!!!!!!!!!!!!!!!local variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: kk1, kk2 !integers we do not want to read from .dm2
integer :: dfact    !double factorial
integer :: ii,iii,sum, ig, ir, summ, np, sm, smm !some integers
double precision :: pot, pot_x, pot_y, pot_z !power for the polynomial V 
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
double precision :: T1, T2, T3, T4, TT1, TT2, TT3, TT4, TT5 !time check
double precision :: Tread, T1screen, T2screen, Tgrid
double precision :: vee, h
double precision, allocatable :: vee_intra(:)
double precision :: intracule_total, intracule_zero

lim=thresh*(dble(nprim)*(dble(nprim)+1.d0)*0.5d0)**(-1.d0) !limit for the 1st integral screening
 
call cpu_time(T1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!Obtain grid points!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (radial_integral) then
    call gridpoints(nradc,nAngc,sfalpha,nquad,cent,Ps) !obtain the grid points
    allocate(r(3,rrgrid))
    r=rrrg
    ngrid=rrgrid
    deallocate(rrrg)
else if (radial_plot .or. vee_flag) then
     call gridpoints2(nblock,tart,stp,n_an_per_part)
     allocate(r(3,rgrid))
     r=rpg
     ngrid=rgrid
else if (cubeintra) then !vectorial plot
     call gridpoints3(center_i,step_i,np_i)    
     allocate(r(3,rgrid))
     r=rg
     ngrid=rgrid
else if (intracule_at_zero) then
     ngrid=1
     allocate(r(3,ngrid))
     r(:,1)=0.d0
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
intracule_zero=0.d0
 do while (.true.)  !loop for primitive quartets.
          call cpu_time(TT1)    
          read(5,end=100) kk1,i,j,k,l,DMval,kk2 !read a line from binary file .dm2
          if (i.eq.0) goto 100 !file is finished            
          call cpu_time(TT2) 
          Tread=Tread+(TT2-TT1)
          trDM2=trDM2+DMval
          smm=smm+1
          if (normalize_dm2p) then
              !normalization of DM2--> Normalize the primitives
              N_prim_i=(2.d0*Alpha(i)/pi)**(0.75d0)&
              *dsqrt(((4.d0*Alpha(i))**(dble(TMN(i,1)+TMN(i,2)+TMN(i,3))))*&
              dble(dfact(2*TMN(i,1)-1)*dfact(2*TMN(i,2)-1)*dfact(2*TMN(i,3)-1))**(-1.d0))  
   
              N_prim_j=(2.d0*Alpha(j)/pi)**(0.75d0)&
              *dsqrt(((4.d0*Alpha(j))**(dble(TMN(j,1)+TMN(j,2)+TMN(j,3))))*&
              dble(dfact(2*TMN(j,1)-1)*dfact(2*TMN(j,2)-1)*dfact(2*TMN(j,3)-1))**(-1.d0)) 
   
              N_prim_k=(2.d0*Alpha(k)/pi)**(0.75d0)&
              *dsqrt(((4.d0*Alpha(k))**(dble(TMN(k,1)+TMN(k,2)+TMN(k,3))))*&
              dble(dfact(2*TMN(k,1)-1)*dfact(2*TMN(k,2)-1)*dfact(2*TMN(k,3)-1))**(-1.d0)) 
         
              N_prim_l=(2.d0*Alpha(l)/pi)**(0.75d0)&
              *dsqrt(((4.d0*Alpha(l))**(dble(TMN(l,1)+TMN(l,2)+TMN(l,3))))*&
              dble(dfact(2*TMN(l,1)-1)*dfact(2*TMN(l,2)-1)*dfact(2*TMN(l,3)-1))**(-1.d0)) 
   
               !compute DMval with the normalization of primitives
               n_prim_t=N_prim_i*N_prim_j*N_prim_k*N_prim_l
               DMval=n_prim_t*DMval 
          end if     
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
         screen1=dabs(DMval)*dsqrt(J_ik(i,k,R_i_K_2)*J_jl(j,l,R_j_l_2))       
         
         call cpu_time(TT3)
         T1screen=T1screen-(TT3-TT2)         
         if (screen1.ge.lim) then    
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 
                 !compute the other variables (eqn. 12)
                 a_ijkl=a_ik+a_jl
                 e_ijkl=(a_ik*a_jl)*a_ijkl**(-1.d0)
                 sqe=dsqrt(e_ijkl)
                 do ii=1,3                   !X_ik,Y_ik,Z_ik,...
                       R_ik(ii)=(Alpha(i)*Cartes(Ra(i),ii)+Alpha(k)*Cartes(Ra(k),ii))*a_ik**(-1.d0)
                       R_jl(ii)=(Alpha(j)*Cartes(Ra(j),ii)+Alpha(l)*Cartes(Ra(l),ii))*a_jl**(-1.d0)
                       R_ijkl(ii)=(a_ik*R_ik(ii)+a_jl*R_jl(ii))*a_ijkl**(-1.d0)
                 end do                      
                 Alf_ijkl=0.5d0*(a_ik-a_jl)*a_ijkl**(-1.d0)

                 !compute Aa_ijkl (the grid independent part)  !eq.18
                 A_ind=(a_ijkl)**(-1.5d0)* dexp(-e_ik*R_i_k_2-e_jl*R_j_l_2)
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
                 screen2=dabs(Aa_ijkl*U_x*U_y*U_z) !eqn.40
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
                                                if (dble(pot).gt.1.d-16) then
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
                                   I_vec(ig)=I_vec(ig)+Aa_ijkl*dexp(-(rp(1)**2.d0+rp(2)**2.d0+rp(3)**2.d0))*V_x*V_y*V_z
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
 close(5)
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!Compute the integrals using the quadrature !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!weights and the vectorial intracule        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call cpu_time(T3)

! Output radial plot using unified system
if (radial_plot) then        
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
            r_intra(i)=(r_intra(i)/(4.0d0*pi))
      end do 
      
      ! Use unified output system for radial plot
      config = create_output_config('Intracule_Radial_Plot', OUTPUT_VECTOR)
      config%description = 'Radial intracule I(r) vs distance'
      config%units = 'atomic units'
      config%method = 'Cioslowski-Liu algorithm with angular integration'
      config%computation_time = T3-T1
      config%notes = 'Electron-electron radial distribution function'
      allocate(config%column_labels(2))
      config%column_labels(1) = 'Distance'
      config%column_labels(2) = 'Intracule'
      
      data%vector_size = nradi
      allocate(data%vector_data(nradi), data%vector_coords(nradi))
      data%vector_data = r_intra
      data%vector_coords = radi
      
      call write_data(config, data, output_status)
      
      ! Cleanup
      deallocate(data%vector_data, data%vector_coords)
      deallocate(config%column_labels)
      deallocate(r_intra)
end if     

! Output Vee calculation using unified system
if (vee_flag) then        
      allocate(r_intra(nradi))
      allocate(vee_intra(nradi))
      
      ig = 0
      sm = 0
      r_intra = 0.0d0
      vee_intra = 0.0d0
      ir = 0
      
      ! Perform angular integration for each radial point
      do i = 1, nradi
            do k = 1, smn(i) ! sum reduced angular points per radial shell
                sm = sm + 1
                r_intra(i) = r_intra(i) + w_ang(sm) * I_vec(sm) ! Perform the angular quadrature 
                vee_intra(i) = vee_intra(i) + w_ang(sm) * I_vec(sm) ! Same angular quadrature
            end do
            
            ! Handle special case at origin
            if (radi(i) .le. 0.0d0) then
                intracule_zero = r_intra(1)
                r_intra(i) = 0.0d0
                vee_intra(i) = 0.0d0
            else
                ! First normalization: without /radi(i) at the end
                r_intra(i) = (r_intra(i)/(4.0d0*pi)) * 4.0d0*pi*radi(i)**2
                
                ! Second normalization: with /radi(i) at the end (original definition)
                vee_intra(i) = (vee_intra(i)/(4.0d0*pi)) * 4.0d0*pi*radi(i)**2/radi(i)
            end if
            
      end do
      
      ! Compute integrations using Simpson's 1/3 rule
      h = radi(2)-radi(1)  ! spacing = 0.1 au
      
      ! First integration for r_intra (intracule)
      intracule_total = simpson_integrate(r_intra, nradi, h)
      
      ! Second integration for vee_intra (Vee)
      vee = simpson_integrate(vee_intra, nradi, h)
     
      ! Output all Vee results through unified system
      
      ! Output intracule total
      config = create_output_config('Intracule_Total', OUTPUT_SCALAR)
      config%description = 'Total integrated intracule value'
      config%units = 'electrons'
      config%method = 'Simpson rule integration'
      config%computation_time = T3-T1
      data%scalar_value = intracule_total
      call write_data(config, data, output_status)
      
      ! Output Vee total  
      config = create_output_config('Vee_Total', OUTPUT_SCALAR)
      config%description = 'Total electron-electron repulsion energy'
      config%units = 'hartree'
      config%method = 'Simpson rule integration of Vee density'
      config%computation_time = T3-T1
      data%scalar_value = vee
      call write_data(config, data, output_status)
      
      ! Output intracule at zero
      config = create_output_config('Intracule_At_Zero', OUTPUT_SCALAR)
      config%description = 'Intracule value at zero electron separation'
      config%units = 'electrons'
      config%method = 'Direct evaluation at coalescence'
      config%computation_time = T3-T1
      data%scalar_value = intracule_zero
      call write_data(config, data, output_status)
      
      ! Output Vee radial distribution
      config = create_output_config('Vee_Radial_Distribution', OUTPUT_VECTOR)
      config%description = 'Electron-electron repulsion energy density vs distance'
      config%units = 'atomic units'
      config%method = 'Normalized intracule with radial weighting'
      config%computation_time = T3-T1
      allocate(config%column_labels(2))
      config%column_labels(1) = 'Distance'
      config%column_labels(2) = 'Vee_Density'
      
      data%vector_size = nradi
      allocate(data%vector_data(nradi), data%vector_coords(nradi))
      data%vector_data = vee_intra
      data%vector_coords = radi
      
      call write_data(config, data, output_status)
      
      ! Cleanup
      deallocate(data%vector_data, data%vector_coords)
      deallocate(config%column_labels)
      deallocate(r_intra)
      deallocate(vee_intra)
end if 

! Output radial integral using unified system
if (radial_integral) then 
     r_integral=0.d0 
     vee=0.d0
     do i=1,rrgrid
         r_integral=r_integral+rweight(i)*I_vec(i)
         vee=vee+rweight_vee(i)*I_vec(i)
     end do 
     deallocate(rweight)
     
     ! Output radial integral as scalar value
     config = create_output_config('Radial_Integral', OUTPUT_SCALAR)
     config%description = 'Total integrated intracule from 0 to infinity'
     config%units = 'electrons'
     config%method = 'Gauss-Legendre quadrature integration'
     config%computation_time = T3-T1
     config%notes = 'Should equal number of electron pairs for exact wavefunction'
     
     data%scalar_value = r_integral
     call write_data(config, data, output_status)
     
     ! Output Vee integral
     config = create_output_config('Vee_Integral', OUTPUT_SCALAR)
     config%description = 'Total electron-electron repulsion energy'
     config%units = 'hartree'
     config%method = 'Gauss-Legendre quadrature integration'
     config%computation_time = T3-T1
     
     data%scalar_value = vee
     call write_data(config, data, output_status)
     
     ! Output intracule at zero
     config = create_output_config('Intracule_At_Zero', OUTPUT_SCALAR)
     config%description = 'Intracule value at zero electron separation'
     config%units = 'electrons'
     config%method = 'Direct evaluation'
     config%computation_time = T3-T1
     
     data%scalar_value = intracule_zero
     call write_data(config, data, output_status)
end if

! Output cube file using unified system
if (cubeintra) then
     ! Use unified output system for cube file
     config = create_output_config('Intracule_3D', OUTPUT_CUBE)
     config%description = 'Three-dimensional intracule distribution'
     config%units = 'electrons per unit volume'
     config%method = 'Cioslowski-Liu algorithm on regular 3D grid'
     config%computation_time = T3-T1
     config%origin = center_i
     config%spacing = step_i
     config%dimensions = np_i
     
     data%grid_size = ngrid
     allocate(data%grid_values(ngrid))
     ! Scale the values as in original code
     data%grid_values = I_vec/2.0d0
     
     call write_data(config, data, output_status)
     
     ! Cleanup
     deallocate(data%grid_values)
end if

! Output intracule at zero using unified system
if (intracule_at_zero) then
     config = create_output_config('Intracule_Zero_Distance', OUTPUT_SCALAR)
     config%description = 'Intracule value at r12 = 0 (electron coalescence)'
     config%units = 'electrons'
     config%method = 'Direct evaluation at zero separation'
     config%computation_time = T3-T1
     config%notes = 'Electron coalescence probability density'
     
     data%scalar_value = intracule_zero
     call write_data(config, data, output_status)
end if

call cpu_time(T4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!Output comprehensive summary using unified system!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Create comprehensive calculation summary using TABLE format
config = create_output_config('Calculation_Summary', OUTPUT_TABLE)
config%description = 'Complete computational statistics and timing analysis'
config%method = 'Cioslowski-Liu algorithm with integral screening'
config%computation_time = T4-T1
config%notes = 'Performance analysis and accuracy validation'

! Prepare summary data table
data%table_rows = 15
data%table_cols = 2
allocate(data%table_data(data%table_rows, data%table_cols))
allocate(config%column_labels(2))
config%column_labels(1) = 'Metric'
config%column_labels(2) = 'Value'

! Fill timing and performance statistics
data%table_data(1, 1) = 1.0d0   ! Total CPU time
data%table_data(1, 2) = T4-T1
data%table_data(2, 1) = 2.0d0   ! Grid setup time
data%table_data(2, 2) = T2-T1
data%table_data(3, 1) = 3.0d0   ! Quartet loop time
data%table_data(3, 2) = T3-T2
data%table_data(4, 1) = 4.0d0   ! Integration time
data%table_data(4, 2) = T4-T3
data%table_data(5, 1) = 5.0d0   ! File reading time
data%table_data(5, 2) = Tread
data%table_data(6, 1) = 6.0d0   ! First screening time
data%table_data(6, 2) = T1screen
data%table_data(7, 1) = 7.0d0   ! Second screening time
data%table_data(7, 2) = T2screen
data%table_data(8, 1) = 8.0d0   ! Grid evaluation time
data%table_data(8, 2) = Tgrid
data%table_data(9, 1) = 9.0d0   ! Grid points used
data%table_data(9, 2) = real(ngrid, 8)
data%table_data(10, 1) = 10.0d0 ! DM2 trace
data%table_data(10, 2) = trDM2
data%table_data(11, 1) = 11.0d0 ! Normalized DM2 trace
data%table_data(11, 2) = trace_DM2prim
data%table_data(12, 1) = 12.0d0 ! Screening threshold
data%table_data(12, 2) = thresh
data%table_data(13, 1) = 13.0d0 ! Quartets skipped (1st)
data%table_data(13, 2) = real(sum, 8)
data%table_data(14, 1) = 14.0d0 ! Quartets skipped (2nd)  
data%table_data(14, 2) = real(summ, 8)
data%table_data(15, 1) = 15.0d0 ! Intracule at zero
data%table_data(15, 2) = intracule_zero

call write_data(config, data, output_status)

! Cleanup
deallocate(data%table_data)
deallocate(config%column_labels)

! Final cleanup
deallocate(I_vec)
if (allocated(r)) deallocate(r)

contains

! Simpson's 1/3 rule integration function
function simpson_integrate(y, n, h) result(integral)
    implicit none
    integer, intent(in) :: n           ! number of points
    real(8), intent(in) :: y(n)        ! function values
    real(8), intent(in) :: h           ! grid spacing
    real(8) :: integral                 ! result
    integer :: i
    
    if (n < 3) then
        integral = 0.0d0
        return
    end if
    
    integral = y(1) + y(n)  ! first and last points
    
    ! Add middle points with alternating weights
    do i = 2, n-1
        if (mod(i, 2) == 0) then
            integral = integral + 4.0d0 * y(i)  ! even indices
        else
            integral = integral + 2.0d0 * y(i)  ! odd indices
        end if
    end do
    
    integral = integral * h / 3.0d0
end function simpson_integrate

end subroutine intracule

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine filewfx(wfxfilename)                                             !
!  1-Reads wfxfile and stores the important values in wfxinfo module        !
!    (nelec, noccmo, nprim, ...)                                            !
!  2-Determines the wavefunction type we have: RHF, UHF, correlated WF      !
!  +relaxed/unrelaxed density                                               !
!  3-Builds TMN matrix (matrix with primitive exponents)                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use wfxinfo
use geninfo
use located
implicit none
character*40, intent(in) :: wfxfilename
character*80 :: line
character*80 :: zaborra  !coses que no volem llegir
integer :: i,j,k, kk
integer :: mon !molecular orbital number  
logical :: rhf, rdens, udens, opsh, clsh 
double precision :: maxim, minim
double precision, parameter :: trsh=1.d-6 !threeshold for orbital occupancies

 corr=.false.  
 uhf=.false.
 rhf=.false.
 udens=.false.
 rdens=.false.  
 opsh=.false.
 clsh=.false.
   
open(unit=1,file=wfxfilename,status='OLD') 

call locate(1,"Number of Nuclei") 
read(1,'(a80)')line ! llegeix el que hi ha a la linia que hi ha despres de la linia que conte "Number of Nuclei"
read(line(1:4),*) natoms ! del contingut de "line", guarda el que hi ha a les 3 primeres posicions en una variable
                        
rewind 1  

call locate(1,"Number of Occupied Molecular Orbitals")
read(1,'(a80)')line
read(line(1:5),*) noccmo

allocate(Occ(noccmo))

rewind 1

call locate(1,"Molecular Orbital Occupation Numbers") !check RHF,UHF, or Correlated
do i=1,noccmo                                         !by     
   read(1,*) Occ(i)                                   !reading the MO occupation values
end do
 
do i=1,noccmo
     if ((Occ(i).ne.1.d0).and.(Occ(i).ne.2.d0)) then
        write(*,*) "Correlated wavefunction"
        corr=.true.
        exit
     end if
end do

maxim=maxval(Occ)                                         !store the maximum occupation value
minim=minval(Occ)   
if (corr) then
    if (maxim.ge. 2.d0) then             !MOs filled by pairs
        if ((maxim-2.d0 .gt. trsh).or.(minim.lt. 0.d0)) then 
           write(*,*) "RELAXED DENSITY (closed shell)"                  
           rdens=.true.
        else     
           write(*,*) "UNRELAXED DENSITY (closed shell)"!all the occupation values between 0 and 2
           udens=.true.
        end if                    
    else                                       !1 electron MOs (Unrestricted)              
        if ((maxim-1.d0 .gt. trsh).or.(minim.lt. 0.d0)) then 
           write(*,*) "RELAXED DENSITY (open shell)"
           rdens=.true.
        else    
           write(*,*) "UNRELAXED DENSITY (open shell)"
           udens=.true.
        end if                    
    end if 
else                    !HF mwavefunction, check if it's RHF or UHF 
   if (minim.eq. 1.d0) then                 !we have 1e- in 1 orbital     
         if (maxim.eq. 2.d0) then           !at least 1 doubly occ orbital
            rhf=.true.                      !open shell RHF
            uhf=.false.
            write(*,*) "Open shell RHF"
         else                                 !all singly occ. orbitals (alpha and beta)
            rhf=.false.                        !UHF   
            uhf=.true. 
            write(*,*) "UHF"
         end if
   else if (minim.eq.2.d0) then           !all doubly occ. orbitals
         rhf=.true.                           !closed shell rhf
         uhf=.false.  
         write(*,*) "Closed shell RHF"
   else 
         write(*,*)   "Error, correlated wf"  
         STOP
   end if    
  deallocate(Occ)                        !with HF wavefunctions we no longer need the occupancy   
end if  
 
call locate(1,"Net Charge")           !store some variables
read(1,*) netch

rewind 1

call locate(1,"Number of Electrons")
read(1,'(a80)')line
read(line(1:5),*) nelec

call locate(1,"Number of Alpha Electrons")
read(1,'(a80)')line
read(line(1:5),*) nalfae

rewind 1

call locate(1,"Number of Beta Electrons")
read(1,'(a80)')line
read(line(1:5),*) nbetae

rewind 1

call locate(1,"Number of Primitives")
read(1,'(a80)')line
read(line(1:5),*) nprim

rewind 1

allocate(cartes(natoms,3))
call locate(1,"Nuclear Cartesian Coordinates")
read(1,*) ((cartes(i,j) , j=1,3 ), i=1,natoms)

rewind 1

allocate(Ra(nprim))
call locate(1,"Primitive Centers")
read(1,*) (Ra(i), i=1,nprim)

rewind 1

allocate(Ptyp(nprim))
call locate(1,"Primitive Types")
read(1,*) (Ptyp(i), i=1,nprim)

rewind 1

allocate(Alpha(nprim))
call locate(1,"Primitive Exponents")
read(1,*) (Alpha(i), i=1,nprim)

rewind 1

allocate(an(natoms))
call locate(1,"Atomic Numbers")
do i=1,natoms
 read(1,*) an(i)
end do

allocate(chrg(natoms))
call locate(1,"Nuclear Charges")
do i=1,natoms
  read(1,*) chrg(i)
end do  

allocate(T(noccmo,nprim))                  !read T matrix (MOs in prims)       
call locate(1,"Molecular Orbital Primitive Coefficients")
do kk=1,noccmo
  read(1,'(a80)')line
  read(1,*) mon                          !molecular orbital number "mo"
 read(1,'(a80)') zaborra  
 read(1,*) (T(mon,i), i=1,nprim)
end do

rewind 1

allocate(TMN(nprim,3))      !create a matrix with t m and n exponents for each primitive (see prim equation)
do i=1,nprim                                                                 !t->lx
   if (ptyp(i).eq.1) then   !s prim                                          !m->ly
     do j=1,3                 						!n->lz
       TMN(i,j)=0
     end do
   else if (ptyp(i).eq.2) then  !px
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=0
   else if (ptyp(i).eq.3) then  !py
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=0   
   else if (ptyp(i).eq.4) then  !pz
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=1   
   else if (ptyp(i).eq.5) then  !dx²
      TMN(i,1)=2
      TMN(i,2)=0
      TMN(i,3)=0   
   else if (ptyp(i).eq.6) then  !dy²
      TMN(i,1)=0
      TMN(i,2)=2
      TMN(i,3)=0  
   else if (ptyp(i).eq.7) then  !dz²
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=2   
   else if (ptyp(i).eq.8) then  !dxy
      TMN(i,1)=1
      TMN(i,2)=1
      TMN(i,3)=0   
   else if (ptyp(i).eq.9) then  !dxz
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=1
   else if (ptyp(i).eq.10) then !dyz 
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=1      
   else 
     write(*,*) "Error, f orbital!"
     stop
   end if
end do
                             !if we have an unrestricted wf-------What about correlated???
if (uhf) then                !split T matrix into T_a and T_b (alpha and beta)
     allocate(T_a(nalfae,nprim))
     allocate(T_b(nbetae,nprim))
     do i=1,nalfae                !for alpha molecular orbitals
          do j=1,nprim    
                T_a(i,j)= T(i,j)  !store the prim. coefficients of alpha orb.   
          end do
     end do
     k=nalfae+1                   !first beta orbital
     do i=1,nbetae     
          do j=1,nprim
                T_b(i,j) = T(k,j)     
          end do
          k=k+1                   !next beta orbital
     end do 
end if 
     
 close(1) 
end subroutine filewfx

!****************************************************************************
subroutine filelog(logfilename)                                             !
!1-reads d (Flg) fixed coeficients of the primitives			       !
!2-Normalizes that coeficient (with respect to the primitives, depending on !
!  primitive type (ptyp)						       !
!3-Builds C matrix (expansion of MOs into AOs)                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use loginfo                    
use geninfo
use wfxinfo
use located
implicit none
integer, parameter :: maxao=100
character*10 :: logfilename
character*80 :: line 
character*80 :: zaborra
character*10 :: tao !atomic orbital type
integer :: i,j,k,l,m
integer :: sm, smm, sm2, sm3, sm4 !sums !smp is the sum of the previous primitives
integer :: smatom
real :: kk !alpha value of the primitives (already got it from wfx) 

open(unit=4 ,file=logfilename, status='OLD') !obre arxiu amb nom name a la unitat 4

allocate(npao(maxao)) !number of primitives of each AO basis. 
allocate(Flg(nprim))  !d coeficient for each primitive
allocate(aotyp(maxao))!atomic orbital type

Flg=0.d0

smm=0 !sum primitives
sm=0  !sum hybrid AOs (S, SP, D, ...)
smatom=1
call locate(4,"AO basis set in the form of general basis input") 
read(4,'(a80)') zaborra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!read fixed coeficients of the AOs in the primitives!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do while (smatom.le.natoms)   
      read(4,'(a80)') line
      read(line(2:3),*) tao 
      if (tao.eq."**") then !new atom finished 
           smatom=smatom+1
           read(4,'(a80)') zaborra   
      else if (tao.eq."S") then
           sm=sm+1  !one more AO
           read(line(5:7),*) npao(sm) !number of primitives of S ao   
           do i=1,npao(sm)
                smm=smm+1 !one more primitive
                read(4,'(a80)') line
                read(line(24:39),*) Flg(smm) !read d coeficient of smm primitive  
           end do
           aotyp(sm)=0  !S type ao
       else if (tao.eq."SP") then 
            sm=sm+1 !one more ao, S type
            aotyp(sm)=0
            read(line(5:7),*) npao(sm)      
            do i=1,npao(sm)
                 smm=smm+1  !one more primitive
                 sm2=sm+1     !px AO
                 aotyp(sm2)=1
                 sm3=sm+2     !py AO
                 aotyp(sm3)=1
                 sm4=sm+3     !pz AO  
                 aotyp(sm4)=1  
                 k=smm+npao(sm) !k is the p type primitives, smm is for the S type primitives     
                 read(4,*) kk, Flg(smm), Flg(k) !read d coeficients for S and for P primitives
                 l=k+npao(sm) !py primitive number                  
                 m=l+npao(sm) !pz primitive number    
                 Flg(l)=Flg(k)     !same value(C_px1=C_py1=C_pz1, 
                                   !(C_px2=C_py2=C_pz2, ...)  
                 Flg(m)=Flg(k)     
            end do   
            smm=smm+3*npao(sm)     !sum all the primitives used in this P type orbital                
            do i=1,3
                 npao(sm+i)=npao(sm) !In an SP hybrid orbital, px, py and pz AOs have the same number of primitives than S
            end do
            sm=sm+3 !sum px, py and pz to the total number of AOs
       else if (tao.eq."D") then
            sm=sm+1  !one more ao
            aotyp=2 !dx, dy, dz  !aotyp is used to the normalization of d coeficient (see norm subroutine)
            read(line(5:7),*) npao(sm)   
            do i=1,npao(sm) !do for the number of primitives of this AO
                 smm=smm+1  
                 read(4,*) kk, Flg(smm)      
                 !check if we have cartesian or spherical coordinates    
                 do j=1,5  !six AOs with same npao and F coeficient (if cartesian)
                        k=smm+j
                        l=sm+j  !current AO 
                        if (l.le.3) then 
                                 aotyp(l)=2  !dx,dy or dz coordinate
                        else
                                 aotyp(l)=3  !square (different normalization)  
                        end if 
                        Flg(k)=Flg(smm)      !same coeficient for each primitive of this AO 
                        npao(sm+j)=npao(sm)    
                 end do
            end do    
            smm=smm+6*npao(sm)  !sum all the primitives and AOs used. 
            sm=sm+5     
       end if
end do
nao=sm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!Normalize for each d coeficient of each primitive!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(N_prim(nprim))  
do i=1,nprim
      if (ptyp(i).eq.1) then 
            N_prim(i)=(pi*(2*Alpha(i))**(-1))**(-3.d0/4.d0)
           write(*,*) N_prim(i)
      else if ((ptyp(i).gt.1).and.(ptyp(i).le.4)) then
            N_prim(i)= (pi*0.5d0)**(-3.d0/4.d0) * (2*(Alpha(i)**(1.25d0)))   
      else if ((ptyp(i).gt.4).and.(ptyp(i).le.7)) then  !dxy,dxz.dyz
            N_prim(i)=(pi*0.5d0)**(-3.d0/4.d0) * (2**2 * Alpha(i)**(1.75d0))
      else if ((ptyp(i).gt.7).and.(ptyp(i).le.10)) then  !dx²,dy².dz²
            N_prim(i)=(pi*0.5d0)**(-3.d0/4.d0) * (2**2 * Alpha(i)**(1.75d0))/(sqrt(3.d0))   
      else
            write(*,*) "f type orbital, insert its normalization equation"
      end if    
end do
do i=1, nprim              !Multiply the normalization of the primitives
      Flg(i)=Flg(i)*N_prim(i)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!Create expansion coef. of MOs in AOs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cmatrix() 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!102 format(30(F4.2,1X))
!103 format(30(I4,1X)) 
 close(4)
end subroutine filelog

!*******************************************************************************
subroutine cmatrix()   
!!calculates C with F and T matrices. C are the contracted coefficients of AOs in
!MOs. F are coefficients of AOs in primitives and T are for MOs in Primirives.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use wfxinfo
use loginfo
use geninfo
implicit none
integer :: i, j, k
     
allocate(Ckalk(noccmo,nao))
write(*,*) nao
write(*,*) noccmo
do i=1,noccmo
 k=1
 do j=1,nao
      Ckalk(i,j)=T(i,k)/(Flg(k)) 
      k=k+npao(j)  
  end do
end do   

! 41 format(10(F7.3,1X))
end subroutine cmatrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program wavefunction !main program
!!!!!reads .wfx and .log files and performs calculations about this info
use inputdat  !information about the calculations we want to do
implicit none
integer :: a !defines to subroutine cubefile what function we want to represent: prim, ao, mo, dens
!double precision :: gradx, grady, gradz !gradient
logical :: normalize_dm2p

readwfx=.false. 
readlog=.false. 

call readinput()  

if (readwfx) call filewfx(wfxfilename)  !reads info from a wfx file 
if (readfchk) call filefchk(fchkfilename)
if (readlog) call filelog(logfilename)  !reads info from a log file

if (cube) then
  if (primcube) then  
      a=1
      call cubefile(a,nameprim) !Generate cubefile with a Primitive
  end if

  if (aocube) then 
      a=2
      call cubefile(a,nameao)  !Generate cubefile with AO
  end if

  if (mocube) then 
      a=3
      call cubefile(a,namemo) !Generate cubefile with a MO 
  end if  

  if (denscube) then
      a=5
      call cubefile(a,namedens) !Generate cubefile with density from MO
  end if

  !if (gradient) then 
  !  call gradient(0,0,0,gradx,grady,gradz) !not impletented
  !end if

  if (laplacian) then 
     a=6     
     call cubefile(a,namelap) !generate a cubefile with laplacian
  end if
end if

if (intracalc) then !compute the intracule
  call readintra()
  write(*,*) "readed intra info"
  normalize_dm2p=.true.
  call intracule(normalize_dm2p)  
  write(*,*) "intracule computed succesfully"
end if

if (c1calc) then
   call c1hole(70,1.d0)
end if
end program wavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   In this subroutine we read an input file where we specify what type of
!   calculations we want to do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readinput() 
use inputdat
use cubeinfo

use locatemod
implicit none
character*80 :: name
integer :: i
call getarg(1,name)  !gets the name of the input file (writen in the prompt)
name=trim(name)      !remove the blank spaces of the string 'name'
open(unit=3,file=name,status='OLD') 
 
!set up logical variables as false
readwfx=.false. 
readfchk=.false. 
readlog=.false.
primcube=.false.
aocube=.false.
MOcube=.false.
denscube=.false.
gradient=.false.
laplacian=.false.
intracalc=.false.


call locate(3,"$wfxfile")
read(3,*) wfxfilename
if (wfxfilename.ne.'no') then
      readwfx=.true.
end if
rewind 3
call locate(3,"$fchkfile")
read(3,*) fchkfilename
if (fchkfilename.ne.'no') then
      readfchk=.true.
end if
rewind 3

call locate(3,"$logfile")
 read(3,*) logfilename
 if (logfilename.ne.'no') then
      readlog=.true.
 end if 

 rewind 3

 call locate(3,"$cubefile")
 read(3,*) cube
 if (cube) then 
   read(3,*) (center(i), i=1,3) !cube centered in (x,y,z)
   read(3,*) (step(i), i=1,3) !distance between points in the axis
   read(3,*) (np(i), i=1,3) !number of points for each axis
  
   rewind 3

   call locate(3,"$primitive")
   read(3,*) nameprim
   if (nameprim.ne.'no') primcube=.true.
   read(3,*) npr !number of primitive we want to print

   rewind 3

   if (readlog) then    
      call locate(3,"$AO")
      read(3,*) nameao
   if (nameao.ne.'no') aocube=.true. 
      read(3,*) cao  !number of ao we want to print
   end if 

   rewind 3

   call locate(3,"$MO")
   read(3,*) namemo
   if (namemo.ne.'no') MOcube=.true.    
   read(3,*) mo  !number of MO we want to print

   rewind 3

   call locate(3,"$density")
   read(3,*) namedens 
   if (namedens.ne.'no') denscube=.true.
       
   rewind 3

   call locate(3,"$gradient")
   read(3,*) namegrad
   if (namegrad.ne.'no') gradient=.true.

   rewind 3

   call locate(3,"$laplacian")
   read(3,*) namelap
   if (namelap.ne.'no') laplacian=.true.
 
   rewind(3)
   close(3)
   
 end if
 
 call locate(3,"$intracule")
 read(3,*) nameint
 if (nameint.ne.'no') then
        intracalc=.true.
 end if

 call locate(3,"$HFhole")
 read(3,*) c1calc
end subroutine readinput 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readintra() 
use intrainfo
use locatemod
implicit none
character*80 :: name
integer :: i, j
character*20 :: ccent
logical :: beckw
call getarg(1,name)  !gets the name of the input file (writen in the prompt)
name=trim(name)      !remove the blank spaces of the string 'name'
open(unit=3,file=name,status='OLD')  
        !read input for intracule
        call locate(3,'$Threshold') !threshold for integral screenings
        read(3,*) thresh  
        rewind(3)
        dm2name = ''
        call locate(3,'$DM2')
        read(3,*) dm2name           !name of the dm2p file
        read(3,*) trsh1, trsh2      !thresholds used in DM2prim
        read(3,*) outname           !name of the output file
        rewind(3)
        
        call locate(3,'$radial_integral')
        read(3,*) radial_integral   
        read(3,*) dif_nodes
        if (radial_integral) then   !read info to perform radial integral
            call locate(3,'$Number of quadrature')
            read(3,*) ccent
            autoc=.false.
            if (ccent.eq."manual") then  !input parameters
                read(3,*) nquad          !number of quadrature centers      
                call locate(3,'Quadrature center')
                allocate(cent(3,nquad))
                do i=1,nquad            !read position of centers
                   read(3,*) (cent(j,i), j=1,3)
                end do
            else !default
                write(*,*) "Computing centers automatically"
                call centercalc()
                write(*,*) "nquad computed =", nquad
                !call pdint() !calculate approximate I_vs_r curve
            end if            
            if (dif_nodes) then       !Different node for each centre
                allocate(nradc(nquad))
                allocate(nangc(nquad))
                allocate(sfalpha(nquad)) !PROBLEM: We can't know nquad beforehand
                call locate(3,'Gauss-Legendre')
                read(3,*) (nradc(i), i=1,nquad)
                read(3,*) sfalpha(:)
                rewind(3)           
                call locate(3,'$Gauss-Lebedev')
                read(3,*) nangc(:)
            else        
                call locate(3,'Gauss-Legendre')
                read(3,*) nrad
                read(3,*) sfalpha
                if (nquad.eq.1) then  !Becke isn't used 
                   read(3,*) a,b       !Choose limits of the radial integral
                end if                 !if 0 the integral is computed from 0 to infty
                rewind(3)
                call locate(3,'$Gauss-Lebedev')
                read(3,*) gpt !write number of angular nodes
                do i=1,nquad
                    nradc(i)=nrad
                    nangc(i)=gpt
                    sfalpha(i)=sfalpha(1)
                end do                
            end if 
            rewind(3)
            call locate(3,'$Center weights')
            allocate(Ps(nquad)) !allocate weight of each centre
            !read(3,*) beckw
            !if (beckw) then
                !call pdint()
                !compute Becke centres weights automatically (TO IMPLEMENT)
            !else
            read(3,*) Ps(:) !the weigth of a positive center must be equal to the neg.
            !end if
            rewind(3)
        end if   
        call locate(3,'$radial_plot')
        read(3,*) radial_plot
        if (radial_plot) then
                read(3,*) r_plot_name
                read(3,*) nblock
                allocate(n_an_per_part(nblock))
                allocate(tart(2,nblock))
                allocate(stp(nblock))
                do i=1,nblock
                  read(3,*) tart(:,i), stp(i), n_an_per_part(i)
                end do
        end if
        rewind(3)
        call locate (3,'$Vectorial_plot')
        read(3,*) cubeintra
        if (cubeintra) then
                read(3,*) cubeintraname
                read(3,*) (center_i(i), i=1,3) !cube centered in (x,y,z)
                read(3,*) (step_i(i), i=1,3) !distance between points in the axis
                read(3,*) (np_i(i), i=1,3) !number of points for each axis
        end if
        rewind(3)
        call locate (3,'$Vee')
        write(*,*) 'yeehaw'
        read(3,*) vee_flag
        if (vee_flag) then
                read(3,*) nblock
                allocate(n_an_per_part(nblock))
                allocate(tart(2,nblock))
                allocate(stp(nblock))
                do i=1,nblock
                  read(3,*) tart(:,i), stp(i), n_an_per_part(i)
                end do
        end if
        rewind(3)
 close(3) 
end subroutine readintra

subroutine autocac()
use geninfo
use radis
implicit none
!character*20 :: param
integer :: i
double precision, dimension(natoms) :: radii
double precision :: maxrad

do i=1,natoms
   radii(i)=BL(an(i)) !assing a radius to each atom of the system. 
end do
maxrad=maxval(radii)


end subroutine autocac






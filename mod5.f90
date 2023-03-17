module wfxinfo !especific data of wfx files 
  implicit none
  logical :: corr   !correlated method
  integer :: noccmo !number of occupied MOs
  double precision, allocatable, dimension(:,:) :: T   !MOs in primitives coeficients (matriu T)
  double precision, allocatable, dimension(:,:) :: T_a !Alpha MOs coeficients
  double precision, allocatable, dimension(:,:) :: T_b !Beta MOs coeficients
  double precision, allocatable, dimension(:) :: Occ   !ocupancies of orbitals
end module wfxinfo
!---------------------------------------------------------------------------------------
module geninfo           !general information, can be obtained from .wfx files.
  implicit none
  logical :: uhf, rhf    !unrestricted-restricted Hartree Fock???
  logical :: rdens,udens !relaxed or unrelaxed density
  logical :: opsh, clsh  !open or closed shell
  integer :: natoms      !number of atoms
  integer :: nelec
  integer :: nprim       !number of primitives
  double precision :: netch !net charge
  integer :: mspin       !electronic spin multiplicity
  integer :: nalfae      !number of alpha electrons
  integer :: nbetae      !number of beta electrons
  integer, allocatable, dimension(:) :: Ra    !atomic centers of primitives
  integer, allocatable, dimension(:) :: Ptyp  !primitive type
  double precision, allocatable, dimension(:) :: Alpha !Primitive exponents
  integer, allocatable, dimension(:,:) :: TMN !matrix with t,m,n coeficients
  double precision, allocatable,  dimension(:,:) :: cartes !nuclear cartesian coordinates
  integer, allocatable, dimension(:) :: an !atomic number
  double precision, allocatable, dimension(:) :: chrg !nuclear charge  
  double precision, parameter :: pi=4.d0*datan(1.d0)
  
end module
!----------------------------------------------------------------------------------------
 module loginfo
   implicit none 
   double precision, allocatable, dimension(:) :: Flg, N_prim !fixed expansion coeficient for the contraction of AOs into primitives, normalization constant for the primitives
   double precision, allocatable, dimension(:,:) :: Ckalk !expansion coeficient of MOs in AOs
   integer, allocatable, dimension(:) :: npao !number of primitives for each atomic orbital
   integer :: nao !number of atomic orbitals in the molecule
   integer, allocatable, dimension(:) :: aotyp !type of atomic orbital: s->0, px, py, pz ->1,  ...   
end module loginfo   

!------------------------------------------------------------------------------------------

module inputdat
implicit none
 logical :: readwfx, readlog, cube, primcube, aocube, MOcube, denscube, gradient, laplacian, intracalc
 logical :: c1calc
 character*40 :: wfxfilename
 character*40 :: logfilename
 character*40 :: nameprim
 character*40 :: nameao
 character*40 :: namemo   !name of the cubefile we will generate
 character*40 :: namedens
 character*40 :: namegrad
 character*40 :: namelap
 character*40 :: nameint
 character*40 :: namec1
end module inputdat

!----------------------------------------------------------------------------------------------------
module intrainfo

!input variables for intracule calculation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
character*20 :: dm2name, outname
double precision :: thresh
double precision :: trsh1, trsh2 !thresholds for DM2prim
!!!!!Radial integral parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical :: autoc                                       !calculate centres automatically
logical :: radial_integral                    !=.true. integral of radial intra from 0 to infty
integer :: nquad !number of quadratures
double precision, allocatable, dimension(:,:) :: cent  !integration centre
double precision, allocatable, dimension(:) :: sfalpha !scaling factor for Gauss-Legendre
double precision, allocatable, dimension(:) :: Ps !Bragg-Slater weight  
logical :: dif_nodes                                   !different nodes per centre
integer, allocatable, dimension(:) :: nradc
integer, allocatable, dimension(:) :: nangc  !angular grid points (diff nodes per centre)
double precision :: a,b                                !limits of the integral
integer :: nrad               !number of Gauss-Legendre quadrature points (same nodes per centre)
double precision :: snglalpha                     !alpha parameter (same for all nodes)

!!!!Radial intracule at several distances!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical :: radial_plot                                !=.true if we want I vs R plot
character*40 :: r_plot_name                           !name of the radial plot
integer :: nblock                                     !number of blocks for the radius
integer, allocatable, dimension(:) :: n_an_per_part   !number of angluar points per block
double precision, allocatable, dimension(:,:) :: tart !space between radius
integer, allocatable, dimension(:) :: stp    !number of points for each block
integer :: gpt, nAng                         !angular grid points (same nodes per centre)
!vectorial plot with a cubefile!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical :: cubeintra
character*40 :: cubeintraname
double precision, dimension(3) :: center_i            !cube centered in (x,y,z)
double precision, dimension(3) :: step_i              !distance between points for each axis
integer, dimension(3) :: np_i                         !number of points in the cube for each axis

end module intrainfo

!--------------------------------------------------------------------------------------------------------

module quadratures
implicit none
!Angular+radial quadrature
double precision, allocatable, dimension(:,:) :: rg, rrrg !grid points for the integral (total, sym+weight)
integer ::  maxgrid, rrgrid !original and final number of (symmetry and weight reduced) grid points
double precision, allocatable, dimension(:) :: rweight !weight of rrrg grid points
!becke multicenter weight

!Angular quadrature (I(s) vs s)
integer :: rgrid   !number of sym reduced grid points
integer, allocatable, dimension(:) :: smn !number of points per quadrature center (Becke)
double precision, allocatable, dimension(:) :: radi !radial points
integer :: nradi !number of radial points
double precision, allocatable, dimension(:) :: w_ang  !angular weights
double precision, allocatable, dimension(:,:) :: rpg  !x,y,z total grid points
integer :: rpgrid !number of grid points for the plot (symmetry reduced)
!Gauss-hermite quadrature
double precision, allocatable, dimension(:) :: rh, w_r !nodes and weights for gauss hermite(coef.)
end module quadratures

!-----------------------------------------------------------------------------------------

module cubeinfo
   integer :: mo !molecular orbital to represent in the cubefile
   integer :: npr  !primitive to represent in the cubefile
   integer :: cao  !atomic orbital to represent in the cubefile
   double precision, dimension(3) :: center !cube centered in (x,y,z)
   double precision, dimension(3) :: step !distance between points for each axis
   integer, dimension(3) :: np !number of points in the cube for each axis
end module cubeinfo

module fractions
   double precision, parameter :: hr=1.d0*(3.d0**(-1.d0))
   double precision, parameter :: boh=5.d0*(3.d0**(-1.d0))
   double precision, parameter :: bs=1.d0*(6.d0**(-1.d0)) 
   double precision, parameter :: bih=2.d0*(3.d0**(-1.d0))
end module fractions        

module radis
   double precision, parameter, dimension(92) :: BL=(/0.327d0, 0.320d0, &
    + 1.219d0, 0.911d0, 0.793d0, 0.766d0, 0.699d0, 0.658d0,&
    + 0.900d0, 0.690d0, 1.545d0, 1.333d0, 1.199d0, 1.123d0, 1.110d0, 1.071d0, 1.039d0,&
    + 0.970d0, 1.978d0, 1.745d0, 1.337d0, 1.274d0, 1.236d0, 1.128d0, 1.180d0, 1.091d0,&
    + 1.089d0, 1.077d0, 1.146d0, 1.187d0, 1.199d0, 1.179d0, 1.209d0, 1.201d0, 1.201d0,&
    + 1.100d0, 2.217d0, 1.928d0, 1.482d0, 1.377d0, 1.353d0, 1.240d0, 1.287d0, 1.212d0,&
    + 1.229d0, 1.240d0, 1.362d0, 1.429d0, 1.385d0, 1.380d0, 1.421d0, 1.400d0, 1.397d0,&
    + 1.300d0, 2.442d0, 2.149d0, 1.653d0, 1.600d0, 1.600d0, 1.600d0, 1.600d0, 1.600d0,&
    + 1.600d0, 1.600d0, 1.500d0, 1.500d0, 1.500d0, 1.500d0, 1.500d0, 1.500d0, 1.500d0,&
    + 1.364d0, 1.346d0, 1.256d0, 1.258d0, 1.222d0, 1.227d0, 1.227d0, 1.273d0, 1.465d0,&
    + 1.531d0, 1.434d0, 1.496d0, 1.500d0, 1.500d0, 1.450d0, 1.500d0, 1.500d0, 1.500d0,&
    + 1.500d0, 1.500d0,1.500d0/)
end module radis

!----------------------------------------------------------------------------------------

module located
        contains      
                subroutine locate(iunit,string)
                integer :: ii        
                integer iunit
                character string*(*)
                character*80 linia

                rewind(iunit)

                ii=0
               do while(ii.eq.0)
                 read(iunit,"(a80)", end=10)linia
                 if(index(linia,string).ne.0) then
                 ii=1
                 return
                !   found=.true.
                 end if
               end do
                10 write(*,*) string, 'section not found'
                return
                end
end module located











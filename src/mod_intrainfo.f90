module intrainfo
!input variables for intracule calculation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
character*40 :: dm2name, outname
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
integer :: nrad               !number of Gauss-Legendre quadrature points (same nodes per centre)
double precision :: snglalpha                     !alpha parameter (same for all nodes)
!!!!Radial intracule at several distances!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical :: radial_plot                                !=.true if we want I vs R plot
character*40 :: r_plot_name                           !name of the radial plot
integer :: nblock                                     !number of blocks for the radius
integer, allocatable, dimension(:) :: n_an_per_part   !number of angluar points per block
double precision, allocatable, dimension(:,:) :: tart !space between radius
double precision, allocatable, dimension(:) :: stp    !step size for each block
integer :: gpt, nAng                         !angular grid points (same nodes per centre)
!vectorial plot with a cubefile!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical :: cubeintra
character*40 :: cubeintraname
double precision, dimension(3) :: center_i            !cube centered in (x,y,z)
double precision, dimension(3) :: step_i              !distance between points for each axis
integer, dimension(3) :: np_i  
logical :: vee_flag !number of points in the cube for each axis
!intracule at zero distance!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical :: intracule_at_zero
end module intrainfo


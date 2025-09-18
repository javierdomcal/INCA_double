module quadratures
implicit none
!Angular+radial quadrature
double precision, allocatable, dimension(:,:) :: rg, rrrg !grid points for the integral (total, sym+weight)
integer ::  maxgrid, rrgrid !original and final number of (symmetry and weight reduced) grid points
double precision, allocatable, dimension(:) :: rweight, rweight_vee !weight of rrrg grid points
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


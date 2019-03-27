module inputdata
real *8 allocatable::pbead(:,:,:),qbead(:,:,:,:),dvdq(:,:,:)!!!the momentum, position, and gradient of each bead of each atom
real *8 allocatable::posi_start(:,:),mass(:)  !!! the initial position of the molecule and the mass of each atom
real *8 allocatable::initq(:,:,:,:),initp(:,:,:,:) !!!initial position and momentum after thermostat
real *8 allocatable::umatrix(:,:),utrans(:,:)  !!! umatrix for transform to Nomal Mode Space
real *8 allocatable::centroidq(:,:),centroidp(:,:)  !!centroid of the ring polymer
!!! the centroid of molecule 
real *8 cm(3),pm(3)              

integer *4 sthermo,
integer *4 Natoms,ngas,nsurf


end module inputdata

module parameters
real *8 parameter:: C1=9.1093856e-31  !!! the mass of the electron
real *8 parameter:: C2=6.02214076e23  !!! NA 
real *8 parameter:: C3=0.03674931d0   !!! eV to hartree times factor
real *8 parameter:: C4=0.529177429d0  !!! bohr to angstrom times factor
real *8 parameter:: C5=1.38064852e-23 !!! kb boltzmann constant
real *8 parameter:: C6=4.359748e-18   !!! hartree to J times factor 
real *8 parameter:: C7=41.341374576d0 !!! fs to Au times factor
real *8 parameter:: C8=1822.888486d0  !!! amu to au (mass) 
real *8 pi
end module parameters

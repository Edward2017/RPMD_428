subroutine RPMD_VERLET()
use inputdata
use parameters
implicit none 
real *8 
integer::i,j,k

if(verlet_style.eq.2) then 
  call get_potential
  
  pbead(i,j,k)=pbead(i,j,k)-0.5d0*dt*dvdq(i,j,k)
    
  if(nbeads.eq.1) then
    do k=1,nbeads
      do j=1,natoms
        do i=1,3
          qbead(i,j,k)=qbead(i,j,k)+pbead(i,j,k)*dt/mass(j)
        enddo
      enddo
    enddo
  else
    call free_ring_polymer()
  endif
  
  call get_potential()
  
  pbead=pbead-0.5d0*dt*dvdq

endif

end subroutine RPMD_VERLET

subroutine free_ring_polymer()
use inputdata
use parameters
implicit none
real *8
integer::i,j,k




end subroutine free_ring_polymer



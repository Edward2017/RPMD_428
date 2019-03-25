subroutine thermostat(p,T,matom,info)
use inputdata
use parameters
implicit none
real *8 p(3,natoms,nbeads),T,matom(natoms)
interger :: i,j,k,nstep,info


do k=1,nbeads
	do i=1,ngas
		qbead(:,i,k)=posi_start(:,i)
	enddo
	do j=1+ngas,natoms
		qbead(:,j,K)=posi_start(:,j)
	enddo
enddo

Nsteps=Nequal+numtraj*step_sample
do nstep=1,Nsteps
	call RPMD_VERLET()
	call get_potential()
	call shift()
!!! for andersen thermostat 
	if(info.eq.0) then
		if(mod(nstep,nnequal).eq.0) then
			call sample_bead_momentum()
		endif
	endif
	
!!!
	if(nstep.ge.Nequal.and.mod(nstep-Nequal,step_sample).eq.0.and.Ntraj.lt.numtraj)then
		Ntraj=Ntraj+1
		initq(:,:,:,ntraj)=qbead(:,:,:)
		initp(:,:,:,ntraj)=pbead(:,:,:)
	endif
	
enddo


end subroutine thermostat

subroutine get_temperature(tnow)
use inputdata
use parameters
implicit none
real *8 tnow,etemp
integer ::i,j,k

Etarget=
etemp=0.0d0
do k=1,nbeads
	do j=1,ngas
		do i=1,3
			etemp=etemp+0.5d0*pbead(i,j,k)**2/mass(j)
		enddo
	enddo
enddo
tnow=etemp/Etarget*TA

end subroutine get_temperature

subroutine get_ZPE(ZPE)
use inputdata
use parameters
implicit none
real *8 
integer :: i,j,k

call get_kinetic_energy
call get_Ering_energy
call get_vtot_energy
ZPE=(Ek-Ering+vtot)/C3/dble(nbeads)
return
end subroutine get_ZPE


subroutine shift
use inputdata
use parameters
implicit none
real *8 center(3) 
integer::i,j,k

center=0.0d0
mm=0.0d0
do k=1,nbeads
	do j=1,natoms
		do i=1,3
			center(i)=center(i)+qbead(i,j,k)*mass(j)
		enddo
	enddo
enddo
do j=1,natoms
	mm=mm+mass(j)
enddo
if(center(3).le.5.0d0/C4)then
do j=1,natoms
	do k=1,nbeads
		qbead(3,j,k)=qbead(3,j,k)+10.0d0/C4
	enddo
enddo

end subroutine shift


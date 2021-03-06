subroutine thermostat(p,T,matom,info)
use inputdata
use parameters
implicit none
real *8 p(3,natoms,nbeads),T,matom(natoms),zpe,tem
interger :: i,j,k,nstep,info


do k=1,nbeads
	do i=1,ngas
		qbead(:,i,k)=posi_start(:,i)/C4
	enddo
	do j=1+ngas,natoms
		qbead(:,j,K)=posi_start(:,j)/C4
	enddo
enddo
open(100,file="samle.dat")
open(200,file="thermostat.dat")
Nsteps=Nequal+numtraj*step_sample
do nstep=1,Nsteps
	call RPMD_VERLET()
	call get_potential()
	call shift()
!!! for andersen thermostat 
	if(info.eq.0) then
		if(mod(nstep,nnequal).eq.0) then
			call sample_bead_momentum()
			call get_ZPE(zpe)
			call get_temperature(tem)
			write(333,'(I4,f18.10)')nstep,tem,zpe
		endif
	endif
	
!!!
	if(nstep.ge.Nequal.and.mod(nstep-Nequal,step_sample).eq.0.and.Ntraj.lt.numtraj)then
		Ntraj=Ntraj+1
		initq(:,:,:,ntraj)=qbead(:,:,:)
		initp(:,:,:,ntraj)=pbead(:,:,:)
		call get_ZPE(zpe)
		call get_temperature(tem)	
		write(100,'(I4,2f18.10)')Ntraj,tem,zpe
	endif
	
enddo
close(100)

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

subroutine autocorrelation(temper,N,coefficient)
implicit none
integer *4 N,i,j,k
real *8 temper(2*N),coefficient，avg,dsig

avg=0.0d0
do i=1,2*N
  avg=avg+temper(i)/(2*N)
enddo
dsig=0.0d0
do i=1,2*N
	dsig=dsig+(temper(i)-avg)**2
enddo
do i=1,N
	coefficient=coefficient+(temper(i)-avg)*(temper(i+N)-avg)/dsig
enddo
return
end subroutine autocorrelation

	
subroutine get_kinetic_energy
use inputdata,only:pbead,mass,ngas,nbeads,Ek
implicit none 
integer *4 i,j,k
Ek = 0.0d0
do i = 1, 3
  do j = 1, Ngas
    do k = 1, Nbeads
      Ek=Ek+0.5d0*pbead(i,j,k)*pbead(i,j,k)/mass(j)
    end do
  end do
end do	

end subroutine get_kinetic_energy


subroutine get_Ering_energy
use inputdata
use parameters
implicit none
real *8 dx,dy,dz,wn
integer *4 i,j,k
Ering = 0.0d0
dx=0.0d0
dy=0.0d0
dz=0.0d0
wn = Nbeads / beta
do j = 1, Ngas
  dx = q(1,j,1) - q(1,j,Nbeads)
  dy = q(2,j,1) - q(2,j,Nbeads)
  dz = q(3,j,1) - q(3,j,Nbeads)
  Ering=Ering+0.5d0*mass(j)*wn*wn*(dx*dx+dy*dy+dz*dz)

  do k = 2, Nbeads
    dx = q(1,j,k-1) - q(1,j,k)
    dy = q(2,j,k-1) - q(2,j,k)
    dz = q(3,j,k-1) - q(3,j,k)
    Ering=Ering+0.5d0*mass(j)*wn*wn*(dx*dx+dy*dy+dz*dz)
  end do
end do

end subroutine get_Ering_energy

subroutine get_total_Ep
use inputdata,only:Ep,nbeads,vtot
implicit none
integer::i,j,k
vtot=0.0d0
do i=1,nbeads
	vtot=vtot+Ep(i)
enddo
end subroutine get_total_Ep

subroutine read_initqp()
use inputdata,only:initq,inito,natoms,nbeads,numtraj
implicit none
integer :: i,j,k,n
open(300,file="initqp.dat")
read(300,*)
do n=1,numtraj
	do j=1,natoms
		do k=1,nbeads
			read(300,'(6f18.10)')initq(:,j,k,n),initp(:,j,k,n)
		enddo
	enddo
	read(300,*)
enddo
close(300)		
end subroutine read_initqp


subroutine init_random_seed()
integer :: i,n,clock
integer*4, allocatable:: seed(:)
call random_seed(size=n)
allocate (seed(n))

call system_clock(count=clock)

seed=clock + 37*(/(i-1,i=1,n)/)

call random_seed(put=seed)

eallocate(seed)
end subroutine init_random_seed





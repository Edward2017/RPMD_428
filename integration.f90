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
real *8 qtemp(3,natoms,nbeads),ptemp(3,natoms,nbeads),p1(3,natoms,nbeads),q1(3,natoms,nbeads)
real *8 poly(2,2,nbeads),wn,wk,wt,
integer::i,j,k,k2

qtemp=qbead
ptemp=pbead
p1=0.0d0
q1=0.0d0
!!!transform to Nomal mode space 
do i = 1, 3
  do j = 1,natoms
    do k=1,nbeads
      do k2=1,nbeads
        p1(i,j,k)=p1(i,j,k)+umatrix(k2,k)*utemp(i,j,k2)
        q1(i,j,k)=q1(i,j,k)+umatrix(k2,K)*qtemp(i,j,k2)
      enddo
    enddo
  end do
end do

!!! intergrate the position
do j=1,natoms
  poly(1,1,1)=1.0d0
  poly(1,2,1)=0.0d0
  poly(2,1,1)=dt/mass(j)
  poly(2,2,1)=1.0d0
  do k=2,nbeads
    beta_n=beta/dble(nbeads)
    wn=1.0/beta_n
    wk=2.0*wn*dsin(dble(k-1)*pi/dble(nbeads) )
    wt=wk*dt
    poly(1,1,k)=dcos(wt)
    poly(1,2,K)=-1.0d0*mass(j)*wk*dsin(wt)
    poly(2,1,K)=dsin(wt)/(mass(j)*wk)
    poly(2,2,k)=dcos(wt)
  enddo
  do k=1,nbeads
    do i=1,3
      ptemp(i,j,k)=poly(1,1,K)*p1(i,j,k)+poly(1,2,k)*q1(i,j,k)
      qtemp(i,j,K)=poly(2,1,K)*p1(i,j,k)+poly(2,2,k)*q1(i,j,k)
    enddo
  enddo
enddo
p1=0.0d0
q1=0.0d0

!!!transform to cartesian space      
do i = 1,3
  do j = 1,natoms
    do k = 1,nbeads
      do k2 = 1,nbeads
        ptemp(i,j,k)=ptemp(i,j,k)+Utrans(k2,k)*P1(i,j,k2)
        qtemp(i,j,k)=qtemp(i,j,k)+Utrans(k2,K)*q1(i,j,k2)
      enddo
    enddo
  end do       
end do

qbead=qtemp
pbead=ptemp

end subroutine free_ring_polymer

subroutine get_umatrix(n)
use inputdata,only:umatrix,utrans
use parameters,only:pi
implicit none 
real *8
integer:: i,j,k,n

ph=2.0*pi/dble(n)
factor1=dsqrt(1.0/dble(n))
if(.not.allocated(umatrix)) allocate(umatrix(n,n),utrans(n,n))
if(n.ne.1)then 
  do j=1,n
    umatrix(j,1)=factor1
  enddo
  do j=1,n
    umatrix(j,n/2+1)=factor1*(-1.0d0)**j
  enddo
  
  factor2=factor1*dsqrt(2.0d0)
  
  do k= 2, n/2
    do j= 1, n
      umatrix(j,k)=factor2*dcos(ph*dble(j*(k-1)))  
    enddo
  enddo
  
  do k= n/2+2, n
    do j= 1, n
      umatrix(j,k)=factor2*dsin(ph*dble(j*(k-1)))
    enddo
  enddo
else
  umatrix(1,1)=1.0d0
endif  

do i=1,n
  do j=1,n
    utrans(i,j)=umatrix(j,i)
  enddo
enddo

end subroutine get_umatrix






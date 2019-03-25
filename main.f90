program main
use inputdatas
use parameters
real *8 
integer *4 

!!! initiate the random_number_seed for random_number
call init_random_seed()

!!! read the parameters in the input_for_rpmd 
call read_inputdata()

!!! for some PES, init_pes is for the parameters of the PES to be read in
if(NINIT_PES.eq.1) call init_pes()

!!! Umatrix is needed for the translation to the normal mode space 
call get_umatrix(nbeads)

!!! define the file output the convergence of each trajectory
open(520,file="convergence.dat")

!!! initial the fundermetal parameters used in the program 
t=0.0d0
beta=C6/(TA*C5)
beta_n=beta/dble(nbeads)
NPATHS=0

!!! generate the initial position and momentum
call sample_bead_momentum()
call thermostat(pbead,TA,mass)


!!! steps means the number of the all trajectories
do nstep=1,steps
  qbead(:,:,:)=initq(:,:,:,nstep)
  pbead(:,:,:)=initp(:,:,:,nstep)
  
  call get_centroid()
  
  !!! shift molecule to the fixed position 
  call shift_molecule()
  
 !!! remove the translational energy of molecule
 if(rm_trans_energy.eq.1)  call remove_translational_energy()
  
  !!! translational energy is added on
  call add_translational_energy()

  !!! each trajectory begins here
  do tstep=1,tsteps
    call RPMD_VERLET()
    call get_centroid()
    call test
    if(NTEST.eq.2) goto 300
    Eall=0.0d0
    call get_kinetic_energy()
    call get_Ering_energy()
    call get_total_ep()
    Eall=(Ek+vtot+Ering)/C3/dble(nbeads)
    if(Eall.ge.HTMAX) HTMAX=Eall
    if(Eall.le.HTMIN) HTMIN=Eall
    if(swrite.eq.1) call gwrite()
    
   enddo
   
300  continue
write(520,'')"traj= ",nstep,"PATHS= ",NPATHS,"convergence= ",(HTMAX-HTMIN)*23.0605d0
close(200)
close(999)
enddo
close(520)

!!! for product analysis
   
end program main   
   

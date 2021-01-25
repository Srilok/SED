!-----------------------------------------------------------------------
! F. Leroy TU-Darmstadt - NOVEMBER 2009
! 
! The program MVACF computes the mass-weighted velocity autocorrelation
! funtion (MVACF).
! Inputs are velocities and masses.
! Outputs are the x,y and z projection of the MVACF and MVACF itself.
!-----------------------------------------------------------------------

program MVACF

  implicit none

  integer :: i,j
  integer :: natom,read_step,corr_length
  integer :: nframe_max,nframe_read
  integer :: pos_per_part
  integer :: corr_nbr
  integer :: shift_data
  character(len=*) :: filename

  real :: aux
  real :: time_interval,time_interval_w
  real, dimension(:,:), allocatable :: mvx,mvy,mvz
  real, dimension(:), allocatable :: acfx,acfy,acfz,acft
  real, dimension(:), allocatable :: mass,id,type

!-----------------------------------------------------------------------

!--- READ INPUT PARAMETERS FILE-----------------------------------------
  call getarg(1,filename)
  open(1, file = "input-mvacf")
  read(1,*) natom         
  read(1,*) nframe_max    
  read(1,*) time_interval 
  read(1,*) read_step     
  read(1,*) corr_length   
  read(1,*) shift_data    
  close(1)

  ! natom         is the number of atoms in the system
  ! nframe_max    is the number of frames in the trajectory
  ! time_interval is the time interval between two frames
  !               for example if the velocities were written every 1 fs
  !               time_interval is 0.001
  ! read_step     is the periodicity with which data are really used
  !               for example the trajectory contains 10,000 time steps
  !               if read_step is 2, time steps 2,4,6,... will be
  !               ignored and the number of frames really used will be
  !               10,000/2=5,000
  ! corr_length   is the correlation length in steps
  !               for example if the number of frames used is 8,500 and
  !               corr_length is 8,192 with time interval 0.001
  !               the correlation function will be computed every 0.001
  !               and will be 8.192 ps long
  ! shift_data    is the shifting register parameter. It corresponds to
  !               the number of steps between two time origins

!--- PRINT TO THE SCREEN FOR VERIFICATION ------------------------------
  print*, " "
  print*, "       --- YOUR PARAMETERS ---"
  print*, "Nbr. of atoms               ",natom
  print*, "Nbr. of frames              ",nframe_max
  print*, "Real time interval          ",time_interval
  print*, "Reading periodicity         ",read_step
  print*, "Correlation length          ",corr_length
  print*, "Shifting register parameter ",shift_data

  if (mod(nframe_max,read_step)/=0) then
     print*,"frame_max should be a multiple of read_step"
     stop
  end if
  nframe_read = nframe_max / read_step

!--- ALLOCATE MEMORY FOR VELOCITY ARRAYS -------------------------------
  !first dimension   number of atoms
  !second dimension  number of frames read
  allocate(x(natom,nframe_read))
  allocate(y(natom,nframe_read))
  allocate(z(natom,nframe_read))
  allocate(vx(natom,nframe_read))
  allocate(vy(natom,nframe_read))
  allocate(vz(natom,nframe_read))  

  allocate(mass(natom))
  allocate(id(natom))
  allocate(type(natom))

!--- READ VELOCITY FILE ------------------------------------------------

  open(2, file = filename)
  print*,">>Reading Trajectory<<"
  do i = 1,nframe_max
  read(2,*)    
  read(2,*)    
  read(2,*)    
  read(2,*)    
  read(2,*)    
  read(2,*)    
  read(2,*)    
  read(2,*)    
  read(2,*)    
     do j = 1, natom
	read(2,*) id(j),type(j),mass(j),x(j,i),y(j,i),z(j,i),vx(j,i),mvy(j,i),mvz(j,i)
     end do
  end do
  close(2)  

  print*,">>Trajectory read<<"

!--- WRITE RESULTS TO THE FILE vacf-tot.dat ----------------------------
  open(101, file = "vacf-tot-1.dat")
  do j = 1, corr_length
     write(101,'(1F15.5,4E14.6)') time_interval_w*(j-1), &
          acfx(j),acfy(j),acfz(j),acft(j)
  end do
  close(101) 
!-----------------------------------------------------------------------
contains
  
  subroutine ACF_CALC(xin,yin,zin,nptot, &
       corrl,corrn,shiftd,posppart, &
       xout,yout,zout,tout)
    
    !performs the average over time origins

    implicit none

    integer :: j,k,l,kk,ll,pp
    integer, intent(in) :: posppart,nptot
    integer, intent(in) :: corrn,corrl,shiftd
    real :: cx,cy,cz
    real :: cox,coy,coz
    real :: xin_aux,yin_aux,zin_aux
    real, intent(inout), dimension(nptot,posppart) :: xin,yin,zin
    real, intent(inout), dimension(corrl) :: xout,yout,zout,tout
    
    xout(:) = 0. ; yout(:) = 0. ; zout(:) = 0.
    do j = 1, corrn
       do k = 1, corrl
          cx = 0.
          cy = 0.
          cz = 0.
          do l = 1, nptot
             cox = xin(l,k) * xin(l,1)
             coy = yin(l,k) * yin(l,1)
             coz = zin(l,k) * zin(l,1)
             cx = cx + cox ; cy = cy + coy ; cz = cz + coz
          end do
          xout(k) = xout(k) + cx
          yout(k) = yout(k) + cy
          zout(k) = zout(k) + cz
       end do
       ! shifting register
       if (j < corrn) then
          do kk = shiftd + 1, posppart
             do ll = 1, nptot
                xin_aux = xin(ll,kk)
                xin(ll,kk-shift_data) = xin_aux
                yin_aux = yin(ll,kk)
                yin(ll,kk-shift_data) = yin_aux
                zin_aux = zin(ll,kk)
                zin(ll,kk-shift_data) = zin_aux
             end do
          end do
       end if
       if (mod(j,50)==0) then
          print*,"correlation ",j," / ",corrn," done"
       end if
    end do
    
    do k = 1, corrl
       xout(k) = xout(k) / corrn
       yout(k) = yout(k) / corrn
       zout(k) = zout(k) / corrn
       tout(k) = xout(k) + yout(k) + zout(k)
    end do
    
  end subroutine ACF_CALC
!-----------------------------------------------------------------------
end program MVACF
!.......................................................................
!=======================================================================


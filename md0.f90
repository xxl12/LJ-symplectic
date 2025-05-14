program main

  !::  MOLECULAR DYNAMICS SIMULATION IN A TWO DIMENSIONAL LENNARD-JONES SYSTEM  ::

  implicit none

  !:: size of the system
  integer, parameter :: nx=32, ny=32, na=2, nmax= nx*ny*na, nd=2

  !:: time period 
  integer, parameter :: nstep=100000

  !:: sample frequency and number of samples
  integer, parameter :: nsf= 20, nsamp= nstep/nsf

  !:: step size for the integrator
  real(8), parameter :: delt= 1D-2
  
  !:: particle position, velocity, force, potential, etc.
  real(8) :: x(nmax,nd), v(nmax,nd), f(nmax,nd), p(nmax), xs(nmax, nd)
  real(8) :: bsize(nd)

  !:: density 
  real(8), parameter :: rho= 0.93179839d0, vol= dble(nmax)/rho
  
  !:: Verlet list
  integer, parameter :: ncf=10 ! check the list every ncf steps 
  real(8), parameter :: rcut = 2.5d0, rlist=3.2d0
  real(8), parameter :: rcut2= rcut*rcut
  real(8), parameter :: rlst2= rlist*rlist 
  integer, parameter :: maxnb= 40
  integer :: point(nmax), list(nmax*maxnb)
  
  !:: adjust the potential energy
  real(8), parameter :: pcut= (rcut**(-12.0) - rcut**(-6.0))

  !:: Quantities being sampled
  real(8), dimension(nsamp) :: tmpr, pr, ekin, epot 
  
  integer :: i, n

  !:: initialize the system ... 
  call INIT_POS

  call INIT_VEL

  !:: time integration
  do n=1, nstep

     call INTEGRATE_RK
!call samp_prssr(1); print*, pr(1)
     call SAMP_DATA

  end do
     
  call OUTPUT
  
contains
  
  subroutine init_pos

    !:: initial the postion the the particle on a triangular lattice

    implicit none

    real(8) :: a0, dx, dy, dr(2), r2
    integer :: i, j, m, n

    a0= dsqrt( 2d0/dsqrt(3d0)/rho )

    dx= a0; dy= a0*dsqrt(3d0)

    bsize= (/ dble(nx)*dx, dble(ny)*dy /)

    x(1, :)= (/0d0, 0d0/) 
    x(2, :)= (/dx/2d0, dy/2d0/)

    n=0
    do i=1, nx
       do j=1, ny
          do m=1, 2 
             x(n+m, :)= x(m, :) + (/dble(i-1)*dx, dble(j-1)*dy/)
          end do
          n= n + 2
       end do
    end do

    call comp_force

    return

  end subroutine init_pos

  subroutine init_vel

    !:: initialize the velocity

    use ran_mod

    integer :: isd(4), i, k

    isd= (/3, 5, 7 ,9/)

    open(11, file='vel0.dat')
    do i=1, nmax
!       do k=1, nd
!          v(i, k)= dlarnd(3, isd)*.4d0
          read(11, *) v(i,1), v(i,2)
!       end do
    end do
    close(11)   

    !:: shift the velocity so that the mean is zero
    do k=1, nd
       v(:,k)=v(:,k) - sum(v(:,k))/dble(nmax)
    end do

    return

  end subroutine init_vel
  

  subroutine comp_force

    implicit none
    
    integer, save :: nf= -1
    logical :: update 

    real(8) :: dr(2), r2, r2i, r6i, ff(2), pp
    integer :: i, j, k, jbeg, jend, jnb, nlist

    nf= nf + 1
    if(mod(nf, ncf) .eq. 0) call check(update)

    !:: update the Verlet's list
    if (update) then

       call savex

       update= .false.

       !:: setting up the Verlet's list
       nlist=0
       do i = 1, nmax-1
          point(i)=nlist+1
          
          do j = i+1, nmax
             dr=x(i,:)-x(j,:)
             
             !:: periodic boundary condition
             r2= 0d0
             do k=1, nd
                dr(k)= dr(k) - dnint(dr(k)/bsize(k))*bsize(k)
                r2   = r2 + dr(k)*dr(k)
             end do

             if (r2 .lt. rlst2) then 
                !:: add the particle to the list
                nlist=nlist+1
                list(nlist)=j
             end if
          end do
       end do
       point(nmax)=nlist+1
    end if
    
    !:: calculate the force and energy
    f= 0d0; p= 0d0
    do i=1, nmax-1
       jbeg = point(i)
       jend = point(i+1)-1
       do jnb= jbeg, jend
          j  = list(jnb)
          dr=x(i,:)-x(j,:)

          !:: periodic boundary condition
          r2= 0d0
          do k=1, nd
             dr(k)= dr(k) - dnint(dr(k)/bsize(k))*bsize(k)
             r2   = r2 + dr(k)*dr(k)
          end do
          
          if(r2 .le. rcut2) then
             r2i= 1D0/r2
             r6i= r2i**3
             ff = r2i*r6i*(r6i-.5D0)*dr
             
             f(i,:)= f(i,:) + ff
             f(j,:)= f(j,:) - ff
             
             pp  = r6i*(r6i - 1d0) - pcut
             p(i)= p(i) + pp
             p(j)= p(j) + pp

          end if
       end do
    end do

    f= f*48d0
    p= p*4D0
    
    return

  end subroutine comp_force

  subroutine savex
    
    !:: save the array x to xs

    xs= x

  end subroutine savex

  subroutine check(update)

    !:: check the Verlet's list

    implicit none

    logical, intent(out) :: update
    
    real(8) :: dispmx, dr(2), r2
    integer :: i, k

    dispmx = 0D0
    do i = 1, nmax
       dr= x(i, :) - xs(i, :)
       
       r2= 0d0
       do k=1, nd
          dr(k)= dr(k) - dnint(dr(k)/bsize(k))*bsize(k)
          r2 = r2 + dr(k)*dr(k)
       end do
       
       dispmx = Dmax1(r2, dispmx)
    end do
    
    dispmx = 2D0 * Dsqrt (dispmx)
    update = (dispmx .gt. (rlist -rcut))
    
    return

  end subroutine check

  subroutine integrate_rk

    implicit none

    real(8), dimension(nmax,nd) :: xo, vo, v1, v2, v3, v4, f1, f2, f3, f4

    xo= x; vo= v

    call adj_pos
    call comp_force
    
    v1=v; f1=f
    
    x= xo + v1*delt/2d0
    v= vo + f1*delt/2d0

    call adj_pos
    call comp_force

    v2=v; f2=f
    
    x= xo + v2*delt/2d0
    v= vo + f2*delt/2d0

    call adj_pos
    call comp_force

    v3=v; f3=f
    
    x= xo + v3*delt
    v= vo + f3*delt

    
    call adj_pos
    call comp_force

    v4=v; f4=f

    x= xo + (v1 + 2d0*v2 + 2d0*v3 + v4)*delt/6d0
    v= vo + (f1 + 2d0*f2 + 2d0*f3 + f4)*delt/6d0

    return

  end subroutine integrate_rk

  subroutine integrate

    v= v + f*delt/2d0
    x= x + v*delt
    
    call adj_pos

    call comp_force
    
    v= v + f*delt/2d0

    return

  end subroutine integrate

  subroutine samp_data

    implicit none

    integer, save :: ntime = 0 
    integer :: nc

!    call samp_vacf
!    call samp_rdf

    ntime= ntime + 1
    if(mod(ntime, nsf) .ne. 0) then 
       return
    else
       nc= ntime/nsf
       call samp_engr (nc)
       call samp_tmpr (nc)
       call samp_prssr(nc)
    end if

    return
       
  end subroutine samp_data

  subroutine samp_engr(nc)

    !:: sample the energy

    implicit none

    integer, intent(in) :: nc

    !:: potential energy
    epot(nc)= sum(p)/dble(nmax*2)
    
    !:: kinetic energy
    ekin(nc)= sum(v*v)/dble(nmax*2)
    
    return

  end subroutine samp_engr

  subroutine samp_tmpr(nc)
    
    !:: sample the temperature 

    implicit none

    integer, intent(in) :: nc
    integer :: i, k

    tmpr(nc)= sum(v*v)/dble(nmax*nd)

    return

  end subroutine samp_tmpr

  subroutine samp_prssr(nc)

    !:: sample the velocity

    implicit none

    integer, intent(in) :: nc

    real(8) :: dr(2), r2, r2i, r6i, ff(2)
    integer :: i, j, k, jbeg, jend, jnb

    !:: calculate the force and energy
    pr(nc)=0d0
    do i=1, nmax-1
       jbeg = point(i)
       jend = point(i+1)-1
       do jnb= jbeg, jend

          j  = list(jnb)
          
          dr=x(i,:)-x(j,:)
          
          !:: periodic boundary condition
          r2= 0d0
          do k=1, nd
             dr(k)= dr(k) - dnint(dr(k)/bsize(k))*bsize(k)
             r2   = r2 + dr(k)*dr(k)
          end do
          
          if(r2 .le. rcut2) then
             r2i=1D0/r2
             r6i=r2i**3
             ff =r2i*r6i*(r6i-.5D0)*dr
          
             do k=1, nd
                pr(nc)= pr(nc) + ff(k)*dr(k)
             end do

          end if
       end do
    end do
    
    pr(nc)= pr(nc)*48d0
    do i=1, nmax
       do k=1, nd
          pr(nc) = pr(nc) + v(i,k)*v(i,k)
       end do
    end do
    
    pr(nc)= pr(nc)/(dble(nd)*vol)
    
    return

  end subroutine samp_prssr

  subroutine output

    implicit none

    integer :: i, n
    
    !:: save the position ...
    open(12, file= 'pos.dat')
    do i=1, nmax
       write(12, '(2F16.7)') x(i, :)
    end do
    close(12)   

    !:: save the velocity ...
    open(12, file= 'vel.dat')
    do i=1, nmax
       write(12, '(2F16.7)') v(i, :)
    end do
    close(12)
    
    !:: save the temperature
    open(13, file= 'tmpr.dat')
    do n=1, nsamp
       write(13, '(F16.7)') tmpr(n)
    end do
    close(13)

    !:: save the pressure
    open(14, file= 'pressure.dat')
    do n=1, nsamp
       write(14, '(F16.7)') pr(n)
    end do
    close(14)

    !:: save the pressure
    open(15, file= 'engr.dat')
    do n=1, nsamp
       write(15, '(2F16.7)') ekin(n), epot(n)
    end do
    close(15)

    return

  end subroutine output

  subroutine adj_pos

    implicit none

    integer :: i, k

    do i=1, nmax
       do k=1, nd
          if(x(i,k) .gt. bsize(k)) x(i,k)= x(i,k) - bsize(k)
          if(x(i,k) .lt. 0d0) x(i,k)= x(i,k) + bsize(k)
       end do
    end do

    return

  end subroutine adj_pos

  subroutine samp_vacf

    !:: sample the velocity auto-correlation ---


    !  0         it0       2*it0      3*it0      4*it0       5*it0     t0max*it0 
    !  |----------|----------|----------|----------|-----------|----------|
    !  


    implicit none

    !:: sample the velocity autocorrelation
    integer, parameter :: it0=200, nsamp=4, neq= 4000
    integer, parameter :: tmax = 2000
    integer, parameter :: t0max= tmax/it0+1
    real(8), save :: vacf(tmax), r2t(tmax), dtime, time_vr(tmax), vrdr(nd)
    real(8), save :: xx0(nmax, nd, t0max), vv0(nmax, nd, t0max)
    integer, save :: ntel, tt0, t1, idelt, time0(t0max), ntime(tmax)

    integer :: i, j, k

    integer, save :: nc= 0
    logical, save :: init_vacf=.true.

    nc= nc + 1
    if(mod(nc, nsamp).ne.0 .or. nc.lt.neq) return

    if(init_vacf) then
       dtime = delt * dble(nsamp)
       vacf = 0D0
       r2t = 0D0
       ntime = 0
       ntel =0
       tt0=0
       
       init_vacf=.false.
    end if

    ntel=ntel+1 

    !:: Define a new {t=0} and save the position and velocity 
    if ( mod(ntel-1, it0) .eq. 0) then
       tt0= tt0+1
       t1= mod(tt0-1, t0max) + 1

       time0(t1)=ntel
       do i=1, nmax
          xx0(i, :, t1) =x(i,:)
          vv0(i, :, t1) =v(i,:)
       end do

    end if

    !:: compute the correlation between the current time
    !:: and all the previous points where (t=0). 
    do j=1, min(tt0, t0max)    
       idelt= ntel- time0(j) +1
       if(idelt .le. tmax) then            
          ntime(idelt) = ntime(idelt) +1    
          do i=1, nmax
             vacf(idelt)= vacf(idelt) + sum(v(i,:) * vv0(i,:,j))

             do k=1, nd
                vrdr(k)= x(i,k)- xx0(i,k,j)
                vrdr(k)= vrdr(k)  - Dnint(vrdr(k)/bsize(k))*bsize(k)
                r2t (idelt)= r2t (idelt) + vrdr(k)**2
             end do

          end do
       end if
    end do
           
    if (nc .eq. nstep) then
       do i=1, tmax
          time_vr(i) = dtime* (i +.5D0)
          vacf(i) = vacf(i) /dble(nmax * ntime(i))
          r2t (i) = r2t (i) /dble(nmax * ntime(i))
       end do
       
       open (12, File='vacf0.dat')
       do j=1, tmax
          write(12, '(3E20.10)') time_vr(j), vacf(j), r2t(j)
       end do
       close(12)   

    end if
    
    return
    
  end subroutine samp_vacf





  
  subroutine samp_rdf

    !:: radial distribution function

    implicit none

    integer, parameter :: nbin= 100, nsamp=50, neq= 4000
    real(8), parameter :: rmax= 10d0, dg= rmax/dble(nbin)
    
    real(8), save :: g(nbin)=0, rg(nbin)
    integer, save :: ngt=0, nc=0

    real(8) :: dr(nd), r, r2, vb

    integer :: i, j, k, nh

    nc= nc + 1
    if(mod(nc,nsamp) .ne. 0 .or. nc.lt.neq) return

    ngt= ngt + 1
    do i=1, nmax-1
       do j=i+1, nmax
          dr= x(i,:)-x(j,:)
          
          r2= 0d0
          do k=1, nd
             dr(k)= dr(k)  - Dnint(dr(k)/bsize(k))*bsize(k)
             r2= r2 + dr(k)*dr(k)
          end do

          r = dsqrt(r2)
          nh= r/dg + 1

          if(nh.le.nbin) g(nh)= g(nh) + 2d0

       end do
    end do

    if(nc .eq. nstep) then
       do i=1, nbin
          rg(i) = dble(i-1)*dg

          vb = Dacos(-1d0)*dble(i**2 - (i-1)**2)*dg**2 !volume of the bin
          g(i)= g(i)/(nmax*rho*vb*ngt)
       end do

       open (12, File='rdf0.dat')
       do i=1, nbin
          write(12, '(2E20.10)') rg(i), g(i)
       end do
       close(12)   

    end if

  end subroutine samp_rdf

end program main

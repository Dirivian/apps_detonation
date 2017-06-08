! This is the non-localized source term for the model equation from Kasimov's 
! paper about chaos in shockwaves
c
c
c =========================================================
      subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
c =========================================================
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
c
c
c     # solve the diffusion equation q_t = q_{xx} using Crank-Nicolson
c     # The LAPACK tridiagonal solver dgtsv is used, which is in tridiag.f
c
      
c     # local storage:
      double precision :: beta
      double precision :: a
      real :: alpha_1
      double precision :: us
      double precision :: sho
      beta = 0.1d0
      pi = 4.d0*datan(1.0d0)
   
      alpha_1 = 4.5
      us =  q(1,mx+1)
      write(46,*) us
      sho=(1.d0/us)**(alpha_1)
      
      
      
      a= 1.d0/(4.d0*sqrt(4.d0*pi*beta)*(1.d0+erf(1.d0/(2*sqrt(beta)))))
      
      do i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         
         q(1,i) = q(1,i)+dt*a*exp(-(xcell+sho)**2/(4*beta))

	 enddo
c
      return
      end


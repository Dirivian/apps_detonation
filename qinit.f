c
c
c =========================================================
       subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
      double precision :: w1
      double precision :: w2
      double precision :: beta
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
c
c
      pi = 4.d0*datan(1.0d0)
      beta= 0.1d0
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         w1=(1+erf((1.0d0)/(2*sqrt(beta))))
         w2=(1+erf((xcell+1.0d0)/(2*sqrt(beta))))
         q(1,i) = 0.5d0*(1+sqrt(w2/w1))
  150    continue
c
      return
      end


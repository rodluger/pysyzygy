      subroutine transit(t,flux,bcirc,rhos,MpMs,esw,ecw,per,
     &           u1,u2,p,exptime,tN,nz,err,nt,ndata)
C  This routine computes the lightcurve for occultation
C  of a quadratically limb-darkened source without microlensing.
C  Please cite Mandel & Agol (2002) if you make use of this routine
C  in your research.  Please report errors or bugs to 
C  agol@tapir.caltech.edu
      implicit none
      integer i,nz,nt,ndata,j,k,istart,istop,s
      double precision exptime,bcirc,rhos,MpMs,esw,ecw,u1,u2,p,
     &       per,t(ndata),flux(ndata),tN(nt),favg,
     &       lambdad(nz),dt0,t0,window,
     &       etad(nz),lambdae(nz),lam,pi,x1,x2,x3,z,omega,kap0,
     &       kap1,q,Kk,Ek,Pk,n,ellec,ellk,rj,dt,G,DAYSEC,aRs,
     &       inc,becc,w,ecc,fi,tperi0,M,EA,f,rRs,fx,fxprime,deltax,
     &       kepler,keplerderiv
      integer err
      integer maxiter
      double precision tol
      integer NOERROR, KEPMAXITER
      parameter (NOERROR = 0, KEPMAXITER = 1)
      parameter (maxiter=20)
      parameter (tol=1.d-14)
      if(abs(p-0.5d0).lt.1.d-3) p=0.5d0
      err=NOERROR
C
C Input:
C
C t        array, length ndata of times (days) to compute the flux at
C flux     array, length ndata. Must be initialized with baseline flux!
C          Typically, this should be set to np.ones(ndata). The final 
C          flux will be stored in this param.
C bcirc    circular impact parameter
C rhos     stellar density (cgs)
C MpMs     Mp/Ms
C esw      e sin omega
C ecw      e cos omega
C per      orbital period (days)
C u1       linear limb-darkening coefficient (gamma_1 in paper)
C u2       quadratic limb-darkening coefficient (gamma_2 in paper)
C p        Rp/Rs
C exptime  The exposure time in days
C tN       An array of length nt with all the times of center-of-transit
C nz       # of points in exposure window over which to average the flux
C nt       The number of transits
C ndata    The number of data points
C
C
C Output:
C
C err      Error flags. 
C          (1) KEPMAXITER: Maximum number of iterations exceeded in
C              the Kepler problem
C
C Limb darkening has the form:
C  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*
C       (1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
C
C Now, compute pure occultation curve:
      omega=1.d0-u1/3.d0-u2/6.d0
      pi=acos(-1.d0)
C Calculate some variables
      G = 6.672d-8
      DAYSEC = 86400.d0
      aRs=((G*rhos*(1.d0+MpMs)*(per*DAYSEC)**2.d0)/
     &     (3.d0*pi))**(1.d0/3.d0)
      inc = acos(bcirc/aRs)
      w = atan2(esw,ecw)
      ecc = sqrt(esw**2.d0 + ecw**2.d0)
      fi = 3.d0*pi/2.d0 - w
      tperi0=per*sqrt(1.d0-ecc**2.d0)/(2.d0*pi)*(ecc*sin(fi)/
     &               (1.d0+ecc*cos(fi)) 
     & -2.d0/sqrt(1.d0-ecc**2.d0)*atan2(sqrt(1.d0-ecc**2)*
     &                           tan(fi/2.d0),1.d0+ecc))
C Equation (7) in Winn review
      becc = bcirc*(1.d0-ecc**2.d0)/(1.d0+ecc*sin(w-pi))
C Equation (14) in Winn review for the transit duration
      window=per/2.d0/pi*asin(((1.d0+p)**2-becc**2)**0.5d0/
     &       (sin(inc)*aRs))
C Correction for eccentricity
      window=window*sqrt(1.d0-ecc**2.d0)/(1.d0+ecc*sin(w-pi))
C Correction for blurring
      window=window+exptime
C Offsets
      dt0 = (exptime + exptime/nz)/(nz+1.d0)
      t0 = -0.5d0*exptime + exptime/(2.d0*nz)

C Now find the impact parameter at the right edge of the window
C corresponding to the first time point.
C If the planet is already occulting the star, we need
C to increase our window size. This is overkill for most
C systems, since the formula above is pretty accurate, but
C you never know, especially at high ecc.
      z = 0.d0      
      window = window/1.1d0
      do while (z.lt.(1.d0 + p))
        window = window*1.1d0
        M = 2.d0*pi/per*(t0 + dt0*(nz-1) - window - tperi0)    
        EA = M  
        do s=1,maxiter
          fx = kepler(EA,ecc,M)
          fxprime = keplerderiv(EA,ecc,M)
          if (abs(fx) < tol) then
            exit
          endif
          deltax = fx/fxprime
          EA = EA - deltax
        enddo
        if (s > maxiter) then
          if (abs(kepler(EA,ecc,M)) > tol) then
            err = KEPMAXITER
            return
          endif
        endif
        f = 2.d0*atan2((1.d0+ecc)**0.5d0*sin(EA/2.d0),
     &                 (1.d0-ecc)**0.5d0*cos(EA/2.d0))
        rRs = aRs*(1.d0-ecc**2.d0)/(1.d0+ecc*cos(f))
        z = rRs*sqrt(1.d0-(sin(w+f)*sin(inc))**2.d0)
      end do 
C And now the impact parameter at the left edge of the window
C corresponding to the last time point.
      z = 0.d0      
      window = window/1.1d0
      do while (z.lt.(1.d0 + p))
        window = window*1.1d0
        M = 2.d0*pi/per*(t0 + window - tperi0)    
        EA = M  
        do s=1,maxiter
          fx = kepler(EA,ecc,M)
          fxprime = keplerderiv(EA,ecc,M)
          if (abs(fx) < tol) then
            exit
          endif
          deltax = fx/fxprime
          EA = EA - deltax
        enddo
        if (s > maxiter) then
          if (abs(kepler(EA,ecc,M)) > tol) then
            err = KEPMAXITER
            return
          endif
        endif
        f = 2.d0*atan2((1.d0+ecc)**0.5d0*sin(EA/2.d0),
     &                 (1.d0-ecc)**0.5d0*cos(EA/2.d0))
        rRs = aRs*(1.d0-ecc**2.d0)/(1.d0+ecc*cos(f))
        z = rRs*sqrt(1.d0-(sin(w+f)*sin(inc))**2.d0)
      end do 
      
C Find istart and istop (indices corresponding to points of first & last
C contact).
      istop = 1
C Loop over all transits      
      do j=1,nt
      do istart=istop,ndata-1
        if(t(istart).ge.(tN(j)-window)) then
          exit
        endif
      enddo
      do istop=istart,ndata-1
        if(t(istop).gt.(tN(j)+window)) then
          exit
        endif
      enddo
C Loop over all indices within the transit
      do k=istart,istop
C It's really sad that this next line is necessary.
      if (istart.eq.istop) then
        exit
      endif      
      dt = t(k)-tN(j)
      favg = 0.d0
C Loop over each point within the exposure:
      do i=1,nz   
C Find impact parameter
        M = 2.d0*pi/per*(t0 + dt0*(i-1) + dt - tperi0) 
C Newton's method ---------------->     
        EA = M  
        do s=1,maxiter
          fx = kepler(EA,ecc,M)
          fxprime = keplerderiv(EA,ecc,M)
          if (abs(fx) < tol) then
            exit
          endif
          deltax = fx/fxprime
          EA = EA - deltax
        enddo
        if (s > maxiter) then
          if (abs(kepler(EA,ecc,M)) > tol) then
            err = KEPMAXITER
            return
          endif
        endif
C <---------------------------------   
        f = 2.d0*atan2((1.d0+ecc)**0.5d0*sin(EA/2.d0),
     &                 (1.d0-ecc)**0.5d0*cos(EA/2.d0))
        rRs = aRs*(1.d0-ecc**2.d0)/(1.d0+ecc*cos(f))
        z = rRs*sqrt(1.d0-(sin(w+f)*sin(inc))**2.d0)
C        
C
        x1=(p-z)**2
        x2=(p+z)**2
        x3=p**2-z**2
C the source is unocculted:
C Table 3, I.
        if(z.ge.1.d0+p) then
          lambdad(i)=0.d0
          etad(i)=0.d0
          lambdae(i)=0.d0
          goto 10
        endif
C the  source is completely occulted:
C Table 3, II.
        if(p.ge.1.d0.and.z.le.p-1.d0) then
          lambdad(i)=1.d0
          etad(i)=1.d0
          lambdae(i)=1.d0
          goto 10
        endif
C the source is partly occulted and the occulting object crosses limb:
C Equation (26):
        if(z.ge.abs(1.d0-p).and.z.le.1.d0+p) then
          kap1=acos(min((1.d0-p*p+z*z)/2.d0/z,1.d0))
          kap0=acos(min((p*p+z*z-1.d0)/2.d0/p/z,1.d0))
          lambdae(i)=p*p*kap0+kap1
          lambdae(i)=(lambdae(i)-0.5d0*sqrt(max(4.d0*z*z-
     &               (1.d0+z*z-p*p)**2,0.d0)))/pi
        endif
C the occulting object transits the source star (but doesn't
C completely cover it):
        if(z.le.1.d0-p) lambdae(i)=p*p
C the edge of the occulting star lies at the origin- special 
C expressions in this case:
        if(abs(z-p).lt.1.d-4*(z+p)) then
C Table 3, Case V.:
          if(z.ge.0.5d0) then
            lam=0.5d0*pi
            q=0.5d0/p
            Kk=ellk(q)
            Ek=ellec(q)
C Equation 34: lambda_3
            lambdad(i)=1.d0/3.d0+16.d0*p/9.d0/pi*(2.d0*p*p-1.d0)*Ek-
     &                 (32.d0*p**4-20.d0*p*p+3.d0)/9.d0/pi/p*Kk
C Equation 34: eta_1
            etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &              (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
            if(p.eq.0.5d0) then
C Case VIII: p=1/2, z=1/2
              lambdad(i)=1.d0/3.d0-4.d0/pi/9.d0
              etad(i)=3.d0/32.d0
            endif
            goto 10
          else
C Table 3, Case VI.:
            lam=0.5d0*pi
            q=2.d0*p
            Kk=ellk(q)
            Ek=ellec(q)
C Equation 34: lambda_4
            lambdad(i)=1.d0/3.d0+2.d0/9.d0/pi*(4.d0*(2.d0*p*p-1.d0)*Ek+
     &                 (1.d0-4.d0*p*p)*Kk)
C Equation 34: eta_2
            etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
            goto 10
          endif
        endif
C the occulting star partly occults the source and crosses the limb:
C Table 3, Case III:
        if((z.gt.0.5d0+abs(p-0.5d0).and.z.lt.1.d0+p).or.(p.gt.0.5d0.
     &      and.z.gt.abs(1.d0-p)*1.0001d0.and.z.lt.p)) then
          lam=0.5d0*pi
          q=sqrt((1.d0-(p-z)**2)/4.d0/z/p)
          Kk=ellk(q)
          Ek=ellec(q)
          n=1.d0/x1-1.d0
          Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
C Equation 34, lambda_1:
          lambdad(i)=1.d0/9.d0/pi/sqrt(p*z)*(((1.d0-x2)*(2.d0*x2+
     &        x1-3.d0)-3.d0*x3*(x2-2.d0))*Kk+4.d0*p*z*(z*z+
     &        7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
C Equation 34, eta_1:
          etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &          (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
          goto 10
        endif
C the occulting star transits the source:
C Table 3, Case IV.:
        if(p.le.1.d0.and.z.le.(1.d0-p)*1.0001d0) then
          lam=0.5d0*pi
          q=sqrt((x2-x1)/(1.d0-x1))
          Kk=ellk(q)
          Ek=ellec(q)
          n=x2/x1-1.d0
          Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
C Equation 34, lambda_2:
          lambdad(i)=2.d0/9.d0/pi/sqrt(1.d0-x1)*((1.d0-5.d0*z*z+p*p+
     &         x3*x3)*Kk+(1.d0-x1)*(z*z+7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
          if(abs(p+z-1.d0).le.1.d-4) then
            lambdad(i)=2/3.d0/pi*acos(1.d0-2.d0*p)-4.d0/9.d0/pi*
     &            sqrt(p*(1.d0-p))*(3.d0+2.d0*p-8.d0*p*p)
          endif
C Equation 34, eta_2:
          etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
        endif
 10     continue
C Now, using equation (33):
        favg=favg+1.d0-((1.d0-u1-2.d0*u2)*lambdae(i)+(u1+2.d0*u2)*
     &      lambdad(i)+u2*etad(i))/omega
      enddo
C Take the mean to get the flux integrated over the exposure window      
      flux(k)=flux(k)*favg/nz
      enddo
      enddo
      return
      end

      FUNCTION kepler(EA, ecc, M)
        REAL*8 kepler, EA, ecc, M
        kepler = EA - ecc*sin(EA) - M
        return
      END

      FUNCTION keplerderiv(EA, ecc, M)
        REAL*8 keplerderiv, EA, ecc, M
        keplerderiv = 1.d0 - ecc*cos(EA)
        return
      END

      FUNCTION rc(x,y)
      REAL*8 rc,x,y,ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2,
     *C3,C4
      PARAMETER (ERRTOL=.04d0,TINY=1.69d-38,SQRTNY=1.3d-19,BIG=3.d37,
     *TNBG=TINY*BIG,COMP1=2.236d0/SQRTNY,COMP2=TNBG*TNBG/25.d0,
     *THIRD=1.d0/3.d0,C1=.3d0,C2=1.d0/7.d0,C3=.375d0,C4=9.d0/22.d0)
      REAL*8 alamb,ave,s,w,xt,yt
      if(x.lt.0..or.y.eq.0..or.(x+abs(y)).lt.TINY.or.(x+
     *abs(y)).gt.BIG.or.(y.lt.-COMP1.and.x.gt.0..and.x.lt.COMP2)) 
     & return
      if(y.gt.0.d0)then
        xt=x
        yt=y
        w=1.
      else
        xt=x-y
        yt=-y
        w=sqrt(x)/sqrt(xt)
      endif
1     continue
        alamb=2.d0*sqrt(xt)*sqrt(yt)+yt
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        ave=THIRD*(xt+yt+yt)
        s=(yt-ave)/ave
      if(abs(s).gt.ERRTOL)goto 1
      rc=w*(1.d0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      FUNCTION rj(x,y,z,p)
      REAL*8 rj,p,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (ERRTOL=.05d0,TINY=2.5d-13,BIG=9.d11,C1=3.d0/14.d0,
     *C2=1.d0/3.d0,C3=3.d0/22.d0,C4=3.d0/26.d0,C5=.75d0*C3,
     *C6=1.5d0*C4,C7=.5d0*C2,C8=C3+C3)
CU    USES rc,rf
      REAL*8 a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,
     *ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt,rc,rf
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY.or.max(x,y,
     *z,abs(p)).gt.BIG) return
      sum=0.d0
      fac=1.d0
      if(p.gt.0.d0)then
        xt=x
        yt=y
        zt=z
        pt=p
      else
        xt=min(x,y,z)
        zt=max(x,y,z)
        yt=x+y+z-xt-zt
        a=1.d0/(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=rc(rho,tau)
      endif
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
        beta=pt*(pt+alamb)**2
        sum=sum+fac*rc(alpha,beta)
        fac=.25d0*fac
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        pt=.25d0*(pt+alamb)
        ave=.2d0*(xt+yt+zt+pt+pt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        delp=(ave-pt)/ave
      if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3.d0*ec
      ee=eb+2.d0*delp*(ea-ec)
      rj=3.d0*sum+fac*(1.d0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*
     *(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
      if (p.le.0.d0) rj=a*(b*rj+3.d0*(rcx-rf(xt,yt,zt)))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      function ellec(k)
      implicit none
      double precision k,m1,a1,a2,a3,a4,b1,b2,b3,b4,ee1,ee2,ellec
C Computes polynomial approximation for the complete elliptic
C integral of the second kind (Hasting's approximation):
      m1=1.d0-k*k
      a1=0.44325141463d0
      a2=0.06260601220d0
      a3=0.04757383546d0
      a4=0.01736506451d0
      b1=0.24998368310d0
      b2=0.09200180037d0
      b3=0.04069697526d0
      b4=0.00526449639d0
      ee1=1.d0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.d0/m1)
      ellec=ee1+ee2
      return
      end

      function ellk(k)
      implicit none
      double precision a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ellk,
     &       ek1,ek2,k,m1
C Computes polynomial approximation for the complete elliptic
C integral of the first kind (Hasting's approximation):
      m1=1.d0-k*k
      a0=1.38629436112d0
      a1=0.09666344259d0
      a2=0.03590092383d0
      a3=0.03742563713d0
      a4=0.01451196212d0
      b0=0.5d0
      b1=0.12498593597d0
      b2=0.06880248576d0
      b3=0.03328355346d0
      b4=0.00441787012d0
      ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1)
      ellk=ek1-ek2
      return
      end

      FUNCTION rf(x,y,z)
      REAL*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08d0,TINY=1.5d-38,BIG=3.d37,THIRD=1.d0/3.d0,
     *C1=1.d0/24.d0,C2=.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0.d0.or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)return
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0NL&WR2.

C TO BULID
C f2py -c transit.f -m transit
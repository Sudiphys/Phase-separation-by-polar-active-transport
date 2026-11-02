cc     Two-dimensional periodic simulation of revised S10a-S10c.

c     Fields: x=psi, cp=c_plus, cm=c_minus, plx=p_x, ply=p_y.
c     gfortran -O2 -ffixed-line-length-none active_phase_separation_final_commented.f -o simulation
c     ./simulation
c     ieee_is_finite detects NaN and infinite field values.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit real*8(a-h,o-z)
      character*8 names(901),name1(901),name2(900),name3(900)
      parameter(k=420,k2=k**2,khalf=k/2,n=1)
      
cc    k refers to the system size

      integer m(k),p(k),s1,is0
c     The last array index stores the current and next time levels.
      DOUBLE PRECISION x(k,k,2),cp(k,k,2),cm(k,k,2)
      DOUBLE PRECISION plx(k,k,2),ply(k,k,2)
      DOUBLE PRECISION pspx(k,k),pspy(k,k),psit(k,k),cdif(k,k,2)
      DOUBLE PRECISION conpx(k,k),conpy(k,k),conmx(k,k),conmy(k,k)
      DOUBLE PRECISION conpxp(k,k),conpyp(k,k),conmxm(k,k),conmym(k,k)
      DOUBLE PRECISION bpsx(k,k),bpsy(k,k),cpt(k,k),cmt(k,k)
      DOUBLE PRECISION amucp(k,k),amucm(k,k),ajcpx(k,k),ajcmx(k,k)
      DOUBLE PRECISION ajcpy(k,k),ajcmy(k,k)
      DOUBLE PRECISION alapsi(k,k),divsp2(k,k),agrp2(k,k),divlps(k,k)
      DOUBLE PRECISION dpcdf(k,k),dpcdfpsi(k,k),apOP(k,k)
      DOUBLE PRECISION aggp2(k,k),aggpx2(k,k),aggpy2(k,k)
      DOUBLE PRECISION divlps1(k,k),divlps2(k,k),dpcdfpsi2(k,k)
      DOUBLE PRECISION divpsx(k,k),divpsy(k,k),dpcdfpsi1(k,k) 
      DOUBLE PRECISION aggpx2s(k,k),aggpy2s(k,k)
      DOUBLE PRECISION agdpx(k,k),agdpy(k,k) 
      DOUBLE PRECISION agdcpxx(k,k),agdcpyy(k,k)
      DOUBLE PRECISION agdcmxx(k,k),agdcmyy(k,k)
      DOUBLE PRECISION da1(k,k),da2(k,k),da3(k,k),da4(k,k),da5(k,k)
      DOUBLE PRECISION da6(k,k),da7(k,k),phi(k,k,2)
      DOUBLE PRECISION cd(k,k),ce(k,k),cf(k,k),cg(k,k),ch(k,k)

      
      REAL(KIND=8) :: Dr, lp, v0, alpha, Le, lsq, d, bet
      REAL(KIND=8) :: anu, Dpsi, Dc, azeta, arho, akap
      REAL(KIND=8) :: Dmt, Dm, anup, azetap, vp, ara

      DOUBLE PRECISION vmfield(k,k),motors
c     First snapshot at step 10000; the interval then doubles.
      integer next_save
      next_save = 10000

       iseed=-137571
       iseed1=-137573
c       v0=-delv0
c       v1=0.03
c       nv0=1

c     Read four lists of 900 filenames, each at most 8 characters.
c     Only names, read from fnamerhol, supplies snapshot filenames.
              open(unit=1,file='fnamerhol')
              open(unit=2,file='fnamerholp')
              open(unit=3,file='fnamerholp5')
              open(unit=4,file='fnamerholp8')

ccc           to save all data in a file
c     data.dat is opened but not written in this version.
c     data_comp_polarity.dat records the derived parameters.
              open(unit=11,file='data.dat')
              open(unit=12,file='data_comp_polarity.dat')

              read(1,*)(names(i),i=1,900)
              read(2,*)(name1(i),i=1,900)
              read(3,*)(name2(i),i=1,900)
              read(4,*)(name3(i),i=1,900)

       pi=2.0d0*asin(1.0d0)
c     Dimensionless timestep and grid spacing.
c     dx2=2*dx for centered gradients; dxx=dx**2 for Laplacians.
       dt=0.000075d0
       dx=0.015d0
       dx2=2*dx
       dxx=dx**2

ccc   average values of different species

c     Mean initial densities and polarization components.
      avgpsi=0.2
      avgcp=0.15
      avgcm=0.15
      avgp=0.0      



       ntime1=0
       ntime=200000000
cccc    defining neighbor latticle points, important for finte
ccc     difference method
c     Periodic neighbors: p is the next index, m the previous index.
          do i=1,k
          p(i)=i+1
          if(i.eq.k)p(i)=1
          m(i)=i-1
          if(i.eq.1)m(i)=k
          enddo
cccc    starting defining initial conditions
        do j=1,k
        do i=1,k

c     Initial fields have uniform random perturbations.
c     Half-widths: 0.025 for psi and p; 0.005 for each motor.
c     den and theta draws are retained exactly as in the source.
        den=(ran2(iseed))
c       psi       
        x(i,j,1)=avgpsi+0.05*(0.5-ran2(iseed))

c       one type of motor       
        cp(i,j,1)=avgcp+0.01*(0.5-ran2(iseed1))
c 0.01*sin(2*pi*50*i/k)*sin(2*pi*50*j/k)

c       other type of motor        
        cm(i,j,1)=avgcm+0.01*(0.5-ran2(iseed))
c cos(2*pi*50*i/k)*cos(2*pi*50*j/k)

c       difference between motors        
        cdif(i,j,1)=cp(i,j,1)-cm(i,j,1)
c       random angle        
        theta=(ran2(iseed)-0.5)*2*pi
c       polarization components        
        plx(i,j,1)=avgp+0.05*(0.5-ran2(iseed))
c cos(theta)
        ply(i,j,1)=avgp+0.05*(0.5-ran2(iseed))
c sin(theta)
        enddo
        enddo


        nrho=100
        npic=100
        ndef=100
        ndef1=100
        npro=100

ccc    parameters
c     Physical parameters and coefficients used for S10.
c     lp=Le=L; lsq=2*L; alpha and bet are alpha and beta.
        Dr    = 0.1D0
        lp    = 5.0D0
        v0    = 0.25D0
        alpha = 1000.0D0
        Le    = 5.0D0
        lsq   = 2.0D0*Le
        d     = 0.02D0
        bet   = 1.0D0
c     anu is the effective dimensional coupling nu.
        anu   = 0.007D0 
        Dpsi  = 0.1D0
        Dc    = 1.0D0
        azeta = 0.0000063D0
        arho  = 1.0D0/(pi*0.02D0*0.02D0)
c       akap = kappa_tilde / psi_0.
        akap   = 0.00005D0

        
c     Dmt=D_mt_tilde and Dm=D_tilde.
        Dmt=(Dpsi)/(Dr*lp*lp)
        Dm=(Dc)/(Dr*lp*lp)
c     anup = nu*rho_0/L = nu_tilde/v0_tilde.
c     With vp=0.5, anup=1.114085 and nu_tilde=0.557042.
c     All three S10 coupling terms use this same anup.
        anup=anu*arho/(lp)
c     azetap = 2*zeta_tilde/(c_plus_0+c_minus_0).
        azetap=azeta*arho/(lp*lp)
        azetapp=azetap*2*avgcp/(avgcp+avgcm)
        azetapm=azetap*2*avgcm/(avgcp+avgcm)
c     vp=v0_tilde; ara=1/R in the active motor transport terms.
        vp=v0/(Dr*lp)
        ara=Dc*DLOG(Le/d)/(pi*d*v0)

        
c        write(*,*)vp,Dmt,Dm,anup,azetap,akap

c     All right-hand sides use iold; Euler updates go into inew.
        iold=1
        inew=2
        itp=0
        do itime=1,ntime

        if(itime.eq.1) then
        write(12,*)vp,Dmt,Dm,anup,azetap,akap,ara
        endif

c     Validate physical denominators; do not clip conserved fields.
        do j=1,k
        do i=1,k
c     atot=psi+c_plus+c_minus; 1-atot is the free fraction.
        atot=x(i,j,iold)+cp(i,j,iold)+cm(i,j,iold)
        if (.not.ieee_is_finite(atot)
     &      .or..not.ieee_is_finite(plx(i,j,iold))
     &      .or..not.ieee_is_finite(ply(i,j,iold))) then
        write(*,*) 'Non-finite fields at step/cell:',itime,i,j
        stop 1
        endif
        if (x(i,j,iold).le.0.0d0.or.atot.ge.1.0d0
     &      .or.cp(i,j,iold).lt.0.0d0
     &      .or.cm(i,j,iold).lt.0.0d0) then
        write(*,*) 'Invalid fields at step/cell:',itime,i,j
        write(*,*) x(i,j,iold),cp(i,j,iold),cm(i,j,iold)
        stop 1
        endif
c     vmfield = vp*(cp-cm)/(cp+cm) = -v_mt_tilde.
c     Use zero when the total motor density is zero.
        motors=cp(i,j,iold)+cm(i,j,iold)
        vmfield(i,j)=0.0d0
        if (motors.gt.0.0d0) then
        vmfield(i,j)=vp*(cp(i,j,iold)-cm(i,j,iold))/motors
        endif
        enddo
        enddo

c     Pass 1: construct local scalars and transport vectors.
c     Complete this pass over the entire grid before pass 2.
c     Full-grid intermediate pass 1
        do j=1,k
        do i=1,k
        pxx=(plx(p(i),j,iold)-plx(m(i),j,iold))/dx2
        pyy=(ply(i,p(j),iold)-ply(i,m(j),iold))/dx2
c     pspx is the scalar psi**2*div(p), for polarization elasticity.
        pspx(i,j)=(pxx+pyy)*x(i,j,iold)*x(i,j,iold)
        polOP=plx(i,j,iold)*plx(i,j,iold)+ply(i,j,iold)*ply(i,j,iold)
        vm=vmfield(i,j)
        cdiffx=(cdif(p(i),j,iold)-cdif(m(i),j,iold))/dx2
        cdiffy=(cdif(i,p(j),iold)-cdif(i,m(j),iold))/dx2
c     S10b active flux: vm*psi*p; its divergence is added.
         conpx(i,j)=vm*x(i,j,iold)*plx(i,j,iold)
         conpy(i,j)=vm*x(i,j,iold)*ply(i,j,iold)
c     A=|p|**2*(1+alpha*|p|**2/2).
c     B=p.grad(cp-cm); stored in dpcdf.
        agrp2(i,j)=polOP*(1.0+(alpha/2)*polOP)
        dpcdf(i,j)=cdiffx*plx(i,j,iold)+cdiffy*ply(i,j,iold)
         apOP(i,j)=sqrt(plx(i,j,iold)**2+ply(i,j,iold)**2)    
        psixx=(x(p(i),j,iold)-2*x(i,j,iold)+x(m(i),j,iold))/dxx
        psiyy=(x(i,p(j),iold)-2*x(i,j,iold)+x(i,m(j),iold))/dxx
        alapsi(i,j)=(psixx+psiyy)
        agdpsx=(x(p(i),j,iold)-x(m(i),j,iold))/dx2
        agdpsy=(x(i,p(j),iold)-x(i,m(j),iold))/dx2
        atot=x(i,j,iold)+cp(i,j,iold)+cm(i,j,iold)
c     S10b crowding vector: (2*L/(pi*d))*psi/(1-atot)*grad(psi).
        agdpx(i,j)=(lsq/(pi*d))*(x(i,j,iold)/(1.0-atot))*(agdpsx)
        agdpy(i,j)=(lsq/(pi*d))*(x(i,j,iold)/(1.0-atot))*(agdpsy)
c     H=psi*(div(p))**2; bpsx,bpsy store psi*p.
        aggp2(i,j)=(pxx+pyy)*(pxx+pyy)*x(i,j,iold)
              bpsx(i,j)=plx(i,j,iold)*x(i,j,iold)
              bpsy(i,j)=ply(i,j,iold)*x(i,j,iold)
c     S10a active vectors: ara*(-vp+vm)*psi*cp*p for plus,
c     and ara*(vp+vm)*psi*cm*p for minus.
c     These vectors enter the density equations with +divergence.
         conpxp(i,j)=ara*(-vp+vm)*x(i,j,iold)*cp(i,j,iold)*plx(i,j,iold)
         conpyp(i,j)=ara*(-vp+vm)*x(i,j,iold)*cp(i,j,iold)*ply(i,j,iold)
         conmxm(i,j)=ara*(vp+vm)*x(i,j,iold)*cm(i,j,iold)*plx(i,j,iold)
         conmym(i,j)=ara*(vp+vm)*x(i,j,iold)*cm(i,j,iold)*ply(i,j,iold)
        agdcpx=(cp(p(i),j,iold)-cp(m(i),j,iold))/dx2
        agdcpy=(cp(i,p(j),iold)-cp(i,m(j),iold))/dx2
c     S10a plus-motor crowding vector: cp/(1-atot)*grad(cp).
        agdcpxx(i,j)=(cp(i,j,iold)/(1.0-atot))*(agdcpx)
        agdcpyy(i,j)=(cp(i,j,iold)/(1.0-atot))*(agdcpy)
         agdcmx=(cm(p(i),j,iold)-cm(m(i),j,iold))/dx2
         agdcmy=(cm(i,p(j),iold)-cm(i,m(j),iold))/dx2
c     S10a minus-motor crowding vector: cm/(1-atot)*grad(cm).
        agdcmxx(i,j)=(cm(i,j,iold)/(1.0-atot))*(agdcmx)
        agdcmyy(i,j)=(cm(i,j,iold)/(1.0-atot))*(agdcmy)
        enddo
        enddo

c     Pass 2: form psi*grad(A), psi*grad(B), and grad(H).
c     Also construct the motor gradient-coupling potentials.
c     Full-grid intermediate pass 2
        do j=1,k
        do i=1,k
        agrpx2=(agrp2(p(i),j)-agrp2(m(i),j))/dx2
        agrpy2=(agrp2(i,p(j))-agrp2(i,m(j)))/dx2
        divpsx(i,j)=x(i,j,iold)*(agrpx2)
        divpsy(i,j)=x(i,j,iold)*(agrpy2)
        dpcdfx=(dpcdf(p(i),j)-dpcdf(m(i),j))/dx2
        dpcdfy=(dpcdf(i,p(j))-dpcdf(i,m(j)))/dx2
        dpcdfpsi1(i,j)=(dpcdfx)*x(i,j,iold)
        dpcdfpsi2(i,j)=(dpcdfy)*x(i,j,iold)
        aggpx2(i,j)=(aggp2(p(i),j)-aggp2(m(i),j))/dx2
        aggpy2(i,j)=(aggp2(i,p(j))-aggp2(i,m(j)))/dx2
         acpxx=(cp(p(i),j,iold)-2*cp(i,j,iold)+cp(m(i),j,iold))/dxx
         acpyy=(cp(i,p(j),iold)-2*cp(i,j,iold)+cp(i,m(j),iold))/dxx
         acmxx=(cm(p(i),j,iold)-2*cm(i,j,iold)+cm(m(i),j,iold))/dxx
         acmyy=(cm(i,p(j),iold)-2*cm(i,j,iold)+cm(i,m(j),iold))/dxx
            bpsxx=(bpsx(p(i),j)-bpsx(m(i),j))/dx2
            bpsyy=(bpsy(i,p(j))-bpsy(i,m(j)))/dx2
c     mu_plus/minus = -azetap*Laplacian(c_plus/minus)
c                      +/- anup*pi*d/(2*L)*div(psi*p).
         amucp(i,j)=-azetap*(acpxx+acpyy)+anup*(pi*d/lsq)*(bpsxx+bpsyy)
         amucm(i,j)=-azetap*(acmxx+acmyy)-anup*(pi*d/lsq)*(bpsxx+bpsyy)
        enddo
        enddo

c     Pass 3: form psi*grad(H) and c_plus/minus*grad(mu).
c     Complete all vectors before taking their divergences.
c     Full-grid intermediate pass 3
        do j=1,k
        do i=1,k
         aggpx2s(i,j)=x(i,j,iold)*aggpx2(i,j)
         aggpy2s(i,j)=x(i,j,iold)*aggpy2(i,j)
          gamucpx=(amucp(p(i),j)-amucp(m(i),j))/dx2
          gamucpy=(amucp(i,p(j))-amucp(i,m(j)))/dx2
          ajcpx(i,j)=cp(i,j,iold)*(gamucpx)
          ajcpy(i,j)=cp(i,j,iold)*(gamucpy)
          gamucmx=(amucm(p(i),j)-amucm(m(i),j))/dx2
          gamucmy=(amucm(i,p(j))-amucm(i,m(j)))/dx2
          ajcmx(i,j)=cm(i,j,iold)*(gamucmx)
          ajcmy(i,j)=cm(i,j,iold)*(gamucmy)
        enddo
        enddo

c     Assemble S10 right-hand sides from the current time level.
c     All intermediate arrays now belong to this timestep.
        do j=1,k
        do i=1,k
        pxx=(plx(p(i),j,iold)-plx(m(i),j,iold))/dx2
        pxxy=(plx(i,p(j),iold)-plx(i,m(j),iold))/dx2
        pyy=(ply(i,p(j),iold)-ply(i,m(j),iold))/dx2
        pyyx=(ply(p(i),j,iold)-ply(m(i),j,iold))/dx2
        pxpx=(pspx(p(i),j)-pspx(m(i),j))/dx2
        pypy=(pspx(i,p(j))-pspx(i,m(j)))/dx2
        polOP=plx(i,j,iold)*plx(i,j,iold)+ply(i,j,iold)*ply(i,j,iold)
        vm=vmfield(i,j)
        cdiffx=(cdif(p(i),j,iold)-cdif(m(i),j,iold))/dx2
        cdiffy=(cdif(i,p(j),iold)-cdif(i,m(j),iold))/dx2
c     S10c x component: local relaxation, motor coupling,
c     elasticity, and advection. Divide their sum by psi.
c     The final coupling is (anup/bet)*grad(cp-cm).
c     The final elastic term is akap/psi*grad(psi**2*div(p)).
        tpbetax=-plx(i,j,iold)*(1.0+alpha*polOP)*x(i,j,iold)
        tnux=anup*cdiffx*x(i,j,iold)
        tpkapx=akap*(pxpx)
c     advpx is negative vm*(p.grad)p_x; subtracting it gives
c     +vm*(p.grad)p_x = -v_mt_tilde*(p.grad)p_x.
        advpx=-vm*plx(i,j,iold)*pxx-vm*ply(i,j,iold)*pxxy
        addpx=tpbetax+(tnux/bet)+tpkapx-advpx*x(i,j,iold)
        addpx=addpx/x(i,j,iold)
c     S10c y component, with the same coefficients and signs.
        tpbetay=-ply(i,j,iold)*(1.0+alpha*polOP)*x(i,j,iold)
        tnuy=anup*cdiffy*x(i,j,iold)
        tpkapy=akap*(pypy)
        advpy=-vm*ply(i,j,iold)*pyy-vm*plx(i,j,iold)*pyyx
        addpy=tpbetay+(tnuy/bet)+tpkapy-advpy*x(i,j,iold)
        addpy=addpy/x(i,j,iold)
         conpxx=(conpx(p(i),j)-conpx(m(i),j))/dx2
         conpyy=(conpy(i,p(j))-conpy(i,m(j)))/dx2
c     S10b active contribution: div(vm*psi*p).
        addcompsi=(conpxx+conpyy)
        agrpx2=(agrp2(p(i),j)-agrp2(m(i),j))/dx2
        agrpy2=(agrp2(i,p(j))-agrp2(i,m(j)))/dx2
        divspx2=(divpsx(p(i),j)-divpsx(m(i),j))/dx2
        divspy2=(divpsy(i,p(j))-divpsy(i,m(j)))/dx2
        dpcdfx=(dpcdf(p(i),j)-dpcdf(m(i),j))/dx2
        dpcdfy=(dpcdf(i,p(j))-dpcdf(i,m(j)))/dx2
        dinupsx=(dpcdfpsi1(p(i),j)-dpcdfpsi1(m(i),j))/dx2
        dinupsy=(dpcdfpsi2(i,p(j))-dpcdfpsi2(i,m(j)))/dx2
        psixx=(x(p(i),j,iold)-2*x(i,j,iold)+x(m(i),j,iold))/dxx
        psiyy=(x(i,p(j),iold)-2*x(i,j,iold)+x(i,m(j),iold))/dxx
        agdpsx=(x(p(i),j,iold)-x(m(i),j,iold))/dx2
        agdpsy=(x(i,p(j),iold)-x(i,m(j),iold))/dx2
        atot=x(i,j,iold)+cp(i,j,iold)+cm(i,j,iold)
        agdpsxx=(agdpx(p(i),j)-agdpx(m(i),j))/dx2
        agdpsyy=(agdpy(i,p(j))-agdpy(i,m(j)))/dx2
        daggpx2=(aggpx2s(p(i),j)-aggpx2s(m(i),j))/dx2
        daggpy2=(aggpy2s(i,p(j))-aggpy2s(i,m(j)))/dx2
c     S10b diffusion plus crowding.
c     tbetpsi=(bet/2)*div(psi*grad(A)).
c     tnupsict=-anup*div(psi*grad(B)).
c     tkapsi=bet*akap*div(psi*grad(H)).
c     Multiply their sum by Dmt, then add active transport.
        tlapsi=(psixx+psiyy+(agdpsxx+agdpsyy))
        tbetpsi=(bet/2)*(divspx2+divspy2)
        tnupsict=-anup*(dinupsx+dinupsy)
        tkapsi=bet*akap*(daggpx2+daggpy2)
        addpsitotal=(tlapsi+tbetpsi+tnupsict+tkapsi)*Dmt
         acpxx=(cp(p(i),j,iold)-2*cp(i,j,iold)+cp(m(i),j,iold))/dxx
         acpyy=(cp(i,p(j),iold)-2*cp(i,j,iold)+cp(i,m(j),iold))/dxx
         acmxx=(cm(p(i),j,iold)-2*cm(i,j,iold)+cm(m(i),j,iold))/dxx
         acmyy=(cm(i,p(j),iold)-2*cm(i,j,iold)+cm(i,m(j),iold))/dxx
            bpsxx=(bpsx(p(i),j)-bpsx(m(i),j))/dx2
            bpsyy=(bpsy(i,p(j))-bpsy(i,m(j)))/dx2
          gamucpx=(amucp(p(i),j)-amucp(m(i),j))/dx2
          gamucpy=(amucp(i,p(j))-amucp(i,m(j)))/dx2
          gamucmx=(amucm(p(i),j)-amucm(m(i),j))/dx2
          gamucmy=(amucm(i,p(j))-amucm(i,m(j)))/dx2
          divcpx=(ajcpx(p(i),j)-ajcpx(m(i),j))/dx2
          divcpy=(ajcpy(i,p(j))-ajcpy(i,m(j)))/dx2 
          divcmx=(ajcmx(p(i),j)-ajcmx(m(i),j))/dx2
          divcmy=(ajcmy(i,p(j))-ajcmy(i,m(j)))/dx2
         conpxxp=(conpxp(p(i),j)-conpxp(m(i),j))/dx2
         conpyyp=(conpyp(i,p(j))-conpyp(i,m(j)))/dx2
         conmxxm=(conmxm(p(i),j)-conmxm(m(i),j))/dx2
         conmyym=(conmym(i,p(j))-conmym(i,m(j)))/dx2
        agdcpx=(cp(p(i),j,iold)-cp(m(i),j,iold))/dx2
        agdcpy=(cp(i,p(j),iold)-cp(i,m(j),iold))/dx2
        agdcpsxxx=(agdcpxx(p(i),j)-agdcpxx(m(i),j))/dx2
        agdcpsyyy=(agdcpyy(i,p(j))-agdcpyy(i,m(j)))/dx2
         agdcmx=(cm(p(i),j,iold)-cm(m(i),j,iold))/dx2
         agdcmy=(cm(i,p(j),iold)-cm(i,m(j),iold))/dx2
        agdcmsxxx=(agdcmxx(p(i),j)-agdcmxx(m(i),j))/dx2
        agdcmsyyy=(agdcmyy(i,p(j))-agdcmyy(i,m(j)))/dx2
c     S10a plus RHS: active divergence plus Dm times diffusion,
c     gradient-coupling divergence, and crowding divergence.
       addtotcp=conpxxp+conpyyp+Dm*(acpxx+acpyy+divcpx+divcpy)
       addtotcp=addtotcp+Dm*(agdcpsxxx+agdcpsyyy)  
c     S10a minus RHS, with its own active and coupling signs.
       addtotcm=conmxxm+conmyym+Dm*(acmxx+acmyy+divcmx+divcmy)
       addtotcm=addtotcm+Dm*(agdcmsxxx+agdcmsyyy)
c     Forward Euler: new field = old field + dt*right-hand side.
c     Store the updated motor difference and output fields.
         x(i,j,inew)=x(i,j,iold)+(addpsitotal+addcompsi)*dt
        cp(i,j,inew)=cp(i,j,iold)+addtotcp*dt
        cm(i,j,inew)=cm(i,j,iold)+addtotcm*dt
        plx(i,j,inew)=plx(i,j,iold)+(addpx)*dt
        ply(i,j,inew)=ply(i,j,iold)+(addpy)*dt 
        cdif(i,j,inew)=cp(i,j,inew)-cm(i,j,inew)
        da1(i,j)=x(i,j,inew)
        da2(i,j)=cp(i,j,inew)
        da3(i,j)=cm(i,j,inew)
        da4(i,j)=plx(i,j,inew)
        da5(i,j)=ply(i,j,inew)
        enddo
        enddo

        
         
c     Save the updated fields at the scheduled timestep.
c     The first filename used is names(101).
          if(itime.eq.next_save) then

          ndef=ndef+1
          ndef1=ndef1+1

          if (ndef1.gt.900) then
          write(*,*) 'Snapshot filename list exhausted'
          stop 1
          endif
          open(unit=22,file=names(ndef1),status='replace')


          do i=1,k
          do j=1,k


c     Output columns: i, j, psi, c_plus, c_minus, p_x, p_y.
          write(22,*)i,j,da1(i,j),da2(i,j),da3(i,j),da4(i,j),da5(i,j)

          enddo

          write(22,*)""
          enddo
          close(22)
         write(*,*) itime

         next_save = next_save*2

         endif


       
c     Swap buffers only after every grid cell has been updated.
        itemp=iold
        iold=inew
        inew=itemp

        
       
        enddo 

        close(1)
        close(2)
        close(3)
        close(4)
        close(11)
        close(12)
        stop
        end 
cccccccccccccccccccccccccccccccccccccccccc
c     Random-number generator used for the initial perturbations.
       FUNCTION ran2(idum)
       INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
       REAL*8 ran2,AM,EPS,RNMX
       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     * NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
       INTEGER idum2,j,k,iv(NTAB),iy
       SAVE iv,iy,idum2
       DATA idum2/123456789/, iv/NTAB*0/, iy/0/
       if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           if (idum.lt.0) idum=idum+IM1
           if (j.le.NTAB) iv(j)=idum
11       continue
         iy=iv(1)
       endif
       k=idum/IQ1
       idum=IA1*(idum-k*IQ1)-k*IR1
       if (idum.lt.0) idum=idum+IM1
       k=idum2/IQ2
       idum2=IA2*(idum2-k*IQ2)-k*IR2
       if (idum2.lt.0) idum2=idum2+IM2
       j=1+iy/NDIV
       iy=iv(j)-idum2
       iv(j)=idum
       if(iy.lt.1)iy=iy+IMM1
       ran2=min(AM*iy,RNMX)

       return
       END
        

       






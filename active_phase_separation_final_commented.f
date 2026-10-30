c     Two-dimensional integration of Eqs. S10a-c.
c     Fields: x=psi, cp=c_plus, cm=c_minus, plx=p_x, ply=p_y.
c     Only explanatory comments have been added.
c     Corrected periodic version: staged intermediates, current vm.
c  first code for one dimensional
c     Import the function used to detect infinite or NaN values.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
       implicit real*8(a-h,o-z)
c     Declare filename arrays, grid size and field arrays.
c     The last field index stores the old and new time levels.
      character*8 names(901),name1(901),name2(900),name3(900)
      parameter(k=340,k2=k**2,khalf=k/2,n=1)
      
cc    k refers to the system size

c     p(i) and m(i) are periodic neighbor indices, not polarity.
      integer m(k),p(k),s1,is0
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

c     Store the local velocity at each site for the current step.
      DOUBLE PRECISION vmfield(k,k),motors
      integer next_save
      next_save = 10000

c     Set the random-number seed variables for initialization.
       iseed=-137571
       iseed1=-137573
c       v0=-delv0
c       v1=0.03
c       nv0=1

c     Open input filename lists and output files.
c     Read 900 filename entries from each input list.
              open(unit=1,file='fnamerhol')
              open(unit=2,file='fnamerholp')
              open(unit=3,file='fnamerholp5')
              open(unit=4,file='fnamerholp8')

ccc           to save all data in a file
              open(unit=11,file='data.dat')
              open(unit=12,file='data_comp_polarity.dat')

              read(1,*)(names(i),i=1,900)
              read(2,*)(name1(i),i=1,900)
              read(3,*)(name2(i),i=1,900)
              read(4,*)(name3(i),i=1,900)

c     Set pi, the time increment and the spatial increments.
c     dx2=2*dx is used for first derivatives; dxx=dx**2 for
c     second derivatives.
       pi=2.0d0*asin(1.0d0)
       dt=0.000075d0
       dx=0.015d0
       dx2=2*dx
       dxx=dx**2

ccc   average values of different species

c     Set nominal initial mean fractions and mean polarity.
      avgpsi=0.2
      avgcp=0.15
      avgcm=0.15
      avgp=0.0      



c     Set the maximum number of timesteps.
       ntime1=0
       ntime=200000000
cccc    defining neighbor latticle points, important for finte
ccc     difference method
c     Construct forward and backward neighbors with periodic
c     wrapping. The same maps are used in both directions.
          do i=1,k
          p(i)=i+1
          if(i.eq.k)p(i)=1
          m(i)=i-1
          if(i.eq.1)m(i)=k
          enddo
cccc    starting defining initial conditions
c     Initialize all fields on the two-dimensional grid.
c     psi: uniform perturbation of amplitude +/-0.025.
c     cp and cm: uniform perturbations of amplitude +/-0.005.
c     Each polarity component: perturbation of +/-0.025.
        do j=1,k
        do i=1,k

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
c     Initialize the motor difference cdif=cp-cm.
        cdif(i,j,1)=cp(i,j,1)-cm(i,j,1)
c       random angle        
c     Generate an angle; the active initialization below uses
c     separate random draws for plx and ply instead.
        theta=(ran2(iseed)-0.5)*2*pi
c       polarization components        
        plx(i,j,1)=avgp+0.05*(0.5-ran2(iseed))
c cos(theta)
        ply(i,j,1)=avgp+0.05*(0.5-ran2(iseed))
c sin(theta)
        enddo
        enddo


c     Initialize counters used for selecting output filenames.
        nrho=100
        npic=100
        ndef=100
        ndef1=100
        npro=100

ccc    parameters
c     Assign diffusion, speed, geometry and coupling parameters.
c     lsq=2*Le is twice the microtubule length.
        Dr    = 0.1D0
        lp    = 5.0D0
        v0    = 0.25D0
        alpha = 1000.0D0
        Le    = 5.0D0
        lsq   = 2.0D0*Le
        d     = 0.02D0
        bet   = 1.0D0
        anu   = 0.007D0 
        Dpsi  = 0.1D0
        Dc    = 1.0D0
        azeta = 0.0000063D0
        arho  = 1.0D0/(pi*0.02D0*0.02D0)
        akap   = 0.00005D0

        
c     Calculate coefficients used in the reduced equations.
c     Dmt and Dm are reduced diffusion coefficients.
c     vp is the reduced motor speed; ara is the inverse duty ratio.
c     anup is the coupling coefficient used in S10.
c     azetap is 2*zeta_tilde/(mean_cp+mean_cm).
c     akap is kappa_tilde/mean_psi.
        Dmt=(Dpsi)/(Dr*lp*lp)
        Dm=(Dc)/(Dr*lp*lp)
        anup=anu*arho/(lp)
        azetap=azeta*arho/(lp*lp)
        azetapp=azetap*2*avgcp/(avgcp+avgcm)
        azetapm=azetap*2*avgcm/(avgcp+avgcm)
        vp=v0/(Dr*lp)
        ara=Dc*DLOG(Le/d)/(pi*d*v0)

        
c        write(*,*)vp,Dmt,Dm,anup,azetap,akap

c     Select the initial old/new buffers and start time evolution.
        iold=1
        inew=2
        itp=0
        do itime=1,ntime

c     Write the parameter row once, at the first timestep.
        if(itime.eq.1) then
        write(12,*)vp,Dmt,Dm,anup,azetap,akap,ara
        endif

c     Validate physical denominators; do not clip conserved fields.
c     STEP 1: Check old fields before calculating derivatives.
c     Stop for invalid fractions or non-finite tested values.
c     No clipping of the fields is performed.
        do j=1,k
        do i=1,k
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
c     Calculate vm=vp*(cp-cm)/(cp+cm) locally.
c     The code convention is vm=-vmt_tilde in Eqs. S10a-c.
c     Use zero when the total motor fraction vanishes.
        motors=cp(i,j,iold)+cm(i,j,iold)
        vmfield(i,j)=0.0d0
        if (motors.gt.0.0d0) then
        vmfield(i,j)=vp*(cp(i,j,iold)-cm(i,j,iold))/motors
        endif
        enddo
        enddo

c     STEP 2: Build the first set of intermediate fields over
c     the full grid, using only the old field values.
c     Full-grid intermediate pass 1
        do j=1,k
        do i=1,k
c     pxx and pyy are first derivatives of polarity.
c     Their sum is q=div(p). Store pspx=psi**2*q.
        pxx=(plx(p(i),j,iold)-plx(m(i),j,iold))/dx2
        pyy=(ply(i,p(j),iold)-ply(i,m(j),iold))/dx2
        pspx(i,j)=(pxx+pyy)*x(i,j,iold)*x(i,j,iold)
        polOP=plx(i,j,iold)*plx(i,j,iold)+ply(i,j,iold)*ply(i,j,iold)
        vm=vmfield(i,j)
        cdiffx=(cdif(p(i),j,iold)-cdif(m(i),j,iold))/dx2
        cdiffy=(cdif(i,p(j),iold)-cdif(i,m(j),iold))/dx2
c     Store the vector vm*psi*p for microtubule transport.
         conpx(i,j)=vm*x(i,j,iold)*plx(i,j,iold)
         conpy(i,j)=vm*x(i,j,iold)*ply(i,j,iold)
c     agrp2=A=|p|**2*(1+alpha*|p|**2/2).
c     dpcdf=B=p dot grad(cp-cm).
        agrp2(i,j)=polOP*(1.0+(alpha/2)*polOP)
        dpcdf(i,j)=cdiffx*plx(i,j,iold)+cdiffy*ply(i,j,iold)
         apOP(i,j)=sqrt(plx(i,j,iold)**2+ply(i,j,iold)**2)    
c     Calculate the Laplacian and gradient of psi.
c     agdpx and agdpy store the excluded-volume vector:
c     (2*Le/(pi*d))*psi/(1-psi-cp-cm)*grad(psi).
        psixx=(x(p(i),j,iold)-2*x(i,j,iold)+x(m(i),j,iold))/dxx
        psiyy=(x(i,p(j),iold)-2*x(i,j,iold)+x(i,m(j),iold))/dxx
        alapsi(i,j)=(psixx+psiyy)
        agdpsx=(x(p(i),j,iold)-x(m(i),j,iold))/dx2
        agdpsy=(x(i,p(j),iold)-x(i,m(j),iold))/dx2
        atot=x(i,j,iold)+cp(i,j,iold)+cm(i,j,iold)
        agdpx(i,j)=(lsq/(pi*d))*(x(i,j,iold)/(1.0-atot))*(agdpsx)
        agdpy(i,j)=(lsq/(pi*d))*(x(i,j,iold)/(1.0-atot))*(agdpsy)
c     aggp2=H=psi*(div(p))**2 is the density bending scalar.
c     bpsx and bpsy are the components of psi*p.
        aggp2(i,j)=(pxx+pyy)*(pxx+pyy)*x(i,j,iold)
              bpsx(i,j)=plx(i,j,iold)*x(i,j,iold)
              bpsy(i,j)=ply(i,j,iold)*x(i,j,iold)
c     Store the active motor vectors inside the RHS divergence:
c     plus: ara*(-vp+vm)*psi*cp*p;
c     minus: ara*(vp+vm)*psi*cm*p.
         conpxp(i,j)=ara*(-vp+vm)*x(i,j,iold)*cp(i,j,iold)*plx(i,j,iold)
         conpyp(i,j)=ara*(-vp+vm)*x(i,j,iold)*cp(i,j,iold)*ply(i,j,iold)
         conmxm(i,j)=ara*(vp+vm)*x(i,j,iold)*cm(i,j,iold)*plx(i,j,iold)
         conmym(i,j)=ara*(vp+vm)*x(i,j,iold)*cm(i,j,iold)*ply(i,j,iold)
c     Calculate motor gradients and excluded-volume vectors:
c     cp/(1-psi-cp-cm)*grad(cp), and similarly for cm.
        agdcpx=(cp(p(i),j,iold)-cp(m(i),j,iold))/dx2
        agdcpy=(cp(i,p(j),iold)-cp(i,m(j),iold))/dx2
        agdcpxx(i,j)=(cp(i,j,iold)/(1.0-atot))*(agdcpx)
        agdcpyy(i,j)=(cp(i,j,iold)/(1.0-atot))*(agdcpy)
         agdcmx=(cm(p(i),j,iold)-cm(m(i),j,iold))/dx2
         agdcmy=(cm(i,p(j),iold)-cm(i,m(j),iold))/dx2
        agdcmxx(i,j)=(cm(i,j,iold)/(1.0-atot))*(agdcmx)
        agdcmyy(i,j)=(cm(i,j,iold)/(1.0-atot))*(agdcmy)
        enddo
        enddo

c     STEP 3: Differentiate the completed first-pass arrays.
c     Full-grid intermediate pass 2
        do j=1,k
        do i=1,k
c     Store psi*grad(A), psi*grad(B), and grad(H).
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
c     Calculate motor Laplacians and div(psi*p).
c     Store the gradient-energy and torque potential parts:
c     amucp=-azetap*lap(cp)+anup*pi*d/(2*Le)*div(psi*p);
c     amucm=-azetap*lap(cm)-anup*pi*d/(2*Le)*div(psi*p).
         acpxx=(cp(p(i),j,iold)-2*cp(i,j,iold)+cp(m(i),j,iold))/dxx
         acpyy=(cp(i,p(j),iold)-2*cp(i,j,iold)+cp(i,m(j),iold))/dxx
         acmxx=(cm(p(i),j,iold)-2*cm(i,j,iold)+cm(m(i),j,iold))/dxx
         acmyy=(cm(i,p(j),iold)-2*cm(i,j,iold)+cm(i,m(j),iold))/dxx
            bpsxx=(bpsx(p(i),j)-bpsx(m(i),j))/dx2
            bpsyy=(bpsy(i,p(j))-bpsy(i,m(j)))/dx2
         amucp(i,j)=-azetap*(acpxx+acpyy)+anup*(pi*d/lsq)*(bpsxx+bpsyy)
         amucm(i,j)=-azetap*(acmxx+acmyy)-anup*(pi*d/lsq)*(bpsxx+bpsyy)
        enddo
        enddo

c     STEP 4: Build psi*grad(H) and the motor potential vectors.
c     ajcp=cp*grad(amucp); ajcm=cm*grad(amucm).
c     Complete these arrays before calculating their divergence.
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

c     STEP 5: Assemble all right-hand sides from old fields
c     and the completed intermediate arrays.
c     All intermediate arrays now belong to this timestep.
        do j=1,k
        do i=1,k
c     Calculate both spatial derivatives of each polarity
c     component and grad(psi**2*div(p)).
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
c     Polarity equation S10c, x component:
c     tpbetax: nonlinear relaxation multiplied by psi;
c     tnux/bet: motor-gradient coupling multiplied by psi;
c     tpkapx: bending contribution before division by psi.
c     advpx=-vm*(p dot grad)*p_x. Subtracting advpx*psi
c     and then dividing by psi gives +vm*(p dot grad)*p_x.
c     Finally divide the assembled expression by psi.
        tpbetax=-plx(i,j,iold)*(1.0+alpha*polOP)*x(i,j,iold)
        tnux=anup*cdiffx*x(i,j,iold)
        tpkapx=akap*(pxpx)
        advpx=-vm*plx(i,j,iold)*pxx-vm*ply(i,j,iold)*pxxy
        addpx=tpbetax+(tnux/bet)+tpkapx-advpx*x(i,j,iold)
        addpx=addpx/x(i,j,iold)
c     Assemble the same polarity terms for the y component.
        tpbetay=-ply(i,j,iold)*(1.0+alpha*polOP)*x(i,j,iold)
        tnuy=anup*cdiffy*x(i,j,iold)
        tpkapy=akap*(pypy)
        advpy=-vm*ply(i,j,iold)*pyy-vm*plx(i,j,iold)*pyyx
        addpy=tpbetay+(tnuy/bet)+tpkapy-advpy*x(i,j,iold)
        addpy=addpy/x(i,j,iold)
c     Microtubule active transport: addcompsi=div(vm*psi*p).
         conpxx=(conpx(p(i),j)-conpx(m(i),j))/dx2
         conpyy=(conpy(i,p(j))-conpy(i,m(j)))/dx2
        addcompsi=(conpxx+conpyy)
c     Calculate divergences of the prepared vectors for psi:
c     div(psi*grad(A)), div(psi*grad(B)), the excluded-volume
c     vector, and div(psi*grad(H)).
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
c     Microtubule equation S10b:
c     tlapsi: ordinary diffusion plus excluded-volume correction;
c     tbetpsi: local polarity energy;
c     tnupsict: negative torque-coupling contribution;
c     tkapsi: bending contribution.
c     Multiply their sum by Dmt to obtain addpsitotal.
        tlapsi=(psixx+psiyy+(agdpsxx+agdpsyy))
        tbetpsi=(bet/2)*(divspx2+divspy2)
        tnupsict=-anup*(dinupsx+dinupsy)
        tkapsi=bet*akap*(daggpx2+daggpy2)
        addpsitotal=(tlapsi+tbetpsi+tnupsict+tkapsi)*Dmt
c     Motor equation S10a: calculate ordinary Laplacians,
c     divergences of ajcp/ajcm, active motor-vector divergences,
c     and excluded-volume-vector divergences.
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
c     Assemble complete motor right-hand sides.
c     The active divergence is outside the Dm factor.
c     Dm multiplies ordinary diffusion, potential-gradient
c     relaxation and the excluded-volume correction.
       addtotcp=conpxxp+conpyyp+Dm*(acpxx+acpyy+divcpx+divcpy)
       addtotcp=addtotcp+Dm*(agdcpsxxx+agdcpsyyy)  
       addtotcm=conmxxm+conmyym+Dm*(acmxx+acmyy+divcmx+divcmy)
       addtotcm=addtotcm+Dm*(agdcmsxxx+agdcmsyyy)
c     STEP 6: Update all five fields by explicit forward Euler.
c     new=old+dt*RHS. Recalculate cdif from the new motor fields.
         x(i,j,inew)=x(i,j,iold)+(addpsitotal+addcompsi)*dt
        cp(i,j,inew)=cp(i,j,iold)+addtotcp*dt
        cm(i,j,inew)=cm(i,j,iold)+addtotcm*dt
        plx(i,j,inew)=plx(i,j,iold)+(addpx)*dt
        ply(i,j,inew)=ply(i,j,iold)+(addpy)*dt 
        cdif(i,j,inew)=cp(i,j,inew)-cm(i,j,inew)
c     Copy updated fields into arrays used for snapshot output.
        da1(i,j)=x(i,j,inew)
        da2(i,j)=cp(i,j,inew)
        da3(i,j)=cm(i,j,inew)
        da4(i,j)=plx(i,j,inew)
        da5(i,j)=ply(i,j,inew)
        enddo
        enddo

        
         
c     STEP 7: Save when itime reaches next_save.
c     Increment the filename counter and check list capacity.
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


c     Write columns i,j,psi,cp,cm,p_x,p_y.
c     Separate fixed-i blocks by a blank line.
          write(22,*)i,j,da1(i,j),da2(i,j),da3(i,j),da4(i,j),da5(i,j)

          enddo

          write(22,*)""
          enddo
          close(22)
         write(*,*) itime

c     Double the next output step: 10000,20000,40000,...
         next_save = next_save*2

         endif


       
c     STEP 8: Exchange old and new buffer indices for the
c     next timestep without copying the field arrays.
        itemp=iold
        iold=inew
        inew=itemp

        
       
        enddo 

c     Close files after the time loop has finished.
        close(1)
        close(2)
        close(3)
        close(4)
        close(11)
        close(12)
        stop
        end 
cccccccccccccccccccccccccccccccccccccccccc
c     Uniform random-number generator used for initial fields.
c     Internal generator state is saved between calls.
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
        

       






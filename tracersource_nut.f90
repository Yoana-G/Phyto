!subroutine tracersource(m,n,dtimel) 
      subroutine tracerinit(stepl)
  !     -----------------------  
  !                                         
      USE header
      integer i,j,k,l,n,step1
      double precision z,mldensity,sig,trseed

  !     tracer sources                                                    
  !   Tracer 1 is nitrate

  !   if this is a continuation of a run, only trinit is to be initialized
  !   Tr 1 is nutrient

  !    phydepth= 30.d0
     trseed=1.d-1    ! Chl =0.08 mg/l, Chl2C=15
mldensity = 1026.5 ! nitracline density

!     TRACERS
!     =========

!     Initialze nutrient
      do j=0,NJ+1
         do i=0,NI+1
            do k=NK,1,-1
               !z = -zc(i,j,k)*DL
               sig = rho(i,j,k)
               Tr(1,i,j,k,0) = (sig-mldensity)*6.0
               if (sig.lt.mldensity) Tr(1,i,j,k,0)=trseed
               if (sig.gt.1027.5) Tr(1,i,j,k,0) = (sig-1027.5)*2+6
               Tr(2,i,j,k,0) = 1.d-5
			   Tr(3,i,j,k,0) = Tr(2,i,j,k,0)
               Tr(4,i,j,k,0) = Tr(2,i,j,k,0)
            end do
         end do
      end do

!     Periodic Boundaries (periodic e-w)
      do n=0,1
         do k=0,NK+1
            do j=0,NJ+1
               do it=1,ntr
                  Tr(it,0,j,k,n)= Tr(it,NI,j,k,n)
                  Tr(it,NI+1,j,k,n)= Tr(it,1,j,k,n)
               end do
            end do
         end do
      end do
      return
      end
                                                        
  integer  i,j,k,n,m 
  REAL(kind=rc_kind) :: dtimel,fac,tau,tauinv,vol,r,GS,GR,GE,kpar,I_0,L_PS,L_PR,L_PE,mPS,mPR,mPE,k_PR,k_PS,k_PE,I
  !                                                                       
  !     consump is the uptake (or if negative, then addition) of tracer   

  tau is the time scale for consump = 1day = 86400 seconds
  tauinv = 1.d0/(1.d0*86400.d0)    ! tauinv is in (per second)
  fac= dtimel/(1.d0*86400.d0*UL/LEN) ! dtime is non-dim time, TL=LEN/UL
  ! fac is non-dim

!note indices --> Tr 1: Nutrients, Tr 2: Pro, Tr 3: Syn, Tr 4: Euks

  kpar = 0.05; ! light decay
  GS = 5;  !max growth rate for Synechococcus
  GR = 1.5; !max growth rate for Prochlorococcus
  GE = 5;  !max growth rate for Picoeukaryotes
  I_0 = 1;   !light intensity
  L_PS = 5; !gamma;light dependence  (Syn)
  L_PR = 250; !light dependence (Pro)
  L_PE = 5;  !light dependence (euks)
  k_PR = 0.25; !half-saturation parameter; 0.01 - 5 "realistic range" / (Pro)
  k_PS = 0.5;  !half-saturation parameter for (Syn)  
  k_PE = 0.5;  !half-saturation parameter for (euks)
  mPS = 0.1; !mortality rate (Syn) 
  mPR = 0.1; !mortality rate (Pro)
  mPE = 0.1; !mortality rate (euks)

  !W = 1; !sinking rate
  !PE0 = 0; ! bottom boundary condition for Euks
                                                             
  do m=0,1 
     do k=0,NK+1 
      !do k=1,NK
        do j=0,NJ+1
         !do j=1,NJ 
           do it=1,ntr
	    !do i=1,NI 
			   z = zc(i,j,k)*DL
			   I = I_0*EXP(kpar*z)
              		Tr(it,0,j,k,m)= Tr(it,NI,j,k,m) 
              		Tr(it,NI+1,j,k,m)= Tr(it,1,j,k,m) 

	       !starting biological reactions
            	Tr(1,i,j,k,n) = Tr(1,i,j,k,m)+fac*(-GR*(EXP(-L_PR*I))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+k_PR)*Tr(2,i,j,k,m)+mPR*Tr(2,i,j,k,m)) + &
				fac*(-GS*(EXP(-L_PS*I))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+k_PS)*Tr(3,i,j,k,m)-mPS*Tr(3,i,j,k,m)) + &
				fac*(-GE*(EXP(-L_PE*I))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+k_PE)*Tr(4,i,j,k,m)-mPE*Tr(4,i,j,k,m))
            	Tr(2,i,j,k,n) = Tr(2,i,j,k,m)+fac*(-GR*(EXP(-L_PR*I))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+k_PR)*Tr(2,i,j,k,m)-mPR*Tr(2,i,j,k,m))
		Tr(3,i,j,k,n) = Tr(3,i,j,k,m)+fac*(-GS*(EXP(-L_PS*I))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+k_PS)*Tr(3,i,j,k,m)-mPS*Tr(3,i,j,k,m))
            	Tr(4,i,j,k,n) = Tr(4,i,j,k,m)+fac*(-GE*(EXP(-L_PE*I))*Tr(1,i,j,k,m)/(Tr(1,i,j,k,m)+k_PE)*Tr(4,i,j,k,m)-mPE*Tr(4,i,j,k,m))

	       
  	       !N = N + dN_dt*dt
  	       !PS = PS + dPS_dt*dt
 	       !PR = PR + dPR_dt*dt
 	       !PE = PE + dPE_dt*dt

           end do
        end do
     end do
  end do
  !                                                                       
  return 
end subroutine tracersource

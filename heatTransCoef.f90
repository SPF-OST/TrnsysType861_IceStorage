!---------------------------------------------------------------------------  
!> @file heatTransCoef.f90
!> @author SPF institut fur SolarTechnik 
!> @author D.Carbonell, W.Lodge
!> @date 10/08/2012
!> @brief convection heat trasnfer coefficients related functions
!--------------------------------------------------------------------------  
    
module heatTransferCoef

    use physProp
    use TrnsysFunctions
    use TrnsysConstants
    
    implicit none
    
    integer,parameter :: LAMINAR=0,&
            TURBULENT=1
    
    integer,parameter :: CIRCULAR=0,&
               RECTANGULAR=1
    
    contains

    double precision function getHydraulicDiameter(thick,lenght)
    
        double precision, intent(in)  :: thick, lenght
        double precision :: area,P
        
        area = thick * lenght                ! : Approx. cross section of one HX
        P = (2 * lenght) + (2 * thick)       ! : Wetted perimeter of HX
        getHydraulicDiameter = 4 * area / P  
        
    end function getHydraulicDiameter
    
    subroutine calculateHeatTransCoefPipeInWong(tBulk,conc,D,L,mIn,alphaIn,thick,crossSectionType,enhancedNu)
        use spfGlobalConst, only: PI
    
        implicit none
        
        double precision, intent(in)  :: tBulk, conc, D, L, mIn, thick, enhancedNu
        integer,intent(in) :: crossSectionType
        double precision, intent(out) :: alphaIn
        double precision :: Nu, Pr, Re, vIn, muFluid,cpFluid,lambdaFluid,rhoFluid, Dh        
        
        cpFluid      = getCpEg(conc,tBulk);
        lambdaFluid  = getLambdaEg(conc,tBulk);
        rhoFluid     = getRhoEg(conc,tBulk);
        muFluid      = getNuEg(conc,tBulk)*rhoFluid; ! be carefull we need dynamic viscosity, not kinematic !!

        ! double lengthChar = diameter*PI;

        if(crossSectionType==RECTANGULAR) then
            if(L/thick<=1.0) then
               Dh = thick                                        
            else if(L/thick<=1.4d0) then
               Dh = 1.17*thick 
            else if(L/thick<=2.0d0) then
               Dh = 1.33*thick 
            else if(L/thick<=3.0d0) then
               Dh = 1.5*thick 
            else if(L/thick<=4.0d0) then
               Dh = 1.6*thick 
            else if(L/thick<=8.0d0) then
               Dh = 1.78*thick 
            else 
               Dh = 2.0*thick 
            endif    
                       
            vIn = mIn/(rhoFluid*thick*L);
            
        else
            Dh = D
            vIn = 4.*mIn/(rhoFluid*pi*Dh*Dh);
        endif
    
        Pr  = getPrandtl (muFluid,cpFluid,lambdaFluid);
        
        Re  = getReynolds (rhoFluid,vIn,Dh,muFluid);
        
                                  
        call NuWongForcedConvection(Re,Pr,Dh,L,muFluid,muFluid,Nu,thick,crossSectionType);
    
        alphaIn = Nu*enhancedNu*lambdaFluid/Dh;
        
        !write(20,*) 'T=',tBulk,' Pr=',Pr,' vIn=',vIn,' mIn=',mIn,' Re=',Re,' alphaIn=',alphaIn,' mu=',muFluid,' cp=',cpFluid,' lambda=',lambdaFluid,' rho^=',rhoFluid
     
    endsubroutine calculateHeatTransCoefPipeInWong
    

    subroutine calculateHeatTransCoefPipeIn(tBulk,conc,thick,lenght,vel,hIn,flatPlate,enhancedNu,iUnit,iType,Re,Nu_i)
            
        use spfGlobalConst, only: PI
              
        implicit none
        
        integer, intent(in) :: iUnit,iType
        double precision, intent(in)  :: tBulk, conc, thick, lenght, vel, enhancedNu 
        double precision, intent(out) :: hIn, Re, Nu_i
        
        double precision :: cp_b, k_b, nu_b, rho_b, alpha_b, Pr_b, P, D_h,  f,  area,tAv       
        logical :: flatPlate
        character (len=maxMessageLength) :: MyMessage
       
        if(conc==0) then
            tAv = min(max(tBulk,1.0),99.)
        else
            tAv = tBulk !DC big error
        endif
        
        ! glycol heat exchanger fluid properties
        
        cp_b    = getCpEg(conc, tAv)
        k_b     = getLambdaEg(conc, tAv)
        nu_b    = getNuEg(conc, tAv)
        rho_b   = getRhoEg(conc, tAv)
        alpha_b = k_b / (rho_b * cp_b)
        Pr_b    = nu_b / alpha_b
      
    !     Internal forced convection:
    
        if(flatPlate==.true.) then            
            area = thick * lenght                  ! : Approx. cross section of one HX
            P = (2 * lenght) + (2 * thick)         ! : Wetted perimeter of HX
            D_h = 4 * area / P                     ! : Dydraulic diameter of HX
        else
            !D_h = pi*thick*thick/4. Why I had this???
            D_h = thick
        endif    
        
        !u_m = (mDot / n_hx) / (rho_hx * A_c)  ! : Mean fluid velocity
        Re = vel * D_h / nu_b                  ! : Reynolds number
        
        if (Re<2300d0) then                ! : Friction coeff. Funadamentals of heat and mass transfer pag 490
            !f = 96.0d0 / Re ! Why 96??? it shold be 64 Eq. 8.19 !!
            f = 64.0d0 / Re !fully developed laminar flow
        elseif (Re<2e4) then
            f = 0.316d0 * Re**(-0.25d0)
        else
            f = 0.184d0 * Re**(-0.2d0)
        endif                
        
        ! Fundamentals of heat and mass transfer pag 515 Nu_D laminar  = 4.36  (qConst) and if Tconst Nu_d laminar = 3.66
      
        if (Re<2300d0) then
            
            if(enhancedNu>=0) then
                Nu_i = 4.36*enhancedNu !8.23*enhancedNu    ! : Nusselt number (internal)
            else 
                Nu_i = 4.36
            endif       
        else
            !TURBULENT FLOW IN CICRULAR TUBES Gnielinski Fundamentals of heat and mass transfer pag 515
            Nu_i = (f/8d0) * (Re - 1000d0) * Pr_b / (1 + (12.7d0 * (f/8d0)**0.5d0 * (Pr_b**(2d0/3d0) - 1)))
        endif                
        
        hIn = Nu_i * k_b / D_h   ! : Convection coeff. (internal)           
        
        if(hIn<1e-10) then
            D_h = D_h + 1e-30
        endif
        
        if(isnan(hIn) ) then                      
            write(MyMessage,'("calculateHeatTransCoefPipeIn: FATAL ERROR : hIn NAN. tAv= = ",f)') tAv            
            call Messages(-1,Trim(MyMessage),'FATAL', iUnit,iType)
        endif     
        
    end subroutine calculateHeatTransCoefPipeIn
  
     subroutine calculateHeatTransCoefFlatPlate(tBulk,conc,thick,width,lenght,vel,hIn,enhancedNu,iUnit,iType,Re,Nu)
            
        use spfGlobalConst, only: PI
        use Trnsysfunctions
        use interpolation
        
        implicit none
        
        integer, intent(in) :: iUnit,iType
        character (len=maxMessageLength) :: MyMessage
        double precision, intent(in)  :: tBulk, conc, thick, width, lenght, vel, enhancedNu 
        double precision, intent(out) :: hIn,Re,Nu        
        double precision :: cpFluid, lambdaFluid, rhoFluid, muFluid, nuFluid, Dh, P, Pr, area,tAv              
        double precision :: K,n,C,m,ReTransLow,ReTransHigh,NuLam,NuTur
        
        
        if(conc==0) then
            tAv = min(max(tBulk,1.0),99.)
        else
            tAv = tBulk 
        endif
        
        ! glycol heat exchanger fluid properties
        
        cpFluid      = getCpEg(conc,tBulk);
        lambdaFluid  = getLambdaEg(conc,tBulk);
        rhoFluid     = getRhoEg(conc,tBulk);
        nuFluid      = getNuEg(conc,tBulk)
        muFluid      = nuFluid*rhoFluid; ! Dynamic viscosity, not kinematic !!        
      
        ! Internal forced convection:
                    
        area = thick * width                  ! : Approx. cross section of one HX. Thick is b corrugated
        P = (2 * width) + (2 * thick)         ! : Wetted perimeter of HX
        Dh = 4 * area / P                     ! : Hydraulic diameter of HX. Approx D_H = 2b
                                
        Pr  = getPrandtl (muFluid,cpFluid,lambdaFluid)        
        Re  = getReynoldsNu (vel,Dh,nuFluid)                
        
        ReTransLow  = 70d0  !10d0
        ReTransHigh = 150d0 !150d0
         
        ! Laminar Nu
        K=1.0
        n = 0.4         
        C = 1.68            
        !Nu = C*(Re*Pr*Dh/lenght)**n; 
        NuLam = C*(Re*Pr*Dh/width)**n; 
        
        ! Turbulent Nu
                
        K=1.0
        n = 0.4      
        C = 0.2   ! 0.1528
        m = 0.67  ! 0.66    
        NuTur = C*(Re**m)*(Pr**n)*K;  
        
        if(Re<=ReTransLow) then
            Nu = NuLam
        elseif(Re>=ReTransHigh) then
            Nu = NuTur
        else
            Nu = linear(ReTransLow,NuLam,ReTransHigh,NuTur,Re)
        endif
                               
        !Nu = max(Nu,15d0)
        
        hIn = enhancedNu * Nu * lambdaFluid  / Dh   ! : Convection coeff. (internal)                          
                 
        if(isnan(hIn) ) then                      
            write(MyMessage,'("calculateHeatTransCoefFlatPlate: FATAL ERROR : hIn NAN. tAv= = ",f)') tAv            
            call Messages(-1,Trim(MyMessage),'FATAL', iUnit,iType)
        endif  
        
        
     end subroutine calculateHeatTransCoefFlatPlate
     
    subroutine calculateHeatTransCoefPipeOut(tFilm,tInf,tWall,lenght,hOut)
    
        use spfGlobalConst, only: g
        implicit none
        
        double precision, intent(in)  :: lenght,tFilm,tInf,tWall
        double precision, intent(out) :: hOut
        double precision :: cp_f, k_f, rho_f, nu_f, beta_f, Pr_f, alpha_f, Ra, Nu_o,&
        mu, myTFilm
    
        myTFilm = max(tFilm,0.1)
        
        cp_f    = getCpWater(myTFilm)            ! : specific heat
        k_f     = getLambdaWater(myTFilm)        ! : conductivity
        rho_f   = getRhoWater(myTFilm)           ! : density
        nu_f    = getNuWater(myTFilm)            ! : kinematic viscosity
        beta_f  = getBetaWater(myTFilm)          ! : volumetric expansion coeff.
        Pr_f    = getPrandtlWater(myTFilm)       ! : Prandtl number
        alpha_f = k_f / (rho_f * cp_f)         ! : thermal diffusivity
        
        Ra = g * beta_f * abs(tWall - tInf) *  lenght**3 / (nu_f * alpha_f)
        
    ! mu = nu_f*rho_f        
    ! Ra = getRayleigh  (beta_f,rho_f,abs(tWall - tInf),lenght,mu,cp_f,k_f)
        
    !   Correlation from Churchill and Chu (see: Incropera and DeWitt)
    
        !this produced inestabilities in the algorithm
        !if(Ra<=1e9) then !LAMINAR
            !Nu_o = 0.68 + ((0.670 * Ra**(1./4.)) / (1 + (0.492/Pr_f)**(9./16.))**(4./9.))
        !else    !TURBULENT
            !Nu_o = (0.825 + ((0.387 * Ra**(1./6.)) / (1 + (0.492/Pr_f)**(9./16.))**(8./27.)))**2. ! pag 571 Eq:9.26
        !endif
        
        !FOR VERTICAL PLATE
        Nu_o = (0.825 + ((0.387 * Ra**(1./6.)) / (1 + (0.492/Pr_f)**(9./16.))**(8./27.)))**2. ! pag 571 Eq:9.26
        !Nu_o = 0.68 + ((0.670 * Ra**(1./4.)) / (1 + (0.492/Pr_f)**(9./16.))**(4./9.)) !laminar Ra<=10^9
        
        hOut = Nu_o * k_f / lenght  ! : Convection coeff. (external)            
    
        if(isnan(hOut)) then
            hOut = Nu_o * k_f / lenght
        end if 
        
    end subroutine calculateHeatTransCoefPipeOut    
    
    !For vertical plate
  
    subroutine calculateHeatTransCoefImmersedPlate(tFilm,tFluid,tWall,lChar,hOut,cUserDefined,nUserDefined,Ra,Nu,iUnit,iType)
    
        use spfGlobalConst, only: PI
        use interpolation
        implicit none
        
        integer, intent(in) :: iUnit,iType
        character (len=maxMessageLength) :: MyMessage
        double precision, intent(in)  :: lChar,tFilm,tFluid,tWall
        double precision, intent(out) :: hOut,Ra,Nu
        double precision, intent(in) :: cUserDefined,nUserDefined
        double precision :: cp, lambda, rho, mu, beta, deltaT, Gr, Pr,&
        c,n, PrFunction,NuLam,NuTur
    
        cp    = getCpWater(tFilm)            ! : specific heat
        lambda = getLambdaWater(tFilm)       ! : conductivity
        rho   = getRhoWater(tFilm)           ! : density
        mu    = getNuWater(tFilm)*rho        ! : dynamic viscosity
        
        beta  = getBetaWater(tFilm)          ! : volumetric expansion coeff.
        !beta  = 1./(tFilm+273.15)           ! : Only valid for ideal gas !!
        deltaT = max(abs(tFluid-tWall),0.1)
        
        Ra = getRayleigh  (beta,rho,deltaT,lChar,mu,cp,lambda,iUnit,iType)               
       
        if(cUserDefined>0.0 .and. nUserDefined>0.0) then
            
            c = cUserDefined
            n = nUserDefined        
            Nu = c*Ra**n
        
        else 
            
            Pr  = getPrandtl (mu,cp,lambda)
            
            !heat exchanger design handbook. chapter 2.5.7 
            PrFunction = 1.d0/(1.0d0 + (0.492/Pr)**(9.d0/16.))**(16.d0/9.)                
            
            !churchill and chu from incroepra Eq. 9.26
            !Nu = (0.825 + ((0.387 * Ra**(1./6.)) / (1 + (0.492/Pr)**(9./16.))**(8./27.)))**2. !generic works quite well for cooling
            
            if(Ra<=1e9) then
            
                !c = 0.59d0
                !n = 0.25d0            
                !from Incropera Eq. 9.27
                !Nu = 0.68 + ((0.670 * Ra**(1./4.)) / (1 + (0.492/Pr)**(9./16.))**(4./9.)) !laminar  
                
                Nu = 0.68 + 0.670 * (Ra*PrFunction)**(0.25d0) !laminar                                
             
            else if (Ra<=1e10) then ! transition
            
                NuLam = 0.68 + 0.670 * (Ra*PrFunction)**(0.25d0)            
                NuTur = 0.15*(Ra*PrFunction)**(1.0d0/3.0d0)
            
                !NuTur = cUserDefined*Ra**nUSerDefined
            
                Nu = linear(1d0*1e9,NuLam,1d0*1e10,NuTur,Ra)
            
            else ! turbulent            
                !c = 0.1d0
                !n = (1.0d0/3.0)
            
                 Nu = 0.15*(Ra*PrFunction)**(1.0d0/3.0d0)                         
             
            endif    
          
        endif        
       
        hOut = Nu * lambda / lChar    ! : Convection coeff. (external)            
    
        !if(Nu<0.75) then
        !    hOut = Nu * lambda / lChar
        !endif 
        
        if(isnan(hOut) ) then                      
            write(MyMessage,'("calculateHeatTransCoefImmersedPlate: FATAL ERROR : hOut NAN. tFilm= = ",f)') tFilm            
            call Messages(-1,Trim(MyMessage),'FATAL', iUnit,iType)
        endif  
               
        
    end subroutine calculateHeatTransCoefImmersedPlate 
    
    subroutine calculateHeatTransCoefImmersedMorgan(tFilm,tFluid,tWall,lChar,hOut,cUserDefined,nUserDefined,Ra,Nu,iUnit,iType)
    
        use spfGlobalConst, only: PI
        use interpolation
        implicit none
        
        integer, intent(in) :: iUnit,iType
        character (len=maxMessageLength) :: MyMessage
        double precision, intent(in)  :: lChar,tFilm,tFluid,tWall
        double precision, intent(out) :: hOut,Ra,Nu
        double precision, intent(in) :: cUserDefined,nUserDefined
        double precision :: cp, lambda, rho, mu, beta, deltaT, Gr, Pr,&
        c,n, PrFunction,NuLam,NuTur
    
        cp    = getCpWater(tFilm)            ! : specific heat
        lambda = getLambdaWater(tFilm)       ! : conductivity
        rho   = getRhoWater(tFilm)           ! : density
        mu    = getNuWater(tFilm)*rho        ! : dynamic viscosity
        
        beta  = getBetaWater(tFilm)          ! : volumetric expansion coeff.
        !beta  = 1./(tFilm+273.15)           ! : Only valid for ideal gas !!
        deltaT = max(abs(tFluid-tWall),0.1)
        
        Ra = getRayleigh  (beta,rho,deltaT,lChar,mu,cp,lambda,iUnit,iType)               
                           
        c = cUserDefined
        n = nUserDefined        
        Nu = c*Ra**n
                           
        hOut = Nu * lambda / lChar    ! : Convection coeff. (external)            
           
        
        if(isnan(hOut) ) then                      
            write(MyMessage,'("calculateHeatTransCoefImmersedMorgan: FATAL ERROR : hOut NAN. tFilm= = ",f)') tFilm            
            call Messages(-1,Trim(MyMessage),'FATAL', iUnit,iType)
        endif  
               
        
    end subroutine calculateHeatTransCoefImmersedMorgan 
    
      subroutine calculateHeatTransCoefImmersedPipe(tFilm,tFluid,tWall,lChar,hOut,cUserDefined,nUserDefined,Ra,Nu,iUnit,iType)
    
        use spfGlobalConst, only: PI
        use interpolation
        implicit none
        
        integer, intent(in) :: iUnit,iType
        character (len=maxMessageLength) :: MyMessage
        double precision, intent(in)  :: lChar,tFilm,tFluid,tWall
        double precision, intent(out) :: hOut,Ra,Nu
        double precision, intent(in) :: cUserDefined,nUserDefined
        double precision :: cp, lambda, rho, mu, beta, deltaT, Gr, Pr,&
        c,n, PrFunction,NuLam,NuTur
    
        cp    = getCpWater(tFilm)            ! : specific heat
        lambda = getLambdaWater(tFilm)       ! : conductivity
        rho   = getRhoWater(tFilm)           ! : density
        mu    = getNuWater(tFilm)*rho        ! : dynamic viscosity
        
        beta  = getBetaWater(tFilm)          ! : volumetric expansion coeff.
        !beta  = 1./(tFilm+273.15)           ! : Only valid for ideal gas !!
        deltaT = max(abs(tFluid-tWall),0.1)
        
        Ra = getRayleigh  (beta,rho,deltaT,lChar,mu,cp,lambda,iUnit,iType)               
       
        if(cUserDefined>0.0 .and. nUserDefined>0.0) then
            
            c = cUserDefined
            n = nUserDefined    
            
            if(c>1.0 .or.n>1.0) then
                if(Ra<=1.0e-2) then
                    c = 0675
                    n=0.058
                else if(Ra<=1.0e2) then
                    c = 1.02
                    n=0.148
                else if(Ra<=1e4) then
                    c = 0.85
                    n=0.188
                else if(Ra<=1e7) then
                    c=0.48
                    n=0.250
                else 
                    c = 0.125
                    n = 0.333
                end if
            end if
            
            Nu = c*Ra**n
        
        else 
            
            Pr  = getPrandtl (mu,cp,lambda)
            
            !Churchill and Chu Int Journal of heat and mass transfer 18, 1049, 1975
            
            Nu = (0.6 + 0.387*Ra**(1.0/6.0) / (1 + (0.559/Pr)**(9.0/16.0) )**(8/27.0) )**2                         
                                       
        endif        
       
        hOut = Nu * lambda / lChar    ! : Convection coeff. (external)            
           
        
        if(isnan(hOut) ) then                      
            write(MyMessage,'("calculateHeatTransCoefImmersedPlate: FATAL ERROR : hOut NAN. tFilm= = ",f)') tFilm            
            call Messages(-1,Trim(MyMessage),'FATAL', iUnit,iType)
        endif  
               
        
      end subroutine calculateHeatTransCoefImmersedPipe

      subroutine calculateHeatTransCoefIceToWater(tFilm,tFluid,tWall,d,hOut,cUserDefined,nUserDefined,Ra,Nu,iUnit,iType)
    
        use spfGlobalConst, only: PI
        use interpolation
        implicit none
        
        integer, intent(in) :: iUnit,iType
        character (len=maxMessageLength) :: MyMessage
        double precision, intent(in)  :: d,tFilm,tFluid,tWall
        double precision, intent(out) :: hOut,Ra,Nu
        double precision, intent(in) :: cUserDefined,nUserDefined
        double precision :: cp, lambda, rho, muKin, muDyn, beta, deltaT, &
        gamma_1,gamma_2,theta,Gr,Pr
        !double precision :: term1,term2,term3
        
        
        cp    = getCpWater(tFilm)            ! : specific heat
        lambda = getLambdaWater(tFilm)       ! : conductivity
        rho   = getRhoWater(tFilm)           ! : density
        muKin = getNuWater(tFilm)        ! : kinematic viscosity
        muDyn = muKin*rho
        beta  = getBetaWater(tFilm)          ! : volumetric expansion coeff.       
        
        gamma_1 = 0.793953e-5
        gamma_2 = -0.655908e-7
        
        deltaT = max(abs(tFluid-0.0),0.1)                
        theta = max(abs((4.0-tFluid)/(0.0-4.0)),0.1)
        
       ! term1 = 9.81*d**3.0/muKin**2
       ! term2 = 2*gamma_1*theta*deltaT**2
        !term3 = 1.0+3*gamma_2*theta*deltaT/2.0/gamma_1
        
        Gr = getGrashof (beta,rho,deltaT,d,muDyn)
        !Gr = (9.81*d**3.0/muKin**2)*(2*gamma_1*theta*deltaT**2)*(1.0+3*gamma_2*theta*deltaT/2.0/gamma_1)              
                               
        !Cheng 1988 
            
        Nu = 0.732*Gr**0.25                 
                                                    
        hOut = Nu * lambda / d    ! : Convection coeff. (external)            
    
        Pr = getPrandtl (muDyn,cp,lambda)
        Ra = Gr * Pr
        
        if(isnan(hOut) ) then                      
            write(MyMessage,'("calculateHeatTransCoefImmersedPlate: FATAL ERROR : hOut NAN. tFilm= = ",f)') tFilm            
            call Messages(-1,Trim(MyMessage),'FATAL', iUnit,iType)
        endif  
               
        
      end subroutine calculateHeatTransCoefIceToWater
      
    subroutine calculateHeatTransCoefImmersedPipe_old(tFilm,tFluid,tWall,diameter,hOut,cUserDefined,nUserDefined,iUnit,iType)
    
        use spfGlobalConst
        implicit none
        
        double precision, intent(in)  :: diameter,tFilm,tFluid,tWall
        double precision, intent(out) :: hOut
        double precision, intent(in) :: cUserDefined,nUserDefined
        integer, intent(in) ::iUnit,iType
        double precision :: cp, lambda, rho, mu, beta, Ra, Nu, deltaT, Gr, Pr, lChar,&
        c,n
    
        cp    = getCpWater(tFilm)            ! : specific heat
        lambda = getLambdaWater(tFilm)       ! : conductivity
        rho   = getRhoWater(tFilm)           ! : density
        mu    = getNuWater(tFilm)*rho        ! : dynamic viscosity
        
        beta  = getBetaWater(tFilm)          ! : volumetric expansion coeff.
        !beta  = 1./(tFilm+273.15)           ! : Only valid for ideal gas !!
        deltaT = abs(tFluid-tWall)       
        
        ! Gr*Pr = Ra
        ! Ra = g * beta_f * abs(tWall - tInf) *  lenght**3 / (nu_f * alpha_f)
          
        ! < where did I get this characteristic lenght??
        lChar= pi*diameter ! IT DEPENS ON THE TANK DIAMETER !!!!!
        
        !Gr = getGrashof (beta,rho,deltaT,lChar,mu)
        !Pr = getPrandtl (mu,cp,lambda) 
        
        Ra = getRayleigh  (beta,rho,deltaT,lChar,mu,cp,lambda,iUnit,iType)           
    
        c = cUserDefined 
        n = nUserDefined               
                
        Nu = c*Ra**n
        
        hOut = Nu * lambda / lChar    ! : Convection coeff. (external)            
    
        !write(20,*) 'Ra= ',Ra,' Gr=',Gr,' Pr=',Pr,' GrPr=',Gr*Pr,' Nu=',Nu,' d=',diameter,' LChar=',lChar,' hOut=',hOut,' deltaT=',deltaT,' cp=',cp,' lambda=',lambda,' Tfilm=',tFilm
        !' cUser=',cUserDefined,' nUserDefined=',nUserDefined
        
    end subroutine calculateHeatTransCoefImmersedPipe_old 
    
    
    subroutine calculateHeatTransferWall(lambdaWall,dIn,dOut,hWall,geometry,iUnit,iType)
              
        integer, intent(in) :: geometry
        double precision, intent(in) ::  lambdaWall,dIn,dOut
        double precision, intent(out) :: hWall
        integer, intent(in) :: iUnit,iType
        character (len=maxMessageLength) :: MyMessage
        
        if(geometry==0)   then     !FLAT PLATE
            ! so dIn is the thick of the plate
           hWall  = lambdaWall / dIn
        else if (geometry==1) then ! CILINDER
           hWall  = 1./(0.5*log(dOut/dIn)*dOut/lambdaWall+1e-30) 
        else                                 
            write(MyMessage,'("calculateHeatTransferWall FATAL: ERROR IN GEOMETRY=",d)') geometry
            call Messages(-1,Trim(MyMessage),'FATAL', iUnit,iType)                                  
           
        endif
    
    end subroutine calculateHeatTransferWall
    
    subroutine NuWongForcedConvection(Re,Pr,d,l,mu,muWall,Nu,thick,crossSectionType)
    
        double precision, intent(in) :: Re, Pr, d,l,mu,muWall,thick
        double precision, intent(out) :: Nu
        integer, intent(in) :: crossSectionType
        double precision :: Gz, C,m,n,K
        integer :: regime
        
        C= 0.0; m = 0.0; n=0.0; K=0.0      

    !>------------------------------------------------------------------------
    !>Circular tube 
    !>------------------------------------------------------------------------
            
        Gz=getGraetz(d,l,Re,Pr);
        
        regime=TURBULENT; 
        
        ! I force turbulent flow
        if (Re<=2000.) regime=LAMINAR;

        select case (crossSectionType) 
            
            case(CIRCULAR) 
            
                select case (regime)
            
                    case (LAMINAR)                                !> Laminar convection 

                        if (Gz>=10.) then
                            C=1.86d0;
                            m=1.d0/3.;
                            n=1.d0/3.;
                            K= (d/(l+1e-30))**(1./3.)*(mu/muWall)**0.14;
                                 
                        else 
                            C=3.66d0;
                            m=0.0d0;
                            n=0.0d0;
                            K=1.0d0;
                                    
                        endif
                                
                    case (TURBULENT)                                !> turbulent convection 

                        if (Pr<0.6) then
                            C=0.023d0;
                            m=0.8d0;
                            n=0.4d0;
                            K=1.d0;
                      
                        else 
                            C=0.027d0;
                            m=0.8d0;
                            n=0.33d0;
                            K=(mu/muWall)**0.14;
                        endif                                  
                
                    case default            
                        !call FoundBadParameter(1,'Fatal','WongForcedConvection unrecognized of flow regim=',regime,' in circular shape') 
                        write(*,*) 'WongForcedConvection unrecognized of flow regim'
                        call exit(0)                
                
                end select !regime
                    
                    
            case (RECTANGULAR)
                
                select case (regime)
            
                    case (LAMINAR)                                !> Laminar convection 
                    
                        m = 0.0
                        n = 0.0
                        k = 1.0
                            
                        if(l/thick<=1.0) then
                            C = 2.98d0;                                                    
                        else if(l/thick<=1.4d0) then
                            C = 3.08d0;
                        else if(l/thick<=2.0d0) then
                            C = 3.39d0;
                        else if(l/thick<=3.0d0) then
                            C = 3.96d0;
                        else if(l/thick<=4.0d0) then
                            C = 4.44d0;
                        else if(l/thick<=8.0d0) then
                            C = 5.95d0;
                        else 
                            C = 7.54d0;
                        endif    
                            
                    case (TURBULENT)                                !> turbulent convection SAME as cilinder                       
                        
                        if (Pr<0.6) then
                            C=0.023d0;
                            m=0.8d0;
                            n=0.4d0;
                            K=1.d0;
                      
                        else 
                            C=0.027d0;
                            m=0.8d0;
                            n=0.33d0;
                            K=(mu/muWall)**0.14;
                        endif                                  
                
                    case default           
                        !call FoundBadParameter(1,'Fatal','WongForcedConvection unrecognized of flow regim=',regime,' in rectangular shape') 
                        write(*,*) 'WongForcedConvection unrecognized of flow regim'
                        call exit(0) 
                    
                    end select !regime
                    
            case default                
                write(*,*) 'FATAL WongForcedConvection unrecognized of crossSectionType=',crossSectionType
                call exit(0)         
                
        end select            
                     
        Nu = C*(Re**m)*(Pr**n)*K;
    
    end subroutine NuWongForcedConvection
            

!---------------------------------------------------------------------------------
!>@brief Prandtl dimensionenless number: Pr=mu*cp/lambda
!>@param  mu : dynamic viscosity 
!>@param  cp : specific heat capacity
!>@param  lambda  : thermal conductivity
!>@return Prandt number
!---------------------------------------------------------------------------------
    double precision function getPrandtl (mu,cp,lambda) 
    
        double precision, intent(in) :: mu,cp,lambda        
        
        getPrandtl = mu*cp/(lambda+1e-30); 

    end function getPrandtl
    
!---------------------------------------------------------------------------------
!>@brief Grashof dimensionenless number:
!> Gr=9.81*beta*(rho**2)*deltaT*(L**3)/(mu**2)
!>@param  beta: thermal expansion coefficient
!>@param  rho: density
!>@param  deltaT: temperature difference
!>@param  L: reference length
!>@param  mu: dynamic viscosity
!>@return Grashof number
!---------------------------------------------------------------------------------

    double precision function getGrashof (beta,rho,deltaT,L,mu) 
    
        double precision, intent(in) :: beta,rho,deltaT,L,mu
        
        !if (deltaT<0.) then
            !write(20,*) 'getGrashof: Gr must be > 0. Check the value of delta_T:',deltaT
            !call exit(0)
        !endif

        getGrashof = 9.81*beta*(rho**2.d0)*deltaT*(L**3.d0)/((mu**2.d0)+1e-30)        

        !write(20,*) 'beta=',beta,'rho=',rho,' deltaT=',deltaT,' L=',L,' mu=',mu
        
    end function getGrashof
    
!---------------------------------------------------------------------------------
!>@brief Reynolds dimensionenless number:
!> Re=rho*V*Dh/mu
!>@param  rho: density           [kg/m^3]
!>@param  V: reference velocity  [m/s]
!>@param  Dh: hydraulic diameter [m]
!>@param  mu: dynamic viscosity  [kg/m.s]
!>@return Reynolds number
!---------------------------------------------------------------------------------

    double precision function getReynolds (rho,V,Dh,mu) 
    
        double precision, intent(in) :: rho,V,Dh,mu        
               
        getReynolds = rho*V*Dh/mu
        
    end function getReynolds
 
!---------------------------------------------------------------------------------
!>@brief Reynolds dimensionenless number:
!> Re=V*Dh/nu
!>@param  V: reference velocity  [m/s]
!>@param  Dh: hydraulic diameter [m]
!>@param  mu: kynematic viscosity  [m^2/s]
!>@return Reynolds number
!---------------------------------------------------------------------------------

    double precision function getReynoldsNu (V,Dh,nu) 
    
        double precision, intent(in) :: V,Dh,nu        
               
        getReynoldsNu = V*Dh/nu
        
    end function getReynoldsNu
    
!---------------------------------------------------------------------------------
!>@brief Rayleigh dimensionenless number:
!> Ra=9.81*beta*(rho**2)*delta_t*(L**3)*cp/(mu*lam) = Gr*Pr
!>@param  beta: thermal expansion coefficient  [1/K]
!>@param  rho: density                         [kg/m^3]
!>@param  deltaT: temperature difference       [K]
!>@param  L: reference length                  [m] 
!>@param  mu: dynamic viscosity                [kg/m.s]      
!>@param  cp: heat capacity                    [J/kg K]
!>@param  lam: thermal conductivity            [W/m K]
!>@return Reynolds number
!---------------------------------------------------------------------------------

    double precision function getRayleigh  (beta,rho,deltaT,L,mu,cp,lam,iUnit,iType) 
    
        double precision, intent(in) :: beta,rho,deltaT,L,mu,cp,lam
        integer, intent(in) :: iUnit,iType      
        character (len=maxMessageLength) :: MyMessage             
         
        if (deltaT<0.) then           
            write(MyMessage,'("getRayleight: Ra must be > 0. Check the value of deltaT:",f)') deltaT
            call Messages(-1,Trim(MyMessage),'FATAL', iUnit,iType)
        endif
        
        getRayleigh = 9.81*beta*(rho**2)*deltaT*(L**3)*cp/(mu*lam)
        
    end function getRayleigh 


!---------------------------------------------------------------------------------
!>@brief Graetz dimensionenless number:
!> Gz=0.25*PI*D*Re*Pr/L
!>@param  D: diameter
!>@param  L: length
!>@param  Re: Reynolds number
!>@param  Pr: Prandtl number    
!>@return Graetz number
!---------------------------------------------------------------------------------

    double precision function getGraetz  (D,L,Re,Pr) 
        use spfGlobalConst, only: PI
    
        double precision, intent(in) :: D,L,Re,Pr        
               
        getGraetz = 0.25*pi*D*Re*Pr/L
        
    end function getGraetz


end module heatTransferCoef

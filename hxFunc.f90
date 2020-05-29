!> --------------------------------------------------------------------------
!> @file: hxFunc.f90
!> @brief: functions to calculate the discretized immersed hx of the storage tank 
!> @author SPF institut fur SolarTechnik 
!> @author D.Carbonell,
!> @date 20/08/2015.
!> @todo 
!> --------------------------------------------------------------------------

module hxFunc
    
    use physProp
    use hxModule
    
    use TrnsysFunctions
     
    implicit none
    
    !Specification. Write here global parameters that can be used inside the module
    
    contains        

!--------------------------------------------------------------------------
!>@brief : set default values of the heat exchanger structure. All equal to 0
!>@param oneImmersedHx : structure of type hxStruct
!--------------------------------------------------------------------------

    subroutine setDefaultHx(oneImmersedHx)

        use hxModule  
        use spfGlobalConst
    
        implicit none
    
        type (hxStruct), intent(inout) :: oneImmersedHX   
    
        oneImmersedHx%dIn  = 0.0d0  !> Inlet diamter of the pipe [m]
        oneImmersedHx%dOut = 0.0d0 !> Outlet diamter of the pipe [m]
        !oneImmersedHx%l    = 0.0d0  !> Total lenght of the pipe [m]
        oneImmersedHx%zIn  = 0.0d0  !> Height of the hx inlet absolute
        oneImmersedHx%zOut = 0.0d0  !> Height of the hx outlet absolute           
        oneImmersedHx%zInDef  = 0.0d0  !> Height of the hx inlet absolute
        oneImmersedHx%zOutDef = 0.0d0  !> Height of the hx outlet absolute           
        oneImmersedHx%glycolConc = 0.0d0  !> Antifreeze percentage [%100]                                
        oneImmersedHx%cpConstant = 0.0d0  !> cp constant from output heat power. To be consisten with other decks
        oneImmersedHx%rhoConstant= 0.0d0  !> rho constant from output heat power. To be consisten with other decks
        oneImmersedHx%volume = 0.0d0      !> volume of the Hx 
        oneImmersedHx%numberOfCv= 0                      
        oneImmersedHx%lambdaWall = 0.0d0 !> the thermal conductivity of the wall hx
        oneImmersedHx%nNuHeat = 0.0d0 ! parameter for Nusselt = C*Ra**n
        oneImmersedHx%nNuCool = 0.0d0 ! parameter for Nusselt = C*Ra**n           
        oneImmersedHx%cNuHeat = 0.0d0 ! parameter for Nusselt = C*Ra**n
        oneImmersedHx%cNuCool = 0.0d0 ! parameter for Nusselt = C*Ra**n   
            
        oneImmersedHx%tFluidIn  = 0.0d0   !> inlet fluid tmeperature [K]
        oneImmersedHx%mDot      = 0.0d0  !> The mass flow rate in [kg/s]           
        oneImmersedHx%tFluidOut = 0.0d0   !> outlet fluid tmeperature [K]
        oneImmersedHx%tFluidAv   = 0.0d0  !> averaged fluid tmeperature [K]
        oneImmersedHx%U      = 0.0d0 !> The global heat transfer coefficient [W/m2k]
        oneImmersedHx%UA     = 0.0d0 !> The global heat transfer coefficient [W/k]
        oneImmersedHx%alphaIn = 0.0d0 !> Inlet fluid heat transfer coefficient [W//m2K]
        oneImmersedHx%alphaOut = 0.0d0 !> Outlet fluid heat transfer coefficient [W//m2K]
        oneImmersedHx%alphaWall = 0.0d0 !> wall heat transfer coefficient [W//m2K] 
        oneImmersedHx%alphaIce = 0.0d0 !> wall heat transfer coefficient [W//m2K] 
        oneImmersedHx%epsilon = 0.0d0 !> efficiency of the heat exchanger epsilon = tIn-tOut/(tIn-tStore)
        oneImmersedHx%lmtd = 0.0d0 !> logaritmic mean temperature                     
        oneImmersedHx%dz    = 0.0d0 !> Height of one hx control volume (zTop-zBot)*zFac
        oneImmersedHx%zTop  = 0.0d0 !> Top height of the hx * zFac 
        oneImmersedHx%zBot  = 0.0d0 !> Bottom height of the hx * zFac                                    
        oneImmersedHx%tStoreAv = 0.0d0     !> weighted-averaged store temperature [K]                                               
        oneImmersedHx%area   = 0.0d0 !> The total external area of the heat exchanger      
     
        oneImmersedHx%goUp = .true.  ! true: the Hx goes from down to up of the store, false: from top to to bottom
        oneImmersedHx%qAcum = 0.0d0 ! Total acumulated heat in the Hx
        oneImmersedHx%qUtil = 0.0d0 ! Total usefull heat qUtil + qAcum = qLoss
        oneImmersedHx%qHxToTnk = 0.0d0 ! Total heat by the Hx. So usefull heat to the storage 
        oneImmersedHx%imbalance = 0.0d0 ! heat balancw
        oneImmersedHx%imbalanceIce = 0.0d0 ! heat balancw

        oneImmersedHx%revertedFlow = .false.                     
        
                       
       
        
        oneImmersedHx%iceThickMelt= 0.0                                   !< melted layer in the heat exchanger [m]
        oneImmersedHx%iceThick = 0.0        !< Ice thickness of all parallel HXs [m]
        oneImmersedHx%iceThickIn = 0.0        !< Ice thickness of all parallel HXs [m]
        oneImmersedHx%iceThickMeltOld =0.0 !< melted layer in the heat exchanger at previous time step [m]   
        oneImmersedHx%iceThickOld = 0.0     !< Ice thickness of all parallel HXs at previous time step [m]                                    
        oneImmersedHx%iceThickInOld = 0.0        !< Ice thickness of all parallel HXs [m]
        oneImmersedHx%qIce = 0.0      !< The latent heat for each time step and heat exchanger                                             
        !oneImmersedHx%iceFormed = 0.0 !< ice formed in one hx including parallel hx's ice = iceFormed+iceOld and iceFormed=ds*nphx          
        oneImmersedHx%iceMass = 0.0     !< the mass of ice in all parallel HXs [kg]
        oneImmersedHx%iceMassOld = 0.0     !< the mass of ice in all parallel HXs [kg]
        oneImmersedHx%iceMassMeltedOutside = 0.0
        oneImmersedHx%addedCapacity = 0.0 ! 4.0e6 ! J/m^3
        
        oneImmersedHx%memoryAllocated = 0  
        oneImmersedHx%geometry = 1
        oneImmersedHx%nParallelHx=1       
        oneImmersedHx%hxMode=0
        oneImmersedHx%reference = 0
        oneImmersedHx%nTubes = 64
        oneImmersedHx%isUsed = 0
        
        oneImmersedHx%dOutIceAvg  = 0
        oneImmersedHx%dInIceAvg   = 0
        oneImmersedHx%dMeltIceAvg = 0
        
        oneImmersedHx%iceThickAvg      = 0
        oneImmersedHx%iceThickInAvg    = 0
        oneImmersedHx%iceThickMeltAvg  = 0
            
        !oneImmersedHx%tWallSurface = 0.0 
        !oneImmersedHx%tWallSurfaceOld = 0.0 
        
    end subroutine setDefaultHx
    
!--------------------------------------------------------------------------
!>@brief : Initialize the temperature of heat exchanger arrays equal to storage temp
!>@param immersedHx : all hx structures of type hxStruct
!>@param iceStore : structure of type iceStore
!--------------------------------------------------------------------------
   
    subroutine initializeTempHx(immersedHx,iceStore)

        use hxModule  
        use iceStoreDef
        use spfGlobalConst
    
        implicit none
    
        type (hxStruct), intent(inout) :: immersedHX(nIHX) 
        type(iceStoreStruct), intent(in), target :: iceStore
        integer :: i,iHx
        
        do iHx=1,iceStore%nHx
            
            if(immersedHx(iHx)%isUsed) then
                call calculateTStoreSeenByHx(immersedHx(iHx),iceStore)
                
                do i=1,immersedHx(iHx)%numberOfCv
                    immersedHx(iHx)%T(i) = immersedHx(iHx)%TStore(i)   
                    immersedHx(iHx)%Tice(i) = iceStore%TFreeze   
                end do
            
                ! The last value
                i = immersedHx(iHx)%numberOfCv            
                immersedHx(iHx)%T(i+1) = immersedHx(iHx)%TStore(i)
            end if
        end do
        
    end subroutine initializeTempHx
    
!--------------------------------------------------------------------------
!>@brief : Dynamic memory allocation for the heat exchanger arrays 
!>@param numberOfHx : number of heat exchangers
!>@param nCv : number of storage Cv
!--------------------------------------------------------------------------

    subroutine setMemoryHx(immersedHx,nCv,numberOfUsedHx)
   
        use hxModule
        use iceStoreConst
    
        implicit none
    
        type(hxStruct), intent(inout) :: immersedHx(nIHX)
        integer, intent(in) :: numberOfUsedHx
        integer, intent(in) :: nCv
        integer :: numberOfCv,error, iHx,i, n
    
        do iHx=1,numberOfUsedHx
        
            if(immersedHx(iHx)%memoryAllocated==0 .and. immersedHx(iHx)%isUsed==1) then
            
                immersedHx(iHx)%memoryAllocated=1
                
                numberOfCv = immersedHx(iHx)%numberOfCv
               
                allocate(immersedHx(iHx)%tStore(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%tStore0(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%UACv(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%qHxCv(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%t0(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%tWall(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%tFilm(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%tWall0(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%tIce(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%tIce0(numberOfCv),stat=error)
                
                !> i+1
                allocate(immersedHx(iHx)%t(numberOfCv+1),stat=error)
                allocate(immersedHx(iHx)%tIte(numberOfCv+1),stat=error)                
                allocate(immersedHx(iHx)%vcly(numberOfCv+1),stat=error)
        
                !> nHxcv,nCv
                allocate(immersedHx(iHx)%factorHxToTank(numberOfCv,nCv),stat=error)
                allocate(immersedHx(iHx)%factorTankToHx(numberOfCv,nCv),stat=error)
                allocate(immersedHx(iHx)%iceMassMeltedOutsideCv(numberOfCv,nCv),stat=error)
                
                !ice related vector
                !allocate(immersedHx(iHx)%qMeltCv(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%qIceCv(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%qAcumIceCv(numberOfCv),stat=error)
               
                allocate(immersedHx(iHx)%iceMassCv(numberOfCv),stat=error) 
                allocate(immersedHx(iHx)%iceMassCvOld(numberOfCv),stat=error) 
                !allocate(immersedHx(iHx)%iceMassCvIte(numberOfCv),stat=error)    
                
                allocate(immersedHx(iHx)%iceMode(numberOfCv),stat=error)  
                allocate(immersedHx(iHx)%indexStore(numberOfCv),stat=error)  
                allocate(immersedHx(iHx)%hxCvBlocked(numberOfCv),stat=error)
                
                if(immersedHx(iHx)%geometry==PLATE) then
               
                    allocate(immersedHx(iHx)%iceThickCv(numberOfCv),stat=error)
                    allocate(immersedHx(iHx)%iceThickInCv(numberOfCv),stat=error)
                    !allocate(immersedHx(iHx)%deltaIceThickCv(numberOfCv),stat=error)                                           
                    allocate(immersedHx(iHx)%iceThickMeltCv(numberOfCv),stat=error)                 
                    allocate(immersedHx(iHx)%iceThickCvOld(numberOfCv),stat=error)
                    allocate(immersedHx(iHx)%iceThickInCvOld(numberOfCv),stat=error)
                    allocate(immersedHx(iHx)%iceThickMeltCvOld(numberOfCv),stat=error)    
                    allocate(immersedHx(iHx)%resetIceThickInCv(numberOfCv),stat=error)
                    allocate(immersedHx(iHx)%resetIceThickInAndMeltCv(numberOfCv),stat=error)
                    
                elseif(immersedHx(iHx)%geometry==COIL) then                     
                                         
                    allocate(immersedHx(iHx)%dOutIce(numberOfCv),stat=error) 
                    allocate(immersedHx(iHx)%dInIce(numberOfCv),stat=error) 
                    allocate(immersedHx(iHx)%dOutMeltIce(numberOfCv),stat=error) 
                    allocate(immersedHx(iHx)%dOutIceOld(numberOfCv),stat=error) 
                    allocate(immersedHx(iHx)%dInIceOld(numberOfCv),stat=error) 
                    allocate(immersedHx(iHx)%dOutMeltIceOld(numberOfCv),stat=error) 
                    allocate(immersedHx(iHx)%dOutIceUsed(numberOfCv),stat=error)      
                    allocate(immersedHx(iHx)%dEqCv(numberOfCv),stat=error)
                    allocate(immersedHx(iHx)%dEqCvOld(numberOfCv),stat=error)   
                    
                    allocate(immersedHx(iHx)%resetDInIceToZero(numberOfCv),stat=error)   
                    
                endif
                
                allocate(immersedHx(iHx)%areaIceUsed(numberOfCv),stat=error) 
                allocate(immersedHx(iHx)%phiOverlap(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%fConst(numberOfCv),stat=error)
                                
                allocate(immersedHx(iHx)%volIceCv(numberOfCv),stat=error)
                allocate(immersedHx(iHx)%volIceCvOld(numberOfCv),stat=error)
                                                              
                do i=1,numberOfCv
                    
                    immersedHx(iHx)%tStore(i) = 0.0d0
                    immersedHx(iHx)%tStore0(i) = 0.0d0
                    immersedHx(iHx)%UACv(i) = 0.0d0
                    immersedHx(iHx)%qHxCv(i) = 0.0d0
                    
                    immersedHx(iHx)%t0(i) = 0.0d0
                    immersedHx(iHx)%tWall(i) = 0.0d0
                    immersedHx(iHx)%tFilm(i) = 0.0d0
                    immersedHx(iHx)%tWall0(i) = 0.0d0
                    
                    immersedHx(iHx)%tIce(i) = 0.0d0
                    immersedHx(iHx)%tIce0(i) = 0.0d0
                    
                    immersedHx(iHx)%t(i) = 0.0d0
                    immersedHx(iHx)%tIte(i) = 0.0d0                   
                    immersedHx(iHx)%vcly(i) = 0.0d0            
                
                    !ice related
                    
                   ! immersedHx(iHx)%qMeltCv(i) = 0.0d0
                    immersedHx(iHx)%qIceCv(i) = 0.0d0
                    immersedHx(iHx)%qAcumIceCv(i) = 0.0d0
                   
                    immersedHx(iHx)%iceMassCv(i) = 0.0d0
                    immersedHx(iHx)%iceMassCvOld(i) = 0.0d0
                    !immersedHx(iHx)%iceMassCvIte(i) = 0.0d0
                   
                     immersedHx(iHx)%iceMode(i) = 0 !integer
                    immersedHx(iHx)%indexStore(i) = 0 !integer
                    immersedHx(iHx)%hxCvBlocked(i) = 0 !integer
                    
                    if(immersedHx(iHx)%geometry==PLATE) then
                   
                        immersedHx(iHx)%iceThickCv(i) = 0.0d0
                        immersedHx(iHx)%iceThickInCv(i) = 0.0d0
                        !immersedHx(iHx)%deltaIceThickCv(i) = 0.0                   
                        immersedHx(iHx)%iceThickMeltCv(i) = 0.0
                        immersedHx(iHx)%iceThickCvOld(i) = 0.0 
                        immersedHx(iHx)%iceThickInCvOld(i) = 0.0d0
                        immersedHx(iHx)%iceThickMeltCvOld(i) = 0.0
                        immersedHx(iHx)%resetIceThickInCv(i) = 0
                        immersedHx(iHx)%resetIceThickInAndMeltCv(i) = 0
                        
                    elseif(immersedHx(iHx)%geometry==COIL) then                                        
                        immersedHx(iHx)%dEqCv(i) = 0.0d0
                        immersedHx(iHx)%dEqCvOld(i) = 0.0d0
                        immersedHx(iHx)%dOutIce(i)        = 0.0
                        immersedHx(iHx)%dOutIceUsed(i)    = 0.0
                        immersedHx(iHx)%dInIce(i)         = 0.0
                        immersedHx(iHx)%dOutMeltIce(i)    = 0.0
                        immersedHx(iHx)%dOutIceOld(i)     = 0.0
                        immersedHx(iHx)%dInIceOld(i)      = 0.0 
                        immersedHx(iHx)%dOutMeltIceOld(i) = 0.0
                        
                        immersedHx(iHx)%resetDInIceToZero(i) = 0
                                        
                    endif
                    
                    immersedHx(iHx)%fConst(i) = 1.0d0
                    
                    if(immersedHx(iHx)%geometry==COIL) then
                        immersedHx(iHx)%phiOverlap(i) = 0.0d0
                    else
                        immersedHx(iHx)%phiOverlap(i) = 1.0d0
                    endif
                    
                    immersedHx(iHx)%volIceCv(i) = 0.0d0
                    immersedHx(iHx)%volIceCvOld(i) = 0.0d0
                    
                    !2D
                    do n=1,nCv
                        immersedHx(iHx)%factorHxToTank(i,n) = 0.0d0
                        immersedHx(iHx)%factorTankToHx(i,n) = 0.0d0
                        immersedHx(iHx)%iceMassMeltedOutsideCv(i,n) = 0.0d0
                    enddo
                    
                    
                    immersedHx(iHx)%areaIceUsed(i)    = 0.0
                    
                enddo
        
                immersedHx(iHx)%t(numberOfCv+1) = 0.0d0
                immersedHx(iHx)%tIte(numberOfCv+1) = 0.0d0
                immersedHx(iHx)%vcly(numberOfCv+1) = 0.0d0 
                
            end if
        enddo       

    end subroutine setMemoryHx

!--------------------------------------------------------------------------
!>@brief : Dynamic memory deallocation  for the heat exchanger arrays 
!>@param numberOfHx : number of heat exchangers
!--------------------------------------------------------------------------

    subroutine setMemoryFreeHx(immersedHx,numberOfUsedHx)

        use hxModule
        use iceStoreConst
    
        implicit none
    
        type(hxStruct), intent(inout) :: immersedHx(nIHX)       
        integer, intent(in) :: numberOfUsedHx
        integer :: error, iHx
    
        do iHx=1,numberOfUsedHx        
            if(immersedHx(iHx)%memoryAllocated==1) then
                
                immersedHx(iHx)%memoryAllocated=0
                
                deallocate(immersedHx(iHx)%tStore,stat=error)
                deallocate(immersedHx(iHx)%tStore0,stat=error)
                deallocate(immersedHx(iHx)%UACv,stat=error)
                deallocate(immersedHx(iHx)%qHxCv,stat=error)
                deallocate(immersedHx(iHx)%t0,stat=error)
                deallocate(immersedHx(iHx)%t,stat=error)
                deallocate(immersedHx(iHx)%tIte,stat=error)
                deallocate(immersedHx(iHx)%tWall,stat=error)
                deallocate(immersedHx(iHx)%tWall0,stat=error)
                deallocate(immersedHx(iHx)%tFilm,stat=error)
                deallocate(immersedHx(iHx)%tIce,stat=error)
                deallocate(immersedHx(iHx)%tIce0,stat=error)
                deallocate(immersedHx(iHx)%vcly,stat=error)
                
                deallocate(immersedHx(iHx)%factorHxToTank,stat=error)
                deallocate(immersedHx(iHx)%factorTankToHx,stat=error)
                deallocate(immersedHx(iHx)%iceMassMeltedOutsideCv,stat=error)  
                
                !======================
                !ice related vectors
                !======================
                
                !deallocate(immersedHx(iHx)%qMeltCv,stat=error)
                deallocate(immersedHx(iHx)%qIceCv,stat=error)
                deallocate(immersedHx(iHx)%qAcumIceCv,stat=error)
               
                deallocate(immersedHx(iHx)%iceMassCv,stat=error) 
                deallocate(immersedHx(iHx)%iceMassCvOld,stat=error) 
                !deallocate(immersedHx(iHx)%iceMassCvIte,stat=error)                
                deallocate(immersedHx(iHx)%iceMode,stat=error)
                deallocate(immersedHx(iHx)%hxCvBlocked,stat=error)
                
                if(immersedHx(iHx)%geometry==PLATE) then
               
                    deallocate(immersedHx(iHx)%iceThickCv,stat=error)
                    deallocate(immersedHx(iHx)%iceThickInCv,stat=error)                               
                    deallocate(immersedHx(iHx)%iceThickMeltCv,stat=error)
                    deallocate(immersedHx(iHx)%iceThickMeltCvOld,stat=error)
                    deallocate(immersedHx(iHx)%iceThickInCvOld,stat=error)
                    deallocate(immersedHx(iHx)%iceThickCvOld,stat=error)
                    deallocate(immersedHx(iHx)%resetIceThickInCv,stat=error)
                    deallocate(immersedHx(iHx)%resetIceThickInAndMeltCv,stat=error)
                    
                elseif(immersedHx(iHx)%geometry==COIL) then
                    
                    deallocate(immersedHx(iHx)%dOutIce,stat=error) 
                    deallocate(immersedHx(iHx)%dOutIceUsed,stat=error) 
                    deallocate(immersedHx(iHx)%dInIce,stat=error) 
                    deallocate(immersedHx(iHx)%dOutMeltIce,stat=error) 
                    deallocate(immersedHx(iHx)%dOutIceOld,stat=error) 
                    deallocate(immersedHx(iHx)%dInIceOld,stat=error) 
                    deallocate(immersedHx(iHx)%dOutMeltIceOld,stat=error) 
                    deallocate(immersedHx(iHx)%dEqCv,stat=error)
                    deallocate(immersedHx(iHx)%dEqCvOld,stat=error)
                    deallocate(immersedHx(iHx)%resetDInIceToZero,stat=error)
                
                endif
                
                deallocate(immersedHx(iHx)%areaIceUsed,stat=error)                 
                deallocate(immersedHx(iHx)%volIceCv,stat=error)
                deallocate(immersedHx(iHx)%volIceCvOld,stat=error)
                              
            end if
            
        enddo       

    end subroutine setMemoryFreeHx

!--------------------------------------------------------------------------
!>@brief : iterated the step by step model.
!>@param : dTime time step in seconds
!>@param : oneImmersedHx : the functions only solves one hx
!>@return : oneImmersedHx with the new calculated T, qhxToTank,qAcum,qUtil
!--------------------------------------------------------------------------

    !subroutine calculateStepByStep(iceStore,oneImmersedHx,dTime)
    subroutine calculateStepByStep(iceStore,immersedHx,iHx,dTime)
    
        use hxModule       
        use iceStoreDef
        use iceStoreConst
        use trnsysConstants
        use spfAlgorithmConst
    
        implicit none
    
        !type (hxStruct), intent(inout), target :: oneImmersedHX     
        type (hxStruct), intent(inout), target :: immersedHX(nIHx)
        type (hxStruct), pointer :: oneImmersedHX 
        
        type (iceStoreStruct), intent(inout), target :: iceStore        
        double precision,intent(in) :: dTime 
        double precision :: errorMax, nMaxIter, error,relax
        integer :: iHx,nHxCv, nIte, status, dummy                          
            
        oneImmersedHX => immersedHX(iHx)
        
        relax = 1.0d0
        status = CONTINUE_ITE
        nIte = 0
        
        errorMax = iceStore%errorMaxStorage                          
        nMaxIter = iceStore%nMaxIteStepByStep                
            
        nHxCv = oneImmersedHx%numberOfCv
            
        oneImmersedHx%tIte(1:nHxCv+1) = oneImmersedHx%t(1:nHxCv+1) !DC last should be updated too
            
        
        if(iceStore%itDtTrnsys==395) then
            dummy=1
        endif
         
            
        do while (status==CONTINUE_ITE .and. nIte<=nMaxIter)                        
                                                
            nIte = nIte+1                        
            
            !call calculateStepByStepOneIte(iceStore,oneImmersedHx,dTime)
            call calculateStepByStepOneIte(iceStore,immersedHx,iHx,dTime)
            
            error = getMaxError(oneImmersedHx%tIte,oneImmersedHx%t,nHxCv)  
            
            if (iceStore%itTrnsys>=20) then !this means iteration problems....
                iceStore%dummy = 0
            endif
            
            if (error<errorMax) then
                status = CONVERGED                                   
            elseif (error>100.0) then
                status = DIVERGED                              
                write(iceStore%MyMessage,'("calculateStepByStep CALCULATE HX DIVERGED")') 
                call Messages(-1,Trim(iceStore%MyMessage),'FATAL', iceStore%iUnit,iceStore%iType)                                           
            elseif (nIte>nMaxIter) then
                status = MAX_NITER_REACHED  
                
                if(iceStore%verboseLevel>=1 .and. nMaxIter > 2) then
                    write(iceStore%MyMessage,*) 'calculateStepByStep MAX NUM ITERATION ERROR-T= ',error                 
                    call Messages(-1,Trim(iceStore%MyMessage),'NOTICE', iceStore%iUnit,iceStore%iType)   
                endif
            else
                status = CONTINUE_ITE                 
            endif          
                
            oneImmersedHx%t(1:nHxCv)    = oneImmersedHx%tIte(1:nHxCv)*(1-relax)+relax*oneImmersedHx%t(1:nHxCv)
            oneImmersedHx%tIte(1:nHxCv+1) = oneImmersedHx%t(1:nHxCv+1)
            
        enddo
        
    end subroutine calculateStepByStep
        
        
!--------------------------------------------------------------------------
!>@brief : solves the heat exchanger using a step by step method from beginning to end.
!>@param : dTime time step in seconds
!>@param : oneImmersedHx : the functions only solves one hx
!>@return : oneImmersedHx with the new calculated T, qhxToTank,qAcum,qUtil
!--------------------------------------------------------------------------

    !subroutine calculateStepByStepOneIte(iceStore,oneImmersedHx,dTime)
    subroutine calculateStepByStepOneIte(iceStore,immersedHx,iHx,dTime)

        use hxModule       
        use iceStoreDef
        use iceStoreConst
        use trnsysConstants
    
        implicit none
    
        !type (hxStruct), intent(inout), target :: oneImmersedHX     
        type (hxStruct), intent(inout), target :: immersedHx(nIHX)     
        type (hxStruct), pointer :: oneImmersedHX 
        integer, intent(in) :: iHx        
        type (iceStoreStruct), intent(inout), target :: iceStore
        
        double precision,intent(in) :: dTime       
    
        double precision :: mDot,cp,rhoCpVOverDt    
        double precision :: ap,aw,bb,acumTerm 
        double precision ::uaStep, tMax, tMin, qUtil,qAcum, qHxToTnk, qHxToTnkCv,&
        tAv,sumUA, tStoreAv, tWall, tCalc, Bterm, tLimit, dThxMin
        double precision :: qUtilCv,alphaIn,alphaWall,alphaOut,alphaIce !AhxCv
        integer :: nIni,nend,nInc,n, nCenter, nCenterCv, nCv   
        logical :: limit = .true.,limitedused=.false.
        character (len=maxMessageLength) :: MyMessage
        double precision :: increasedThermalInertia, tSink, tWall2, dummy, myUA,nParallelHx
        integer :: useHxCapacity = 1
        
        oneImmersedHX => immersedHx(iHx)
        
        nCv = oneImmersedHx%numberOfCv
        
        if(oneImmersedHX%goUp) then
            nIni = 1
            nend = nCv+1
            nInc = 1
            nCenter  = 0
        else
            nIni = nCv+1
            nend = 1
            nInc = -1
            nCenter = -1
        endif    
        
        if(useHxCapacity) then
            !capacity term can cause overloading of ice since even that the Cv is blocked the acum term can produce ice.
            !but is it so high to be relevant or it is not well implemented?
            rhoCpVOverDt = (oneImmersedHx%cpConstant*oneImmersedHx%rhoConstant+oneImmersedHx%addedCapacity)*oneImmersedHx%volume/dTime    
            acumTerm = rhoCpVOverDt/oneImmersedHx%numberOfCv                  
        else
            acumTerm = 0.
        end if
        
        
        ! Lowest delta-T for heat exchange calculation without mass flow
        ! Below this value we assumne there is no heat exchange
    
        dThxMin = 0.01   ! Used to be 1.0.
    
        sumUA = 0.0
        tStoreAv = 0.0
        qUtil = 0.0
        qAcum = 0.0
        qHxToTnk = 0.0
        
        !AhxCv     = oneImmersedHx%area/oneImmersedHx%numberOfCv !external heat exchanger area per Cv
        
        cp = oneImmersedHx%cpConstant
        mDot = oneImmersedHx%mDot      
    
        oneImmersedHx%t(nIni) = oneImmersedHx%tFluidIn    
        oneImmersedHx%alphaIn = 0.0
        oneImmersedHx%alphaWall = 0.0
        oneImmersedHx%alphaOut = 0.0        
        oneImmersedHx%alphaIce = 0.0     
        oneImmersedHx%UaIn = 0.0
        oneImmersedHx%UaWall = 0.0
        oneImmersedHx%UaOut = 0.0        
        oneImmersedHx%UaIce = 0.0  
        
        oneImmersedHx%fConstrainedAvg = 0.0
        oneImmersedHx%phiOverlapAvg = 0.0
        
        nParallelHx = oneImmersedHx%nParallelHx 
        
        if(mDot>zeroMassFlowInKgs .and. mDot<=minMassFlowInKgs) then
            mDot=minMassFlowInKgs
            nParallelHx=nParallelHx*mDot/minMassFlowInKgs
        endif
        
        if(mDot<zeroMassFlowInKgs) then
                        
            do n=nIni,nend-nInc,nInc
            
                nCenterCv = n+nCenter
                
                if(oneImmersedHx%iceMode(nCenterCv)>0.) then
                !if(oneImmersedHx%iceMassCv(nCenterCv)>0.) then
                    tSink = iceStore%Tfreeze
                else
                    tSink = oneImmersedHx%tStore(nCenterCv)
                endif
                
                if(abs(tSink-oneImmersedHx%t0(nCenterCv))>dThxMin .and. acumTerm>1e-20) then                                                                                       
                                     
                    tWall = 0.5*(tSink+oneImmersedHx%t(n)) 
                    !tWall = oneImmersedHx%tWall0(nCenterCv)                    
                    !!TO DO
                    
                    !call calculateUAPhysicalOld(iceStore,oneImmersedHx,oneImmersedHx%t(n),oneImmersedHx%tStore(nCenterCv),tWall,oneImmersedHx%UACv(nCenterCv),oneImmersedHx%iceMode(nCenterCv),alphaIn,alphaWall,alphaOut)                                                   
                    
                    !call calculateUAPhysical(iceStore,oneImmersedHx,alphaIn,alphaWall,alphaOut,alphaIce,n,nCenterCv)
                    call calculateUAPhysical(iceStore,immersedHx,iHx,alphaIn,alphaWall,alphaOut,alphaIce,n,nCenterCv)
                    
                    !if the Cv of strage is full, then this Cv can not build ice if ice is in and t<0.
                    ! of t>0 UA!= 0 because otherwise we will never melt.
                    !if(oneImmersedHx%hxCvBlocked(nCenterCv)==1 .and. oneImmersedHx%iceMode(nCenterCv)==1 .and. oneImmersedHx%t(n)<0.0) then
                        !oneImmersedHx%UACv(nCenterCv)=0.0                             
                    !endif
                    
                    ! Use UA to evaluate the losses                    
                    
                    Bterm = oneImmersedHx%UACv(nCenterCv)/acumTerm                                      
                    tCalc = (Bterm*tSink + oneImmersedHx%t0(nCenterCv))/(1+Bterm)   
                        
                    tLimit = tSink
                    
                    ! We limit in a way that the face temperature is not out of limits!!
                    
                    if(n /= nIni) then
                        tLimit = 0.5*(tSink+oneImmersedHx%t(n))
                    endif    
                    
                    !if qHxToTnkCv >0 means that t0 > tStore, then t can not be lower than tStore
                    
                    if(oneImmersedHx%t0(nCenterCv) > tSink) then                             
                        tCalc = max(tCalc,tLimit)                                                                                   
                    else ! if qHxToTnk <0 tStore > t0, then t can no be higher than tStore                                       
                        tCalc = min(tCalc,tLimit) 
                    endif                
                    
                    qHxToTnkCv = -acumTerm * (tCalc - oneImmersedHx%t0(nCenterCv))                                             
                    
                                   
                else
                    ! Do nothing
                    qHxToTnkCv = 0.0
                    tCalc = oneImmersedHx%t0(nCenterCv)                    
                endif    
                
                qHxToTnk = qHxToTnk + qHxToTnkCv   
                oneImmersedHx%qHxCv(nCenterCv) = -qHxToTnkCv*nParallelHx                               
                                          
                tStoreAv = tStoreAv + tSink/oneImmersedHx%numberOfCv
                                    
                if(n==nIni) then
                    oneImmersedHx%t(n+nInc) = tCalc
                    oneImmersedHx%t(n) = tCalc
                else
                    oneImmersedHx%t(n+nInc) = 2*tCalc - oneImmersedHx%t(n)
                endif                
                             
                oneImmersedHx%tWall(nCenterCv) = 0.5*(tStoreAv+tCalc)                            
            
            enddo  
                
            qUtil = 0.0              
            qAcum =-qHxToTnk               

!> =====================================
!> CALCULATION WHEN THERE IS MASS FLOW
!> =====================================

        else
            do n=nIni,nend-nInc,nInc                       
        
                nCenterCv = n+nCenter                        
               
                !if we use tIn of Cv instead of tAv we do not need to iterate !!!
                !tWall = 0.5*(oneImmersedHx%tStore(nCenterCv)+oneImmersedHx%t(n)) !oneImmersedHx%tWall0(nCenterCv)
                !tWall = oneImmersedHx%tWall0(nCenterCv)
                 
                !call calculateUAPhysical(iceStore,oneImmersedHx,alphaIn,alphaWall,alphaOut,alphaIce,n,nCenterCv)
                call calculateUAPhysical(iceStore,immersedHx,iHx,alphaIn,alphaWall,alphaOut,alphaIce,n,nCenterCv)
                 
                !ap Tp = aw * Tw + b            
                !if the Cv of storage is full, then this Cv can not build ice if ice is in and t<0.
                ! of t>0 UA!= 0 because otherwise we will never melt.                                                              
                
                myUA = oneImmersedHx%UACv(nCenterCv)
                
                !If Cv is blocked we do not allow to produce ice.
                if(oneImmersedHx%hxCvBlocked(nCenterCv)==1 .and. (oneImmersedHx%t(n)<0.0 .or. oneImmersedHx%t0(nCenterCv)<0.0)) then                   
                    myUA = 0.0 
                    acumTerm=0.0
                    oneImmersedHx%UACv(nCenterCv) = myUA                                   
                    
                endif
                
                ap = mDot*cp+0.5*(myUA+acumTerm)                       
                aw = mDot*cp-0.5*(myUA+acumTerm)                                
                
                if(oneImmersedHx%iceMode(nCenterCv)==1) then
                    if(oneImmersedHx%geometry==PLATE) then
                        ! Using actual values can cause convergence problems within TRNSYS
                        if(oneImmersedHx%iceThickMeltCvOld(nCenterCv)<1e-10 .or. oneImmersedHx%iceThickInCvOld(nCenterCv)>0.) then
                            tSink = iceStore%tFreeze
                        else
                            tSink = oneImmersedHx%tStore(nCenterCv)
                        endif
                    else if (oneImmersedHx%geometry==COIL) then
                        if(oneImmersedHx%dOutMeltIceOld(nCenterCv)<oneImmersedHx%dOut .or. oneImmersedHx%dInIceOld(nCenterCv)>oneImmersedHx%dOut) then
                            tSink = iceStore%tFreeze    
                        else
                            tSink = oneImmersedHx%tStore(nCenterCv)
                        endif                                        
                    endif
                else
                    tSink = oneImmersedHx%tStore(nCenterCv)
                endif        
                              
                bb = tSink*myUA+ oneImmersedHx%t0(nCenterCv)*acumTerm
                
                oneImmersedHx%t(n+nInc) = (aw * oneImmersedHx%t(n)+bb)/ap                                                                                            
                                
                !We must use the limit function
           
                if(limit) then ! otherwise we might have unphysical values. Check when, low velocities? when conduction may be relevant?
                
                    tMax = max(oneImmersedHx%t(n),tSink)
                    tMax = max(tMax,oneImmersedHx%t0(nCenterCv))
                    tMin = min(oneImmersedHx%t(n),tSink)
                    tMin = min(tMin,oneImmersedHx%t0(nCenterCv))
            
                    if((oneImmersedHx%t(n+nInc)>tMax .or. oneImmersedHx%t(n+nInc)<tMin)) then
                        limitedused = .true.                           
                        oneImmersedHx%t(n+nInc)= max(oneImmersedHx%t(n+nInc),tMin)
                        oneImmersedHx%t(n+nInc)= min(oneImmersedHx%t(n+nInc),tMax)
                    else
                        limitedused = .false.
                    endif                                                                       
                endif               
            
                tAv  = 0.5*(oneImmersedHx%t(n+nInc)+oneImmersedHx%t(n))
                 
                ! From hx view (out negative) < 0 cooling the hx and therefore heating the storage
                qUtilCv = mDot*cp*(oneImmersedHx%t(n+nInc)-oneImmersedHx%t(n))  
                qUtil   = qUtil +  qUtilCv                                                                                                               
                qAcum   = qAcum + acumTerm*(tAv-oneImmersedHx%t0(nCenterCv))                                                                     
                                
            
                if(limitedused==.true.) then
                    qHxToTnkCv = - (mDot*cp*(oneImmersedHx%t(n+nInc)-oneImmersedHx%t(n)) +acumTerm*(tAv-oneImmersedHx%t0(nCenterCv)))
                else    
                    qHxToTnkCv = myUA*(tAv-tSink)
                endif
            
                qHxToTnk = qHxToTnk + qHxToTnkCv
            
                ! including all parallel HXs
                
                oneImmersedHx%qHxCv(nCenterCv) = -qHxToTnkCv*nParallelHx 
                                
                if(oneImmersedHx%iceMassCv(nCenterCv)>0.) then                    
                    if(iceStore%fixTwallInIceMode==1) then
                        tWall   = (tAv+iceStore%TFreeze)*0.5                        
                    else
                        tWall2   = tAv + (qUtilCv/oneImmersedHx%areaIceUsed(nCenterCv)) * (1.0 / alphaIn + 1.0 /alphaWall)      
                        tWall   = tAv - (qHxToTnkCv/oneImmersedHx%areaIceUsed(nCenterCv)) * (1.0 / alphaIn + 1.0 /alphaWall) 
                    endif
                else
                    tWall2   = tAv + (qUtilCv/oneImmersedHx%areaIceUsed(nCenterCv)) * (1.0 / alphaIn + 1.0 /alphaWall) 
                    tWall   = tAv - (qHxToTnkCv/oneImmersedHx%areaIceUsed(nCenterCv)) * (1.0 / alphaIn + 1.0 /alphaWall) 
                endif
          
               
                if(tWall<oneImmersedHx%t(n) .and.  tWall<oneImmersedHx%tStore(nCenterCv) .and. tWall<iceStore%TFreeze) then                   
                    write(iceStore%MyMessage,*) 'tWall Out of limits Cv=',n,' tWall=',tWall,' tHx=',oneImmersedHx%t(n),' tStore=',oneImmersedHx%tStore(nCenterCv)
                    !call messages(-1,trim(iceStore%MyMessage),'WARNING',iceStore%iUnit,iceStore%iType)
                endif
              
                oneImmersedHx%tWall(nCenterCv) = tWall
                
                sumUA = sumUA + oneImmersedHx%UACv(nCenterCv)    
                tStoreAv = tStoreAv + oneImmersedHx%tStore(nCenterCv)/oneImmersedHx%numberOfCv
                
                oneImmersedHx%alphaIn   = oneImmersedHx%alphaIn + alphaIn/oneImmersedHx%numberOfCv
                oneImmersedHx%alphaWall = oneImmersedHx%alphaWall + alphaWall/oneImmersedHx%numberOfCv
                oneImmersedHx%alphaOut  = oneImmersedHx%alphaOut + alphaOut/oneImmersedHx%numberOfCv
                oneImmersedHx%alphaIce  = oneImmersedHx%alphaIce + alphaIce/oneImmersedHx%numberOfCv
        
                oneImmersedHx%UaIn   = oneImmersedHx%UaIn   + alphaIn*oneImmersedHx%areaIceUsed(nCenterCv)
                oneImmersedHx%UaWall = oneImmersedHx%UaWall + alphaWall*oneImmersedHx%areaIceUsed(nCenterCv)
                oneImmersedHx%UaOut  = oneImmersedHx%UaOut  + alphaOut*oneImmersedHx%areaIceUsed(nCenterCv)
                oneImmersedHx%UaIce  = oneImmersedHx%UaIce  + alphaIce*oneImmersedHx%areaIceUsed(nCenterCv)
        
                !oneImmersedHx%fConstrainedAvg = oneImmersedHx%fConstrainedAvg + iceStore%fConstrained/oneImmersedHx%numberOfCv  
                oneImmersedHx%fConstrainedAvg = oneImmersedHx%fConstrainedAvg + oneImmersedHx%fConst(nCenterCv)/oneImmersedHx%numberOfCv   
                oneImmersedHx%phiOverlapAvg   = oneImmersedHx%phiOverlapAvg   + oneImmersedHx%phiOverlap(nCenterCv)/oneImmersedHx%numberOfCv               
                
            enddo                             
        endif          
                         
        oneImmersedHx%tFluidOut =  oneImmersedHx%t(nend)
        oneImmersedHx%qAcum = qAcum*nParallelHx
        oneImmersedHx%qUtil = mDot*cp*(oneImmersedHx%tFluidOut-oneImmersedHx%tFluidIn)*nParallelHx
        oneImmersedHx%qHxToTnk = qHxToTnk*nParallelHx        
        oneImmersedHx%imbalance = oneImmersedHx%qAcum + oneImmersedHx%qUtil + oneImmersedHx%qHxToTnk
        !< Ua is per Hx, not considering all parallel Hx's !!
        oneImmersedHx%UA = sumUA
        oneImmersedHx%U  = oneImmersedHx%UA/oneImmersedHx%area        
        oneImmersedHx%tStoreAv = tStoreAv                                     
        
        !if(oneImmersedHx%hxCvBlocked(nCenterCv)==1 .and. oneImmersedHx%qHxToTnk>0.) then
        !    dummy=1
        !end if
        
    end subroutine calculateStepByStepOneIte
    
  subroutine heatExchangerStatus(iceStore,oneImmersedHx)    
    
    use hxModule            
    use iceStoreDef
    use iceStoreConst
    
    implicit none
    
    type (hxStruct), intent(inout), target :: oneImmersedHX  
    type (iceStoreStruct), intent(inout), target :: iceStore
    integer :: n
        
    do n=1,oneImmersedHx%numberOfCv 
                    
       ! oneImmersedHx%tFilm(n) = 0.75*oneImmersedHx%tWall(n)+0.25*oneImmersedHx%tStore(n)
        oneImmersedHx%tFilm(n) = 0.75*oneImmersedHx%tWall0(n)+0.25*oneImmersedHx%tStore0(n)
            
        if(oneImmersedHx%iceMassOld>0. .and. (oneImmersedHx%tFilm(n)<0. .or. oneImmersedHx%iceMassCvOld(n)>1e-8) ) then                               
            oneImmersedHx%iceMode(n) = 1
        else if(oneImmersedHx%iceMassOld<=0 .and. (oneImmersedHx%tFilm(n) <= min(iceStore%tSubcool,-0.1) .or. oneImmersedHx%iceMode(n) == 1)) then             
            oneImmersedHx%iceMode(n) = 1            
        else
            oneImmersedHx%iceMode(n) = 0            
        endif        
    enddo
        
    
    !for stability
    do n=2,oneImmersedHx%numberOfCv-1 
        if(oneImmersedHx%iceMode(n) == 1 .and. (oneImmersedHx%iceMode(n-1) == 0 .and. oneImmersedHx%iceMode(n+1)==0)) then
            oneImmersedHx%iceMode(n) = 0
        endif
        
        if(oneImmersedHx%iceMode(n) == 0 .and. (oneImmersedHx%iceMode(n-1) == 1 .and. oneImmersedHx%iceMode(n+1)==1)) then
            oneImmersedHx%iceMode(n) = 1
        endif
  
    enddo

  end subroutine heatExchangerStatus
    
  
  subroutine checkHeatExchangerMass(iceStore,oneImmersedHx,ds,n,myCase)    
    
    use hxModule            
    use iceStoreDef
    use iceStoreConst
    use trnsysConstants
        
    implicit none
    
    type (hxStruct), intent(inout), target :: oneImmersedHX      
    type (iceStoreStruct), intent(inout), target :: iceStore    
    integer, intent(in) :: n,myCase
    double precision, intent (in) :: ds
    double precision :: massIn,massOut,diffMass,addMass,AhxCv,mass,mass0,calcMass,calcMassFromVol
  
    if(iceStore%verboseLevel>=3) then 
                   
        AhxCv = oneImmersedHx%area/oneImmersedHx%numberOfCv
        
        if(oneImmersedHx%geometry==PLATE) then
            massIn    = oneImmersedHx%iceThickInCv(n)*AhxCv*iceStore%rhoIce
            massOut   = oneImmersedHx%iceThickCv(n)*AhxCv*iceStore%rhoIce
            addMass   = ds*oneImmersedHx%nParallelHx*AhxCv*iceStore%rhoIce                                  
        endif
        mass  = oneImmersedHx%iceMassCv(n)
        mass0 = oneImmersedHx%iceMassCvOld(n)
        calcMass = mass0 + addMass
        calcMassFromVol = oneImmersedHx%volIceCv(n)*iceStore%rhoIce
        
        diffMass  = mass-massIn-massOut
                        
        if( abs(diffMass) > 1e-10  ) then
            write(iceStore%MyMessage,*) 'calculateIceHxThickness Cv=',n,'MyCase',myCase,'diffMass= ',diffMass,' iceMass=',oneImmersedHx%iceMassCv(n),' addMassInThisDT=',addMass ,'iceMassHX=',massIn+massOut,'iceMassHX-In=',massIn,'iceMassHX-Out=',massOut                            
            call messages(-1,trim(iceStore%MyMessage),'NOTICE',iceStore%iUnit,iceStore%iType)
        endif
    endif
  
  end subroutine checkHeatExchangerMass
  
  
  subroutine calculateIceHxThicknessPipeCv(iceStore,oneImmersedHx,n)    
    
    use hxModule       
    use iceStoreDef
    use iceStoreConst
    use trnsysConstants
    
    implicit none
    
    type (hxStruct), intent(inout), target :: oneImmersedHX
    type (iceStoreStruct), intent(inout), target :: iceStore        
    integer, intent (in) :: n
    
    double precision :: dMeltEq,qIce,r1,r2,r1Total,r2Total,volIceOutAdded,dVolIce,&
                        volOfIceTotal,volOld1,volOld2,vDiff,volIcedCv,volOfIceFromOld,&
                        volOfIceOld, qLatentCalc, qDiff,qAcumIce,iceThickInSum,newIceIn, dl
    double precision :: diffVolError    
    integer ::  severalLayer,myCase,limiting,notpossible
        
    dl = oneImmersedHx%dL        
     
    oneImmersedHx%resetDinIceToZero(n) = 0
    
    ! It could be that we switch from one condition to another during Trnsys iteration and therefore all ifs that are used to decide the case type should be
    ! linked to old values and not to actual ones. Otherwise convergence problems and even divegence may occur
    
    if(oneImmersedHx%dInIceOld(n) > oneImmersedHx%dOut) then ! several layers  DC-FEB-17 DON'T CHANGE !!!!!         
        severalLayer = 1                                                                               
        if(oneImmersedHx%qHxCv(n)>=0.) then !ICING 
            r1 = 0.5*oneImmersedHx%dInIceOld(n)                            
            myCase = 1                              
        else ! MELTING                       
             ! If there is some ice that has not reach yet dMetlIce, then we reorganize the mass
             ! such that dInIce=0 and dMelt is decreased to add the mass of the ice included in dInIce
           
             !To calculate the existing volume we need ot use old time steps !! DC-MAR-2017
             volOld1 = pi/4.0*(oneImmersedHx%dOutIceOld(n)**2-oneImmersedHx%dOutMeltIceOld(n)**2)*dl                        
             volOld2 = pi/4.0*(oneImmersedHx%dInIceOld(n)**2-oneImmersedHx%dOut**2)*dl             
                    
             ! Calculate the equivalent diameter that conserve mass.            
            
             dMeltEq= sqrt(oneImmersedHx%dOutIceOld(n)**2-4d0*(volOld1+volOld2)/pi/dl)
             
             if(isnan(dMeltEq)) then
                 iceStore%dummy=1
             endif
             oneImmersedHx%resetDinIceToZero(n) = 1            
             !oneImmersedHx%dInIce(n) = oneImmersedHx%dOut changing this can create inestabilities DC-FEB-17                                                                  
             r1 = 0.5*dMeltEq                            
             myCase = 2
                            
        endif                                                                        
    else      ! only one layer
        severalLayer = 0
                        
        if(oneImmersedHx%qHxCv(n)>=0.) then !ICING
            ! dIceIn GROW If melted ice then we grow from dOut
            ! In this condition dMelt should not change so we use the Old value 
                                    
            if(oneImmersedHx%dOutMeltIceOld(n) > oneImmersedHx%dOut) then !it grows from rPipe as DiceIn
                r1 = oneImmersedHx%dOut*0.5 
                myCase = 3
            ! dIceOut GROW If not we grow from dOutIce
            else
                r1 = 0.5*oneImmersedHx%dOutIceOld(n) ! it grows from dOutIce
                myCase = 4
            endif
        else ! Melting one layer
            r1 = 0.5*oneImmersedHx%dOutMeltIceOld(n) ! it melts from dOutMelt
            myCase = 5
        endif
    endif                                                                                                      
                    
    volIcedCv = abs(oneImmersedHx%qHxCv(n))*iceStore%dtsec/(iceStore%rhoIce*iceStore%H_pc*oneImmersedHx%nParallelHx)                                      
    r2 = sqrt(volIcedCv/pi/dl+r1**2)                                       
                    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !VOLUME ADDED-MELTED CALCULATION!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                        
    dVolIce = 0.0
    volOfIceTotal = oneImmersedHx%volIceCvOld(n)
    
    if(myCase==1 .or. myCase==3) then ! ICING  : myCase=1 several layer, it grows from dInIce and its limited by dMeltOut
                                      ! ICING  : myCase=3 one layer it grows from dOut as dInIce and its limited by dMeltOut         
        oneImmersedHx%dInIce(n) = r2*2.0
        oneImmersedHx%dOutMeltIce(n) = oneImmersedHx%dOutMeltIceOld(n)
        
        dVolIce      = pi/4.0*(oneImmersedHx%dInIce(n)**2-oneImmersedHx%dInIceOld(n)**2)*dl
        volOfIceTotal = pi/4.0*(oneImmersedHx%dInIce(n)**2-oneImmersedHx%dOut**2)*dl
        !We use the current dOutIce because it may be reduced by external melting !!!
        volOfIceTotal = volOfIceTotal + pi/4.0*(oneImmersedHx%dOutIce(n)**2-oneImmersedHx%dOutMeltIceOld(n)**2)*dl !DC-FEB-17
        !volOfIceTotal = volOfIceTotal + pi/4.0*(oneImmersedHx%dOutIceOld(n)**2-oneImmersedHx%dOutMeltIceOld(n)**2)*dl
        
        limiting = 0
        ! Limiting. The inner ice reached the melted ice layer. The rest of ice is built on the outer surface
        ! DC-FEB-17 I set it like this because it can change from melting to icing during different iterations and then it can happen that dOutMeltIce< dInIceOld !!! 
        ! This will cause big amount of ice produced, negative water volume and divergence
        if(oneImmersedHx%dInIce(n)>oneImmersedHx%dOutMeltIceOld(n)) then                    
            limiting = 1
            vDiff = pi/4.0*(oneImmersedHx%dInIce(n)**2-oneImmersedHx%dOutMeltIceOld(n)**2)*dl !DC-FEB-17                                                              
            oneImmersedHx%dInIce(n) = oneImmersedHx%dOutMeltIceOld(n)                                                                                                     
            !After ice grows in the inner layer, it grows in the outer layer
            !I move the vDiff ice volume to the outer surface
            oneImmersedHx%dOutIce(n) = sqrt(4.0*vDiff/pi/dl+oneImmersedHx%dOutIceOld(n)**2)
                        
            if(oneImmersedHx%dInIce(n)>oneImmersedHx%dOutIce(n)) then
                iceStore%dummy=1
            endif
                
            volOfIceTotal = pi/4.0*(oneImmersedHx%dOutIce(n)**2-oneImmersedHx%dOut**2)*dl
            
            if(abs(dVolIce-pi/4.0*(oneImmersedHx%dInIce(n)**2-oneImmersedHx%dInIceOld(n)**2)*dl-vDiff)>1e-10) then
                iceStore%dummy=1
            endif
            !volOIce = pi/4.0*(oneImmersedHx%dInIce(n)**2-oneImmersedHx%dInIceOld(n)**2)*dl
            !volOfIce = volOfIce + vDiff
                                                                           
        endif  
            
    else if(myCase==4) then !ICING : one layer, it grows from dOutIce
    
        oneImmersedHx%dOutIce(n) = r2*2.0                                  
        
        ! Total mass created in this time step. 
                            
        dVolIce = pi*(r2**2-r1**2)*dl       
        
        ! total mass in HX (all in one block)
        
        r1Total = oneImmersedHx%dOut*0.5
        r2Total = r2
        volOfIceTotal = pi*(r2Total**2-r1Total**2)*dl  
        
    else if(myCase==2 .or. myCase==5) then ! myCase=2 MELTING : several layer. Melting from dMeltIceOld. We move dInIce to dOutIce and calculate dEq.
                                           ! myCase=5 MELTING : one layer. Melting from dOutMeltIceOld
                                            
        limiting = 0                
        oneImmersedHx%dOutMeltIce(n) = r2*2.0             
        ! It can not be that in one iteration both grow, but it could happen in different iterations so I need to reset them DC-MAR-2017
        oneImmersedHx%dInIce(n) = oneImmersedHx%dInIceOld(n)
        
        if(isnan(oneImmersedHx%dOutMeltIce(n))) then
            iceStore%dummy=1
        endif
                                       
        !Limiting
        if(oneImmersedHx%dOutMeltIce(n)>oneImmersedHx%dOutIce(n)) then
            r2 = oneImmersedHx%dOutIce(n)/2.0
            oneImmersedHx%dOutMeltIce(n) = oneImmersedHx%dOutIce(n)                            
            limiting = 1
        endif                                               
        
         if(oneImmersedHx%dInIce(n)>oneImmersedHx%dOutMeltIce(n) .and. oneImmersedHx%resetDInIceToZero(n)==0) then
                iceStore%dummy=1
         endif
         
        !r1 = oneImmersedHx%dOutMeltIceOld(n)/2.0 DC-FEB-17  ERROR r1 is calculated before, for myCase2=dEq for myCase=5 is dOutMeltOld
        r2 = oneImmersedHx%dOutMeltIce(n)/2.0
                
        ! total mass left in HX
        r2Total = 0.5*oneImmersedHx%dOutIce(n)
        r1Total = r2                                       
                
        ! Total mass melted in this time step. We simplify and linealize for one time step                    
        dVolIce = pi*(r2**2-r1**2)*dl  
        ! total mass in HX
        volOfIceTotal = pi*(r2Total**2-r1Total**2)*dl            
        
    endif
          
    qIce = dVolIce*iceStore%rhoIce*iceStore%H_pc*oneImmersedHx%nParallelHx/iceStore%dtsec    
    oneImmersedHx%Tice(n) = 0.5*(oneImmersedHx%T(n)+iceStore%Tfreeze)                    
    !qAcumIce is = 0. Not implemented yet
    qAcumIce = dVolIce*iceStore%rhoIce*iceStore%cpIce*(oneImmersedHx%Tice(n)-oneImmersedHx%Tice0(n))*oneImmersedHx%nParallelHx/iceStore%dtsec                    
    oneImmersedHx%qAcumIceCv(n) = qAcumIce
                                                            
    if(oneImmersedHx%qHxCv(n)>=0) then 
        !ice qAcumIce <0 
        oneImmersedHx%qIceCv(n) = max(qIce + qAcumIce,0.)
        volOfIceFromOld = oneImmersedHx%volIceCvOld(n) + dVolIce
    else
        !melt qAcumIce>0
        oneImmersedHx%qIceCv(n) = min(-qIce + qAcumIce,0.)
        volOfIceFromOld  = oneImmersedHx%volIceCvOld(n) - dVolIce
    endif   
    
    if(abs(volOfIceFromOld-volOfIceTotal)>1e-8) then
        diffVolError = volOfIceFromOld-volOfIceTotal
        iceStore%dummy=1
    endif        
                    
    oneImmersedHx%volIceCv(n)  = volOfIceTotal ! volOfIceFromOld 
    oneImmersedHx%iceMassCv(n) = volOfIceTotal*iceStore%rhoIce*oneImmersedHx%nParallelHx             
    oneImmersedHx%dEqCv(n)   = sqrt(4.0*oneImmersedHx%volIceCv(n)/pi/dl+oneImmersedHx%dOut**2)                 
    
  end subroutine calculateIceHxThicknessPipeCv
   
  ! There is no control inside this function that avoids built more ice than the maxIceFrac fixed 
  
  subroutine calculateIceHxThicknessPlateCv(iceStore,oneImmersedHx,n)    
    
    use hxModule       
    use iceStoreDef
    use iceStoreConst
    use trnsysConstants
    
    implicit none
    
    type (hxStruct), intent(inout), target :: oneImmersedHX
    type (iceStoreStruct), intent(inout), target :: iceStore        
    integer, intent (in) :: n
    double precision :: AhxCv, ds, newIceIn, hSensible, addMass, meltedThick, oldMass
    integer :: myCase, notpossible, dummy
            
           
    hSensible = 0.0
    AhxCv = oneImmersedHx%area/oneImmersedHx%numberOfCv               
        
    !ds of one Hx (not of all in parallel)
        
    ds = iceStore%dtsec*oneImmersedHx%qHxCv(n) / ((iceStore%H_pc+hSensible)*iceStore%rhoIce*AhxCv*oneImmersedHx%nParallelHx)                      
                    
    ! We only allow to melt the present ice
    if(ds<0.d0) then
        oldMass = oneImmersedHx%iceThickCvOld(n)+oneImmersedHx%iceThickInCvOld(n) !Error ice was never completely melted because iceThickInCv was not considered
        if(abs(ds)>oldMass/oneImmersedHx%nParallelHx) then !complete melting
            ds = -oldMass/oneImmersedHx%nParallelHx !make sure that if this is limited the sensiebl heat is adapted accordinly
        endif            
    endif 
                                                                          
    !qIce for all HX's
    oneImmersedHx%qIceCv(n) = ds*AhxCv*iceStore%rhoIce*oneImmersedHx%nParallelHx*(iceStore%H_pc+hSensible)/iceStore%dtsec                                              
    oneImmersedHx%resetIceThickInCv(n) = 0
    oneImmersedHx%resetIceThickInAndMeltCv(n) = 0
        
    if(ds<0.) then ! melting   
                                               
        !A fluid film forms between ice sheet and HX surface                              
            
        !if(oneImmersedHx%iceThickInCvOld(n)>0.) then ITERATION PROBLEMS ????
        if(oneImmersedHx%iceThickInCv(n)>0.) then
                            
            ! if a 2nd layer exist but we melt, then we add the two ice layers in and out to the hx surface !!
            !oneImmersedHx%iceThickCv(n)   = oneImmersedHx%iceThickCvOld(n) + oneImmersedHx%iceThickInCvOld(n) + ds*oneImmersedHx%nParallelHx !! DANYK
            oneImmersedHx%iceThickCv(n)   = oneImmersedHx%iceThickCvOld(n) + oneImmersedHx%iceThickInCv(n) + ds*oneImmersedHx%nParallelHx         
                
            !oneImmersedHx%iceThickMeltCv(n) = oneImmersedHx%iceThickMeltCvOld(n) - ds*oneImmersedHx%nParallelHx - oneImmersedHx%iceThickInCv(n)
            oneImmersedHx%iceThickMeltCv(n) = oneImmersedHx%iceThickMeltCvOld(n) - ds*oneImmersedHx%nParallelHx
                
            ! we then need to reset the 2nd layer of ice to zero RESET SHOULD BE DONE AT UPDATE TIME STEP
            oneImmersedHx%resetIceThickInCv(n) = 1
                             
            myCase = 1
        else
            oneImmersedHx%iceThickMeltCv(n) = oneImmersedHx%iceThickMeltCvOld(n) - ds*oneImmersedHx%nParallelHx    
            oneImmersedHx%iceThickCv(n)     = oneImmersedHx%iceThickCvOld(n)     + ds*oneImmersedHx%nParallelHx  
            myCase = 2
                            
            if(oneImmersedHx%iceThickCv(n)==0) then !complete melting 
                myCase=21
            endif
        endif                                                
                        
    else if(ds>0) then !icing
                                
        !meltedThick = oneImmersedHx%iceThickMeltCvOld(n)
        meltedThick = oneImmersedHx%iceThickMeltCv(n)                                   
                                       
        if(meltedThick > 0.) then !a second layer of ice growing is formed 
                
            newIceIn = oneImmersedHx%iceThickInCvOld(n) + ds*oneImmersedHx%nParallelHx               
                
            if(newIceIn<=meltedThick) then !all formed ice on internal 2nd layer
                oneImmersedHx%iceThickInCv(n) = newIceIn
                myCase =  4
                oneImmersedHx%iceThickMeltCv(n) = max(oneImmersedHx%iceThickMeltCvOld(n) - ds*oneImmersedHx%nParallelHx,0.)
            else !part formed on internal and part on external
                oneImmersedHx%iceThickInCv(n) = meltedThick !then this will be set to 0, but if I do it here we can have convergence problems.
                oneImmersedHx%iceThickCv(n)   = oneImmersedHx%iceThickCvOld(n)+ (newIceIn-meltedThick)
                !oneImmersedHx%iceThickCv(n)   = oneImmersedHx%iceThickCvOld(n)+ newIceIn
                oneImmersedHx%resetIceThickInAndMeltCv(n) = 1
                myCase =  5                    
            endif                                
                            
        else !only one layer exist
            myCase = 6
            oneImmersedHx%iceThickCv(n) = oneImmersedHx%iceThickCvOld(n) + ds*oneImmersedHx%nParallelHx
        endif
    else !ds = 0
        oneImmersedHx%iceThickCv(n) = 0.
    endif                                            
                    
    addMass = ds*oneImmersedHx%nParallelHx*AhxCv*iceStore%rhoIce
    
    !if(addMass>0. .and. oneImmersedHx%hxCvBlocked(n)==1 ) then !DC-OCT.2018
    !    dummy=1
    !endif
    
    oneImmersedHx%iceMassCv(n)  = max(oneImmersedHx%iceMassCvOld(n) + addMass,0.0d0)              
    oneImmersedHx%volIceCv(n)  = (oneImmersedHx%iceThickCv(n)+oneImmersedHx%iceThickInCv(n))*AhxCv                             
        
    
  end subroutine calculateIceHxThicknessPlateCv
    
!--------------------------------------------------------------------------
!>@brief : calculates heat exchanger ice thickness melt or .
!>@param : iceStore struct 
!>@param : oneImmersedHx : the functions only solves one hx
!>@return : oneImmersedHx with the new calculated ice thickness 
!--------------------------------------------------------------------------

    subroutine calculateIceHxThickness(iceStore,oneImmersedHx)    
    
        use hxModule       
        use iceStoreDef
        use iceStoreConst
        use trnsysConstants
    
        implicit none
    
        type (hxStruct), intent(inout), target :: oneImmersedHX
        type (iceStoreStruct), intent(inout), target :: iceStore        
        integer :: n, dummy
       
        double precision :: qIceSum, qAcumIceSum, iceMassSum, iceThickSum, &
                            iceThickInSum, iceThickMeltSum !, rhoCpVOverDt, acumTerm 
                                    
        !hSensible = 0.0  
        qAcumIceSum     = 0.0
        qIceSum     = 0.0
        iceMassSum  = 0.0 
        iceThickSum = 0.0
        iceThickInSum = 0.0
        iceThickMeltSum = 0.0
        !volOfIceTotal = 0.0d0
        !dl = oneImmersedHx%dL                                                                
                 
        
        do n=1,oneImmersedHx%numberOfCv                                                                                
                           
            !hSensible = oneImmersedHx%cpConstant*iceStore%tSubcool
                        
            if(oneImmersedHx%iceMode(n)==0 .or. oneImmersedHx%qHxCv(n)==0) then                                
                oneImmersedHx%qIceCv(n) = 0.0                                               
            else
                
                if(oneImmersedHx%hxCvBlocked(n) .and. oneImmersedHx%qHxCv(n)>0.) then
                    dummy=1
                end if
                
                if(oneImmersedHx%geometry==PLATE) then
                    
                   call calculateIceHxThicknessPlateCv(iceStore,oneImmersedHx,n)
                    
                else if(oneImmersedHx%geometry==COIL) then
                                                                                     
                  call calculateIceHxThicknessPipeCv(iceStore,oneImmersedHx,n)                                                                     
                
                endif    ! hx type PLATE or COIL                                                                     
             
            endif                    
                
            qAcumIceSum     = qAcumIceSum + oneImmersedHx%qAcumIceCv(n)
            qIceSum         = qIceSum + oneImmersedHx%qIceCv(n)
            iceMassSum      = iceMassSum + oneImmersedHx%iceMassCv(n)
            
            if(oneImmersedHx%geometry==PLATE) then
                iceThickSum      = iceThickSum     + oneImmersedHx%iceThickCv(n)
                iceThickInSum    = iceThickInSum   + oneImmersedHx%iceThickInCv(n)
                iceThickMeltSum  = iceThickMeltSum + oneImmersedHx%iceThickMeltCv(n)                           
            endif                      
             
            
        enddo     ! n=1,nCv                                                                       
        
        oneImmersedHx%qIce     = qIceSum
        oneImmersedHx%qAcumIce = qAcumIceSum
        oneImmersedHx%iceMass  = iceMassSum
                
        if(oneImmersedHx%geometry==PLATE) then   
            ! The following are averaged values over the Cvs including all Hx's
            oneImmersedHx%iceThick      = icethickSum/oneImmersedHx%numberOfCv     ! *oneImmersedHx%nParallelHx
            oneImmersedHx%iceThickIn    = icethickInSum/oneImmersedHx%numberOfCv   ! *oneImmersedHx%nParallelHx
            oneImmersedHx%iceThickMelt  = iceThickMeltSum/oneImmersedHx%numberOfCv ! *oneImmersedHx%nParallelHx  
            !the following a re per hx plate
            oneImmersedHx%iceThickAvg      = oneImmersedHx%iceThick/oneImmersedHx%nParallelHx
            oneImmersedHx%iceThickInAvg    = oneImmersedHx%iceThickIn/oneImmersedHx%nParallelHx
            oneImmersedHx%iceThickMeltAvg  = oneImmersedHx%iceThickMelt/oneImmersedHx%nParallelHx 
        else
            !The following are per pipe without all hx's or pipes
            oneImmersedHx%dOutIceAvg  = sum(oneImmersedHx%dOutIce(1:oneImmersedHx%numberOfCv))/oneImmersedHx%numberOfCv
            oneImmersedHx%dInIceAvg   = sum(oneImmersedHx%dInIce(1:oneImmersedHx%numberOfCv))/oneImmersedHx%numberOfCv
            oneImmersedHx%dMeltIceAvg = sum(oneImmersedHx%dOutMeltIce(1:oneImmersedHx%numberOfCv))/oneImmersedHx%numberOfCv
                       
        end if       
        
end subroutine calculateIceHxThickness
    
!--------------------------------------------------------------------------
!>@brief : calculates the internal heat trasnfer coef of the heat exchanger for one Cv
!>@param : iceStore: ice storage structure
!>@param : oneImmersedHx : hx structure for one hx
!>@param : tBulk: the averaged water temperature
!>@return : hIn: the internal heat transfer coef in W/m2 of external hx area !!
!--------------------------------------------------------------------------

    subroutine calculateInternalHeatTransferCoef(iceStore,oneImmersedHx,tFluid,hIn)

        use iceStoreDef
        use hxModule
        use heatTransferCoef        
        use TrnsysFunctions

        implicit none     
        
        type(iceStoreStruct), intent(inout), target :: iceStore      
        type(hxStruct), intent(inout), target :: oneImmersedHx 
        
        !integer, intent(in) :: i
        double precision, intent(in) :: tFluid
        double precision, intent(out) :: hIn
        double precision :: Ac, heightHx, thickHx,&
        vel,lChar,Lhx,dOut,b,ReIn,NuIn
        
        heightHx = oneImmersedHx%Hhx   ! : HX characteristic height
        thickHx  = oneImmersedHx%Whx   ! : HX characteristic thickness                                                    
        
        if(isnan(tFluid)) then
            write(iceStore%MyMessage,'("calculateInternalHeatTransferCoef TFluid is NAN. This is typically cause by some error in the inputs of the Hx")') tFluid
            call Messages(-1,Trim(iceStore%MyMessage),'FATAL', iceStore%iUnit,iceStore%iType)     
        endif
            
        if(oneImmersedHx%geometry==PLATE) then
            if(iceStore%useCorrugatedPlate==1) then
                b  = 0.5*thickHx-oneImmersedHx%dxWallHx  !corrugated parameter  
                Ac = b * heightHx
            else
                Ac  = (thickHx-2.0*oneImmersedHx%dxWallHx) * heightHx                 ! : Approx. cross section of one HX 
            endif
        else            
            Ac    = pi*oneImmersedHx%dIn**2/4.0
        endif    
                
        vel = oneImmersedHx%mDot / (oneImmersedHx%rhoConstant * Ac)    ! : Mean fluid velocity
        
        !Dh = getHydraulicDiameter(thick_hx,height_hx) ! L_hx must be total
       
        if(oneImmersedHx%geometry==PLATE) then !PLATE HX                     
            
            if(iceStore%useCorrugatedPlate==1) then                    
                lHx = oneImmersedHx%LHx !oneImmersedHx%area/(2.0*heightHx) !area is total 2 faces !!!                
                call calculateHeatTransCoefFlatPlate(tFluid,oneImmersedHx%glycolConc,b,heightHx,lHx,vel,hIn,oneImmersedHx%nEnhanceNu,iceStore%iUnit,iceStore%iType,oneImmersedHx%ReIn,oneImmersedHx%NuIn)                  
            else
                call calculateHeatTransCoefPipeIn(tFluid,oneImmersedHx%glycolConc,thickHx,heightHx,vel,hIn,.true.,oneImmersedHx%nEnhanceNu,iceStore%iUnit,iceStore%iType,oneImmersedHx%ReIn,oneImmersedHx%NuIn)                 
            endif
            
        else if(oneImmersedHx%geometry==COIL) then !CILINDER HX                        
            !call calculateHeatTransCoefPipeInWong(tBulk,oneImmersedHx%glycolConc,oneImmersedHx%dIn,oneImmersedHx%Lhx,oneImmersedHx%mDot,hIn,0.0d0,CIRCULAR,oneImmersedHx%nEnhanceNu)              
            !vel = vel/oneImmersedHx%nTubes !! DC ERROR becasue the mDot/nHy and nHx includes all parallel pipes !!!!
            call calculateHeatTransCoefPipeIn(tFluid,oneImmersedHx%glycolConc,oneImmersedHx%dIn,oneImmersedHx%Lhx,vel,hIn,.false.,oneImmersedHx%nEnhanceNu,iceStore%iUnit,iceStore%iType,oneImmersedHx%ReIn,oneImmersedHx%NuIn)
            
        endif    
        
    end subroutine calculateInternalHeatTransferCoef
    
    subroutine calculateExternalHeatTransferCoef(iceStore,oneImmersedHx,tFilm,tStorage,tWall,hOut,withIce)
    
        use iceStoreDef
        use hxModule
        
        use heatTransferCoef
     
        implicit none     
        
        type(iceStoreStruct), intent(inout), target :: iceStore     
        type(hxStruct), intent(inout), target :: oneImmersedHx
        integer, intent(in) :: withIce
        double precision, intent(out):: hOut
        double precision, intent(in) :: tFilm,tStorage,tWall
        double precision :: dOut,lChar,cNuMelting,nNuMelting
        

        if(oneImmersedHx%geometry==0) then !PLATE HX
            !The expressions are for vertical plates, therefore the characteristic lenght is height not L_hx
            lChar = oneImmersedHx%Lhx  
            !Lhx  = Ahx/(2.0d0*heightHx)    
        else !CILINDER HX
            dOut  = oneImmersedHx%dOut 
            lChar = dOut                               
        endif    
        
        !Be carefull. we need height_hx as L_hx in RECTANGULAR, but then Graetz number is calculated with height and not with L!!!
                          
        if(oneImmersedHx%tFluidIn<=(tStorage-0d0)) then ! cooling
                       
             if(oneImmersedHx%geometry==PLATE) then
                 call calculateHeatTransCoefImmersedPlate(tFilm,tStorage,&
                    tWall,lChar,hOut,oneImmersedHx%cNuCool,oneImmersedHx%nNuCool,oneImmersedHx%Ra,oneImmersedHx%Nu,iceStore%iUnit,iceStore%iType) 
             else
                 
                 if(withIce==1) then !icing 
                                        
                    !this produces too low alpha and a too low icing process when modified GRashof is used
                    !call calculateHeatTransCoefIceToWater(tFilm,tStorage,tWall,lChar,hOut,oneImmersedHx%cNuCool,oneImmersedHx%nNuCool,oneImmersedHx%Ra,oneImmersedHx%Nu,iceStore%iUnit,iceStore%iType)
                   
                    call calculateHeatTransCoefImmersedPipe(tFilm,tStorage,&
                    tWall,lChar,hOut,oneImmersedHx%cNuCool,oneImmersedHx%nNuCool,oneImmersedHx%Ra,oneImmersedHx%Nu,iceStore%iUnit,iceStore%iType) 
                 else
                    call calculateHeatTransCoefImmersedPipe(tFilm,tStorage,&
                        tWall,lChar,hOut,oneImmersedHx%cNuCool,oneImmersedHx%nNuCool,oneImmersedHx%Ra,oneImmersedHx%Nu,iceStore%iUnit,iceStore%iType) 
                endif
             endif
                                              
        else !heating
            
            if(withIce==1) then ! melting
                
                if(oneImmersedHx%geometry==PLATE) then
                
                    !if we set to zero we get the generic function 
                    !cNuMelting = 0.5  ! 0.5
                    !nNuMelting = 0.0  ! 0.28
                 
                    call calculateHeatTransCoefImmersedPlate(tFilm,tStorage,&
                    tWall,lChar,hOut,oneImmersedHx%cNuHeat,oneImmersedHx%nNuHeat,oneImmersedHx%Ra,oneImmersedHx%Nu,iceStore%iUnit,iceStore%iType)
                    
                else
                 !from Koller et al. 2012 IJR
                 cNuMelting = 0.3
                 nNuMelting = 0.208                
                 
                 call calculateHeatTransCoefImmersedMorgan(tFilm,tStorage,&
                 tWall,lChar,hOut,cNuMelting,nNuMelting,oneImmersedHx%Ra,oneImmersedHx%Nu,iceStore%iUnit,iceStore%iType)       
                    
                endif
                                                          
                 
            elseif(withIce==0) then
                
                 if(oneImmersedHx%geometry==PLATE) then
                    call calculateHeatTransCoefImmersedPlate(tFilm,tStorage,&
                    tWall,lChar,hOut,oneImmersedHx%cNuHeat,oneImmersedHx%nNuHeat,oneImmersedHx%Ra,oneImmersedHx%Nu,iceStore%iUnit,iceStore%iType)                
                 else
                     call calculateHeatTransCoefImmersedPipe(tFilm,tStorage,&
                    tWall,lChar,hOut,oneImmersedHx%cNuHeat,oneImmersedHx%nNuHeat,oneImmersedHx%Ra,oneImmersedHx%Nu,iceStore%iUnit,iceStore%iType) 
                     
                endif
            endif
            
        endif
        
        
            
    end subroutine calculateExternalHeatTransferCoef
    
!--------------------------------------------------------------------
!>@brief Sets the values of oneImmersedHx%tStore for each control 
! volume of one heat exchanger to the temperature this control volume 
! "sees" from the storage tank bulk elements 
!--------------------------------------------------------------------

    subroutine preCalculateTankDataToHEx(immersedHx,iceStore)
   
       use hxModule    
       use iceStoreDef
       
       implicit none
    
       type (hxStruct), intent(inout) :: immersedHX(nIHX)
       type(iceStoreStruct), intent(inout), target :: iceStore        
       double precision :: dhTank,x 
       double precision :: zHxTop, zHxBot,hBotCv,hTopCv,mySum,volHx, view
       integer :: iHex,n,jHx,i,j
            
       do jHx=1,iceStore%nHx
           
           if(immersedHx(jHx)%isUsed) then
               
           do iHex=1,immersedHx(jHx)%numberOfCv                                   
                                  
                zHxTop = immersedHx(jHx)%vcly(iHex+1) ! relative [0,zFac]
                zHxBot = immersedHx(jHx)%vcly(iHex)                    
        
		        do n=1,iceStore%nCv
                     
                    ! I use height of the store to obtain absolute height, zDown is relative height [0-1]!!
            
                    hTopCv = iceStore%vcly(n+1)          ! relative [0,zFac]
                    hBotCv = iceStore%vcly(n)            ! relative [0,zFac]                            

                    dhTank= hTopCv-hBotCv                                 
            
                    if (hBotCv>=zHxTop) then ! out of limits
                
                        exit  !out of loop do n=1,nCv
                
                    else if(zHxTop>hBotCv .and. zHxBot<hTopCv) then !inside of limits
                                
                        x = 0.  ! lenght usefull of each cv of the tank                    
                                                                        
                        if(zHxTop<=hTopCv .and. zHxBot>=hBotCv) then ! all hx is inside 1 storage tank cv
                            x = immersedHx(jHx)%dz                                        
                        else if(zHxTop<=hTopCv .and. zHxBot<=hBotCv) then ! Top part of cv of the storage in higher 
                            x = (zHxTop-hBotCv)*1.d0                             
                        else if(zHxTop>=hTopCv .and. zHxBot<=hBotCv) then ! All storage Cv is considered
                            x = (hTopCv-hBotCv)*1.0d0                    
                        else if (zHxTop>=hTopCv .and. zHxBot>=hBotCv) then ! Bottom part of storage tank is lower 
                            x = (hTopCv-zHxBot)*1.0d0                  
                        else                   
                            call messages(-1,'copyTankDataToHEx() heat exchanger cv <= tank cv. Not controlled case','Fatal',iceStore%iUnit,iceStore%iType)                       
                        endif    
                
                        immersedHx(jHx)%factorHxToTank(iHex,n) = x/immersedHx(jHx)%dz
                        
                        immersedHx(jHx)%indexStore(iHex)=n
                        
                    endif                                  
                enddo             
           enddo
           
           endif
           
       end do
        
       ! facorHxtoTank . We assume that usedHx are sharing the view equally.
       
       do n=1,iceStore%nCv
           do jHx=1,iceStore%nHx
                if(immersedHx(jHx)%isUsed) then
               
                    view = 0.0
                    do iHex=1,immersedHx(jHx)%numberOfCv   
                        view = view + immersedHx(jHx)%factorHxToTank(iHex,n)
                    enddo
                
                    do iHex=1,immersedHx(jHx)%numberOfCv   
                        if(view == 0.0) then
                            immersedHx(jHx)%factorTankToHx(iHex,n) = 0.0
                        else
                            immersedHx(jHx)%factorTankToHx(iHex,n) = immersedHx(jHx)%factorHxToTank(iHex,n)/view
                        endif
                    enddo
                endif    
           enddo
       enddo
       
        
        do i=1,iceStore%nCv
            
            volHx = 0.0
            
            do iHex=1,iceStore%nHx
                if(immersedHx(iHex)%isUsed) then
                    do j=1,immersedHx(iHex)%numberOfCv   
                        volHx      = volHx  + immersedHx(iHex)%volume*immersedHx(iHex)%factorHxToTank(j,i)/immersedHx(iHex)%numberOfCv
                    enddo
                endif
            enddo
            
            iceStore%vHxCv(i)   = volHx
            iceStore%vTankMinusHxCv(i) =  iceStore%vTankCv(i) - iceStore%vHxCv(i)
            
        enddo
        
        iceStore%vTankMinusHx = sum(iceStore%vTankMinusHxCv(1:iceStore%nCv))
    
        if(1) then                  
           
            do jHx=1,iceStore%nHx  
                if(immersedHx(jHx)%isUsed) then
                do iHex=1,immersedHx(jHx)%numberOfCv                                                                                                                                            		                                 
                    mySum =  sum(immersedHx(jHx)%factorHxToTank(iHex,1:iceStore%nCv))
                    
                    if(abs(mySum-1.0d0)>1e-10) then
                        write(iceStore%MyMessage,'("ERROR IN FactorToTankHx in Hx",i" for Hx Cv",i" sum=",e)') jHx,iHex,mySum
                        call Messages(-1,Trim(iceStore%MyMessage),'WARNING', iceStore%iUnit,iceStore%iType)
                        !return
                    end if
                    
                end do
                endif
                
            end do
           
       end if
       
    end subroutine preCalculateTankDataToHEx
  
 !--------------------------------------------------------------------
!>@brief calculates the T bulk water that is seen for each control 
! volume of one heat exchanger 
!--------------------------------------------------------------------

    subroutine calculateTStoreSeenByHx(oneImmersedHx,iceStore)
   
        use hxModule        
        use iceStoreDef
        
        implicit none
        
        type (hxStruct), intent(inout) :: oneImmersedHX    
        type(iceStoreStruct), intent(in), target :: iceStore
        double precision :: dhTank,x 
        integer :: zHxTop, zHxBot,hBotCv,hTopCv
        integer ::i,n
        
        ! This loop should be done in the other way. First n=1,nCv and then iHex. It will be faster
        ! To speed up I can precalculate the factor to be multiplied for each Hx and
        ! then only loop for all and multiply for this factor.
        
        do i=1,oneImmersedHx%numberOfCv         
            
            oneImmersedHx%tStore(i) = 0.      
            
            do n=1,iceStore%nCv            
                !factorHxToTank goes from [0-1]
                oneImmersedHx%tStore(i) = oneImmersedHx%tStore(i) + iceStore%T(n)*oneImmersedHx%factorHxToTank(i,n)      
            end do
            
        end do
       
    end subroutine calculateTStoreSeenByHx
 
   
    
    subroutine calculateGeoHx(immersedHx,iceStore)
    
        use iceStoreDef
        
        type(hxStruct), intent(inout) :: immersedHx(nIHX)
        type(iceStoreStruct), intent(inout) :: iceStore                
        integer :: iHx, i
        double precision :: volOfIce,r1,r2
                                   
        do iHx=1,iceStore%nHx                           
            
            if(immersedHx(iHx)%isUsed==1) then
                
                immersedHx(iHx)%dl = immersedHx(iHx)%Lhx/immersedHx(iHx)%numberOfCv
                
                immersedHx(iHx)%dz = (immersedHx(iHx)%zTop-immersedHx(iHx)%zBot)/immersedHx(iHx)%numberOfCv  
            
                immersedHx(iHx)%vcly(1) =  immersedHx(iHx)%zBot      
            
                do i=1,immersedHx(iHx)%numberOfCv
            
                    immersedHx(iHx)%vcly(i+1) = immersedHx(iHx)%vcly(i)+immersedHx(iHx)%dz ! relative [0-zFac]
                                
                enddo    
            
                immersedHx(iHx)%vcly(immersedHx(iHx)%numberOfCv+1) = immersedHx(iHx)%zTop ! I enforce it to avoid rounding problems    
            
            
                if(immersedHx(iHx)%geometry==PLATE) then !flat plate                                
                                              
                    immersedHx(iHx)%volume = immersedHx(iHx)%Hhx*(immersedHx(iHx)%Whx-2.0*immersedHx(iHx)%dxWallHx)*immersedHx(iHx)%LHx
                    immersedHx(iHx)%area   = immersedHx(iHx)%Hhx*immersedHx(iHx)%LHx*2.0 !including two faces                               
                    
                else if(immersedHx(iHx)%geometry==COIL) then !pipes, coils
                                          
                    immersedHx(iHx)%volume = pi*(immersedHx(iHx)%dIn**2/4.0)*immersedHx(iHx)%LHx 
                    immersedHx(iHx)%area   = pi*immersedHx(iHx)%dOut*immersedHx(iHx)%LHx                             
                   
                end if
                
                immersedHx(iHx)%dA = immersedHx(iHx)%area/immersedHx(iHx)%numberOfCv
                
                do i=1,immersedHx(iHx)%numberOfCv
                    
                    immersedHx(iHx)%volIceCv(i) = immersedHx(iHx)%iceMassCv(i)/iceStore%rhoIce/immersedHx(iHx)%nParallelHx
                    
                    if(immersedHx(iHx)%geometry==COIL) then
                                                
                        immersedHx(iHx)%dOutIce(i)  = sqrt(4.0*immersedHx(iHx)%volIceCv(i)/pi/immersedHx(iHx)%dL+immersedHx(iHx)%dOut**2)                      
                        immersedHx(iHx)%dInIce(i)      = immersedHx(iHx)%dOut
                        immersedHx(iHx)%dOutMeltIce(i) = immersedHx(iHx)%dOut
                    else
                        immersedHx(iHx)%iceThickCv(i)    = immersedHx(iHx)%iceMassCv(i)/(immersedHx(iHx)%area/immersedHx(iHx)%numberOfCv * iceStore%rhoIce)                 
                        immersedHx(iHx)%iceThickCvOld(i) = immersedHx(iHx)%iceThickCv(i)
                    endif  
                        
                    immersedHx(iHx)%areaIceUsed(i) = immersedHx(iHx)%area/immersedHx(iHx)%numberOfCv
                enddo
                    
            endif
            
        enddo
        
        
        call preCalculateTankDataToHEx(immersedHx,iceStore)
        
    end subroutine calculateGeoHx
    
    subroutine calculateUAPhysical(iceStore,immersedHx,iHx,alphaIn,alphaWall,alphaOut,alphaIce,iFace,i)

        use hxModule
        use iceStoreDef
        use spfGlobalConst
        use heatTransferCoef
        use interpolation
        
        implicit none
    
        !type (hxStruct), intent(inout) :: hxData
        type(hxStruct), intent(inout), target :: immersedHx(nIHX)
        type (hxStruct), pointer :: hxData ,hxData2
        type (iceStoreStruct), intent(inout) :: iceStore
        integer, intent(in) :: iFace,i,iHx
        
        double precision, intent(out)  :: alphaIn,alphaOut,alphaWall,alphaIce        
        double precision :: U,notused, iceForOneHx,cNuHeat,nNuHeat
        double precision :: tFluidAv,tStoreAv,tWall,tFilm
        double precision :: lambdaWater,dUsed, phiOverlap,iceThickBetweenTwoPlates
        double precision :: magicFactor,lowLimit,upLimit,alphaIceCond,alphaIceConv, dsRange
        
        !double precision,dimension(immersedHx(iHx)%numberOfCv),intent(out) :: iceThickCvSecondHx
        !double precision, pointer :: iceThickCvSecondHx(:)      
        integer          :: iSecondHx,ii,ntubesX,ntubesY
        double precision :: myPhiOverlap,f
        
        hxData =>  immersedHx(iHx)
        
        if(iHx==1 .and. iceStore%nUsedHx==2) then           
            iSecondHx = 2
        else if(iHx==2 .and. iceStore%nUsedHx==2) then           
            iSecondHx = 1
        else            
            iSecondHx = iHx
        endif
        
        hxData2 => immersedHx(iSecondHx)
        !iceThickCvSecondHx => immersedHx(iSecondHx)%iceThickCv
        
        
        !i = iCenter of CV
        !tFluidAv = hxData%t(iFace)   ! changed to Old 14.11.2016 for stability reasons
        tFluidAv = hxData%t0(i)       ! changed to Old 14.11.2016 for stability reasons         
        tStoreAv = hxData%tStore0(i)  ! changed for stability
        tWall    = hxData%tWall0(i)   ! changed to Old 14.11.2016 for stability reasons  ... added 31.08.2016  
        
        !tWall    = 0.5*(hxData%tStore(i)+hxData%t(iFace))   ! erased  31.08.2016  
        !Changed it was only for no ice 31.08.2016
        tFilm = max(1.0d0, min(90.d0, 0.5*(tWall + tStoreAv))) ! be carefull if using this in other cases !!!   
        
        if(isnan(tFluidAv) .or. isnan(tStoreAv) .or. isnan(tWall) ) then
            write(iceStore%myMessage,*) 'calculateUAPhysical Nan temperatures tFluidAv=',tFluidAv,' tStoreAv=',tStoreAv,' tWall=',tWall
            call messages(-1,trim(iceStore%myMessage),'fatal',iceStore%iUnit,iceStore%iType)              
        end if               
        
        notUsed = 0.0
                       
        if(hxData%iceMode(i)==0) then
                                                            
            !tFilm = tWall 
            
            call calculateExternalHeatTransferCoef(iceStore,hxData,tFilm,tStoreAv,tWall,alphaOut,0) !withoutIce           
            
            if(hxData%geometry==COIL) dUsed = hxData%dOut
            
            alphaIce= 1e5
            
            iceStore%fConstrained = 1.0d0
            hxData%fConst(i)     = 1.0
        else
                                                               
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  U MELT                                         !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
            if(hxData%qHxCv(i) < (-1e-8) .or. hxData%t(iFace)>tStoreAv) then ! Sometimes we get q=-1e-12 but melting is not important and makes fConst to jump to 1.
                        
                alphaOut = 1e5 !not used
                
                ! Water is encapsulated with ice so the surrounding temp is the ice temperature                                                 
                ! When we melt water in between the pipe and the ice layer.                            
                
                iceStore%fConstrained = 1.0d0
                lambdaWater = getLambdaWater(tStoreAv)  !DANYK this should be 0.1 oC??  
                        
                if(hxData%geometry==PLATE) then                                       
                    
                    hxData%phiOverlap(i) = 1.0d0                                                                             
               
                    !lambdaWaterEq = max(lambdaWater,lambdaWater*hxData%Nu)                                       
                    !iceThickBetweenTwoPlates = (hxData%iceThickMeltCvOld(i)+hxData2%iceThickMeltCvOld(i))/hxData%nParallelHx
                    !
                    !if(iceThickBetweenTwoPlates >= (iceStore%xBetweenPipes-1e-4)) then                                           
                    !    hxData%phiOverlap(i) = 0.0
                    !    hxData%fConst(i)     = 0.0
                    !else                           
                    !    hxData%phiOverlap(i) = 1.0d0 
                    !    hxData%fConst(i)     = 1.0d0
                    !endif   
                                            
                    if(hxData%iceThickMeltCv(i)>0.) then     !1cm                        
                        
                        iceForOneHX =  hxData%iceThickMeltCv(i)/hxData%nParallelHx
                        
                        !dsRange  = 1e-3             ! 5 mm
                        !lowLimit = dsRange          !  
                        !upLimit  = lowLimit+dsRange ! 
                                                                      
                        !alphaIce =  lambdaWaterEq*hxData%nParallelHx/hxData%iceThickMeltCv(i) !W/m2  
                        !call calculateExternalHeatTransferCoef(iceStore,hxData,tFilm,tStoreAv,tWall,alphaIceUp,0) !it was set useIceModel=0
                       
                        if(0) then
                            call calculateExternalHeatTransferCoef(iceStore,hxData,tFilm,tStoreAv,tWall,notused,1) !it was set                  
                            alphaIceCond =  lambdaWater*max(hxData%Nu,1.0)/iceForOneHx 
                        else
                            alphaIceCond =  lambdaWater/iceForOneHx ! Pure conduction in water
                        endif
                        
                        !convection 
                        cNuHeat = hxData%cNuHeat
                        nNuHeat = hxData%nNuHeat
                        if(cNuHeat<100 .and. nNuHeat<100) then ! if we use very large values we simulate a slurry with inf alpha                            
                            if(iceStore%xBetweenPipes >= 1/7.) then
                                ! for large distances we use convection Nu = C*Ra^n
                                hxData%cNuHeat = 0.5
                                hxData%nNuHeat = 0.28
                            else
                                ! for low distances we use churchill standard equation
                                hxData%cNuHeat = 0.0
                                hxData%nNuHeat = 0.0
                            endif
                        end if
                        
                        call calculateExternalHeatTransferCoef(iceStore,hxData,tFilm,tStoreAv,tWall,alphaIceConv,1)
                        
                        hxData%cNuHeat = cNuHeat
                        hxData%nNuHeat = nNuHeat
                        
                        !alphaIceLow = alphaIceLow*max(hxData%Nu,1.0)
                        
                        !if(iceForOneHx<lowLimit) then
                        !    alphaIce = alphaIceCond 
                        !else if(iceForOneHx>=lowLimit .and. iceForOneHx<=upLimit) then
                        !    alphaIce = linear(lowLimit,alphaIceCond,upLimit,alphaIceConv,iceForOneHX)
                        !else
                        !   alphaIce = alphaIceConv
                        !endif                                                                        
                        
                        !if(alphaIceHigh>alphaIceLow) then
                        !    iceStore%dummy=1
                        !else
                        !    iceStore%dummy=1
                        !endif                        
                        
                        alphaIce = max(alphaIceCond,alphaIceConv)
                        
                        if(iceStore%timeInHours>=5.) then
                             iceStore%dummy=1
                        endif
        
                        !alphaIceCond is HIGHER than alphaIceConv for low thickness!
                        
                    else
                        alphaIce = 1000                                   
                    endif
                    
                else if(hxData%geometry==COIL) then   !cilinder
                    
                    !This values account that not all tubes are exposed to a 4 surrounding scheme.
                    
                    !ntubesX = (2.0*hxData%nTubes-2.)/hxData%nTubes
                    !ntubesY = (2.0*iceStore%nRealHx-2.)/iceStore%nRealHx
                    
                    ntubesX = hxData%nTubes
                    ntubesY = iceStore%nRealHx
                    
                    dUsed = max(hxData%dOutMeltIce(i),hxData%dOut)
                    
                    if(0) then
                        call getOverlappingAngle(iceStore,iceStore%xBetweenPipes,iceStore%  yBetweenPipes,ntubesX,ntubesY,0.5*dUsed,0.5*dUsed,hxData%dOut,hxData%phiOverlap(i))                                                        
                        
                    else                
                        hxData%phiOverlap(i) = 0.
                    endif 
                    
                    if(abs(dUsed-hxData%dOut)<1e-10) then 
                        alphaIce = 1e5                       
                    else                                                                           
                        f = 1-hxData%phiOverlap(i)/2.0/pi
                        
                        if(f<=0) then   
                            !like if there was no ice. Never used, just in case we may need it
                            call calculateExternalHeatTransferCoef(iceStore,hxData,tFilm,tStoreAv,tWall,alphaIce,1)        
                        else
                            call calculateExternalHeatTransferCoef(iceStore,hxData,tFilm,tStoreAv,tWall,notUsed,1)                                                                                                                     
                            myPhiOverlap = hxData%phiOverlap(i)
                                                                                    
                            if(1) then
                                !myPhiOverlap = 0.0                               
                                alphaIce =  (2.0*pi-myPhiOverlap)*lambdaWater*max(hxData%Nu,1.0)/(pi*dUsed*log(dUsed/hxData%dOut)) !W/m 
                                !alphaIce =  (2.0*pi-myPhiOverlap)*lambdaWater*max(hxData%Nu/dUsed,1.0)/(pi*dUsed*log(dUsed/hxData%dOut)) !W/m 
                            else                                
                                alphaIce =  (2.0*pi-myPhiOverlap)*lambdaWater/(pi*dUsed*log(dUsed/hxData%dOut)) !W/m2 
                            endif
                            
                        endif
                    endif 
                    
                    iceStore%fConstrained = (2.0*pi-hxData%phiOverlap(i))/2.0/pi
                    hxData%fConst(i)      = (2.0*pi-hxData%phiOverlap(i))/2.0/pi
                endif    
                                                       
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  U ICE
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            else
                
                tFilm = max(0.1d0,0.5*(iceStore%TFreeze + tStoreAv)) 
                !tFilm = max(0.1d0,0.5*(tWall + tStoreAv))
                 
                if(hxData%geometry==PLATE) then                                             
            
                    if(0) then
                        call calculateExternalHeatTransferCoef(iceStore,hxData,tFilm,tStoreAv,tWall,alphaOut,1)    !withIce=1
                    else
                        alphaOut=1e5 ! when icing this is irrelevant and does not work for flat plates
                    endif
                    ! for flat plate this is the thickness of only one side because we are using 
                    ! the total area for the calculation including both sides. 
                    ! Its like having the plate 2 times bigger but insolated in the back
                    ! Therefore when we print the thickness from one side is the same as this one
                    
                    if(hxData%iceThickInCv(i)>0.) then !icing the free space heat exchanger and iced blocked due to internal melting
                        iceForOneHX =  hxData%iceThickInCv(i)/hxData%nParallelHx 
                        hxData%phiOverlap(i) = 1.0d0
                    else
                        iceForOneHX =  hxData%iceThickCv(i)/hxData%nParallelHx
                        
                        !to avoid convergence problems
                        iceThickBetweenTwoPlates = (hxData%iceThickCvOld(i)+hxData2%iceThickCvOld(i))/hxData%nParallelHx
                    
                        if(iceThickBetweenTwoPlates >= (iceStore%xBetweenPipes-1e-4)) then                    
                            !iceStore%fConstrained=1e-15
                            hxData%phiOverlap(i) = 0.0
                            hxData%fConst(i)     = 0.0
                        else
                            !iceStore%fConstrained=1.0d0       
                            hxData%phiOverlap(i) = 1.0d0 ! =fConstrained
                            hxData%fConst(i)     = 1.0d0
                        endif   
                    
                    endif                                            
                    
                    alphaIce  = iceStore%kIce*hxData%phiOverlap(i)/(iceForOneHx+1e-15)      ! W/m2K                                                                     
                            
                else !cilinder                                                                                                       
                    
                    if(1) then ! for cilinders we need this alpha, for plates we dont. Why is that?
                        call calculateExternalHeatTransferCoef(iceStore,hxData,tFilm,tStoreAv,tWall,alphaOut,1) !withIce
                    else
                        alphaOut=1e5
                    endif
                    
                    if(hxData%dInIce(i) > hxData%dOut) then   ! partial icing with a ice layer around
                        dUsed = hxData%dInIce(i)
                        hxData%phiOverlap(i) = 0.0
                        
                    else ! continous icing
                        dUsed = hxData%dOutIce(i)
                        
                       ! ntubesX = (2.0*hxData%nTubes-2.)/hxData%nTubes
                        !ntubesY = (2.0*iceStore%nRealHx-2.)/iceStore%nRealHx                             
                        
                        ntubesX = hxData%nTubes
                        ntubesY = iceStore%nRealHx
                        
                         if(iceStore%timeInHours>20) then
                            iceStore%dummy=1
                         endif
                         
                        !hxData%phiOverlap(i) = 0. 
                        call getOverlappingAngle(iceStore,iceStore%xBetweenPipes,iceStore%yBetweenPipes,ntubesX,ntubesY,0.5*hxData%dOutIce(i),0.5*hxData2%dOutIce(i),hxData%dOut,hxData%phiOverlap(i))

                    endif                                                                                              
                           
                    if(abs(dUsed-hxData%dOut)<1e-10) then
                       alphaIce = 1e5
                    else                                  
                       alphaIce =  (2.0*pi-hxData%phiOverlap(i))*iceStore%kIce/(pi*dUsed*log(dUsed/hxData%dOut)) !W/m2          
                    endif                                                       
                    
                    iceStore%fConstrained = (2.0*pi-hxData%phiOverlap(i))/2.0/pi                                       
                    hxData%fConst(i)      = (2.0*pi-hxData%phiOverlap(i))/2.0/pi 
                    
                            
                endif                                                                      
            endif                                      
        endif                                         
        
         !> calculation of RWall
        
        if(hxData%geometry==PLATE) then ! flat plate
            
            !call calculateHeatTransferWall(hxData%lambdaWall,hxData%dxWallHx,notUsed,alphaWall,hxData%geometry,iceStore%iUnit,iceStore%iType)                       
            alphaWall  = hxData%lambdaWall / hxData%dxWallHx !W/m2            
                            
        else !pipe
                                 
            !alphaWall =  2.0*hxData%lambdaWall/(dUsed*log(hxData%dOut/hxData%dIn)*iceStore%fConstrained) !W/m2 of ice or of external diameter
            alphaWall =  2.0*hxData%lambdaWall/(dUsed*log(hxData%dOut/hxData%dIn)*hxData%fConst(i)) !W/m2 of ice or of external diameter
            
        end if                                                             
        
        call calculateInternalHeatTransferCoef(iceStore,hxData,tFluidAv,alphaIn)
        
        ! To convert it to /m2 of Aext.  
        if(hxData%geometry==COIL) then 
            alphaIn = alphaIn*hxData%dIn/dUsed/iceStore%fConstrained
        end if
                         
        ! Overall heat transfer coefficient [W/m2]:        
        
        if(alphaIn<=0.0 .or. alphaWall<=0.0 .or. alphaOut<=0.0 .or. alphaIce<=0.0) then
            U = 0.0
        else    
            U = 1.0d0 / ((1.0 / alphaIn) + (1.0/alphaWall)+ (1.0/alphaIce) + (1.0 / alphaOut))          
        endif                
        
        if(iceStore%useInfinitUA .or. hxData%nNuHeat>100. .or. hxData%nNuCool>100.)  then ! parameter for Nusselt = C*Ra**n 
            U = 1e4
        endif
        
        if (iceStore%itTrnsys>=20 .and. i==2) then ! This means iteration problems....
            iceStore%dummy = 0
        endif
                
        if(hxData%geometry==COIL) then 
            hxData%dOutIceUsed(i) = dUsed
            hxData%areaIceUsed(i) = (pi-0.5*hxData%phiOverlap(i))*dUsed*hxData%dL            
            hxData%UACv(i)        = U*hxData%areaIceUsed(i)
        else
            hxData%areaIceUsed(i) = hxData%area/hxData%numberOfCv 
            hxData%UACv(i)        = U*hxData%areaIceUsed(i)                                  
        endif   
               
        !if(iceStore%timeInHours>=35 .and. hxData%fConst(i)>0.9) then
        !    iceStore%dummy=1
        !endif
        
    end subroutine calculateUAPhysical
         
    !---------------------------------------------------------------------------     
    !> @brief calculation of the storage tank one dimensional heat conduction
    !> equation including a source term of the heat exchanger
    !> @param x1 : distance between hx
    !> @param y1 : distance between pipes
    !> @param ntubesX : number of tubes in x direction (between hx)
    !> @param ntubesY : number of tubes in y direction (for one hx)
    !> @param rIce1 : ice radious for hx 1
    !> @param rIce2 : ice radious for hx 2
    !> @return phiOverlap
    !---------------------------------------------------------------------------
    
    subroutine getOverlappingAngle(iceStore,x1,y1,ntubesX,ntubesY,rIce1,rIce2,rPipe,phiOverlap)
    
        use spfGlobalConst 
        use iceStoreConst               
        use iceStoreDef        
        
        type (iceStoreStruct), intent(inout) :: iceStore       
        double precision, intent(in)  :: x1,y1,rIce1,rIce2,rPipe        
        integer, intent(in)           :: ntubesX,nTubesY ! Because we us an equivalent value to copensate that not all pipes are surrounded by 4
        double precision, intent(out) :: phiOverlap
        double precision :: f
        double precision :: alpha, beta, gamma, small, long, g
        double precision ::v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,rIce1Used,rIce2Used
        
        integer :: Neto = 0
        
        small = 0.5*min(x1,y1) 
        long  = 0.5*max(x1,y1)
                
        if(Neto==0) then
            !if(abs(x1-y1)<1e-10) then
            !    call getOverlappingAngleNewton(iceStore,x1,y1,ntubesX,ntubesY,rIce1,rIce2,f)
            !else
            !     call getOverlappingAngleMattia(iceStore,x1,y1,ntubesX,ntubesY,rIce1,rIce2,rPipe,f)
            !endif
        
            if(small/long>=(1/6.0)-1e-8) then
                call getOverlappingAngleMattia(iceStore,x1,y1,rIce1,rIce2,rPipe,ntubesX,ntubesY,f)
            else
                !call getOverlappingAngleMattia(iceStore,x1,y1,rIce1,rIce2,rPipe,ntubesX,ntubesY,f)
                call getOverlappingAngleMattiaLowArea(iceStore,x1,y1,rIce1,rIce2,rPipe,ntubesX,ntubesY,f)
            endif
            
            !call getOverlappingAngleMattiaOld(iceStore,x1,y1,rIce1,rIce2,rPipe,ntubesX,ntubesY,f)
            !call getOverlappingAngleJekel(iceStore,x1,y1,rIce1,rIce2,f)
       
            ! Modification to consider that the last two pipes in each distance is only affected by one side and not by two.                    
            f = 0.5*((f*(nTubesX-2.) + (f+1.))/nTubesX + (f*(nTubesY-2.) + (f+1.))/nTubesY)               
            phiOverlap = (1.0-f)*2.0*pi  
        
        else
            if(1) then
                v1=-0.24882159
	            v2=3.21603436
	            v3=-0.36384021
	            v4=1.98830385
	            v5=0.38896584
	            v6=0.53800171
	            v7=2.14297332
	            v8=-3.20396299
	            v9=1.41740863
	            v10=-2.38774657                                    
                                                             
                alpha = small/long 
                beta  = rPipe/long 
                
                gamma = rIce1/small 
                
                g = v1*gamma**2+v2*beta**2+v3*alpha**2+v4*gamma*beta+v5*gamma*alpha+v6*alpha*beta+v7*gamma+v8*beta+v9*alpha+v10
                rIce1Used = max(1.0,g)*rIce1
                
                gamma = rIce2/small
                g = v1*gamma**2+v2*beta**2+v3*alpha**2+v4*gamma*beta+v5*gamma*alpha+v6*alpha*beta+v7*gamma+v8*beta+v9*alpha+v10
                rIce2Used = max(1.0,g)*rIce2
            else
                rIce1Used=rIce1
                rIce2Used=rIce2
            endif
            
            call getOverlappingAngleNeto(iceStore,x1,y1,ntubesX,ntubesY,rIce1Used,rIce2Used,phiOverlap)        
        endif
               
        
    end subroutine getOverlappingAngle
    
    subroutine getOverlappingAngleMattia(iceStore,x1,y1,rIce1,rIce2,rPipe,ntubesX,ntubesY,f)
                 
        use spfGlobalConst 
        use iceStoreConst               
        use iceStoreDef
        
        ! implicit none 
        type (iceStoreStruct), intent(inout) :: iceStore
        double precision, intent(in)  :: x1,y1,rIce1,rIce2,rPipe
        integer, intent(in)           :: ntubesX,nTubesY
        double precision, intent(out) :: f
        double precision :: rIceReal,acosSmall,aCosLarge
        double precision :: small,Ar,a1,a2,a3,a4,a5,long,rIce,b1,b3,b4
        integer :: method , method_Rreal
        double precision :: a,b,c,d,e,f1,g,h,i,j,k
        double precision :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
        double precision :: w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,alpha,beta,gamma,alpha_f,beta_f
                                                       
        !x1 = x between pipes
        !y1 = x between hx                
               
        rIce = 0.5*(rIce1+rIce2)                       
        small = 0.5*min(x1,y1) 
        long  = 0.5*max(x1,y1)

        gamma = rIce/small
        
        method_Rreal = 2
        method       = 2
               
        if(rIce >= small) then                                                                       
           
            alpha = small/long 
            beta  = rPipe/long 
            
            ! This is not working 
            ! alpha = small/min(rIce,long) 
            ! beta  = rPipe/min(rIce,long) 
                    
            if(method_Rreal==1) then
                b1 =1.10201842 
                b3 =1.078907    
                b4 =-1.27731665               
                k = rIce/long 
                g = (b3*k)+b4+b1*k/alpha                                             
            else
                v1=-0.24882159
	            v2=3.21603436
	            v3=-0.36384021
	            v4=1.98830385
	            v5=0.38896584
	            v6=0.53800171
	            v7=2.14297332
	            v8=-3.20396299
	            v9=1.41740863
	            v10=-2.38774657    
                !g = v1*gamma**2+v2*alpha**2+v3*beta**2+v4*gamma*alpha+v5*gamma*beta+v6*alpha*beta+v7*gamma+v8*alpha+v9*beta+v10     
                g = v1*gamma**2+v2*beta**2+v3*alpha**2+v4*gamma*beta+v5*gamma*alpha+v6*alpha*beta+v7*gamma+v8*beta+v9*alpha+v10 
            endif
        
            rIceReal = max(1.0,g)*rIce             
            acosLarge  = acos(min(long/rIceReal,1.))
            acosSmall = acos(min(small/rIceReal,1.))          
            Ar = max(1.-2.0*(acosLarge+acosSmall)/pi,0.)                                                                           
                
            alpha_f = small/min(rIceReal,long)
            beta_f = rPipe/min(rIceReal,long)                                                                      
                
            if(method==1) then
                                 
                a1 = 5.22737543e-01   
                a2 = 9.63108062e-01   
                a3 = 1.65792512e-04  
                a4 =-9.51896537e+02
                a5 = 9.51991794e+02                                                                                                  
                f = Ar**((a1*beta_f**a2+a4*alpha_f**a3)+a5)                                
                
            else if(method==2) then                               	                                            
                    
                w1=0.18142704
	            w2=0.12806988
	            w3=0.59817959
	            w4=0.38678919
	            w5=-1.63787425
	            w6=-0.47262632
	            w7=1.43174888
	            w8=0.79350135
	            w9=-0.80224704
	            w10=0.41533001                               
                    
                f = AR**(w1*AR**2+w2*beta_f**2+w3*alpha_f**2+w4*AR*beta_f+w5*AR*alpha_f+w6*beta_f*alpha_f+w7*AR+w8*beta_f+w9*alpha_f+w10)                                                                                                
            endif
            
            f = min(max(f,0.01),0.99)
        else
            f= 1.0d0
        endif
        
    end subroutine getOverlappingAngleMattia
    
    subroutine getOverlappingAngleMattiaLowArea(iceStore,x1,y1,rIce1,rIce2,rPipe,ntubesX,ntubesY,f)
                 
        use spfGlobalConst 
        use iceStoreConst               
        use iceStoreDef
        
        ! implicit none 
        type (iceStoreStruct), intent(inout) :: iceStore
        double precision, intent(in)  :: x1,y1,rIce1,rIce2,rPipe
        integer, intent(in)           :: ntubesX,nTubesY
        double precision, intent(out) :: f
        double precision :: rIceReal,acosSmall,aCosLarge
        double precision :: small,Ar,a1,a2,a3,a4,a5,long,rIce,b1,b3,b4
        integer :: method , method_Rreal
        double precision :: a,b,c,d,e,f1,g,h,i,j,k
        double precision :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
        double precision :: w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,alpha,beta,gamma,alpha_f,beta_f
                                                       
        !x1 = x between pipes
        !y1 = x between hx                
               
        rIce = 0.5*(rIce1+rIce2)                       
        small = 0.5*min(x1,y1) 
        long  = 0.5*max(x1,y1)

        gamma = rIce/small
        
        method_Rreal = 1       
               
        if(rIce >= small) then                                                                       
           
            alpha = small/long 
            beta  = rPipe/long 
                       
            if(method_Rreal==1) then
                beta = ((((rIce**2-rPipe**2)*3.14)/4./small)**2+small**2)**(1./2)
                alpha = (rIce-(small))/((small*long/3.14*4)**(1./2)-(small))
                if (alpha<0.) then
                    rIceReal = rIce
                else
                    rIceReal = min(alpha*beta+(1-alpha)*rIce,beta)
                endif                
            else
                v1=-0.24882159
	            v2=3.21603436
	            v3=-0.36384021
	            v4=1.98830385
	            v5=0.38896584
	            v6=0.53800171
	            v7=2.14297332
	            v8=-3.20396299
	            v9=1.41740863
	            v10=-2.38774657    
                !g = v1*gamma**2+v2*alpha**2+v3*beta**2+v4*gamma*alpha+v5*gamma*beta+v6*alpha*beta+v7*gamma+v8*alpha+v9*beta+v10     
                g = v1*gamma**2+v2*beta**2+v3*alpha**2+v4*gamma*beta+v5*gamma*alpha+v6*alpha*beta+v7*gamma+v8*beta+v9*alpha+v10 
                rIceReal = max(1.0,g)*rIce 
            endif
                                
            acosLarge  = acos(min(long/rIceReal,1.))
            acosSmall = acos(min(small/rIceReal,1.))          
            Ar = max(1.-2.0*(acosLarge+acosSmall)/pi,0.)                                                                           
                
            alpha_f = small/min(rIceReal,long)
            beta_f = rPipe/min(rIceReal,long)                                                                      
            
            w1 = 0.23940494 
            w2 = -0.04949191  
            w3 = 0.60486184  
            w4 = 0.11707023 
            w5 = -1.4871727  
            w6 = -0.56806277
            w7 = 1.32557062  
            w8 = 0.94878816 
            w9 = -0.81571666  
            w10 = 0.4157662                               
                    
            f = AR**(w1*AR**2+w2*beta_f**2+w3*alpha_f**2+w4*AR*beta_f+w5*AR*alpha_f+w6*beta_f*alpha_f+w7*AR+w8*beta_f+w9*alpha_f+w10)                              
            f = min(max(f,0.01),0.99)
        else
            f= 1.0d0
        endif
        
    end subroutine getOverlappingAngleMattiaLowArea
    
     subroutine getOverlappingAngleMattiaOld(iceStore,x1,y1,rIce1,rIce2,rPipe,ntubesX,ntubesY,f)
                 
        use spfGlobalConst 
        use iceStoreConst               
        use iceStoreDef
        
        ! implicit none 
        type (iceStoreStruct), intent(inout) :: iceStore
        double precision, intent(in)  :: x1,y1,rIce1,rIce2,rPipe
        integer, intent(in)           :: ntubesX,nTubesY
        double precision, intent(out) :: f
        double precision :: rFrac,phiOverlap1,phiOverlap2, rIceReal,acosSmall,aCosLarge
        double precision :: rCrit,Ar,r,a1,a2,a3,a4,a5,long,rIce,myX, b1,b3,b4,WD,beta
        integer :: method = 2, useRIceReal=1
        double precision :: a,b,c,d,e,f1,g,h,i,j,DD
        double precision :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
        double precision :: w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,alpha,gamma
                                                       
        !x1 = x between pipes
        !y1 = x between hx                
               
        rIce = 0.5*(rIce1+rIce2)                       
        rCrit = 0.5*min(x1,y1) 
        long  = 0.5*max(x1,y1)
                
        if(rIce >= rCrit) then                                                           
            
            if(method==1) then
                
                !Ar = 1.0 - 2*(acos(min(long/rIce,1.0)) + acos(min(rCrit/rIce,1.)))/pi
                 if(useRIceReal==1) then
                     
                    b1 =1.10201842 
                    b3 =1.078907    
                    b4 =-1.27731665
                    WD = rCrit/long
                    beta = rIce/long    
                    rIceReal = max(1.0,((b3*beta)+b4+b1*beta/WD))*rIce                     
                    
                    a1 = 5.22737543e-01   
                    a2 = 9.63108062e-01   
                    a3 = 1.65792512e-04  
                    a4 =-9.51896537e+02
                    a5 = 9.51991794e+02
                    
                    acosLarge  = acos(min(long/rIceReal,1.))
                    acosSmall = acos(min(rCrit/rIceReal,1.))          
                    Ar = max(1.-2.0*(acosLarge+acosSmall)/pi,0.)
                                     
                 else
                    rIceReal = rIce                    
                    a1 = 3.38430075e-01   
                    a2 = 9.91253378e-01  
                    a3 = -2.24162364e+00   
                    a4 = 2.17633222e-03 
                    a5 = 1.03929081e-01
                    Ar = max(1.-(rIceReal-rCrit)/(2.*sqrt(long*rCrit/pi)-rCrit),0.)
                 endif                                                      
                                                   
                 if(rIceReal<long) then
                     f = Ar**((a1*(rPipe/rIceReal)**a2+a4*(rCrit/rIceReal)**a3)+a5)
                 else
                     f = Ar**((a1*(rPipe/long)**a2+a4*(rCrit/long)**a3)+a5)
                 endif  
                 
            else if(method==2) then                               
	                            
                b1 =1.10201842 
                b3 =1.078907    
                b4 =-1.27731665
                WD = rCrit/long
                beta = rIce/long    
                rIceReal = max(1.0,((b3*beta)+b4+b1*beta/WD))*rIce
                    
                a=0.18142704
	            b=0.12806988
	            c=0.59817959
	            d=0.38678919
	            e=-1.63787425
	            f1=-0.47262632
	            g=1.43174888
	            h=0.79350135
	            i=-0.80224704
	            j=0.41533001
                
                acosLarge  = acos(min(long/rIceReal,1.))
                acosSmall = acos(min(rCrit/rIceReal,1.))          
                Ar = max(1.-2.0*(acosLarge+acosSmall)/pi,0.)                
                
	            if(rIceReal<rCrit) then
		            f = 1
	            else if(rIceReal<long) then		
		            WD = rCrit/rIceReal
		            DD = rPipe/rIceReal
		            f = AR**(a*AR**2+b*DD**2+c*WD**2+d*AR*DD+e*AR*WD+f1*DD*WD+g*AR+h*DD+i*WD+j)
	            else
		            WD = rCrit/long
		            DD = rPipe/long
                    f = AR**(a*AR**2+b*DD**2+c*WD**2+d*AR*DD+e*AR*WD+f1*DD*WD+g*AR+h*DD+i*WD+j)
                end if       
            else if(method==3) then
                
                gamma = rIce/rCrit
                if(gamma<1) then
                    gamma=1.
                end if
                !alpha = small/min(rIce,long)
                !beta  = rPipe/min(rIce,long)
                alpha = rCrit/long
                beta  = rPipe/long
                
                v1=-0.24882159
	            v2=3.21603436
	            v3=-0.36384021
	            v4=1.98830385
	            v5=0.38896584
	            v6=0.53800171
	            v7=2.14297332
	            v8=-3.20396299
	            v9=1.41740863
	            v10=-2.38774657                   
                
                g = v1*gamma**2+v2*alpha**2+v3*beta**2+v4*gamma*alpha+v5*gamma*beta+v6*alpha*beta+v7*gamma+v8*alpha+v9*beta+v10  
                
                rIceReal = max(1.0,g)*rIce
                
                w1=0.18142704
	            w2=0.12806988
	            w3=0.59817959
	            w4=0.38678919
	            w5=-1.63787425
	            w6=-0.47262632
	            w7=1.43174888
	            w8=0.79350135
	            w9=-0.80224704
	            w10=0.41533001
                               
                acosLarge  = acos(min(long/rIceReal,1.))
                acosSmall = acos(min(rCrit/rIceReal,1.))   
                
                Ar = max(1.-2.0*(acosLarge+acosSmall)/pi,0.)                                
                alpha = rCrit/min(rIceReal,long)
                beta  = rPipe/min(rIceReal,long)
                
                f = AR**(w1*AR**2+w2*alpha**2+w3*beta**2+w4*AR*alpha+w5*AR*beta+w6*alpha*beta+w7*AR+w8*alpha+w9*beta+w10)
            else
                
                !icing/largeDistance
                !starting from highest order
                !G16: [ -147.55566018 686.9555387  -1195.55390495  920.69894783 -263.54153282]
                !S16: [ -9.87908828   23.89746704  -19.5122458     5.38881242   0.6112868 ]
                !G8:  [-21.80156826   63.181515    -66.59318098    29.56300956  -3.65442805]
                !S8:  [  4.1265486   -11.6583505   12.33448801     -6.17840486  1.74148702]
                
                myX = 2.0*rIce/sqrt(x1**2+y1**2)
                
                if(nTubesX==64) then !G-type   
                    
                    if(nTubesY==16) then ! G16 Hx = 0.03 
                        a1=24.65570037  
                        a2=-5461.33333432  
                        a3=12120.59640216  
                        a4=-9005.66398349 
                        a5=2233.48857733
                    else if(nTubesY==8) then ! G8
                        a1=-1974.52258502  
                        a2=4465.53890017 
                        a3=-3752.68619431  
                        a4=1387.32709196  
                        a5=-189.31048724
                    endif                    
                else if(nTubesX==48) then !S-type
                    
                    if(nTubesY==16) then ! S-16 Hx = 0.06
                        a1=-321.24525115  
                        a2=568.71209859 
                        a3=-368.10322572  
                        a4=101.45424488   
                        a5=-9.0525556                       
                    else if(nTubesY==8) then ! S-8hx = 0.12
                        a1=-2.64716725e+01   
                        a2=2.28437917e+01   
                        a3=1.21645091e-02  
                        a4=-5.45798929e+00 
                        a5=1.80026313e+00                             
                    endif 
                else
                    iceStore%dummy=1
                endif
                
                f = a5+a4*myX+a3*myX**2+a2*myX**3+a1*myX**4
                
            endif
            
            f = min(max(f,0.01),0.99)
        else
            f= 1.0d0
        endif
        
    end subroutine getOverlappingAngleMattiaOld
     subroutine getOverlappingAngleNeto(iceStore,x1,y1,ntubesX,ntubesY,rIce1,rIce2,phiOverlap)
                 
        use spfGlobalConst 
        use iceStoreConst               
        use iceStoreDef
        
        ! implicit none 
        type (iceStoreStruct), intent(inout) :: iceStore
        integer, intent(in) :: ntubesX,nTubesY ! Because we us an equivalent value to copensate that not all pipes are surrounded by 4
        double precision, intent(in) :: x1,y1,rIce1,rIce2       
        double precision, intent(out) :: phiOverlap   
        double precision :: rFrac,phiOverlap1,phiOverlap2
        integer :: ellipse = 1
        double precision :: rCrit,areaIce,rIceUsedX,rIceUsedY,factor,rIce,nx,ny
        
        
        phiOverlap1 = 0.0d0       
        phiOverlap2 = 0.0d0           
                    
        nx = (2.0*nTubesX-2.)/nTubesX
        ny = (2.0*nTubesY-2.)/ntubesY 
                      
        if(ellipse==0) then
                
            !> see Neto and Krarti Deterministic model ASHRAE Transactions Vol 103 1997.
                 
            if((rIce1+rIce2)>x1 ) then           
                phiOverlap1 =  2.0d0*acos(min((x1**2+rIce1**2-rIce2**2)/(2*rIce1*x1),1.0d0))*nx
            endif                                                
    
            if(rIce1>y1*0.5) then !I'm assuming that both hx have the same ice thickness      
                phiOverlap2 =  2.0d0*acos(min(y1/(2*rIce1),1.0d0))*ny
            endif
        else                                 
            areaIce   = pi*rIce1**2
            rIceUsedY = rIce1
            rIceUsedX = rIce1
            
            if(x1<=y1) then
                
                areaIce  = pi*rIce1**2
                factor   = 1 ! Why did I use this factor ? 2.5/4.
                    
                if(rIce1 >= x1*factor .and. abs(x1-y1)>1e-4) then
                    rIceUsedX = x1*factor
                    rIceUsedY = areaIce/(pi*rIceUsedX) !I use the area of an ellipse (A=pi*b*a) to correct the second ice radious 
                endif
                                        
                if(rIceUsedX>x1*0.5) then !I'm assuming that both hx have the same ice thickness                               
                    phiOverlap1 =  2.0d0*acos(min(x1/(2*rIceUsedX),1.0d0))*nx
                endif
                    
                if(rIceUsedY>y1*0.5) then !I'm assuming that both hx have the same ice thickness                               
                    phiOverlap2 =  2.0d0*acos(min(y1/(2*rIceUsedY),1.0d0))*ny 
                endif
            else                                
                if(rIce1 >= y1) then 
                    rIceUsedX = areaIce/(pi*y1) !I use the area of an ellipse (A=pi*b*a) to correct the second ice radious 
                    rIceUsedY = y1
                endif
                     
                if(rIceUsedY>y1*0.5) then !I'm assuming that both hx have the same ice thickness                               
                    phiOverlap2 =  2.0d0*acos(min(y1/(2*rIceUsedY),1.0d0))*ny
                endif
                    
                if(rIceUsedX>x1*0.5) then !I'm assuming that both hx have the same ice thickness                               
                    phiOverlap1 =  2.0d0*acos(min(x1/(2*rIceUsedX),1.0d0))*nx 
                endif                                                   
               
            endif
            
        endif
            
        phiOverlap = phiOverlap1 + phiOverlap2                                           
       
        !if((rIce1)>y1*0.5 .and. (rIce1+rIce2)>x1) then             
        !    phiOverlap =phiOverlap*magicFactor
        !endif
            
        !elseif((rIce1+rIce2)/2.0>=x1 .or. rIce1>=y1) then
        !    phiOverlap =phiOverlap*1.2
        !end if
        
        !phiOverlap = min(phiOverlap,1.95*pi)            
        

     end subroutine getOverlappingAngleNeto
     
     subroutine getOverlappingAngleJekel(iceStore,x1,y1,rIce1,rIce2,f)
                 
        use spfGlobalConst 
        use iceStoreConst               
        use iceStoreDef
        
        ! implicit none 
        type (iceStoreStruct), intent(inout) :: iceStore        
        double precision, intent(in) :: x1,y1,rIce1,rIce2       
        double precision, intent(out) :: f             
        double precision :: rCrit,rIce,Ar
        integer :: useAveraged=0
                
        !x1 = x between pipes
        !y1 = x between hx        
                             
        rIce  = 0.5*(rIce1+rIce2)         
        !We can't consider a two domensional. Newton assumed a symmetric geometry.
        if(useAveraged==1) then
            rCrit = 0.5*(x1+y1)/2.0
        else
            rCrit = min(x1,y1)/2.0
        endif
                   
        if(rIce >= rCrit) then                            
            Ar = 1.0 - 4.0*acos(min(rCrit/rIce1,1.0))/pi                        
            f = -1.441*Ar+2.455*Ar**0.5+(3.116*Ar-3.158*Ar**0.5)*rIce1/rCrit                                                                
            f = min(max(f,0.01),0.99)
        else
            f= 1.0d0
        endif                    

      end subroutine getOverlappingAngleJekel
      
end module hxFunc
    
    
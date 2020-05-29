!---------------------------------------------------------------------------  
!> @file iceStoreDefinition.f90
!> @author SPF institut fur SolarTechnik 
!> @author D.Carbonell, W.Lodge
!> @date 10/08/2012
!> @brief The type iceStorageStruct inside the iceStorageDef module defines 
!> all the variables needed for the ice storage tank.
!--------------------------------------------------------------------------      

module iceStoreConst
            
    use spfGlobalConst
    
    integer, parameter ::   nMax=100,& ! maximum number of tank nodes allowed
                            nIHX=4     !  Maximum number of immersed heat exchangers                        
                                                          
    integer, parameter :: ne = 10,& ! number of equations
                          na = 22,& ! number of arguments
                          nr = 20   ! number of Newton Rapson solutions
    
    integer, parameter :: COIL  = 1,&
                          PLATE = 0  
    
    double precision, parameter :: tol = 1e-6 ! DANYK some imbalances                                                       
    
    double precision :: zeroMassFlowInKgs = 0.5/3600d0  ! 0.5 kg/h
    double precision :: minMassFlowInKgs  = 2/3600d0   ! 2   kg/h 
    
end module iceStoreConst
          
module hxModule
            
        implicit none        
    
        ! the discretization of the heat exchanger is always from bottom to top.
        ! it does not matter the direction of the mass flow        
        
        type hxStruct
            
            !---------------------------------------------------------------------------
            !> Input parameters of the model. 15 input parameters
            !---------------------------------------------------------------------------
            
            !> 8 parameters for the 3 hx cases
            
            integer :: geometry,&  !< 0 flat plate, 1 cilinder, 2 capilar mat
                       nParallelHx,& !< Number of parallel hx            
                       memoryAllocated,&
                       hxMode,&!< 0 no ice  1 ice 
                       reference,&!> the code I use to know which one it is 1,2,3,4
                       nTubes,& !> in case of callira mats how many tubes per hx.
                       isUsed 
            
            !> common to all geometries
            
            double precision :: Lhx,&    !< HX characteristic Lenght (before was not an input, but AHX)
                                dl,&     !>L/nCv
                                dA       !>Area/nCv
            !< PIPE 
            
            double precision :: dIn,&   !> Inlet diamter of the pipe [m]
                                dOut  !> Outlet diamter of the pipe [m]                                 
            
            !< FLAT PLATE
            
            double precision :: Hhx,&    !< HX height [m] (before was Lhx)
                                Whx,&    !< HX characteristic thickness [m] (before was dHx)                                
                                dxWallHx    !< HX wall thickness [m]
            
            !< capilar mats (PIPE)+
            
            !double precision:: xBetweenManifolds,& !< distance between manifolds                        
                               !xBetweenPipes       !> distance between pipes along the manifold 
            
            integer :: numberOfPipes !> Number of pipes per capilar mat
            
            !< Generic            
         
            double precision :: zIn,zInDef,&    !> Height of the hx inlet absolute. Def means readed from parameters. zIn may change becasue of mass flow direction
                                zOut,zOutDef,&  !> Height of the hx outlet absolute.Def means readed from parameters. zOut may change becasue of mass flow direction           
                                glycolConc,&  !> Antifreeze percentage [%100]                                
                                cpConstant,&  !> cp constant from output heat power. To be consisten with other decks
                                rhoConstant,&!> rho constant from output heat power. To be consisten with other decks
                                addedCapacity !> capacity added to the water capacity of the Hx to include solid structure and piping 
                                
            
            integer :: numberOfCv,&
                       orderHx !> if we have a value from 1 to 4 it means tha the outlet T goes to the following hx indicated by the index
                                      !> 0 or -1 means that the outlet does not go into any hx inlet
                       !allHxIsBlocked !> if all Cv are blocked for ice production =1, else = 0 
            
            double precision :: lambdaWall,& !> the thermal conductivity of the wall hx
                                nNuHeat,& ! parameter for Nusselt = C*Ra**n
                                cNuHeat,& ! parameter for Nusselt = C*Ra**n                       
                                nNuCool,& ! parameter for Nusselt = C*Ra**n
                                cNuCool,& ! parameter for Nusselt = C*Ra**n   
                                nEnhanceNu ! enhanced Nu for intrernal heat trasnfer coefficient Nu = Nu*enhanceNu
           
            !---------------------------------------------------------------------------
            !> Inputs of the model
            !---------------------------------------------------------------------------
            
            double precision :: tFluidIn,&
                                tFluidInRev,& !> inlet fluid tmeperature and reverted fluid temperature [K]
                                tFluidInIt,&   !> inlet from previous iteration
                                mDot,&          !> The mass flow rate in per hx [kg/s]
                                mDotAllHxInKgh     !> The mass flow rate for all hx [kg/h]
            
            !---------------------------------------------------------------------------
            !> Outputs of the model (some may be not used but they can be)
            !---------------------------------------------------------------------------
            
            double precision :: tFluidOut,&    !> outlet fluid tmeperature [K]
                                tFluidAv,&     !> averaged fluid tmeperature [K]
                                tWallAv,&      !> averaged wall tmeperature. Used for icing conditions [K]
                                U,&      !> The global heat transfer coefficient [W/m2k]
                                UA,&     !> The global heat transfer coefficient [W/k]
                                alphaIn,& !> Inlet fluid heat transfer coefficient [W//m2K]
                                alphaOut,& !> Outlet fluid heat transfer coefficient [W//m2K]
                                alphaIce,& !> Outlet fluid heat transfer coefficient [W//m2K]
                                alphaWall,& !> wall heat transfer coefficient [W//m2K] 
                                UaIn,& !> Inlet fluid heat transfer coefficient [W/K]
                                UaOut,& !> Outlet fluid heat transfer coefficient [W/K]
                                UaIce,& !> Outlet fluid heat transfer coefficient [W/K]
                                UaWall,& !> wall heat transfer coefficient [W/K] 
                                epsilon,& !> efficiency of the heat exchanger epsilon = tIn-tOut/(tIn-tStore)
                                fConstrainedAvg,&
                                phiOverlapAvg,&
                                lmtd !> logaritmic mean temperature
            
            !---------------------------------------------------------------------------
            !< Internal parameters
            !---------------------------------------------------------------------------                                                      
            
            double precision :: tStoreAv,& !> weighted-averaged store temperature [K]                                               
                                area,&     !> The total external area of the heat exchanger                                
                                volume,&   !> volume of the Hx 
                                dz,&    !> Height of the hx (zTop-zBot)*zFac
                                zTop,&  !> Top height of the hx in absolute terms
                                zBot  !> Bottom height of the hx in absolute terms                                 
            
            double precision, allocatable :: tStore(:),&   !> Temperature corresponding Cv of the store  
                                             tStore0(:),&
                                             UACv(:),&            
                                             qHxCv(:),& !> heat output for each Hx Cv
                                             t0(:),&   !> temperatures at previosu time step (1,nCv+1) 
                                             t(:),&   !> temperatures (1,nCv+1) 
                                             tIte(:),&   !> temperatures previius iteration (1,nCv+1) 
                                             tWall(:),&
                                             tWall0(:),&   !> Film temperatures at wall (1,nCV) at previous time step
                                             tFilm(:),&    !> Film temperatures at wall (1,nCV)
                                             vcly(:),& !> relative positions of each CV in the storage size (nCv+1) zPosition(1)=hBot zPosition(nCv+1)=zTop
                                             tIce(:),& !> ice temperature in each layer. When all ice is build I allow to extract sensible heat from ice.
                                             tIce0(:)
            
            double precision, allocatable :: factorHxToTank(:,:),& !> (nHxCv,nCv) the weithed factor [0-1] To pass from Hx to Tank
                                             factorTankToHx(:,:),&   !> (nHxCv,nCv) the weithed factor [0-1] To pass from Tank To Hx 
                                             iceMassMeltedOutsideCv(:,:)
            
            !> needed to calculate the T of the storage that is seen for each Hx cv.
            !> 1 all Hx Cv is inside Cv storage, 0 the hx cv don't see the Cv of the storage. 
            !> The sum for 1:nCv = 1 (the Hx control volume is fully considered for the whole storage)
                    
            logical :: goUp  ! true: the Hx goes from down to up of the store, false: from to to bottom
            logical :: revertedFlow  ! true: the flow is negative so the in out heights are changed.
            !> heat flusxes include all parallel Hx
            double precision :: qAcum,& ! [W] Total acumulated heat in the Hx 
                                        ! (>0 the new T is higher than the previous one T0 ==> storage accumulates sensible heat)
                                qUtil,& ! [W] Total usefull heat in the fluid (>0 (cooling storage) Tout > Tin)
                                qHxToTnk,& ! [W] Total heat by the Hx. So usefull heat to the storage (>0 goes to the storage) -qUtil+qSensibleReleased
                                imbalance,& ! [W] heat balance in the Hx balance = qUtil + qAcum + qHxToTnk
                                imbalanceIce ! [W] heat balance in the Hx with ice balance = qUtil + qAcumHx + qHxToTnk + qIce + qAcumIce 
            
 !> Ice related variables
 
            double precision :: iceThickMelt,&    !< melted layer in the heat exchanger [m]
                                iceThick,&        !< Ice thickness of all parallel HXs [m] 
                                iceThickIn,&      !< 2nd ice thickness of all parallel HXs [m] 
                                iceThickMeltOld,& !< melted layer in the heat exchanger at previous time step [m]   
                                iceThickOld,&     !< Ice thickness of all parallel HXs at previous time step [m]                                    
                                iceThickInOld,&   !< 2nd ice thickness of all parallel HXs at previous time step [m]
                                iceMass,&     !< the mass of ice in all parallel HXs [kg]
                                iceMassOld,&     !< the mass of ice in all parallel HXs [kg]
                                iceMassMeltedOutside,& !< the mass of ice in all parallel HXs melted from bulk water [kg]
                                qIce,&        !< The latent heat for each time step and heat exchanger      !< qIce = ds*H_pc*rho_ice / dt => qIce = H_pc*dMice/dt
                                qAcumIce
            
            integer, allocatable :: resetIceThickInCv(:),resetIceThickInAndMeltCv(:)
            
             
            double precision, allocatable :: qIceCv(:),&  ! formed ice o each cv
                                             !qMeltCv(:),& ! melted ice on hx for each cv   
                                             qAcumIceCv(:),&  ! Acumulated sensible heat from the ice
                                             iceThickCv(:),&   ! ice thickness on each cv. Total thickness on one face in the flat plate
                                                               ! radial thickness for pipes (distance from dOut to dIce)
                                             iceThickInCv(:),&   ! 2nd ice layer thickness on each cv before it reaches the iceThickMelt layer. 
                                             iceThickInCvOld(:),&   ! 2nd ice layer thickness on each cv before it reaches the iceThickMelt layer. 
                                             iceThickCvOld(:),&  ! previous time step
                                             iceMassCv(:),&      ! mass of ice for each cv in kg
                                       !      iceMassCvIte(:),&      ! mass of ice for each cv in kg
                                             iceMassCvOld(:),&      ! mass of ice for each cv in kg
                                             volIceCv(:),&    !> Mass of ice for each hx in each Cv
                                             volIceCvOld(:),& !> Mass of ice for each hx in each Cv 
                                             dEqCv(:),& !> Mass of ice for each hx in each Cv
                                             dEqCvOld(:),& !> Mass of ice for each hx in each Cv                                                             
                                             deltaIceThickCv(:),&   ! ice thickness formed in one time step by all hx
                                             iceThickMeltCv(:),&    ! ice thickness formed in one time step by all hx
                                             iceThickMeltCvOld(:),&   ! ice thickness formed in one time step by all hx
                                             indexStore(:),& !> to which index correspnond the Cv. If two are onvolved we have 2.5 (it means 1/2 the HxCv is in Cv=2 and 1/2 in 3.
                                             dOutIce(:),&      ! outlet diameter of ice growing in a pipe
                                             phiOverlap(:),&   ! overlapping angle= real surface = (2pi-phiOverlap)*r*dl
                                             fConst(:),&       ! overlapping function fConstrained = (2.0*pi-phiOverlap(i))/2.0/pi
                                             dOutIceUsed(:),&  ! outlet diameter of ice growing in a pipe used for UA depending on icing melting or partially icing
                                             areaIceUsed(:),&  !> The total external area of growing ice on the heat exchanger for one CV.
                                             dInIce(:),&  ! inside diameter of ice growing in a pipe ==> dInIce + water = dMeltIce +ice thick = dOutIce
                                             dOutMeltIce(:),&  ! outlet diameter of melted ice 
                                             dOutIceOld(:),&  ! outlet diameter of ice growing in a pipe
                                             dInIceOld(:),&       ! inside diameter of ice growing in a pipe ==> dInIce + water = dMeltIce +ice thick = dOutIce
                                             dOutMeltIceOld(:)  ! outlet diameter of melted ice                                    
                                            
            integer, allocatable :: iceMode(:),&  !> 1 icing mode 0 no icing mode
                                    hxCvBlocked(:) !> ice growing is block 1, 0 not blocked
                                   
                                    
 !> variables to be printed in a file           
            double precision :: Ra,& 
                                Nu,&
                                ReIn,&
                                NuIn,&
                                dOutIceAvg,&
                                dInIceAvg,&
                                dMeltIceAvg,&
                                iceThickAvg,&
                                iceThickInAvg,&
                                iceThickMeltAvg
            
            integer, allocatable :: resetDInIceToZero(:)
            
        end type hxStruct
   
        !type (hxStruct) :: immersedHx(nIHX) ! We must define it globally, otherwise some values are lost for each iteration               
        !logical :: needToSolveHx
        
 end module hxModule 
    
module iceStoreDef

    use iceStoreConst
    use TrnsysConstants
    use TrnsysFunctions
    use hxModule
    
    implicit none        
     
    type iceStoreStruct
        
        ! All parameters must be fixed, but the rest can have a dynamic memory.
                          
        !integer np,ni,nout,nd,luw        
             
        !> Storage Tank variables 
        
        integer :: nCv,&   !> number of storage tank control volumes !!!  
                   iType,& !> TRNSYS type number
                   iUnit,&   !> TRNSYS unit number
                   noticeFound 
        
        character (len=maxMessageLength) :: MyMessage
        
        double precision, allocatable :: Told(:),&   !< Temperature at old time step [K]
                            T(:),&      !< Tank temperature[K]                                                                            
                            tIteGlobal(:),&
                            tIteStore(:),&
                            Tenv(:),&   !< Enviromental temperature  [K]                            
                            vcly(:),&     !< Height of each Cv face (vcly) (vcly(1)=0, vcly(nCv+1)=HTank) [m]
                            nody(:),&     !< Height of each T (nody(1)=(vcly(2)+vcly(1))/2) [m]
                            M(:),&     !< Mass for each control volume (M = volume*rho = [kg])
                            H(:),&     !< Same as dy because Cv heights are always constant [m]
                            A(:),&     !< Exchange area with surrounding  
                            volWater(:),& !< the volume of water for each CV (excluding ice )
                            volWaterOld(:),& !< the volume of water for each CV (excluding ice ) at previous time
                            vTankCv(:),& !< the volume of the tank for each cv                           
                            Uloss(:),&  !< Tank heat transfer losses [W/m^2.K]   
                            UAloss(:),&  !< Tank heat transfer losses [W/K]
                            areaBot(:),& !< Bottom area for each Cv with previous [m^2]. Only change for horizontal cilinder     
                            areaTop(:),& !< Top area for each Cv with previous [m^2]. Only change for horizontal cilinder 
                            qloss(:),&  !< Heat transfer losses for each Cv [W]                            
                            iceFrac(:),&!< Mass percentage of ice fraction [%]                            
                            iceFloatMass(:),& !> Mass of deatached ice in each CV (From top to bottom due to gravity effect)   
                            iceMassHxCv(:),& !> Mass of ice for each hx in each Cv                                                                  
                            iceMassHxCvOld(:),& !> Mass of ice for each hx in each Cv                            
                            iceTotalMassCv(:),&
                            iceFloatMassOld(:),& !> Mass of deatached ice in each CV at previous time step(From top to bottom due to gravity effect)
                            qFusedFloatCv(:),& !> heat of fusion of deatached ice for each CV  
                            qFusedHxCv(:),& !> heat of fusion of ice on Hx for each CV  DC Jan 2015                            
                            !QvIte(:),& !> Source term in energy equation from previous iteration, for relaxation process.
                            !QvHx(:),& !> Source term in energy equation from sensible heat from Hx 
                            isHxLayerRatio(:),& !> [0,1] if hx is inside Cv = 1 else =0, partially = ratio Used to calculate the floating ice on Hx CV's
                            vHxCv(:),&
                            vTankMinusHxCv(:),&
                            storeCvBlocked(:)
        
        double precision :: dy,&           !< Height of each control volume (dy = HTank/nCv) [m]                            
                            VTank,&        !< Tank volume [m3]
                            VTankEff,&     !< Tank volume - heat exchanger volume - ice formed [m3]
                            HTank,&        !< Tank height [m] 
                            WTank,&        !< Tank width [m] or diameter
                            LTank,&        !< Tank lenght [m]
                            !UTank,&        !< Tank heat transfer coefficient U [W/m^2.K]
                            tEnvTop,&      !< Enviromental temperature at the top surface [°C] 
                            tEnvBot,&      !< Enviromental temperature at the bottom surface [°C] 
                            UlossTop,&     !< Additional losses of the top surface [W/m2] 
                            UlossBot,&       !< Additional losses of the bottom surface [W/m2] 
                            areaTopSurface,&  !< the area of the top surface [m2]         
                            areaBotSurface,&  !< the area of the bottom surface [m2]
                            vTankMinusHx      ! the maximum freezable water Tank volume - heat exchanger volume
        
        double precision :: kEff,&         !< Eff. thermal cond. tank [W/m.K]
                            kEffDefined, &
                            alphaEff
        
         double precision :: zSensor1,&         !< Sensor height RELATIVE positions 
                             zSensor2,&         !< Sensor height RELATIVE positions 
                             zSensor3,&         !< Sensor height RELATIVE positions 
                             zSensor4,&         !< Sensor height RELATIVE positions 
                             zSensor5,&         !< Sensor height RELATIVE positions 
                             tSensor1,&         !< T Sensor height RELATIVE positions 
                             tSensor2,&         !< T Sensor height RELATIVE positions 
                             tSensor3,&         !< T Sensor height RELATIVE positions 
                             tSensor4,&         !< T Sensor height RELATIVE positions 
                             tSensor5           !< T Sensor height RELATIVE positions 
         
        !< Immersed heat exchanger related variables 
        
        integer :: nHx,& ! Number of heat exchangers Fixed to 4 !!!
                   nRealHx,& !this will be used for capillar mats to calculate nTubes !
                   mechanicalDeIce,& ! 1: active 0:inactive  
                   iceIsReleased,& ! 1 ice is release in this time step 0 ice remains on the hx surface                   
                   checkParameters,& ! 1 check parameters 0 not                  
                   order(nIHX),& 
                   solutionOrderPositiveMDot(nIHX),&
                   solutionOrderNegativeMDot(nIHX),&
                   solutionOrder(nIHX),&
                   nUsedHx,&
                   nHxStore,&
                   iceBlockGrowingMode
        
        double precision :: tStoreAv,& !> average temperature of the store [oC]
                            qAcumStore,&      !> acumulated heat for a time step [W]                                
                            sumQHx,&
                            sumQLoss,&
                            qFused,&
                            sumQIce,&   !> The sum of all heat exchangers     
                            sumQIceAcum,&   !> The sum of all heat exchangers     
                            qAcumHx,&   !> acumulated heat in heat exchangers
                            fusedIceFloatMass,&
                            imbalance   !> imbalance as imb = sumQHx - qAcumStore - isumQLoss  - qFused + sumQIce   
                           ! nCvInHxs   !> number of Storage Cv that belong to Hxs
                            
        double precision ::  meltCrit,&   !< Film critical melting thickness [m]
                             deltaT,&     !
                             deltaTsub,&
                             maxIceFrac,& !< Maximum allowed fraction of ice formation for each Cv. V_ice/V_available,water
                             iceFloatingIni,&  !
                             iceFloatingOld,&
                             iceFloating,&   !< Mass of ice floating in the store [kg] 
                             iceHx,& 
                             iceTotalMass,& ! !< Mass of ice floating in the store + mass of ice in the Hx = total mass of Ice [kg] 
                             massIceFrac,&   ! >Mass of ice/ mass of (water+ice)
                             massIceFracHx,& ! >Mass of ice/ mass of (water+ice-hx)
                             iceTotalMassOld,& ! !< Mass of ice floating in the store + mass of ice in the Hx = total mass of Ice [kg] 
                             e0,&
                             e1,&
                             H_pc,&   !< ice enthalpy for phase change [J/kg]
                             Tfreeze,&
                             TSubcool,&
                             TBoil,&
                             A_ice,&
                             P_ice,&
                             L_ice,&
                             Ra_ice,&
                             Nu_ice,&
                             rhoIce,&
                             cpIce,&
                             kIce,&
                             h_ice,&
                             deltaEnergy,& !< Change of internal energy due to S [J]
                             !sumOneQHx(nIHX),&
                             qLossTop,&
                             qLossBot, &
                             maxIceFracIceLayer ! maximum ice floating ratio in layers were ice is produced    
        
        integer :: tankGeometry,& !> 0 is a box, 1 vertical cilinder (width = radius), 2 horizontal cylinder (not implemented, only modify calculateGeometry)
                   useTwallOld,& ! 1 for system simulations to improve convergence
                   deIceIsPossible ! de-ice capability or not. Used only on rearrange ice
      
        double precision :: xBetweenPipes,& !> S 0.02 G 0.03
                            yBetweenPipes 
        
        ! Combined ice store and hx variables         
        
        double precision, allocatable :: qhx(:,:),& !< Total HX heat transfer for each CV store (positive is out of the heat exchanger?)
                                         qIceCv(:,:),& !< Energy used for icing or melting in the heat exchangers for each Cv.    
                                         qIceAcumCv(:,:),& !< Sensible ice energy in the heat exchangers for each Cv.  
                                         iceMassHx(:,:) !> Mass of ice for each hx in each CV        
        ! Miscellaneous
        double precision :: RaMeltingFloat,&
                            NuMeltingFloat
        
        double precision :: rhoWater,cpWater                               
        double precision :: timeInHours,&
                            dtsec !> it was an integer before !!!
        
        integer :: indexIceHx,&     ! the index from 1 to 4 of the heat exchangers used to produce ice. From bottom to a height below HTank.
                   storeFull,&      ! if kg ice = 70% of total capacity (ehre 70% is user defined)
                   storeTotallyFull ! (the ice storage is at 95% of capacity, the hx stops working in ice mode)
        
        double precision :: iceThickMax,&   !< Maximum ice thickness allowed before constrained---> storeTotallyFull = true.
                                            ! calculated as ice storage wide / (nhx*2)        
                            volumeTop,&     !< Volume above the heat exllcharchangers used to say storefull = true and then release of ice is not possible           
                            volumeBot,&   !< Volume of layers were the heat exchangers are present (volume=volumeTop+volumeBottom)
                            heightIceHx,&
                            fConstrained    !< constrained factor for touching growing layers
        
        double precision :: errorMaxStorage    = 1e-6 !1e-6  ! lower causes some imbalances
        double precision :: errorMaxIceStorage = 1e-5 !1e-5  ! this used qImb as a criteria
        double precision :: errorMaxHx = 1e-8          ! 1e-8 created no imbalance
        double precision :: tFlip = 4.0d0 ! temperature when inversion occurs           
    !    double precision :: mDotMin = 5./3600d0    !120d0/3600  !  minimum mass flow allowed in heat exchangers, below this number 0 is assumed 120 Kgh = 2 l/min ; 5kg/h / 3600 = 0.0055 kg/s
        double precision :: cNuMelting = 0.15d0   !0.55d0
        double precision :: nNuMelting = 1d0/3d0  !0.33d0
        double precision :: errorStorageLoop,errorIceStorageLoop 
        integer :: meltingPlate  = 1 ! if 1 the melting of floating ice layers is calculated with an equation.
        
        integer :: fatalError = 0    ! 0: no fatal error, 1: ice storage diverged, 2: storage diverged (inside ice storage), 3: hx divered with ice, 4:hx diverged wiout ice
        integer :: nMaxIterHx = 20
        integer :: nMaxIteStorage = 20
        integer :: nMaxIteIceStorage = 20
        integer :: nMaxIteStepByStep = 20 ! I activate to find an error ! 1 !iteration here brings lttle and cost a lot
        integer :: myDebugUnit = 20
        integer :: useInfinitUA = 0
        integer :: printEachIteWarning = 0
        integer :: iUnitFree
        integer :: recalculateInDti = 1
        integer :: useOldTimeStepWhenDiverged = 0 !this is for trnsys system simulation to avoid a case where the model diverge. Old values are used.
        
        !>=================== CHANGE ============!
        integer :: printFileOut  = 0
        integer :: printUAForIce = 1
       
        integer :: verboseLevel = 0 ! verbose level 0 no messatges, 1 only important messages, 3 all messatges, 4 debug mode ((print where is in the code to check if it hangs out) 
        
        integer :: itIceStore,itStore,itTrnsys,itDtTrnsys, dummy
        
        !integer :: allowFloatingIceInIceHx = 1  allow and change parameters if not desired      
        integer :: useAlphaOut = 1 !This needs to be 1 always, otherwise we dont excahnger energy between ice and water
        integer :: fixTwallInIceMode = 0
        double precision :: factortWall=1.0d0               
        
        integer :: constantPhysicalProp,useCorrugatedPlate
                
        
    end type iceStoreStruct
    
    type oneDStruct
                  
         integer :: nCv
         double precision, allocatable :: U(:), &
                                          RhoCpDxOverDt(:),&
                                          Qv(:),&
                                          !QvFromHx(:),&
                                          !QvFromIce(:),&
                                          !QvFromMeltHx(:),&
                                          !QvFromMeltFloating(:),&
                                          QvOld(:),&                  
                                          LambdaOverDx(:),&
                                          UFuse(:),&
                                          ScDx(:),&
                                          SpDx(:)                     
    
    end type oneDStruct
         
    !global variables !!
    
    !type (hxStruct)  :: immersedHx(nIHX)
         
end module iceStoreDef


    !---------------------------------------------------------------------------  
    !> @file calculateIceStore.f90
    !> @author SPF institut fur SolarTechnik 
    !> @author D.Carbonell, W.Lodge
    !> @date 10/08/2012
    !> @brief This file includes the main functions to solve the ice storage tank
    !--------------------------------------------------------------------------  
    module iceStoreCalculate


    use iceStoreDef  
    use hxModule
    use iceStoreConst   
    use Trnsysfunctions

    implicit none

    contains

    !---------------------------------------------------------------------------     
    !> @brief Ice storage calculation
    !> for one time step 
    !> @param[in] iceStore : iceStorage structure
    !> @param[out]  iceStore
    !> @return
    !--------------------------------------------------------------------------- 
    
    subroutine calculateIceStore(iceStore,immersedHx,oneD)

    use spfAlgorithmConst
    use iceStoreFunc
    use hxFunc
    use solvers
    
    implicit none                

    type(iceStoreStruct), intent(inout), target :: iceStore  
    type(oneDStruct), intent(inout), target :: oneD
    type(hxStruct), intent(inout), target :: immersedHx(nIHX)
    
    double precision :: iceFloating, iceFloatingOld        
    double precision , pointer ::   T(:), Told(:), M(:)                     
    integer :: nCv, nHx, conv, iIni,iEnd,iInc,iOrigin
    double precision :: tfreeze, alphaEff, dy
    
    ! internal variables

    integer ::  i, iHx, nMaxIter, status, nIte, dtsec,j !hxMode,hxModeOld
    double precision :: dum, error, errorMax,imb,relax
    logical :: existMassFlow, iceFormation   
    integer :: jHx,jHxBefore, dummy

    
   ! if(iceStore%itDtTRnsys==403) then
   if(iceStore%timeInHours > 8075.83) then
         dummy = 1    
    endif
    
    !< Assignments to iceStore structure

    nCv = iceStore%nCv;  nHx = iceStore%nHx
    iceFloatingOld  = iceStore%iceFloatingOld;   iceFloating  = iceStore%iceFloating
    T => iceStore%T;  Told => iceStore%Told;        
    M => iceStore%M; 
    TFreeze= iceStore%TFreeze; dtsec = iceStore%dtsec
    alphaEff = iceStore%alphaEff; dy = iceStore%dy

    ! End assignments to type iceStorageStruct

    conv = 0
    dum = 0.d0
    nMaxIter = iceStore%nMaxIteIceStorage                       
   
    status = CONTINUE_ITE
    nIte = 0

    iceStore%itIceStore = 0    
    relax = 0.8d0
    
    ! Iterative loop . Main calculation the iceStorage

    existMassFlow = .true.                                

    ! > Initialize T for new iterations 

    iceFloating = iceFloatingOld               
                
    do i=1,nCv
       !T(i) = Told(i) !necessary and correct? To me this is making more difficult the iterative process over several TRNSYS ite
       do j=1,nHx
           iceStore%qIceCv(j,i) = 0.0d0 
           iceStore%qhx(j,i) = 0.0d0
       enddo
    enddo
          
    !Is this necessary? 
    do i=1,nHx
        immersedHx(i)%iceThickMelt= immersedHx(i)%iceThickMeltOld
        immersedHx(i)%iceThickIn  = immersedHx(i)%iceThickInOld
        immersedHx(i)%iceThick    = immersedHx(i)%iceThickOld
    end do        
    !    
    errorMax = iceStore%errorMaxIceStorage 
    
    iIni = 1
    iEnd = nHx
    iInc = 1
    
   ! if(iceStore%revertedFlow(1)) then
   !     iIni = nHx
   !     iInc = -1
   !     iEnd = 1
   ! endif     
    
    if(iceStore%ItDtTrnsys>=5172) then
        iceStore%dummy=1
    endif
    
    do while (status==CONTINUE_ITE .and. nIte<=nMaxIter)      
        
        nIte = nIte+1
        
        iceStore%itIceStore = iceStore%itIceStore+1
        
        iceStore%tIteGlobal(1:nCv) = T(1:nCv)                                           
           
        do iHx=1,iceStore%nHx             
            
            ! We simulate hx in order so that outlet on one can be used as inlet of another.            
            
            jHx = iceStore%solutionOrder(iHx)                      
                        
            !This is the way to control that the Hx is solved
            
            if(immersedHx(jHx)%isUsed) then
                
                if(iHx>1) then
                    jHxBefore = iceStore%solutionOrder(iHx-1)
                    if(immersedHx(jHx)%orderHx>=1 .and. immersedHx(jHxBefore)%orderHx>=1 .and. jHxBefore>0) then                                
                    endif     
                endif
            
                call calculateTStoreSeenByHx(immersedHx(jHx),iceStore)
                ! Make sure that iteration inside StepByStep are not necessary            
                !call calculateStepByStep(iceStore,immersedHx(jHx),iceStore%dtsec) 
                call calculateStepByStep(iceStore,immersedHx,jHx,iceStore%dtsec)
                call heatExchangerStatus(iceStore,immersedHx(jHx)) !are we in iceMode or not?
                call calculateIceHxThickness(iceStore,immersedHx(jHx))
            endif
            
        end do
            
        call calculateQFromHxToStorage(iceStore,immersedHx)               
        
        ! With internal iterations
                        
        call calculateStorage(iceStore,oneD,immersedHx,1)    
        
        if(errorFound()) return
        
        call postProcessIte(iceStore,immersedHx)                              
        
        ! Convergence of this iteration?    

        error = getMaxError(iceStore%tIteGlobal,T,nCv)                      
        imb   = abs(iceStore%sumQHx - iceStore%qAcumStore - iceStore%sumQLoss  - iceStore%qFused + iceStore%sumQIce)
        
        iceStore%errorIceStorageLoop = error
        
        if (iceStore%itTrnsys>=20) then !this means iteration problems....
            iceStore%dummy = 0
        endif
         
        T(1:nCv) = relax*T(1:nCv)+(1-relax)*iceStore%tIteGlobal(1:nCv)                          
         
        iceStore%fatalError = 0
        
        !if (error<errorMax .and. imb<errorMax) then
        if (error<errorMax) then
        !if (error<errorMax) then
            status = CONVERGED                
            
            !print *, 'CONVERGED with ITE ',nIte,'ERROR ',error,' Imb ',imb 
            !write (myScreenUnit,*) 'TIME ',iceStore%timeInSeconds/3600.,' CONVERGED GLOBAL with ite = ',nIte,' ERROR= ',error,' hxMode= ',hxMode
            
        elseif (error>100.) then
            status = DIVERGED           
            
            if(iceStore%useOldTimeStepWhenDiverged==1) then
                call useOldTimeStep(iceStore,immersedHx)          
            end if
            
            iceStore%fatalError = 2
            write(iceStore%MyMessage,'("calculateIceStore GLOBAL ALGORITHM DIVERGED")') 
            call Messages(-1,Trim(iceStore%MyMessage),'NOTICE', iceStore%iUnit,iceStore%iType)
                    
           
        elseif (nIte>nMaxIter) then
            status = MAX_NITER_REACHED    
            
            write(iceStore%MyMessage,*) 'calculateIceStorage MAX NUM ITERATION GLOBAL with ite = ',nIte,' ERROR-T= ',error,' ERROR-Qimb= ',imb    
            if(iceStore%printEachIteWarning==1) then
                call Messages(-1,Trim(iceStore%MyMessage),'NOTICE', iceStore%iUnit,iceStore%iType)                       
            endif
            
        else
            status = CONTINUE_ITE  
            !write (myScreenUnit,*) 'TIME ',iceStore%timeInSeconds/3600.,' CONTINUE ITE ',nIte,'ERROR ',error ,' HXMODE ',hxMode 
        endif                                              
       ! if(iceStore%verboseLevel==4) write (myScreenUnit,*) 'DB(ICE): GLOBAL LOOP nIte',nIte,'ERROR',error,'IMB',imb
    enddo 
        
    
    call reversionEliminationAlgorithm(iceStore) !DC
   
        
    !-----------------------------------------------------------------------
    !                             END MAIN LOOP
    !-----------------------------------------------------------------------

    end subroutine calculateIceStore
        
    !---------------------------------------------------------------------------     
    !> @brief calculation of the storage tank one dimensional heat conduction
    !> equation including a source term of the heat exchanger
    !> @param iceStore : iceStorageSrtruct type
    !> @param oneD : oneD struct      
    !> @return iceStore with calculated data
    !---------------------------------------------------------------------------
     subroutine calculateStorage(iceStore,oneD,immersedHx,recalculateMass)
     
        use solvers
        use spfGlobalConst
        use iceStoreFunc
        use hxFunc
        use spfAlgorithmConst
        use heatTransferCoef        
       
        implicit none
        
        type(iceStoreStruct), intent(inout), target :: iceStore
        type(oneDStruct), intent(inout), target :: oneD
        type(hxStruct), intent(inout), target :: immersedHx(nIHX)
              
        double precision :: relax
        integer :: nCv, nHx, i, nIte,nMaxIter, status, recalculateMass
        double precision :: areaRef,lambdaEff      
        double precision :: error, errorMax        
        double precision :: qSourceHx, relaxQ, lambda, cp, rho       
        double precision :: iceMassOneCv,qHxSensible,qHxLatent !,iceVolumeRatioCv              
        integer :: j,k !DC 14.01.15
               
        
        nCv = iceStore%nCv; nHx = iceStore%nHx
              
        ! We assume all cross area is equal and therefore all coefficients are calculated as W/m2k
        ! being the m2 the cross sectional area. This is important for UA because the external area differs from
        ! the cross sectional area. if we set some areaTop=0 in some Cv it will not affect here
    
        areaRef =  iceStore%areaTop(1)                        
        lambdaEff = iceStore%kEff
        
        status = CONTINUE_ITE
        nIte = 0        
        iceStore%itStore = 0
        
        errorMax = iceStore%errorMaxStorage 
              
        relaxQ = 1.0d0
        relax  = 1.0d0
               
        nMaxIter = iceStore%nMaxIteStorage                
              
        !> This iteration is only necessary if we melt ice because water volume changes or when physical properties of water change.
        !> When no ice exist nMaxIter can be set to 1
        
        if(iceStore%ItDtTrnsys>=5178) then
            iceStore%dummy=1
        endif
        
        do while (status==CONTINUE_ITE .and. nIte<=nMaxIter)                        
                                                
            nIte = nIte+1                                        
            iceStore%itStore = iceStore%itStore + 1
            
           ! iceStore%fusedKgIceFromBulkWater(1:nHx) = 0.0d0
            
            do i=1,nCv
            
                ! Keep this here because the water mass is recalculated for each time step
                if(iceStore%constantPhysicalProp) then
                    cp = iceStore%cpWater                    
                    rho = iceStore%rhoWater
                else    
                    cp = getCpWater(max(iceStore%T(i),0.1))
                    rho = getRhoWater(max(iceStore%T(i),0.1))                        
                endif
               
                oneD%RhoCpDxOverDt(i)  = iceStore%volWater(i)*rho*cp/(iceStore%dtsec*areaRef)  ! W/Km2
                
                if(i /= nCv) then             
                    if(iceStore%constantPhysicalProp) then
                        oneD%LambdaOverDx(i)  = lambdaEff/iceStore%H(i)
                    else
                        lambda = getLambdaWater(0.5*(iceStore%T(i)+iceStore%T(i+1)))
                        oneD%LambdaOverDx(i)  = lambda*lambdaEff/iceStore%H(i)                                             
                    endif
                endif      
                
                iceStore%tIteStore(i) = iceStore%T(i)                                                                                                                    
                qHxSensible = sum(iceStore%qhx(1:nHx,i))
                
                if(sum(immersedHx(1:nHx)%iceMass)>1e-15) then
                    qHxLatent = sum(iceStore%qIceCv(1:nHx,i)) 
                else
                    qHxLatent = 0.0
                end if                                                
                                
                call calculateQvFromMeltHXCv(immersedHx,iceStore,oneD,i,areaRef)   !!calculate iceStore%qFusedHxCv(i)
                
                ! We only melt the old mass of floating ice.
                                                    
                call calculateQMeltFloatingCv(immersedHx,iceStore,oneD,i,areaRef) !calculate iceStore%qFusedFloatCv(i)
                
                oneD%Qv(i) = qHxSensible + qHxLatent - iceStore%qFusedFloatCv(i) -iceStore%qFusedHxCv(i)          
                oneD%Qv(i) = oneD%Qv(i)/areaRef !! W/m2k                                                                                        
                oneD%U(i) = iceStore%UAloss(i)/areaRef                                                                        
                
                !oneD%ScDx(i) = max(oneD%Qv(i),0.0) + oneD%U(i)*iceStore%Tenv(i)    
                !oneD%SpDx(i) = -oneD%U(i) + min(oneD%Qv(i),0.0)/(iceStore%T(i)+1e-10)   
                            
                oneD%ScDx(i) = oneD%Qv(i) + oneD%U(i)*iceStore%Tenv(i)    
                oneD%SpDx(i) = -oneD%U(i)
                
                iceMassOneCv = iceStore%iceFloatMass(i)+sum(iceStore%iceMassHx(1:nHx,i))
                !iceVolumeRatioCv = iceMassOneCv/(iceStore%rhoIce*iceStore%vTankCv(i))
                                                                                           
            enddo                                
                                
            call oneDimensionalConductionPatankar(nCv,iceStore%Told,oneD%ScDx,oneD%SpDx,oneD%LambdaOverDx,oneD%RhoCpDxOverDt,1.0d0,iceStore%T,1)                                                
            !iceStore%qFused = sum(iceStore%qFusedFloatCv(1:nCv))+sum(iceStore%qFusedHxCv(1:nCv)) !DC 15.01.14
                                
            if(recalculateMass) then
                call calculateMassOfIce(iceStore)                    
                call calculateMassOfWater(iceStore,immersedHx)
            endif
                                    
            ! Convergence of this iteration?    

            error = getMaxError(iceStore%tIteStore,iceStore%T,nCv)                                                                                                 
            iceStore%errorStorageLoop = error
            
            !if(iceStore%itTrnsys>=20) then !this means iteration problems....
            !    iceStore%dummy = 0
            !endif
            
            iceStore%fatalError = 0
            
            if (error<errorMax) then
                status = CONVERGED                
                !print *, 'CONVERGED with ITE ',nIte,'ERROR ',error
                !write (myScreenUnit,*) 'TIME ',iceStore%timeInSeconds/3600.,'STORAGE CONVERGED with ite = ',nIte,' ERROR= ',error,' relax= ',relax
            elseif (error>100.0) then
                status = DIVERGED
                iceStore%fatalError = 2
                
                write(iceStore%MyMessage,'("calculateStore CALCULATE STORAGE DIVERGED")') 
                !call Messages(-1,Trim(iceStore%MyMessage),'FATAL', iceStore%iUnit,iceStore%iType)                           
                
            elseif (nIte>nMaxIter) then
                status = MAX_NITER_REACHED   
                
                write(iceStore%MyMessage,*) 'calculateStorage MAX NUM ITERATION ERROR-T= ',error    
                if(iceStore%printEachIteWarning==1) then
                    call Messages(-1,Trim(iceStore%MyMessage),'NOTICE', iceStore%iUnit,iceStore%iType) 
                endif                
            else
                status = CONTINUE_ITE                 
            endif          
                    
            do i=1,nCv
                if(isnan(iceStore%T(i))) then
                    iceStore%dummy=1
                end if
            
                if(iceStore%T(i)<-1) then
                    iceStore%errorStorageLoop = error
                else if(iceStore%T(i)>90.) then
                     iceStore%errorStorageLoop = error
                
                endif
            enddo
            
            iceStore%T(1:nCv) = relax*iceStore%T(1:nCv)+(1-relax)*iceStore%tIteStore(1:nCv)                     
            
            
        enddo       
       
               
     end subroutine calculateStorage
           
     
   
  
  subroutine calculateQvFromMeltHXCv(immersedHx,iceStore,oneD,i,areaRef) ! qv from melting in ice all Hx is calculated
        
    use hxModule
    use iceStoreDef
    use TrnsysFunctions
    use TrnsysConstants
    use heatTransferCoef
    
    implicit none
        
    type(hxStruct), intent(inout) :: immersedHx(nIHX)
    type(oneDStruct), intent(in), target :: oneD
    type(iceStoreStruct), intent(inout) :: iceStore  
    double precision, intent(in) :: areaRef
    integer, intent(in) :: i 
    double precision :: areaIce,alphaOutHx,qUsedToFuseHx,qUsedToFuseHxMax, kgIceInAllParallelHx,fusedKgIce,&
                        Ra, Nu,lChar,FourVolOverPiL,tUsed,possibleMeltedIce,volOuter,qFusedHxCv,&
                        volInner, dum, iceMassHxToCvStore
    integer :: j,k
           
    tUsed = iceStore%TOld(i)
    iceStore%qFusedHxCv(i) =0.0  
           
    if(iceStore%iceMassHxCv(i)>1e-10 .and. tUsed>0.1  .and. iceStore%useAlphaOut==1) then
        
        ! Melting of ice layer on the hx's !                
        ! write(iceStore%MyMessage,*) 'Melting from outside DtIte=',iceStore%itDtTrnsys,' iCv=',i,' itIceStor=',iceStore%itIceStore,'itStore=',iceStore%itStore,' itTrnsy=',iceStore%itTrnsys                
        ! call Messages(-1,Trim(iceStore%MyMessage),'NOTICE', iceStore%iUnit,iceStore%iType) 
        
        do j=1,iceStore%nHx                         
            
            if(immersedHx(j)%isUsed) then                                                               
                do k=1,immersedHx(j)%numberOfCv 
                         
                    immersedhx(j)%iceMassMeltedOutsideCv(k,i) = 0.0d0
                    
                    !If I'm melting from inside I don't melt from outside, so when qIceCv<0. I'm melting !!
                    if(immersedHx(j)%iceMassCv(k)*immersedHx(j)%factorHxToTank(k,i)>1e-10 ) then                                            
                        
                        if(immersedHx(j)%geometry==PLATE) then     
                            
                            areaIce = immersedHx(j)%areaIceUsed(k)*immersedHx(j)%factorHxToTank(k,i)*immersedHx(j)%nParallelHx
                            
                            call calculateHeatTransCoefImmersedPlate(0.5*(iceStore%Tfreeze+tUsed),tUsed,iceStore%Tfreeze,immersedHx(j)%Lhx,&
                            alphaOutHx,iceStore%cNuMelting,iceStore%nNuMelting,Ra,Nu,iceStore%iUnit,iceStore%iType)  
                            
                        else                            
                            ! WE ARE ONLY MELTING THE OUTSIDE ICE LAYER !!                         
                            areaIce = pi*immersedHx(j)%dOutIce(k)*(immersedHx(j)%dL)*immersedHx(j)%factorHxToTank(k,i)*immersedHx(j)%nParallelHx                                                       
                            !lChar   = immersedHx(j)%LHx/immersedHx(j)%numberOfCv
                            !call calculateHeatTransCoefPipeOut(immersedHx(j)%tFilm(k),iceStore%Told(i),immersedHx(j)%tWall(k),lChar,alphaOutHx)                                                          
                            lChar   = immersedHx(j)%dOutIce(k)
                            call calculateHeatTransCoefImmersedPlate(0.5*(iceStore%Tfreeze+tUsed),tUsed,iceStore%Tfreeze,immersedHx(j)%Lhx,&
                            alphaOutHx,iceStore%cNuMelting,iceStore%nNuMelting,Ra,Nu,iceStore%iUnit,iceStore%iType)  
                            
                        endif
                        
                        qUsedToFuseHx = alphaOutHx*areaIce*max(tUsed-iceStore%Tfreeze,0.0)  ! W                                                         
                        qUsedToFuseHxMax = oneD%RhoCpDxOverDt(i)*areaRef*(iceStore%T(i)-iceStore%Tfreeze)*immersedHx(j)%factorTankToHx(k,i)-iceStore%qFusedHxCv(i)
                        
                        qUsedToFuseHx = min(qUsedToFuseHx,qUsedToFuseHxMax)     
                                      
                        ! The kg melted from existent kg
                        ! kgIceInAllParallelHx = immersedHx(j)%iceMassCv(k)*immersedHx(j)%factorHxToTank(k,i)
                        
                        if(immersedHx(j)%geometry==COIL) then  
                            volOuter =  pi*immersedHx(j)%dL*(immersedHx(j)%dOutIce(k)**2-immersedHx(j)%dOutMeltIce(k)**2)/4.0
                            kgIceInAllParallelHx = volOuter*immersedHx(j)%factorHxToTank(k,i)*immersedHx(j)%nParallelHx*iceStore%rhoIce
                        else
                            kgIceInAllParallelHx = immersedHx(j)%iceThickCv(k)*immersedHx(j)%areaIceUsed(k)*immersedHx(j)%factorHxToTank(k,i)*iceStore%rhoIce
                        endif                                                
                                                                                                                             
                        possibleMeltedIce = max(qUsedToFuseHx*iceStore%dtsec/iceStore%H_pc,0.0)                            
                        fusedKgIce = min(possibleMeltedIce,kgIceInAllParallelHx)  ! W*s/J/kg = kg   
                        
                        if(kgIceInAllParallelHx > possibleMeltedIce) then
                            !limited !!!
                            iceStore%dummy=1
                        endif
                                              
                        ! Real fusion heat considering all hx's
                        
                        iceStore%qFusedHxCv(i) = iceStore%qFusedHxCv(i) + fusedKgIce * iceStore%H_pc / iceStore%dtsec                                                        
                        
                        immersedhx(j)%iceMassMeltedOutsideCv(k,i) = fusedKgIce                                                  
                                                                    
                    else
                        
                    endif                                      
                    
                enddo !HX loop i=1 nHxCv                             
                
            endif
        enddo !for hx=1 nHX
    else                                       
        iceStore%qFusedHxCv(i) = 0.0   
        do j=1,iceStore%nHx   
            if(immersedHx(j)%isUsed) then                                                               
                do k=1,immersedHx(j)%numberOfCv 
                    immersedhx(j)%iceMassMeltedOutsideCv(k,i) = 0.0d0
                enddo
            endif
        enddo       
    end if
       
   
    
  end subroutine calculateQvFromMeltHXCv

! Called only once at end of time Step

                        
   subroutine calculateQMeltFloatingCv(immersedHx,iceStore,oneD,i,areaRef) ! qv from melting floating ice is calculated
        
        use hxModule
        use iceStoreDef
        use heatTransferCoef
    
        implicit none
        
        type(hxStruct), intent(in) :: immersedHx(nIHX)
        type(iceStoreStruct), intent(inout) :: iceStore    
        type(oneDStruct), intent(inout), target :: oneD
        double precision, intent(in) :: areaRef
        integer, intent(in) :: i
        !double precision, intent(out) :: qMeltFloating(iceStore%nCv)                              
        double precision :: iceThickAv, areaIcePlate,areaExchange,alphaOut,qUsedToFuse,qUsedToFuseMax,&
                            fusedKgIce,fAreaBroken,myTStorage
            
        fAreaBroken = 1.0
        qUsedToFuse = 0.0                                             
        qUsedToFuseMax = 0.0                        
        
        if(iceStore%iceFloatMassOld(i)>0.0) then  
        
            !> We need to calculate each time the volume of water is calculated, but not every iteration becasue Told is used
            
            myTStorage = max(iceStore%Told(i),0.0d0)
            
            if(iceStore%meltingPlate) then                   
                ! I assume that all layers are of 1cm thick. From Incropera 2006 cites in Tamasauskas solar energy 2012.
                ! AreaExchange = iceStore%iceFloatMass(i)*fAreaBroken*2.0d0/(iceStore%rhoIce*0.01)                        
                iceThickAv   = 0.01 ! 1cm
                areaIcePlate = 0.5  ! m2                                                
                areaExchange = (iceStore%iceFloatMass(i)*fAreaBroken*areaIcePlate/(iceStore%rhoIce*iceThickAv))*(2.0*areaIcePlate+4.0*sqrt(areaIcePlate)*iceThickAv)                        
                call calculateHeatTransCoefImmersedPlate(0.5*(iceStore%Tfreeze+myTStorage),myTStorage,iceStore%Tfreeze,immersedHx(1)%Lhx,alphaOut,iceStore%cNuMelting,iceStore%nNuMelting,iceStore%RaMeltingFloat,iceStore%NuMeltingFloat,iceStore%iUnit,iceStore%iType)                                                                                                    
                qUsedToFuse    = alphaOut*areaExchange*max(myTStorage-iceStore%Tfreeze,0.0d0)   ! W                                              
                qUsedToFuseMax = oneD%RhoCpDxOverDt(i)*areaRef*(myTStorage-iceStore%Tfreeze)                                                                      
                qUsedToFuse = min(qUsedToFuse,qUsedToFuseMax)                       
            else
                qUsedToFuse = oneD%RhoCpDxOverDt(i)*areaRef*(myTStorage-iceStore%Tfreeze) ! W   CHANGED 05.08.2014    I ONLY METL WITH TOLD !!!!!!! 
            endif
            
            ! The kg melted f of kg existent
            fusedKgIce = max(min(qUsedToFuse*(iceStore%dtsec/iceStore%H_pc),iceStore%iceFloatMass(i)),0.0)  ! W*s/J/kg = kg                                           
            ! Real fusion heat 
            iceStore%qFusedFloatCv(i)= fusedKgIce * iceStore%H_pc / iceStore%dtsec                 
        else            
            iceStore%qFusedFloatCv(i)= 0.0d0
        end if            
            
        
        
   end subroutine calculateQMeltFloatingCv
                                     
                
    end module iceStoreCalculate
    
       
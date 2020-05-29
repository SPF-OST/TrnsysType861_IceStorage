
    
subroutine TYPE861 
!***********************************************************************
!                                              ____    ____    ______  *
!                                             /  _ \  / __ \  / ____/  *
!                                             \ \ \/ / /_/ / / /___    *
!                                              \ \  / ____/ / ____/    *
! Institut für                             /\__/ / / /     / /         *
! Solartechnik                             \____/ /_/     /_/          *
! Copyright © 2015                                                     *
! @author D.Carbonell                                                  *
!***********************************************************************
!
! DESCRIPTION
!
!     Stratification storage tank wihout direct ports.
!     4 Immersed heat exchangers.
!     Heat exchangers can form ice. 
!
! REVISION HISTORY
!     version v30 20-08-2015. Type860 moved to Type861 with immersedHx that are discretized and solved with step by step model
!
! TODO 
!   - Ice sensible capacity should be considered when ice grows. 
!   - Maximum ice fraction should consider that water moves when ice grows, then only the volume (and not the mass) covered by the hx can be iced
!   - Now it is simplified from inputs such that rho,ice = rho,water, but a coherent formulation should de done
!   - This mode can be easily extended to other PCM's
! DETECTED PROBLEMS
!  - In some cases a negative water volume was detected until divergence occurs. This should never happen when the maximum ice fraction is set below 0.9 or so, but it still happens some times.
!  - A possible solution could be to avoid adding ice to a CV which is mostly covered by a blocked CV.
!  
! DIMENSIONS
!     All calculations inside TYPE861 are done in SI units.
!
! REFERENCE
! The use of this model should be referenced as :
! 
! for capillary mats formulation and validation
!
! D. Carbonell, M. Battaglia, D. Philippen, M.Y. Haller 
! "Numerical and experimental evaluation of ice storages with ice on capillary mat heat exchangers for solar-ice systems"
! International Journal of Refrigeration, Volume 88, April 2018, Pages 383-401.
! 
! and the following for the flat plate formulation and validation
! 
! Carbonell, et al., 2017. Ice-Ex “Heat Exchanger Analyses for Ice Storages in Solar and Heat Pump Applications”. 
! Institut für Solartechnik for Swiss Federal Office of Energy (SFOE), Research Programme Solar Heat and Heat Storage, CH-3003 Bern.
! 
!***********************************************************************
!DEC$ATTRIBUTES DLLEXPORT :: TYPE861

!     USER STATEMENTS
      use TrnsysConstants
      use TrnsysFunctions

      use spfGlobalConst
      use iceStoreConst
      use iceStoreDef
      use hxModule
      use hxFunc
      
      use iceStoreFunc
      use iceStoreCalculate
     
      use interpolation      

      implicit none
             
      character (len=12) iStr
      character (len=maxMessageLength) TypeMsg(10),TypeMsgFatal(5)  
      integer :: iunit,itype,np,ni,nout,nd,nCv,i,timeStepIt 
      double precision :: simulationTime 
      double precision :: time0, dt, tInit(100)
      
      integer :: errorCheck,iteration,iHx,k
      
      type(iceStoreStruct) :: iceStore
      type(oneDStruct) :: oneD
      type (hxStruct)  :: immersedHx(nIHX)    
      
!-----------------------------------------------------------------------
!                RETRIEVE VALUES OF GLOBAL VARIABLES
!-----------------------------------------------------------------------      

      simulationTime = getSimulationTime()
      time0 = getSimulationStartTime()
      dt = getSimulationTimeStep()
      timeStepIt = getTimestepIteration()
      
      iceStore%itTrnsys = timeStepIt                      
      iceStore%timeInHours = simulationTime
                
      if(iceStore%verboseLevel==4) then 
          write (iceStore%myDebugUnit,*) '========================'
          write (iceStore%myDebugUnit,*) 'DB(ICE): INSIDE TYPE 861 ',iceStore%timeInHours
          write (iceStore%myDebugUnit,*) '========================'
      endif
      
!-----------------------------------------------------------------------
!     SET THE VERSION INFORMATION FOR TRNSYS
!-----------------------------------------------------------------------     

      if(getIsVersionSigningTime()) then
          call SetTypeVersion(17)
          return 
      endif
     
!-----------------------------------------------------------------------
!                  PERFORM LAST CALL MANIPULATIONS 
!-----------------------------------------------------------------------
    
      if(getIsLastCallofSimulation()) then         
          
        call setMemoryFreeIceStorage(iceStore) 
        call setMemoryFreeOneDStruct(oneD)
        call setMemoryFreeHx(immersedHx,iceStore%nHx)
        
         if(iceStore%printFileOut==1) then
            CLOSE(iceStore%iUnitFree)  
         endif
         
        return 
        
      endif

!------------------------------------------------------------------------------
! PERFORM POST CONVERGENCE MANIPULATIONS. 
! Perform Any "End of Timestep" Manipulations That May Be Required
!------------------------------------------------------------------------------       
      
      if(getIsEndOfTimestep()) then      
                                                                  
        if(iceStore%verboseLevel==4) write (iceStore%myDebugUnit,*) 'DB(ICE):BEGIN IS END OF TIME STEP'
                 
        if(iceStore%timeInHours>=6) then
            iceStore%dummy=1
        endif 
        
        call releaseIce(iceStore,immersedHx)
        call reDistributeIce(iceStore)        
        call checkStorageStatus(iceStore,immersedHx)               
        call calculateMassOfWater(iceStore,immersedHx)     
        call postProcessDti(iceStore,immersedHx)                    
        call printMessagesLastTimeStep(iceStore,immersedHx)        
        call updateTimeStep(iceStore,immersedHx)          
        call printMyOutFile(iceStore,immersedHx)               
        
        if(iceStore%verboseLevel==4) then 
          write (iceStore%myDebugUnit,*) '========================'
          write (iceStore%myDebugUnit,*) 'DB(ICE): IS END OF TIME STEP TYPE 861 ',iceStore%timeInHours
          write (iceStore%myDebugUnit,*) '========================'
        endif
        
        !if(iceStore%timeInHours>5842.) then
        if(immersedHx(1)%mDotAllHxInKgh>10) then
            !immersedHx(1)%fConstrainedAvg>=0.5
            iceStore%dummy=0
        endif
        
        call setOutputs(iceStore,immersedHx)
        
        iceStore%itDtTrnsys = iceStore%itDtTrnsys+1                                
        
      endif
 
!-----------------------------------------------------------------------
!                  PERFORM FIRST CALL MANIPULATIONS
!-----------------------------------------------------------------------
     If(getIsFirstCallofSimulation()) then
          
!       Retrieve unit and type number

        iceStore%iUnit = getCurrentUnit()
        iceStore%iType = getCurrentType()
        iceStore%iUnitFree = getNextAvailableLogicalUnit()
        
!       Retrieve the number of derivatives (Tinit)

        nCv = getNumberOfDerivatives()      
        nd = nCv
        np = 22+(nIHX*19)+2+4+5+2 !+4 (top,bottom topCv and botCv)
        ni = nCv+2+(nIHX*3)+1    ! +1 mechanicalDeIce +2 TTop and tbot
        nout = 40+2*nCv+2+5+5+10+1+4        
        
        call SetNumberofParameters(np)          
        call SetNumberofInputs(ni)	
        call SetNumberofDerivatives(nd)        
        call SetNumberofOutputs(nout)          
        call SetIterationMode(1)  
                  
        if(nCv<1 .or. nCv>nMax) then     
          write(iceStore%MyMessage,'("The number of storage Cv =",d"  must be between 1 and =",d)') nCv,nMax
          call Messages(-1,Trim(iceStore%MyMessage),'FATAL', iceStore%iUnit,iceStore%iType)                    
        endif
        
        !> Set dynamic memory 
        
        iceStore%nCv =nCv                
        
        call setMemoryIceStorage(iceStore)
        call setMemoryOneDStruct(oneD,nCv)        
        
        do iHx=1,nIHX
            call setDefaultHx(immersedHx(iHx))
        end do
        
        immersedHx(1)%numberOfCv = getParameterValue(40)
        immersedHx(2)%numberOfCv = getParameterValue(40+19)
        immersedHx(3)%numberOfCv = getParameterValue(40+19*2)
        immersedHx(4)%numberOfCv = getParameterValue(40+19*3)     
        
        immersedHx(1)%geometry = getParameterValue(27)
        immersedHx(2)%geometry = getParameterValue(27+19)
        immersedHx(3)%geometry = getParameterValue(27+19*2)
        immersedHx(4)%geometry = getParameterValue(27+19*3)     
        
        
        iceStore%nUsedHx = 0
        
        
        do iHx=1,nIHX 
            if(immersedHx(iHx)%numberOfCv>0) then
                immersedHx(iHx)%isUsed = 1
                iceStore%nUsedHx = iceStore%nUsedHx+1
               
            else 
                immersedHx(iHx)%isUsed = 0
            endif
        enddo
        
        
        
        iceStore%nHx = getParameterValue(1) 
        
        call setMemoryHx(immersedHx,nCv,nIHX)
        
         if(iceStore%printFileOut==1) then
            OPEN (iceStore%iUnitFree,File = "UA.dat")     
            if(iceStore%printUAForIce==1) then
                if(immersedHx(1)%geometry==COIL) then
                     write(iceStore%iUnitFree,*) 'time UAhx hIn hWall hout hIce fcons phi dIceIn dIceOut dMelt Nu AR dice/small'
                else
                     write(iceStore%iUnitFree,*) 'time UAhx hIn hWall hout hIce fcons phi sIceInAvg sIceAvg sMeltAvg Nu - -'
                endif
               
                write(iceStore%iUnitFree,*) 'Index  W/K   W/m2k   W/m2k W/m2k [] [rad] [cm] [cm] [cm] [-] [-] [-]'                                 
            else
                write(iceStore%iUnitFree,*) 'time   UAhx hIn hWall hout hIce Ra Nu Rei Nui'
                write(iceStore%iUnitFree,*) 'Index  W/K   W/m2k   W/m2k W/m2k [] [] [] [] [] [] []'                                 
            endif
         endif
         
                 
        return 
        
      endif
                  
      if(getIsReReadParameters()) then     
          write(iceStore%MyMessage,'("More that one type is not allowed. You can compile two dll with different name instead")') 
          call Messages(-1,Trim(iceStore%MyMessage),'FATAL', iceStore%iUnit,iceStore%iType)   
          return          
      endif
      
!-----------------------------------------------------------------------
!               PERFORM INITIAL TIMESTEP MANIPULATIONS             
!-----------------------------------------------------------------------
     
      if (getIsStartTime()) then                                          
        
        do i=1,nd        
            tInit(i) = getNumericalSolution(i)
        enddo
        
        call initialize(iceStore,tInit) 
        
        !-----------------------------------------------------------------------
        !                     READ PARAMETERS (serially):
        !-----------------------------------------------------------------------               
        
        call readParameters(iceStore,immersedHx)
        
        if(iceStore%checkParameters) then
            call checkParameters(iceStore,immersedHx)
        endif 
        
        if(errorFound()) then
            return                                              
        end if
        
        call calculateGeo(iceStore,immersedHx)                        
        call calculateGeoHx(immersedHx,iceStore) 
        
        call initializeMassOfIce(iceStore,immersedHx)
        
        call initializeTempHx(immersedHx,iceStore) !Hx temp = tStorage
        
        call updateTimeStep(iceStore,immersedHx)                 
          
      endif
      
!-----------------------------------------------------------------------
!                    MAIN CALCULATION/ITERATION LOOP
!-----------------------------------------------------------------------

    !>-----------------------------------------------------------------------
    !>                            READ INPUTS
    !>-----------------------------------------------------------------------
                                                          
      call readInputs(iceStore,immersedHx)
      
      if(iceStore%ItDtTrnsys==9126) then
          iceStore%dummy=1
      endif
       
      if(iceStore%TimeInHours>=8320.233) then
          iceStore%dummy=1
      endif
      
      if(errorFound()) return
      
    !>-----------------------------------------------------------------------
    !>                        TIME RELATED THINGS
    !>-----------------------------------------------------------------------
    !>     Time-step in seconds

      if (simulationTime <= (time0 + 1e-5)) then
         iceStore%dtsec = 0d0
         do i=1,nout
            call setOutputValue(i,0.0d0)
         enddo             
         return          
      else
         iceStore%dtsec = dt*3600d0 
      endif                

      call calculateIceStore(iceStore,immersedHx,oneD) 

      if(errorFound()) return
    
      if(iceStore%fatalError /= 0 .and. iceStore%useOldTimeStepWhenDiverged==0) then                        
        call Messages(-1,iceStore%fatalError,'FATAL', iceStore%iUnit,iceStore%iType)   
        return            
        !call FoundBadParameter(1,'Fatal',TypeMsgFatal(iceStore%fatalError))  
      endif
!-----------------------------------------------------------------------
!                           SET OUTPUTS
!-----------------------------------------------------------------------

     call setOutputs(iceStore,immersedHx)       
      
     do i=1,nIHX
        if(immersedHx(i)%isUsed==0) then
            do k=1,immersedHx(i)%numberOfCv 
                immersedHx(i)%dInIce(k)      = immersedHx(i)%dInIceOld(k)
                immersedHx(i)%dOutMeltIce(k) = immersedHx(i)%dOutMeltIceOld(k)
            enddo
        endif
     enddo
    
     
     if(iceStore%itTrnsys>=20) then ! This means iteration problems....
          iceStore%dummy = 0
     endif
      
     if(iceStore%timeInHours>=24.466) then ! This means iteration problems....
          iceStore%dummy = 0
     endif
     
     if(iceStore%verboseLevel==4) then 
          write (iceStore%myDebugUnit,*) '========================'
          write (iceStore%myDebugUnit,*) 'DB(ICE): END OF TYPE 861 ',iceStore%timeInHours
          write (iceStore%myDebugUnit,*) '========================'
       endif
     return

end SUBROUTINE TYPE861

subroutine printMyOutFile(iceStore,immersedHx)
    use iceStoreDef
    use hxModule
    use interpolation
    use TrnsysFunctions
     
    implicit none
   
    type(iceStoreStruct), intent(inout), target :: iceStore  
    type(hxStruct), intent(inout), target :: immersedHx(nIHX)
     
    integer i

    double precision :: dIceIn,dIceOut,dIceMelt,sumIceThick
    double precision :: Ra, Ar, long, small, dIceOverSmall 
    
    if(iceStore%printFileOut==1) then          
  
        i = 1   
                        
        if(immersedHx(i)%geometry==COIL) then
            dIceIn   = 100*immersedHx(i)%dInIceAvg   ! cm
            dIceOut  = 100*immersedHx(i)%dOutIceAvg  ! cm    
            dIceMelt = 100*immersedHx(i)%dMeltIceAvg ! cm
            
            long  = max(iceStore%xBetweenPipes,iceStore%yBetweenPipes)*100 ! cm
            small = min(iceStore%xBetweenPipes,iceStore%yBetweenPipes)*100 ! cm
            Ar = 1.0 - 2*(acos(min(long/dIceOut,1.0)) + acos(min(small/dIceOut,1.)))/pi            
            dIceOverSmall = dIceOut/small
        else
            dIceIn   = 100*immersedHx(i)%iceThickInAvg
            dIceOut  = 100*immersedHx(i)%iceThickAvg
            dIceMelt = 100*immersedHx(i)%iceThickMeltAvg
            Ar = 1
            dIceOverSmall = 1
        endif
        
        if(iceStore%printUAForIce==1) then
            if(iceStore%timeInHours==1) then
                iceStore%dummy =1
            endif
            
        !write(iceStore%iUnitFree,11) iceStore%timeInHours,immersedHx(i)%UA,immersedHx(i)%alphaIn,immersedHx(i)%alphaWall,immersedHx(i)%alphaOut
            write(iceStore%iUnitFree,11) iceStore%timeInHours,immersedHx(i)%UA*immersedHx(i)%nParallelHx,immersedHx(i)%UaIn*immersedHx(i)%nParallelHx,immersedHx(i)%UaWall*immersedHx(i)%nParallelHx,min(immersedHx(i)%UaOut*immersedHx(i)%nParallelHx,1e5),&
            min(immersedHx(i)%UaIce*immersedHx(i)%nParallelHx,1e5),immersedHx(i)%fConstrainedAvg,immersedHx(i)%phiOverlapAvg,dIceIn,dIceOut,dIceMelt,immersedHx(i)%Nu,Ar,dIceOverSmall
     
             11   format(F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.8,1x,F12.8)
        else
            write(iceStore%iUnitFree,12) iceStore%timeInHours,immersedHx(i)%UA*immersedHx(i)%nParallelHx,immersedHx(i)%UaIn*immersedHx(i)%nParallelHx,immersedHx(i)%UaWall*immersedHx(i)%nParallelHx,immersedHx(i)%UaOut*immersedHx(i)%nParallelHx,&
            min(immersedHx(i)%UaIce*immersedHx(i)%nParallelHx,10000.),immersedHx(i)%Ra,immersedHx(i)%Nu,immersedHx(i)%ReIn,immersedHx(i)%NuIn  
            Ra = immersedHx(i)%Ra
          !  12   format(F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4)
          12   format(F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,F12.4,1x,ES,1x,ES,1x,F12.4,1x,F12.4,1x,F12.4)
                
        endif
        
        !11   format(E5.6,1x,E5.6,1x,E5.6,1x,E5.6,1x,E5.6,1x,E5.6,1x,E5.6,1x,E5.6,1x,E5.6,1x,E5.6,1x,E5.6)
       
             
    endif
end subroutine printMyOutFile

    
subroutine setOutputs(iceStore,immersedHx)
          
    use iceStoreDef
    use hxModule
    use interpolation
    use TrnsysFunctions
     
    implicit none
   
    type(iceStoreStruct), intent(inout), target :: iceStore  
    type(hxStruct), intent(inout), target :: immersedHx(nIHX)
     
    integer i,k,oInc,nCv
    double precision :: QLoss, USideAv,myQLoss
    double precision :: tZ01,tZ02,tZ03,tZ04,tZ05,tZ06,tZ07,tZ08,tZ09,tZ10,sumIceThick
    !double precision :: dIceIn,dIceOut,dIceMelt
    !double precision :: Ra 
    
    nCv = iceStore%nCv        
    
    call setOutputValue(1,iceStore%tStoreAv)
             
    ! Storage timestep internal energy change (due to S!) [W]
    ! iceStore%imbalance = iceStore%sumQHx - iceStore%qAcumStore - iceStore%sumQLoss  - iceStore%qFused + iceStore%sumQIce
    ! iceStore%sumQHx,iceStore%qAcumStore,iceStore%sumQLoss,iceStore%qFused,iceStore%sumQIce,iceStore%imbalance
      
    call setOutputValue(2,iceStore%sumQHx)
    call setOutputValue(3,iceStore%qAcumStore)
    call setOutputValue(4,iceStore%sumQLoss)
      
    if(iceStore%qFused>0.0 .and. iceStore%qFused<1e-20) then
        call setOutputValue(5,0.0d0)
    else
        call setOutputValue(5,iceStore%qFused)
    endif
     
    call setOutputValue(6,iceStore%sumQIce)
    call setOutputValue(7,iceStore%imbalance)
                 
    ! Amount of floating ice in storage tank [kg]
    
    if(iceStore%iceFloating<1e-20) then
        call setOutputValue(8,0.0d0)
    else
        call setOutputValue(8,iceStore%iceFloating)
    endif      
       
        
    if(sum(immersedHx(1:iceStore%nHx)%iceThick)<1e-20) then                    
        sumIceThick = 0.0d0 ! thickness of ice layer in total [m]
    else        
        sumIceThick = sum(immersedHx(1:iceStore%nHx)%iceThick)       
    endif

    if(iceStore%timeInHours>=10.0) then
        nCv = iceStore%nCv
    endif
    
     call setOutputValue(9,sumIceThick)
     
    !total amount of ice  
    if(iceStore%iceTotalMass<1e-20) then
        call setOutputValue(10,0.0d0)
    else
        call setOutputValue(10,iceStore%iceTotalMass)
    endif
                
    !     The storage is full with ice to the limit set by the user    [0:1]
    call setOutputValue(11,iceStore%storeFull*1d0)
            
    oInc=12
    !     Output for all IHXs, even those inactive:
    do i=1,nIHX
        if(immersedHx(i)%isUsed==0) then
            do k=oInc,oInc+7
                call setOutputValue(k,0.d0)
            enddo
        else    
            
        !       HX inlet temperature:                                      [°C]
            call setOutputValue(oInc,immersedHx(i)%tFluidIn); oInc=oInc+1
        !       HX outlet temperature:                                      [°C]                
            call setOutputValue(oInc,immersedHx(i)%tFluidOut); oInc=oInc+1                                               
        !       HX wall temperature:                                      [°C]
            !call setOutputValue(oInc,immersedHx(i)%tWallAv); oInc=oInc+1
            call setOutputValue(oInc,immersedHx(i)%fConstrainedAvg); oInc=oInc+1
        !!      HX mass flow for all parallel hx: [kg/h]
        !!    call setOutputValue(oInc,immersedHx(i)%mDot*3600d0*immersedHx(i)%nParallelHx); oInc=oInc+1
        !       HX energy exchanged:                                        [W]
            !!call setOutputValue(oInc,immersedHx(i)%qHxToTnk); oInc=oInc+1   
            call setOutputValue(oInc,immersedHx(i)%qUtil); oInc=oInc+1 ! 27.01.2016 From the heat exchanger perspective. NOT from the storage !!!!!
            
            if(immersedHx(i)%geometry==PLATE) then
                !       HX ice thickness including all parallel. Since we use 2 times the area the value represents the thickness in one side in flat plate and radious thickness in cilinder [m]
                call setOutputValue(oInc,immersedHx(i)%iceThick); oInc=oInc+1             
                !       HX ice thickness including all parallel. Since we use 2 times the area the value represents the thickness in one side in flat plate and radious thickness in cilinder [m]
                call setOutputValue(oInc,immersedHx(i)%iceThickMelt); oInc=oInc+1                                              
            else                                
                call setOutputValue(oInc,sum(immersedHx(i)%dOutIce(1:immersedHx(i)%numberOfCv))/immersedHx(i)%numberOfCv); oInc=oInc+1            
                call setOutputValue(oInc,sum(immersedHx(i)%dOutMeltIce(1:immersedHx(i)%numberOfCv))/immersedHx(i)%numberOfCv); oInc=oInc+1                                
            endif            
            !       UA HX  including all parallel      [W/mK]
            call setOutputValue(oInc,immersedHx(i)%UA*immersedHx(i)%nParallelHx); oInc=oInc+1   
            
            if(iceStore%itDtTrnsys>=5965) then
                iceStore%dummy=1
            endif
        endif                
    enddo

    call setOutputValue(40,iceStore%storeTotallyFull*1d0)
      
    oInc=40+1
            
    ! % 5 sensors
      
    iceStore%Tsensor1 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,iceStore%zSensor1,0)
    iceStore%Tsensor2 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,iceStore%zSensor2,0)
    iceStore%Tsensor3 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,iceStore%zSensor3,0)
    iceStore%Tsensor4 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,iceStore%zSensor4,0)
    iceStore%Tsensor5 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,iceStore%zSensor5,0)
      
    call setOutputValue(oInc,iceStore%Tsensor1);oInc = oInc + 1
    call setOutputValue(oInc,iceStore%Tsensor2);oInc = oInc + 1
    call setOutputValue(oInc,iceStore%Tsensor3);oInc = oInc + 1
    call setOutputValue(oInc,iceStore%Tsensor4);oInc = oInc + 1
    call setOutputValue(oInc,iceStore%Tsensor5);oInc = oInc + 1
        
    tZ01 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.05*iceStore%HTank,0)
    tZ02 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.1*iceStore%HTank,0)
    tZ03 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.2*iceStore%HTank,0)
    tZ04 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.3*iceStore%HTank,0)
    tZ05 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.4*iceStore%HTank,0)
    tZ06 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.5*iceStore%HTank,0)
    tZ07 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.6*iceStore%HTank,0)
    tZ08 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.7*iceStore%HTank,0)
    tZ09 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.8*iceStore%HTank,0)
    tZ10 = vectorLinearInterpolationLimited(iceStore%nCv,iceStore%nody,iceStore%T,0.95*iceStore%HTank,0)
      
    call setOutputValue(oInc,tZ01);oInc = oInc + 1
    call setOutputValue(oInc,tZ02);oInc = oInc + 1
    call setOutputValue(oInc,tZ03);oInc = oInc + 1
    call setOutputValue(oInc,tZ04);oInc = oInc + 1
    call setOutputValue(oInc,tZ05);oInc = oInc + 1
    call setOutputValue(oInc,tZ06);oInc = oInc + 1
    call setOutputValue(oInc,tZ07);oInc = oInc + 1
    call setOutputValue(oInc,tZ08);oInc = oInc + 1
    call setOutputValue(oInc,tZ09);oInc = oInc + 1
    call setOutputValue(oInc,tZ10);oInc = oInc + 1      
      
    i = 1
    iceStore%qLossBot = iceStore%qLoss(i)-iceStore%ULoss(i)*iceStore%A(i)*(iceStore%T(i) - iceStore%TEnv(i))      
    i = iceStore%nCv
    iceStore%qLossTop = iceStore%qLoss(i)-iceStore%ULoss(i)*iceStore%A(i)*(iceStore%T(i) - iceStore%TEnv(i))
      
    call setOutputValue(oInc,iceStore%qLossBot);oInc = oInc + 1
      
    do i=1,iceStore%nCv
        if(i==1) then  
            myQLoss = iceStore%qLoss(i)-iceStore%qLossBot
        else if (i==iceStore%nCv) then
            myQLoss = iceStore%qLoss(i)-iceStore%qLossTop
        else
            myQLoss = iceStore%qLoss(i)
        endif        
        call setOutputValue(oInc,myQLoss)        
        oInc = oInc + 1
    enddo
      
    call setOutputValue(oInc,iceStore%qLossTop);oInc = oInc + 1           
    call setOutputValue(oInc,iceStore%Tenv(1))           ;oInc = oInc + 1
    call setOutputValue(oInc,iceStore%Tenv(iceStore%nCv));oInc = oInc + 1     
    call setOutputValue(oInc,iceStore%ULossBot);oInc = oInc + 1
      
    USideAv = 0.0
      
    do i=1,iceStore%nCv		
        USideAv = USideAv + iceStore%ULoss(i)/iceStore%nCv
    enddo           
      
    call setOutputValue(oInc,USideAv)          ;oInc = oInc + 1
    call setOutputValue(oInc,iceStore%ULossTop);oInc = oInc + 1
       
   do i=1,iceStore%nCv
       call setOutputValue(oInc,iceStore%T(i))
       oInc = oInc + 1
   enddo
           
   call setOutputValue(oInc,iceStore%qAcumHx);oInc = oInc + 1                  
   call setOutputValue(oInc,iceStore%sumQIceAcum);oInc = oInc + 1
   
   if(iceStore%timeInHours>=12.6) then
       iceStore%dummy=1
    end if 
   
   
end subroutine setOutputs

subroutine checkParameters(iceStore,immersedHx)
      
    use spfGlobalConst
    use iceStoreDef
    use hxModule
    use TrnsysFunctions
     
    implicit none
   
    type(iceStoreStruct), intent(inout), target :: iceStore  
    type(hxStruct), intent(inout) :: immersedHx(nIHX)
    
    integer :: iUnit,iType,i,j
    character (len=nMaxCharLength) :: myMessage
         
    iUnit = iceStore%iUnit
    iType = iceStore%iType       
                    
    if(iceStore%nHx <1 .or. iceStore%nHx>nIHX) then ! number of heat exchangers 
        write(myMessage,*) 'The number of heat exchangers must be between 1 &
        and 4'
        call messages(-1,trim(myMessage),'fatal',iUnit,iType)              
    endif
        
    if(iceStore%VTank<=0.0) then
        write(myMessage,*) 'Tank volume can not be zero or negative'
        call messages(-1,trim(myMessage),'fatal',iUnit,iType)          
    
    endif    
   
    if(iceStore%HTank<=0.0) then
        write(myMessage,*) 'Tank height can not be zero or negative'
        call messages(-1,trim(myMessage),'fatal',iUnit,iType)          
       
    endif 
    
    if(iceStore%WTank<=0.0 .and. iceStore%tankGeometry==0) then
        write(myMessage,*) 'Tank width can not be zero or negative in rectangular shape'
        call messages(-1,trim(myMessage),'fatal',iUnit,iType)          
        
    endif 
    
    if(iceStore%tankGeometry /= 0 .and. iceStore%tankGeometry /= 1) then
        write(myMessage,*) 'Tank geometry must be or 0 (box) or 1 (cilinder)'
        call messages(-1,trim(myMessage),'fatal',iUnit,iType)          
        
    endif 

   ! if(iceStore%UTank<0.0) then
   !     write(myMessage,'("Tank losses (U-value) can not be negative")')
   !     call messages(-1,trim(myMessage),'fatal',iUnit,iType)          
   ! endif 
    
    if(iceStore%kEff<0.0) then
        write(myMessage,*) 'Tank effective conductivity can not be negative'
        call messages(-1,trim(myMessage),'fatal',iUnit,iType)                
    endif 
   
    
    if(iceStore%rhoWater<700.0) then
        write(myMessage,*) 'Water density is too low <700. It must be in [kg/m3]'
        call messages(-1,trim(myMessage),'warning',iUnit,iType)          
       
    endif 
   
     if(iceStore%cpWater<3800.0) then
        write(myMessage,*)'Water thermal capacity is too low <3800. It must be in [K/KgK]'
        call messages(-1,trim(myMessage),'warning',iUnit,iType)          
      
     endif 
     
     if(iceStore%rhoIce<700.0) then
        write(myMessage,*)'Ice density is too low <700. It must be in [kg/m3]'
        call messages(-1,trim(myMessage),'warning',iUnit,iType)          
        
     endif 
     
     if(iceStore%kIce>3 .or. iceStore%kIce<0.5) then
         
        write(myMessage,*) 'Ice conductivity =',iceStore%kIce,' is out of range [0.5-3].It must be in [W/mK]'               
        call messages(-1,trim(myMessage),'warning',iUnit,iType)                  
     endif 
     
    if(iceStore%H_pc<3.0e5) then        
        write(myMessage,*) 'Ice heat of fusion =',iceStore%H_pc,' is too low (<3e5).It must be in [J/kgK]'
        call messages(-1,trim(myMessage),'warning',iUnit,iType)                  
    endif 
    
    if(iceStore%TSubcool>1.0) then
        write(myMessage,*) 'Subcooling T =',iceStore%TSubcool,' is too high (>1 oC)'
        call messages(-1,trim(myMessage),'warning',iUnit,iType)                  
    endif   
        
   if(iceStore%TFreeze>5.0 .or. iceStore%TFreeze<0.) then
        write(myMessage,*) 'Freezing T =',iceStore%TFreeze,' is out of range (0 < TFreeze < 5 oC)'
        call messages(-1,trim(myMessage),'warning',iUnit,iType)                  
   endif    

   if(iceStore%iceFloatingIni<0.0 .or. (iceStore%iceFloatingIni/iceStore%rhoIce)>iceStore%VTank) then
        write(myMessage,*) 'Initial block Kg in the store iceFloatingIni =',iceStore%iceFloatingIni,'&
         is out of range (0 < iceFloatingIni/rhoIce < volume of tank)'
        call messages(-1,trim(myMessage),'fatal',iUnit,iType)                  
       
   endif
       
    if(iceStore%meltCrit<0.0) then
        write(myMessage,*) 'Critical melting distance =',iceStore%meltCrit,'&
         must be >= 0'
        call messages(-1,trim(myMessage),'fatal',iUnit,iType)                               
    endif
      
    if(iceStore%maxIceFrac<0.0 .or. iceStore%maxIceFrac>=1.0) then
        write(myMessage,*) 'Maximum ice mass fraction =',iceStore%maxIceFrac,'&
         is out of range [0-0.99]. Some water must be left in each layer !!'
        call messages(-1,trim(myMessage),'fatal',iUnit,iType)                  
           
    endif
            
    do i=1,iceStore%nHx
        if(immersedHx(i)%orderHx>4 .or. immersedHx(i)%orderHx<0) then
            write(myMessage,*) 'Wrong orderHx =',immersedHx(i)%orderHx,' (4 maximum) for Hx=',i,' if Not connected use 0'
            call messages(-1,trim(myMessage),'fatal',iUnit,iType)
        end if
               
        !Check that we dont have two times the same order Hx (it would mean that we need to divide mass flow
        ! It can be done, but its not implemented yet.
        do j=1,iceStore%nHx
            if(i .ne. j) then
                if(immersedHx(i)%orderHx==immersedHx(j)%orderHx .and. immersedHx(j)%orderHx>0) then
                    write(myMessage,*) 'Wrong orderHx =',immersedHx(i)%orderHx,' for Hx=',i,' Same index as Hx',j
                    call messages(-1,trim(myMessage),'fatal',iUnit,iType)
                end if
            end if
        end do                                                      
                     
    end do
    
    !do i=1,iceStore%nHx
    !    
    !    do j=1,iceStore%nHx
    !        
    !        if(iceStore%solutionOrder(j)==i) then                                      
    !        
    !            if(immersedHx(j)%orderHx .ne. (immersedHx(iceStore%solutionOrder(i-1))%orderHx+1) ) then
    !                write(myMessage,*) 'Wrong orderHx =',immersedHx(i)%orderHx,' for Hx=',i,' It should follow Hx',iceStore%solutionOrder(i-1)
    !                call messages(-1,trim(myMessage),'fatal',iUnit,iType)                
    !            end if 
    !        
    !            end if
    !        end if 
    !        
    !    enddo
    !    
    !enddo
    
        do i=1,iceStore%nHx                         
            
             if(immersedHx(i)%numberOfCv>0) then
                 
                 if(immersedHx(i)%geometry>2 .or. immersedHx(i)%geometry<0) then
                     write(myMessage,*) 'Wrong geometry [0,1] =',immersedHx(i)%geometry,'for Hx=',i
                     call messages(-1,trim(myMessage),'fatal',iUnit,iType)
                 end if
            
                 if(immersedHx(i)%zInDef<0.0 .or. immersedHx(i)%zInDef>iceStore%HTank) then               
                    write(myMessage,*) 'Inlet height are out of range for Hx=',i
                    call messages(-1,trim(myMessage),'fatal',iUnit,iType)   
                 end if            
             
                 if(immersedHx(i)%zOut<0.0 .or. immersedHx(i)%zOut>iceStore%HTank) then
                    write(myMessage,*) 'Outlet height are out of range for Hx=',i
                    call messages(-1,trim(myMessage),'fatal',iUnit,iType)      
                 end if
            
                 if(abs(immersedHx(i)%zOut-immersedHx(i)%zIn)<0.001) then
                    write(myMessage,*) 'Outlet height is too close to inlet height for Hx=',i,'. Use a wider range of at least 0.001'
                    call messages(-1,trim(myMessage),'fatal',iUnit,iType)
                 end if  
             
                 if(immersedHx(i)%rhoConstant<500.0) then
                     write(myMessage,*) 'Fluid density =',immersedHx(i)%rhoConstant,' of HX=',i,' is too low (<500).&
                    Units must be in [kg/m3]'                       
                    call messages(-1,trim(myMessage),'warning',iUnit,iType)                         
                 endif 
   
                 if(immersedHx(i)%cpConstant<3000.0) then
                    write(myMessage,*) 'Fluid thermal capacity =',immersedHx(i)%cpConstant,' of HX=',i,' &
                    is too low <3000. It must be in [J/KgK]'
                    call messages(-1,trim(myMessage),'warning',iUnit,iType)                         
                 endif 
     
                 if(immersedHx(i)%glycolConc<0.0 .or. immersedHx(i)%glycolConc>100.0) then
                    write(myMessage,*) 'Fraction of antifreeze =',immersedHx(i)%glycolConc,' &
                    of HX=',i,' is out of range [0-100]'  
                    call messages(-1,trim(myMessage),'fatal',iUnit,iType)          
             
                 endif    
           
                 if(immersedHx(i)%glycolConc>0.0 .and. immersedHx(i)%glycolConc<=1.0) then
                    write(myMessage,*) 'Fraction of antifreeze =',immersedHx(i)%glycolConc,' &
                    of HX=',i,' must be in [%] so from 0-100. It seems too low !'  
                    call messages(-1,trim(myMessage),'warning',iUnit,iType)                      
                 endif                                                    
                 if(immersedHx(i)%Lhx<0.0) then
                   write(myMessage,*) 'Height =',immersedHx(i)%Lhx,' &
                   of HX=',i,'can not be negative'  
                   call messages(-1,trim(myMessage),'fatal',iUnit,iType)                  
                 endif   
                 
            endif
             
        enddo                      
       
                 
        do i=1,iceStore%nCv                                                          
            if(iceStore%Uloss(i)<0.) then
                write(myMessage,*) 'ULoss calculated  =',iceStore%Uloss(i),' &
                from control volume =',i,' can not be negative'          
                call messages(-1,trim(myMessage),'fatal',iUnit,iType)          
                
            endif    
               
        enddo
                       
                  
end subroutine checkParameters

subroutine printMessagesLastTimeStep(iceStore,immersedHx)
    
    use spfGlobalConst
    use iceStoreDef
    use TrnsysFunctions
    use hxModule
    
    implicit none
   
    type(iceStoreStruct), intent(inout), target :: iceStore  
    type(hxStruct), intent(inout) :: immersedHx(nIHX)
     
    integer i
    
    if(iceStore%verboseLevel>=1) then
            
            do i=1,iceStore%nCv
            if(iceStore%T(i)<-1) then                    
                write(iceStore%MyMessage,'("NEGATIVE T=",f" in STORAGE Cv=",i)') iceStore%T(i),i
                call Messages(-1,Trim(iceStore%MyMessage),'NOTICE', iceStore%iUnit,iceStore%iType)                    
            endif
        enddo
             
        if(abs(iceStore%imbalance)>1e-3) then                         
            write(iceStore%MyMessage,'("IMB Res Balance=",f)') iceStore%imbalance
            call Messages(-1,Trim(iceStore%MyMessage),'NOTICE', iceStore%iUnit,iceStore%iType)               
                
        endif
            
        if(iceStore%errorStorageLoop>iceStore%errorMaxStorage) then
            write(iceStore%MyMessage,'("MAX NUM of storage fluid ITERATION with ERROR= =",f)') iceStore%errorStorageLoop
            call Messages(-1,Trim(iceStore%MyMessage),'NOTICE', iceStore%iUnit,iceStore%iType)
        endif
            
        if(iceStore%errorIceStorageLoop>iceStore%errorMaxIceStorage) then           
            write(iceStore%MyMessage,'("MAX NUM of global ITERATION with ERROR= =",f)') iceStore%errorIceStorageLoop  
            call Messages(-1,Trim(iceStore%MyMessage),'NOTICE', iceStore%iUnit,iceStore%iType)
        endif
            
    endif

end subroutine printMessagesLastTimeStep

subroutine readInputs(iceStore,immersedHx)
     
    use iceStoreConst
    use iceStoreDef
    use iceStoreFunc
    use hxModule
    
    use TrnsysFunctions
    use Trnsysconstants
    
    implicit none
    
    type(iceStoreStruct), intent(inout), target :: iceStore  
    type(hxStruct), intent(inout),target :: immersedHx(nIHX)
    
    integer :: i, j, iHx, iInc, revertedFlow, iOrigin, oneFlowIsReverted
    double precision :: areaTop, areaSide, uTop,uSide,areaBot,uBot 
    logical :: needToSolveHx
        
    !-----------------------------------------------------------------------
    !           Converted to DIMENSIONS from standard DECK units: 
    !-----------------------------------------------------------------------
        
    iInc=1             
    needToSolveHx = .false.
    
    !  Read inputs from used Hx
    
    oneFlowIsReverted = 0 
    
    do iHx=1,iceStore%nHx    
        
        ! used to check if this was the problem
        !if(iceStore%itTrnsys<25) then ! This means iteration problems....         
        !    immersedHx(iHx)%tFluidIn    = getInputValue(iInc) ! HX Inlet temperature [°C]                                          
        !endif
        
        immersedHx(iHx)%tFluidInIt = immersedHx(iHx)%tFluidIn
        immersedHx(iHx)%tFluidIn   = getInputValue(iInc)        
        
        if(isnan(immersedHx(iHx)%tFluidIn)) then    
            call FoundBadInput(iInc,'Fatal','Inlet Fluid temperature in the HX is Not a Number (NAN)')                         
        else if((immersedHx(iHx)%tFluidIn>0 .and. immersedHx(iHx)%tFluidInIt<0.) .or. (immersedHx(iHx)%tFluidIn<0 .and. immersedHx(iHx)%tFluidInIt>0.) ) then      
            
            if(iceStore%verboseLevel>=2 .and. iceStore%itTrnsys>=5) then
                write(iceStore%myMessage,*) 'Inlet temperature of heat exchanger = ',iHx,' is switching from positive to negative in ItTrnsys=',iceStore%itTrnsys,' This will make convergence very difficult (if itTrnsys<5 this message is not printed !!)'
                call messages(-1,trim(iceStore%myMessage),'warning',iceStore%iUnit,iceStore%iType)
            end if
        endif   
        
        immersedHx(iHx)%mDotAllHxInKgh = getInputValue(iInc+1)
        immersedHx(iHx)%mDot        = immersedHx(iHx)%mDotAllHxInKgh/(3600d0*immersedHx(iHx)%nParallelHx) ! HX Inlet massflow [kg/h] -> [kg/s]        
        immersedHx(iHx)%tFluidInRev = getInputValue(iInc+2)                                      ! HX Inlet temperature of reverted flow [°C]  
        
        iInc = iInc+3
        
        immersedHx(iHx)%zIn  = immersedHx(iHx)%zInDef
        immersedHx(iHx)%zOut = immersedHx(iHx)%zOutDef
        
        revertedFlow = .false.
        
        if(abs(immersedHx(iHx)%mDot)<zeroMassFlowInKgs) then
            
            if(iceStore%verboseLevel>=1) then
                write(iceStore%myMessage,*) 'Mass flow in HX below the limit = ',immersedHx(iHx)%mDot*3600.,' kg/h'
                call messages(-1,trim(iceStore%myMessage),'warning',iceStore%iUnit,iceStore%iType)
            endif
            
            immersedHx(iHx)%mDot = 0.0            
            revertedFlow = immersedHx(iHx)%revertedFlow
            
        else            
            needToSolveHx = .true.                        
            if(immersedHx(iHx)%mDot<0.0) then
                immersedHx(iHx)%mDot = abs(immersedHx(iHx)%mDot)
                revertedFlow = .true.
                immersedHx(iHx)%zIn  = immersedHx(iHx)%zOutDef
                immersedHx(iHx)%zOut = immersedHx(iHx)%zInDef
                immersedHx(iHx)%tFluidIn = immersedHx(iHx)%tFluidInRev
            endif       
        endif
                
        if(immersedHx(iHx)%zOut>immersedHx(iHx)%zIn) then
            immersedHx(iHx)%goUp=.true.
        else
            immersedHx(iHx)%goUp=.false.
        endif 
        
        if(revertedFlow) oneFlowIsReverted = 1            
            
        immersedHx(iHx)%revertedFlow = revertedFlow
                             
    enddo   
                 
   
    if(oneFlowIsReverted) then   
        iceStore%solutionOrder(1:nIHX) = iceStore%solutionOrderNegativeMDot(1:nIHX)
    else
        iceStore%solutionOrder(1:nIHX) = iceStore%solutionOrderPositiveMDot(1:nIHX)
    endif    
               
    
    iInc = 13
    
    iceStore%mechanicalDeIce =  getInputValue(iInc)
        
    if(iceStore%mechanicalDeIce==1) then
        iInc=iInc+1
        iInc=iInc-1
    endif
    iInc    = iInc+1
        
    do i=1,iceStore%nCv
        iceStore%Tenv(i) = getInputValue(iInc)           ! Temp. of environment        [°C]            
        iInc    = iInc+1        
    enddo
        
    iceStore%tEnvBot = getInputValue(iInc)
    iceStore%tEnvTop = getInputValue(iInc+1)
      
     if(iceStore%nCv>=3) then
        ! rewritte TEnv with TTop and TTbot        
        
        areaTop  = iceStore%areaTopSurface
        areaSide = iceStore%A(iceStore%nCv)
        uTop     = iceStore%ULossTop
        uSide    = iceStore%ULoss(iceStore%nCv)
        
        iceStore%Tenv(iceStore%nCv) = (uTop*areaTop*iceStore%tEnvTop+uSide*areaSide*iceStore%Tenv(iceStore%nCv))/iceStore%UALoss(iceStore%nCv)
        
        areaBot  = iceStore%areaBotSurface 
        areaSide = iceStore%A(1)
        uBot     = iceStore%ULossBot
        uSide    = iceStore%ULoss(1)      
        
        iceStore%Tenv(1) = (uBot*areaBot*iceStore%tEnvBot+uSide*areaSide*iceStore%Tenv(1))/iceStore%UALoss(1)
     end if
     
 end subroutine readInputs
 
subroutine readParameters(iceStore,immersedHx)

!-----------------------------------------------------------------------
!                     READ PARAMETERS (serially):
!-----------------------------------------------------------------------                                         
       
        use iceStoreConst
        use iceStoreDef 
        use iceStoreFunc
        use TrnsysFunctions
        use hxModule
        use util
    
        implicit none                
        
        type(iceStoreStruct), intent(inout), target :: iceStore           
        type(hxStruct),dimension(nIHX),intent(inout) :: immersedHx
    
        double precision :: uAddBotCv,uAddTopCv,uAddBot,uAddTop,notUsed
        integer :: pInc 
        integer :: i=0,iHx,j,count
                  
        !iceStore%printFileOut = getParameterValue(1)
        
        iceStore%nHx        = 4 ! getParameterValue(1)    !: Number of heat exchangers int[1,4]
        
        iceStore%VTank      = getParameterValue(2)        !: Tank Volume [m^3]
        iceStore%HTank      = getParameterValue(3)        !: Tank height [m]
        iceStore%WTank      = getParameterValue(4)        !: Tank width [m]
        iceStore%tankGeometry  = getParameterValue(5)     !: Tank geometry 0:box 1:cilinder
       
        iceStore%xBetweenPipes = getParameterValue(6)     !: Distance between pipes in one hx for ice constrained        
        iceStore%yBetweenPipes = getParameterValue(7)     !: Distance between pipes between hx's for ice constrained  
        
        iceStore%cpIce = 0.0 !4223*0.0 !4223.0
        
        !iceStore%addNodalLossCoef = getParameterValue(7)  !: Add. nodal loss coeff. int[0,1]
        
        iceStore%kEff        = getParameterValue(8)    !: Eff. thermal cond. tank [W/m.K]
        iceStore%kEffDefined = iceStore%kEff 
        !print *,' nHX=',iceStore%nHx," VTank=",iceStore%VTank
        
        iceStore%rhoWater   = getParameterValue(9)          !: Water density [kg/m^3]
        iceStore%CpWater    = getParameterValue(10)      !: Water specific heat [J/kg.K]
        iceStore%rhoIce     = iceStore%rhoWater    ! DC !getParameterValue(11)         !: Ice density [kg/m^3]
        iceStore%kIce       = getParameterValue(12)     !: Ice therm. cond. [W/m.K]
        
        iceStore%H_pc       = getParameterValue(13)     !: Water<->Ice enthalpy [J/kg]
        iceStore%TSubcool   = getParameterValue(14)         !: Water<->Ice temperature [°C]
        iceStore%TFreeze    = getParameterValue(15)         !: Water<->Freezing temperature [°C]
        iceStore%iceFloatingIni = getParameterValue(16)     !: initial ice kg [kg]
        iceStore%meltCrit   = getParameterValue(17)         !: Film critical melting thickness [m]
        
        if(iceStore%meltCrit>=1.) then
            iceStore%deIceIsPossible=0
        endif        
        
        iceStore%maxIceFrac = getParameterValue(18)         !: The maximum ice fraction [%]
        
        iceStore%checkParameters = getParameterValue(19)                       
        
        !CHANGE OF PARAMETERS AND ORDER
!     iceStore%arrangmentHx = getParameterValue(20)        
        iceStore%heightIceHx  = getParameterValue(21)*iceStore%HTank !is realtive to total height
        iceStore%maxIceFracIceLayer = getParameterValue(22)
        iceStore%useTwallOld = getParameterValue(23) !use TWallOld for UA values. Convergence improvement       
        
        iceStore%nRealHx = getParameterValue(24)
        iceStore%constantPhysicalProp = getParameterValue(25)
        iceStore%useCorrugatedPlate   = getParameterValue(26)
        
    !    Heat exchanger parameters (read all IHX, even unused ones):
        !do i=1,iceStore%nHx
                
        pInc = 27
 
        iceStore%nHxStore = 0
         
        do iHx=1,nIHX ! 19 inputs for each HX                         
                            
            ! general (2) 
            
            immersedHx(iHx)%reference = iHx          
            immersedHx(iHx)%geometry    = getParameterValue(pInc); pInc=pInc+1 
            immersedHx(iHx)%nParallelHx = getParameterValue(pInc); pInc=pInc+1 
            
            if(immersedHx(iHx)%isUsed==1) iceStore%nHxStore = iceStore%nHxStore+immersedHx(iHx)%nParallelHx
                       
            !for each hx type (5 + 1 lambdaWall  =6)
            
            if(immersedHx(iHx)%geometry==PLATE) then !flat plate                                
                
                immersedHx(iHx)%Hhx   = getParameterValue(pInc); pInc=pInc+1       !  HX characteristic height           [m]
                immersedHx(iHx)%Whx   = getParameterValue(pInc); pInc=pInc+1       !  Heat exchanger characteristic thickness                  [m]
                immersedHx(iHx)%LHx     = getParameterValue(pInc); pInc=pInc+1     !  Lenght of Hx           [m^2]
                immersedHx(iHx)%dxWallHx  = getParameterValue(pInc); pInc=pInc+1   !  Heat exchanger wall thickness                
                                  
            else if(immersedHx(iHx)%geometry==COIL) then ! pipes, coils
                
                immersedHx(iHx)%dIn  = getParameterValue(pInc); pInc=pInc+1   !> Inlet diamter of the pipe [m]
                immersedHx(iHx)%dOut = getParameterValue(pInc); pInc=pInc+1   !> Outlet diamter of the pipe [m]
                immersedHx(iHx)%LHx  = getParameterValue(pInc); pInc=pInc+1   !< lenght of pipe
                immersedHx(iHx)%addedCapacity = getParameterValue(pInc); pInc=pInc+1   !< J/m3                    
                immersedHx(iHx)%nTubes = immersedHx(iHx)%nParallelHx/iceStore%nRealHx !                
                
                !one not used             
                !pInc=pInc+1                       
            else
                pInc = pInc+4
            end if
            
            immersedHx(iHx)%orderHx= getParameterValue(pInc); pInc=pInc+1
            
            if(iHx>iceStore%nHx) then
                !This will never be used, but I want to have it outside limits to order the vector of solution in a proper way
                immersedHx(iHx)%orderHx=nIHX+1    
            endif
            
            immersedHx(iHx)%lambdaWall = getParameterValue(pInc); pInc=pInc+1
            
            !> height relative positions (2)
            
            immersedHx(iHx)%zInDef = getParameterValue(pInc); pInc=pInc+1     ! use relative height in the deck                       
            immersedHx(iHx)%zOutDef = getParameterValue(pInc); pInc=pInc+1   ! use relative height in the deck                        
                               
            !Convert to absolute height.
            
            immersedHx(iHx)%zInDef  = immersedHx(iHx)%zInDef*iceStore%HTank
            immersedHx(iHx)%zOutDef = immersedHx(iHx)%zOutDef*iceStore%HTank
            
            immersedHx(iHx)%zIn  = immersedHx(iHx)%zInDef
            immersedHx(iHx)%zOut = immersedHx(iHx)%zOutDef                        
                        
            !> fluid properties (3)
            
            immersedHx(iHx)%rhoConstant = getParameterValue(pInc); pInc=pInc+1  ! kg/m3               
            immersedHx(iHx)%cpConstant  = getParameterValue(pInc); pInc=pInc+1  ! kJ/kgK -> J/KgK
            immersedHx(iHx)%glycolConc  = getParameterValue(pInc); pInc=pInc+1
                                                                  
            !> numerical parameters (6) 
            
            immersedHx(iHx)%numberOfCv = getParameterValue(pInc); pInc=pInc+1 ! integer   
            
            
            
            immersedHx(iHx)%cNuHeat = getParameterValue(pInc); pInc=pInc+1
            immersedHx(iHx)%nNuHeat = getParameterValue(pInc); pInc=pInc+1            
            immersedHx(iHx)%cNuCool = getParameterValue(pInc); pInc=pInc+1
            immersedHx(iHx)%nNuCool = getParameterValue(pInc); pInc=pInc+1  
            
            immersedHx(iHx)%nEnhanceNu   = getParameterValue(pInc); pInc=pInc+1    ! Increase Nu
            
            !< Process some data    
        
            immersedHx(iHx)%zBot = min(immersedHx(iHx)%zOut,immersedHx(iHx)%zIn) !
            immersedHx(iHx)%zTop = max(immersedHx(iHx)%zOut,immersedHx(iHx)%zIn) !                                                                                 		    
			
            if(immersedHx(iHx)%zOut>immersedHx(iHx)%zIn) then
                immersedHx(iHx)%goUp=.true.
            else
                immersedHx(iHx)%goUp=.false.
            endif
       
        enddo    
              
        ! pInc = 20 + 18*4 = 92           
               
        uAddBotCv = getParameterValue(pInc) ; pInc=pInc+1
        uAddTopCv = getParameterValue(pInc) ; pInc=pInc+1
        uAddBot   = getParameterValue(pInc) ; pInc=pInc+1
        uAddTop   = getParameterValue(pInc) ; pInc=pInc+1
        
        ! Relative z !!
        
        iceStore%zSensor1 = getParameterValue(pInc)*iceStore%HTank ; pInc=pInc+1
        iceStore%zSensor2 = getParameterValue(pInc)*iceStore%HTank ; pInc=pInc+1
        iceStore%zSensor3 = getParameterValue(pInc)*iceStore%HTank ; pInc=pInc+1
        iceStore%zSensor4 = getParameterValue(pInc)*iceStore%HTank ; pInc=pInc+1
        iceStore%zSensor5 = getParameterValue(pInc)*iceStore%HTank ; pInc=pInc+1
                
        do i=1,iceStore%nCv 
            if(i>int(iceStore%nCv/2)) then                                    
                iceStore%Uloss(i) = uAddTopCv                                   
            else                
                iceStore%Uloss(i) = uAddBotCv                                   
            endif                
        enddo
        
        iceStore%UlossTop = uAddTop
        iceStore%UlossBot = uAddBot                                  
                            
        iceStore%order(1:nIHX) = immersedHx(1:nIHX)%orderHx                                     
        
        call sortIndex(nIHX,iceStore%order,iceStore%solutionOrderPositiveMDot)
        
        iceStore%solutionOrderNegativeMDot(1:nIHX) = nIHX+1
                
        do iHx=1,nIHX            
            iceStore%solutionOrderNegativeMDot(iHx) = nIHX-iceStore%solutionOrderPositiveMDot(iHx)+1            
        enddo                              
        
        iceStore%solutionOrder(1:nIHX) = iceStore%solutionOrderPositiveMDot(1:nIHX)
         
    end subroutine readParameters
   
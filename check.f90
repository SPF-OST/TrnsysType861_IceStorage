!    ! 861 !!!
 if(hxData%geometry==PLATE) then        
				   
	   ! For flat plate this is the thickness of only one side because we are using 
	   ! the total area for the calculation including both sides. 
	   ! Its like having the plate 2 times bigger but insolated in the back
   
	   iceForOneHX =  hxData%iceThickCv(iCenter)/hxData%nParallelHx
   else
	   if(hxData%dInIce(iCenter) > hxData%dOutHx(iCenter)) then   ! partial icing with a ice layer around
		   dUsed = hxData%dInIce(iCenter)
	   else ! continous icing
		   dUsed = hxData%dOutHx(iCenter)
	   endif
				   
	   !lChar =  hxData%dOutHx(iCenter)                                                                         
				  
	   if(abs(dUsed-hxData%dOutHx(iCenter))<1e-10) then
		   h_ice = 1e10
		   R_ice = 0.0
	   else                                  
		   ntubesX=2
		   ntubesY=2                             
		   call getOverlappingAngle(iceStore%xBetweenPipes,iceStore%yBetweenPipes,ntubesX,ntubesY,0.5*dUsed,phiOverlap)
		   h_ice =  (2.0*pi-phiOverlap)*iceStore%kIce/(pi*dUsed*log(dUsed/hxData%dOutHx(iCenter))) !W/m2
		   R_ice = 1./(h_ice)
		   iceStore%fConstrained = (2.0*pi-phiOverlap)/2.0/pi                            
	   endif
   
   end if
   
			 
   !MELTING 
   if(iceForOneHx<1e-15 .or. hxData%t(i)>tStoreAv) then                                               
	   if(iceStore%storeFull) then
		   !If storage is full we are to good becasue we neglect the water layer that is formed after a small melting     
		   !so I consider 1/2 of the iceThickness as an average number
		   
		   alphaOut = iceStore%kIce/(iceForOneHx*0.5)   ! W/m2K                                            
	   else
		   !DCAR 03.05.2013 if melting, then the ice layer does not play any role!!!
		   alphaOut = 1e20 ! W/m2K                                                               
	   endif 
   ! ICING
   else
	   if(hxData%geometry==PLATE) then                            
		   alphaOut = iceStore%kIce/iceForOneHx     ! W/m2K                                          
	   else !cilinder 
			
		   !calculateHeatTransferWall(lambdaWall,dIn,dOut,hWall,geometry)
		   !iceStore%dOutIce(i) = thick_hx+2*iceStore%ice(i)
		   !call calculateHeatTransferWall(iceStore%kIce,thick_hx,iceStore%dOutIce(i),h_o,iceStore%hxGeometry(i),iceStore%iUnit,iceStore%iType)
		   !Ro = 1./h_o
	   endif
			   
   endif
else
   call calculateExternalHeatTransferCoef(iceStore,hxData,tFilm,tStoreAv,tWall,alphaOut)

end if
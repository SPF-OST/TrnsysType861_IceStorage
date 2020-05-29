!---------------------------------------------------------------------------  
!> @file interpolation.f90
!> @author SPF institut fur SolarTechnik 
!> @author D.Carbonell
!> @date 10/10/2012
!> @brief interpolation functions
!> @todo quadratic and liner vector interpolation
!--------------------------------------------------------------------------  



    
    module interpolation



    !use spfGlobalConst
    
    !use ! other modules
    implicit none
    
!#define LIN_I(x1,y1,x2,y2,x)      ((y1)+(((x)-(x1))*((y2)-(y1))/((x2)-(x1)))) 

    !Specification. Writte here global parameters that can be used inside the module
    
    contains
    
  
    
    double precision function linear(x1,y1,x2,y2,x) 
    
       double precision, intent(in) :: x1,y1,x2,y2,x
              
       linear = y1+((x-x1)*(y2-y1)/(x2-x1))
    
    end function linear 
    
    double precision function quadratic(x1,y1,x2,y2,x3,y3,x) 
    
       double precision, intent(in) :: x1,y1,x2,y2,x3,y3,x
       
       quadratic = ((y1)*((x)-(x2))*((x)-(x3))/(((x1)-(x2))*((x1)-(x3)))+&
            (y2)*((x)-(x1))*((x)-(x3))/(((x2)-(x1))*((x2)-(x3)))+&
            (y3)*((x)-(x1))*((x)-(x2))/(((x3)-(x1))*((x3)-(x2))))
    
    end function quadratic
    
    subroutine findIndex(sizeArray,array,xIn,indexUp,indexDown)    
    
     integer, intent(in) :: sizeArray
     integer, intent(out):: indexDown, indexUp
     integer :: getOut, iinf,isup,im
     double precision, intent(in) :: array(1:sizeArray) 
     double precision :: xIn
     
     iinf=1;
     isup=sizeArray;
     getOut = 0
     
     if(sizeArray==1) then
        indexUp   = 1;
        indexDown = 1;
     else if(array(1)>xIn) then
        indexUp   = 2;
        indexDown = 1;
     else if(array(sizeArray)<xIn) then
        indexUp   = sizeArray;
        indexDown = sizeArray-1;         
     else
     
         !im = int((isup-iinf)*0.5)  
         ! if only 1 or two Cv we must control im to void 0 
         im = max(int((isup-iinf)*0.5),iinf)
         im = min(im,isup)

         do while (getOut==0)
     
              if(array(im)>xIn) then
                isup=im
              else                
                iinf=im
              endif     
	 
              im=int(iinf+(isup-iinf)*0.5);

              if((isup-iinf)==1) getOut=1;

              !print *,'isup=',isup,' iinf=',iinf,' xIn=',xIn,' ArrayUp=',array(isup),' ArrayDown=',array(iinf)
          
          enddo

          indexUp   = isup;
          indexDown = iinf;
          
     endif
     
    end subroutine findIndex

  
    double precision function vectorLinearInterpolation(sizeArray,x,y,xActual,info)

        double precision, intent(in) :: xActual
        integer, intent(in) :: sizeArray, info
        double precision  :: x(1:sizeArray),y(1:sizeArray)
        integer :: indexUP,indexDown        	    
            
		call findIndex(sizeArray,x,xActual,indexUp,indexDown)
		   
		vectorLinearInterpolation = linear(x(indexDown),y(indexDown),x(indexUp),y(indexUp),xActual);            
		
		!if(info) then
		!	write(unitDebug,*) ' vectorLinearInterpolation: valueInterpolated',vectorLinearInterpolation,' yUp( ',indexUp,')= ',y(indexUp),' yDown( ',indexDown,')= ',&
  !          y(indexDown),' xDown ',x(indexDown),' xUp ',x(indexUp),' x=',xActual
  !      endif            		                
		
    end function vectorLinearInterpolation
    
    double precision function vectorLinearInterpolationLimited(sizeArray,x,y,xActual,info)

        double precision, intent(in) :: xActual
        integer, intent(in) :: sizeArray, info
        double precision :: x(1:sizeArray),y(1:sizeArray)
        integer :: indexUP,indexDown        
	
	    if(xActual>x(sizeArray)) then
            !if(info) then
            !    write(unitDebug,*),' vectorLinearInterpolation: x=',xActual,' > x(sizeArray)=',x(sizeArray),' y(sizeArray)= ',y(sizeArray)
            !endif

            vectorLinearInterpolationLimited = y(sizeArray);
        
	    else if (xActual<x(1)) then
            if(info) then
                print*,' vectorLinearInterpolation: x=',xActual,' < x(1)=',x(1),' y(1)= ',y(1)
            endif

            vectorLinearInterpolationLimited = y(1)            
        
        else 
            
		    call findIndex(sizeArray,x,xActual,indexUp,indexDown)
		   

		    vectorLinearInterpolationLimited = linear(x(indexDown),y(indexDown),x(indexUp),y(indexUp),xActual);
            
		
		    !if(info) then
			   ! write(unitDebug,*)' vectorLinearInterpolationLimited: valueInterpolated',vectorLinearInterpolationLimited,' yUp( ',indexUp,')= ',y(indexUp),' yDown( ',indexDown,')= ',&
      !          y(indexDown),' xDown ',x(indexDown),' xUp ',x(indexUp),' x=',xActual
      !      endif            		    
            
	    endif
		
    end function vectorLinearInterpolationLimited
 
    !be carefull with index. Used for Trnsys Type708
    double precision function matrixLinearInterpolation(sizeArrayX,sizeArrayY,nodx,nody,var,xActual,yActual)
    
        double precision, intent(in) :: xActual,yActual
        integer, intent(in) :: sizeArrayX,sizeArrayY
        double precision ,intent(in) :: nodx(1:sizeArrayX),nody(1:sizeArrayY),var(0:sizeArrayX-1,0:sizeArrayY-1)      
        integer :: jIndexUp,jIndexDown,iIndexUp,iIndexDown
        double precision:: varDown,varUp         
        integer :: index 
        
        index = 0
        
        call findIndex(sizeArrayY,nody,yActual,jIndexUp,jIndexDown)
        call findIndex(sizeArrayX,nodx,xActual,iIndexUp,iIndexDown)
             
                
        varDown=linear(nodx(iIndexDown),var(iIndexDown+index,jIndexDown+index),nodx(iIndexUp),var(iIndexUp+index,jIndexDown+index),xActual)
        
        varUp=linear(nodx(iIndexDown),var(iIndexDown+index,jIndexUp+index),nodx(iIndexUp),var(iIndexUp+index,jIndexUp+index),xActual)
  
        matrixLinearInterpolation=linear(nody(jIndexDown),varDown,nody(jIndexUp),varUp,yActual)
    
        if (abs(matrixLinearInterpolation)<1e-15) then
           matrixLinearInterpolation=0.
        endif            
        
    end function matrixLinearInterpolation
    
     double precision function matrixLinearInterpolationNew(sizeArrayX,sizeArrayY,nodx,nody,var,xActual,yActual)
    
        double precision, intent(in) :: xActual,yActual
        integer, intent(in) :: sizeArrayX,sizeArrayY
        double precision ,intent(in) :: nodx(1:sizeArrayX),nody(1:sizeArrayY),var(1:sizeArrayX,1:sizeArrayY)      
        integer :: jIndexUp,jIndexDown,iIndexUp,iIndexDown
        double precision:: varDown,varUp         
       
        call findIndex(sizeArrayY,nody,yActual,jIndexUp,jIndexDown)
        call findIndex(sizeArrayX,nodx,xActual,iIndexUp,iIndexDown)
                             
        varDown=linear(nodx(iIndexDown),var(iIndexDown,jIndexDown),nodx(iIndexUp),var(iIndexUp,jIndexDown),xActual)        
        varUp=linear(nodx(iIndexDown),var(iIndexDown,jIndexUp),nodx(iIndexUp),var(iIndexUp,jIndexUp),xActual)
  
        matrixLinearInterpolationNew=linear(nody(jIndexDown),varDown,nody(jIndexUp),varUp,yActual)
    
        if (abs(matrixLinearInterpolationNew)<1e-15) then
           matrixLinearInterpolationNew=0.
        endif            
        
    end function matrixLinearInterpolationNew
       

end module interpolation

!---------------------------------------------------------------------------  
!> @file utils.f90
!> @author SPF institut fur SolarTechnik 
!> @author D.Carbonell, W.Lodge
!> @date 10/08/2012
!> @brief utility functions
!> @todo in the future we should have mathematical utility functions 
!> and miscellaneous utilitys
!--------------------------------------------------------------------------  
    module util

    !use ! other modules
    implicit none
    
    !Specification. Writte here global parameters that can be used inside the module
    
    contains


    
    double precision function evalcratl(order, logx, logy, x, y, c)
    
!     Evaluate a Chebyshev X,Y Rational Order 5/6
        integer :: order,logx,logy
        double precision :: x,y,c(*),z
        integer :: tcnt,j,m
        double precision :: tx(12),ty(12),num,den
        
        if(logx.NE.1) then
            x=(x-(50.00000000000000d0))/(50.00000000000000d0)
        else
            x=(dlog(x)-(0.000000000000000d0))/(0.000000000000000d0)
        end if
        if(logy.NE.1) then
            y=(y-(40.00000000000000d0))/(80.00000000000000d0)
        else
            y=(dlog(y)-(0.000000000000000d0))/(0.000000000000000d0)
        end if
        
        select case (order)
            case (6)
                tcnt=3
            case (10)
                tcnt=4
            case (14)
                tcnt=5
            case (18)
                tcnt=6
            case (22)
                tcnt=7
            case (26)
                tcnt=8
            case (30)
                tcnt=9
            case (34)
                tcnt=10
            case (38)
                tcnt=11
            case (42)
                tcnt=12
            case default
                z=0.0
            return
        end select
            
        if(tcnt.GT.7) then
            if(x.LT.-1.d0) x=-1.d0
            if(x.GT. 1.d0) x= 1.d0
            if(y.LT.-1.d0) y=-1.d0
            if(y.GT. 1.d0) y= 1.d0
        end if
        
        tx(1)=1.d0
        ty(1)=1.d0
        tx(2)=x
        ty(2)=y
        
        do j=3,tcnt
            tx(j)=2*x*tx(j-1)-tx(j-2)
            ty(j)=2*y*ty(j-1)-ty(j-2)
        enddo
        
        m=2
        num=c(1)
        den=1.0+c(2)*tx(m)+c(3)*ty(m)
        do j=4,order,4
            num=num+c(j)*tx(m)
            num=num+c(j+1)*ty(m)
            m=m+1
            den=den+c(j+2)*tx(m)
            den=den+c(j+3)*ty(m)
        enddo
        
        if(den.EQ.0.0) then
            evalcratl=0.0
        else
            evalcratl=(num/den)*(0.1065000000000000d0)+(1.049500000000000d0)
        end if
        
        !return
    end function evalcratl
         
  double precision function getMaxError(T,Tnew,nCv)
        integer :: nCv, i
        double precision :: T(nCv),Tnew(nCv),maxError,error
        
        maxError = 0.
        error = 0.
        
        do i=1,nCv        
            error  = error+ abs(T(i)-Tnew(i))            
            maxError  = max(abs(T(i)-Tnew(i)),maxError) 
            !write(20,*) 'i=',i,' maxError ',maxError,' error ',abs(T(i)-Tnew(i))           
        enddo
        
        error = error/maxval(abs(T(1:nCv)))
        
        getMaxError = maxError
        !getMaxError = error
        
        !if (sqrt(dot_product(errT(1:n), errT(1:n))) < (tol * maxval(abs(T(1:n)))))                 
        
  end function getMaxError
  
  double precision function getLMTD(tIn, tOut, tSource)
   
        double precision, intent(in) :: tIn, tOut, tSource
        double precision :: dTi, dTo
        
        dTi = (tIn-tSource)
        dTo = (tOut-tSource)
        
        if(dTi==dTo) then
            getLMTD = 0.0d0
        else if(dTi<1e-15 .or. dTo<1e-15 .or. (dTi*dTo)<1e-15 ) then
            getLMTD = 0.0d0
        else
            getLMTD = (dTi-dTo)/(log(dTi/dTo))            
        endif    
            
  end function getLMTD
  
  subroutine sortIndex(n,vector,indexSorted)
   
    integer, intent(in) :: n
    integer, intent(in) :: vector(n)
    integer, intent(out) :: indexSorted(n)     
    integer :: min,i,j,jj,notUsedYet
    
    indexSorted(1:n)=-1     
   
    
    do i=1,n ! loop for remaining indexSorted         
        min=n+1                    
        do j=1,n  ! loop for all n vectors        
            notUsedYet=1
            do jj=1,n
                if(j==indexSorted(jj)) notUsedYet=0
            enddo
            
            if(vector(j)<min .and. notUsedYet==1) then                                    
                min=vector(j)
                indexSorted(i)=j
            endif    
            
        enddo
    enddo
    
  end subroutine sortIndex
  
end module util

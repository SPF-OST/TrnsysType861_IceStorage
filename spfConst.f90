!---------------------------------------------------------------------------  
!> @file spfConst.f90
!> @author SPF institut fur SolarTechnik 
!> @author D.Carbonell
!> @date 10/08/2012
!> @brief constant values
!> @todo 
!--------------------------------------------------------------------------  

module spfGlobalConst

    implicit none                

    double precision :: PI = 3.14109265358979323
    double precision :: g = 9.81       
    double precision :: sigmaBolzman = 5.670373d-8 !W/(m2*K4)
    integer , parameter :: nMaxCharLength = 200
    !integer :: myScreenUnit = 20
      
end module spfGlobalConst

module spfAlgorithmConst

    integer :: CONTINUE_ITE = 0, &
               CONVERGED = 1, &                            
               MAX_NITER_REACHED = 2, &
               DIVERGED = 3  
    
end module spfAlgorithmConst


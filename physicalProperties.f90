!---------------------------------------------------------------------------  
!> @file physicalProperties.f90
!> @author SPF institut fur SolarTechnik 
!> @author W.Lodge, D.Carbonell
!> @date 10/08/2012
!> @brief physical properties for water and antifreeze solutions
!> @todo this functions should be used in an structure of physical data were all programs
!> can access to it.
!--------------------------------------------------------------------------  
module physProp

    use util
    
    implicit none
    
    !Specification. Writte here global parameters that can be used inside the module
    
    contains

!------------------------------------------------------------------------
    double precision function getAlphaDestratify(T)
        double precision :: T
        if (T.gt.4.0d0) then
            !       "relatively" linear up until the boiling temperature
            getAlphaDestratify = 3.729d-10*T + 1.348d-7
        else
            !       dramatically destratifying (using high alpha -> 0°C)
            getAlphaDestratify = 5.642d-6*T**2 - 4.753d-5*T + 1d-4
        endif
        ! return in [m^2/s]
    end function getAlphaDestratify
!!-----------------------------------------------------------------------
    double precision function getBetaWater(T)
        double precision :: T
        if (T.gt.8.0d0) then
            !"relatively" linear up until the boiling temperature
            getBetaWater = ((0.8d0 * T**0.5348d0) - 1.9114d0) * 1e-4
        else if (T < 1.0d0) then
            !write(20,*) 'WARNING getBetaWater : Water at freezing temperature T',T
            getBetaWater = 5.22d-5
        else        
            !parabolic around 4°C
            getBetaWater = 3.257d-6*T**2 - 2.606d-5*T + 5.22d-5
        endif
    ! return                                                     ![K^-1]
    end function getBetaWater
!!-----------------------------------------------------------------------
    double precision function getCpWater(T)
        double precision :: T
        if (T <= 0.5d0) then
            !write(unitDebug,*) 'WARNING getCpWater : Water at freezing temperature T',T        
            getCpWater = 4209.1d0
        else 
            getCpWater = 4209.1d0 - (132.8e-2 * T) + (143.2e-4 * T**2)
        endif
        !return                                                 ![J/(kg.K)]
    end function getCpWater
!!-----------------------------------------------------------------------
    double precision function getLambdaWater(T)
        double precision :: T
        if(T <= 0.5d0) then
            getLambdaWater = 0.520d0 
        else
            getLambdaWater = 0.520d0 + (0.0198d0 * T**0.46d0)
        endif
        !return                                                  ![W/(m.K)]
    end function getLambdaWater

!!-----------------------------------------------------------------------
    double precision function getNuWater(T)
        double precision :: T
        if (T <= 0.5d0) then
            !write(20,*) 'WARNING getNuWater : Water at freezing temperature T',T            
            getNuWater = 1.477e-6
        else    
            getNuWater = 1.477e-6 * exp(-1.747e-2 * T)
        endif    
        !return                                                    ![m^2/s]
    end function getNuWater
    !
!!-----------------------------------------------------------------------
    double precision function getRhoWater(T)
        double precision :: T
        if (T > 8.0d0) then
            getRhoWater = 1000d0 - (0.0128d0 * T**1.76d0)
        else if (T < 1.0d0) then
            !write(20,*) 'WARNING getRhoWater : Water at freezing temperature T',T
            getRhoWater = 999.5d0
        else
    !        parabolic around 4°C
             getRhoWater = -0.01546d0*T**2 + 0.1237d0*T + 999.5d0
        endif
    !      return                                                   ![kg/m^3]
    end function getRhoWater
!      
!!-----------------------------------------------------------------------
    double precision function getPrandtlWater(T)
        double precision :: T
        getPrandtlWater = (39.5345d0 * T**-0.144d0) - 18.8396d0
        !return                                                        ![-]
    end function getPrandtlWater
!
!!-----------------------------------------------------------------------
    double precision function getCpEg(x,y)
!     returns specific heat of Clariant Antifrogen N
!     x : Concentration [v/v]
!     y : Temperature   [°C]
        double precision x,y,z
        double precision z1,z2,z3,z4
        
        if(x<0.001) then 
            getCpEg = getCpWater(y)
        else    
            z1=4.187215724620146d0+x*(-0.1284734006185055d0+x*(0.001721285576589111d0+x*(-7.254731103592794d-06)))
            z2=y*(0.005560804328893099d0+y*(5.708359419409403d-06))
            z3=1.000000000000000d0+x*(-0.02809958303149450d0+x*(0.0003711762756828782d0+x*(-1.336458939136372d-06)))
            z4=y*(0.001368223833593263d0)
            z=(z1+z2)/(z3+z4)
            getCpEg = z * 1e3 ! convert kJ to J
        endif    
        !return                                                 ![J/(kg.K)]
    end function getCpEg
!
!!-----------------------------------------------------------------------
    double precision function getLambdaEg(x,y)
    !     returns thermal conductivity of Clariant Antifrogen N
    !     x : Concentration [v/v]
    !     y : Temperature   [°C]
        double precision x,y,z
        double precision z1,z2,z3,z4
        
        if(x<0.001) then 
            getLambdaEg = getLambdaWater(y)
        else               
            z1=0.5788960571804519d0+x*(-0.01193439180735611d0+&
            x*(0.0001645555229344245d0+x*(-6.445468236604122d-07)))
            z2=y*(-0.003159998140240872d0+&
            y*(1.264728802532987d-05))
            z3=1.000000000000000d0+x*(-0.01306568812884102d0+&
            x*(0.0001594368726120601d0))
            z4=y*(-0.007228670055182665d0+&
            y*(2.924554017266508d-05))
            z=(z1+z2)/(z3+z4)
            getLambdaEg = z
        endif
        
        ! return                                                  ![W/(m.K)]
    end function getLambdaEg

!!-----------------------------------------------------------------------
    double precision function getNuEg(x,y)
    !     returns kinematic viscosity of Clariant Antifrogen N
    !     x : Concentration [v/v]
    !     y : Temperature   [°C]
        double precision, intent(in) :: x,y
        double precision :: z= 0.d0
        
        if(x<0.001) then 
            getNuEg = getNuWater(y)
        else   
            z=(1.976480821894529d0+&
            x*(0.02058411691044437d0+x*(5.080655484562748d-05))+&
            y*(-0.01762174452170111d0+y*(0.0001174060019474829d0))+&
            x*y*(-1.106550029261898d-06))/&
            (1+x*(-0.01694716459943308d0+x*(7.468375992873698d-05))+&
            y*(0.02535986151154673d0+y*(0.0001346786752652095d0))+&
            x*y*(-0.0002352776033071207d0))
            getNuEg = z * 1e-6 ! convert mm^2 to m^2
        endif
        
        !return                                                    ![m^2/s]
    end function getNuEg
!      
!!-----------------------------------------------------------------------

    double precision function getRhoEg(x,y)
    !     returns density of Clariant Antifrogen N
    !     x : Concentration [v/v]
    !     y : Temperature   [°C]
        double precision, intent(in) :: x,y
        double precision :: z,x1,y1
        double precision :: c(22+1)
        
        if(x<0.001) then 
            getRhoEg = getRhoWater(y)
        else   
            data c(1)/0.01953055969786295d0/
            data c(2)/-0.1285671415225388d0/
            data c(3)/0.1294231581454231d0/
            data c(4)/0.6816495959420304d0/
            data c(5)/-0.3575538097164653d0/
            data c(6)/0.04965691136446271d0/
            data c(7)/-0.08477006297533489d0/
            data c(8)/-0.08677012028585240d0/
            data c(9)/-0.08706374513717987d0/
            data c(10)/0.02510117459658974d0/
            data c(11)/-0.008050075207446592d0/
            data c(12)/0.009451106482821699d0/
            data c(13)/0.02001072910627018d0/
            data c(14)/-0.01580223956530836d0/
            data c(15)/0.002123046330383432d0/
            data c(16)/0.01225430487063658d0/
            data c(17)/0.005153019143803391d0/
            data c(18)/-0.001863831885494377d0/
            data c(19)/-0.003652859624057242d0/
            data c(20)/-0.003561961931131665d0/
            data c(21)/-3.396406970700756d-05/
            data c(22)/-0.008822474424373947d0/
            data c(23)/-0.003196774313307572d0/
            x1 = x
            y1 = y
            z = evalcratl(22,0,0,x1,y1,c)
            getRhoEg = z * 1e3 ! convert litres to m^3
            !return                                                   ![kg/m^3]
        endif
    end function getRhoEg
      
end module physProp
    

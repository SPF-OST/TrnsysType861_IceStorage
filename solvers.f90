!---------------------------------------------------------------------------  
!> @file solvers.f90
!> @author SPF institut fur SolarTechnik 
!> @author W.Lodge, D.Carbonell
!> @date 10/08/2012
!> @brief functions to solve systems of equations
!--------------------------------------------------------------------------  
module solvers

contains

!-----------------------------------------------
!>@brief SOLVE TRI-DIAGONAL MATRIX | A x = d |\n
!> b*t(i) + a*t(i-1) + c*(i+1) = d \n
!> A : \n
!>       | b(1) c(1)  0    0    ...    0    |\n
!>       | a(2) b(2) c(3)  0    ...    0    |\n
!>       |  0   a(2) b(3) c(4)  ...    0    |\n
!>       |  :    :    :    :     :     :    |\n
!>       |  0    0    0   a(n)  b(n)  c(n)  |\n
!>       |  0    0    0    0   a(n+1)b(n+1) |\n
!>x : \n
!>       |  U(1)  |\n
!>       |   :    |\n
!>       | U(n+1) |\n
!> d: \n
!>       |  d(1)  |\n
!>       |   :    |\n
!>       | d(n+1) |\n
!>@param n : Number of control volumes
!>@param a : Sub-diagonal of A
!>@param b : Main-diagonal of A
!>@param c : Super-diagonal of A
!>@param d : Source term 
!>@param U : The solution to the independent variable
!>@return U : The solution to the independent variable
!---------------------------------------------------

     subroutine solveTDMA(n,a,b,c,d,U)
     
      implicit none
      integer,intent(in) :: n
      double precision,dimension(n),intent(in) :: a,b,c,d
      double precision,dimension(n),intent(out) :: U
     
      double precision,dimension(n) :: bp,dp
      double precision :: m, res
      integer i
    
      
!     Make copies of b and v (keep originals!)
      bp(:) = b(:)
      dp(:) = d(:)
 
!     The first pass (setting the coefficients):
      do i = 2,n
        m = a(i) / bp(i-1)
        bp(i) = b(i) - m * c(i-1)
        dp(i) = d(i) - m * dp(i-1)
      enddo
      U(n) = dp(n) / bp(n)
      
!     The second pass (backwards-substition)
      do i = n-1, 1, -1
        U(i) = (dp(i) - c(i) * U(i+1)) / bp(i)                     
      enddo
 
      !> Check that it solves correctly
      !do i=1,n
      !  if(i==1) then
      !      res = b(i)*U(i)+c(i)*U(i+1)-d(i)
      !  elseif (i==n) then
      !      res = b(n)*U(n)+a(n)*U(n-1)-d(n)              
      !  else
      !  res = b(i)*U(i)+a(i)*U(i-1)+c(i)*U(i+1)-d(i)           
      !  endif
      !    
      !  print *,' TDMA i ',i,' RES= ',res,' T(i)=',U(i)
      !  
      !  if(abs(res)>1e-5) then
      !      call exit(0)
      !  endif                               
      !  
      !enddo
      
     
      
    end subroutine solveTDMA
    
    subroutine solveBDTBMA(n,a,b,c,d,U)
     
      implicit none
      integer,intent(in) :: n
      double precision,dimension(n,2,2),intent(in) :: a,b,c
      double precision, dimension(n,2),intent(in) :: d
      double precision,dimension(100,2),intent(inout) :: U
     
      double precision,dimension(n,2,2) :: bp  , a_inv
      double precision,dimension(n,2) :: Y,dp
      double precision :: m, res
      integer i
    
      
!     Make copies of b and v (keep originals!)
      bp(:,:,:) = b(:,:,:)
      dp(:,:) = d(:,:)
      
!     calculate inverse
      do i = 1,n
         a_inv(i,1,1)=a(i,2,2)/(-a(i,2,1)*a(i,1,2)+a(i,1,1)*a(i,2,2))
         a_inv(i,1,2)=-a(i,1,2)/(-a(i,2,1)*a(i,1,2)+a(i,1,1)*a(i,2,2))
         a_inv(i,2,1)=-a(i,2,1)/(-a(i,2,1)*a(i,1,2)+a(i,1,1)*a(i,2,2))
         a_inv(i,2,2)=a(i,1,1)/(-a(i,2,1)*a(i,1,2)+a(i,1,1)*a(i,2,2))
      enddo
    
!     i=1
      Y(1,:) = matmul(a_inv(1,:,:),d(1,:))
!     The first pass (setting the coefficients):
      U(1,:)=Y(1,:)
      do i = 2,n
        Y(i,:) = matmul(a_inv(i,:,:),d(i,:)-matmul(b(i,:,:),Y(i-1,:)))
        U(i,:) = Y(i,:)
      enddo
      
      
!     
 
      !> Check that it solves correctly
      !do i=1,n
      !  if(i==1) then
      !      res = b(i)*U(i)+c(i)*U(i+1)-d(i)
      !  elseif (i==n) then
      !      res = b(n)*U(n)+a(n)*U(n-1)-d(n)              
      !  else
      !  res = b(i)*U(i)+a(i)*U(i-1)+c(i)*U(i+1)-d(i)           
      !  endif
      !    
      !  print *,' TDMA i ',i,' RES= ',res,' T(i)=',U(i)
      !  
      !  if(abs(res)>1e-5) then
      !      call exit(0)
      !  endif                               
      !  
      !enddo
      
     
      
    end subroutine solveBDTBMA
    
    subroutine solveTDMA2D(fi,ap,ae,aw,an,as,bb,TX,TY,errorMax,relax,myError,nIte)
      
        implicit none
                     
        double precision, intent(in) :: errorMax,relax
        integer, intent(in) :: TX,TY
        integer, intent(out) :: nIte
        double precision, intent(out) :: myError
        integer:: i,j,horizontalPass,verticalPass,nIteMax
        double precision,dimension(TX) :: aex,awx,apx,fix,bx
        double precision,dimension(TY) :: any,asy,apy,fiy,by
        double precision,dimension(TX,TY) :: fi,fiOld,ap,ae,aw,an,as,bb
       
        horizontalPass = 1
        verticalPass = 0
        myError = 1e10
        nIte = 0
        nIteMax = 100
        
        do while(myError>errorMax .and. nIte<nIteMax)
                        
            fiOld(:,:) = fi(:,:)

            nIte = nIte+1
                                
            ! horizontal line  
            if(horizontalPass==1) then  
                
                do j=1,TY 
                    do i=1,TX                                             
                        if(j==1) then
                            bx(i) = bb(i,j) - an(i,j)*fi(i,j+1) 
                        else if(j==TY) then
                            bx(i) = bb(i,j) - as(i,j)*fi(i,j-1)
                        else                
                            bx(i) = bb(i,j) - an(i,j)*fi(i,j+1) - as(i,j)*fi(i,j-1)
                        endif
                    enddo    
                    aex = ae(:,j)
                    apx = ap(:,j)
                    awx = aw(:,j)                    
                                                                                                    
                    call solveTDMA(TX,awx,apx,aex,bx,fix)
  
                    do i=1,TX 
                        fi(i,j) = fix(i)
                    enddo                    
               enddo
            endif
              
            ! vertical line
            if(verticalPass==1) then
                do i=1,TX                
                    do j=1,TY            
                        if(i==1) then
                            by(j) = bb(i,j) - ae(i,j)*fi(i+1,j) 
                        else if(i==TX) then
                            by(j) = bb(i,j) - aw(i,j)*fi(i-1,j)
                        else                
                            by(j) = bb(i,j) - ae(i,j)*fi(i+1,j) - aw(i,j)*fi(i-1,j)
                        endif
                    enddo
                    
                    any = an(i,:) 
                    apy = ap(i,:)
                    asy = as(i,:)                   
                    
                    call solveTDMA(TY,asy,apy,any,by,fiy)
                                        
                    
                    do j=1,TY 
                        fi(i,j) = fiy(j)
                    enddo   
                   
                enddo
            endif

            myError = 0.0d0                       
   
             do j=1,TY 
                do i=1,TX  
                    myError = max(myError,abs(fi(i,j)-fiOld(i,j)))                    
                    fi(i,j) = fiOld(i,j)*(1.-relax)+ fi(i,j)*relax
                enddo
             enddo
             
    enddo
        
   
    
 end subroutine solveTDMA2D   
    
    
!---------------------------------------------------
!>@brief     1D HEAT EQUATION SOLVED WITH CRANK-NICOLSON SCHEME (NEUMANN BCs)
!>@param  n     : Number of control volumes
!>@param  dt    : Time-step                                        [seconds]
!>@param  dx    : Constant mesh discretizer                              [m]
!>@param  alpha : Diffusion coefficient                              [m^2/s]
!>@param  S     : The Sink/Source term                                 [K/s]
!>@param  T     : Temperature vector for each control volume
!>@param relax  : relaxation factor
!>@return  T    : Temperature (dependent variable at intervals dx)      [°C]
!---------------------------------------------------

    subroutine marchTime(n, dt, dx, alpha, T, S,relax)                   

      implicit none
      
      integer,intent(in) :: n, dt
      double precision, intent(in) :: dx, alpha, relax
      double precision, dimension(n), intent(inout) :: T
      double precision, dimension(n), intent(inout) :: S
      double precision, dimension(n) :: a,b,c,d
      double precision :: r, Told(n)      
      integer :: i

!     Compute Fourier number
      r = (alpha * dt) / (dx**2+1e-20)
      
      !print *, 'r alpha dt dx ', r, alpha, dt,dx
      
!     Create Tri-diagonal matrix | A x = d |
!     Initialise coeffs. for the left-hand-side of the PDE system (A)
      a(1:n) = -r
      c(1:n) = -r
      b(1:n) = (2.0d0 + (2.0d0 * r))
      d(1:n) = 0.0d0
    
      !do i=1,n
          !print *, 'MarchTime a',i,a(i)
      !enddo
      
!     Neumann boundary conditions
      b(1) = 2.0d0 + r
      b(n) = 2.0d0 + r

!     The right-hand-side (d) includes time derivative (in r)
!     and the crank-nicolson terms from the previous time-step:
      d(2:n-1) = (r * T(1:n-2)) + ((2.0d0 - 2.0d0*r) * T(2:n-1)) + (r * T(3:n))
      
!     Enforce Neumann boundary conditions at mesh ends
      d(1) = ((2.0d0 - r) * T(1)) + (r * T(2))
      d(n) = (r * T(n-1)) + ((2.0d0 - r) * T(n))
      
!     Add sink/source term to d (twice; for old and new times-steps):
!     (In essence the power densities are created outside marchTime,
!     whereas in the documentation they are shown inside)
!     Units of S -> [K/s]
      d(1:n) = d(1:n) + (2d0 * dt * S(1:n))
      
!     Find solution with Tri-Diagonal Matric Algorithm
      Told(1:n) = T(1:n)
      
      call solveTDMA(n,a,b,c,d,T)            
      
      T(1:n) = relax*T(1:n)+(1-relax)*Told(1:n)
      
    end subroutine marchTime
    
    subroutine gaussPivot(n, a, b, x)
    
!     GAUSSIAN ELIMINATION WITH SCALED ROW PIVOTING

      implicit none
      integer :: n, i, j, k, z, p(n), pk
      double precision :: sum, smax, r, rmax
      double precision :: a(n,n), b(n), s(n), x(n)
      
      b(:) = -b(:)
      do i=1,n
        p(i) = i
        smax = 0.0 
        do j=1,n
          smax = max(smax,abs(a(i,j)))
        enddo
        s(i) = smax
      enddo

!     Row interchange, if needed
      do k=1,n-1
        rmax = 0.0
        do i=k,n
          r = abs(a(p(i),k))/s(p(i))
          if (r.gt.rmax) then
            j = i
            rmax = r
          endif
        enddo

!       For some reason, j sometimes appears outside the limits
!       of 1 and n, hence this hack job:
        if (j.gt.n) then
          j = n
        endif
        
        pk = p(j)
        p(j) = p(k)
        p(k) = pk

!       Elimination
        do i=k+1,n      
          z = a(p(i),k)/a(p(k),k)       
          a(p(i),k) = z
          do j=k+1,n    
            a(p(i),j) = a(p(i),j) - z*a(p(k),j)      
         enddo
        enddo  
      enddo    
      do k=1,n-1
        do i=k+1,n
          b(p(i)) = b(p(i)) - a(p(i),k)*b(p(k))
        enddo
      enddo

!     Backwards substitution
      do i=n,1,-1
        sum = b(p(i))
        do j=i+1,n
          sum = sum - a(p(i),j)*x(j)
        enddo
        x(i) = sum/a(p(i),i)
      enddo
      
    end subroutine gaussPivot
    
    subroutine newtonRaphson (f, ne, na, nr, x, args, res, tol, conv)
        
!     SOLVES SIMULTANEOUS EQUATIONS WITH THE NEWTON-RAPHSON METHOD:

      external :: f
      integer :: ne, na, nr, i, maxIter
      integer, intent(inout) :: conv
      double precision, intent(inout) :: args(na), res(nr), x(ne)
      double precision :: jac(ne,ne), f0(ne), dx(ne), tol
      logical :: converged
      parameter (maxIter=30)
      
      f0(:) = 0.0d0
      dx(:) = 0.0d0
      jac(:,:) = 0.0d0
      
      do i=1,maxIter
        call jacobian(f, ne, na, nr, x, args, res, jac, f0)
        if (sqrt(dot_product(f0, f0) / real(ne)) < tol) then
          goto 1 ! : convergence criteria reached
        endif
        
        call gaussPivot(ne, jac, f0, dx)
        
        x(:) = x(:) + dx(:)                
        
        if (sqrt(dot_product(dx, dx)) < (tol * maxval(abs(x)))) then
          goto 1 ! : convergence criteria reached
        endif
      enddo
      
!     if the loop over maxIter finished without convergence:
      conv = -1
      
1     continue      
                    
      end subroutine newtonRaphson

!-----------------------------------------------------------------------
      subroutine jacobian(f, ne, na, nr, x, args, res, jac, f0)

!     COMPUTES THE JACOBIAN MATRIX FOR f(x)
      
      integer :: ne, na, nr, i
      double precision, intent(inout) :: args(na), res(nr), x(ne) 
      double precision, intent(inout) :: jac(ne,ne), f0(ne)
      double precision :: h, temp, f1(ne)
      interface
        subroutine f(ne, na, nr, x, args, res, zero)
          integer, intent(in) :: ne, na, nr
          double precision, intent(inout) :: args(na), res(nr)
          double precision, intent(inout) :: x(ne), zero(ne)
        end subroutine f
      end interface
      
      f1(:) = 0.0d0
      h = 1e-4
      call f(ne, na, nr, x, args, res, f0)
      do i=1,ne
        temp = x(i)
        x(i) = temp + h
        call f(ne, na, nr, x, args, res, f1)
        x(i) = temp
        jac(:,i) = (f1(:) - f0(:)) / h
      enddo
        
      end subroutine jacobian
      
!---------------------------------------------------------------------------------
!>@brief 1D HEAT EQUATION SOLVED WITH GENERIC SCHEME (NEUMANN BCs)
!> It is assumed that top and bottom boundaries are adiabatic, 
!> so that lambdaOverDx has (nCv-1) size and lambdaOverDx(1) is between Cv=1 and Cv=2
!> The equations are formulated such that ap*tp + ae*Te + ae*Tw = b
!>@param  nCv     : Number of control volumes
!>@param  T0(nCv) : Temperature vector for each control volume in previous time step
!>@param  U(nCv)  : Heat loss coefficient [W/m2*k]
!>                  in order to consider the external area we should do U = U*externalArea/crossArea
!>                  being the crossArea the same for all Cv.
!>@param  Tamb(nCv)           : surrounding or ambient temperature for heat losses/gains [k]
!>@param  lambdaOverDx(nCv-1)  : lambda/dx                 [W/m2*k]
!>@param  rhoCpDxOverDt(nCv)   : rho*cp*dx/dt           [W/m2*k]
!>@param  qv(nCv)              : source term for each Cv (+ rise the temperature) [W/m2]
!>                               include here the heat losses as qv = qv + U*Tambient(i)
!>@param  fTime(constant)      : time relaxation factor           between [0-1]
!>                             (fTime=0 Explicit method, f=1, Implicit method (stable), f=0.5 Crank-Nicolson) 
!>@param relax(constant)  : relaxation factor   [0-1] 
!>                          use 1 by default and only in cases were convergence is difficult use 0.8 or 0.5
!>@return  T(nCv)         : Temperature (dependent variable at intervals dx)      [°C]
!-----------------------------------------------------------------------------------
      subroutine oneDimensionalConduction(nCv,T0,Qv,U,Tamb,LambdaOverDx,RhoCpDxOverDt,fTime,relax,T) !,Tmin,Tmax)                                        
      
        implicit none
        integer, intent(in) :: nCv
        double precision :: Ap(nCv), Atop(nCv), Abot(nCv), B(nCv), Tite(nCv),schemeTerm
        double precision, intent(in) ::  T0(nCv),Qv(nCv),U(nCv),Tamb(nCv),LambdaOverDx(nCv-1),&
        fTime, relax, RhoCpDxOverDt(nCv) ! fTime=0 Explicit method, f=1, Implicit method, f=0.5 Crank-Nicolson
      !  double precision, intent(in) :: Tmin,Tmax ! the minimum and maximum allowed temperatures.
        double precision, intent(inout) ::  T(nCv)
    
        integer :: i                
        
        Tite(1:nCv) = T(1:nCv)
        
        schemeTerm = min((1.d0-fTime),1.0d0)
        
        if(schemeTerm<0.001) schemeTerm = 0.0d0 ! to avoid rounding problems
        
        ! We assume adiabatic at i=1 and i=nCv, 
        
        i = 1
        
        if(nCv==1) then
        
            Atop(i) = 0.0d0
            Abot(i) = 0.0d0
            B(i)    = RhoCpDxOverDt(i) * T0(i)  + U(i)*Tamb(i) + Qv(i)
            Ap(i)   = RhoCpDxOverDt(i) + U(i)
        
        else     
            Atop(i) = -LambdaOverDx(i)*fTime
            Abot(i) = 0.0d0
            B(i)    = (RhoCpDxOverDt(i) + schemeTerm*Atop(i))* T0(i) - Atop(i)*schemeTerm*T0(i+1) + U(i)*Tamb(i)+Qv(i)
            !Ap(i)   = -fTime*Atop(i)+RhoCpDxOverDt(i)+U(i)
            Ap(i) = -Atop(i)+RhoCpDxOverDt(i)+U(i)
        
            !write (20,*) 'INSIDE i=',i,' Qv=',Qv(i),' UA*T=',U(i)*Tamb(i),' ap0=',RhoCpDxOverDt(i),' T0=',T0(i)
        
            do i=2,nCv-1
            
                Atop(i) =  -LambdaOverDx(i)*fTime              
                Abot(i) =  -LambdaOverDx(i-1)*fTime  
                B(i)    = (RhoCpDxOverDt(i) + schemeTerm*Atop(i)+ schemeTerm*Abot(i) )* T0(i) - Atop(i)*schemeTerm*T0(i+1) - Abot(i)*schemeTerm*T(i-1) + U(i)*Tamb(i)+Qv(i)
                !Ap(i)   = -(fTime*Atop(i)+fTime*Abot(i))+RhoCpDxOverDt(i)+U(i)
                 Ap(i)   = -(Atop(i)+Abot(i))+RhoCpDxOverDt(i)+U(i)
            
            enddo    
      
            i = nCv
        
            Atop(i) = 0.0d0
            Abot(i) = -LambdaOverDx(i-1)*fTime
            B(i)    = (RhoCpDxOverDt(i) + schemeTerm*Abot(i))* T0(i) - Abot(i)*schemeTerm*T0(i-1)+ U(i)*Tamb(i) + Qv(i)                
            !Ap(i)   = -fTime*Abot(i)+RhoCpDxOverDt(i)+U(i)
            Ap(i)   = -Abot(i)+RhoCpDxOverDt(i)+U(i)
            
        endif
        
        call solveTDMA(nCv,Abot,Ap,Atop,B,T)                 
                        
        T(1:nCv) = relax*T(1:nCv)+(1-relax)*Tite(1:nCv)   
            
      end subroutine oneDimensionalConduction
      
    subroutine oneDimensionalConductionDV(nCv,T0,Qv,U,Tamb,LambdaOverDx,RhoCpDxOverDt,fTime,T) 
      
        implicit none
        integer, intent(in) :: nCv
        double precision :: Ap(nCv), Atop(nCv), Abot(nCv), B(nCv), schemeTerm
        double precision, intent(in) ::  T0(nCv),Qv(nCv),U(nCv),Tamb(nCv),LambdaOverDx(nCv-1),&
        fTime, RhoCpDxOverDt(nCv) ! fTime=0 Explicit method, f=1, Implicit method, f=0.5 Crank-Nicolson      
        double precision, intent(inout) ::  T(nCv)
    
        integer :: i                               
        
        schemeTerm = min((1.d0-fTime),1.0d0)
        
        if(schemeTerm<0.001) schemeTerm = 0.0d0 ! to avoid rounding problems
        
        ! We assume adiabatic at i=1 and i=nCv, 
        
        i = 1
        
        if(nCv==1) then
        
            Atop(i) = 0.0d0
            Abot(i) = 0.0d0
            B(i)    = RhoCpDxOverDt(i) * T0(i)  + U(i)*Tamb(i) + Qv(i)
            Ap(i)   = RhoCpDxOverDt(i) + U(i)
                    
            
           
            !Ap0(i)  = RhoCpDxOverDt(i)
            !B(i)    = Ap0(i)*T0(i)  + ScDx(i)
            !Ap(i)   = Ap0(i) - SpDx(i)
            
        else     
            Atop(i) = -LambdaOverDx(i)*fTime
            Abot(i) = 0.0d0
            B(i)    = (RhoCpDxOverDt(i) + schemeTerm*Atop(i))* T0(i) - Atop(i)*schemeTerm*T0(i+1) + U(i)*Tamb(i)+Qv(i)           
            Ap(i) = -Atop(i)+RhoCpDxOverDt(i)+U(i)
            
           
            !Ap0(i)  = RhoCpDxOverDt(i)
            !B(i)    = (Ap0(i) + schemeTerm*Atop(i))* T0(i) - Atop(i)*schemeTerm*T0(i+1) + ScDx(i)         
            !Ap(i)   = -Atop(i)+Ap0(i)-SpDx(i)
            
            do i=2,nCv-1
            
                Atop(i) =  -LambdaOverDx(i)*fTime              
                Abot(i) =  -LambdaOverDx(i-1)*fTime  
                B(i)    = (RhoCpDxOverDt(i) + schemeTerm*Atop(i)+ schemeTerm*Abot(i) )* T0(i) - Atop(i)*schemeTerm*T0(i+1) - Abot(i)*schemeTerm*T(i-1) + U(i)*Tamb(i)+Qv(i)            
                Ap(i)   = -(Atop(i)+Abot(i))+RhoCpDxOverDt(i)+U(i)
         
                !B(i)    = (Ap0(i) + schemeTerm*Atop(i)+ schemeTerm*Abot(i) )* T0(i) - Atop(i)*schemeTerm*T0(i+1) - Abot(i)*schemeTerm*T(i-1) + ScDx(i)            
                !Ap(i)   = -(Atop(i)+Abot(i)) + Ap0(i) - SpDx(i)
            
            enddo    
      
            i = nCv
        
            Atop(i) = 0.0d0
            Abot(i) = -LambdaOverDx(i-1)*fTime
            B(i)    = (RhoCpDxOverDt(i) + schemeTerm*Abot(i))* T0(i) - Abot(i)*schemeTerm*T0(i-1)+ U(i)*Tamb(i) + Qv(i)                           
            Ap(i)   = -Abot(i)+RhoCpDxOverDt(i)+U(i)
            
            
            !Ap0(i)  = RhoCpDxOverDt(i)
            !B(i)    = (Ap0(i) + schemeTerm*Abot(i))* T0(i) - Abot(i)*schemeTerm*T0(i-1)+ ScDx(i)                     
            !Ap(i)   = -Abot(i)+ Ap0(i) - SpDx(i)
            
        endif
        
        call solveTDMA(nCv,Abot,Ap,Atop,B,T)                 
                        
        !T(1:nCv) = relax*T(1:nCv)+(1-relax)*Tite(1:nCv)   
            
    end subroutine oneDimensionalConductionDV
     
    subroutine oneDimensionalConductionPatankar(nCv,T0,ScDx,SpDx,LambdaOverDx,RhoCpDxOverDt,fTime,T,alwaysPositive) 
      
        implicit none
        integer, intent(in) :: nCv,alwaysPositive
        double precision :: Ap(nCv), Ap0(nCv), Atop(nCv), Abot(nCv), B(nCv), schemeTerm
        double precision, intent(in) ::  T0(nCv),LambdaOverDx(nCv-1),&
        fTime, RhoCpDxOverDt(nCv), ScDx(nCv), SpDx(nCv) ! fTime=0 Explicit method, f=1, Implicit method, f=0.5 Crank-Nicolson      
        double precision, intent(inout) ::  T(nCv)
    
        integer :: i                               
        
        schemeTerm = min((1.d0-fTime),1.0d0)
        
        if(schemeTerm<0.001) schemeTerm = 0.0d0 ! to avoid rounding problems
        
        ! We assume adiabatic at i=1 and i=nCv, 
        
        i = 1
        
        if(nCv==1) then
        
            Atop(i) = 0.0d0
            Abot(i) = 0.0d0
            Ap0(i)  = RhoCpDxOverDt(i)
            B(i)    = Ap0(i)*T0(i)  + ScDx(i)
            Ap(i)   = Ap0(i) - SpDx(i)
        
        else     
            Atop(i) = -LambdaOverDx( i)*fTime
            Abot(i) = 0.0d0
            Ap0(i)  = RhoCpDxOverDt(i)
            B(i)    = (Ap0(i) + schemeTerm*Atop(i))* T0(i) - Atop(i)*schemeTerm*T0(i+1) + ScDx(i)         
            Ap(i)   = -Atop(i)+Ap0(i)-SpDx(i)
            
            do i=2,nCv-1
                Ap0(i)  = RhoCpDxOverDt(i)
                Atop(i) =  -LambdaOverDx(i)*fTime              
                Abot(i) =  -LambdaOverDx(i-1)*fTime  
                B(i)    = (Ap0(i) + schemeTerm*Atop(i)+ schemeTerm*Abot(i) )* T0(i) - Atop(i)*schemeTerm*T0(i+1) - Abot(i)*schemeTerm*T(i-1) + ScDx(i)            
                Ap(i)   = -(Atop(i)+Abot(i)) + Ap0(i) - SpDx(i)
            
            enddo    
      
            i = nCv
        
            Atop(i) = 0.0d0
            Abot(i) = -LambdaOverDx(i-1)*fTime
            Ap0(i)  = RhoCpDxOverDt(i)
            B(i)    = (Ap0(i) + schemeTerm*Abot(i))* T0(i) - Abot(i)*schemeTerm*T0(i-1)+ ScDx(i)                     
            Ap(i)   = -Abot(i)+ Ap0(i) - SpDx(i)
            
        endif
       
        call solveTDMA(nCv,Abot,Ap,Atop,B,T)                 
                        
        !T(1:nCv) = relax*T(1:nCv)+(1-relax)*Tite(1:nCv)   
            
    end subroutine oneDimensionalConductionPatankar
    

       
    Function BrentRoots(AFunction,x1, x2,n,arguments, Tolerance,maxIterations,valueAtRoot,niter, error )  
    
    !*****************************************************
    !*              Brent Method Function                *
    !* Reference:  BORLAND MATHEMATICAL LIBRARY          *
    !*                                                   *
    !*                F90 version by J-P Moreau, Paris.  *
    !*                       (www.jpmoreau.fr)
    !*                                                   *
    !* Small changes by Mattia Battaglia SPF             *
    !*                                                   *
    !*                                                   *
    !* ------------------------------------------------- *
    !* The purpose is to find a real root of a real      *
    !* function f(x) using Brent method.                 *
    !*                                                   *
    !* INPUTS:  AFunction : non-linear input function    *
    !*          n         : number of additional (constant) input arguments *
    !*          arguments : additional arguments
    !*          x1,x2     : interval of root             *
    !*          Tolerance : desired accuracy for root    *
    !*          maxIter   : maximum number of iterations *
    !*                                                   *
    !* OUTPUTS: The function returns the root value      *
    !*          ValueAtRoot : value of f(root)           *
    !*          niter    : number of done iterations     *
    !*          error    : =0, all OK                    *
    !*                   : =1, no root found in interval *
    !*                   : =2, no more iterations !      *
    !*****************************************************
        
        parameter (FPP = 1.d-11, nearzero = 1.d-20)
        interface
        subroutine AFunction(T,n,arguments,name)
          integer, intent(in) :: n
          double precision, intent(in) :: T
          double precision, intent(in),dimension(:) :: arguments
          double precision, intent(inout) :: name
          
        end subroutine AFunction
        end interface
        double precision ::  x1,x2,Tolerance,valueAtRoot
        integer,intent(inout) :: error,nIter
        double precision,dimension(:) :: arguments
        double precision :: resultat, AA, BB, CC, DD, EE, FA, FB, FC, Tol1, PP, QQ, RR, SS, xm
        integer :: i, done, n

        i = 0
        done = 0
        error = 0
        AA = x1  
        BB = x2 
        call AFunction(AA,n,arguments,FA) 
        call AFunction(BB,n,arguments,FB)
        if (RootBracketed(FA,FB).eq.0) THEN
            error = 1
        else 
            FC = FB
            do while (done.eq.0.and.i < maxIterations)
                if (RootBracketed(FC,FB).eq.0) then
                    CC = AA 
                    FC = FA
                    DD = BB - AA
                    EE = DD
                endif
                if (dabs(FC) < dabs(FB)) then
                    AA = BB 
                    BB = CC 
                    CC = AA
                    FA = FB
                    FB = FC
                    FC = FA
                endif
                Tol1 = 2.0 * FPP * dabs(BB) + 0.5 * Tolerance
                xm = 0.5 * (CC-BB)
                if ((dabs(xm) <= Tol1).or.(dabs(FA) < nearzero)) then
                ! A root has been found
                resultat = BB
                done = 1
                call AFunction(resultat,n,arguments,valueAtRoot)
              else 
                if ((dabs(EE) >= Tol1).and.(dabs(FA) > dabs(FB))) then
              SS = FB/ FA
              if (dabs(AA - CC) < nearzero) then
                PP = 2.0 * xm * SS
                QQ = 1.0 - SS
              else 
                QQ = FA/FC
                RR = FB /FC
                PP = SS * (2.0 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1.0))
                QQ = (QQ - 1.0) * (RR - 1.0) * (SS - 1.0)
              endif
              if (PP > nearzero) QQ = -QQ
              PP = dabs(PP);
              if ((2.0 * PP) < Minimum(3.0*xm *QQ-dabs(Tol1 * QQ), dabs(EE * QQ))) then
                EE = DD
                DD = PP/QQ
              else 
                DD = xm
                EE = DD
              endif
            else 
              DD = xm
              EE = DD
            endif
            AA = BB
            FA = FB
            if (dabs(DD) > Tol1) then 
              BB = BB + DD
            else 
              if (xm > 0) then 
	        BB = BB + dabs(Tol1)
              else 
	        BB = BB - dabs(Tol1)
              endif
            endif
            call AFunction(BB,n,arguments,FB)
            i=i+1
          endif
	    end do
        if (i >= maxIterations) error = 2
      endif
      niter = i
      BrentRoots = resultat

    
    end function BrentRoots ! BrentRoots()
        
    Function RootBracketed(x1,x2)
      real*8 x1,x2 
      integer resultat
      if ((x1 > 0.and.x2 > 0).or.(x1 < 0.and.x2 < 0)) then 
        resultat = 0
      else
        resultat = 1
      endif
      RootBracketed = resultat
    end Function RootBracketed

    ! returns the minimum of two real numbers
    Function Minimum(x1,x2) 
      real*8 x1,x2,resultat
      if (x1 < x2) then
        resultat = x1
      else 
        resultat = x2
      endif
      Minimum = resultat
    end Function Minimum
    
    
    
! Newtons method for one dimensional function using analytic form of derivative
! Input arguments:
!       - f         : non linear function f(x,na,arguments(na),name) (name: name of argument x in calling function)
!       - derivative: derivative of f
!       - na        : number of additional arguments
!       - arguments : vector containing additional arguments
!       - maxIter   : maximum number of iterations before abort
!       - iterations: in out variable to track number of iterations for debugging
!       - Tolerance : 0+/-Tolerance determince how close to the root the algorithm goes
!       - x0        : starting value for the iterations
!       - error     : in out variable 1 if the algorithm did not converge
    Function Newton(f,derivative,na,arguments,maxIter,iterations,Tolerance,x0,error)  
        double precision :: e,x0,yy,y1,Y,Newton
        interface
        subroutine f(T,n,arguments,name)
          integer, intent(in) :: n
          double precision, intent(in) :: T
          double precision, intent(in),dimension(:) :: arguments
          double precision, intent(inout) :: name
          
        end subroutine f
        end interface
        interface
        subroutine derivative(T,n,arguments,name)
          integer, intent(in) :: n
          double precision, intent(in) :: T
          double precision, intent(in),dimension(:) :: arguments
          double precision, intent(inout) :: name
          
        end subroutine derivative
        end interface
        integer,intent(in) :: maxIter,na
        integer,intent(inout) ::iterations,error
        double precision, intent(in),dimension(na) :: arguments
        double precision, intent(in) :: Tolerance
        iterations=0
        error=0
        yy=1
        ! Get y and y1
        do while (abs(yy)>Tolerance)
            call f(x0,na,arguments,yy)
            call derivative(x0,na,arguments,y1)
            ! Update estimate
            x0=x0-(yy/y1)
            iterations=iterations+1
            if (iterations>=maxIter) THEN
                error = 1
                exit
            endif
            
        
        end do
        Newton = x0    
        return
    end Function Newton

end module solvers
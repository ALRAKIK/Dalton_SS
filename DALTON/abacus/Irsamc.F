C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

      module turin_irsamc
      
      ! logical !

      LOGICAL    toura 
      
      ! character ! 
      
      CHARACTER*8 ANGO , space 

      ! integer !

      ! real ! 
      
      REAL*8     LX, LY, LZ
      REAL*8     ax, ay, az

      end module turin_irsamc
      
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ******************************************************************** ! 
C           
C                       One electron integrals
C
! ******************************************************************** !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

!     it is realted to the overlap matrix from type S 


      SUBROUTINE DISTANA(x1,x2,X,y1,y2,Y,z1,z2,Z,DIST)
      
        use turin_irsamc
  
        implicit none 
  
#include "priunit.h"

! ---------------------------------------------------------------------- !

        ! input !
  
        REAL*8, intent(in)  :: x1 , x2    ! the gaussian coordinates on X 
        REAL*8, intent(in)  :: y1 , y2    ! the gaussian coordinates on Y
        REAL*8, intent(in)  :: z1 , z2    ! the gaussian coordinates on Z 
  
        ! local !
  
        REAL*8 :: X
        REAL*8 :: Y
        REAL*8 :: Z
  
        ! output !
  
        REAL*8, intent(out) :: DIST      ! the distance 

! ---------------------------------------------------------------------- !

C ********************************************************************** C        

!       from her1car.F ! 
  
!       DIFABX = CORAX - CORBX
!       DIFABY = CORAY - CORBY
!       DIFABZ = CORAZ - CORBZ
!       DISTAB = DIFABX*DIFABX + DIFABY*DIFABY + DIFABZ*DIFABZ
  
C ********************************************************************** C

      if (toura) then
  
        if (space == "POLYMER") then
  
          X = x1 - x2                       ! calculate the distance in normal way 
              
          if (abs(X) > 0.5d0*Lx) then       ! PBC 
            X = Lx - abs(X)
          end if

          Y = y1 - y2
          Z = z1 - z2  

          end if
        
      else
  
      X = x1 - x2
      Y = y1 - y2
      Z = z1 - z2
  
      end if
  
        DIST = X*X+Y*Y+Z*Z                  ! the distance 
  
  
!       write(lupri,'(a,3f12.6)')  'CORAX',    x1  
!       write(lupri,'(a,3f12.6)') 'x2,y2,z2  :',x2,y2,z2
!       write(lupri,'(a,3f12.6)') 'X , Lx    :',X
!       write(lupri,'(a,3f12.6)') 'D         :',DIST
!       write(lupri,'(a,3f12.6)')

        End  

C ---------------------------------------------------------------------- C
C ---------------------------------------------------------------------- C
C ---------------------------------------------------------------------- C

!     calculate the coordinates of the barycenter using the equations
!     xp = e1*x1+e2*x2/(e1+e2)


      SUBROUTINE    BarCENT(e1,e2,e1p,e2p,
     &                         gx1,gx2,Xp,
     &                         gy1,gy2,Yp,
     &                         gz1,gz2,Zp)
      
      use turin_irsamc
      
      implicit none
      
#include "priunit.h"

! ----------------------------------------------------------------------!
      
        ! input !
      
         real*8,intent(in)       :: e1  , e2     ! e1/(e1+e2) , e2/(e1+e2)
         real*8,intent(in)       :: e1p , e2p    ! this just to print as the program calculate e1/(e1+e2) and e2/(e1+e2) before
         real*8,intent(in)       :: gx1 , gx2    ! the gaussian coordinates on X 
         real*8,intent(in)       :: gy1 , gy2    ! the gaussian coordinates on Y 
         real*8,intent(in)       :: gz1 , gz2    ! the gaussian coordinates on Z 
      
      
        ! output ! 

         real*8,intent(out)      :: Xp , Yp , Zp ! the barycenter coordinates 

! ----------------------------------------------------------------------!

C ********************************************************************** C

!       from her1car.F !

!      CORPX  = EXPAPI*CORAX + EXPBPI*CORBX
!      CORPY  = EXPAPI*CORAY + EXPBPI*CORBY
!      CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ

C ********************************************************************** C

      if (toura) then
      
        if (space == 'POLYMER') then 

          if (abs(gx1-gx2) > 0.5d0*Lx) then 

              if (gx1 > gx2) then 

                Xp = e1*gx1 + e2*(gx2+Lx)

              else if (gx2 > gx1) then

                Xp = e1*(gx1+Lx)+e2*gx2
              
              else 

                Xp = e1*gx1+e2*gx2 
              
              end if 

          else 

          Xp = e1*gx1+e2*gx2   

          end if 

!       X  = e1*x1 + e2*x2 
        Yp  = e1*gy1 + e2*gy2
        Zp  = e1*gz1 + e2*gz2

        end if
      
      else
      
      Xp  = e1*gx1 + e2*gx2
      Yp  = e1*gy1 + e2*gy2
      Zp  = e1*gz1 + e2*gz2
      
      end if      

      if (gx1 == gx2 ) then 
        Xp = gx1 
      end if 

!        write(lupri,'(a,3f12.6)')  '###################################'
!        write(lupri,'(a,3f12.6)')  ''
!        write(lupri,'(a,3f12.6)')  ''
!        write(LUPRI,'(a,3f12.6)')  "***********************************"
!        write(LUPRI,'(a,3f12.6)')  '               Barycenter'
!        write(LUPRI,'(a,3f12.6)')  "***********************************"
!        write(LUPRI,'(a,3f12.6)')  ''
!        write(lupri,'(a,3f12.6)')  '          e1                 e2'
!        write(lupri,'(a,3f12.6)')  ''
!        write(lupri,'(4x,f12.6,a,f12.6)')     e1p ,'      ', e2p 
!        write(lupri,'(a,3f12.6)')  ''
!        write(lupri,'(a,3f12.6)')  '       gx1         gx2         Xp'
!        write(lupri,'(a,3f12.6)')  ''
!        write(lupri,'(3(f12.6))')     gx1  ,  gx2 ,   Xp
!        write(lupri,'(a,3f12.6)')  ''

        End
C ---------------------------------------------------------------------- C
C ---------------------------------------------------------------------- C
C ---------------------------------------------------------------------- C        

!     calculate the distance between the gaussian and the barycenter 
!     important for P,d,f,... functions and to get the ODC 
!     
!     if these distance is correct, then the overlap and the kinetic 
!     integrals should be ok 

      SUBROUTINE DISTANAPC(Px,Py,Pz,gx1,gy1,gz1,
     &                              gx2,gy2,gz2,
     &                           Pxg1,Pyg1,Pzg1,
     &                           Pxg2,Pyg2,Pzg2)

      use turin_irsamc

      implicit none

#include "pi.h"
#include "priunit.h"      

! ---------------------------------------------------------------------- !

      ! input !

      REAL*8, intent(in) :: Px , Py , Pz
      REAL*8, intent(in) :: gx1 , gy1 , gz1
      REAL*8, intent(in) :: gx2 , gy2 , gz2

      ! local !

      REAL*8             :: ngx1 , ngx2 

      ! output !

      REAL*8, intent(out) :: Pxg1 , Pyg1 , Pzg1
      REAL*8, intent(out) :: Pxg2 , Pyg2 , Pzg2 

! ---------------------------------------------------------------------- !


C ********************************************************************** C

!     from her1odc.F !

!     PAX = CORPX - CORAX
!     PAY = CORPY - CORAY
!     PAZ = CORPZ - CORAZ
!     PBX = CORPX - CORBX
!     PBY = CORPY - CORBY
!     PBZ = CORPZ - CORBZ

C ********************************************************************** C

      if (toura) then

        if (space == "POLYMER") then 
        
          if (abs(gx1-gx2) > 0.5*Lx) then 

            if (gx1 > gx2) then 
              ngx1 = gx1 
              ngx2 = gx2 + Lx 
            end if 

            if (gx2 > gx1) then 
              ngx2 = gx2 
              ngx1 = gx1+Lx 
            end if 

          else 
            ngx1 = gx1 
            ngx2 = gx2 
          end if 
        
        Pxg1 = Px - ngx1 
        Pxg2 = Px - ngx2 

          
        Pyg1 = Py - gy1
        Pzg1 = Pz - gz1
        
        Pyg2 = Py - gy2
        Pzg2 = Pz - gz2  
        
        end if   
      
      else  
      
      Pxg1 = Px - gx1  
      Pyg1 = Py - gy1
      Pzg1 = Pz - gz1
      
      Pxg2 = Px - gx2
      Pyg2 = Py - gy2
      Pzg2 = Pz - gz2 
      
      end if 
       
!       write(LUPRI,'(a,3f12.6)') "____________________________________"
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)')  '                  ODC'
!       write(LUPRI,'(a,3f12.6)') "____________________________________"
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)') '              Bary_centers  ' 
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(3(4x,f12.6))')         Px , Lx - Px 
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)') '              Gaussian locations '
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(3(4x,f12.6))')         gx1 , gx2  
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)') '               Bary-Gaus loc     ' 
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(3(4x,f12.6))')         Pxg1 , Pxg2
!       write(LUPRI,'(a,3f12.6)') "" 

      End 

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C  

!     Calculate the distance between the atoms and the barycenter to get 
!     the nuclear attraction, (xa - xp), we need this distance and the 
!     ODC from the subroutine DISTANAPC (2nd one in one elctron integrals)
!     to get the correct integrals

      SUBROUTINE DISTARCP(xa,xp,Xap,ya,yp,Yap,za,zp,Zap)

      use turin_irsamc

      implicit none

#include "pi.h"
#include "priunit.h"

!----------------------------------------------------------------------!

      ! input !

      REAL*8, intent(in)  :: xa , xp
      REAL*8, intent(in)  :: ya , yp
      REAL*8, intent(in)  :: za , zp

      ! local !

      integer             :: sign
      REAL*8              :: nxa , nxp 

      ! output !

      REAL*8, intent(out) :: Xap
      REAL*8, intent(out) :: Yap
      REAL*8, intent(out) :: Zap

!----------------------------------------------------------------------!

C**********************************************************************C

!     from her1car.F !

!     DIFCPX = COORC(1,IATOMC) - CORPX
!     DIFCPY = COORC(2,IATOMC) - CORPY
!     DIFCPZ = COORC(3,IATOMC) - CORPZ

C**********************************************************************C


      if (toura) then 
      
        if (space == "POLYMER") then
          
          if (abs(xa-xp) > 0.5*Lx) then 

            if (xa > xp) then 
              Xap = xa - (xp+Lx)                
            end if 

            if (xp > xa) then 
              Xap = (xa+Lx) - xp 
            end if 

          else 

            Xap = xa - xp 

          end if 

        
        if (Xap >= 0) then 
          sign =   1.d0 
        else 
          sign = - 1.d0 
        end if 
      
!      Xap = sign * abs(Xap) 

        Xap = sign*(ax**(-2.d0)*(2.d0-2.d0*cos(ax*Xap)))**(0.5d0)


        Yap  = ya - yp
        Zap  = za - zp
      
        end if 
      
      else 
      
      Xap  = xa - xp
      Yap  = ya - yp
      Zap  = za - zp
      
      end if

!     write(LUPRI,'(a,3f12.6)') "-------------------------------------"
!     write(LUPRI,'(a,3f12.6)') '          Nuclear attraction'
!     write(LUPRI,'(a,3f12.6)') "-------------------------------------"
!     write(LUPRI,'(a,3f12.6)') '             Bary_centers'
!     write(LUPRI,'(a,3f12.6)')
!     write(LUPRI,'(3(4x,f12.6))')         xp , Lx - xp
!     write(LUPRI,'(a,3f12.6)') 
!     write(LUPRI,'(a,3f12.6)') '             Nuclei coord'
!     write(LUPRI,'(a,3f12.6)')
!     write(LUPRI,'((11x,f12.6))')         xa    
!     write(LUPRI,'(a,3f12.6)')
!     write(LUPRI,'(a,3f12.6)') '             The distance :' 
!     write(LUPRI,'(a,3f12.6)')
!     write(LUPRI,'((11x,f12.6))')         Xap 
!     write(LUPRI,'(a,3f12.6)')
!     write(LUPRI,'(a,3f12.6)')

      END 

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************************! 
C           
C                       Two electron integrals
C
!**********************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     it is realted to the overlap matrix from type S 

      SUBROUTINE DISTANA2(x1,x2,X,y1,y2,Y,z1,z2,Z,DIST)
      
      use turin_irsamc
  
      implicit none 
  
#include "pi.h"
#include "priunit.h"

!----------------------------------------------------------------------!

        ! input !
  
        REAL*8, intent(in)  :: x1 , x2   ! the gaussian coordinates on X
        REAL*8, intent(in)  :: y1 , y2   ! the gaussian coordinates on Y  
        REAL*8, intent(in)  :: z1 , z2   ! the gaussian coordinates on Z
  
        ! local !
  
        REAL*8 :: X
        REAL*8 :: Y
        REAL*8 :: Z  
  
  
        ! output !
  
        REAL*8, intent(out) :: DIST

!----------------------------------------------------------------------!

C**********************************************************************C

  !      from her2odc.F !
  
  !     DIFX = CRX1 - CRX2                       
  !     DIFY = CRY1 - CRY2                       
  !     DIFZ = CRZ1 - CRZ2                       
  !     DIST12 = DIFX*DIFX + DIFY*DIFY + DIFZ*DIFZ
  
C**********************************************************************C

      if (toura) then
  
        if (space == "POLYMER") then
  
              X = x1 - x2
  
              if (abs(X) > 0.5d0*Lx) then
                  X = Lx - abs(X)
              end if
              
        end if
  
        Y = y1 - y2
        Z = z1 - z2
  
      else
  
      X = x1 - x2
      Y = y1 - y2
      Z = z1 - z2
  
      end if
  
      DIST = X*X+Y*Y+Z*Z

      End  

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

!     calculate the coordinates of the barycenter using the equations
!     xp = e1*x1+e2*x2/(e1+e2)


      SUBROUTINE    BarCENT2(e1,e2,e1p,e2p,
     &                          gx1,gx2,Xp,
     &                          gy1,gy2,Yp,
     &                          gz1,gz2,Zp)

      use turin_irsamc

      implicit none

#include "priunit.h"
      
!----------------------------------------------------------------------!

      ! input !

      real*8,intent(in)       :: e1  , e2     ! e1/(e1+e2) , e2/(e1+e2)
      real*8,intent(in)       :: e1p , e2p    ! this just to print as the program calculate e1/(e1+e2) and e2/(e1+e2) before
      real*8,intent(in)       :: gx1 , gx2    ! the gaussian coordinates on X 
      real*8,intent(in)       :: gy1 , gy2    ! the gaussian coordinates on Y 
      real*8,intent(in)       :: gz1 , gz2    ! the gaussian coordinates on Z 


      ! output ! 
      
      real*8,intent(out)      :: Xp , Yp , Zp ! the barycenter coordinates 
      
!----------------------------------------------------------------------!
      
C**********************************************************************C
      
!       from her1car.F !
      
!      CORPX  = EXPAPI*CORAX + EXPBPI*CORBX
!      CORPY  = EXPAPI*CORAY + EXPBPI*CORBY
!      CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ
      
C**********************************************************************C
      
      if (toura) then
      
        if (space == 'POLYMER') then 
        
          if (abs(gx1-gx2) > 0.5d0*Lx) then 
          
              if (gx1 > gx2) then 
              
                Xp = e1*gx1 + e2*(gx2+Lx)
              
              else if (gx2 > gx1) then
              
                Xp = e1*(gx1+Lx)+e2*gx2
              
              else 
              
                Xp = e1*gx1+e2*gx2 
              
              end if 
             
          else 
          
          Xp = e1*gx1+e2*gx2   
          
          end if 
         
!       Xp  = e1*gx1 + e2*gx2 
        Yp  = e1*gy1 + e2*gy2
        Zp  = e1*gz1 + e2*gz2
         
        end if
       
      else
      
      Xp  = e1*gx1 + e2*gx2
      Yp  = e1*gy1 + e2*gy2
      Zp  = e1*gz1 + e2*gz2
      
      end if      
       
      if (gx1 == gx2 ) then 
        Xp = gx1 
      end if 
       
      End
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

!     calculate the distance between the gaussian and the barycenter 
!     important for P,d,f,... functions and to get the ODC 
!     
!     if these distance is correct, then the overlap and the kinetic 
!     integrals should be ok 

      SUBROUTINE DISTANAPC2(Px,Py,Pz,gx1,gy1,gz1, 
     &                      gx2,gy2,gz2,
     &                      Pxg1,Pyg1,Pzg1,Pxg2,Pyg2,Pzg2)

      use turin_irsamc

      implicit none

#include "pi.h"
#include "priunit.h"      

!----------------------------------------------------------------------!

      ! input !

      REAL*8, intent(in) :: Px , Py , Pz
      REAL*8, intent(in) :: gx1 , gy1 , gz1
      REAL*8, intent(in) :: gx2 , gy2 , gz2

      ! local !

      REAL*8             :: ngx1 , ngx2

      ! output !

      REAL*8, intent(out) :: Pxg1 , Pyg1 , Pzg1
      REAL*8, intent(out) :: Pxg2 , Pyg2 , Pzg2 

!----------------------------------------------------------------------!


C**********************************************************************C

!     from her1odc.F !

!     PAX = CORPX - CORAX
!     PAY = CORPY - CORAY
!     PAZ = CORPZ - CORAZ
!     PBX = CORPX - CORBX
!     PBY = CORPY - CORBY
!     PBZ = CORPZ - CORBZ

C**********************************************************************C

      if (toura) then
     
        if (space == "POLYMER") then 
       
          if (abs(gx1-gx2) > 0.5*Lx) then 

            if (gx1 > gx2) then 
              ngx1 = gx1 
              ngx2 = gx2 + Lx 
            end if 

            if (gx2 > gx1) then 
              ngx2 = gx2 
              ngx1 = gx1+Lx 
            end if 

          else 
            ngx1 = gx1 
            ngx2 = gx2 
          end if 

        Pxg1 = Px - ngx1 
        Pxg2 = Px - ngx2  
          
        Pyg1 = Py - gy1
        Pzg1 = Pz - gz1
         
        Pyg2 = Py - gy2
        Pzg2 = Pz - gz2  
         
        end if   
       
      else  
     
      Pxg1 = Px - gx1  
      Pyg1 = Py - gy1
      Pzg1 = Pz - gz1
     
      Pxg2 = Px - gx2
      Pyg2 = Py - gy2
      Pzg2 = Pz - gz2 
     
      end if 

!       write(LUPRI,'(a,3f12.6)') "____________________________________"
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)')  '                  ODC'
!       write(LUPRI,'(a,3f12.6)') "____________________________________"
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)') '              Bary_centers  ' 
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(3(4x,f12.6))')         Px  , Lx - Px 
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)') '              Gaussian locations '
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(3(4x,f12.6))')         gx1 , gx2  
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)') '               Bary-Gaus loc     ' 
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(3(4x,f12.6))')         Pxg1 , Pxg2
!       write(LUPRI,'(a,3f12.6)') "" 
    
      End 

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C 

!     Calculate the distance between the two barycenters to get 
!     the nuclear attraction, (xp1 - xp2), we need this distance and the 
!     ODC from the subroutine DISTANAPC (2nd one in two elctron integrals)
!     to get the correct integrals

      SUBROUTINE DISTARCPP(xp1,xp2,Xpp,yp1,yp2,Ypp,zp1,zp2,Zpp)

      use turin_irsamc

      implicit none

#include "pi.h"
#include "priunit.h"

!----------------------------------------------------------------------!

! input !

      REAL*8, intent(in)  :: xp1 , xp2
      REAL*8, intent(in)  :: yp1 , yp2
      REAL*8, intent(in)  :: zp1 , zp2

! local !

      integer             :: sign 

! output !

      REAL*8, intent(out) :: Xpp
      REAL*8, intent(out) :: Ypp
      REAL*8, intent(out) :: Zpp

!----------------------------------------------------------------------!

C**********************************************************************C

!     from her2her.F !

!          PQXI = PX - SGN34X*COOR34(IOD34,1)
!          PQYI = PY - SGN34Y*COOR34(IOD34,2)
!          PQZI = PZ - SGN34Z*COOR34(IOD34,3)

C**********************************************************************C


      if (toura) then 
  
        if (space == "POLYMER") then
            
           Xpp = xp1 - xp2  

          if (Xpp >= 0) then 
            sign =   1.d0 
          else 
            sign = - 1.d0 
          end if

!          Xpp = sign * abs(Xpp) 

          Xpp = sign*(ax**(-2.d0)*(2.d0-2.d0*cos(ax*Xpp)))**(0.5d0)

          Ypp  = yp1 - yp2
          Zpp  = zp1 - zp2
  
        end if 
  
      else 
  
      Xpp  = xp1 - xp2
      Ypp  = yp1 - yp2
      Zpp  = zp1 - zp2
  
      end if

!       write(LUPRI,'(a,3f12.6)') "-------------------------------------"
!       write(LUPRI,'(a,3f12.6)') '          Electron repultions'
!       write(LUPRI,'(a,3f12.6)') "-------------------------------------"
!       write(LUPRI,'(a,3f12.6)') '             Bary_centers'
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(3(4x,f12.6))')         xp1 , xp2
!       write(LUPRI,'(a,3f12.6)') 
!       write(LUPRI,'(a,3f12.6)') '             The distance :' 
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'((11x,f12.6))')         Xpp  
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)')

      END 
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************************! 
C           
C                       Nuclear repulsion
C
!**********************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE DISTNUC(x1,x2,y1,y2,z1,z2,D)
          
      use turin_irsamc

      implicit none 

#include "pi.h"
#include "priunit.h"

      REAL*8, intent(in)  :: x1 , x2         ! the coordinates on X
      REAL*8, intent(in)  :: y1 , y2         ! the coordinates on Y
      REAL*8, intent(in)  :: z1 , z2         ! the coordinates on Z 


      REAL*8, intent(out) :: D               ! the distance

      ! local ! 

      REAL*8              :: X   , Y   , Z 


      if (toura) then 

        if (space == "POLYMER") then
        
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)


!          X = x1 - x2
          Y = y1 - y2
          Z = z1 - z2

        end if
        
      else 
              
      X = x1 - x2
      Y = y1 - y2
      Z = z1 - z2 
      
      end if
    
      D = sqrt(X*X+Y*Y+Z*Z)

!          write(LUPRI,'(a,3f12.6)') 'x1,x2,X   :', x1,x2 , X 
!          write(LUPRI,'(a,3f12.6)') 'y1,y2,Y   :', y1,y2 , Y 
!          write(LUPRI,'(a,3f12.6)') 'D         :', D  
!          write(LUPRI,'(a,3f12.6)') '*******************'

      END
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
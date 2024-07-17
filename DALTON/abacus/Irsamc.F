C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

      module turin_irsamc

      implicit none 
      
      ! logical !

      LOGICAL    toura 
      
      LOGICAL    same_center_X , same_center_Y ,same_center_Z
      
      ! character ! 
      
      CHARACTER*8 ANGO , space 

      ! integer !

      INTEGER(KIND=8),parameter ::   decimals = 12 
      
      ! real ! 
      
      REAL*8     LX, LY, LZ
      REAL*8     ax, ay, az

      ! matrix real ! 

      REAL*8,allocatable   ::  over(:,:)

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
!#include "turin.h"

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

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            1D              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "POLYMER") then
  
        X = x1 - x2                       ! calculate the distance in normal way 
              
        if (abs(X) > 0.5d0*Lx) then       ! PBC 
          X = Lx - abs(X)
        end if
    
        Y = y1 - y2
        Z = z1 - z2  

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            2D              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "SLAB") then
  
        X = x1 - x2                       ! calculate the distance in normal way 
                
        if (abs(X) > 0.5d0*Lx) then       ! PBC 
          X = Lx - abs(X)
        end if
      
        Y = y1 - y2

        if (abs(Y) > 0.5d0*Ly) then       ! PBC 
          Y = Ly - abs(Y)
        end if
    
        Z = z1 - z2  
  
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            3D              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "CRYSTAL") then
  
          X = x1 - x2                       ! calculate the distance in normal way 
                  
          if (abs(X) > 0.5d0*Lx) then       ! PBC 
            X = Lx - abs(X)
          end if
        
          Y = y1 - y2
  
          if (abs(Y) > 0.5d0*Ly) then       ! PBC 
            Y = Ly - abs(Y)
          end if
      
          Z = z1 - z2  

          if (abs(Z) > 0.5d0*Lz) then       ! PBC 
            Z = Lz - abs(Z)
          end if
      
          end if  
        
      else
  
      X = x1 - x2
      Y = y1 - y2
      Z = z1 - z2
  
      end if
  
        DIST = X*X+Y*Y+Z*Z                  ! the distance 
  
  
!       write(lupri,'(a,3f12.6)')  'CORAX',    x1 , x2
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
     &                         gz1,gz2,Zp,
     &                         NUMO1 , NUMO2)
      
      use turin_irsamc
      
      implicit none
      
#include "priunit.h"
#include "turin.h"

! ----------------------------------------------------------------------!
      
        ! input !
      
         real*8,intent(in)       :: e1    , e2     ! e1/(e1+e2) , e2/(e1+e2)
         real*8,intent(in)       :: e1p   , e2p    ! this just to print as the program calculate e1/(e1+e2) and e2/(e1+e2) before
         real*8,intent(in)       :: gx1   , gx2    ! the gaussian coordinates on X 
         real*8,intent(in)       :: gy1   , gy2    ! the gaussian coordinates on Y 
         real*8,intent(in)       :: gz1   , gz2    ! the gaussian coordinates on Z 
         integer, intent(in)     :: NUMO1 , NUMO2
      
      
        ! output ! 

         real*8,intent(out)      :: Xp , Yp , Zp ! the barycenter coordinates 

! ----------------------------------------------------------------------!

C ********************************************************************** C

!       from her1car.F !

!      CORPX  = EXPAPI*CORAX + EXPBPI*CORBX
!      CORPY  = EXPAPI*CORAY + EXPBPI*CORBY
!      CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ

C ********************************************************************** C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
!      New implimentation after changing the distance 
!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (toura) then
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            1D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

          if (gx1 == gx2) then 
            Xp = gx1 
            if (NUMO1 .ne. NUMO2) then 
!            write(lupri,'(I3,f12.6,I3,f12.6)')NUMO1 ,e1p ,NUMO2 ,e2p
            same_center_X = .TRUE. 
            else
            same_center_X = .FALSE.  
            end if 
          end if 

          if (e1 == e2) then 
            if (abs(gx1-gx2) <= 0.5d0*Lx) then 
              Xp = 0.5*(gx1 + gx2)
            else 
              Xp = 0.5*(gx1 + gx2 + Lx)
            end if 
          end if 

          if (XP >= Lx) then 
            Xp = Xp - Lx
          end if 

!       Xp  = e1*gx1 + e2*gx2 
        Yp  = e1*gy1 + e2*gy2
        Zp  = e1*gz1 + e2*gz2

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            2D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == 'SLAB') then 

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

        if (gx1 == gx2) then 
          Xp = gx1 
          if (NUMO1 .ne. NUMO2) then 
!            write(lupri,'(I3,f12.6,I3,f12.6)')NUMO1 ,e1p ,NUMO2 ,e2p
          same_center_X = .TRUE. 
          else
          same_center_X = .FALSE.  
          end if 
        end if 

        if (e1 == e2) then 
          if (abs(gx1-gx2) <= 0.5d0*Lx) then 
            Xp = 0.5*(gx1 + gx2)
          else 
            Xp = 0.5*(gx1 + gx2 + Lx)
          end if 
        end if 

        if (XP >= Lx) then 
          Xp = Xp - Lx
        end if 

        if (abs(gy1-gy2) > 0.5d0*Ly) then 

          if (gy1 > gy2) then 

            Yp = e1*gy1 + e2*(gy2+Ly)

          else if (gy2 > gy1) then

            Yp = e1*(gy1+Ly)+e2*gy2
          
          else 

            Yp = e1*gy1+e2*gy2 
          
          end if 

      else 

      Yp = e1*gy1+e2*gy2   

      end if 

      if (gy1 == gy2) then 
        Yp = gy1 
        if (NUMO1 .ne. NUMO2) then 
!            write(lupri,'(I3,f12.6,I3,f12.6)')NUMO1 ,e1p ,NUMO2 ,e2p
        same_center_Y = .TRUE. 
        else
        same_center_Y = .FALSE.  
        end if 
      end if 

      if (e1 == e2) then 
        if (abs(gy1-gy2) <= 0.5d0*Ly) then 
          Yp = 0.5*(gy1 + gy2)
        else 
          Yp = 0.5*(gy1 + gy2 + Lx)
        end if 
      end if 

      if (Yp >= Ly) then 
        Yp = Yp - Ly
      end if 

!       Xp  = e1*gx1 + e2*gx2 
!       Yp  = e1*gy1 + e2*gy2
        Zp  = e1*gz1 + e2*gz2

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            3D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == 'CRYSTAL') then 

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

        if (gx1 == gx2) then 
          Xp = gx1 
          if (NUMO1 .ne. NUMO2) then 
!            write(lupri,'(I3,f12.6,I3,f12.6)')NUMO1 ,e1p ,NUMO2 ,e2p
          same_center_X = .TRUE. 
          else
          same_center_X = .FALSE.  
          end if 
        end if 

        if (e1 == e2) then 
          if (abs(gx1-gx2) <= 0.5d0*Lx) then 
            Xp = 0.5*(gx1 + gx2)
          else 
            Xp = 0.5*(gx1 + gx2 + Lx)
          end if 
        end if 

        if (XP >= Lx) then 
          Xp = Xp - Lx
        end if 

        if (abs(gy1-gy2) > 0.5d0*Ly) then 

          if (gy1 > gy2) then 

            Yp = e1*gy1 + e2*(gy2+Ly)

          else if (gy2 > gy1) then

            Yp = e1*(gy1+Ly)+e2*gy2
          
          else 

            Yp = e1*gy1+e2*gy2 
          
          end if 

      else 

      Yp = e1*gy1+e2*gy2   

      end if 

      if (gy1 == gy2) then 
        Yp = gy1 
        if (NUMO1 .ne. NUMO2) then 
!            write(lupri,'(I3,f12.6,I3,f12.6)')NUMO1 ,e1p ,NUMO2 ,e2p
        same_center_Y = .TRUE. 
        else
        same_center_Y = .FALSE.  
        end if 
      end if 

      if (e1 == e2) then 
        if (abs(gy1-gy2) <= 0.5d0*Ly) then 
          Yp = 0.5*(gy1 + gy2)
        else 
          Yp = 0.5*(gy1 + gy2 + Ly)
        end if 
      end if 

      if (Yp >= Ly) then 
        Yp = Yp - Ly
      end if 

      if (abs(gz1-gz2) > 0.5d0*Lz) then 

      if (gz1 > gz2) then 

        Zp = e1*gz1 + e2*(gz2+Lz)

      else if (gz2 > gz1) then

        Zp = e1*(gz1+Lz)+e2*gz2
        
      else 

        Zp = e1*gz1+e2*gz2 
        
      end if 

      else 
        
      Zp = e1*gz1+e2*gz2   
        
      end if 

      if (gz1 == gz2) then 
        Zp = gz1 
        if (NUMO1 .ne. NUMO2) then 
!          write(lupri,'(I3,f12.6,I3,f12.6)')NUMO1 ,e1p ,NUMO2 ,e2p
        same_center_Z = .TRUE. 
        else
        same_center_Z = .FALSE.  
        end if 
      end if 

      if (e1 == e2) then 
        if (abs(gz1-gz2) <= 0.5d0*Lz) then 
          Zp = 0.5*(gz1 + gz2)
        else 
          Zp = 0.5*(gz1 + gz2 + Lz)
        end if 
      end if 

      if (Zp >= Lz) then 
        Zp = Zp - Lz
      end if 

!       Xp  = e1*gx1 + e2*gx2 
!       Yp  = e1*gy1 + e2*gy2
!       Zp  = e1*gz1 + e2*gz2

        end if
      
      else
      
      Xp  = e1*gx1 + e2*gx2
      Yp  = e1*gy1 + e2*gy2
      Zp  = e1*gz1 + e2*gz2
      
      end if      

  !      write(lupri,'(a,3f12.6)')  '###################################'
  !      write(lupri,'(a,3f12.6)')  ''
  !      write(lupri,'(a,3f12.6)')  ''
  !      write(LUPRI,'(a,3f12.6)')  "***********************************"
  !      write(LUPRI,'(a,3f12.6)')  '               Barycenter'
  !      write(LUPRI,'(a,3f12.6)')  "***********************************"
  !      write(LUPRI,'(a,3f12.6)')  ''
  !      write(lupri,'(a,3f24.16)')  '          e1                 e2'
  !      write(lupri,'(a,3f12.6)')  ''
  !      write(lupri,'(4x,f24.16,a,f24.16)')     e1p ,'      ', e2p 
  !      write(lupri,'(a,3f12.6)')  ''
  !      write(lupri,'(a,3f12.6)')  '       gx1         gx2         Xp'
  !      write(lupri,'(a,3f12.6)')  ''
  !      write(lupri,'(3(f24.16))')     gx1  ,  gx2 ,   Xp
  !      write(lupri,'(a,3f12.6)')  ''

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
     &                           Pxg2,Pyg2,Pzg2,J)

      use turin_irsamc

      implicit none

#include "pi.h"
#include "priunit.h"
!#include "turin.h"      

! ---------------------------------------------------------------------- !

      ! input !

      REAL*8, intent(in) :: Px , Py , Pz
      REAL*8, intent(in) :: gx1 , gy1 , gz1
      REAL*8, intent(in) :: gx2 , gy2 , gz2
      INTEGER,intent(in) :: J

      ! local !

      REAL*8             :: ngx1 , ngx2 
      INTEGER            :: sign 

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
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           
!      New implimentation after changing the distance
!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (toura) then

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           1D             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "POLYMER") then 
        
        call ssd(Px,gx1,Pxg1)
  
        call ssd(Px,gx2,Pxg2)

!!EUC!! 
!        call essd(Px,gx1,Pxg1)
!!EUC!!   
!        call essd(Px,gx2,Pxg2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! force implementation  !   

!        if (Pxg1 > 0) then 
!          Pxg1 = - Pxg1 
!        end if 
! 
!        if (Pxg2 < 0) then 
!          Pxg2 = - Pxg2 
!        end if 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

        Pyg1 = Py - gy1
        Pzg1 = Pz - gz1
        
        Pyg2 = Py - gy2
        Pzg2 = Pz - gz2  

        end if   

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           2D             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "SLAB") then 
        
          call ssd(Px,gx1,Pxg1)
    
          call ssd(Px,gx2,Pxg2)
          
          call ssd(Py,gy1,Pyg1)
    
          call ssd(Py,gy2,Pyg2)
          
          Pzg1 = Pz - gz1          
          Pzg2 = Pz - gz2  
  
          end if        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           3D             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (space == "CRYSTAL") then 
        
            call ssd(Px,gx1,Pxg1)
      
            call ssd(Px,gx2,Pxg2)
            
            call ssd(Py,gy1,Pyg1)
      
            call ssd(Py,gy2,Pyg2)
            
            call ssd(Pz,gz1,Pzg1)
      
            call ssd(Pz,gz2,Pzg2)
    
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
!       write(LUPRI,'(3(4x,f24.16))')         Px 
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)') '              Gaussian locations '
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(3(4x,f24.16))')         gx1 , gx2  
!       write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'(a,3f12.6)') '               Bary-Gaus loc     ' 
!       write(LUPRI,'(a,3f12.6)') 
!       write(LUPRI,'(3(4x,f24.16))')         Pxg1 , Pxg2
!       write(LUPRI,'(a,3f12.6)') ""

      End 

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C  

!     Calculate the distance between the atoms and the barycenter to get 
!     the nuclear attraction, (xa - xp), we need this distance and the 
!     ODC from the subroutine DISTANAPC (2nd one in one elctron integrals)
!     to get the correct integrals

      SUBROUTINE DISTARCP(xa,xp,Xap,ya,yp,Yap,za,zp,Zap,J)

      use turin_irsamc

      implicit none

#include "pi.h"
#include "priunit.h"
!#include "turin.h"

!----------------------------------------------------------------------!

      ! input !

      REAL*8, intent(in)  :: xa , xp
      REAL*8, intent(in)  :: ya , yp
      REAL*8, intent(in)  :: za , zp
      integer,intent(in)  ::  J 
      ! local !

      integer             :: sign
      REAL*8              :: nxp 

      ! output !

      REAL*8, intent(out) :: Xap
      REAL*8, intent(out) :: Yap
      REAL*8, intent(out) :: Zap

!----------------------------------------------------------------------!

C**********************************************************************C

!     from her1car.F !

!     DIFCPX = COORC(1,IATOMC) - CORPX      ! 
!     DIFCPY = COORC(2,IATOMC) - CORPY      !
!     DIFCPZ = COORC(3,IATOMC) - CORPZ      !

C**********************************************************************C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!      New implimentation after changing the distance                  !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (toura) then 
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           1D             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "POLYMER") then
          
        call SSD(xa,xp,Xap)
        
          if (Xap >=  0) then 
            sign =     1 
          else 
            sign =  -   1 
          end if 

          if (xa == xp ) Xap = 0.d0 

          Xap = sign*abs(Xap)

          Xap = sign*(ax**(-2.d0)*(2.d0-2.d0*cos(ax*Xap)))**(0.5d0)

!          if (same_center_X .and. abs(xa-xp) == 0.5d0*Lx) Xap = 0.d0 
                
          Yap = ya - yp 
          Zap = za - zp

        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           2D             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "SLAB") then
          
          call SSD(xa,xp,Xap)
          
            if (Xap >=  0) then 
              sign =     1 
            else 
              sign =  -   1 
            end if 
  
            if (xa == xp ) Xap = 0.d0 
  
            Xap = sign*abs(Xap)
  
            Xap = sign*(ax**(-2.d0)*(2.d0-2.d0*cos(ax*Xap)))**(0.5d0)
  
            if (same_center_X) Xap = 0.d0 
                  
            call SSD(ya,yp,Yap)
          
            if (Yap >=  0) then 
              sign =     1 
            else 
              sign =  -   1 
            end if 
  
            if (ya == yp ) Yap = 0.d0 
  
            Yap = sign*abs(Yap)
  
            Yap = sign*(ay**(-2.d0)*(2.d0-2.d0*cos(ay*Yap)))**(0.5d0)
  
            if (same_center_Y) Yap = 0.d0 
                  
            Zap = za - zp
  
          end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           3D             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (space == "CRYSTAL") then
          
            call SSD(xa,xp,Xap)
            
              if (Xap >=  0) then 
                sign =     1 
              else 
                sign =  -   1 
              end if 
    
              if (xa == xp ) Xap = 0.d0 
    
              Xap = sign*abs(Xap)
    
              Xap = sign*(ax**(-2.d0)*(2.d0-2.d0*cos(ax*Xap)))**(0.5d0)
    
              if (same_center_X) Xap = 0.d0 
                    
              call SSD(ya,yp,Yap)
            
              if (Yap >=  0) then 
                sign =     1 
              else 
                sign =  -   1 
              end if 
    
              if (ya == yp ) Yap = 0.d0 
    
              Yap = sign*abs(Yap)
    
              Yap = sign*(ay**(-2.d0)*(2.d0-2.d0*cos(ay*Yap)))**(0.5d0)
    
              if (same_center_Y) Yap = 0.d0 
                    
              call SSD(za,zp,Zap)
            
              if (Zap >=  0) then 
                sign =     1 
              else 
                sign =  -   1 
              end if 
    
              if (za == zp ) Zap = 0.d0 
    
              Zap = sign*abs(Zap)
    
              Zap = sign*(az**(-2.d0)*(2.d0-2.d0*cos(az*Zap)))**(0.5d0)
    
              if (same_center_Z) Zap = 0.d0 
    
            end if

      else 
      
      Xap  = xa - xp
      Yap  = ya - yp
      Zap  = za - zp

      end if

!      write(LUPRI,'(a,3f12.6)') "-------------------------------------"
!      write(LUPRI,'(a,3f12.6)') '          Nuclear attraction'
!      write(LUPRI,'(a,3f12.6)') "-------------------------------------"
!      write(LUPRI,'(a,3f12.6)') '             Bary_centers'
!      write(LUPRI,'(a,3f12.6)')
!      write(LUPRI,'(11x,2f12.6)')           xp 
!      write(LUPRI,'(a,3f12.6)') 
!      write(LUPRI,'(a,3f12.6)') '             Nuclei coord'
!      write(LUPRI,'(a,3f12.6)')
!      write(LUPRI,'((11x,2f12.6))')         xa  
!      write(LUPRI,'(a,3f12.6)')
!      write(LUPRI,'(a,3f12.6)') '             The distance :' 
!      write(LUPRI,'(a,3f12.6)')
!       write(LUPRI,'((11x,4f24.16))') xa,xp, Xap
!      write(LUPRI,'((11x,4f24.8))') xa,xp, Xap 
!      write(LUPRI,'(a,3f12.6)')
!      write(LUPRI,'(a,3f12.6)')


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

      SUBROUTINE DISTANA2(x1,x2,y1,y2,z1,z2,DIST)
      
      use turin_irsamc
  
      implicit none 
  
#include "pi.h"
#include "priunit.h"
!#include "turin.h"

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

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            1D           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "POLYMER") then
  
          X = x1 - x2
  
          if (abs(X) > 0.5d0*Lx) then
            X = Lx - abs(X)
          end if

          Y = y1 - y2
          Z = z1 - z2
  
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            2D           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "SLAB") then
  
          X = x1 - x2
  
          if (abs(X) > 0.5d0*Lx) then
            X = Lx - abs(X)
          end if

          Y = y1 - y2

          if (abs(Y) > 0.5d0*Ly) then
            Y = Ly - abs(Y)
          end if

          Z = z1 - z2
  
        end if        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            3D           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "CRYSTAL") then
  
          X = x1 - x2
  
          if (abs(X) > 0.5d0*Lx) then
            X = Lx - abs(X)
          end if

          Y = y1 - y2

          if (abs(Y) > 0.5d0*Ly) then
            Y = Ly - abs(Y)
          end if

          Z = z1 - z2

          if (abs(Z) > 0.5d0*Lz) then
            Z = Lz - abs(Z)
          end if
  
        end if
  
  
      else
  
      X = x1 - x2
      Y = y1 - y2
      Z = z1 - z2
  
      end if
  
      DIST = X*X+Y*Y+Z*Z

!      write(lupri,'(a,3f12.6)') 'D         :',DIST

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
!#include "turin.h"
      
!----------------------------------------------------------------------!

      ! input !

      real*8,intent(in)       :: e1  , e2     ! e1/(e1+e2) , e2/(e1+e2)
      real*8,intent(in)       :: e1p , e2p    ! this just to print as the program calculate e1/(e1+e2) and e2/(e1+e2) before
      real*8,intent(in)       :: gx1 , gx2    ! the gaussian coordinates on X 
      real*8,intent(in)       :: gy1 , gy2    ! the gaussian coordinates on Y 
      real*8,intent(in)       :: gz1 , gz2    ! the gaussian coordinates on Z 

      ! local ! 

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
      
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            1D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

          if (gx1 == gx2 ) then 
            Xp = gx1 
          end if 


          if (e1 == e2) then 
            if (abs(gx1-gx2) <= 0.5d0*Lx) then 
              Xp = 0.5*(gx1 + gx2)
            else 
              Xp = 0.5*(gx1 + gx2 + Lx)
            end if 
          end if 


          if (Xp >= Lx) then 
            Xp = Xp  - Lx 
          end if 

!       Xp  = e1*gx1 + e2*gx2 
        Yp  = e1*gy1 + e2*gy2
        Zp  = e1*gz1 + e2*gz2

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            2D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == 'SLAB') then 
        
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

          if (gx1 == gx2 ) then 
            Xp = gx1 
          end if 


          if (e1 == e2) then 
            if (abs(gx1-gx2) <= 0.5d0*Lx) then 
              Xp = 0.5*(gx1 + gx2)
            else 
              Xp = 0.5*(gx1 + gx2 + Lx)
            end if 
          end if 

          if (Xp >= Lx) then 
            Xp = Xp  - Lx 
          end if 

          if (abs(gy1-gy2) > 0.5d0*Ly) then 
          
            if (gy1 > gy2) then 
            
              Yp = e1*gy1 + e2*(gy2+Ly)
            
            else if (gy2 > gy1) then
            
              Yp = e1*(gy1+Ly)+e2*gy2
            
            else 
            
              Yp = e1*gy1+e2*gy2 
            
            end if 
           
        else 
        
        Yp = e1*gy1+e2*gy2   
        
        end if 

        if (gy1 == gy2 ) then 
          Yp = gy1 
        end if 


        if (e1 == e2) then 
          if (abs(gy1-gy2) <= 0.5d0*Ly) then 
            Yp = 0.5*(gy1 + gy2)
          else 
            Yp = 0.5*(gy1 + gy2 + Ly)
          end if 
        end if 

        if (Yp >= Ly) then 
          Yp = Yp  - Ly 
        end if

!       Xp  = e1*gx1 + e2*gx2 
!       Yp  = e1*gy1 + e2*gy2
        Zp  = e1*gz1 + e2*gz2

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            3D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == 'CRYSTAL') then 
        
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

          if (gx1 == gx2 ) then 
            Xp = gx1 
          end if 


          if (e1 == e2) then 
            if (abs(gx1-gx2) <= 0.5d0*Lx) then 
              Xp = 0.5*(gx1 + gx2)
            else 
              Xp = 0.5*(gx1 + gx2 + Lx)
            end if 
          end if 

          if (Xp >= Lx) then 
            Xp = Xp  - Lx 
          end if 

          if (abs(gy1-gy2) > 0.5d0*Ly) then 
          
            if (gy1 > gy2) then 
            
              Yp = e1*gy1 + e2*(gy2+Ly)
            
            else if (gy2 > gy1) then
            
              Yp = e1*(gy1+Ly)+e2*gy2
            
            else 
            
              Yp = e1*gy1+e2*gy2 
            
            end if 
           
        else 
        
        Yp = e1*gy1+e2*gy2   
        
        end if 

        if (gy1 == gy2 ) then 
          Yp = gy1 
        end if 


        if (e1 == e2) then 
          if (abs(gy1-gy2) <= 0.5d0*Ly) then 
            Yp = 0.5*(gy1 + gy2)
          else 
            Yp = 0.5*(gy1 + gy2 + Ly)
          end if 
        end if 

        if (Xp >= Ly) then 
          Yp = Yp  - Ly 
        end if

        if (abs(gz1-gz2) > 0.5d0*Lz) then 
          
          if (gz1 > gz2) then 
          
            Zp = e1*gz1 + e2*(gz2+Lz)
          
          else if (gz2 > gz1) then
          
            Zp = e1*(gz1+Lz)+e2*gz2
          
          else 
          
            Zp = e1*gz1+e2*gz2 
          
          end if 
         
        else 
      
        Zp = e1*gz1+e2*gz2   
      
        end if 

        if (gz1 == gz2 ) then 
          Zp = gz1 
        end if 


        if (e1 == e2) then 
          if (abs(gz1-gz2) <= 0.5d0*Lz) then 
            Zp = 0.5*(gz1 + gz2)
          else 
            Zp = 0.5*(gz1 + gz2 + Lz)
          end if 
        end if 

        if (Zp >= Lz) then 
          Zp = Zp  - Lz 
        end if

!       Xp  = e1*gx1 + e2*gx2 
!       Yp  = e1*gy1 + e2*gy2
!       Zp  = e1*gz1 + e2*gz2

        end if        
       
      else
      
      Xp  = e1*gx1 + e2*gx2
      Yp  = e1*gy1 + e2*gy2
      Zp  = e1*gz1 + e2*gz2
      
      end if      
       
!      write(lupri,'(a,3f12.6)')  '###################################'
!      write(lupri,'(a,3f12.6)')  ''
!      write(lupri,'(a,3f12.6)')  ''
!      write(LUPRI,'(a,3f12.6)')  "***********************************"
!      write(LUPRI,'(a,3f12.6)')  '               Barycenter'
!      write(LUPRI,'(a,3f12.6)')  "***********************************"
!      write(LUPRI,'(a,3f12.6)')  ''
!      write(lupri,'(a,3f12.6)')  '          e1                 e2'
!      write(lupri,'(a,3f12.6)')  ''
!      write(lupri,'(4x,f24.16,a,f24.16)')     e1p ,'      ', e2p 
!      write(lupri,'(a,3f12.6)')  ''
!      write(lupri,'(a,3f24.16)')  '       gx1         gx2         Xp'
!      write(lupri,'(a,3f12.6)')  ''
!       write(lupri,'(3(f24.16))')     gx1  ,  gx2 ,   Xp
!      write(lupri,'(a,3f12.6)')  ''

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
!#include "turin.h"      

!----------------------------------------------------------------------!

      ! input !

      REAL*8, intent(in) :: Px , Py , Pz
      REAL*8, intent(in) :: gx1 , gy1 , gz1
      REAL*8, intent(in) :: gx2 , gy2 , gz2

      ! local !


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!      New implimentation after changing the distance                  !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (toura) then
     
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           1D             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "POLYMER") then 
       
        call ssd(Px,gx1,Pxg1)
        call ssd(Px,gx2,Pxg2)

        Pyg1 = Py - gy1
        Pzg1 = Pz - gz1
         
        Pyg2 = Py - gy2
        Pzg2 = Pz - gz2  

        end if   

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           2D             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (space == "SLAB") then 
       
          call ssd(Px,gx1,Pxg1)
          call ssd(Px,gx2,Pxg2)
  
          call ssd(Py,gy1,Pyg1)
          call ssd(Py,gy2,Pyg2)

          Pzg1 = Pz - gz1
          Pzg2 = Pz - gz2  
  
          end if   
  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           3D             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
          if (space == "CRYSTAL") then 
       
            call ssd(Px,gx1,Pxg1)
            call ssd(Px,gx2,Pxg2)
    
            call ssd(Py,gy1,Pyg1)
            call ssd(Py,gy2,Pyg2)
  
            call ssd(Pz,gz1,Pzg1)
            call ssd(Pz,gz2,Pzg2)
    
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
!       write(LUPRI,'(3(4x,f12.6))')         Px  !, Lx - Px 
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
!#include "turin.h"

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
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            1D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "POLYMER") then

       if (xp1 > xp2) then 
         Xpp = xp1 - xp2
         Xpp = - Xpp 
         if (abs(Xpp) > 0.5d0*Lx) then 
           Xpp = Lx - abs(Xpp)
         end if 
       end if 

        if (xp2 > xp1) then 
          Xpp = xp1 - xp2
          Xpp = - Xpp 
          if (abs(Xpp) > 0.5d0*Lx) then 
            Xpp = Lx - abs(Xpp)
            Xpp = - Xpp 
          end if 
        end if 

        if   (abs(xp1-xp2) > 0.5d0*Lx - 0.00000000001 ) then 
          if (abs(xp1-xp2) < 0.5d0*Lx + 0.00000000001 ) then 
            if (xp1 > xp2) then 
              Xpp = - 0.5d0*Lx 
            else 
              Xpp =   0.5d0*Lx
            end if 
          end if 
        end if 

      if (Xpp >  0) then 
        sign = -  1.d0 
      else  
        sign =    1.d0 
      end if       

      if ((xp1 == xp2))   Xpp = 0.d0

      Xpp = sign*abs(Xpp)

      Xpp = sign*(ax**(-2.d0)*(2.d0-2.d0*cos(ax*Xpp)))**(0.5d0)

        Ypp  = yp1 - yp2
        Zpp  = zp1 - zp2
  
        end if 
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            2D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (space == "SLAB") then

          if (xp1 > xp2) then 
            Xpp = xp1 - xp2
            Xpp = - Xpp 
            if (abs(Xpp) > 0.5d0*Lx) then 
              Xpp = Lx - abs(Xpp)
            end if 
          end if 
   
           if (xp2 > xp1) then 
             Xpp = xp1 - xp2
             Xpp = - Xpp 
             if (abs(Xpp) > 0.5d0*Lx) then 
               Xpp = Lx - abs(Xpp)
               Xpp = - Xpp 
             end if 
           end if 
   
           if   (abs(xp1-xp2) > 0.5d0*Lx - 0.00000000001 ) then 
             if (abs(xp1-xp2) < 0.5d0*Lx + 0.00000000001 ) then 
               if (xp1 > xp2) then 
                 Xpp = - 0.5d0*Lx 
               else 
                 Xpp =   0.5d0*Lx
               end if 
             end if 
           end if 
   
         if (Xpp >  0) then 
           sign = -  1.d0 
         else  
           sign =    1.d0 
         end if       
   
        if ((xp1 == xp2))   Xpp = 0.d0
   
         Xpp = sign*abs(Xpp)
   
         Xpp = sign*(ax**(-2.d0)*(2.d0-2.d0*cos(ax*Xpp)))**(0.5d0)

        if (yp1 > yp2) then 
          Ypp = yp1 - yp2
          Ypp = - Ypp 
          if (abs(Ypp) > 0.5d0*Ly) then 
            Ypp = Ly - abs(Ypp)
          end if 
        end if 
 
         if (yp2 > yp1) then 
           Ypp = yp1 - yp2
           Ypp = - Ypp 
           if (abs(Ypp) > 0.5d0*Ly) then 
             Ypp = Ly - abs(Ypp)
             Ypp = - Ypp 
           end if 
         end if 
 
         if   (abs(yp1-yp2) > 0.5d0*Ly - 0.00000000001 ) then 
           if (abs(yp1-yp2) < 0.5d0*Ly + 0.00000000001 ) then 
             if (yp1 > yp2) then 
               Ypp = - 0.5d0*Ly 
             else 
               Ypp =   0.5d0*Ly
             end if 
           end if 
         end if 
 
       if (Ypp >  0) then 
         sign = -  1.d0 
       else  
         sign =    1.d0 
       end if       
 
      if ((yp1 == yp2))   Ypp = 0.d0
 
       Ypp = sign*abs(Ypp)
 
       Ypp = sign*(ay**(-2.d0)*(2.d0-2.d0*cos(ay*Ypp)))**(0.5d0)

       Zpp  = zp1 - zp2
     
      end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            3D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
      if (space == "CRYSTAL") then

        if (xp1 > xp2) then 
          Xpp = xp1 - xp2
          Xpp = - Xpp 
          if (abs(Xpp) > 0.5d0*Lx) then 
            Xpp = Lx - abs(Xpp)
          end if 
        end if 
 
        if (xp2 > xp1) then 
          Xpp = xp1 - xp2
          Xpp = - Xpp 
          if (abs(Xpp) > 0.5d0*Lx) then 
            Xpp = Lx - abs(Xpp)
            Xpp = - Xpp 
          end if 
        end if 
 
        if   (abs(xp1-xp2) > 0.5d0*Lx - 0.00000000001 ) then 
          if (abs(xp1-xp2) < 0.5d0*Lx + 0.00000000001 ) then 
            if (xp1 > xp2) then 
              Xpp = - 0.5d0*Lx 
            else 
              Xpp =   0.5d0*Lx
            end if 
          end if 
        end if 
 
       if (Xpp >  0) then 
         sign = -  1.d0 
       else  
         sign =    1.d0 
       end if       
 
      if ((xp1 == xp2))   Xpp = 0.d0
 
       Xpp = sign*abs(Xpp)
 
       Xpp = sign*(ax**(-2.d0)*(2.d0-2.d0*cos(ax*Xpp)))**(0.5d0)

        if (yp1 > yp2) then 
          Ypp = yp1 - yp2
          Ypp = - Ypp 
          if (abs(Ypp) > 0.5d0*Ly) then 
            Ypp = Ly - abs(Ypp)
          end if 
        end if 

         if (yp2 > yp1) then 
           Ypp = yp1 - yp2
           Ypp = - Ypp 
           if (abs(Ypp) > 0.5d0*Ly) then 
             Ypp = Ly - abs(Ypp)
             Ypp = - Ypp 
           end if 
         end if 

         if   (abs(yp1-yp2) > 0.5d0*Ly - 0.00000000001 ) then 
           if (abs(yp1-yp2) < 0.5d0*Ly + 0.00000000001 ) then 
             if (yp1 > yp2) then 
               Ypp = - 0.5d0*Ly 
             else 
               Ypp =   0.5d0*Ly
             end if 
           end if 
         end if 

       if (Ypp >  0) then 
         sign = -  1.d0 
       else  
         sign =    1.d0 
       end if       

      if ((yp1 == yp2))   Ypp = 0.d0

       Ypp = sign*abs(Ypp)

       Ypp = sign*(ay**(-2.d0)*(2.d0-2.d0*cos(ay*Ypp)))**(0.5d0)

       if (zp1 > zp2) then 
        Zpp = zp1 - zp2
        Zpp = - Zpp 
        if (abs(Zpp) > 0.5d0*Lz) then 
          Zpp = Lz - abs(Zpp)
        end if 
      end if 

       if (zp2 > zp1) then 
         Zpp = zp1 - zp2
         Zpp = - Zpp 
         if (abs(Zpp) > 0.5d0*Lz) then 
           Zpp = Lz - abs(Zpp)
           Zpp = - Zpp 
         end if 
       end if 

       if   (abs(zp1-zp2) > 0.5d0*Lz - 0.00000000001 ) then 
         if (abs(zp1-zp2) < 0.5d0*Lz + 0.00000000001 ) then 
           if (zp1 > zp2) then 
             Zpp = - 0.5d0*Lz
           else 
             Zpp =   0.5d0*Lz
           end if 
         end if 
       end if 

       if (Zpp >  0) then 
         sign = -  1.d0 
       else  
         sign =    1.d0 
       end if       

      if ((zp1 == zp2))   Zpp = 0.d0

       Zpp = sign*abs(Zpp)

       Zpp = sign*(az**(-2.d0)*(2.d0-2.d0*cos(az*Zpp)))**(0.5d0)
     
      end if

      else 
  
      Xpp  = xp1 - xp2
      Ypp  = yp1 - yp2
      Zpp  = zp1 - zp2
  
      end if
             
!      write(LUPRI,'(a,3f12.6)') "-------------------------------------"
!      write(LUPRI,'(a,3f12.6)') '          Electron repultions'
!      write(LUPRI,'(a,3f12.6)') "-------------------------------------"
!      write(LUPRI,'(a,3f12.6)') '             Bary_centers'
!      write(LUPRI,'(a,3f12.6)')
!      write(LUPRI,'(3(4x,f12.6))')         xp1 , xp2
!      write(LUPRI,'(a,3f12.6)') 
!      write(LUPRI,'(a,3f12.6)') '             The distance :' 
!      write(LUPRI,'(a,3f12.6)')
!      write(LUPRI,'(3(f24.14))')         xp1,xp2,Xpp
!      write(LUPRI,'(3(f24.10))')         Xpp
!      write(LUPRI,'(f24.14,a,L,f24.14,a,L,f24.14)') xp1,' ', center1 ,
!     &           xp2,' ', center2, Xpp
!      if (abs(xp1-xp2) > 0.5d0*Lx - 0.0001 ) then 
!        if (abs(xp1-xp2) < 0.5d0*Lx + 0.0001) then 
!      write(LUPRI,'(f24.14,a,L,f24.14,a,L,f24.14)') xp1,' ', center1 ,
!     &           xp2,' ', center2, Xpp
!        end if 
!      end if 
!      write(LUPRI,'(a,3f12.6)')
!      write(LUPRI,'(a,3f12.6)') 


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
#include "turin.h"

      REAL*8, intent(in)  :: x1 , x2         ! the coordinates on X
      REAL*8, intent(in)  :: y1 , y2         ! the coordinates on Y
      REAL*8, intent(in)  :: z1 , z2         ! the coordinates on Z 


      REAL*8, intent(out) :: D               ! the distance

      ! local ! 

      REAL*8              :: X   , Y   , Z 


      if (toura) then 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            1D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (space == "POLYMER") then
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
          Y = y1 - y2
          Z = z1 - z2

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            2D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (space == "SLAB") then
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
          Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
          Z = z1 - z2

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !            3D            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (space == "CRYSTAL") then
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
          Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
          Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)

        end if

      else 
              
      X = x1 - x2
      Y = y1 - y2
      Z = z1 - z2 
      
      end if
    
      D = sqrt(X*X+Y*Y+Z*Z)

!          write(LUPRI,'(a,L)') "toura  " , toura 
!          write(LUPRI,'(a,3f12.6)') 'x1,x2,X   :', x1,x2 , X 
!          write(LUPRI,'(a,3f12.6)') 'y1,y2,Y   :', y1,y2 , Y 
!          write(LUPRI,'(a,3f12.6)') 'D         :', D  
!          write(LUPRI,'(a,3f12.6)') '*******************'

      END
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************************! 
C           
C                       Short signed distance 
C
!**********************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE SSD(x1,x2,X)

        use turin_irsamc
      
        implicit none

!#include "turin.h"
        
        REAL*8, intent(in)  :: x1 , x2         ! the coordinates on X 
      
        REAL*8, intent(out) :: X               ! the distance
      
        ! local ! 
      
        INTEGER              :: sign
      
      
        X = x1 - x2 

        if (X >= 0) then 
          sign =  1 
        else
          sign = -1 
        end if 
        
        if (abs(X) > 0.5d0*Lx) then 
          X = (-1)*sign*(Lx-abs(X))
        end if 

        END

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

      subroutine truncate(number, decimals)
        implicit none
        real(kind=8), intent(inout) :: number
        integer(kind=8), intent(in) :: decimals
        real(kind=8) :: factor
      
        ! Calculate the factor for truncation
        factor = 10.0d0 ** decimals
      
        ! Truncate the number
        number = int(number * factor,8) / factor
      end subroutine truncate

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

      subroutine para_torus(core)

        use turin_irsamc

        implicit none

        INTEGER, intent(in) :: core 

        ! local ! 

        LOGICAL        topo , POLY 
        CHARACTER*8    pspace 
        REAL*8         pax , pay , paz 
        REAL*8         pLx , pLy , pLz 

!#include "turin.h"
      
!        if (core == 0) then 
!          print*, "i am the core and toura " ,core,toura,ax 
!          if (toura) then 
!            topo   = toura
!            pspace = space 
!            if (space == "POLYMER") then 
!              POLY = .TRUE. 
!            end if 
!          end if 
!          print*, "i am the core and toura " ,core,topo, ax
!        end if 

!        if (core .ne. 0) then 
!          print*, "i am the core and toura " ,core,toura,ax
!          print*, ""
!            toura = topo 
!            space = pspace 
!            if (POLY) then 
!              space = "POLYMER"
!            end if 
!            space = "POLYMER"
!          print*, "i am the core and toura " ,core,toura,ax
!          print*, "----------"
!        end if

      end subroutine para_torus 

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
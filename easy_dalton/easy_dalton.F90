program easy_dalton

      implicit none 
      

      integer               :: i
      integer               :: num_atom , charge , chargea , chargeb 
      integer               :: Froom , Too
      character(len=20)     :: BASIS
      character(len=3)      :: ATOM
      character(len=10)     :: num_atom_str , step_str
      character(len=100)    :: command 
      character(len=20)     :: arg

      logical               :: check_length      = .FALSE.
      logical               :: check_num         = .FALSE.
      logical               :: check_dis         = .FALSE.
      logical               :: check_torus       = .FALSE.
      logical               :: check_ring        = .FALSE.
      logical               :: check_infinite    = .FALSE.
      logical               :: check_all         = .FALSE.
      logical               :: check_move        = .FALSE.
      logical               :: check_AO          = .FALSE.
      logical               :: check_GHOST       = .FALSE.
      logical               :: check_CC          = .FALSE.
      logical               :: check_DC          = .FALSE.




      real*8,parameter      :: pi = DACOS(-1.D0)
      real*8                :: radius
      real*8                :: theta 
      real*8                :: step 
      real*8                :: step_h 
      real*8                :: step_increment
      real*8                :: max_length 
      real*8                :: dis_atom
      real*8                :: AO
      real*8                :: length
      real*8                :: G_D

 
      BASIS    = "STO-3G"
      ATOM     = "H"
      charge   = 1 
      Length   = 1.d0  

      call system("clear")


      if (COMMAND_ARGUMENT_COUNT().eq.0) then
            write(*,*) ""
            write(*,*) "Bienvenue          /)"
            write(*,*) "          /\___/\ (("
            write(*,*) "          \`@_@'/  ))"
            write(*,*) "   Amer   {_:Y:.}_// Alrakik  "
            write(*,*) "----------{_}^-'{_}---------- "
            
            write(*,*) ""
            write(*,*) "This a kind of Manual for your program!"
            write(*,*) "write an input file and put the keywords you want!"
                write(*,*) ""
                write(*,*) "OPTIONS: "
                write(*,*) ""
                write(*,*) " N        Write the number of Atoms or use all"
                write(*,*) " D        Write the distance between the atoms"
                write(*,*) " L        Write the length of the TORUS "
                write(*,*) " Basis    Write the Basis set"
                write(*,*) " Name     Write the name of the atoms"
                write(*,*) " C        Write the charge of the atoms"
                write(*,*) " T        To calculate the Torus    space"
                write(*,*) " R        To calculate the Ring     space"
                write(*,*) " I        To calculate the Infinite space"
                write(*,*) " all      To calculate the Torus, Ring and the Infinite space"
                write(*,*) " Move     Move the atoms around the Torus"
                write(*,*) " AO       Limit to delete the MO's based on the overlap of the AO's"
                write(*,*) " FF       Read the geometry from file and modify it"
                write(*,*) " GHOST    add a Ghost atoms before and after (P function)"
                write(*,*) " CCSDT    calculate using CCSDT instead of HF "
		write(*,*) " Double   Two charge "
                write(*,*) ""
                stop
              else
              call get_command_argument(1,arg)
        endif

            call readff(arg,num_atom,dis_atom,length,charge,check_torus,&
     &                  check_ring,check_infinite,check_all,Froom,Too,  &
     &                  check_Move,step_increment,check_AO,AO,check_Ghost,&
     &                  check_CC,Basis,check_dis,check_num,check_DC,&
     &                  chargea,chargeb)

     
        if (.not.(check_dis .and. (check_num .or. check_all))) then
            write(*,*) ""
            write(*,*) "Au revoir          /)"
            write(*,*) "          /\___/\ (("
            write(*,*) "          \`@_@'/  ))"
            write(*,*) "   Amer   {_:Y:.}_// Alrakik  "
            write(*,*) "----------{_}^-'{_}---------- "
            write(*,*) "You need at least :"
            write(*,*) "                    1 - the Number of atoms"
            write(*,*) "                    2 - the Distance between them"
            write(*,*) ""
            write(*,*) "Default basis set is 'STO-3G'"
            write(*,*) "Default charge  is Hydrogen"
            stop
        end if 

!     -------------------------------------------------------------     !      
!                               Torus 
!     -------------------------------------------------------------     !      

      if (check_torus) then

      call system('mkdir Torus')

      if (check_all) then 

      call system('mkdir Torus/from_2_to_16/')

      num_atom = Froom

      do while (num_atom <= Too)

      open(1,file="Torus.dal")
      write(1,'(a)') "**DALTON INPUT"
      write(1,'(a)') ".RUN WAVE FUNCTIONS"
      write(1,'(a)') ".TORUS"
      write(1,'(f16.12,a,5x,f16.12,a,5x,f16.12,a,5x,a)') num_atom*dis_atom,"d0",&
     &                                          num_atom*dis_atom,"d0",&
     &                                          num_atom*dis_atom,"d0 As ",&
     &                                          "POLYMER"
      write(1,'(a)') "**INTEGRALS"
      write(1,'(a)') "**WAVE FUNCTIONS"
      if (check_CC) then
      write(1,'(a)') ".CC"
      write(1,'(a)') "*CC INPUT"
      write(1,'(a)') ".MP2"
      write(1,'(a)') ".CCSD"
      write(1,'(a)') ".CC(T)"
      else
      write(1,'(a)') ".HF"
      end if
      write(1,'(a)') "*SCF IN"
      write(1,'(a)') ".FOCK ITERATIONS"
      write(1,'(a)') "4"
!      write(1,'(a)') ".NODIIS"
      write(1,'(a)') "*ORBITAL"
      if (check_AO) then
        write(1,'(a)') ".AO DELETE"
        write(1,'(f16.12)') AO
      end if
      write(1,'(a)') ".MOSTART"
      write(1,'(a)') "H1DIAG"
      write(1,'(a)') "**END OF DALTON INPUT"

      close(1)

      call itoa(num_atom,num_atom_str)

      open(2,file="geo_" // trim(adjustl(num_atom_str)) // ".mol")

      write(2,'(a)') "BASIS"
      write(2,'(a)') BASIS
      write(2,'(a,a)') ATOM , "POLYMER"
      write(2,'(a,a,a)') "using the ",  BASIS , "basis"
      if (check_GHOST) then 
      write(2,'(a)') "Atomtypes=3 NoSymmetry"
      else
        if (check_DC) then  
          write(2,'(a)') "Atomtypes=2 NoSymmetry"
        else 
          write(2,'(a)') "Atomtypes=1 NoSymmetry"
        end if   
      end if 
      

      if (check_DC) then 

      write(2,'(a,I3,a,I5)') "Charge=",chargea,".0  Atoms=",num_atom/2

      do i = 0 , num_atom-1 , 2 
        theta = i*dis_atom
        write(2,'(a,f16.12,a,a,a)') ATOM , theta ,"d0","  0.000000d0" ,"  0.000000d0"
      end do

      write(2,'(a,I3,a,I5)') "Charge=",chargeb,".0  Atoms=",num_atom/2

      do i = 1 , num_atom-1 , 2 
        theta = i*dis_atom
        write(2,'(a,f16.12,a,a,a)') ATOM , theta ,"d0","  0.000000d0" ,"  0.000000d0"
      end do

      else 

      write(2,'(a,I3,a,I5)') "Charge=",charge,".0  Atoms=",num_atom

      do i = 0 , num_atom-1
        theta = i*dis_atom
        write(2,'(a,f16.12,a,a,a)') ATOM , theta ,"d0","  0.000000d0" ,"  0.000000d0"
      end do

      end if 


      if (check_GHOST) then 
        write(2,'(a,I2,a,I2,a)') "Charge=",1,".0  Atoms=",num_atom, " GHOST"
        do i = 0 , num_atom-1
        theta = i*dis_atom - G_D 
        if (theta < 0.d0) theta = theta + num_atom*dis_atom
        write(2,'(a,f16.12,a,a,a)') ATOM , theta ,"d0","  0.000000d0" ,"  0.000000d0"
        end do
      end if 

      if (check_GHOST) then 
        write(2,'(a,I2,a,I2,a)') "Charge=",2,".0  Atoms=",num_atom, " GHOST"
        do i = 0 , num_atom-1
        theta = i*dis_atom + G_D
        if (theta >= num_atom*dis_atom) theta = theta - num_atom*dis_atom
        write(2,'(a,f16.12,a,a,a)') ATOM , theta ,"d0","  0.000000d0" ,"  0.000000d0"
        end do
      end if 
    
      close(2)
    
      command = 'dalton -noarch -gb 12 Torus.dal geo_'// trim(adjustl(num_atom_str)) // '.mol'


      call system(command)

      call system('mv *.out Torus/from_2_to_16/')
      call system('rm -rf *.dal *.mol')

      num_atom = num_atom + 1

      end do 

      else 

      open(1,file="Torus.dal")

      write(1,'(a)') "**DALTON INPUT"
      write(1,'(a)') ".RUN WAVE FUNCTIONS"
      write(1,'(a)') ".TORUS"
      if (check_length) then 
      write(1,'(f16.12,a,5x,f16.12,a,5x,f16.12,a,5x,a)') length,"d0",&
     &                              length,"d0",&
     &                              length,"d0  As ",&
     &                              "POLYMER"
      else
      write(1,'(f16.12,a,5x,f16.12,a,5x,f16.12,a,5x,a)') num_atom*dis_atom,"d0",&
     &                              num_atom*dis_atom,"d0",&
     &                              num_atom*dis_atom,"d0  As ",&
     &                              "POLYMER"
      end if 
      write(1,'(a)') "**INTEGRALS"
      write(1,'(a)') "**WAVE FUNCTIONS"
      if (check_CC) then
      write(1,'(a)') ".CC"
      write(1,'(a)') "*CC INPUT"
      write(1,'(a)') ".MP2"
      write(1,'(a)') ".CCSD"
      write(1,'(a)') ".CC(T)"
      else
      write(1,'(a)') ".HF"
      end if      
      write(1,'(a)') "*SCF IN"
      write(1,'(a)') ".FOCK ITERATIONS"
      write(1,'(a)') "4"
!      write(1,'(a)') ".NODIIS"
      write(1,'(a)') "*ORBITAL"
      if (check_AO) then
        write(1,'(a)') ".AO DELETE"
        write(1,'(f16.12)') AO
      end if
      write(1,'(a)') ".MOSTART"
      write(1,'(a)') "H1DIAG"
      write(1,'(a)') "**END OF DALTON INPUT"

      close(1)

      if (check_move) then 
      
      step = step_increment
      if (length >= num_atom*dis_atom) then 
        max_length = length 
      else 
        max_length = num_atom*dis_atom
      end if 

      do while (step  < max_length)
 
      call itoa(num_atom,num_atom_str)
      call itoa2(step,step_str)

      open(2,file="geo_" // trim(adjustl(num_atom_str)) //"_"// trim(adjustl(step_str)) // ".mol")
  
      write(2,'(a)') "BASIS"
      write(2,'(a)') BASIS
      write(2,'(a,a)') ATOM , "POLYMER"
      write(2,'(a,a,a)') "using the ",  BASIS , "basis"
      if (check_GHOST) then 
        write(2,'(a)') "Atomtypes=3 NoSymmetry"
      else
        if (check_DC) then  
          write(2,'(a)') "Atomtypes=2 NoSymmetry"
        else 
          write(2,'(a)') "Atomtypes=1 NoSymmetry"
        end if 
      end if 


      if (check_DC) then 

      write(2,'(a,I3,a,I5)') "Charge=",chargea,".0  Atoms=",num_atom/2

      do i = 0 , num_atom-1 ,2 
        step_h = real(i*dis_atom,8)+step
        if (step_h >= num_atom*dis_atom) step_h = step_h - num_atom*dis_atom
        write(2,'(a,f16.8,a,2x,a,a)') ATOM , step_h ,"d0","  0.000000d0" ,"  0.000000d0"
      end do

      write(2,'(a,I3,a,I5)') "Charge=",chargeb,".0  Atoms=",num_atom/2

      do i = 1 , num_atom-1 ,2 
        step_h = real(i*dis_atom,8)+step
        if (step_h >= num_atom*dis_atom) step_h = step_h - num_atom*dis_atom
        write(2,'(a,f16.8,a,2x,a,a)') ATOM , step_h ,"d0","  0.000000d0" ,"  0.000000d0"
      end do

      else 

      write(2,'(a,I3,a,I5)') "Charge=",charge,".0  Atoms=",num_atom
    
      do i = 0 , num_atom-1
        step_h = real(i*dis_atom,8)+step
        if (step_h >= num_atom*dis_atom) step_h = step_h - num_atom*dis_atom
        write(2,'(a,f16.8,a,2x,a,a)') ATOM , step_h ,"d0","  0.000000d0" ,"  0.000000d0"
      end do

      end if 

      if (check_GHOST) then 
        write(2,'(a,I2,a,I2,a)') "Charge=",1,".0  Atoms=",num_atom , " GHOST"
        do i = 0 , num_atom-1
          step_h = real(i*dis_atom,8)+step - G_D 
          if (step_h >= num_atom*dis_atom) step_h = step_h - num_atom*dis_atom
          write(2,'(a,f16.8,a,2x,a,a)') ATOM , step_h ,"d0","  0.000000d0" ,"  0.000000d0"
        end do

      end if 

      if (check_GHOST) then 
        write(2,'(a,I2,a,I2,a)') "Charge=",2,".0  Atoms=",num_atom , " GHOST"
        do i = 0 , num_atom-1
          step_h = real(i*dis_atom,8)+step + G_D 
          if (step_h >= num_atom*dis_atom) step_h = step_h - num_atom*dis_atom
          write(2,'(a,f16.8,a,2x,a,a)') ATOM , step_h ,"d0","  0.000000d0" ,"  0.000000d0"
        end do

      end if 

      close(2)
  
      command = 'dalton -noarch -gb 12 Torus.dal geo_'// trim(adjustl(num_atom_str)) //"_"// trim(adjustl(step_str)) //  '.mol'

      call system(command)
      
      call system('mv *.out Torus/')
      call system('rm -rf *.mol')

      step  = step + step_increment

      end do 

      call system('rm -rf *.dal')

      else 

      call itoa(num_atom,num_atom_str)

      open(2,file="geo_" // trim(adjustl(num_atom_str)) // ".mol")

      write(2,'(a)') "BASIS"
      write(2,'(a)') BASIS
      write(2,'(a,a)') ATOM , "POLYMER"
      write(2,'(a,a,a)') "using the ",  BASIS , "basis"
      if (check_GHOST) then 
      write(2,'(a)') "Atomtypes=3 NoSymmetry"
      else
        if (check_DC) then 
          write(2,'(a)') "Atomtypes=2 NoSymmetry"
        else 
          write(2,'(a)') "Atomtypes=1 NoSymmetry"
        end if 
      end if

      if (check_DC) then 
        
        write(2,'(a,I3,a,I5)') "Charge=",chargea,".0  Atoms=",num_atom/2
        do i = 0 , num_atom-1 ,2 
          theta = i*dis_atom
          write(2,'(a,f16.12,a,a,a)') ATOM , theta   ,"d0","  0.000000d0" ,"  0.000000d0"
        end do 
        
        write(2,'(a,I3,a,I5)') "Charge=",chargeb,".0  Atoms=",num_atom/2
        do i = 1 , num_atom-1 ,2 
          theta = i*dis_atom
          write(2,'(a,f16.12,a,a,a)') ATOM , theta   ,"d0","  0.000000d0" ,"  0.000000d0"
        end do

      else

        write(2,'(a,I3,a,I5)') "Charge=",charge,".0  Atoms=",num_atom
        do i = 0 , num_atom-1
          theta = i*dis_atom
          write(2,'(a,f16.12,a,a,a)') ATOM , theta   ,"d0","  0.000000d0" ,"  0.000000d0"
        end do 

      end if 

      if (check_GHOST) then 
        write(2,'(a,I2,a,I2,a)') "Charge=",1,".0  Atoms=",num_atom, " GHOST"
        do i = 0 , num_atom-1
          theta = i*dis_atom - G_D 
          if (theta < 0) theta = theta + num_atom*dis_atom
        write(2,'(a,f16.12,a,a,a)') ATOM , theta  ,"d0","  0.000000d0" ,"  0.000000d0"
        end do 
      end if 

      if (check_GHOST) then 
        write(2,'(a,I2,a,I2,a)') "Charge=",2,".0  Atoms=",num_atom, " GHOST"
        do i = 0 , num_atom-1
          theta = i*dis_atom + G_D 
          if (theta < 0) theta = theta - num_atom*dis_atom
        write(2,'(a,f16.12,a,a,a)') ATOM , theta  ,"d0","  0.000000d0" ,"  0.000000d0"
        end do 
      end if 

      close(2)

      command = 'dalton -noarch -gb 12 Torus.dal geo_'// trim(adjustl(num_atom_str)) // '.mol'

      call system(command)
      
      call system('mv *.out Torus/')
      call system('rm -rf *.dal *.mol')
    
      end if

      end if 

      end if 

!     -------------------------------------------------------------     !      
!                               Ring  
!     -------------------------------------------------------------     !      

      if (check_ring) then

      call system('mkdir Ring')

      if (check_all) then 

      call system('mkdir Ring/from_2_to_16/')
      
      num_atom = Froom
      
      do while (num_atom <= Too)
      
      open(1,file="Ring.dal")
      write(1,'(a)') "**DALTON INPUT"
      write(1,'(a)') ".RUN WAVE FUNCTIONS"
      write(1,'(a)') "**INTEGRALS"
      write(1,'(a)') "**WAVE FUNCTIONS"
      if (check_CC) then
      write(1,'(a)') ".CC"
      write(1,'(a)') "*CC INPUT"
      write(1,'(a)') ".MP2"
      write(1,'(a)') ".CCSD"
      write(1,'(a)') ".CC(T)"
      else
      write(1,'(a)') ".HF"
      end if
      write(1,'(a)') "*SCF IN"
      write(1,'(a)') ".FOCK ITERATIONS"
      write(1,'(a)') "4"
!      write(1,'(a)') ".NODIIS"
      write(1,'(a)') "*ORBITAL"
      write(1,'(a)') ".MOSTART"
      write(1,'(a)') "H1DIAG"
      write(1,'(a)') "**END OF DALTON INPUT"
      
      close(1)
      
      call itoa(num_atom,num_atom_str)
      
      open(2,file="geo_" // trim(adjustl(num_atom_str)) // ".mol")
      
      write(2,'(a)') "BASIS"
      write(2,'(a)') BASIS
      write(2,'(a,a)') ATOM , "POLYMER"
      write(2,'(a,a,a)') "using the ",  BASIS , "basis"
      if (check_GHOST) then 
      write(2,'(a)') "Atomtypes=3 NoSymmetry"
      else
        if (check_DC) then 
          write(2,'(a)') "Atomtypes=2 NoSymmetry"
        else 
          write(2,'(a)') "Atomtypes=1 NoSymmetry"
        end if 
      end if 


      if (check_DC) then 

      radius = num_atom*dis_atom/(2*pi)

      write(2,'(a,I3,a,I5)') "Charge=",chargea,".0  Atoms=",num_atom/2
          
      do i = 0 , num_atom-1 , 2 
        theta = (real(i*dis_atom,8)*2*pi)/(num_atom*dis_atom)
        write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000000d0"
      end do 

      write(2,'(a,I3,a,I5)') "Charge=",chargeb,".0  Atoms=",num_atom/2
          
      do i = 1 , num_atom-1 , 2 
        theta = (real(i*dis_atom,8)*2*pi)/(num_atom*dis_atom)
        write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000000d0"
      end do 

      else 

      radius = num_atom*dis_atom/(2*pi)

      write(2,'(a,I3,a,I5)') "Charge=",charge,".0  Atoms=",num_atom
        
      do i = 0 , num_atom-1
        theta = (real(i*dis_atom,8)*2*pi)/(num_atom*dis_atom)
        write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000000d0"
      end do 

      end if 


      if (check_GHOST) then 
      write(2,'(a,I2,a,I2,a)') "Charge=",1,".0  Atoms=",num_atom, " GHOST"
      do i = 0 , num_atom-1
        theta = (real(i*dis_atom-G_D,8)*2*pi)/(num_atom*dis_atom)
        write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000000d0"
      end do 
      end if 

      if (check_GHOST) then 
        write(2,'(a,I2,a,I2,a)') "Charge=",2,".0  Atoms=",num_atom, " GHOST"
        do i = 0 , num_atom-1
          theta = (real(i*dis_atom+G_D,8)*2*pi)/(num_atom*dis_atom)
          write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000000d0"
        end do 
        end if 
  

      close(2)
          
      command = 'dalton -noarch -gb 12 Ring.dal geo_'// trim(adjustl(num_atom_str)) // '.mol'
      
      call system(command)
            
      call system('mv *.out Ring/from_2_to_16/')
      call system('rm -rf *.dal *.mol')
      
      num_atom = num_atom + 1
      
      end do 
      
      else 

      open(1,file="Ring.dal")
    
      write(1,'(a)') "**DALTON INPUT"
      write(1,'(a)') ".RUN WAVE FUNCTIONS"
      write(1,'(a)') "**INTEGRALS"
      write(1,'(a)') "**WAVE FUNCTIONS"
      if (check_CC) then
      write(1,'(a)') ".CC"
      write(1,'(a)') "*CC INPUT"
      write(1,'(a)') ".MP2"
      write(1,'(a)') ".CCSD"
      write(1,'(a)') ".CC(T)"
      else
      write(1,'(a)') ".HF"
      end if
      write(1,'(a)') "*SCF IN"
      write(1,'(a)') ".FOCK ITERATIONS"
      write(1,'(a)') "4"
!      write(1,'(a)') ".NODIIS"
      write(1,'(a)') "*ORBITAL"
      write(1,'(a)') ".MOSTART"
      write(1,'(a)') "H1DIAG"
      write(1,'(a)') "**END OF DALTON INPUT"
      
      close(1)

      if (check_Move) then 

      step = step_increment
      if (length >= num_atom*dis_atom) then 
        max_length = length 
      else 
        max_length = num_atom*dis_atom
      end if 

      do while (step  < max_length)
 
      call itoa(num_atom,num_atom_str)
      call itoa2(step,step_str)

      open(2,file="geo_" // trim(adjustl(num_atom_str)) //"_"// trim(adjustl(step_str)) // ".mol")
  
      write(2,'(a)') "BASIS"
      write(2,'(a)') BASIS
      write(2,'(a,a)') ATOM , "POLYMER"
      write(2,'(a,a,a)') "using the ",  BASIS , "basis"
      if (check_GHOST) then 
      write(2,'(a)') "Atomtypes=3 NoSymmetry"
      else 
      write(2,'(a)') "Atomtypes=1 NoSymmetry"
      end if 
      write(2,'(a,I3,a,I5)') "Charge=",charge,".0  Atoms=",num_atom
    
      radius = num_atom*dis_atom/(2*pi)

      do i = 0 , num_atom-1
        theta = (real(i*dis_atom+step,8)*2*pi)/(num_atom*dis_atom)
        write(2,'(a,f16.12,5x,f16.12,a)') ATOM,radius*cos(theta),radius*sin(theta),"  0.000000d0"
      end do 

      close(2)
  
      command = 'dalton -noarch -gb 12 Ring.dal geo_'// trim(adjustl(num_atom_str)) //"_"// trim(adjustl(step_str)) //  '.mol'

      call system(command)
      
      call system('mv *.out Ring/')
      call system('rm -rf *.mol')

      step  = step + step_increment

      end do 

      else 
            
      call itoa(num_atom,num_atom_str)
      
      open(2,file="geo_" // trim(adjustl(num_atom_str)) // ".mol")
      
      write(2,'(a)') "BASIS"
      write(2,'(a)') BASIS
      write(2,'(a,a)') ATOM , "POLYMER"
      write(2,'(a,a,a)') "using the ",  BASIS , "basis"
      if (check_GHOST) then 
        write(2,'(a)') "Atomtypes=3 NoSymmetry"
      else
        if (check_DC) then 
        write(2,'(a)') "Atomtypes=2 NoSymmetry"
        else 
        write(2,'(a)') "Atomtypes=1 NoSymmetry"
        end if 
      end if 


      if (check_DC) then 

      radius = num_atom*dis_atom/(2*pi)

      write(2,'(a,I3,a,I5)') "Charge=",chargea,".0  Atoms=",num_atom/2
      
      do i = 0 , num_atom-1,2
        theta = (real(i*dis_atom,8)*2*pi)/(num_atom*dis_atom)
        write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000d0"
      end do 
  
      write(2,'(a,I3,a,I5)') "Charge=",chargeb,".0  Atoms=",num_atom/2
      
      do i = 1 , num_atom-1,2
        theta = (real(i*dis_atom,8)*2*pi)/(num_atom*dis_atom)
        write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000d0"
      end do

      else 

      write(2,'(a,I3,a,I5)') "Charge=",charge,".0  Atoms=",num_atom
      
      radius = num_atom*dis_atom/(2*pi)

      do i = 0 , num_atom-1
        theta = (real(i*dis_atom,8)*2*pi)/(num_atom*dis_atom)
        write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000d0"
      end do 

      end if 

      if (check_GHOST) then 
      write(2,'(a,I2,a,I2,a)') "Charge=",1,".0  Atoms=",num_atom, " GHOST"  
      do i = 0 , num_atom-1
      theta = (real(i*dis_atom-G_D,8)*2*pi)/(num_atom*dis_atom)
      write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000d0"
      end do
      end if 

      if (check_GHOST) then 
        write(2,'(a,I2,a,I2,a)') "Charge=",2,".0  Atoms=",num_atom, " GHOST"  
        do i = 0 , num_atom-1
        theta = (real(i*dis_atom+G_D,8)*2*pi)/(num_atom*dis_atom)
        write(2,'(a,f16.12,a,f16.12,a)') ATOM,radius*cos(theta),"  ",radius*sin(theta),"  0.000000d0"
        end do
      end if 
      
      close(2)
      
      command = 'dalton -noarch -gb 12 Ring.dal geo_'// trim(adjustl(num_atom_str)) // '.mol'
      
      call system(command)
            
      call system('mv *.out Ring/')
      call system('rm -rf *.dal *.mol')
            
      end if       

      end if 

      end if 
!     -------------------------------------------------------------     !      
!                               infinite  
!     -------------------------------------------------------------     ! 

      if (check_infinite) then

      call system('mkdir Infinite ')

      if (check_all) then 

      call system('mkdir Infinite/from_2_to_16/')
      
      num_atom = Froom
      
      do while (num_atom <= Too)
      
      open(1,file="Infinite.dal")
      write(1,'(a)') "**DALTON INPUT"
      write(1,'(a)') ".RUN WAVE FUNCTIONS"
      write(1,'(a)') "**INTEGRALS"
      write(1,'(a)') "**WAVE FUNCTIONS"
      if (check_CC) then
      write(1,'(a)') ".CC"
      write(1,'(a)') "*CC INPUT"
      write(1,'(a)') ".MP2"
      write(1,'(a)') ".CCSD"
      write(1,'(a)') ".CC(T)"
      else
      write(1,'(a)') ".HF"
      end if
      write(1,'(a)') "*SCF IN"
      write(1,'(a)') ".FOCK ITERATIONS"
      write(1,'(a)') "4"
!      write(1,'(a)') ".NODIIS"
      write(1,'(a)') "*ORBITAL"
      write(1,'(a)') ".MOSTART"
      write(1,'(a)') "H1DIAG"
      write(1,'(a)') "**END OF DALTON INPUT"
      
      close(1)
      
      call itoa(num_atom,num_atom_str)
      
      open(2,file="geo_" // trim(adjustl(num_atom_str)) // ".mol")

      write(2,'(a)') "BASIS"
      write(2,'(a)') BASIS
      write(2,'(a,a)') ATOM , "POLYMER"
      write(2,'(a,a,a)') "using the ",  BASIS , "basis"
      write(2,'(a)') "Atomtypes=1 NoSymmetry"
      write(2,'(a,I3,a,I5)') "Charge=",charge,".0  Atoms=",num_atom
      
      do i = 0 , num_atom-1
      theta = i*dis_atom
      write(2,'(a,f16.12,a,a,a)') ATOM ,theta  ,"d0","  0.000000d0" ,"  0.000000d0"
          
      end do 
          
      close(2)
          
      command = 'dalton -noarch -gb 12 Infinite.dal geo_'// trim(adjustl(num_atom_str)) // '.mol'
      
      call system(command)
            
      call system('mv *.out Infinite/from_2_to_16/')
      call system('rm -rf *.dal *.mol')
      
      num_atom = num_atom + 1 
      
      end do 
      
      else 
      
      open(1,file="Infinite.dal")
      
      write(1,'(a)') "**DALTON INPUT"
      write(1,'(a)') ".RUN WAVE FUNCTIONS"
      write(1,'(a)') "**INTEGRALS"
      write(1,'(a)') "**WAVE FUNCTIONS"
      write(1,'(a)') ".HF"
      write(1,'(a)') "*ORBITAL"
      write(1,'(a)') ".MOSTART"
      write(1,'(a)') "H1DIAG"
      write(1,'(a)') "**END OF DALTON INPUT"
      
      close(1)
      
      call itoa(num_atom,num_atom_str)
      
      open(2,file="geo_" // trim(adjustl(num_atom_str)) // ".mol")
      
      write(2,'(a)') "BASIS"
      write(2,'(a)') BASIS
      write(2,'(a,a)') ATOM , "POLYMER"
      write(2,'(a,a,a)') "using the ",  BASIS , "basis"
      write(2,'(a)') "Atomtypes=1 NoSymmetry"
      write(2,'(a,I2,a,I2)') "Charge=",charge,".0  Atoms=",num_atom
      
      do i = 0 , num_atom-1
        theta = i*dis_atom
      write(2,'(a,f16.12,a,a,a)') ATOM , theta  ,"d0","  0.000000d0" ,"  0.000000d0"
      
      end do 
      
      close(2)
      
      command = 'dalton -noarch -gb 12 Infinite.dal geo_'// trim(adjustl(num_atom_str)) // '.mol'
      
      call system(command)
            
      call system('mv *.out Infinite/')
      call system('rm -rf *.dal *.mol')
            
      end if 

      end if 

end program easy_dalton

subroutine itoa(i,s)
      implicit none
      integer, intent(in)           :: i
      character(len=10),intent(out) :: s
  
      write(s, '(I0)') i
end subroutine itoa

subroutine itoa2(i,s)
  implicit none
  real*8, intent(in)            :: i
  character(len=10),intent(out) :: s

  write(s, '(f8.4)') i
end subroutine itoa2



subroutine readff(arg,num_atom,dis_atom,Length,charge,&
     &            check_torus,check_ring,check_infinite,&
     &            check_all,Froom,Too,check_Move,step_increment,&
     &            check_AO,AO,check_GHOST,check_CC,Basis,check_dis,&
     &            check_num,check_DC,chargea,chargeb)

      implicit none

      ! ---- ! input  ! ---- ! 
      character (len = 20 ) , intent(in)         :: arg
      ! -------------------- ! 


      ! ---- ! local  ! ---- ! 

      integer                                    :: i , nlines
      character (len = 20 )                      :: s_line 




      ! -------------------- ! 


      ! ---- ! output ! ---- !

      integer , intent(out)                      :: num_atom
      integer , intent(out)                      :: charge 
      real*8  , intent(out)                      :: dis_atom
      real*8  , intent(out)                      :: length
      logical , intent(out)                      :: check_torus      
      logical , intent(out)                      :: check_ring       
      logical , intent(out)                      :: check_infinite  
      logical , intent(out)                      :: check_all 
      integer , intent(out)                      :: Froom ,Too 
      logical , intent(out)                      :: check_Move
      real*8  , intent(out)                      :: step_increment
      logical , intent(out)                      :: check_AO
      real*8  , intent(out)                      :: AO
      logical , intent(out)                      :: check_Ghost
      logical , intent(out)                      :: check_CC
      character(len=20) , intent(out)            :: BASIS
      logical , intent(out)                      :: check_dis
      logical , intent(out)                      :: check_num
      logical , intent(out)                      :: check_DC
      integer , intent(out)                      :: chargea,chargeb

      ! -------------------- !


      ! ---- ! CODE ! ---- !


      nlines = 0

      open (1, file = arg)
      do
        read (1,*, end=1)
          nlines = nlines + 1
        end do 
1     close (1)

      open(1,file = arg)
      do i = 1,nlines
        
      read(1,*,end=2) s_line

      if (s_line == "N") then 
        read(1,*) num_atom
        check_num = .TRUE.
      end if 

      if (s_line == "D") then 
        read(1,*) dis_atom
        check_dis = .TRUE.
      end if

      if (s_line == "L") then 
        read(1,*) Length 
      end if

      if (s_line == "C") then 
        read(1,*) charge 
      end if

      if (s_line == "Torus") then 
        check_torus = .TRUE. 
      end if

      if (s_line == "Ring") then 
        check_ring = .TRUE. 
      end if

      if (s_line == "Infinite") then 
        check_infinite = .TRUE. 
      end if

      if (s_line == "All") then 
        check_all = .TRUE.
        read(1,*) Froom, Too 
      end if

      if (s_line == "Move") then 
        check_Move = .TRUE.
        read(1,*) step_increment
      end if

      if (s_line == "AO") then 
        check_AO = .TRUE.
        read(1,*) AO
      end if

      if (s_line == "Ghost") then 
        check_GHOST = .TRUE.
      end if

      if (s_line == "CCSDT") then 
        check_CC    = .TRUE.
      end if

      if (s_line == "Basis") then 
        read(1,*) BASIS
      end if

      if (s_line == "Double") then 
        read(1,*) chargea , chargeb 
        check_DC    = .TRUE.
      end if

      end do 
2     close(1)

end subroutine readff

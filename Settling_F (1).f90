program Settling
    ! Include the modules
    use iso_fortran_env
    
    implicit none

    ! Initial Variables
    
    character(len=100)                     :: typeKeyword
    character(len=256)                     :: inputFilePath, DensityFilePath, SalandTempFilePath, ReynoldsFilePath, outputFileName        
    integer                                :: i, j, LP  
    integer                                :: iostat
    integer                                :: last_occurence
    character(len=2048)                    :: path
    character(len=1)                       :: char_to_find              = "\"
    integer                                :: position_of_last
    real                                   :: Z, v1
    real                                   :: MP_DENSITY, WATER_DENSITY, DIAMETER, SALINITY, TEMPERATURE
    integer :: Correl
    real, dimension(:), allocatable :: POS, DENS, SAL, TEMP, VELPROPERTY 

    !Find the complied executable to retrieve the folder path
    call getcwd(path)
    path = trim(path) // "\"
    print *, "===============================================" 
    print *, "Finding the .exe file ", trim(path), " check 1"
    print *, "Retrieving the folder path"
    print *, "Folder path imposed, files must be in the same folder as .exe"
    inputFilepath = trim(path) // "input.dat"
    DensityFilePath = trim(path) // "density.dat"
    SalandTempFilePath = trim(path) // "SalandTemp.dat"
    ReynoldsFilePath = trim(path) // "Reynolds.dat"
    OutputFileName = trim(path) // "output.dat"
    
    print *, "================================================="
    print *, "Creating path"
    
    print *, trim(inputFilepath)
    print *, trim(DensityFilePath)
    print *, trim(SalandTempFilePath)
    print *, trim(ReynoldsFilePath)
    print *, trim(OutputFileName)
    
    open(unit=10, file=inputFilepath, status='old', action='read', iostat=iostat)
    if (iostat /= 0) then
        print *, "Error opening input file. The program cannot continue. Check for input.dat in the same folder as .exe"
        close(10)
        stop
    else
        print *, "Input file working properly"
    end if
    
     open(unit=11, file=DensityFilePath, status='old', action='read', iostat=iostat)
     if (iostat /= 0) then
        open(unit=12, file=SalandTempFilePath, status='old', action='read', iostat=iostat)
            if (iostat /= 0) then
                print *, "Error opening density files. Check for Density.dat or SalandTemp in the same folder as .exe"
            stop
            close(11)
            close(12)
            else
                print *, "SalandTemp file was found and accepted to calculate density"
            endif
    else
         print *, "Density.dat file found and accepted to supply the density values"
            
    endif

     open(unit=13, file=ReynoldsFilepath, status='old', action='read', iostat=iostat)
    if (iostat /= 0) then
        print *, "Reynolds file not found. Reynolds table cannot be used as a method (use CORREL : 1 or 2)."
        close(13)
    else
        print *, "Reynolds file was found."
    end if
    print *, "Passed 1"
    print *, "==============================================="
    
    
    ! Print program control panel
    print *, "================================================"
    print *, "=                                              ="
    print *, "=                                              ="
    print *, "=                                              ="
    print *, "=====+======= Program Started =================="
    print *, "=                                              ="
    print *, "=                                              ="
    print *, "=                                              ="
    print *, "================================================"
    
    print *, " "
    
    call Loader(InputFilePath, MP_DENSITY, WATER_DENSITY, DIAMETER, SALINITY, TEMPERATURE, Z, CORREL, POS, SAL, TEMP, ReynoldsFilePath, DENS)
    print *, "Finish input"

    call Calculator(Z, MP_DENSITY, WATER_DENSITY, DIAMETER, SALINITY, TEMPERATURE, CORREL, POS, SAL, TEMP, ReynoldsFilePath, DENS, LP, v1, VelProperty)
    
   
    
  if (Z /= 0) then
    open(30, file=OutputfileName, status='replace', action='write')
    write(30, fmt='(*(g0,:,","))') '     Position', '     Velocity'
    do i = 1, LP
      write(30, *) POS(i), VelProperty(i)
    enddo
    close(30)
    print *, "Position", "      Velocity"
     do i = 1, LP
      print *, POS(i), VelProperty(i)
     enddo
 !   write(30, *) POS(i), VelProperty(i), SAL(i), TEMP(i)!, CD e Reynolds
    ! colocar reynolds
    !, WDProperty(i), &
     !                                SALProperty(i), TempProperty(i), CDProperty(i)
 !   enddo
    
  !  close(30)
    else
        print *, WATER_DENSITY, v1, MP_DENSITY
        open(30, file=OutputfileName, status='replace', action='write')
        write(30, fmt='(*(g0,:,","))') '     Velocity'
        write(30, *) v1
        close(30)
    endif
    
    
    print *, "==================================="
    print *, "The program reached the end"
    print *, "==================================="
    
    !deallocate (POS)
    !deallocate (SAL)
    !deallocate (TEMP)


    contains

    
    
  subroutine Loader(inputFileName, MP_DENSITY, WATER_DENSITY, DIAMETER, SALINITY, TEMPERATURE, Z, CORREL, POS, SAL, TEMP, ReynoldsFilePath, DENS)
    character(len=100), intent(in) :: inputFileName
    character(len=256), intent(in) :: ReynoldsFilePath
    character(len=100) :: line, value, line1
    real :: CDS
    integer :: outputFileUnit
    real :: MP_DENSITY, WATER_DENSITY, DIAMETER, SALINITY, TEMPERATURE
    integer :: colon_position
    character(len=100) :: keyword
    integer :: keyword_pos
    character(len=100), dimension(:), allocatable :: lines
    character(len=100), dimension(:), allocatable :: lines1
    real, dimension(:), allocatable :: POS, DENS, SAL, TEMP 
    integer :: i, j, numLines, numlines1, iostat, n
    real:: u, CDnew
    real :: calculateCdVel
    real :: Z
    integer :: Correl
    real :: p1, t, salT, d1, MAXCOLUMN, POSCHECKDIF, POSCHECK, dp
    logical :: insideBlock
    real, parameter :: gravitational_constant = 9.8        ! m / s^2
    real, parameter :: water_dynamic_viscosity = 8.90E-4   ! Pa s

   open(unit=10, file=InputFilepath, status='old', action='read', iostat=iostat)
  if (iostat /= 0) then
      print *, "Error #01 Reading Input"
        stop
    end if

    ! Count the number of lines in the file
    numLines = 0
print *, "Checkpoint 1"
    do
        read(10, '(A)', iostat=iostat)
        if (iostat /= 0) then
            exit ! Exit the loop if end of file is reached
        else
            numLines = numLines + 1
        end if
    end do

    ! Allocate memory for the array to store lines
    allocate(lines(numLines))

    ! Rewind the file to read from the beginning
    rewind(10)
    
        !Guarantee that the values are zeroed
        MP_DENSITY = 0
        WATER_DENSITY = 0
        DIAMETER = 0
        SALINITY = 0
        TEMPERATURE = 0
        Z = 0 
        MAXCOLUMN = 0
        CORREL = 0
        
    ! Read each line and store it into the array
    insideBlock = .false.    
    do i = 1, numLines
        read(10, '(A)') line

        ! Check if the line indicates the start of the block
        if (index(line, "<BegincalculateSettling>") /= 0) then
            insideBlock = .true.
        end if

        ! Check if the line indicates the end of the block
        if (index(line, "<EndcalculateSettling>") /= 0) then
            insideBlock = .false.
        end if
        
        ! Continue processing only if inside the block
        if (insideBlock) then
            lines(i) = line
            colon_position = index(line, ":")
            if (index(line, "MP_DENSITY:") /= 0) then
                read(line(colon_position+1:), *) MP_DENSITY
            endif    
            if (index(line, "WATER_DENSITY:") /= 0) then
                read(line(colon_position+1:), *) WATER_DENSITY
            endif
            if (index(line, "DIAMETER:") /= 0) then
                read(line(colon_position+1:), *) DIAMETER
            endif
            if (index(line, "SALINITY:") /= 0) then
                read(line(colon_position+1:), *) SALINITY
            endif
            if (index(line, "TEMPERATURE:") /= 0) then
                read(line(colon_position+1:), *) TEMPERATURE
            endif
            if (index(line, "Z:") /= 0) then
                read(line(colon_position+1:), *) Z
            endif
            if (index(line, "CORREL:") /= 0) then
                read(line(colon_position+1:), *) Correl
            endif
            
        endif

    end do
    print *, "Data available:"
    if (MP_DENSITY /= 0) then
        print *, "MP_DENSITY:", MP_DENSITY
    endif
    if (WATER_DENSITY /= 0) then
        print *, "WATER_DENSITY:", WATER_DENSITY
        print *, "Only applicable to instantaneaous calculation without Z"
        print *, "Otherwise, use Density.dat or SalandTemp.dat"
    endif
    if (DIAMETER /= 0) then
        print *, "DIAMETER:", DIAMETER
    endif
        if (Z /= 0) then
        print *, "Z:", Z
        print *, "Water column accepted"
        print *, "Mandatory the use of Density.dat or SalandTemp.dat"
        else
            if (Z == 0) then
                print *, "Z:", Z
                print *, "No water column informed"
                print *, "Use the input.dat input for calculation. Reynolds.dat is also applicable."
                endif
    endif
      

    if (SALINITY /= 0) then
        print *, "SALINITY:", SALINITY
        print *, "Only applicable to instantaneaous calculation without Z"
        print *, "Otherwise, use Density.dat or SalandTemp.dat"
    endif
    if (TEMPERATURE /= 0) then
        print *, "TEMPERATURE:", TEMPERATURE
        print *, "Only applicable to instantaneaous calculation without Z"
        print *, "Otherwise, use Density.dat or SalandTemp.dat"
    endif
    if (CORREL /= 0) then
        print *, "CORREL:", CORREL
        print *, "Correlation method, 1- Correlation using Morrison(2003), 2- Correlation using interpolation, 3- Correlation using the Reynolds.dat table"
    else
        print *, "Correl is mandatory, please decide the method: 0- Correlation using Morrison(2003), 1- Correlation using interpolation, 2- Correlation using the Reynolds.dat table"
    endif
    
    open(unit=12, file=SalandTempFilePath, status='old', action='read', iostat=iostat)
    if (iostat /= 0) then
        print *, "SalandTemp file not used"
    else
    print *, "SalandTemp read"
    ! First, count the number of lines
    numLines = 0
    do
        read(12, '(A)', iostat=iostat) line
        if (iostat /= 0) exit
        numLines = numLines + 1
    end do
    numLines = numLines - 1 !THe header has to come out
    rewind(12)  ! Rewind after counting
    read (12, *)
    ! Allocate arrays based on the number of lines
       allocate(POS(numLines), SAL(numLines), TEMP(numLines))
    ! Now, read each line and extract the densities
    do i = 1, numLines
        read(12, *) p1, salT, t
        print *, "Leitura", p1, salT, t
        POS(i) = p1
        SAL(i) = salT
        TEMP(i) = t

    end do
    !print *, POS, "Check Position"
    !print *, SAL, "Check Salinity"
    !print *, TEMP, "Check Temperature"
    !print *, DENS, "Check Density"
    close(12)
    
    endif
    
    
    !deallocate(POS, DENS)
        
    
    open(unit=11, file=DensityFilePath, status='old', action='read', iostat=iostat)
     if (iostat /= 0) then
        print *, "Passing to Sal and Temp.dat"
     else
         print *, "Reading Density.dat file"
    ! First, count the number of lines
    numLines = 0
    do
        read(11, '(A)', iostat=iostat) line
        if (iostat /= 0) exit
        numLines = numLines + 1
    end do
    rewind(11)  ! Rewind after counting
!    print *, "OK"
    !read (11, *)
    ! Allocate arrays based on the number of lines
    allocate(POS(numLines), DENS(numLines))
    
    ! Now, read each line and extract the densities
   ! rewind(11)
    do i = 1, numLines
        read(11, *) p1, d1  ! Directly read two real numbers separated by space
        !print *, p1, d1
        POS(i) = p1
        DENS(i) = d1
    end do

    !print *, "=============================================================="
    ! Check dT consistency
    POSCHECK = POS(2) - POS(1)
    do i = 2, n - 1
        POSCHECKDIF = POS(i+1) - POS(i)
        if (POSCHECKDIF /= POSCHECK) then 
        print *, "Error in consistency"
        stop
        endif
    enddo
    
    close(11)
    endif
    !deallocate(POS, DENS)
            
    !print *, "results ", CDS, MP_DENSITY, WATER_DENSITY, DIAMETER, SALINITY, TEMPERATURE
    !print *, "===== ALL NECESSARY VARIABLES WERE READ ====="
    end subroutine 
    
    
    !Ler arquivo Reynolds.dat
    
subroutine ReadReynolds(RRey, CDNew, ReynoldsFilePath)
character(len=256) :: ReynoldsFilePath
integer :: numDataPoints, numLines
character(len=100) :: line
real :: reynoldsNumber, CDS, RRey, CDNew
real, dimension(:), allocatable :: Reynoldst, CDt
integer :: i, j

!print *, ReynoldsFilePath
   open(unit=13, file=ReynoldsFilePath, status='old', iostat = iostat)
   if (iostat /= 0) then
       print *, "Error opening the file Reynolds.dat"
       stop
   endif
   
   numLines = 0
   print *, "Checkpoint 1"
       do
        read (13, '(A)', iostat=iostat)
        if (iostat /= 0) then
        exit
       else 
        numLines = numLines + 1
       endif
       enddo       
   print *, "Reynolds opened and read"
        
        
    ! Allocate memory for arrays to store Reynolds and CD data
    allocate(Reynoldst(numLines))
    allocate(CDt(numLines))
 !   print *, "==========================================="
 !   print *, "REYNOLDS WERE READ"
 !   print *, Reynoldst
 !   print *, "==========================================="
 !   print *, "==========================================="
 !   print *, "CD VALUES WERE READ"
 !   print *, CDt
!    print *, "==========================================="
    ! Rewind the file
    rewind(13)
    
    ! Skip the header line
    read(13, *)
    
    ! Read Reynolds and CD data from the file
    do i = 1, numLines-1
        read(13, *) Reynoldst(i), CDt(i)
      !  print *, "Reading Reynolds File ", Reynoldst(i), CDt(i)
        
    end do
    print *, "Matching Reynolds =", RRey
    ! Close the file
    close(11)
   ! print *, "Start loop for Reynolds"
    do j = 1, numLines-1
!        print *, "Check 1"
        if (Reynoldst(j) <= RRey .and. Reynoldst(j+1) >= RRey) then
            CDnew = interpolateCD(reynoldsNumber, Reynoldst(j), Reynoldst(j+1), CDt(j), CDt(j+1))
            print *, "In the range"
            print *, "Reynolds Number:", RRey
            print *, "Corresponding CD:", CDnew
                    
        else if (RRey < 0.2 .and. j == 1) then
            print *, "CD adopted will be 24/Reynolds"
            CDnew = 24/RRey
       !     exit
            
        else if (RRey > 1000 .and. j == numLines) THEN
            print *, "CD adopted will be 0.44"
            CDnew = 0.44
       !     exit
        endif
            
        
    end do

    end subroutine ReadReynolds

 
function interpolateCD(x, x1, x2, y1, y2) result(y)
        real, intent(in) :: x, x1, x2, y1, y2
        real :: y
        
        y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
end function interpolateCD

function polyval(coeffs, x)
    real, dimension(:), intent(in) :: coeffs
    real, intent(in) :: x
    real :: polyval
    
    integer :: i, n
    real(8) :: xpowers
    
    n = size(coeffs)
    polyval = coeffs(n)
    xpowers = 1.0
    
    do i = n-1, 1, -1
        xpowers = xpowers * x
        polyval = polyval + coeffs(i) * xpowers
    end do
    
    return
    
end function polyval


!
!Subrotina para calcular o CD baseado na forma Liu e Evans
!
subroutine cd_sphere(Re, CDNew) 
    real :: Re
    real :: CDNew
    real, dimension(4) :: p
    real :: x1
    
    if (Re <= 0.0) then
        CDNew = 0.0
    else if (Re > 8.0e6) then
        CDNew = 0.2
    else if (Re > 0.0 .and. Re <= 0.5) then
        CDNew = 24.0 / Re
    else if (Re > 0.5 .and. Re <= 100.0) then
        p = [4.22, -14.05, 34.87, 0.658]
        CDNew = polyval(p, 1.0 / Re)
    else if (Re > 100.0 .and. Re <= 1.0e4) then
        p = [-30.41, 43.72, -17.08, 2.41]
        CDNew = polyval(p, 1.0 / log10(Re))
    else if (Re > 1.0e4 .and. Re <= 3.35e5) then
        p = [-0.1584, 2.031, -8.472, 11.932]
        CDNew = polyval(p, log10(Re))
    else if (Re > 3.35e5 .and. Re <= 5.0e5) then
        x1 = log10(Re / 4.5e5)
        CDNew = 91.08 * x1**4 + 0.0764
    else
        p = [-0.06338, 1.1905, -7.332, 14.93]
        CDNew = polyval(p, log10(Re))
    end if
end subroutine cd_sphere


   !rotina para calcular reynolds
    subroutine CalcRey (DIAMETER, v, v0, water_density, water_dynamic_viscosity)
    real :: rei, DIAMETER, v, v0, water_density, water_dynamic_viscosity
    
    Rei = DIAMETER * v0 * water_density / water_dynamic_viscosity
    
    
    end subroutine CalcRey
    
! Rotina para calcular o erro
subroutine Error(v0, v1, C)
real :: C, v0, v1

C = abs(v1-v0)/v1 * 100

end subroutine Error

!Calcular CD por Morisson (2003)
function calculateCd(Re, Cd)
real :: cd, Re, calculateCd
   
        Cd = 24.0 / Re + 2.6 * (Re / 5.0) / (1.0 + (Re / 5.0)**1.52) &
        + 0.411 * (Re / 263.0)**(-7.94) / (1.0 + (Re / 263.0)**(-8.0)) &
        + 0.25 * (Re / 1.0E6) / (1.0 + (Re / 1.0E6))
   return
end function calculateCd


!Routine to calculate Stokes velocity
    subroutine CalcVelStokes(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, water_dynamic_viscosity, v0)
    real :: gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, water_dynamic_viscosity, v0

    v0 = gravitational_constant * DIAMETER **2 * (MP_DENSITY - WATER_DENSITY) / (18*water_dynamic_viscosity)
    !calcular numero de reynolds
    end subroutine CalcVelStokes
    
!Rotina para calculo de velocidade com CD
    function CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, CD, v1)
    
    real, intent(inout) :: v1
    real:: CalcVelCD 
    real :: gravitational_constant 
    real :: DIAMETER 
    real :: MP_DENSITY
    real :: WATER_DENSITY
    real :: CD
    real :: dT
    !print *, MP_Density, WATER_DENSITY, DIAMETER, gravitational_constant, CD !check input global
      v1 = ((4.0 / 3.0) * (MP_Density - WATER_DENSITY ) * DIAMETER * gravitational_constant / WATER_DENSITY * CD)**0.5
           
      return 
!    print *, v1 !check calculation
     !print *, 'Velocity obtained check', v1 
    end function CalcVelCD
    
!Rotina para calcular a densidade baseado em salinidade e temperatura
    subroutine DensWater(Sl, Tp, WATER_DENSITY)
    
    real, parameter :: b0 = 8.24493e-1
    real, parameter :: b1 = -4.0899e-3
    real, parameter :: b2 = 7.6438e-5
    real, parameter :: b3 = -8.2467e-7
    real, parameter :: b4 = 5.3872e-9
    real, parameter :: c0 = -5.72466e-3
    real, parameter :: c1 = 1.0227e-4
    real, parameter :: c2 = -1.6546e-6
    real, parameter :: a0 = 999.842594
    real, parameter :: a1 = 6.793952e-2
    real, parameter :: a2 = -9.095290e-3
    real, parameter :: a3 = 1.001685e-4
    real, parameter :: a4 = -1.120083e-6
    real, parameter :: a5 = 6.536332e-9
    real, parameter :: d0 = 4.8314e-4
    real :: pw, WATER_DENSITY 
    real :: Sl, Tp
        
    pw = a0 + a1*Tp + a2*Tp**2 + a3*Tp**3 + a4*Tp**4 + a5*Tp**5
    
    WATER_DENSITY = pw + (b0 + b1*Tp + b2*Tp**2 + b3*Tp**3 + b4*Tp**4)*Sl + (c0 + c1*Tp + c2*Tp**2)*Sl**(3/2)+d0*Sl**2
    
    end subroutine

    
    !---------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------
    

!Calculator this is the one
subroutine Calculator(Z, MP_DENSITY, WATER_DENSITY, DIAMETER, SALINITY, TEMPERATURE, Correl, POS, SAL, TEMP, REYNOLDSFILEPATH, DENS, LP, v1, VelProperty)

character(len=1) :: userInput
logical :: inputReceived
character(len=256) :: ReynoldsFilePath
real :: DIAMETER, MP_DENSITY, WATER_DENSITY, CD, C, RREY1, RREY, RV0, RV1, Z, MAXCOLUMN, POSCHECK, CDNEW, v0, SALINITY, TEMPERATURE
real, parameter :: gravitational_constant = 9.8        ! m / s^2
real, parameter :: water_dynamic_viscosity = 8.90E-4   ! Pa s
real :: Re, v1, Re1, Sl, Tp, LS, PS, dS, d
integer :: i, Correl, k, j, l, LP
real, dimension(:), allocatable :: Posproperty, VelProperty, WDProperty, SALProperty, TempProperty, CDProperty
real, dimension(:) :: POS, SAL, TEMP, DENS
integer :: STEP
!print *, "Calculator Started"
 !   print *, 'Checking inputs:'
  !  print *, gravitational_constant, " = g"
   ! print *, DIAMETER, " = diameter"
!   print *, MP_DENSITY, " = MP_DENSITY"
 !   print *, water_dynamic_viscosity, ' = viscosity'
!    print *, ReynoldsfilePath
!---------------------------------------------------------------------------------------------------------
! Option 1: Z = 0
!---------------------------------------------------------------------------------------------------------
    if (Z == 0) then
    print *, "-----------------------------------------------------"
    print *, "Calculating Instantaneous Velocity"
    print *, "Z = 0"
    print *, "-----------------------------------------------------"
    
!R = Resultado     

    call CalcVelStokes(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, water_dynamic_viscosity, v0) 
 !   print *, "Attempt 1: Considering Stokes:"
  !  print *, "Velocity Stokes = ", v0
   ! print *, Diameter, V0, Water_density, Water_dynamic_viscosity
    Re = DIAMETER * v0 * WATER_DENSITY/water_dynamic_viscosity
    !print *, "Initial Reynolds =", Re
    
   if (Correl == 1) then
     !  print *, "Calculating Cd using Morrisson"
       RRey = calculateCd(Re, Cd)
      ! print *, "CD = ", CD
       
       Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, CD, v1)
    !   print *, gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, CD, v1
    !   print *, "Calculating new velocity using the CD"
    !   print *, "Velocity 1 = ", v1
        Re1 = DIAMETER * v1 * WATER_DENSITY/water_dynamic_viscosity
    !    print *, "Calculating New Reynolds using V1 =", v1
    !    print *, "Reynolds = ", Re1
        call Error(v0, v1, C)
    !    print *, "Error = ", C
    i = 0
    do while (C > 0.01)
        i = i + 1
     !   print *, "Seeking Error < Tolerance"    
     !   print *, "New Reynolds", Re1, "now will be old Reynlods"
        Re = Re1
        v0 = v1
      !  print *, "Old reynolds = ", Re
            RRey = calculateCd(Re, Cd)
       !     print *, "Calculating Cd using Morrisson, Cd = ", CD
           Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, Cd, v1)
       !    print *, "Velocity 1 = ", v1
          Re1 = DIAMETER * v1 * WATER_DENSITY/water_dynamic_viscosity
        !    print *, "Calculating New Reynolds using V1 =", v1
            call Error(v0, v1, C)
         !   print *, "Error = ", C
    enddo
   print *, "velocity = ", v1
   print *, "It took ", i , "interations to finish"
   
   endif
   if (CORREL == 2) then
    !                        print *, "Correl == 2 OK"
                            call cd_sphere(Re, Cd)
     !                       print *, "CD = ", CD
                            Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, CD, v1)
      !                      print *, "Velocity in Stokes = ", v1
                            Re1 = DIAMETER * v1 * WATER_DENSITY/water_dynamic_viscosity
1                            print *, "New Reynolds obtained = ", Re1
                             print *, v1
                            call Error(v0, v1, C)
 !                           print *, "Tolerance = ", C
                            j = 0
                            do while (C > 0.01)
                                 j = j + 1
  !                                      print *, "Seeking Tolerance"    
   !                                     print *, "=================================================="
!                                        print *, "New Reynolds", Re1, "now will be old Reynlods"
                                        Re = Re1
                                        Re1 = 0
                                        v0 = v1
                                        v1 = 0
 !                                       print *, "+++++++++++++++++++++"
  !                                      print *, "Old reynolds = ", Re
                                        Re1 = DIAMETER * v0 * WATER_DENSITY/water_dynamic_viscosity
   !                                     print *, "New Reynolds = ", Re1
    !                                    print *, "Old CD = ", CD
                                        RRey = calculateCd(Re1, Cd)
!                                        print *, "New Cd =",  CD
 !                                       print *, "+++++++++++++++++++++"
                                        Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, Cd, v1)
  !                                      print *, "Velocity 0 = ", v0
   !                                     print *, "Velocity 1 = ", v1
                                        call Error(v0, v1, C)
    !                                    print *, "Tolerance = ", C
     !                                   print *, "Attempt number", j
                                         
                            enddo
                            
                          print *, "velocity = ", v1
                          print *, "It took ", j , "interations to finish"
    
   endif
   if (CORREL == 3) then
     !  print *, "CORREL == 3"
                            call ReadReynolds(Re, CD, ReynoldsFilePath)
      !                      print *, "CD = " ,CD
                            Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, CD, v1)
       !                     print *, "Velocity in Stokes = ", v1
                            Re1 = DIAMETER * v1 * WATER_DENSITY/water_dynamic_viscosity
        !                    print *, "New Reynolds obtained = ", Re1
                            call Error(v0, v1, C)
         !                   print *, "Error = ", C
                            j = 0
                            do while (C > 0.01)
                                 j = j + 1
!                                        print *, "Seeking Tolerance"    
 !                                       print *, "=================================================="
  !                                      print *, "New Reynolds", Re1, "now will be old Reynlods"
                                        Re = Re1
                                        Re1 = 0
                                        v0 = v1
                                        v1 = 0
   !                                     print *, "+++++++++++++++++++++"
    !                                    print *, "Old reynolds = ", Re
                                        Re1 = DIAMETER * v0 * WATER_DENSITY/water_dynamic_viscosity
     !                                   print *, "New Reynolds = ", Re1
      !                                  print *, "Old CD = ", CD
                                        RRey = calculateCd(Re1, Cd)
       !                                 print *, "New Cd =",  CD
        !                                print *, "+++++++++++++++++++++"
                                      !  print *, POSCHECK 
                                        Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, Cd, v1)
         !                               print *, "Velocity 0 = ", v0
          !                              print *, "Velocity 1 = ", v1
                                        call Error(v0, v1, C)
                                        print *, "Error = ", C
                                        print *, "Attempt number", j
                                         
                            enddo
                        
                        endif
    endif
    
            
!---------------------------------------------------------------------------------------------------------
! Option 2: Z > 0
!---------------------------------------------------------------------------------------------------------
    if (Z > 0) then
        print *, "-----------------------------------------------------"
        print *, "Calculating instantaneous velocity Z > 0"
        print *, "-----------------------------------------------------"
        
!        print *, "Positions to be used as steps: ",POS
   !     v1 = 0.0
   !     v0 = 0.0

        LP = size(POS)
        allocate(VelProperty(LP))
       ! print *, "Allocation PosProperty ok"
        print *, VelProperty
        print *, LP
 
            k = 0
            l = 1
            ! CALCULAR A VELOCIDADE EM STOKES
                
                do while (l <= LP) 
        !         print *, "Calculating at Z = ", Z
        !        print *, "Maximum of ", LP
                !print *, 'POSPROPERTY', PosProperty(i)
                !S = size(SAL)
                !D = size(DEN)

    
                    if (size(DENS) /= 0) then
                     WATER_DENSITY = DENS(l)
                    elseif (size(SAL) /= 0) then
                     call DensWater(SAL(l), TEMP(l), WATER_DENSITY)
                    endif
                    
         !           print *, "Parameters :", gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, water_dynamic_viscosity
                     call CalcVelStokes(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, water_dynamic_viscosity, v0) 
          !                  print *, "Attempt ", k, "Stokes Velocity calculated: ", v0
           !                 print *, "Position ", Z
            !                print *, "Calculate Reynolds for Stokes velocity"
                            Re = DIAMETER * v0 * WATER_DENSITY/water_dynamic_viscosity
                            
                    ! CONSIDERAR MORRISSON
                            if (Correl == 1) then
             !                       print *, "-----------------------------------------------------------------------------"
!                                    print *, "Process of calculation for the next position begun . It may take few minute
            !                        print *, "-----------------------------------------------------------------------------"
 !                                   print *, "Calculating Cd using Morrisson"
                                    RRey = calculateCd(Re, Cd)
  !                                  print *, "CD = ", CD                                    
                                    Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, CD, v1)
   !                                 print *, "Calculating new velocity using the CD"
    !                                print *, "Velocity 1 = ", v1
                                    Re1 = DIAMETER * v1 * WATER_DENSITY/water_dynamic_viscosity
      !                              print *, "Calculating New Reynolds using V1 =", v1
     !                               print *, "Reynolds = ", Re1
                                    call Error(v0, v1, C)
       !                             print *, "Error = ", C
                                    j = 0
                                        do while (C > 1)
                                            j = j + 1
        !                                    print *, "=================================================="
         !                                   print *, "                 Seeking Tolerance                 "    
          !                                  print *, "=================================================="
           !                                 print *, "New Reynolds", Re1, "now will be old Reynlods"
                                            Re = Re1
                                            Re1 = 0
                                            v0 = v1
                                            v1 = 0
            !                                print *, "+++++++++++++++++++++"
             !                               print *, "Old reynolds = ", Re
                                            Re1 = DIAMETER * v0 * WATER_DENSITY/water_dynamic_viscosity
              !                              print *, "New Reynolds = ", Re1
               !                             print *, "Old CD = ", CD
                                            RRey = calculateCd(Re1, Cd)
                !                            print *, "New Cd =",  CD
                 !                           print *, "+++++++++++++++++++++"
                                            !print *, "CHECK", PS
                                            Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, Cd, v1)
                  !                          print *, "Velocity 0 = ", v0
                   !                         print *, "Velocity 1 = ", v1
                                            call Error(v0, v1, C)
                                            print *, "velocity = ", v1
                                            print *, "Error = ", C
                                            print *, "Attempt number", j
                                        enddo
                                Z = POS(l)
                                !print *, "Next stop ", Z
                                VelProperty(l) = v1
                            endif
                            if (CORREL == 2) then
                      !          print *, "Correl == 2 OK"
                                call cd_sphere(Re, Cd)
                       !         print *, "CD = ", CD
                                Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, CD, v1)
                        !        print *, "Velocity in Stokes = ", v1
                                Re1 = DIAMETER * v1 * WATER_DENSITY/water_dynamic_viscosity
                         !       print *, "New Reynolds obtained = ", Re1
                                call Error(v0, v1, C)
                          !      print *, "Error = ", C
                                j = 0
                                do while (C > 0.01)
                                     j = j + 1
               !                             print *, "Seeking Tolerance"    
                !                            print *, "=================================================="
                 !                           print *, "New Reynolds", Re1, "now will be old Reynlods"
                                            Re = Re1
                                            Re1 = 0
                                            v0 = v1
                                            v1 = 0
                  !                          print *, "+++++++++++++++++++++"
                   !                         print *, "Old reynolds = ", Re
                                            Re1 = DIAMETER * v0 * WATER_DENSITY/water_dynamic_viscosity
                    !                        print *, "New Reynolds = ", Re1
                     !                       print *, "Old CD = ", CD
                                            RRey = calculateCd(Re1, Cd)
                      !                      print *, "New Cd =",  CD
                       !                     print *, "+++++++++++++++++++++"
                                            Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, Cd, v1)
                        !                    print *, "Velocity 0 = ", v0
                         !                   print *, "Velocity 1 = ", v1
                                            call Error(v0, v1, C)
                                            print *, "velocity =", v1
                                            print *, "Error = ", C
                                            print *, "Attempt number", j
                                        enddo
                                Z = POS(l)
                                print *, POS(l)
                                !print *, "Next stop ", Z
                                VelProperty(l) = v1
                            endif
                            if (CORREL == 3) then
  !                              print *, "CORREL == 3"
                                call ReadReynolds(Re, CD, ReynoldsFilePath)
   !                             print *, "CD = " ,CD
                                Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, CD, v1)
    !                            print *, "Velocity in Stokes = ", v1
                                Re1 = DIAMETER * v0 * WATER_DENSITY/water_dynamic_viscosity
!                                print *, "New Reynolds obtained = ", Re1
                                call Error(v0, v1, C)
 !                               print *, "Error = ", C
                                j = 0
                                
                                do while (C > 0.01) 
                                         j = j + 1
  !                                          print *, "Seeking Tolerance"    
   !                                         print *, "=================================================="
    !                                        print *, "New Reynolds", Re1, "now will be old Reynlods"
                                            Re = Re1
                                            Re1 = 0
                                            v0 = v1
                                            v1 = 0
 !                                           print *, "+++++++++++++++++++++"
  !                                          print *, "Old reynolds = ", Re
                                            Re1 = DIAMETER * v0 * WATER_DENSITY/water_dynamic_viscosity
   !                                         print *, "New Reynolds = ", Re1
    !                                        print *, "Old CD = ", CD
                                            RRey = calculateCd(Re1, Cd)
     !                                       print *, "New Cd =",  CD
      !                                      print *, "+++++++++++++++++++++"
                                          !  print *, POSCHECK 
                                            Rv1 = CalcVelCD(gravitational_constant, DIAMETER, MP_DENSITY, WATER_DENSITY, Cd, v1)
       !                                     print *, "Velocity 0 = ", v0
        !                                    print *, "Velocity 1 = ", v1
                                            call Error(v0, v1, C)
                                            print *, "velocity = ", v1
                                            print *, "Error = ", C  
                                            print *, "Attempt number", j
                                enddo
                                Z = POS(l)
                                print *, POS(l)
                                !print *, "Next stop ", Z
                                Velproperty(l) = v1
                            endif
                     l = l + 1                       
                    enddo
                
    endif
   
!---------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------
!    endif
     !print *, "Saving"
    !---------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------
  !  if (Z > 0) then

end subroutine Calculator




end program Settling

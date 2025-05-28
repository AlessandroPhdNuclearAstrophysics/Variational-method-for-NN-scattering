PROGRAM EVALUATE_LOW_ENERGY_SCATTERING_OBSERVABLES
  IMPLICIT NONE
  
  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=10) :: channel1, channel2
  INTEGER :: L, L1, L2
  INTEGER :: ios, i, npoints
  LOGICAL :: is_coupled
  REAL, ALLOCATABLE :: K2(:), KCOTD(:)
  
  ! Get filename from command line
  IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
    WRITE(*,*) "Error: Please provide input filename"
    STOP
  END IF
  
  CALL GET_COMMAND_ARGUMENT(1, filename)
  
  ! Parse filename to determine channel type and extract L values
  CALL parse_filename(filename, is_coupled, channel1, channel2, L, L1, L2)
  
  ! Read data from file
  CALL read_data(filename, K2, KCOTD, npoints)
  
  ! Display results
  IF (is_coupled) THEN
    WRITE(*,*) "Coupled channel: ", TRIM(channel1), "-", TRIM(channel2)
    WRITE(*,*) "L1 = ", L1, ", L2 = ", L2
  ELSE
    WRITE(*,*) "Single channel: ", TRIM(channel1)
    WRITE(*,*) "L = ", L
  END IF
  
  WRITE(*,*) "Read ", npoints, " data points"
  WRITE(*,*) "K2 range: ", K2(1), " to ", K2(npoints)
  WRITE(*,*) "KCOTD range: ", KCOTD(1), " to ", KCOTD(npoints)
  
  ! Clean up
  DEALLOCATE(K2, KCOTD)
  
CONTAINS

  SUBROUTINE parse_filename(fname, is_coup, ch1, ch2, ang_mom, ang_mom1, ang_mom2)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    LOGICAL, INTENT(OUT) :: is_coup
    CHARACTER(LEN=*), INTENT(OUT) :: ch1, ch2
    INTEGER, INTENT(OUT) :: ang_mom, ang_mom1, ang_mom2
    
    INTEGER :: i, j, k
    CHARACTER(LEN=1) :: first_letter
    
    ! Check if it's a coupled channel (contains "epsilon")
    is_coup = INDEX(fname, "epsilon") > 0
    
    IF (is_coup) THEN
      ! Format: k2_kcotd_epsilon_XXX-YYY.dat
      i = INDEX(fname, "epsilon_") + 9
      j = INDEX(fname(i:), "-") + i - 1
      k = INDEX(fname(j:), ".") + j - 2
      
      ch1 = fname(i:j)
      ch2 = fname(j+2:k+1)
      
      ! Convert spectroscopic notation to L values
      first_letter = ch1(1:1)
      SELECT CASE(first_letter)
        CASE('S', 's'); ang_mom1 = 0
        CASE('P', 'p'); ang_mom1 = 1
        CASE('D', 'd'); ang_mom1 = 2
        CASE('F', 'f'); ang_mom1 = 3
        CASE('G', 'g'); ang_mom1 = 4
        CASE('H', 'h'); ang_mom1 = 5
        CASE DEFAULT; ang_mom1 = -1  ! Error
      END SELECT
      
      first_letter = ch2(1:1)
      SELECT CASE(first_letter)
        CASE('S', 's'); ang_mom2 = 0
        CASE('P', 'p'); ang_mom2 = 1
        CASE('D', 'd'); ang_mom2 = 2
        CASE('F', 'f'); ang_mom2 = 3
        CASE('G', 'g'); ang_mom2 = 4
        CASE('H', 'h'); ang_mom2 = 5
        CASE DEFAULT; ang_mom2 = -1  ! Error
      END SELECT
      
      ang_mom = -1  ! Not applicable for coupled channels
    ELSE
      ! Format: k2_kcotd_XXX.dat
      i = INDEX(fname, "kcotd_") + 6
      j = INDEX(fname(i:), ".") + i - 2
      
      ch1 = fname(i:j)
      ch2 = ""
      
      ! Convert spectroscopic notation to L value
      first_letter = ch1(2:2)
      SELECT CASE(first_letter)
        CASE('S', 's'); ang_mom = 0
        CASE('P', 'p'); ang_mom = 1
        CASE('D', 'd'); ang_mom = 2
        CASE('F', 'f'); ang_mom = 3
        CASE('G', 'g'); ang_mom = 4
        CASE('H', 'h'); ang_mom = 5
        CASE DEFAULT; ang_mom = -1  ! Error
      END SELECT
      
      ang_mom1 = -1  ! Not applicable for single channels
      ang_mom2 = -1  ! Not applicable for single channels
    END IF
  END SUBROUTINE parse_filename
  
  SUBROUTINE read_data(fname, k2_values, kcotd_values, num_points)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    REAL, ALLOCATABLE, INTENT(OUT) :: k2_values(:), kcotd_values(:)
    INTEGER, INTENT(OUT) :: num_points
    
    INTEGER :: unit, io_stat, i, temp_size
    REAL :: k2_temp, kcotd_temp
    REAL, ALLOCATABLE :: temp_k2(:), temp_kcotd(:)
    
    ! First count the number of data points
    OPEN(NEWUNIT=unit, FILE=fname, STATUS='old', ACTION='read', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
      WRITE(*,*) "Error opening file: ", TRIM(fname)
      STOP
    END IF
    
    num_points = 0
    READ(unit, *, IOSTAT=io_stat) 
    DO
      READ(unit, *, IOSTAT=io_stat) k2_temp, kcotd_temp
      IF (io_stat /= 0) EXIT
      num_points = num_points + 1
    END DO
    
    ! Rewind and read the data
    REWIND(unit)
    
    ALLOCATE(k2_values(num_points), kcotd_values(num_points))
    
    READ(unit, *, IOSTAT=io_stat) 
    DO i = 1, num_points
      READ(unit, *) k2_values(i), kcotd_values(i)
    END DO
    
    CLOSE(unit)
  END SUBROUTINE read_data
  
END PROGRAM EVALUATE_LOW_ENERGY_SCATTERING_OBSERVABLES
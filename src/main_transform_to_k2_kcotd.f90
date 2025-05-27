PROGRAM TRANSFORM_FROM_PHASE_SHIFTS_TO_KCOTD
  USE SCATTERING_NN_VARIATIONAL
  USE OPERATING_SYSTEM_LINUX
  USE REALLOCATE_UTILS
  USE QUANTUM_NUMBERS
  IMPLICIT NONE
  INTEGER :: J, L, S, TZ
  
  ! Variables for file handling
  CHARACTER(LEN=255) :: folder_path, file_path, output_file
  CHARACTER(LEN=255) :: file_list(1000), current_file
  INTEGER :: num_files, i, io_status, unit_in, unit_out
  
  ! Variables for data processing
  DOUBLE PRECISION :: energy, delta1, delta2, epsilon, k, k2, kcotd1, kcotd2
  DOUBLE PRECISION :: HTM, pi
  DOUBLE PRECISION, ALLOCATABLE :: energies(:), deltas1(:), deltas2(:), epsilons(:)
  DOUBLE PRECISION, ALLOCATABLE :: k_vals(:), k2_vals(:), kcotd1_vals(:), kcotd2_vals(:)
  INTEGER :: num_lines, j_line, num_cols
  LOGICAL :: is_coupled
  
  ! Constants
  pi = 4.0D0 * ATAN(1.0D0)
  
  ! Get command line argument for folder path
  IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
    WRITE(*,*) "Usage: ./main_transform_to_k2_kcotd folder_path"
    STOP
  END IF
  CALL GET_COMMAND_ARGUMENT(1, folder_path)
  
  ! Initialize scattering module to get HTM
  ! Default parameters for np (TZ=0) system
  TZ = 0
  
  ! Set variational parameters to get HTM
  CALL SET_VARIATIONAL_PARAMETERS(0, 0, 0, TZ, 18)
  HTM = GET_HTM()
  WRITE(*,*) "HTM = ", HTM
  
  ! Get list of files in the folder using the module's function
  CALL LIST_FILES_IN_DIRECTORY(TRIM(folder_path), file_list, num_files)
  WRITE(*,*) "Found ", num_files, " files in folder"
  
  ! Process each file
  DO i = 1, num_files
    current_file = TRIM(file_list(i))
    
    ! Check if it's a phase shift file
    IF (INDEX(current_file, "delta_") /= 1) CYCLE
    
    ! Parse filename to extract L, S, J
    CALL PARSE_FILENAME(TRIM(current_file), S, L, J)
    WRITE(*,*) "Processing file: ", TRIM(current_file), " with L=", L, " S=", S, " J=", J
    
    ! Construct full file path
    file_path = TRIM(folder_path) // "/" // TRIM(current_file)
    
    ! Count number of lines and columns in the file
    CALL COUNT_LINES_AND_COLUMNS(TRIM(file_path), num_lines, num_cols)
    
    ! Determine if this is a coupled channel
    is_coupled = (num_cols > 2)
    
    ! Allocate arrays using REALLOCATE_UTILS
    CALL REALLOCATE_1D_1(energies, num_lines)
    CALL REALLOCATE_1D_1(deltas1, num_lines)
    CALL REALLOCATE_1D_3(k_vals, k2_vals, kcotd1_vals, num_lines)
    
    IF (is_coupled) THEN
      CALL REALLOCATE_1D_2(deltas2, epsilons, num_lines)
      CALL REALLOCATE_1D_1(kcotd2_vals, num_lines)
    END IF
    
    ! Read data from file
    OPEN(NEWUNIT=unit_in, FILE=TRIM(file_path), STATUS='OLD', ACTION='READ')
    
    DO j_line = 1, num_lines
      IF (is_coupled) THEN
        READ(unit_in, *, IOSTAT=io_status) energies(j_line), deltas1(j_line), deltas2(j_line), epsilons(j_line)
      ELSE
        READ(unit_in, *, IOSTAT=io_status) energies(j_line), deltas1(j_line)
      END IF
      
      IF (io_status /= 0) EXIT
      
      ! Calculate k, k^2, and k^(2L+1)*cot(delta)
      k = SQRT(energies(j_line)/HTM)
      k2 = k*k
      
      ! Convert degrees to radians for calculations
      deltas1(j_line) = deltas1(j_line) * pi / 180.0D0
      
      ! Calculate k^(2L+1)*cot(delta)
      kcotd1 = k**(2*L+1) * COS(deltas1(j_line)) / SIN(deltas1(j_line))
      
      k_vals(j_line) = k
      k2_vals(j_line) = k2
      kcotd1_vals(j_line) = kcotd1
      
      IF (is_coupled) THEN
        deltas2(j_line) = deltas2(j_line) * pi / 180.0D0
        kcotd2 = k**(2*L+1) * COS(deltas2(j_line)) / SIN(deltas2(j_line))
        kcotd2_vals(j_line) = kcotd2
      END IF
    END DO
    
    CLOSE(unit_in)
    
    ! Create output file name
    output_file = TRIM(folder_path) // "/k2_kcotd_" // TRIM(current_file(7:))
    
    ! Write results to output file
    OPEN(NEWUNIT=unit_out, FILE=TRIM(output_file), STATUS='REPLACE', ACTION='WRITE')
    
    ! Write header
    IF (is_coupled) THEN
      WRITE(unit_out, '(A)') "# k^2    k^(2L+1)*cot(delta1)    k^(2L+1)*cot(delta2)    epsilon"
    ELSE
      WRITE(unit_out, '(A)') "# k^2    k^(2L+1)*cot(delta)"
    END IF
    
    ! Write data
    DO j_line = 1, num_lines
      IF (is_coupled) THEN
        WRITE(unit_out, '(4E18.10)') k2_vals(j_line), kcotd1_vals(j_line), &
                                    kcotd2_vals(j_line), epsilons(j_line)
      ELSE
        WRITE(unit_out, '(2E18.10)') k2_vals(j_line), kcotd1_vals(j_line)
      END IF
    END DO
    
    CLOSE(unit_out)
    WRITE(*,*) "Created output file: ", TRIM(output_file)
    
    ! Arrays will be deallocated automatically at next allocation or program end
  END DO
  
  WRITE(*,*) "All files processed successfully"

CONTAINS

  SUBROUTINE PARSE_FILENAME(filename, S_out, L_out, J_out)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: S_out, L_out, J_out
    
    CHARACTER(LEN=1) :: L_char
    INTEGER :: s_val, pos_start, pos_end
    
    ! Extract the part after "delta_"
    pos_start = INDEX(filename, "delta_") + 6
    pos_end = INDEX(filename, ".dat") - 1
    
    IF (pos_start <= 6 .OR. pos_end <= 0) THEN
      WRITE(*,*) "Error parsing filename: ", TRIM(filename)
      S_out = 0
      L_out = 0
      J_out = 0
      RETURN
    END IF
    
    ! Parse 2S+1 from first character
    READ(filename(pos_start:pos_start), '(I1)') s_val
    S_out = (s_val - 1) / 2
    
    ! Parse L from second character
    L_char = filename(pos_start+1:pos_start+1)
    SELECT CASE(L_char)
      CASE('S')
        L_out = 0
      CASE('P')
        L_out = 1
      CASE('D')
        L_out = 2
      CASE('F')
        L_out = 3
      CASE('G')
        L_out = 4
      CASE DEFAULT
        WRITE(*,*) "Unknown angular momentum: ", L_char
        L_out = -1
    END SELECT
    
    ! Parse J from last character
    READ(filename(pos_start+2:pos_end), '(I1)') J_out
  END SUBROUTINE PARSE_FILENAME
  
  SUBROUTINE COUNT_LINES_AND_COLUMNS(file_path, num_lines, num_cols)
    CHARACTER(LEN=*), INTENT(IN) :: file_path
    INTEGER, INTENT(OUT) :: num_lines, num_cols
    
    INTEGER :: unit_num, io_stat
    CHARACTER(LEN=1000) :: line
    CHARACTER(LEN=20) :: temp
    INTEGER :: i, pos
    
    num_lines = 0
    num_cols = 0
    
    OPEN(NEWUNIT=unit_num, FILE=TRIM(file_path), STATUS='OLD', ACTION='READ')
    
    ! Read first line to determine number of columns
    READ(unit_num, '(A)', IOSTAT=io_stat) line
    IF (io_stat == 0) THEN
      line = ADJUSTL(line)
      IF (line(1:1) /= '#') THEN
        ! Count columns in first line
        i = 1
        pos = 1
        DO WHILE (pos <= LEN_TRIM(line))
          READ(line(pos:), *, IOSTAT=io_stat) temp
          IF (io_stat /= 0) EXIT
          num_cols = num_cols + 1
          
          ! Find position after this number
          DO WHILE (pos <= LEN_TRIM(line) .AND. line(pos:pos) /= ' ')
            pos = pos + 1
          END DO
          
          ! Skip spaces
          DO WHILE (pos <= LEN_TRIM(line) .AND. line(pos:pos) == ' ')
            pos = pos + 1
          END DO
        END DO
        
        num_lines = 1
      END IF
    END IF
    
    ! Count remaining lines
    DO
      READ(unit_num, *, IOSTAT=io_stat)
      IF (io_stat /= 0) EXIT
      num_lines = num_lines + 1
    END DO
    
    CLOSE(unit_num)
  END SUBROUTINE COUNT_LINES_AND_COLUMNS

END PROGRAM TRANSFORM_FROM_PHASE_SHIFTS_TO_KCOTD
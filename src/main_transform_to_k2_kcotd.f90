!> \file main_transform_to_k2_kcotd.f90
!! \brief Transform phase shift files to k^2 and k^{2L+1}cot(delta) representation.
!!
!! This program reads phase shift files (delta_XXX.dat) for nucleon-nucleon scattering channels,
!! computes the corresponding k^2 and k^{2L+1}cot(delta) values, and writes them to output files.
!! For coupled channels, it also processes the second phase shift and mixing angle (epsilon).
!!
!! \details
!! - Input: Folder containing phase shift files named as delta_XXX.dat.
!! - Output: Files with k^2 and kcotd for each channel, and epsilon for coupled channels.
!! - Uses modules: SCATTERING_NN_VARIATIONAL, OPERATING_SYSTEM_LINUX, REALLOCATE_UTILS, QUANTUM_NUMBERS.
!!
!! \author Alessandro
!! \date 2025

PROGRAM TRANSFORM_FROM_PHASE_SHIFTS_TO_KCOTD
  USE SCATTERING_NN_VARIATIONAL
  USE OPERATING_SYSTEM_LINUX
  USE REALLOCATE_UTILS
  USE QUANTUM_NUMBERS
  IMPLICIT NONE
  INTEGER :: J, L, S, TZ

  ! Variables for file handling
  CHARACTER(LEN=255) :: folder_path, file_path, output_file, output_file1, output_file2
  CHARACTER(LEN=255) :: file_list(1000), current_file
  INTEGER :: num_files, i, io_status, unit_in = 32412, unit_out = 23421

  ! Variables for data processing
  DOUBLE PRECISION :: k, k2, kcotd1, kcotd2
  DOUBLE PRECISION :: HTM, PI
  DOUBLE PRECISION, ALLOCATABLE :: energies(:), deltas1(:), deltas2(:), epsilons(:)
  DOUBLE PRECISION, ALLOCATABLE :: k_vals(:), k2_vals(:), kcotd1_vals(:), kcotd2_vals(:)
  INTEGER :: num_lines, j_line, NCH
  LOGICAL :: is_coupled
  TYPE(SCATTERING_CHANNEL) :: TMP

  ! Constants
  pi = 4.0D0 * DATAN(1.0D0)

  ! Get command line argument for folder path
  IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
    WRITE(*,*) "Usage: ./main_transform_to_k2_kcotd folder_path"
    STOP
  END IF
  CALL GET_COMMAND_ARGUMENT(1, folder_path)
  ! Remove trailing '/' from folder_path if present
  IF (LEN_TRIM(folder_path) > 0) THEN
    IF (folder_path(LEN_TRIM(folder_path):LEN_TRIM(folder_path)) == '/') THEN
      folder_path = folder_path(:LEN_TRIM(folder_path)-1)
    END IF
  END IF

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

  !> \brief Main loop: process each phase shift file in the folder.
  !> - Detects single/coupled channels.
  !> - Reads energies and phase shifts.
  !> - Computes k, k^2, and k^{2L+1}cot(delta).
  !> - Writes results to output files.
  ! Process each file
  DO i = 1, num_files
    current_file = TRIM(file_list(i))

    ! Check if it's a phase shift file
    IF (INDEX(current_file, "delta_") /= 1) CYCLE

    ! Parse filename to extract L, S, J
    IF (LEN_TRIM(PARSE_FILENAME(TRIM(current_file))) > 7) CYCLE
    TMP = GET_CHANNEL_FROM_NAME(PARSE_FILENAME(TRIM(current_file)))
    NCH = GET_CHANNEL_NCH(TMP)

    WRITE(*,*) "Processing file: ", TRIM(current_file), " with"
    L = GET_CHANNEL_L(TMP, 1)
    S = GET_CHANNEL_S(TMP, 1)
    J = GET_CHANNEL_J(TMP)
    WRITE(*,*) "L=", L, " S=", S, " J=", J
    IF (NCH > 1) THEN
      L = GET_CHANNEL_L(TMP, 2)
      S = GET_CHANNEL_S(TMP, 2)
      J = GET_CHANNEL_J(TMP)
      WRITE(*,*) "L=", L, " S=", S, " J=", J
    END IF

    ! Construct full file path
    file_path = TRIM(folder_path) // "/" // TRIM(current_file)

    ! Determine if this is a coupled channel
    IF (IS_CHANNEL_COUPLED(TMP)) THEN
      is_coupled = .TRUE.
      WRITE(*,*) "Coupled channel detected"
    ELSE
      is_coupled = .FALSE.
      WRITE(*,*) "Single channel detected"
    END IF

    ! Open the file to count the number of data lines
    OPEN(NEWUNIT=unit_in, FILE=TRIM(file_path), STATUS='OLD', ACTION='READ')
    num_lines = 0
    DO
      READ(unit_in, *, IOSTAT=io_status)
      IF (io_status /= 0) EXIT
      num_lines = num_lines + 1
    END DO
    REWIND(unit_in)

    ! Allocate arrays using REALLOCATE_UTILS
    CALL REALLOCATE_1D_1(energies, num_lines)
    CALL REALLOCATE_1D_1(deltas1, num_lines)
    CALL REALLOCATE_1D_3(k_vals, k2_vals, kcotd1_vals, num_lines)

    IF (is_coupled) THEN
      CALL REALLOCATE_1D_2(deltas2, epsilons, num_lines)
      CALL REALLOCATE_1D_1(kcotd2_vals, num_lines)
    END IF

    ! Read data from file
    ! OPEN(NEWUNIT=unit_in, FILE=TRIM(file_path), STATUS='OLD', ACTION='READ')

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
      L = GET_CHANNEL_L(TMP, 1)
      kcotd1 = k**(2*L+1) * COS(deltas1(j_line)) / SIN(deltas1(j_line))

      k_vals(j_line) = k
      k2_vals(j_line) = k2
      kcotd1_vals(j_line) = kcotd1

      IF (is_coupled) THEN
        L = GET_CHANNEL_L(TMP, 2)
        deltas2(j_line) = deltas2(j_line) * pi / 180.0D0
        if(deltas2(j_line) /= 0.D0) then
          kcotd2 = k**(2*L+1) * COS(deltas2(j_line)) / SIN(deltas2(j_line))
        else
          kcotd2 = HUGE(1.0D0)  ! Handle zero phase shift case
        end if
        kcotd2_vals(j_line) = kcotd2
      END IF
    END DO

    CLOSE(unit_in)

    ! Write results to output file
    if (is_coupled) THEN
      output_file = TRIM(folder_path) // "/k2_kcotd_" // TRIM(current_file(7:9)) // ".dat"
      output_file1= TRIM(folder_path) // "/k2_kcotd_" // TRIM(current_file(11:13)) // ".dat"
      output_file2= TRIM(folder_path) // "/k2_kcotd_epsilon_" // TRIM(current_file(7:))
    ELSE
      output_file = TRIM(folder_path) // "/k2_kcotd_" // TRIM(current_file(7:9)) // ".dat"
    END IF

    ! Write header
    OPEN(unit_out, FILE=TRIM(output_file), STATUS='REPLACE', ACTION='WRITE')
    IF (is_coupled) OPEN(unit_out+1, FILE=TRIM(output_file1), STATUS='REPLACE', ACTION='WRITE')
    IF (is_coupled) OPEN(unit_out+2, FILE=TRIM(output_file2), STATUS='REPLACE', ACTION='WRITE')

    WRITE(unit_out, '(A)') "# k^2    k^(2L+1)*cot(delta1)"
    IF (is_coupled) THEN
      WRITE(unit_out+1, '(A)') "# k^2    k^(2L+1)*cot(delta2)"
      WRITE(unit_out+2, '(A)') "# k^2    epsilon"
    END IF

    ! Write data
    DO j_line = 1, num_lines
      IF (is_coupled) THEN
        IF (IS_FINITE(kcotd1_vals(j_line))) WRITE(unit_out  , '(4E30.16)') k2_vals(j_line), kcotd1_vals(j_line)
        IF (IS_FINITE(kcotd2_vals(j_line))) WRITE(unit_out+1, '(4E30.16)') k2_vals(j_line), kcotd2_vals(j_line)
        IF (IS_FINITE(epsilons(j_line)))    WRITE(unit_out+2, '(4E30.16)') k2_vals(j_line), epsilons(j_line)
      ELSE
        IF (IS_FINITE(kcotd1_vals(j_line))) WRITE(unit_out, '(2E30.16)') k2_vals(j_line), kcotd1_vals(j_line)
      END IF
    END DO

    CLOSE(unit_out)
    if (is_coupled) CLOSE(unit_out+1)
    if (is_coupled) CLOSE(unit_out+2)
    ! Output file created
    IF (is_coupled) THEN
      WRITE(*,*) "Created output files: ", TRIM(output_file), " and ", TRIM(output_file1), " and ", TRIM(output_file2)
    ELSE
      WRITE(*,*) "Created output file: ", TRIM(output_file)
    END IF

    ! Arrays will be deallocated automatically at next allocation or program end
  END DO

  WRITE(*,*) "All files processed successfully"

CONTAINS

  !> \brief Check if a floating-point number is finite.
  !! \param[in] NUMBER Value to check
  !! \return .TRUE. if finite, .FALSE. otherwise
  FUNCTION IS_FINITE(NUMBER) RESULT(FINITE)
    DOUBLE PRECISION, INTENT(IN) :: NUMBER
    LOGICAL :: FINITE
    FINITE = (ABS(NUMBER) < HUGE(1.0D0))
  END FUNCTION IS_FINITE

  !> \brief Parse the filename to extract the channel name.
  !! \param[in] FILENAME Input filename (e.g., "delta_XXX.dat")
  !! \return Channel name string (e.g., "XXX")
  FUNCTION PARSE_FILENAME(FILENAME) RESULT(PARSED)
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    CHARACTER(LEN=8) :: PARSED

    INTEGER :: POS_START, POS_END

    ! EXTRACT THE PART AFTER "DELTA_"
    POS_START = INDEX(FILENAME, "delta_") + 6
    POS_END = INDEX(FILENAME, ".dat") - 1

    PARSED = FILENAME(POS_START:POS_END)
  END function PARSE_FILENAME


END PROGRAM TRANSFORM_FROM_PHASE_SHIFTS_TO_KCOTD
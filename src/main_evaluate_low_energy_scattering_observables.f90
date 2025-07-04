!> \file main_evaluate_low_energy_scattering_observables.f90
!! \brief Evaluate low-energy scattering observables from k^2-kcotd data.
!!
!! This program reads k^2 and k^{2L+1}cot(delta) data files, performs polynomial
!! regression fits (up to 5th order), and prints the fit results and predictions.
!! Handles both single and coupled channels (epsilon mixing).
!!
!! \details
!! - Input: k2_kcotd_XXX.dat or k2_kcotd_epsilon_XXX-YYY.dat
!! - Output: Fit coefficients and predicted values for each order.
!! - Uses: FIT_MODULE for regression routines.
!!
!! \author Alessandro
!! \date 2025

PROGRAM EVALUATE_LOW_ENERGY_SCATTERING_OBSERVABLES
  USE FIT_MODULE
  IMPLICIT NONE

  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=10) :: channel1, channel2
  INTEGER :: L, L1, L2
  INTEGER :: i, npoints
  LOGICAL :: is_coupled
  DOUBLE PRECISION, ALLOCATABLE :: K2(:), KCOTD(:)
  DOUBLE PRECISION :: C, B, A, M, Q, Y0
  DOUBLE PRECISION :: COEFF0(1)
  DOUBLE PRECISION :: COEFF1(2)
  DOUBLE PRECISION :: COEFF2(3)
  DOUBLE PRECISION :: COEFF3(4)
  DOUBLE PRECISION :: COEFF4(5)
  DOUBLE PRECISION :: COEFF5(6)

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

  WRITE(*,*) "#Read ", npoints, " data points"
  IF(.NOT.is_coupled) THEN
    WRITE(*,*) "#Single channel: ", TRIM(channel1)
    WRITE(*,*) "#L = ", L
    WRITE(*,*) "#K2 range: ", K2(1), " fm^2 to ", K2(npoints), " fm^2"
    WRITE(*,*) "#K^(2L+1)COTdelta range: ", KCOTD(1), " to ", KCOTD(npoints)

    CALL LINEAR_REGRESSION(KCOTD, K2, 1, 10, M, Q)
    Y0 = Q
    WRITE(*,'(A,F10.6)') "#y =", Y0

    CALL LINEAR_REGRESSION(KCOTD, K2, 1, npoints/10, M, Q)
    WRITE(*,'(A,F10.6,A,F10.6,A)') "#y =", Q, " +", M, "*x"

    CALL QUADRATIC_REGRESSION(KCOTD, K2, 1, npoints, A, B, C)
    WRITE(*,'(A,F10.6,A,F10.6,A,F10.6,A)') "#y =", C, " +", B, "*x +", A, "*x^2"

    CALL POLYNOMIAL_REGRESSION(KCOTD, K2, 3, npoints, coeff3)
    WRITE(*,'(A,F10.6,A,F10.6,A,F10.6,A,F10.6,A)') "#y =", coeff3(1), " +", coeff3(2), "*x +", &
      coeff3(3), "*x^2 +", coeff3(4), "*x^3"

    CALL POLYNOMIAL_REGRESSION(KCOTD, K2, 4, npoints, coeff4)
    WRITE(*,'(A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F13.6,A)') "#y =", coeff4(1), " +", coeff4(2), &
       "*x +", coeff4(3), "*x^2 +", coeff4(4), "*x^3 +", coeff4(5), "*x^4"

    CALL POLYNOMIAL_REGRESSION(KCOTD, K2, 5, npoints, coeff5)
    WRITE(*,'(A,F10.6,A,F10.6,A,F10.6,A,F13.6,A,F13.6,A,F13.6,A)') "#y =", coeff5(1), " +", &
      coeff5(2), "*x +", coeff5(3), "*x^2 +", coeff5(4), "*x^3 +", coeff5(5), "*x^4 +", coeff5(6), "*x^5"

    WRITE(*,*)
    WRITE(*,*) "# K2 (fm^2)  KCOTD  Oth-order  1st-order  2nd-order  3rd-order  4th-order  5th-order"

    DO I=1, npoints
      WRITE(*, *) K2(I), KCOTD(I), &
            Y0, &
            Q + M*K2(I), &
            C + B*K2(I) + A*K2(I)**2, &
            coeff3(1) + coeff3(2)*K2(I) + coeff3(3)*K2(I)**2 + coeff3(4)*K2(I)**3, &
            coeff4(1) + coeff4(2)*K2(I) + coeff4(3)*K2(I)**2 + coeff4(4)*K2(I)**3 + coeff4(5)*K2(I)**4, &
            coeff5(1) + coeff5(2)*K2(I) + coeff5(3)*K2(I)**2 + coeff5(4)*K2(I)**3 + coeff5(5)*K2(I)**4 + coeff5(6)*K2(I)**5
    ENDDO
  ELSE
    WRITE(*,*) "#Coupled channel mixing-angles: ", TRIM(channel1), "-", TRIM(channel2)
    WRITE(*,*) "#L1 = ", L1, ", L2 = ", L2

    WRITE(*,*) "#K2 range: ", K2(1), " fm^2 to ", K2(npoints), " fm^2"
    WRITE(*,*) "#epsilon range: ", KCOTD(1), " to ", KCOTD(npoints)

    Y0 = 0.0D0
    WRITE(*,'(A,F10.6)') "#y =", Y0

    CALL POLYNOMIAL_REGRESSION_NO_CONSTANT(KCOTD, K2, 1, 10, COEFF0)
    WRITE(*,'(A,F10.6,A)') "#y =", COEFF0(1), "*x"

    CALL POLYNOMIAL_REGRESSION_NO_CONSTANT(KCOTD, K2, 2, npoints, COEFF1)
    WRITE(*,'(A,F10.6,A,F10.6,A)') "#y =", COEFF1(1), "*x + ", COEFF1(2), "*x^2"

    CALL POLYNOMIAL_REGRESSION_NO_CONSTANT(KCOTD, K2, 3, npoints, COEFF2)
    WRITE(*,'(A,F10.6,A,F10.6,A,F10.6,A)') "#y =", COEFF2(1), "*x + ", COEFF2(2), "*x^2 + ", COEFF2(3), "*x^3"

    CALL POLYNOMIAL_REGRESSION_NO_CONSTANT(KCOTD, K2, 4, npoints, COEFF3)
    WRITE(*,'(A,F10.6,A,F10.6,A,F10.6,A,F10.6,A)') "#y =", COEFF3(1), "*x + ", COEFF3(2), "*x^2 + ", COEFF3(3), &
     "*x^3 + ", COEFF3(4), "*x^4"

    CALL POLYNOMIAL_REGRESSION_NO_CONSTANT(KCOTD, K2, 5, npoints, COEFF4)
    WRITE(*,'(A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A)') "#y =", COEFF4(1), "*x + ", COEFF4(2), &
     "*x^2 + ", COEFF4(3), "*x^3 + ", COEFF4(4), "*x^4 + ", COEFF4(5), "*x^5"

    WRITE(*,*)
    WRITE(*,*) "# K2 (fm^2)  KCOTD  Oth-order  1st-order  2nd-order  3rd-order  4th-order  5th-order"

    DO I=1, npoints
      WRITE(*, *) K2(I), KCOTD(I), &
            Y0, &
            COEFF0(1)*K2(I), &
            COEFF1(1)*K2(I) + COEFF1(2)*K2(I)**2, &
            COEFF2(1)*K2(I) + COEFF2(2)*K2(I)**2 + COEFF2(3)*K2(I)**3, &
            COEFF3(1)*K2(I) + COEFF3(2)*K2(I)**2 + COEFF3(3)*K2(I)**3 + COEFF3(4)*K2(I)**4, &
            COEFF4(1)*K2(I) + COEFF4(2)*K2(I)**2 + COEFF4(3)*K2(I)**3 + COEFF4(4)*K2(I)**4 + COEFF4(5)*K2(I)**5
    ENDDO
  ENDIF

  ! DO I=NPOINTS/10, npoints, 10
  !   CALL POLYNOMIAL_REGRESSION(KCOTD, K2, 1, i, COEFF1)
  !   CALL POLYNOMIAL_REGRESSION(KCOTD, K2, 2, i, COEFF2)
  !   CALL POLYNOMIAL_REGRESSION(KCOTD, K2, 3, i, COEFF3)
  !   CALL POLYNOMIAL_REGRESSION(KCOTD, K2, 4, i, COEFF4)
  !   WRITE(120,*) K2(I), COEFF1(1), COEFF2(1), COEFF3(1), COEFF4(1)
  !   WRITE(121,*) K2(I), COEFF1(2), COEFF2(2), COEFF3(2), COEFF4(2)
  !   WRITE(122,*) K2(I), COEFF2(3), COEFF3(3), COEFF4(3)
  !   WRITE(123,*) K2(I), COEFF3(4), COEFF4(4)
  ! ENDDO

  ! Clean up
  DEALLOCATE(K2, KCOTD)

CONTAINS

  !> \brief Parse the filename to determine channel type and extract L values.
  !! \param[in]  fname   Input filename
  !! \param[out] is_coup .TRUE. if coupled channel (epsilon), .FALSE. otherwise
  !! \param[out] ch1     First channel name
  !! \param[out] ch2     Second channel name (if coupled)
  !! \param[out] ang_mom, ang_mom1, ang_mom2  Orbital angular momenta
  SUBROUTINE parse_filename(fname, is_coup, ch1, ch2, ang_mom, ang_mom1, ang_mom2)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    LOGICAL, INTENT(OUT) :: is_coup
    CHARACTER(LEN=*), INTENT(OUT) :: ch1, ch2
    INTEGER, INTENT(OUT) :: ang_mom, ang_mom1, ang_mom2

    INTEGER :: idx_i, j, k
    CHARACTER(LEN=1) :: first_letter

    ! Check if it's a coupled channel (contains "epsilon")
    is_coup = INDEX(fname, "epsilon") > 0

    IF (is_coup) THEN
      ! Format: k2_kcotd_epsilon_XXX-YYY.dat
      idx_i = INDEX(fname, "epsilon_") + 9
      j = INDEX(fname(idx_i:), "-") + idx_i - 1
      k = INDEX(fname(j:), ".") + j - 2

      ch1 = fname(idx_i:j)
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
      idx_i = INDEX(fname, "kcotd_") + 6
      j = INDEX(fname(idx_i:), ".") + idx_i - 2

      ch1 = fname(idx_i:j)
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

  !> \brief Read k^2 and kcotd data from file.
  !! \param[in]  fname         Input filename
  !! \param[out] k2_values     Array of k^2 values
  !! \param[out] kcotd_values  Array of kcotd values
  !! \param[out] num_points    Number of data points read
  SUBROUTINE read_data(fname, k2_values, kcotd_values, num_points)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT) :: k2_values(:), kcotd_values(:)
    INTEGER, INTENT(OUT) :: num_points

    INTEGER :: unit, io_stat, j
    DOUBLE PRECISION :: k2_temp, kcotd_temp

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
    DO j = 1, num_points
      READ(unit, *) k2_values(j), kcotd_values(j)
    END DO

    CLOSE(unit)
  END SUBROUTINE read_data

END PROGRAM EVALUATE_LOW_ENERGY_SCATTERING_OBSERVABLES
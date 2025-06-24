!> @brief Computes the double factorial of an integer N.
!>
!> @details
!> The double factorial of N, denoted as N!!, is defined as the product of all the integers from N down to 1 that have the same parity (odd or even) as N.
!> For example, 7!! = 7 * 5 * 3 * 1, 8!! = 8 * 6 * 4 * 2.
!> By definition, (-1)!! = 1 and 0!! = 1.
!> If N < -1, the function stops with an error.
!> \ingroup math
!>
!> @param[in] N Integer input value. Must be greater than or equal to -1.
!> @return DFACT The double factorial of N as an integer.
!> @note Only defined for N >= -1. If N < -1, the function will terminate execution.
FUNCTION DOUBLE_FACTORIAL(N) RESULT(DFACT)
    INTEGER, INTENT(IN) :: N
    INTEGER :: I, DFACT

    IF (N < -1) THEN
        STOP 'DOUBLE_FACTORIAL: N must be greater than or equal to -1'
    ELSE IF (N == -1 .OR. N == 0) THEN
        DFACT = 1
    ELSE
        DFACT = 1
        DO I = N, 1, -2
            DFACT = DFACT * I
        END DO
    END IF
END FUNCTION DOUBLE_FACTORIAL

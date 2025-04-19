PROGRAM TEST_AV18
    IMPLICIT NONE

    ! DECLARE VARIABLES
    DOUBLE PRECISION :: R
    DOUBLE PRECISION, DIMENSION(2, 2) :: POTENTIAL_F90, POTENTIAL_F, PERCENTAGE_DIFF
    INTEGER :: L, S, J, T, T1Z, T2Z, I, LPOT, A, B

    ! INITIALIZE TEST PARAMETERS
    LPOT = 1

    ! LOOP OVER R VALUES FROM 0 TO 90
    DO I = 1, 90
        R = REAL(I, KIND=8)
        ! LOOP OVER COMBINATIONS OF QUANTUM NUMBERS (EXAMPLE: L, S, J, T, T1Z, T2Z)
        DO L = 0, 2
            DO S = 0, 1
                DO J = MAX(0, ABS(L - S)), L + S
                    DO T = 0, 1
                        IF (MOD(L + S + T, 2) == 1) THEN
                            DO T1Z = -T, T, 2
                                DO T2Z = -T, T, 2
                                    ! CALL THE POTENTIAL SUBROUTINES
                                    CALL AV18PW  (LPOT, L, S, J, T, T1Z, T2Z, R, POTENTIAL_F)
                                    CALL AV18PW90(LPOT, L, S, J, T, T1Z, T2Z, R, POTENTIAL_F90)

                                    ! CALCULATE PERCENTAGE DIFFERENCE
                                    ! CALCULATE PERCENTAGE DIFFERENCES FOR EACH ELEMENT
                                    DO A = 1, 2
                                        DO B = 1, 2
                                            IF (POTENTIAL_F90(A, B) /= 0.0D0) THEN
                                                PERCENTAGE_DIFF(A, B) = ABS(POTENTIAL_F90(A, B) - POTENTIAL_F(A, B)) &
                                                    / ABS(POTENTIAL_F90(A, B)) * 100.0D0
                                            ELSE
                                                PERCENTAGE_DIFF(A, B) = 0.0D0
                                            END IF
                                        END DO
                                    END DO

                                    ! STOP IF PERCENTAGE DIFFERENCE EXCEEDS TOLERANCE
                                    IF (ANY(PERCENTAGE_DIFF > 1.D-10)) THEN
                                        PRINT *, "AV18PW90 and AV18PW differ by more than 10^-10% at R =", R
                                        PRINT *, "L =", L, "S =", S, "J =", J, "T =", T, "T1Z =", T1Z, "T2Z =", T2Z
                                        PRINT *, "PERCENTAGE DIFFERENCE:"
                                        DO A = 1, 2
                                            DO B = 1, 2
                                                PRINT *, PERCENTAGE_DIFF(A, B)
                                            END DO
                                        END DO
                                    END IF
                                END DO
                            END DO
                        END IF
                    END DO
                END DO
            END DO
        END DO
    END DO
    PRINT *, "All tests passed. AV18PW and AV18PW90 are consistent within 10^10%."
END PROGRAM TEST_AV18
! FILEPATH: SRC/FIT.F90
PROGRAM FIT
    USE SCATTERING_NN_VARIATIONAL
    USE QUANTUM_NUMBERS
    USE OPERATING_SYSTEM_LINUX
    USE EFT_PLESS
    IMPLICIT NONE
    INTEGER :: NE = 200
    INTEGER :: LMAX = 2
    INTEGER :: JMAX = 2
    DOUBLE PRECISION :: EMAX = 1.0D0
    TYPE(SCATTERING_CHANNEL), ALLOCATABLE :: CHANNELS(:)
    DOUBLE PRECISION, ALLOCATABLE :: ENERGIES(:)
    INTEGER :: I, NCHANNELS
    TYPE(PHASE_SHIFT_RESULT) :: PS
    INTEGER :: J, L, S, TZ, IPOT, ILB, LEMP
    INTEGER :: IE
    DOUBLE PRECISION :: E
    CHARACTER(LEN=256) :: OUTPUT_DIR, CHANNEL_NAME
    NAMELIST /INPUT/ IPOT, ILB, LEMP, EMAX, NE, LMAX, JMAX, OUTPUT_DIR

    IPOT = 18

    IF (COMMAND_ARGUMENT_COUNT() == 1) THEN
        CALL GET_COMMAND_ARGUMENT(1, VALUE=CHANNEL_NAME)
        READ(CHANNEL_NAME, *) IPOT
    END IF


    LEMP = 0
    OUTPUT_DIR = 'tmp_output_library/'
    SELECT CASE(IPOT)
    CASE(18)
        ! Set parameters for the NN potential
        ILB = 1
        LEMP = 0
        OUTPUT_DIR = TRIM(OUTPUT_DIR)//'AV18/'
    CASE(19)
        ! Set parameters for the NN potential
        ILB = 15
        LEMP = 0
        OUTPUT_DIR = TRIM(OUTPUT_DIR)//'EFT_pless_15/'
    CASE DEFAULT
        PRINT *, 'Unknown potential type. Exiting.'
        STOP
    END SELECT

    WRITE(*, INPUT)

    CALL CREATE_DIRECTORY(OUTPUT_DIR)

    ALLOCATE(ENERGIES(NE))
    
    ENERGIES = [(EMAX/NE*I, I=1, NE)]

    CALL PREPARE_CHANNELS(LMAX, JMAX, 0, CHANNELS)
    NCHANNELS = SIZE(CHANNELS)

    CALL SET_ENERGIES(ENERGIES)
    CALL SET_CHANNELS(CHANNELS)
    
    DO I=1, NCHANNELS
        CALL PRINT_SCATTERING_CHANNEL(CHANNELS(I))
        J = GET_CHANNEL_J(CHANNELS(I))
        L = GET_CHANNEL_L(CHANNELS(I), 1)
        S = GET_CHANNEL_S(CHANNELS(I), 1)
        TZ = GET_CHANNEL_TZ(CHANNELS(I))
        CHANNEL_NAME = GET_CHANNEL_NAME_FROM_OBJECT(CHANNELS(I))
        
        OPEN(UNIT=15, FILE=TRIM(OUTPUT_DIR)//"delta_"//TRIM(CHANNEL_NAME)//'.dat', STATUS='REPLACE', FORM='FORMATTED')
        WRITE(15, '(A)') '# Scattering channel: ' // TRIM(CHANNEL_NAME)
        WRITE(15, '(A)') '# J = ' // TRIM(ADJUSTL(TO_STRING(J)))
        WRITE(15, '(A)') '# L = ' // TRIM(ADJUSTL(TO_STRING(L)))
        WRITE(15, '(A)') '# S = ' // TRIM(ADJUSTL(TO_STRING(S)))
        WRITE(15, '(A)') '# TZ = ' // TRIM(ADJUSTL(TO_STRING(TZ)))
        WRITE(15, '(A)') '# IPOT = ' // TRIM(ADJUSTL(TO_STRING(IPOT)))
        WRITE(15, '(A)') '# ILB = ' // TRIM(ADJUSTL(TO_STRING(ILB)))
        WRITE(15, '(A)') '# LEMP = ' // TRIM(ADJUSTL(TO_STRING(LEMP)))
        WRITE(15, '(A)') '# Date and time: ' // TRIM(ADJUSTL(DATE_AND_TIME_STRING()))
        WRITE(15, '(A)') '# Convention: Blatt-Biedenharn'
        WRITE(15, '(A)') '# Energy (MeV), delta_1(deg), delta_2(deg), mixing angle(deg)'
        DO IE=1, NE
            E = ENERGIES(IE)
            CALL NN_SCATTERING_VARIATIONAL(E, J, L, S, TZ, IPOT, ILB, LEMP, PS)
            WRITE(15,*) E, PS%delta1_BB, PS%delta2_BB, PS%epsilon_BB
        ENDDO
        CLOSE(15)
    END DO
    
    CONTAINS
    FUNCTION TO_STRING(NUM_I) RESULT(NUM)
        INTEGER, INTENT(IN) :: NUM_I
        CHARACTER(LEN=32) :: NUM
        WRITE(NUM, '(I0)') NUM_I
    END FUNCTION TO_STRING

    FUNCTION DATE_AND_TIME_STRING() RESULT(str)
        CHARACTER(LEN=60) :: str
        INTEGER :: values(8)
        CALL DATE_AND_TIME(values=values)
        WRITE(str, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
            values(1), values(2), values(3), values(5), values(6), values(7)
    END FUNCTION DATE_AND_TIME_STRING

END PROGRAM FIT
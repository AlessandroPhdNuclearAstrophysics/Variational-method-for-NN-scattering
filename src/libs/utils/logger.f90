!> \file logger.f90
!! \defgroup logger Logger
!! \ingroup utils
!! \brief Logger utilities and LOGGER type for colored, leveled, and file/terminal logging.
!!
!! This module defines the LOGGER type and provides routines for configuring, writing,
!! and formatting log messages with levels and colors. It supports file and terminal output,
!! colored messages, and log level filtering. Type-bound procedures allow flexible logging
!! for scientific and general applications.
!!
!! \author Alessandro
!! \date 2025
MODULE LOG
  IMPLICIT NONE 
  PUBLIC :: LOGGER
  PUBLIC :: LOG_MESSAGE
  PUBLIC :: COLOR_RESET, COLOR_RED, COLOR_YELLOW, COLOR_CYAN, COLOR_GREEN

  PRIVATE
  !> \brief Logger type for managing log output, levels, formatting, and file state.
  !! \ingroup logger
  !! This type stores configuration and state for logging, including log file name, logger name,
  !! maximum log level, color output flag, file output flag, log unit, and file open state.
  !! It provides type-bound procedures for setting configuration, writing messages, and controlling output.
  TYPE :: LOGGER
  PRIVATE
    CHARACTER(LEN=256) :: LOG_FILE
    CHARACTER(LEN=256) :: LOGGER_NAME = ""
    INTEGER :: MAX_LOG_LEVEL = 1 ! 0=ERROR, 1=WARNING, 2=INFO, 3=DEBUG
    LOGICAL :: COLORED_OUTPUT = .TRUE. ! Whether to use colored output (terminal syntax style)
    LOGICAL :: LOG_TO_FILE = .FALSE. 
    INTEGER :: LOG_UNIT
    LOGICAL :: LOG_UNIT_OPEN = .FALSE. ! Whether the log unit is open
  CONTAINS
    PROCEDURE, PUBLIC :: SET_LOG_FILE
    PROCEDURE, PUBLIC :: SET_LOG_MAX_LEVEL
    PROCEDURE, PUBLIC :: OPEN_LOG_FILE
    PROCEDURE, PUBLIC :: CLOSE_LOG_FILE
    PROCEDURE, PUBLIC :: SET_UNIT
    PROCEDURE, PUBLIC :: CHANGE_UNIT
    PROCEDURE, PUBLIC :: SET_LOGGER_NAME
    PROCEDURE, PUBLIC :: LOG_ERR
    PROCEDURE, PUBLIC :: LOG_WARNING
    PROCEDURE, PUBLIC :: LOG_INFO
    PROCEDURE, PUBLIC :: LOG_DEBUG
    PROCEDURE, PUBLIC :: LOG_MESSAGE => LOG_MESSAGE_CLASS
    PROCEDURE, PUBLIC :: SET_COLORED_OUTPUT
  END TYPE LOGGER

  CHARACTER(LEN=5), PARAMETER :: COLOR_RESET  = CHAR(27)//'[0m'
  CHARACTER(LEN=5), PARAMETER :: COLOR_RED    = CHAR(27)//'[31m'
  CHARACTER(LEN=5), PARAMETER :: COLOR_YELLOW = CHAR(27)//'[33m'
  CHARACTER(LEN=5), PARAMETER :: COLOR_CYAN   = CHAR(27)//'[36m'
  CHARACTER(LEN=5), PARAMETER :: COLOR_GREEN  = CHAR(27)//'[32m'

CONTAINS
  !> \brief Set the log unit number for file operations. Validates the range.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to set unit for
  !! \param[in] UNIT Log unit number (0-99)
  SUBROUTINE SET_UNIT(LOG_OBJ, UNIT)
    IMPLICIT NONE 
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    INTEGER, INTENT(IN) :: UNIT
    IF (UNIT < 0 .OR. UNIT > 99) THEN
      STOP "ERROR: Log unit must be between 0 and 99"
    END IF
    LOG_OBJ%LOG_UNIT = UNIT
  END SUBROUTINE SET_UNIT

  !> \brief Set the logger name, used as a prefix in log messages.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to set name for
  !! \param[in] NAME Logger name string
  SUBROUTINE SET_LOGGER_NAME(LOG_OBJ, NAME)
    IMPLICIT NONE 
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    IF (LEN(NAME) > 256) THEN
      STOP "ERROR: Logger name exceeds 256 characters"
    END IF
    LOG_OBJ%LOGGER_NAME = TRIM(NAME)
  END SUBROUTINE SET_LOGGER_NAME
  
  !> \brief Enable or disable colored output for terminal messages.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to set color output for
  !! \param[in] COLORED_OUTPUT Logical flag to enable/disable color
  SUBROUTINE SET_COLORED_OUTPUT(LOG_OBJ, COLORED_OUTPUT)
    IMPLICIT NONE 
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    LOGICAL, INTENT(IN) :: COLORED_OUTPUT
    LOG_OBJ%COLORED_OUTPUT = COLORED_OUTPUT
  END SUBROUTINE SET_COLORED_OUTPUT

  !> \brief Change the log file unit, closing the previous one if open.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to change unit for
  !! \param[in] NEW_UNIT New log unit number
  SUBROUTINE CHANGE_UNIT(LOG_OBJ, NEW_UNIT)
    IMPLICIT NONE 
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    INTEGER, INTENT(IN) :: NEW_UNIT
    IF (LOG_OBJ%LOG_UNIT_OPEN) THEN
      CLOSE(LOG_OBJ%LOG_UNIT)
      LOG_OBJ%LOG_UNIT_OPEN = .FALSE.
    END IF
    CALL SET_UNIT(LOG_OBJ, NEW_UNIT)
  END SUBROUTINE CHANGE_UNIT

  !> \brief Open the log file for writing if not already open.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to open file for
  SUBROUTINE OPEN_LOG_FILE(LOG_OBJ)
    IMPLICIT NONE 
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    INTEGER :: IOSTAT
    IF (.NOT. LOG_OBJ%LOG_UNIT_OPEN) THEN
      OPEN(UNIT=LOG_OBJ%LOG_UNIT, FILE=LOG_OBJ%LOG_FILE, STATUS='UNKNOWN', ACTION='WRITE', FORM='FORMATTED', IOSTAT=IOSTAT)
      IF (IOSTAT /= 0) THEN
        STOP "ERROR: Unable to open log file"
      END IF
      LOG_OBJ%LOG_UNIT_OPEN = .TRUE.
    END IF
  END SUBROUTINE OPEN_LOG_FILE

  !> \brief Close the log file if open and reset LOG_UNIT_OPEN flag.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to close file for
  SUBROUTINE CLOSE_LOG_FILE(LOG_OBJ)
    IMPLICIT NONE 
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    IF (LOG_OBJ%LOG_UNIT_OPEN) THEN
      CLOSE(LOG_OBJ%LOG_UNIT)
      LOG_OBJ%LOG_UNIT_OPEN = .FALSE.
    END IF
  END SUBROUTINE CLOSE_LOG_FILE

  !> \brief Configure the logger to use a file, set unit, file name, log level, and color option.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to configure
  !! \param[in] UNIT Log unit number
  !! \param[in] FILE_NAME Log file name
  !! \param[in] LOG_MAX_LEVEL Optional maximum log level
  !! \param[in] COLORED_OUTPUT Optional flag for colored output
  SUBROUTINE SET_LOG_FILE(LOG_OBJ, UNIT, FILE_NAME, LOG_MAX_LEVEL, COLORED_OUTPUT)
    IMPLICIT NONE 
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    INTEGER, INTENT(IN) :: UNIT
    CHARACTER(LEN=*), INTENT(IN) :: FILE_NAME
    INTEGER, INTENT(IN), OPTIONAL :: LOG_MAX_LEVEL
    LOGICAL, INTENT(IN), OPTIONAL :: COLORED_OUTPUT
    IF (LEN(FILE_NAME) > 255) THEN
      STOP "ERROR: Log file name exceeds 255 characters"
    END IF
    IF (FILE_NAME == "") THEN
      STOP "ERROR: Log file name cannot be empty"
    END IF
    LOG_OBJ%LOG_UNIT = UNIT
    LOG_OBJ%LOG_FILE = TRIM(FILE_NAME)
    IF (PRESENT(LOG_MAX_LEVEL)) THEN
      CALL SET_LOG_MAX_LEVEL(LOG_OBJ, LOG_MAX_LEVEL)
    END IF
    IF (PRESENT(COLORED_OUTPUT)) THEN
      LOG_OBJ%COLORED_OUTPUT = COLORED_OUTPUT
    END IF
    LOG_OBJ%LOG_TO_FILE = .TRUE.
  END SUBROUTINE SET_LOG_FILE

  !> \brief Set the maximum log level for filtering messages (0=ERROR, 1=WARNING, 2=INFO, 3=DEBUG).
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to set level for
  !! \param[in] LOG_MAX_LEVEL Maximum log level
  SUBROUTINE SET_LOG_MAX_LEVEL(LOG_OBJ, LOG_MAX_LEVEL)
    IMPLICIT NONE
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    INTEGER, INTENT(IN) :: LOG_MAX_LEVEL
    IF (LOG_MAX_LEVEL < 0 .OR. LOG_MAX_LEVEL > 3) THEN
      STOP "ERROR: Log max level must be between 0 and 3"
    END IF
    LOG_OBJ%MAX_LOG_LEVEL = LOG_MAX_LEVEL
  END SUBROUTINE SET_LOG_MAX_LEVEL

  !> \brief Write a log message to the specified unit if the level is within the max log level.
  !! \ingroup logger
  !! \param[in] LOG_UNIT Log unit number
  !! \param[in] MESSAGE Message to write
  !! \param[in] LEVEL Log level of the message
  !! \param[in] MAX_LOG_LEVEL Maximum log level allowed
  SUBROUTINE LOG_MESSAGE(LOG_UNIT, MESSAGE, LEVEL, MAX_LOG_LEVEL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LOG_UNIT
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
    INTEGER, INTENT(IN) :: LEVEL
    INTEGER, INTENT(IN) :: MAX_LOG_LEVEL

    IF (LEVEL < 0 .OR. LEVEL > MAX_LOG_LEVEL) THEN
      RETURN ! Ignore messages below the max log level
    END IF

    WRITE(LOG_UNIT,*) TRIM(MESSAGE)
  END SUBROUTINE LOG_MESSAGE

  !> \brief Write a log message with optional color and logger name. Handles file or terminal output.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to write message for
  !! \param[in] MESSAGE Message to write
  !! \param[in] LEVEL Log level of the message
  !! \param[in] COLOR Optional color code
  !! \param[in] NUMBER Optional integer to append to the message
  SUBROUTINE LOG_MESSAGE_CLASS(LOG_OBJ, MESSAGE, LEVEL, COLOR, NUMBER)
    USE STRINGS_UTILS
    IMPLICIT NONE
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
    INTEGER, INTENT(IN) :: LEVEL
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: COLOR
    INTEGER, OPTIONAL, INTENT(IN) :: NUMBER
    CHARACTER(LEN=1024) :: FULL_MESSAGE
    CHARACTER(LEN=2) :: SEPARATOR
    IF (MESSAGE(:1) /= ' ') THEN
      SEPARATOR = ""
    ELSE
      SEPARATOR = " "
    END IF

    IF (PRESENT(COLOR) .AND. LOG_OBJ%COLORED_OUTPUT) THEN
      FULL_MESSAGE = TRIM(COLOR)//TRIM(LOG_OBJ%LOGGER_NAME) // TRIM(SEPARATOR) // TRIM(MESSAGE)//TRIM(COLOR_RESET)
    ELSE
      FULL_MESSAGE = TRIM(LOG_OBJ%LOGGER_NAME) // TRIM(SEPARATOR) // TRIM(MESSAGE)
    END IF

    IF (PRESENT(NUMBER)) THEN
      FULL_MESSAGE = TRIM(FULL_MESSAGE) // TRIM(TO_STRING(NUMBER))
    END IF

    IF (LOG_OBJ%LOG_TO_FILE) THEN
      IF (.NOT. LOG_OBJ%LOG_UNIT_OPEN) THEN
        OPEN(UNIT=LOG_OBJ%LOG_UNIT, FILE=LOG_OBJ%LOG_FILE, STATUS='UNKNOWN', ACTION='WRITE', FORM='FORMATTED')
        LOG_OBJ%LOG_UNIT_OPEN = .TRUE.
      END IF
      IF (LOG_OBJ%COLORED_OUTPUT .AND. PRESENT(COLOR)) THEN
        CALL LOG_MESSAGE_COLOR_IMPL(LOG_OBJ%LOG_UNIT, FULL_MESSAGE, COLOR)
      ELSE
        CALL LOG_MESSAGE_COLOR_IMPL(LOG_OBJ%LOG_UNIT, FULL_MESSAGE, "")
      END IF
    ELSE
      IF (LOG_OBJ%COLORED_OUTPUT .AND. PRESENT(COLOR)) THEN
        WRITE(LOG_OBJ%LOG_UNIT,*) TRIM(COLOR)//TRIM(FULL_MESSAGE)//TRIM(COLOR_RESET)
      ELSE
        WRITE(LOG_OBJ%LOG_UNIT,*) TRIM(FULL_MESSAGE)
      END IF
    END IF
  END SUBROUTINE LOG_MESSAGE_CLASS

  !> \brief Log an error message (level 0) in red, with logger name. Only if max log level allows.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to log error for
  !! \param[in] MESSAGE Error message to log
  !! \param[in] NUMBER Optional integer to append to the message
  SUBROUTINE LOG_ERR(LOG_OBJ, MESSAGE, NUMBER)
    IMPLICIT NONE
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
    INTEGER, INTENT(IN), OPTIONAL :: NUMBER
    CHARACTER(LEN=1024) :: FULL_MESSAGE

    IF (LOG_OBJ%MAX_LOG_LEVEL >= 0) THEN
      CALL LOG_MESSAGE_CLASS(LOG_OBJ, MESSAGE, 0, COLOR = COLOR_RED,  NUMBER = NUMBER)
    END IF
  END SUBROUTINE LOG_ERR

  !> \brief Log a warning message (level 1) in yellow, with logger name. Only if max log level allows.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to log warning for
  !! \param[in] MESSAGE Warning message to log
  !! \param[in] NUMBER Optional integer to append to the message
  SUBROUTINE LOG_WARNING(LOG_OBJ, MESSAGE, NUMBER)
    IMPLICIT NONE
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
    INTEGER, INTENT(IN), OPTIONAL :: NUMBER

    IF (LOG_OBJ%MAX_LOG_LEVEL >= 1) THEN
      CALL LOG_MESSAGE_CLASS(LOG_OBJ, MESSAGE, 1, COLOR = COLOR_YELLOW,  NUMBER = NUMBER)
    END IF
  END SUBROUTINE LOG_WARNING

  !> \brief Log an info message (level 2) in cyan, with logger name. Only if max log level allows.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to log info for
  !! \param[in] MESSAGE Info message to log
  !! \param[in] NUMBER Optional integer to append to the message
  SUBROUTINE LOG_INFO(LOG_OBJ, MESSAGE, NUMBER)
    IMPLICIT NONE
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
    INTEGER, INTENT(IN), OPTIONAL :: NUMBER

    IF (LOG_OBJ%MAX_LOG_LEVEL >= 2) THEN
      CALL LOG_MESSAGE_CLASS(LOG_OBJ, MESSAGE, 2, COLOR = COLOR_CYAN,  NUMBER = NUMBER)
    END IF
  END SUBROUTINE LOG_INFO

  !> \brief Log a debug message (level 3) with logger name. Only if max log level allows.
  !! \ingroup logger
  !! \param[inout] LOG_OBJ Logger object to log debug for
  !! \param[in] MESSAGE Debug message to log
  !! \param[in] NUMBER Optional integer to append to the message
  SUBROUTINE LOG_DEBUG(LOG_OBJ, MESSAGE, NUMBER)
    IMPLICIT NONE
    CLASS(LOGGER), INTENT(INOUT) :: LOG_OBJ
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
    INTEGER, INTENT(IN), OPTIONAL :: NUMBER

    IF (LOG_OBJ%MAX_LOG_LEVEL >= 3) THEN
      CALL LOG_MESSAGE_CLASS(LOG_OBJ, MESSAGE, 3, NUMBER = NUMBER)
    END IF
  END SUBROUTINE LOG_DEBUG

  !> \brief Write a colored message to the specified log unit, using ANSI escape codes.
  !! \ingroup logger
  !! \param[in] LOG_UNIT Log unit number
  !! \param[in] MESSAGE Message to write
  !! \param[in] COLOR Color code string
  SUBROUTINE LOG_MESSAGE_COLOR_IMPL(LOG_UNIT, MESSAGE, COLOR)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LOG_UNIT
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
    CHARACTER(LEN=*), INTENT(IN) :: COLOR
    IF (COLOR /= "") THEN
      WRITE(LOG_UNIT, '(XA)', ADVANCE='NO') TRIM(COLOR)//TRIM(MESSAGE)//TRIM(COLOR_RESET)
    ELSE
      WRITE(LOG_UNIT, '(XA)', ADVANCE='NO') TRIM(MESSAGE)
    END IF
  END SUBROUTINE LOG_MESSAGE_COLOR_IMPL

END MODULE LOG
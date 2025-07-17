!> \file test_logger.f90
!! \brief Minimal test program for logger.f90 module
!! \ingroup logger
PROGRAM TEST_LOGGER
  USE LOG
  IMPLICIT NONE
  TYPE(LOGGER) :: log_obj

  ! Set logger name and max level
  CALL log_obj%SET_LOGGER_NAME("TEST_LOGGER")
  CALL log_obj%SET_LOG_MAX_LEVEL(3)

  ! Test error, warning, info, debug
  CALL log_obj%LOG_ERR("This is an error message")
  CALL log_obj%LOG_WARNING("This is a warning message")
  CALL log_obj%LOG_INFO("This is an info message")
  CALL log_obj%LOG_DEBUG("This is a debug message")

END PROGRAM TEST_LOGGER

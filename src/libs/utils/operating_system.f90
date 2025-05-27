MODULE OPERATING_SYSTEM_LINUX
  IMPLICIT NONE

  PUBLIC :: CREATE_DIRECTORY

CONTAINS

  !> \brief Create a directory if it does not exist.
  SUBROUTINE CREATE_DIRECTORY(OUT_DIR)
    CHARACTER(LEN=*), INTENT(IN) :: OUT_DIR
    CHARACTER(LEN=256) :: command

    command = 'mkdir -p "' // TRIM(OUT_DIR) // '"'
    CALL SYSTEM(command)
  END SUBROUTINE CREATE_DIRECTORY

  !> \brief Check if a file exists.
  LOGICAL FUNCTION FILE_EXISTS(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: ios

    OPEN(UNIT=10, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios == 0) THEN
      FILE_EXISTS = .TRUE.
      CLOSE(10)
    ELSE
      FILE_EXISTS = .FALSE.
    END IF
  END FUNCTION FILE_EXISTS

  !> \brief Remove a file if it exists.
  SUBROUTINE REMOVE_FILE(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=256) :: command

    IF (FILE_EXISTS(filename)) THEN
      command = 'rm -f "' // TRIM(filename) // '"'
      CALL SYSTEM(command)
    END IF
  END SUBROUTINE REMOVE_FILE

  !> \brief List files in a directory. Optionally filter by extension.
  SUBROUTINE LIST_FILES_IN_DIRECTORY(dir_path, files, count, extension)
    CHARACTER(LEN=*), INTENT(IN) :: dir_path
    CHARACTER(LEN=255), INTENT(OUT) :: files(*)
    INTEGER, INTENT(OUT) :: count
    CHARACTER(LEN=*), OPTIONAL :: extension

    CHARACTER(LEN=256) :: command
    INTEGER :: ios, i
    CHARACTER(LEN=255) :: temp_file

    temp_file = 'temp_file_list.txt'
    command = 'ls "' // TRIM(dir_path) // '" > "' // TRIM(temp_file) // '"'
    CALL SYSTEM(command)

    OPEN(UNIT=10, FILE=TRIM(temp_file), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(*,*) "Error reading directory contents"
      count = 0
      RETURN
    END IF

    count = 0
    DO
      READ(10, '(A)', IOSTAT=ios) files(count + 1)
      IF (ios /= 0) EXIT
      IF (PRESENT(extension)) THEN
        IF (INDEX(files(count + 1), TRIM(extension)) == 0) CYCLE
      END IF
      count = count + 1
    END DO

    CLOSE(10)
    command = 'rm -f "' // TRIM(temp_file) // '"'
    CALL SYSTEM(command)
  END SUBROUTINE LIST_FILES_IN_DIRECTORY

  !> \brief Get the current working directory.
  SUBROUTINE GET_CURRENT_DIRECTORY(current_dir)
    CHARACTER(LEN=*), INTENT(OUT) :: current_dir
    CHARACTER(LEN=256) :: command
    INTEGER :: ios

    command = 'pwd'
    OPEN(UNIT=10, FILE='current_dir.txt', STATUS='REPLACE', ACTION='WRITE')
    WRITE(10, '(A)') command
    CLOSE(10)

    CALL SYSTEM(command)
    OPEN(UNIT=10, FILE='current_dir.txt', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios == 0) THEN
      READ(10, '(A)') current_dir
      CLOSE(10)
      CALL REMOVE_FILE('current_dir.txt')
    ELSE
      WRITE(*,*) "Error getting current directory"
      current_dir = ''
    END IF
  END SUBROUTINE GET_CURRENT_DIRECTORY

  !> \brief Change the current working directory.
  SUBROUTINE CHANGE_DIRECTORY(new_dir)
    CHARACTER(LEN=*), INTENT(IN) :: new_dir
    CHARACTER(LEN=256) :: command

    command = 'cd "' // TRIM(new_dir) // '"'
    CALL SYSTEM(command)
  END SUBROUTINE CHANGE_DIRECTORY

  !> \brief Get the operating system name.
  CHARACTER(LEN=256) FUNCTION GET_OS_NAME()
    CHARACTER(LEN=256) :: command
    INTEGER :: ios

    command = 'uname -s'
    OPEN(UNIT=10, FILE='os_name.txt', STATUS='REPLACE', ACTION='WRITE')
    WRITE(10, '(A)') command
    CLOSE(10)

    CALL SYSTEM(command)
    OPEN(UNIT=10, FILE='os_name.txt', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios == 0) THEN
      READ(10, '(A)') GET_OS_NAME
      CLOSE(10)
      CALL REMOVE_FILE('os_name.txt')
    ELSE
      WRITE(*,*) "Error getting OS name"
      GET_OS_NAME = 'Unknown'
    END IF
  END FUNCTION GET_OS_NAME

  

END MODULE OPERATING_SYSTEM_LINUX
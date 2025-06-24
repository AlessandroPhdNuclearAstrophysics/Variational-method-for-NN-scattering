!> \file operating_system.f90
!! \brief Utilities for interacting with the operating system (Linux).
!! \defgroup operating_system Operating System
!! \ingroup utils
!!
!! Provides routines for directory and file management, listing files, and
!! basic string operations, designed for use in scientific Fortran codes.
!!
!! \author Alessandro
!! \date 2025

MODULE OPERATING_SYSTEM_LINUX
  IMPLICIT NONE 

CONTAINS
  !> \brief Create a directory if it does not exist.
  !! \ingroup operating_system
  !! \param[in] OUT_DIR Directory path to create
  SUBROUTINE CREATE_DIRECTORY(OUT_DIR)
    CHARACTER(LEN=*), INTENT(IN) :: OUT_DIR
    CHARACTER(LEN=256) :: command

    command = 'mkdir -p "' // TRIM(OUT_DIR) // '"'
    CALL SYSTEM(command)
  END SUBROUTINE CREATE_DIRECTORY

  !> \brief Delete a directory and its contents.
  !! \ingroup operating_system
  !! \param[in] OUT_DIR Directory path to delete
  !! \note This will remove the directory and all its contents recursively.
  !!       Use with caution as it cannot be undone.
  SUBROUTINE DELETE_DIRECTORY(OUT_DIR)
    CHARACTER(LEN=*), INTENT(IN) :: OUT_DIR
    CHARACTER(LEN=256) :: command
    INTEGER :: ios

    command = 'rm -rvf "' // TRIM(OUT_DIR) // '"'
    CALL SYSTEM(command, STATUS=ios)
    IF (ios /= 0) THEN
      WRITE(*,*) "Error deleting directory"
    END IF
  END SUBROUTINE DELETE_DIRECTORY

  !> \brief Check if a file exists.
  !! \ingroup operating_system
  !! \param[in] filename File path to check
  !! \return .TRUE. if file exists, .FALSE. otherwise
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
  !! \ingroup operating_system
  !! \param[in] filename File path to remove
  SUBROUTINE REMOVE_FILE(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=256) :: command

    IF (FILE_EXISTS(filename)) THEN
      command = 'rm -f "' // TRIM(filename) // '"'
      CALL SYSTEM(command)
    END IF
  END SUBROUTINE REMOVE_FILE

  !> \brief List files in a directory. Optionally filter by extension.
  !! \ingroup operating_system
  !! \param[in]  dir_path  Directory to list
  !! \param[out] files     Array of file names
  !! \param[out] count     Number of files found
  !! \param[in]  extension (optional) Filter by file extension
  SUBROUTINE LIST_FILES_IN_DIRECTORY(dir_path, files, count, extension)
    CHARACTER(LEN=*), INTENT(IN) :: dir_path
    CHARACTER(LEN=255), INTENT(OUT) :: files(*)
    INTEGER, INTENT(OUT) :: count
    CHARACTER(LEN=*), OPTIONAL :: extension

    CHARACTER(LEN=256) :: command
    INTEGER :: ios
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
  !! \ingroup operating_system
  !! \param[out] current_dir Path of current directory
  SUBROUTINE GET_CURRENT_DIRECTORY(current_dir)
    CHARACTER(LEN=*), INTENT(OUT) :: current_dir
    CHARACTER(LEN=256) :: command
    INTEGER :: ios

    command = 'pwd > current_dir.txt'
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

  !> \brief Get the operating system name.
  !! \ingroup operating_system
  !! \return OS name string
  CHARACTER(LEN=256) FUNCTION GET_OS_NAME()
    CHARACTER(LEN=256) :: command
    INTEGER :: ios

    command = 'uname -s > os_name.txt'
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

  !> \brief Find all case-insensitive substring matches in a string.
  !! \ingroup operating_system
  !! \param[in]  input_string String to search
  !! \param[in]  pattern      Pattern to match
  !! \param[out] matches      Array of matches found
  !! \param[out] n_matches    Number of matches found
  SUBROUTINE FIND_STRING_CASE_INSENSITIVE(input_string, pattern, matches, n_matches)
    USE, INTRINSIC :: ISO_C_BINDING
    CHARACTER(LEN=*), INTENT(IN) :: input_string
    CHARACTER(LEN=*), INTENT(IN) :: pattern
    CHARACTER(LEN=*), INTENT(OUT) :: matches(:)
    INTEGER, INTENT(OUT) :: n_matches

    INTEGER :: pos, next_pos, pattern_len, input_len
    CHARACTER(LEN=:), ALLOCATABLE :: lowercase_input, lowercase_pattern

    ! NOTE: This is a simplified implementation that only performs case-insensitive
    ! substring matching, not full regex pattern matching. A real implementation
    ! would require C/C++ regex libraries linked to Fortran.

    n_matches = 0
    pattern_len = LEN_TRIM(pattern)
    input_len = LEN_TRIM(input_string)

    IF (pattern_len == 0 .OR. input_len == 0) RETURN

    ! Convert to lowercase for case-insensitive matching
    ALLOCATE(CHARACTER(LEN=input_len) :: lowercase_input)
    ALLOCATE(CHARACTER(LEN=pattern_len) :: lowercase_pattern)

    lowercase_input = input_string(1:input_len)
    lowercase_pattern = pattern(1:pattern_len)

    CALL TO_LOWERCASE(lowercase_input)
    CALL TO_LOWERCASE(lowercase_pattern)

    ! Find all occurrences
    pos = 1
    DO WHILE (pos <= input_len - pattern_len + 1)
      next_pos = INDEX(lowercase_input(pos:), TRIM(lowercase_pattern))
      IF (next_pos == 0) EXIT

      ! Found a match
      n_matches = n_matches + 1
      IF (n_matches <= SIZE(matches)) THEN
        pos = pos + next_pos - 1
        matches(n_matches) = input_string(pos:pos+pattern_len-1)
        pos = pos + 1
      ELSE
        EXIT  ! Array is full
      END IF
    END DO

    DEALLOCATE(lowercase_input)
    DEALLOCATE(lowercase_pattern)
  END SUBROUTINE FIND_STRING_CASE_INSENSITIVE

  !> \brief Convert a string to lowercase (in-place).
  !! \ingroup operating_system
  !! \param[inout] string String to convert
  SUBROUTINE TO_LOWERCASE(string)
    CHARACTER(LEN=*), INTENT(INOUT) :: string
    INTEGER :: i, char_code

    DO i = 1, LEN_TRIM(string)
      char_code = IACHAR(string(i:i))
      IF (char_code >= IACHAR('A') .AND. char_code <= IACHAR('Z')) THEN
        string(i:i) = ACHAR(char_code + 32)  ! Convert to lowercase
      END IF
    END DO
  END SUBROUTINE TO_LOWERCASE

  !> \brief Find a file by name within a directory (recursively).
  !! \ingroup operating_system
  !! \param[in] FILENAME Name of the file to search for
  !! \param[in] DIRECTORY Directory path to search in
  !! \return FILE_NAME_WITH_RELATIVE_PATH Relative path to the found file, or empty string if not found
  FUNCTION FIND_FILE(FILENAME, DIRECTORY) RESULT(FILE_NAME_WITH_RELATIVE_PATH)
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    CHARACTER(LEN=*), INTENT(IN) :: DIRECTORY
    CHARACTER(LEN=255) :: FILE_NAME_WITH_RELATIVE_PATH
    CHARACTER(LEN=256) :: command
    INTEGER :: ios

    ! Create a temporary file to store the results
    command = 'find "' // TRIM(DIRECTORY) // '" -type f -name "' // TRIM(FILENAME) // '" > temp_file_list.txt'
    CALL SYSTEM(command)

    ! Read the first file found
    OPEN(UNIT=10, FILE='temp_file_list.txt', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios == 0) THEN
      READ(10, '(A)', IOSTAT=ios) FILE_NAME_WITH_RELATIVE_PATH
      CLOSE(10)
      CALL REMOVE_FILE('temp_file_list.txt')
      IF (ios /= 0 .OR. LEN_TRIM(FILE_NAME_WITH_RELATIVE_PATH) == 0) THEN
        FILE_NAME_WITH_RELATIVE_PATH = ''
        WRITE(*,*) "No files found."
      END IF
    ELSE
      WRITE(*,*) "Error finding files in directory"
      FILE_NAME_WITH_RELATIVE_PATH = ''
      WRITE(*,*) "No files found."
    END IF
  END FUNCTION FIND_FILE

  !> \brief Get the current working directory as a string.
  !! \ingroup operating_system
  !! \return current_dir The current working directory path
  FUNCTION GET_CURRENT_WORKING_DIRECTORY() RESULT(current_dir)
    CHARACTER(LEN=256) :: current_dir
    CHARACTER(LEN=256) :: command
    INTEGER :: ios

    ! Get the current working directory
    command = 'pwd > current_dir.txt'
    CALL SYSTEM(command)

    OPEN(UNIT=10, FILE='current_dir.txt', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios == 0) THEN
      READ(10, '(A)') current_dir
      CLOSE(10)
      CALL REMOVE_FILE('current_dir.txt')
      current_dir = TRIM(current_dir)
    ELSE
      WRITE(*,*) "Error getting current directory"
      current_dir = ''
    END IF
  END FUNCTION GET_CURRENT_WORKING_DIRECTORY

END MODULE OPERATING_SYSTEM_LINUX
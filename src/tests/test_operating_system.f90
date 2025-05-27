PROGRAM test_operating_system
  USE OPERATING_SYSTEM_LINUX
  IMPLICIT NONE
  
  CHARACTER(LEN=255) :: current_dir, test_dir, temp_file, os_name
  CHARACTER(LEN=255) :: files(100)
  INTEGER :: i, count, stat
  LOGICAL :: exists
  
  WRITE(*,*) "=== Testing OPERATING_SYSTEM_LINUX Module ==="
  
  ! Test 1: Get OS name
  os_name = GET_OS_NAME()
  WRITE(*,*) "Current OS: ", TRIM(os_name)
  
  ! Test 2: Get current directory
  CALL GET_CURRENT_DIRECTORY(current_dir)
  WRITE(*,*) "Current directory: ", TRIM(current_dir)
  
  ! Test 3: Create a test directory
  test_dir = TRIM(current_dir) // "/test_os_dir"
  CALL CREATE_DIRECTORY(test_dir)
  WRITE(*,*) "Created test directory: ", TRIM(test_dir)
  
  ! Verify directory was created
  exists = FILE_EXISTS(test_dir)
  IF (.NOT. exists) THEN
    WRITE(*,*) "ERROR: Failed to create test directory"
    STOP 1
  END IF
  
  ! Test 4: Create test files
  DO i = 1, 5
    WRITE(temp_file, '(A,A,I1,A)') TRIM(test_dir), "/test_file_", i, ".txt"
    OPEN(UNIT=10, FILE=TRIM(temp_file), STATUS='REPLACE')
    WRITE(10,*) "This is test file ", i
    CLOSE(10)
    WRITE(*,*) "Created test file: ", TRIM(temp_file)
    
    ! Verify file was created
    exists = FILE_EXISTS(temp_file)
    IF (.NOT. exists) THEN
      WRITE(*,*) "ERROR: Failed to create test file ", i
      STOP 1
    END IF
  END DO
  
  ! Test 5: List files in directory
  CALL LIST_FILES_IN_DIRECTORY(test_dir, files, count)
  WRITE(*,*) "Found ", count, " files in test directory:"
  DO i = 1, count
    WRITE(*,*) "  ", TRIM(files(i))
  END DO
  
  ! Test 6: List files with extension filter
  CALL LIST_FILES_IN_DIRECTORY(test_dir, files, count, ".txt")
  WRITE(*,*) "Found ", count, " .txt files in test directory:"
  DO i = 1, count
    WRITE(*,*) "  ", TRIM(files(i))
  END DO
  
  ! Test 7: Remove test files
  DO i = 1, 5
    WRITE(temp_file, '(A,A,I1,A)') TRIM(test_dir), "/test_file_", i, ".txt"
    CALL REMOVE_FILE(temp_file)
    
    ! Verify file was removed
    exists = FILE_EXISTS(temp_file)
    IF (exists) THEN
      WRITE(*,*) "ERROR: Failed to remove test file ", i
      STOP 1
    END IF
    WRITE(*,*) "Removed test file: ", TRIM(temp_file)
  END DO
  
  ! Test 8: Remove test directory (using system command)
  CALL SYSTEM("rmdir " // TRIM(test_dir), stat)
  IF (stat /= 0) THEN
    WRITE(*,*) "ERROR: Failed to remove test directory"
    STOP 1
  END IF
  WRITE(*,*) "Removed test directory"
  
  ! Verify directory was removed
  exists = FILE_EXISTS(test_dir)
  IF (exists) THEN
    WRITE(*,*) "ERROR: Test directory still exists after removal"
    STOP 1
  END IF
  
  ! Test 9: Try changing directory
  CALL CHANGE_DIRECTORY("/tmp")
  CALL GET_CURRENT_DIRECTORY(current_dir)
  WRITE(*,*) "Changed directory to: ", TRIM(current_dir)
  
  ! Test 10: Test regex pattern matching
  WRITE(*,*) "Testing regex pattern matching..."
  
  ! Test 10.1: Simple pattern match
  CALL TEST_REGEX_SIMPLE_MATCH()
  
  ! Test 10.2: Multiple matches
  CALL TEST_REGEX_MULTIPLE_MATCHES()
  
  ! Test 10.3: No matches
  CALL TEST_REGEX_NO_MATCHES()

  WRITE(*,*) "All tests passed successfully!"

CONTAINS

  SUBROUTINE TEST_REGEX_SIMPLE_MATCH()
    CHARACTER(LEN=100) :: input_string
    CHARACTER(LEN=20) :: pattern
    CHARACTER(LEN=50) :: matches(5)
    INTEGER :: n_matches
    
    input_string = "This is a test string."
    pattern = "test"
    
    CALL FIND_STRING_CASE_INSENSITIVE(input_string, pattern, matches, n_matches)
    
    WRITE(*,*) "Simple pattern test: found ", n_matches, " matches for '", TRIM(pattern), "'"
    IF (n_matches > 0) THEN
      DO i = 1, n_matches
        WRITE(*,*) "  Match ", i, ": ", TRIM(matches(i))
      END DO
    END IF
  IF (n_matches /= 1) THEN
      WRITE(*,*) "ERROR: Expected 1 match but found ", n_matches
      STOP 1
    ELSE
      WRITE(*,*) "  Test passed: Found expected match"
    END IF
  END SUBROUTINE TEST_REGEX_SIMPLE_MATCH
  
  SUBROUTINE TEST_REGEX_MULTIPLE_MATCHES()
    CHARACTER(LEN=100) :: input_string
    CHARACTER(LEN=20) :: pattern
    CHARACTER(LEN=50) :: matches(10)
    INTEGER :: n_matches
    
    input_string = "The quick brown fox jumps over the lazy dog."
    pattern = "the"
    
    CALL FIND_STRING_CASE_INSENSITIVE(input_string, pattern, matches, n_matches)
    
    WRITE(*,*) "Multiple pattern test: found ", n_matches, " matches for '", TRIM(pattern), "'"
    IF (n_matches > 0) THEN
      DO i = 1, n_matches
        WRITE(*,*) "  Match ", i, ": ", TRIM(matches(i))
      END DO
    END IF
    IF (n_matches /= 2) THEN
      WRITE(*,*) "ERROR: Expected 2 matches but found ", n_matches
      STOP 1
    ELSE
      WRITE(*,*) "  Test passed: Found expected matches"
    END IF
  END SUBROUTINE TEST_REGEX_MULTIPLE_MATCHES
  
  SUBROUTINE TEST_REGEX_NO_MATCHES()
    CHARACTER(LEN=100) :: input_string
    CHARACTER(LEN=20) :: pattern
    CHARACTER(LEN=50) :: matches(5)
    INTEGER :: n_matches
    
    input_string = "This is a test string."
    pattern = "nonexistent"
    
    CALL FIND_STRING_CASE_INSENSITIVE(input_string, pattern, matches, n_matches)
    
    WRITE(*,*) "No matches test: found ", n_matches, " matches for '", TRIM(pattern), "'"
    IF (n_matches /= 0) THEN
      WRITE(*,*) "ERROR: Expected 0 matches but found ", n_matches
      DO i = 1, n_matches
        WRITE(*,*) "  Match ", i, ": ", TRIM(matches(i))
      END DO
      STOP 1
    ELSE
      WRITE(*,*) "  Test passed: No matches found as expected"
    END IF
  END SUBROUTINE TEST_REGEX_NO_MATCHES
  
END PROGRAM test_operating_system

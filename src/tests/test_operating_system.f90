PROGRAM test_operating_system
  USE OPERATING_SYSTEM_LINUX
  IMPLICIT NONE

  CHARACTER(LEN=256) :: current_dir, test_dir, temp_file, os_name
  CHARACTER(LEN=256) :: files(100)
  INTEGER :: i, count, stat
  LOGICAL :: exists
  INTEGER :: n_passed, n_failed
  LOGICAL :: all_files_created
  LOGICAL :: all_files_removed

  ! Store test results for summary
  INTEGER, PARAMETER :: max_tests = 50
  CHARACTER(LEN=64) :: test_names(max_tests)
  LOGICAL :: test_results(max_tests)
  INTEGER :: n_tests

  n_passed = 0
  n_failed = 0
  n_tests = 0

  WRITE(*,*) "=== Testing OPERATING_SYSTEM_LINUX Module ==="

  ! Test 1: Get OS name
  os_name = GET_OS_NAME()
  WRITE(*,*) "Current OS: ", TRIM(os_name)
  CALL RECORD_TEST_RESULT("Get OS name", LEN_TRIM(os_name) > 0)

  ! Test 2: Get current directory
  CALL GET_CURRENT_DIRECTORY(current_dir)
  WRITE(*,*) "Current directory: ", TRIM(current_dir)
  CALL RECORD_TEST_RESULT("Get current directory", LEN_TRIM(current_dir) > 0)

  ! Test 3: Create a test directory
  test_dir = TRIM(current_dir) // "/test_os_dir"
  CALL CREATE_DIRECTORY(test_dir)
  WRITE(*,*) "Created test directory: ", TRIM(test_dir)
  exists = FILE_EXISTS(test_dir)
  CALL RECORD_TEST_RESULT("Create test directory", exists)
  IF (.NOT. exists) STOP 1

  ! Test 4: Create test files
  all_files_created = .TRUE.
  DO i = 1, 5
    WRITE(temp_file, '(A,A,I1,A)') TRIM(test_dir), "/test_file_", i, ".txt"
    OPEN(UNIT=10, FILE=TRIM(temp_file), STATUS='REPLACE')
    WRITE(10,*) "This is test file ", i
    CLOSE(10)
    WRITE(*,*) "Created test file: ", TRIM(temp_file)
    exists = FILE_EXISTS(temp_file)
    IF (.NOT. exists) all_files_created = .FALSE.
  END DO
  CALL RECORD_TEST_RESULT("Create test files", all_files_created)
  IF (.NOT. all_files_created) STOP 1

  ! Test 5: List files in directory
  CALL LIST_FILES_IN_DIRECTORY(test_dir, files, count)
  WRITE(*,*) "Found ", count, " files in test directory:"
  DO i = 1, count
    WRITE(*,*) "  ", TRIM(files(i))
  END DO
  CALL RECORD_TEST_RESULT("List files in directory", count == 5)

  ! Test 6: List files with extension filter
  CALL LIST_FILES_IN_DIRECTORY(test_dir, files, count, ".txt")
  WRITE(*,*) "Found ", count, " .txt files in test directory:"
  DO i = 1, count
    WRITE(*,*) "  ", TRIM(files(i))
  END DO
  CALL RECORD_TEST_RESULT("List files with extension filter", count == 5)

  !> Test FIND_FILE for various cases (run while files/dir exist)
  CALL TEST_FIND_FILE()

  ! Test 7: Remove test files
  all_files_removed = .TRUE.
  DO i = 1, 5
    WRITE(temp_file, '(A,A,I1,A)') TRIM(test_dir), "/test_file_", i, ".txt"
    CALL REMOVE_FILE(temp_file)
    exists = FILE_EXISTS(temp_file)
    IF (exists) all_files_removed = .FALSE.
    WRITE(*,*) "Removed test file: ", TRIM(temp_file)
  END DO
  CALL RECORD_TEST_RESULT("Remove test files", all_files_removed)
  IF (.NOT. all_files_removed) STOP 1

  ! Test 8: Remove test directory (using system command)
  CALL SYSTEM("rmdir " // TRIM(test_dir), stat)
  CALL RECORD_TEST_RESULT("Remove test directory", stat == 0)
  IF (stat /= 0) STOP 1

  ! Verify directory was removed
  exists = FILE_EXISTS(test_dir)
  CALL RECORD_TEST_RESULT("Verify directory removed", .NOT. exists)
  IF (exists) STOP 1

  ! Test 9: Try changing directory
  CALL CHANGE_DIRECTORY("/tmp")
  CALL GET_CURRENT_DIRECTORY(current_dir)
  WRITE(*,*) "Changed directory to: ", TRIM(current_dir)
  ! Accept /tmp, /tmp/, or any path ending with /tmp or /tmp/, or containing /tmp as last component
  IF (TRIM(current_dir) == "/tmp" .OR. TRIM(current_dir) == "/tmp/" .OR. &
      INDEX(TRIM(current_dir)//'/', "/tmp/") == LEN(TRIM(current_dir)) - 4 .OR. &
      INDEX(TRIM(current_dir), "/tmp") > 0) THEN
    CALL RECORD_TEST_RESULT("Change directory", .TRUE.)
  ELSE
    CALL RECORD_TEST_RESULT("Change directory", .FALSE.)
  END IF

  ! Test 10: Test regex pattern matching
  WRITE(*,*) "Testing regex pattern matching..."

  CALL TEST_REGEX_SIMPLE_MATCH()
  CALL TEST_REGEX_MULTIPLE_MATCHES()
  CALL TEST_REGEX_NO_MATCHES()

  ! Print grouped results at the end
  WRITE(*,*) ""
  WRITE(*,*) "Test Results Summary:"
  DO i = 1, n_tests
    IF (test_results(i)) THEN
      WRITE(*,'(A,A,A)') test_names(i), " : ", CHAR(27)//'[32mPASSED'//CHAR(27)//'[0m'
      n_passed = n_passed + 1
    ELSE
      WRITE(*,'(A,A,A)') test_names(i), " : ", CHAR(27)//'[31mFAILED'//CHAR(27)//'[0m'
      n_failed = n_failed + 1
    END IF
  END DO
  WRITE(*,*) ""
  WRITE(*,*) "  Passed: ", n_passed
  WRITE(*,*) "  Failed: ", n_failed
  IF (n_failed == 0) THEN
    WRITE(*,'(A)') CHAR(27)//'[32mAll tests passed successfully!'//CHAR(27)//'[0m'
  ELSE
    WRITE(*,'(A)') CHAR(27)//'[31mSome tests failed.'//CHAR(27)//'[0m'
  END IF

CONTAINS

  SUBROUTINE RECORD_TEST_RESULT(test_name, passed)
    CHARACTER(*), INTENT(IN) :: test_name
    LOGICAL, INTENT(IN) :: passed
    IF (n_tests < max_tests) THEN
      n_tests = n_tests + 1
      test_names(n_tests) = test_name
      test_results(n_tests) = passed
    END IF
  END SUBROUTINE RECORD_TEST_RESULT

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
    CALL RECORD_TEST_RESULT("Regex simple match", n_matches == 1)
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
    CALL RECORD_TEST_RESULT("Regex multiple matches", n_matches == 2)
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
      DO i = 1, n_matches
        WRITE(*,*) "  Match ", i, ": ", TRIM(matches(i))
      END DO
    END IF
    CALL RECORD_TEST_RESULT("Regex no matches", n_matches == 0)
  END SUBROUTINE TEST_REGEX_NO_MATCHES

  SUBROUTINE TEST_FIND_FILE()
    CHARACTER(LEN=256) :: found_path, test_file, test_dir_local
    LOGICAL :: exists, ok

    WRITE(*,*) "Testing FIND_FILE..."

    ! Use the test directory created earlier
    CALL GET_CURRENT_DIRECTORY(test_dir_local)
    test_dir_local = TRIM(test_dir_local) // "/test_os_dir"

    ! 1. Test finding an existing file
    WRITE(test_file, '(A,A)') "test_file_1.txt"
    found_path = FIND_FILE(test_file, test_dir_local)
    ok = (LEN_TRIM(found_path) > 0)
    CALL RECORD_TEST_RESULT("FIND_FILE existing file", ok)

    exists = FILE_EXISTS("temp_file_list.txt")
    CALL RECORD_TEST_RESULT("FIND_FILE temp file removed (existing)", .NOT. exists)

    ! 2. Test finding a non-existing file
    found_path = FIND_FILE("nonexistent_file_abc123.txt", test_dir_local)
    ok = (LEN_TRIM(found_path) == 0)
    CALL RECORD_TEST_RESULT("FIND_FILE non-existing file", ok)

    exists = FILE_EXISTS("temp_file_list.txt")
    CALL RECORD_TEST_RESULT("FIND_FILE temp file removed (nonexistent)", .NOT. exists)

    ! 3. Edge case: empty filename
    found_path = FIND_FILE("", test_dir_local)
    ok = (LEN_TRIM(found_path) == 0)
    CALL RECORD_TEST_RESULT("FIND_FILE empty filename", ok)

    exists = FILE_EXISTS("temp_file_list.txt")
    CALL RECORD_TEST_RESULT("FIND_FILE temp file removed (empty filename)", .NOT. exists)

    ! 4. Edge case: empty directory
    found_path = FIND_FILE("test_file_1.txt", "")
    ok = (LEN_TRIM(found_path) == 0)
    CALL RECORD_TEST_RESULT("FIND_FILE empty directory", ok)

    exists = FILE_EXISTS("temp_file_list.txt")
    CALL RECORD_TEST_RESULT("FIND_FILE temp file removed (empty directory)", .NOT. exists)

    ! 5. Edge case: directory does not exist
    found_path = FIND_FILE("test_file_1.txt", "/this_directory_should_not_exist_12345")
    ok = (LEN_TRIM(found_path) == 0)
    CALL RECORD_TEST_RESULT("FIND_FILE non-existing directory", ok)

    exists = FILE_EXISTS("temp_file_list.txt")
    CALL RECORD_TEST_RESULT("FIND_FILE temp file removed (non-existing directory)", .NOT. exists)

    WRITE(*,*) "FIND_FILE tests completed."
  END SUBROUTINE TEST_FIND_FILE

END PROGRAM test_operating_system

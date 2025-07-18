!> \file strings.f90
!! \brief String and character manipulation utilities: classification, case conversion, search, split, and replace.
!! \defgroup strings Strings
!! \ingroup utils
!! This module defines utility functions for character classification (e.g., is_upper, is_digit),
!! and string manipulation (e.g., to_upper, replace, split) for use in Fortran programs.
!!
!! \author Alessandro Grassi
!! \date 2025
module strings_utils
  implicit none
  private

  public :: is_upper, is_lower, is_letter, is_digit
  public :: is_alphanumeric, is_whitespace, is_special
  public :: to_upper, to_lower, trim_str
  public :: replace, contains_str, index_of, split
  public :: to_string  ! New function to convert integers, reals, and doubles to strings

  interface to_string
    module procedure to_string_int, to_string_double, to_string_complex
    module procedure to_string_array_double, to_string_array_int, to_string_array_2d
  end interface to_string

contains

  !> @brief Returns true if the character is an uppercase letter (A-Z).
  !> \ingroup strings
  !>
  !> @param[in] c Character to check.
  !> @return true if the character is uppercase, false otherwise.
  logical function is_upper(c)
    character, intent(in) :: c
    is_upper = (c >= 'A' .and. c <= 'Z')
  end function is_upper

  !> @brief Returns true if the character is a lowercase letter (a-z).
  !> \ingroup strings
  !>
  !> @param[in] c Character to check.
  !> @return true if the character is lowercase, false otherwise.
  logical function is_lower(c)
    character, intent(in) :: c
    is_lower = (c >= 'a' .and. c <= 'z')
  end function is_lower

  !> @brief Returns true if the character is a letter (A-Z or a-z).
  !> \ingroup strings
  !>
  !> @param[in] c Character to check.
  !> @return true if the character is a letter, false otherwise.
  logical function is_letter(c)
    character, intent(in) :: c
    is_letter = is_upper(c) .or. is_lower(c)
  end function is_letter

  !> @brief Returns true if the character is a digit (0-9).
  !> \ingroup strings
  !>
  !> @param[in] c Character to check.
  !> @return true if the character is a digit, false otherwise.
  logical function is_digit(c)
    character, intent(in) :: c
    is_digit = (c >= '0' .and. c <= '9')
  end function is_digit

  !> @brief Returns true if the character is a letter or digit.
  !> \ingroup strings
  !>
  !> @param[in] c Character to check.
  !> @return true if the character is alphanumeric, false otherwise.
  logical function is_alphanumeric(c)
    character, intent(in) :: c
    is_alphanumeric = is_letter(c) .or. is_digit(c)
  end function is_alphanumeric

  !> @brief Returns true if the character is whitespace (space, tab, newline, carriage return).
  !> \ingroup strings
  !>
  !> @param[in] c Character to check.
  !> @return true if the character is whitespace, false otherwise.
  logical function is_whitespace(c)
    character, intent(in) :: c
    is_whitespace = (c == ' ' .or. c == '\t' .or. c == '\n' .or. c == '\r')
  end function is_whitespace

  !> @brief Returns true if the character is not alphanumeric or whitespace.
  !> \ingroup strings
  !>
  !> @param[in] c Character to check.
  !> @return true if the character is special (not alphanumeric or whitespace), false otherwise.
  logical function is_special(c)
    character, intent(in) :: c
    is_special = .not. (is_alphanumeric(c) .or. is_whitespace(c))
  end function is_special

  !> @brief Converts a string to uppercase.
  !> \ingroup strings
  !>
  !> @details Each lowercase letter is converted to its uppercase equivalent. The output string has the same length as the input.
  !> @param[in] str Input string.
  !> @return String in uppercase.
  function to_upper(str) result(upper_str)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: upper_str
    integer :: i
    upper_str = str
    do i = 1, len(str)
      if (is_lower(str(i:i))) then
        upper_str(i:i) = achar(iachar(str(i:i)) - 32)
      end if
    end do
  end function to_upper

  !> @brief Converts a string to lowercase.
  !> \ingroup strings
  !>
  !> @details Each uppercase letter is converted to its lowercase equivalent. The output string has the same length as the input.
  !> @param[in] str Input string.
  !> @return String in lowercase.
  function to_lower(str) result(lower_str)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: lower_str
    integer :: i
    lower_str = str
    do i = 1, len(str)
      if (is_upper(str(i:i))) then
        lower_str(i:i) = achar(iachar(str(i:i)) + 32)
      end if
    end do
  end function to_lower

  !> @brief Trims leading and trailing whitespace from a string.
  !> \ingroup strings
  !>
  !> @details Uses Fortran's trim and adjustl to remove whitespace.
  !> @param[in] str Input string.
  !> @return Trimmed string.
  function trim_str(str) result(trimmed)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: trimmed
    trimmed = trim(adjustl(str))
  end function trim_str

  !> @brief Replaces all occurrences of a substring with another substring.
  !> \ingroup strings
  !>
  !> @details The output string is allocated to twice the input length for safety. If the replacement string is longer than the search string, truncation may occur.
  !> @param[in] str Input string.
  !> @param[in] search Substring to search for.
  !> @param[in] replace_with Substring to replace with.
  !> @return String with replacements.
  function replace(str, search, replace_with) result(result_str)
    character(len=*), intent(in) :: str, search, replace_with
    character(len=len(str)*2) :: result_str
    integer :: pos, curr_pos
    result_str = ""
    curr_pos = 1
    do
      pos = index(str(curr_pos:), search)
      if (pos == 0) then
        result_str = trim(result_str) // str(curr_pos:)
        exit
      end if
      result_str = trim(result_str) // str(curr_pos:curr_pos+pos-2) // replace_with
      curr_pos = curr_pos + pos + len(search) - 1
      if (curr_pos > len(str)) exit
    end do
  end function replace

  !> @brief Returns true if a string contains a given substring.
  !> \ingroup strings
  !>
  !> @param[in] str Input string.
  !> @param[in] substr Substring to search for.
  !> @return Logical: .true. if found, false otherwise.
  logical function contains_str(str, substr) result(found)
    character(len=*), intent(in) :: str, substr
    found = index(str, substr) > 0
  end function contains_str

  !> @brief Returns the index of the first occurrence of a substring in a string.
  !> \ingroup strings
  !>
  !> @param[in] str Input string.
  !> @param[in] substr Substring to search for.
  !> @return Index of first occurrence (0 if not found).
  integer function index_of(str, substr) result(pos)
    character(len=*), intent(in) :: str, substr
    pos = index(str, substr)
  end function index_of

  !> @brief Splits a string into an array of substrings using a delimiter.
  !> \ingroup strings
  !>
  !> @details The output array is allocated. Each part is a substring between delimiters. If allocation fails, the function returns immediately.
  !> @param[in] str Input string.
  !> @param[in] delimiter Delimiter string.
  !> @return Allocatable array of substrings.
  function split(str, delimiter) result(parts)
    character(len=*), intent(in) :: str, delimiter
    character(len=len(str)), allocatable :: parts(:)
    integer :: i, count, prev_pos, alloc_stat
    ! Count number of parts
    count = 1
    do i = 1, len_trim(str) - len(delimiter) + 1
      if (str(i:i+len(delimiter)-1) == delimiter) count = count + 1
    end do
    ! Allocate array for parts
    allocate(parts(count), stat=alloc_stat)
    if (alloc_stat /= 0) then
      ! Handle allocation error
      return
    end if
    ! Extract parts
    prev_pos = 1
    count = 0
    do i = 1, len_trim(str) - len(delimiter) + 1
      if (str(i:i+len(delimiter)-1) == delimiter) then
        count = count + 1
        parts(count) = str(prev_pos:i-1)
        prev_pos = i + len(delimiter)
      end if
    end do
    ! Add the last part
    count = count + 1
    parts(count) = str(prev_pos:len_trim(str))
  end function split

  !> @brief Converts an integer to a string.
  !>
  !> @details This function takes an integer input and converts it to a string representation.
  !>
  !> @param[in] num The integer number to convert.
  !> @return The string representation of the integer.
  function to_string_int(num) result(str)
    integer, intent(in) :: num
    character(len=20) :: str
    write(str, *) num  ! Convert integer to string
    str=adjustl(trim(str))
  end function to_string_int

  !> @brief Converts a double precision number to a string.
  !>
  !> @details This function takes a double precision number finput and converts it to a string representation.
  !>
  !> @param[in] num The double precision number to convert.
  !> @return The string representation of the double precision number.
  function to_string_double(num) result(str)
    double precision, intent(in) :: num
    character(len=32) :: str
    write(str, *) num  ! Convert double precision to string
    str=adjustl(trim(str))
  end function to_string_double 


  !> @brief Converts a complex number to a string.
  !>
  !> @details This function takes a complex input and converts it to a string representation.
  !> @param[in] num The complex number to convert.
  !> @return The string representation of the complex number.
  function to_string_complex(num) result(str)
    complex, intent(in) :: num
    character(len=256) :: str
    write(str, *) num  ! Convert complex to string
    str=adjustl(trim(str))
  end function to_string_complex


  !> @brief Converts an array of double precision numbers to an array of strings.
  !>
  !> @details Each element is converted individually using to_string_double.
  !> @param[in] arr Array of double precision numbers to convert.
  !> @return Array of string representations.
  function to_string_array_double(arr) result(str)
    double precision, intent(in), allocatable :: arr(:)
    character(len=32), allocatable :: str(:)
    integer :: i
    allocate(str(size(arr)))
    do i = 1, size(arr)
      write(str(i), *) trim(to_string_double(arr(i)))  ! Convert each element to string
      str(i)=adjustl(trim(str(i)))
    end do
  end function to_string_array_double

  !> @brief Converts an array of integers to an array of strings.
  !>
  !> @details Each element is converted individually using to_string_int.
  !> @param[in] arr Array of integers to convert.
  !> @return Array of string representations.
  function to_string_array_int(arr) result(str)
    integer, intent(in), allocatable :: arr(:)
    character(len=32), allocatable :: str(:)
    integer :: i
    allocate(str(size(arr)))
    do i = 1, size(arr)
      write(str(i), *) trim(to_string_int(arr(i)))  ! Convert each element to string
    end do
  end function to_string_array_int

  !> @brief Converts a 2D array of double precision numbers to a 2D array of strings.
  !>
  !> @details Each element is converted individually using to_string_double.
  !> @param[in] arr 2D array of double precision numbers to convert.
  !> @return 2D array of string representations.
  function to_string_array_2d(arr) result(str)
    double precision, intent(in), allocatable :: arr(:,:)
    character(len=32), allocatable :: str(:,:)
    integer :: i, j
    allocate(str(size(arr, 1), size(arr, 2)))
    do i = 1, size(arr, 1)
      do j = 1, size(arr, 2)
        write(str(i,j), *) trim(to_string_double(arr(i,j)))  ! Convert each element to string
      end do
    end do
  end function to_string_array_2d

end module strings_utils

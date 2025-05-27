module strings_utils
  implicit none
  private

  ! Public functions for character classification
  public :: is_upper, is_lower, is_letter, is_digit
  public :: is_alphanumeric, is_whitespace, is_special
  
  ! Public functions for string manipulation
  public :: to_upper, to_lower, trim_str
  public :: replace, contains_str, index_of, split

contains

  ! Character classification functions
  logical function is_upper(c)
    character, intent(in) :: c
    is_upper = (c >= 'A' .and. c <= 'Z')
  end function is_upper

  logical function is_lower(c)
    character, intent(in) :: c
    is_lower = (c >= 'a' .and. c <= 'z')
  end function is_lower

  logical function is_letter(c)
    character, intent(in) :: c
    is_letter = is_upper(c) .or. is_lower(c)
  end function is_letter

  logical function is_digit(c)
    character, intent(in) :: c
    is_digit = (c >= '0' .and. c <= '9')
  end function is_digit

  logical function is_alphanumeric(c)
    character, intent(in) :: c
    is_alphanumeric = is_letter(c) .or. is_digit(c)
  end function is_alphanumeric

  logical function is_whitespace(c)
    character, intent(in) :: c
    is_whitespace = (c == ' ' .or. c == '\t' .or. c == '\n' .or. c == '\r')
  end function is_whitespace

  logical function is_special(c)
    character, intent(in) :: c
    is_special = .not. (is_alphanumeric(c) .or. is_whitespace(c))
  end function is_special

  ! String manipulation functions
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

  function trim_str(str) result(trimmed)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: trimmed
    
    trimmed = trim(adjustl(str))
  end function trim_str

  function replace(str, search, replace_with) result(result_str)
    character(len=*), intent(in) :: str, search, replace_with
    character(len=len(str)*2) :: result_str
    integer :: i, pos, curr_pos
    
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

  logical function contains_str(str, substr) result(found)
    character(len=*), intent(in) :: str, substr
    
    found = index(str, substr) > 0
  end function contains_str

  integer function index_of(str, substr) result(pos)
    character(len=*), intent(in) :: str, substr
    
    pos = index(str, substr)
  end function index_of

  function split(str, delimiter) result(parts)
    character(len=*), intent(in) :: str, delimiter
    character(len=len(str)), allocatable :: parts(:)
    integer :: i, count, pos, prev_pos, alloc_stat
    
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

end module strings_utils

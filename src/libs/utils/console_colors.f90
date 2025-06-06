module console_colors
  implicit none

  ! Make all public since this is a utility module
  public

  ! ANSI color code constants
  character(len=*), parameter :: COLOR_RED = char(27)//'[31m'
  character(len=*), parameter :: COLOR_GREEN = char(27)//'[32m'
  character(len=*), parameter :: COLOR_YELLOW = char(27)//'[33m'
  character(len=*), parameter :: COLOR_BLUE = char(27)//'[34m'
  character(len=*), parameter :: COLOR_MAGENTA = char(27)//'[35m'
  character(len=*), parameter :: COLOR_CYAN = char(27)//'[36m'
  character(len=*), parameter :: COLOR_WHITE = char(27)//'[37m'
  character(len=*), parameter :: COLOR_RESET = char(27)//'[0m'

contains
  !> Print text in specified color to standard output
  subroutine print_colored(text, color, advance)
    character(len=*), intent(in) :: text
    character(len=*), intent(in) :: color
    logical, intent(in), optional :: advance
    
    logical :: adv
    
    adv = .true.
    if (present(advance)) adv = advance
    
    if (adv) then
      write(*, '(A)') color//trim(text)//COLOR_RESET
    else
      write(*, '(A)', advance='no') color//trim(text)//COLOR_RESET
    end if
  end subroutine print_colored

  !> Get a string wrapped with color codes
  function colored_string(text, color) result(output)
    character(len=*), intent(in) :: text
    character(len=*), intent(in) :: color
    character(len=len_trim(text)+20) :: output
    
    output = color//trim(text)//COLOR_RESET
  end function colored_string

end module console_colors

module file_utils
  implicit none

  public
  
contains

  subroutine get_unused_unit (unit)
    implicit none
    integer, intent (out) :: unit
    character(20) :: read, write, readwrite
    unit = 50
    do
       inquire (unit=unit, read=read, write=write, readwrite=readwrite)
       if (read == "UNKNOWN" .and. write == "UNKNOWN" &
            .and. readwrite == "UNKNOWN") &
            then
          return
       end if
       unit = unit + 1
    end do
    return
  end subroutine get_unused_unit

  subroutine open_old_bin_file(funit,filen,istatus)
    integer, intent(out) :: funit, istatus
    character(*) :: filen

    call get_unused_unit(funit)
    open(unit=funit, file=filen, status='old', position='rewind', iostat=istatus, form='unformatted')
    if (istatus == 0) then
       write(*,*) '*Open Old Bin:', funit, filen
    else
       write(*,*) '*Failed to Open Old Bin:', funit, filen
    end if
    return
  end subroutine open_old_bin_file

  subroutine open_new_bin_file(funit,filen,istatus)
    integer, intent(out) :: funit, istatus
    character(*) :: filen

    call get_unused_unit(funit)
    open(unit=funit, file=filen, status='replace', position='rewind', iostat=istatus, form='unformatted')
    write(*,*) '*Open New Bin:', funit, filen
    return
  end subroutine open_new_bin_file

  subroutine open_new_asc_file(funit,filen,istatus)
    integer, intent(out) :: funit, istatus
    character(*) :: filen

    call get_unused_unit(funit)
    open(unit=funit, file=filen, status='replace', position='rewind', iostat=istatus)
    write(*,*) '*Open New Asc:', funit, filen
    return
  end subroutine open_new_asc_file

  subroutine open_old_asc_file(funit,filen,istatus)
    integer, intent(out) :: funit, istatus
    character(*) :: filen

    call get_unused_unit(funit)
    open(unit=funit, file=filen, status='old', position='rewind', iostat=istatus)
    write(*,*) '*Open Old Asc:', funit, filen, istatus
    
    return
  end subroutine open_old_asc_file
end module file_utils
  
  

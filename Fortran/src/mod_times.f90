!==============================================================================
!> Module which provides access to globaly available named timers
module timer

  implicit none

  type tTimer

     Character(len=256)    :: desc = ""

     Real(kind=8)          :: ut_start = -1._8, ut_end = -1._8
     Integer, Dimension(8) :: rt_start = -1   , rt_end = -1

     Type(tTimer), pointer :: next => null()

  end type tTimer

  Type(tTimer), target , private :: base
  Type(tTimer), pointer, private :: timelist=>null(), end=>null()

contains

  !============================================================================
  Subroutine start_timer(desc,reset)

    Character(len=*) , intent(in)            :: desc
    Logical          , intent(in), optional  :: reset

    Real(kind=8)                   :: cput
    Integer, Dimension(8)          :: realt

    Logical                        :: exist, loc_reset

    !--------------------------------------------------------------------------

    if (present(reset)) then
       loc_reset = reset
    Else
       loc_reset = .TRUE.
    End if

    Call date_and_time(values=realt)
    Call cpu_time(cput)

    exist = .FALSE.

    timelist => base

    Do while (associated(timelist))

       If (trim(desc) == trim(timelist%desc) ) then

          if (loc_reset) then

             exist    = .TRUE.
             timelist%ut_start = cput
             timelist%rt_start = realt

             !write(*,*)"Reseted timer : ",trim(desc)

          else

             exist    = .TRUE.
             
             timelist%ut_start = cput  - (timelist%ut_end-timelist%ut_start)
             timelist%rt_start = realt - diff_realtimes(timelist%rt_end , timelist%rt_start)

             !write(*,*)"Restarted timer : ",trim(desc)

          End if
          
          exit

       End If
       
       timelist=>timelist%next
       
    End Do

    If (.NOT. exist) then
       
       If (.NOT. associated(end)) end => base

       end%desc = desc
       end%ut_start = cput
       end%rt_start = realt
       end%ut_end   = -1._8
       end%rt_end   = -1

       allocate(end%next)
       end => end%next

    end If
       
  end Subroutine start_timer

  !============================================================================
  Subroutine end_timer(desc)

    Character(len=*) , intent(in)  :: desc

    Real(kind=8)                   :: cput
    Integer, Dimension(8)          :: realt

    Logical                        :: exist

    !--------------------------------------------------------------------------

    Call date_and_time(values=realt)
    Call cpu_time(cput)

    exist = .FALSE.

    timelist => base

    Do while (associated(timelist))

       If (trim(desc) == trim(timelist%desc) ) then

          exist    = .TRUE.
          timelist%ut_end   = cput
          timelist%rt_end   = realt
         
          exit

       End If
       
       timelist=>timelist%next
       
    End Do
    
    If (.NOT. exist) then
       Write(*,*)"Timer with desc = ",trim(desc)," was not found in list"
    End If

  End Subroutine end_timer
 
  !============================================================================
  Function get_timer(desc) result(timer)

    Character(len=*) , intent(in)  :: desc

    Type(tTimer), pointer          :: timer

    Integer                        :: ii
    Logical                        :: timer_found

    timer_found = .FALSE.

    timelist => base

    Do while (associated(timelist%next))

       if ( trim(timelist%desc) == trim(desc) ) then
          timer => timelist
          timer_found = .TRUE.
          exit
       End if

       timelist=>timelist%next            

    End Do

    If ( .not. timer_found ) then
       Write(*,*)"!! Warning !! Timer with description "//trim(desc)//" was not found."
       allocate(timer)
    End If

  End Function get_timer

  !============================================================================
  Subroutine write_timelist(timers,exclude,unit)

    character(len=*),Dimension(:), intent(in), optional :: timers
    Logical                      , intent(in), optional :: exclude
    Integer                      , intent(in), optional :: unit
    Integer                                             :: ii, loc_u
    Character(len=4)                                    :: fmt
    integer                                             :: max_dl
    Logical                                             :: loc_exclude

    loc_u = 6
    if (present(unit)) loc_u = unit

    loc_exclude=.false.
    if (present(exclude))loc_exclude=exclude

    timelist => base

    max_dl = len_trim("Timer description")

    Do while (associated(timelist%next))
       if ( len_trim(timelist%desc) > max_dl) &
            max_dl = len_trim(timelist%desc)
       timelist=>timelist%next
    End Do
          
    timelist => base 
    
    write(fmt,"(I0)")max_dl+2
    
    write(loc_u,"("//trim(fmt)//"('-'),"//&
          "'+---------------------------+---------------------------+--------------+ ')")
    write(loc_u,"(1X,'                 ',T"//fmt//","//&
         "' |         User time         |        Elapsed time       | User/Elapsed | ')")
    write(loc_u,"(1X,'Timer description',T"//fmt//","//&
         "' |   [sec]    | hh:mm:ss:ms  |    [sec]   | hh:mm:ss:ms  |      [%]     | ')")
    write(loc_u,"("//trim(fmt)//"('-'),"//&
          "'+------------+--------------+------------+--------------+--------------+ ')")

    !--------------------------------------------------------------------------
    if (present(timers)) then

       if (loc_exclude) then

          list_loop:Do while (associated(timelist%next))

             Do ii = 1, size(timers)
                if ( trim(timelist%desc) == timers(ii) ) then
                   timelist=>timelist%next
                   cycle list_loop
                End if
             End Do
             
             call calc_times(timelist,max_dl,loc_u)
             
             timelist=>timelist%next
             
          End Do list_loop

       else
          Do while (associated(timelist%next))
             Do ii = 1, size(timers)
                
                if ( trim(timelist%desc) == timers(ii) ) then
                   call calc_times(timelist,max_dl,loc_u)                
                End if
                
             End Do
             timelist=>timelist%next
          End Do
       End if
    !--------------------------------------------------------------------------
    Else
       
       Do while (associated(timelist%next))

          call calc_times(timelist,max_dl,loc_u)

          timelist=>timelist%next
       
       End Do

    End if

    write(loc_u,"("//trim(fmt)//"('-'),"//&
         "'+---------------------------+---------------------------+--------------+ ')")
    write(loc_u,*)

  End Subroutine write_timelist

  !============================================================================
  subroutine calc_times(times,max_dl,unit)

    Type(tTimer), intent(in)          :: times
    Integer, intent(in)               :: max_dl
    Integer, intent(in), Optional     :: unit
    Integer                           :: hh, mm, ss, msec, loc_u
    Integer(Kind=8)                   :: msec_s, msec_e, rt_msec_diff, ut_msec_diff
    Character(len=12)                 :: ut_elapsed,rt_elapsed
    Character(len=4)                  :: fmt

    loc_u = 6
    if (present(unit)) loc_u = unit

    write(fmt,'(I0)')max_dl

    msec_s = 0
    msec_s = msec_s + times%rt_start(8) + times%rt_start(7) * 1000 + &
             times%rt_start(6) * 60 * 1000 + times%rt_start(5) * 60 * 60 * 1000

    msec_e = 0
    msec_e = msec_e + times%rt_end(8) + times%rt_end(7) * 1000 + &
             times%rt_end(6) * 60 * 1000 + times%rt_end(5) * 60 * 60 * 1000

    if (msec_e < msec_s) then
       msec_e = msec_e + 24*60*60*1000
    End if

    rt_msec_diff = msec_e - msec_s

    hh = rt_msec_diff/(60 * 60 * 1000)
    mm = (rt_msec_diff - hh * (60 * 60 * 1000)) / (60 * 1000)
    ss = (rt_msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000) / 1000
    msec = (rt_msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000 - ss * 1000)

    write(rt_elapsed,"(I2.2,':',I2.2,':',I2.2,'.',I3.3)") hh,mm,ss,msec         

    ut_msec_diff = INT((times%ut_end-times%ut_start)*1000._8)

    hh = ut_msec_diff/(60 * 60 * 1000)
    mm = (ut_msec_diff - hh * (60 * 60 * 1000)) / (60 * 1000)
    ss = (ut_msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000) / 1000
    msec = (ut_msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000 - ss * 1000)
    
    write(ut_elapsed,"(I2.2,':',I2.2,':',I2.2,'.',I3.3)") hh,mm,ss,msec         

    Write(loc_u,"(1X,A"//fmt//",2(' | ',F10.3,' | ',A12),' | ',F10.3,'   | ')")&
         times%desc,&
         (times%ut_end-times%ut_start),ut_elapsed,&
         real(rt_msec_diff)/1000.,rt_elapsed,&
         ((times%ut_end-times%ut_start)/(real(rt_msec_diff)/1000.))*100.

  End subroutine calc_times

  !============================================================================
  Function diff_realtimes(b,a) result(c)

    Integer, Dimension(8), intent(in) :: a, b
    Integer, Dimension(8)             :: c

    Integer(Kind=8)       :: msec_s, msec_e, rt_msec_diff

    msec_s = 0
    msec_s = msec_s + a(8) + a(7) * 1000 + &
             a(6) * 60 * 1000 + a(5) * 60 * 60 * 1000

    msec_e = 0
    msec_e = msec_e + b(8) + b(7) * 1000 + &
             b(6) * 60 * 1000 + b(5) * 60 * 60 * 1000

    if (msec_e < msec_s) then
       msec_e = msec_e + 24*60*60*1000
    End if

    rt_msec_diff = msec_e - msec_s

    c = 0

    c(5) = rt_msec_diff/(60 * 60 * 1000)
    c(6) = (rt_msec_diff - c(5) * (60 * 60 * 1000)) / (60 * 1000)
    c(7) = (rt_msec_diff - c(5) * (60 * 60 * 1000) - c(6) * 60 * 1000) / 1000
    c(8) = (rt_msec_diff - c(5) * (60 * 60 * 1000) - c(6) * 60 * 1000 - c(7) * 1000)

  End Function diff_realtimes

  !============================================================================
  Function frac_realtime(a,n) result(c)

    Integer, Dimension(8), intent(in) :: a
    Integer              , intent(in) :: n
    Integer, Dimension(8)             :: c

    Integer(Kind=8)       :: msec_s, rt_msec_frac

    msec_s = 0
    msec_s = msec_s + a(8) + a(7) * 1000 + &
             a(6) * 60 * 1000 + a(5) * 60 * 60 * 1000

    rt_msec_frac =  msec_s/n

    c = 0

    c(5) = rt_msec_frac/(60 * 60 * 1000)
    c(6) = (rt_msec_frac - c(5) * (60 * 60 * 1000)) / (60 * 1000)
    c(7) = (rt_msec_frac - c(5) * (60 * 60 * 1000) - c(6) * 60 * 1000) / 1000
    c(8) = (rt_msec_frac - c(5) * (60 * 60 * 1000) - c(6) * 60 * 1000 - c(7) * 1000)

  End Function frac_realtime

  !============================================================================
  subroutine write_realtime(rt,unit)

    Integer, Dimension(8), intent(in)           :: rt
    Integer              , intent(in), optional :: unit

    Integer                                     :: loc_u
    
    loc_u = 6
    if (present(unit)) loc_u = unit
    
    write(loc_u,"(I2.2,':',I2.2,':',I2.2,'.',I3.3)")rt(5:8)

  End subroutine write_realtime
  
end module timer

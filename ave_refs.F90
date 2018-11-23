program main
  integer, parameter :: ave_num = 30, row_num = 1, d = 1396
  real(8) :: x(ave_num, row_num, d, d), ave(row_num, d, d) = 0.0e0
  character filename*128, ofilename*128
  integer :: unit_num, row, j, k, l

! Average corref.XX
  do j = 1, ave_num
    do row = 1, row_num
      unit_num = 30*row_num+j
      write (filename, '("ERmod_", i4.4, "/refs/corref.", i2.2)') j,row
      open (unit = unit_num, file = filename, form = "UNFORMATTED", action ="read", status = "OLD")
      read (unit = unit_num) x(j, row, :, :)
    end do
  end do

  do k = 1, d
    do l = 1, d
      do row = 1, row_num
        do j = 1, ave_num
          ave(row, k, l) = ave(row, k, l) + x(j, row, k, l)
        end do
        ave(row, k, l) = ave(row, k, l) / dble(ave_num)
        open(unit = cor_io, file =
        print *, k, l, ave(row, k, l)
      end do
    end do
  end do

do row = 1, row_num
  unit_num = 13*row_num
  write (filename, '("ave_ERmod01/refs/corref.", i2.2)') row
  do k = 1, d
    do l = 1, d
    open (unit = unit_num, file = ofilename, form = "UNFORMATTED", action ="write", status = "OLD")
    write (unit = unit_num) ave(row, :, :)
    end do
  end do   
end do

! Average weight_refs
 do j = 1, ave_num
    do row = 1, row_num
      unit_num = 30*row_num+j
      write (filename, '("ERmod_", i4.4, "/refs/corref.", i2.2)') j,row
      open (unit = unit_num, file = filename, form = "UNFORMATTED", action ="read", status = "OLD")
      read (unit = unit_num) x(j, row, :, :)
    end do
  end do

  do k = 1, d
    do l = 1, d
      do row = 1, row_num
        do j = 1, ave_num
          ave(row, k, l) = ave(row, k, l) + x(j, row, k, l)
        end do
        ave(row, k, l) = ave(row, k, l) / dble(ave_num)
        open(unit = cor_io, file =
        print *, k, l, ave(row, k, l)
      end do
    end do
  end do

do row = 1, row_num
  unit_num = 13*row_num
  write (filename, '("ave_ERmod01/refs/corref.", i2.2)') row
  do k = 1, d
    do l = 1, d
    open (unit = unit_num, file = ofilename, form = "UNFORMATTED", action ="write", status = "OLD")
    write (unit = unit_num) ave(row, :, :)
    end do
  end do   
end do


end program main

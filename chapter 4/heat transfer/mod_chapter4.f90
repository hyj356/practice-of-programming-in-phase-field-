module chapter4
  use iso_fortran_env, only: WP => real64, stdout => output_unit
  implicit none
  private
  public :: WP, stdout, initial_temp, write_temperature
  
  real(WP), parameter, public, dimension(3) :: three_point_stencil = [1.0_WP, -2.0_WP, 1.0_WP]
  real(WP), parameter, public :: mu = 1.0_WP  !< 如果密度ρ, 比热容cp和热传导系数lamda在仿真过程都保持定制，那么可以使用一个常数来描述它们, 其值等于lamda/(ρ*cp) 1 mm^2 /s
  real(WP), parameter, public :: dt = 0.2_WP  !< 单步步长, 0.2 s
  real(WP), parameter, public :: dx = 1.0 !< 相邻节点的距离为1mm
  integer, parameter, public :: Nx = 128  !< 分成了128个节点

contains
  pure function initial_temp(N, xlo, xhi, target_temp) result(ret)
    !! 分配好内存之后, 将指定区间[xlo, xhi]的温度初始化为目标温度target
    !! 其余区间默认初始温度为0.0
    integer, intent(in) :: N
    integer, intent(in) :: xlo
    integer, intent(in) :: xhi
    real(WP), intent(in) :: target_temp
    real(WP), allocatable, dimension(:) :: ret 
    integer :: i

    allocate(ret(N), source = 0.0_WP)

    do i = xlo, xhi
      ret(i) = target_temp
    end do
  end function initial_temp

  subroutine write_temperature(filename, temp)
    character(len=*), intent(in) :: filename
    real(WP), intent(in) :: temp(:)
    integer :: fileid
    integer :: N
    integer :: i

    N = size(temp)
    open(newunit=fileid, file=filename, action='write')

    do i = 1, N
      write(fileid, *) temp(i)
    end do

    close(fileid)
  end subroutine write_temperature

end module
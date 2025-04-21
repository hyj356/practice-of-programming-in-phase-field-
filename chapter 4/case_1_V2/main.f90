program main
  ! gfortran .\mod_chapter4.f90 .\main.f90 -O2; .\a.exe
  use case1, only: write_con, micro_ch_pre, WP, stdout,&
                   dt, mobility, dx, dy, free_energ_ch => free_energ_ch_v2, &
                   grad_coeff, write_vtk_file, calculate_energy
  implicit none
  integer, parameter :: Nx = 64  !< x方向上有64个节点
  integer, parameter :: Ny = 64  !< y方向上有64个节点
  integer, parameter :: nstep = 10000 !< 总共迭代nstep步
  integer, parameter :: nprint = 1000 !< 每nprint步输出一次结果
  real(WP), parameter :: c0 = 0.4_WP  !< 目标初始平均浓度
  real(WP), allocatable :: con0(:, :) !< 当前时间步的浓度(concentration)数组
  real(WP), allocatable :: dummy(:, :)!< 中间计算的临时数组, 形状与con0一致
  integer :: i, j
  integer :: k

  ! 初始化浓度矩阵
  call micro_ch_pre(Nx, Ny, c0, con0)
  allocate(dummy(Nx, Ny))
  k = 0
  ! 将浓度矩阵导出到vtk文件中以进行可视化
  call write_vtk_file("./VTK/con_", Nx, Ny, con0, k)
  write(stdout, *) k, calculate_energy(Nx, Ny, con0)
  do k = 1, 10000
    ! 首先对浓度数组con0进行laplace算符, 然后计算出中间变量填入数组dummy中
    dummy =   free_energ_ch(con0) - grad_coeff *                 &
             (cshift(con0, 1, dim=1) + cshift(con0, -1, dim=1) + &
              cshift(con0, 1, dim=2) + cshift(con0, -1, dim=2) - &
              4 * con0) / (dx*dy)
    ! 接下来我们需要对dummy进行laplace算符运算
    con0 =    con0 + dt * mobility * &
             (cshift(dummy, 1, dim=1) + cshift(dummy, -1, dim=1) + &
              cshift(dummy, 1, dim=2) + cshift(dummy, -1, dim=2) - &
              4 * dummy) / (dx*dy)
    ! 对一些数值误差进行校正
    do j = 1, Ny
      do i = 1, Nx
        if(con0(i, j) > 0.9999_WP) con0(i, j) = 0.9999_WP
        if(con0(i, j) < 0.0001_WP) con0(i, j) = 0.0001_WP
      end do
    end do
    ! 写出vtk文件
    if(mod(k, nprint) == 0) then
      call write_vtk_file("./VTK/con_", Nx, Ny, con0, k)
      write(stdout, *) k, calculate_energy(Nx, Ny, con0)
    end if
  end do
end program 


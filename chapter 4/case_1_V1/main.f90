program main
  ! gfortran mod_chapter4.f90 main.f90 -O2; ./a.exe
  use case1, only: write_con, micro_ch_pre, WP, stdout,&
                   dt, mobility, dx, dy, free_energ_ch, grad_coeff, &
                   write_vtk_file, calculate_energy
  implicit none
  integer, parameter :: Nx = 64  !< x方向上有64个节点
  integer, parameter :: Ny = 64  !< y方向上有64个节点
  integer, parameter :: nstep = 10000 !< 总共迭代nstep步
  integer, parameter :: nprint = 200 !< 每nprint步输出一次结果
  real(WP), parameter :: c0 = 0.4_WP  !< 目标初始平均浓度
  real(WP), allocatable :: con0(:, :) !< 当前时间步的浓度(concentration)数组
  real(WP), allocatable :: dummy(:, :)!< 中间计算的临时数组, 形状与con0和con1一致
  real(WP) :: lap_con !< 中间计算的临时变量, 其值等于con的梯度
  real(WP), allocatable :: con1(:, :) !< 下一个时间步的浓度(concentration)数组
  real(WP) :: dfdcon  !< f对浓度c的导数
  integer :: i, j
  integer :: k
  integer :: ip, im, jp, jm !< 分别对应i+1, i-1, j+1, j-1


  ! 初始化浓度矩阵
  call micro_ch_pre(Nx, Ny, c0, con0)
  allocate(con1(Nx, Ny), dummy(Nx, Ny))
  k = 0
  ! 将浓度矩阵导出到vtk文件中以进行可视化
  call write_vtk_file("./VTK/con_", Nx, Ny, con0, k)
  write(stdout, *) k, calculate_energy(Nx, Ny, con0)
  do k = 1, 10000
    ! 首先对浓度数组con0进行laplace算符, 然后计算出中间变量填入数组dummy中
    do j = 1, Ny
      do i = 1, Nx
        ip = i + 1; im = i - 1
        jp = j + 1; jm = j - 1
        if(ip == Nx + 1) ip = 1
        if(im == 0) im = Nx
        if(jp == Ny + 1) jp = 1
        if(jm == 0) jm = Ny
        lap_con = (con0(ip, j) + con0(im, j) + con0(i, jp) + con0(i, jm) &
                  - 4.0_WP * con0(i, j)) /(dx*dy)
        dfdcon = free_energ_ch(con0, i, j)
        dummy(i, j) = dfdcon - grad_coeff * lap_con
      end do
    end do
    ! 接下来我们需要对dummy进行laplace算符运算
    do j = 1, Ny
      do i = 1, Nx
        ip = i + 1; im = i - 1
        jp = j + 1; jm = j - 1
        if(ip == Nx + 1) ip = 1
        if(im == 0) im = Nx
        if(jp == Ny + 1) jp = 1
        if(jm == 0) jm = Ny
        ! 先将临时结果存入con1中
        con1(i, j) = (dummy(ip, j) + dummy(im, j) + dummy(i, jp) + dummy(i, jm) &
                    - 4.0_WP * dummy(i, j)) /(dx*dy)
        ! 接下来进行时域积分
        con1(i, j) = con0(i, j) + dt * mobility * con1(i, j)
        ! 对一些数值误差进行校正
        if(con1(i, j) > 0.9999_WP) con1(i, j) = 0.9999_WP
        if(con1(i, j) < 0.0001_WP) con1(i, j) = 0.0001_WP
      end do
    end do
    ! 将con1复制到con0中
    con0 = con1
    ! 写出vtk文件
    if(mod(k, nprint) == 0) then
      call write_vtk_file("./VTK/con_", Nx, Ny, con0, k)
      write(stdout, *) k, calculate_energy(Nx, Ny, con0)
    end if
  end do
end program 


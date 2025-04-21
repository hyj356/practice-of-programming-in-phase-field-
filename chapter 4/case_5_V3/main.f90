program main
  ! gfortran mod_case_5_V3.f90 main.f90 -O2 -fopenmp; ./a.exe
  use case_5_V3, only: write_vtk_file, free_energy_v2, &
                       micro_poly_cell, laplace_2, gradient_mat, &
                       ! 仿真用到的不变量
                       lamda, mu, kisa, WP, stdout
  use omp_lib, only: omp_get_wtime, omp_set_num_threads
  use iso_fortran_env, only: DP => real64
  implicit none
  ! 建模大小的相关参数
  integer, parameter :: Nx = 200
  integer, parameter :: Ny = 200
  integer, parameter :: NxNy = Nx*Ny
  real(WP), parameter :: dx = 1.0_WP
  real(WP), parameter :: dy = 1.0_WP
  real(WP), parameter :: idx = 1.0_WP / dx
  real(WP), parameter :: idy = 1.0_WP / dy
  real(WP), parameter :: idxdy = idx*idy
  ! 仿真时长与步长的相关参数
  integer, parameter :: nstep = 3000
  integer, parameter :: nprint = 50
  real(WP) :: dt = 5.0e-3_WP      !< 注意本次仿真中, 仿真的时间步长不再是一个常量
  ! 决定初始细胞数量和大小的参数
  real(WP), parameter :: R = 12.0_WP  !< ncell为2时, R为25, ncell为80的时候, R为12
  ! 参与仿真计算的不变量
  real(WP), parameter :: gamma_ = 5.0_WP
  real(WP), parameter :: kappa = 60.0_WP
  real(WP), parameter :: pi = 4.0_WP*atan(1.0_WP)
  real(WP), parameter :: const1 = 30.0_WP/(lamda*lamda)
  real(WP), parameter :: const2 = 2.0_WP*mu/(pi*R*R)
  real(WP), parameter :: const3 = 60.0_WP*kappa/(lamda*lamda*kisa)
  real(WP), parameter :: piR2 = pi*R*R
  ! 仿真过程中需要计算的中间变量
  real(WP) :: gamma_cell
  real(WP) :: vinteg    
  real(WP) :: vintegx
  real(WP) :: vintegy
  real(WP) :: sum_phi
  real(WP) :: vnx 
  real(WP) :: vny 
  real(WP), allocatable :: vnphi(:, :)
  real(WP), allocatable :: term2(:, :)
  real(WP), allocatable :: term3(:, :)
  real(WP), allocatable :: lap_phi(:, :)
  real(WP), allocatable :: dfdphi(:, :)
  real(WP) :: vac_cell
  real(DP) :: start_time, end_time
  ! 仿真数据数组及决定其维度大小的变量
  integer :: nccel      !< 仿真模型中的软细胞的数量
  integer :: ncell = 80  !< 模型中细胞的数量, 这个值在当前的模拟中只能为80或者2
  integer, allocatable :: ccell(:)
  real(WP), allocatable, target :: phis(:, :, :)
  real(WP), allocatable :: phi_dx(:, :)
  real(WP), allocatable :: phi_dy(:, :)
  real(WP), allocatable :: phi2(:, :)
  real(WP), allocatable :: vac(:)
  real(WP), pointer :: phi(:, :)
  ! 仿真过程用到的循环变量
  integer :: istep, icell, iccel, jcell

  ! 分配内存
  call micro_poly_cell(Nx, Ny, ncell, R, phis, vac, nccel, ccell)
  write(stdout, "(A, *(I0, 1x))") "Number of soft cell: ", nccel
  write(stdout, "(A, *(I0, 1x))") "Soft cell id: ", ccell
  allocate(phi_dx(Nx, Ny), phi_dy(Nx, Ny), lap_phi(Nx, Ny), vnphi(Nx, Ny))
  allocate(phi2(Nx, Ny), dfdphi(Nx, Ny), term2(Nx, Ny), term3(Nx, Ny))
  ! call write_eta('./initial_phi.txt', sum(phis, dim=3))
  ! 进行并行环境的设置
  if (ncell == 2) then
    call omp_set_num_threads(2)
    write(stdout, *) 2, "of threads are using."
  else 
    call omp_set_num_threads(12)
    write(stdout, *) 12, "of threads are using."
  end if
  ! 记录开始运行的时间
  start_time = omp_get_wtime()
  ! 开始演化过程
  !$omp parallel default(shared) private(gamma_cell, phi, lap_phi, phi_dx, phi_dy, vinteg, vintegx, vintegy) &
  !$omp private(sum_phi, dfdphi, term2, term3, vac_cell, vnx, vny, vnphi, jcell)
  do istep = 1, nstep
    ! 当仿真步数大于500的时候, 单步时间步长变为0.01
    if ( istep >= 500 ) dt = 1.0e-2_WP
    ! 遍历每个cell, 进行时间演化
    !$omp do
    do icell = 1, ncell
      gamma_cell = gamma_
      ! 计算软细胞(soft cell)的弹性参数
      do iccel = 1, nccel
        if (icell == ccell(iccel)) gamma_cell = 0.5_WP*gamma_
      end do
      ! 取出第icell个细胞的序参数数组
      phi => phis(:, :, icell)
      ! 计算拉普拉斯算子项和方向梯度项
      call laplace_2(Nx, Ny, idxdy, phi, lap_phi)
      call gradient_mat(Nx, Ny, dx, dy, phi, phi_dy, phi_dx)
      ! 接下来对每一个节点计算体积积分项, 当然由于当前的仿真
      ! 是在二维空间进行的, 所以实际上是面积分, 而且由于我们网格
      ! 的面积就是1, 所以原代码省去了乘以dA的部分
      vinteg = sum(phi(:, :)*phi(:, :))
      sum_phi = 0.0_WP
      do jcell = 1, ncell
        if (icell /= jcell) then
          sum_phi  = sum_phi + sum(phis(:, :, jcell)*phis(:, :, jcell))
        end if
      end do
      vintegx = sum(phi(:, :)*phi_dx(:, :))*sum_phi
      vintegy = sum(phi(:, :)*phi_dy(:, :))*sum_phi
      ! 计算剩余系数, 逐个节点进行时间积分
      call free_energy_v2(Nx, Ny, dfdphi, icell, ncell, gamma_, kappa, phi, phis)
      ! 计算第二项, 即自由能导数项乘以常数
      term2(:, :) = -const1 * dfdphi(:, :)
      ! 计算第三项
      term3(:, :) = -const2 * (vinteg - piR2) * phi(:, :)
      ! 更新细胞的速度中的常量部分
      vac_cell = merge(0.0, vac(icell), istep <= 200)
      ! 计算细胞速度的2个分量
      vnx = vac_cell + const3*vintegx; vny = vac_cell + const3*vintegy
      ! 当然对于只有两个细胞的时候, 取消竖直方向上的速度
      if (ncell == 2) then
        vnx = vac_cell; vny = 0.0_WP
      end if
      ! 计算第四项, 进行两个向量的点乘
      vnphi(:, :) = vnx*phi_dx(:, :) + vny*phi_dy(:, :)
      ! 进行时间积分
      phi(:, :) = phi(:, :) + dt * (gamma_cell * lap_phi(:, :) + term2(:, :) + term3(:, :) - vnphi(:, :))
      ! 数值误差校正
      where (phi > 0.9999999_WP) 
        phi = 0.9999999_WP
      else where (phi < 0.0_WP) 
        phi = 0.0_WP
      end where
    end do
    !$omp end do
    ! 输出轨迹文件
    !$omp single
    if ( mod(istep, nprint) == 0) then
      write(stdout, *) "Done step: ", istep
      phi2 = sum(phis(:, :, :)*phis(:, :, :), dim=3)
      call write_vtk_file('./VTK/dump_', Nx, Ny, dx, dy, phi2, istep)
    end if
    !$omp end single
  end do
  !$omp end parallel
  ! 记录结束运行的时间
  end_time = omp_get_wtime()
  write(stdout, *) "The program cost ", (end_time-start_time)," seconds in total."
end program main
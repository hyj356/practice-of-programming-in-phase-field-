program main
  ! gfortran mod_case_5_V4.f90 main.f90 -O2 -fopenmp; ./a.exe
  use case_5_V4, only: write_vtk_file, free_energy_v1, &
                       micro_poly_cell, &
                       ! 仿真用到的不变量
                       lamda, mu, kisa, WP, stdout
  use iso_fortran_env, only: DP => real64
  use omp_lib, only: omp_get_wtime, omp_set_num_threads
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
  real(WP), allocatable :: lap_phi(:, :)
  real(WP) :: vinteg    
  real(WP) :: vintegx
  real(WP) :: vintegy
  real(WP) :: sum_phi
  real(WP) :: vnphi
  real(WP) :: term2
  real(WP) :: term3
  real(WP) :: vnx 
  real(WP) :: vny 
  real(WP) :: dfdphi
  real(WP) :: vac_cell
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
  ! openMP并行相关变量
  real(DP) :: start_time, end_time
  ! 仿真过程用到的循环变量
  integer :: istep, i, j, icell, iccel, jcell
  integer :: ip, im, jp, jm

  ! 分配内存
  call micro_poly_cell(Nx, Ny, ncell, R, phis, vac, nccel, ccell)
  write(stdout, "(A, *(I0, 1x))") "Number of soft cell: ", nccel
  write(stdout, "(A, *(I0, 1x))") "Soft cell id: ", ccell
  allocate(phi_dx(Nx, Ny), phi_dy(Nx, Ny), lap_phi(Nx, Ny))
  allocate(phi2(Nx, Ny))
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
  !$omp parallel default(shared) private(gamma_cell, phi, lap_phi, phi_dx, phi_dy) &
  !$omp private(sum_phi, dfdphi, term2, term3, vac_cell, vnx, vny, vnphi)          &
  !$omp private(ip, im, jp, jm, vinteg, vintegx, vintegy, i, j, jcell)
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
      do j = 1, Ny
        do i = 1, Nx
          ip = i + 1; im = i - 1
          jp = j + 1; jm = j - 1
          if(ip == Nx + 1) ip = 1
          if(im == 0) im = Nx
          if(jp == Ny + 1) jp = 1
          if(jm == 0) jm = Ny
          lap_phi(i, j) = (phi(im, j) + phi(ip, j) + phi(i, jm) + phi(i, jp) - &
                           4.0_WP*phi(i, j)) * idxdy
          phi_dx(i, j) = ( phi(ip, j) - phi(im, j) ) * idx
          phi_dy(i, j) = ( phi(i, jp) - phi(i, jm) ) * idy
        end do
      end do
      ! 接下来对每一个节点计算体积积分项, 当然由于当前的仿真
      ! 是在二维空间进行的, 所以实际上是面积分, 而且由于我们网格
      ! 的面积就是1, 所以原代码省去了乘以dA的部分
      vinteg = 0.0_WP; vintegx = 0.0_WP; vintegy = 0.0_WP 
      do j = 1, Ny
        do i = 1, Nx
          vinteg = vinteg + phi(i, j)*phi(i, j) 
          sum_phi = 0.0_WP
          do jcell = 1, ncell
            if (icell /= jcell) then
              sum_phi  = sum_phi + phis(i, j, jcell)*phis(i, j, jcell)
            end if
          end do
          vintegx = vintegx + phi(i, j)*phi_dx(i, j)*sum_phi
          vintegy = vintegy + phi(i, j)*phi_dy(i, j)*sum_phi
        end do
      end do
      ! 计算剩余系数, 逐个节点进行时间积分
      do j = 1, Ny
        do i = 1, Nx
          ! 计算第二项, 即自由能导数项
          dfdphi = free_energy_v1(i, j, icell, ncell, gamma_, kappa, phi, phis)
          term2 = -const1 * dfdphi
          ! 计算第三项
          term3 = -const2 * (vinteg - piR2) * phi(i, j)
          ! 更新细胞的速度中的常量部分
          vac_cell = merge(0.0, vac(icell), istep <= 200)
          ! 计算细胞速度的2个分量
          vnx = vac_cell + const3*vintegx; vny = vac_cell + const3*vintegy
          ! 当然对于只有两个细胞的时候, 取消竖直方向上的速度
          if (ncell == 2) then
            vnx = vac_cell; vny = 0.0_WP
          end if
          ! 计算第四项, 进行两个向量的点乘
          vnphi = vnx*phi_dx(i, j) + vny*phi_dy(i, j)
          ! 进行时间积分
          phi(i, j) = phi(i, j) + dt * (gamma_cell * lap_phi(i, j) + term2 + term3 - vnphi)
        end do
      end do
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
      end_time = omp_get_wtime()
      write(stdout, *) "Done step: ", istep, ". Time used: ", (end_time-start_time), " seconds"
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
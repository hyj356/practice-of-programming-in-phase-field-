program main
  ! gfortran mod_case_4_V2.f90 main.f90 -O2; ./a.exe
  use case_4_V2, only: WP, stdout, gradient_mat, nucleus, laplacian_2, tau, epsilonb, mu, &
  kappa, delta, aniso, alpha, gamma_, Teq, theta0, ipi, itau, write_vtk_file, write_eta,  &
  gradient_matdx, gradient_matdy, laplace_2
  implicit none
  ! 设定模型求解区域的参数
  integer, parameter :: Nx = 300
  integer, parameter :: Ny = 300
  integer, parameter :: NxNy = Nx*Ny
  real(WP), parameter :: dx = 0.03_WP
  real(WP), parameter :: dy = 0.03_WP
  real(WP), parameter :: idx = 1.0_WP / dx  !< dx的倒数
  real(WP), parameter :: idy = 1.0_WP / dy  !< dy的倒数
  real(WP), parameter :: idxdy = 1.0_WP / (dx*dy) !< 其值等于dx*dy的倒数
  real(WP), parameter :: seed = 5.0_WP  !< 初始化成核区域的大小
  ! 设定模型求解时间和步长的相关参数
  integer, parameter :: nstep = 4000
  integer, parameter :: nprint = 100
  real(WP), parameter :: dt = 1.0e-4_WP
  ! 设置需要求解的场变量数组
  real(WP), allocatable :: phi(:, :)      !< 序参数数组, 1表示被固相占据
  real(WP), allocatable :: tempr(:, :)    !< 无量纲温度数组, 通过单位处理使得平衡温度正好为1
  ! 其他求解需要用到的变量
  integer :: istep      !< 循环变量
  real(WP), allocatable :: lap_phi(:, :), lap_tempr(:, :)
  real(WP), allocatable :: term1(:, :)
  real(WP), allocatable :: term2(:, :)
  real(WP), allocatable :: m(:, :)
  real(WP), allocatable :: dummy(:, :)
  real(WP), allocatable :: theta(:, :)
  real(WP), allocatable :: phi_dx(:, :), phi_dy(:, :)
  real(WP), allocatable :: phi_old(:, :)
  real(WP), allocatable :: epsilon_(:, :) !< 各向异性梯度能系数
  real(WP), allocatable :: epsilon_deriv(:, :)  !< 各向异性梯度能系数的导数, 对theta求导
  ! 初始化场变量数组
  call nucleus(Nx, Ny, seed, phi, tempr)
  allocate(theta(Nx, Ny), epsilon_(Nx, Ny), epsilon_deriv(Nx, Ny))
  allocate(lap_phi(Nx, Ny), lap_tempr(Nx, Ny), phi_old(Nx, Ny))
  allocate(term1(Nx, Ny), term2(Nx, Ny))
  allocate(m(Nx, Ny), dummy(Nx, Ny))
  call write_eta('./inital_phi.txt', phi)
  ! 开始迭代演化
  do istep = 1, nstep
    ! 计算方向导数
    call gradient_mat(Nx, Ny, dx, dy, matx=phi, matdx=phi_dx, matdy=phi_dy)
    ! 计算实际偏转角
    theta = atan2(phi_dy, phi_dx)
    ! 计算各向异性梯度能系数
    epsilon_ = epsilonb * ( 1 + delta * cos( aniso*( theta-theta0 ) ) )
    epsilon_deriv = -epsilonb * delta * aniso * sin( aniso*( theta-theta0 ) )
    ! 对phi和tempr进行拉普拉斯算子操作
    call laplace_2(Nx, Ny, idxdy, phi, lap_phi)
    call laplace_2(Nx, Ny, idxdy, tempr, lap_tempr)
    ! 存储好中间变量
    phi_old(:, :) = phi(:, :)
    dummy(:, :) = epsilon_(:, :) * epsilon_deriv(:, :)
    ! 计算控制方程的第1项, 即∂/(∂y)[ε*(∂ε/∂θ)*(∂φ/∂x)]
    call gradient_matdy(Nx, Ny, dx, dy, matx=dummy(:, :)*phi_dx(:, :), matdy=term1)
    ! 计算控制方程的第2项, 即-∂/(∂x)[ε*(∂ε/∂θ)*(∂φ/∂y)]
    call gradient_matdx(Nx, Ny, dx, dy, matx=dummy(:, :)*phi_dy(:, :), matdx=term2)
    ! 注意term2前有一个负号
    term2 = -term2
    ! 计算系数m
    m(:, :) = (alpha*ipi)*atan(gamma_*(Teq - tempr(:, :)))
    ! 对序参数phi进行时间积分
    phi(:, :) = phi(:, :) + (dt*itau)*( term1(:, :) + term2(:, :) + &
                epsilon_(:, :)*epsilon_(:, :)*lap_phi(:, :) + &     
                phi_old(:, :) * (1.0_WP - phi_old(:, :)) * (phi_old(:, :) - 0.5_WP + m(:, :)) )
    ! 接着对无量纲温度tempr进行时间积分
    tempr(:, :) = tempr(:, :) + dt*lap_tempr(:, :) + &
                  kappa*(phi(:, :) - phi_old(:, :))
    ! 打印信息
    if(mod(istep, nprint) == 0) then
      write(stdout, *) "Done step: ", istep
      call write_vtk_file('./VTK_phi/dump_', Nx, Ny, dx, dy, phi, istep)
      call write_vtk_file('./VTK_tempr/dump_', Nx, Ny, dx, dy, tempr, istep)
    end if
  end do 
  
end program main
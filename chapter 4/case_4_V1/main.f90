program main
  ! gfortran mod_case_4_V1.f90 main.f90 -O2; ./a.exe
  use case_4_V1, only: WP, stdout, gradient_mat, nucleus, laplacian_2, tau, epsilonb, mu, &
  kappa, delta, aniso, alpha, gamma_, Teq, theta0, ipi, itau, write_vtk_file, write_eta
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
  integer :: istep, i, j      !< 循环变量
  integer :: im, ip, jm, jp
  real(WP), allocatable :: lap_phi(:, :), lap_tempr(:, :)
  real(WP) :: phi_old
  real(WP) :: term1
  real(WP) :: term2
  real(WP) :: m
  real(WP), allocatable :: theta(:, :)
  real(WP), allocatable :: phi_dx(:, :), phi_dy(:, :)
  real(WP), allocatable :: epsilon_(:, :) !< 各向异性梯度能系数
  real(WP), allocatable :: epsilon_deriv(:, :)  !< 各向异性梯度能系数的导数, 对theta求导
  ! 初始化场变量数组
  call nucleus(Nx, Ny, seed, phi, tempr)
  allocate(theta(Nx, Ny), epsilon_(Nx, Ny), epsilon_deriv(Nx, Ny))
  allocate(lap_phi(Nx, Ny), lap_tempr(Nx, Ny))
  call write_eta('./inital_phi.txt', phi)
  ! 开始迭代演化
  do istep = 1, nstep
    ! 计算方向导数
    call gradient_mat(Nx, Ny, dx, dy, matx=phi, matdx=phi_dx, matdy=phi_dy)
    ! 计算实际偏转角
    theta = atan2(phi_dy, phi_dx)
    ! 计算各向异性梯度能系数
    epsilon_ = epsilonb * ( 1 + delta * cos( aniso*( theta-theta0 ) ) )
    epsilon_deriv = -epsilonb*delta*aniso*sin( aniso*( theta-theta0 ) )
    ! 逐一对phi和tempr进行拉普拉斯算子操作
    do j = 1, Ny
      do i = 1, Nx
        lap_phi(i, j) = laplacian_2(phi, i, j, Nx, Ny, idxdy)
        lap_tempr(i, j) = laplacian_2(tempr, i, j, Nx, Ny, idxdy)
      end do
    end do
    ! 计算控制方程的剩余项, 以用于时间积分
    do j = 1, Ny
      do i = 1, Nx
        ip = i + 1; im = i - 1
        jp = j + 1; jm = j - 1
        if(ip == Nx + 1) ip = 1
        if(im == 0) im = Nx
        if(jp == Ny + 1) jp = 1
        if(jm == 0) jm = Ny
        phi_old = phi(i, j)
        ! 计算控制方程的第一项, 即-∂/(∂x)[ε*(∂ε/∂θ)*(∂φ/∂y)]
        term1 = -(epsilon_(i, jp)* epsilon_deriv(i, jp)*phi_dy(i, jp) - &
                 epsilon_(i, jm)* epsilon_deriv(i, jm)*phi_dy(i, jm)) * idx
        ! 计算控制方程第二项, 即∂/(∂y)[ε*(∂ε/∂θ)*(∂φ/∂x)]
        term2 = (epsilon_(ip, j)* epsilon_deriv(ip, j)*phi_dx(ip, j) - &
                 epsilon_(im, j)* epsilon_deriv(im, j)*phi_dx(im, j)) * idy
        ! 计算系数m
        m = (alpha*ipi)*atan(gamma_*(Teq - tempr(i, j)))
        ! 对序参数phi进行时间积分
        phi(i, j) = phi(i, j) + (dt*itau)*( term1 + term2 + &
                    epsilon_(i, j)*epsilon_(i, j)*lap_phi(i, j) + &     
                    phi_old*(1.0_WP - phi_old)*(phi_old - 0.5_WP + m) )
        ! 接着对无量纲温度tempr进行时间积分
        tempr(i, j) = tempr(i, j) + dt*lap_tempr(i, j) + &
                      kappa*(phi(i, j) - phi_old)
      end do
    end do
    ! 打印信息
    if(mod(istep, nprint) == 0) then
      write(stdout, *) "Done step: ", istep
      call write_vtk_file('./VTK_phi/dump_', Nx, Ny, dx, dy, phi, istep)
      call write_vtk_file('./VTK_tempr/dump_', Nx, Ny, dx, dy, tempr, istep)
    end if
  end do 
  
end program main
program main
  ! gfortran .\mod_case2_v1.f90 .\main.f90 -O2; ./a.exe
  use case_2_v1, only: WP, stdout, init_grain_micro, write_vtk_file,&
                       free_energ_fd_ca_v1, mobil, grcoeff, write_eta
  implicit none
  integer, parameter :: Nx = 64
  integer, parameter :: Ny = 64
  integer, parameter :: NxNy = Nx*Ny
  integer, parameter :: nstep = 10000
  integer, parameter :: nprint = 100
  real(WP), parameter :: dx = 0.5_WP !< x方向上相邻节点的距离为dx
  real(WP), parameter :: dy = 0.5_WP !< y方向上相邻节点的距离为dx
  real(WP), parameter :: dt = 0.005_WP  !< 时间步长
  real(WP), parameter :: iNxNy = 1.0_WP / (Nx*Ny) !< Nx*Ny的导数
  real(WP), allocatable, target :: etas(:, :, :)
  real(WP), allocatable :: eta2(:, :)
  real(WP), pointer :: eta(:, :)  !< 计算第i个晶粒的数组eta 
  real(WP) :: lap_eta
  real(WP) :: dfdeta
  real(WP) :: grain_sum
  integer, allocatable :: glist(:)  !< 记录晶粒是否存活的数组, 1表示存活, 0表示湮灭
  integer :: ngrain       !< 总共有几个晶粒
  integer :: i, j, k, igrain
  integer :: ip, im, jp, jm

  call init_grain_micro(etas, ngrain, glist, Nx, Ny, dx, dy, iflag=2)
  allocate(eta2(Nx, Ny), source=0.0_WP)
  do k = 1, nstep
    ! 遍历所有晶粒的序参数
    do igrain = 1, ngrain
      ! 如果第i个晶粒还存在的话
      if(glist(igrain) == 1) then
        ! 令eta指向第igrain个晶粒的序参数数组
        eta => etas(:, :, igrain)
        do j = 1, Ny
          do i = 1, Nx
            ! 施加周期性边界条件
            ip = i + 1; im = i - 1
            jp = j + 1; jm = j - 1
            if(ip == Nx + 1) ip = 1
            if(im == 0) im = Nx
            if(jp == Ny + 1) jp = 1
            if(jm == 0) jm = Ny
            ! 计算拉普拉斯算子的值
            lap_eta = (eta(ip, j) + eta(im, j) + eta(i, jp) + eta(i, jm) - &
                        4.0_WP * eta(i, j)) / (dx*dy)
            ! 计算导数
            dfdeta = free_energ_fd_ca_v1(i, j, ngrain, etas, igrain)
            ! 进行时域积分
            eta(i, j) = eta(i, j) - dt*mobil*(dfdeta - grcoeff*lap_eta)
            ! 修正数值误差
            if(eta(i, j) >= 0.9999_WP) eta(i, j) = 0.9999_WP
            if(eta(i, j) <= 0.0001_WP) eta(i, j) = 0.0001_WP
          end do
        end do
        ! 计算第igrain个晶粒的体积分数
        grain_sum = sum(eta) * iNxNy
        if(grain_sum < 0.001_WP) glist(igrain) = 0
      end if
    end do
    if(mod(k, nprint) == 0) then
      eta2 = 0.0_WP
      do i = 1, ngrain
        eta2 = eta2 + etas(:, :, i)*etas(:, :, i) 
      end do
      call write_vtk_file('./VTK_polycrystal/eta_', Nx, Ny, dx, dy, eta2, k)
      write(stdout, *) "Done step: ", k
    end if
  end do

end program main
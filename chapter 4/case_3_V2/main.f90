program main
  ! gfortran mod_case_3_V2.f90 main.f90 -O2; ./a.exe
  use case_3_V2, only: WP, stdout, micro_sint_pre, free_energ_sint_v1, &
  write_eta, write_vtk_file, coefk, coefl, coefm, Dgrb, Dsur, Dvap, Dvol, &
  laplacian_2, laplacian_4, free_energ_sint_v2
  implicit none
  integer, parameter :: npart = 9
  integer, parameter :: Nx = 100
  integer, parameter :: Ny = 100
  real(WP), parameter :: dx = 0.5_WP
  real(WP), parameter :: dy = 0.5_WP
  real(WP), parameter :: idxdy = 1.0_WP / (dx*dy)
  real(WP), parameter :: dt = 1.0e-4_WP
  integer :: nstep = 5000
  integer :: nprint = 100
  integer :: istep
  integer :: i, j
  integer :: ipart, jpart
  integer :: iflag
  real(WP), allocatable, target :: dfunc(:, :, :)
  real(WP), pointer :: dfdcon(:, :), dfdeta(:, :)
  real(WP), pointer :: eta(:, :)
  real(WP), allocatable, target :: etas(:, :, :)
  real(WP), allocatable :: eta2(:, :)
  real(WP), allocatable :: con(:, :)
  real(WP), allocatable :: dummy(:, :)
  real(WP) :: phi
  real(WP) :: lap_dummy, lap_eta
  real(WP) :: sum_
  real(WP) :: mobil

  
  allocate(eta(Nx, Ny), source=0.0_WP)
  allocate(dfunc(Nx, Ny, 2))
  dfdcon => dfunc(:, :, 1); dfdeta => dfunc(:, :, 2)
  allocate(eta2(Nx, Ny))
  allocate(dummy(Nx, Ny))
  ! 初始化序参数数组和密度数组
  call micro_sint_pre(Nx, Ny, npart, 1, etas, con)
  ! 将结果导出以进行可视化
  call write_eta('./data/con_initial.txt', con)
  call write_eta("./data/eta_grain1.txt", etas(:, :, 1))
  call write_eta("./data/eta_grain2.txt", etas(:, :, 2))
  ! 开始进行迭代
  do istep = 1, nstep
    ! 首先演化浓度数组con, 对con进行第1次拉普拉斯算子操作
    iflag = 1
    ! 计算自由能对浓度concetration的导数
    call free_energ_sint_v2(Nx, Ny, con, eta, etas, dfunc, iflag)
    ! 对con进行第2次拉普拉斯算子操作
    do j = 1, Ny
      do i = 1, Nx
        lap_dummy = laplacian_2(dfdcon, i, j, Nx, Ny, idxdy) - &
             0.5_WP*coefm*laplacian_4(con, i, j, Nx, Ny, idxdy)
        ! 计算流动性系数phi
        phi = con(i, j)**3 *(10.0_WP - 15.0_WP*con(i, j) + 6*con(i, j)**2)
        sum_ = 0.0_WP
        do ipart = 1, npart
          do jpart = 1, npart
            if (ipart /= jpart) then
              sum_ = sum_ + etas(i, j, ipart)*etas(i, j, jpart)
            end if
          end do
        end do
        ! 计算与微结构相关的扩散系数D, 对应公式(4.40)
        mobil = Dvol*phi + Dvap*(1.0_WP - phi) + &
                Dsur*con(i, j)*(1.0_WP - con(i, j)) + &
                2*Dgrb*sum_ ! 疑似需要除以2
        ! 进行时间积分
        con(i, j) = con(i, j) + dt*mobil*lap_dummy
      end do
    end do
    ! 数值误差校正
    where (con > 0.9999999_WP) 
      con = 0.9999999_WP
    else where (con < 0.000001_WP) 
      con = 0.0000001_WP
    end where
    ! 现在开始对晶粒的序参数进行演化 
    iflag = 2
    do ipart = 1, npart
      eta => etas(:, :, ipart)
      call free_energ_sint_v2(Nx, Ny, con, eta, etas, dfunc, iflag)
      do j = 1, Ny
        do i = 1, Nx
          lap_eta = laplacian_2(eta, i, j, Nx, Ny, idxdy)
          ! 进行时间积分
          eta(i, j) = eta(i, j) - dt*coefl*(dfdeta(i, j) - 0.5_WP*coefl*lap_eta)
        end do
      end do
      ! 数值误差校正
      where (eta > 0.9999999_WP) 
        eta = 0.9999999_WP
      else where (eta < 0.000001_WP) 
        eta = 0.000001_WP
      end where
    end do
    ! 输出信息
    if ( mod(istep, nprint) == 0 ) then
      write(stdout, '(A, I0)') "Done step: ", istep
      eta2 = sum(etas(:, :, :)*etas(:, :, :), dim=3)
      call write_vtk_file('./VTK_con/con_', Nx, Ny, dx, dy, con, istep)
      call write_vtk_file('./VTK_eta2/eta2_', Nx, Ny, dx, dy, eta2, istep)
    end if
  end do
  write(stdout, '(A)') "-----------------------------All done!-----------------------------"
end program main
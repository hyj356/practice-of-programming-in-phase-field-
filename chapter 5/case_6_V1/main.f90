program hello
    use case_6_V1, only: wp, cwp, stdout, fft_2d, ifft_2d, micro_ch_pre,&
                         prepare_fft, write_vtk_file, free_energ_ch_v1, coefA
    implicit none
    integer, parameter :: Nx = 64
    integer, parameter :: Ny = 64
    real(wp), parameter :: iNxNy = 1.0_wp / (Nx*Ny)
    integer, parameter :: nstep = 20000
    integer, parameter :: nprint = 200
    real(wp), parameter :: dx = 1.0_wp
    real(wp), parameter :: dy = 1.0_wp
    real(wp), parameter :: dt = 1.0e-2_wp
    real(wp), parameter :: c0 = 0.4_wp
    real(wp), parameter :: mobility = 1.0_wp
    real(wp), parameter :: grad_coef = 0.5_wp
    real(wp), allocatable :: con(:, :)
    complex(cwp), allocatable :: conk(:, :)
    real(wp), allocatable :: dfdcon(:, :)
    complex(cwp), allocatable :: dfdconk(:, :)
    real(wp), allocatable :: k2(:, :), k4(:, :), kx(:), ky(:)
    complex(cwp) :: numer
    real(wp) :: denom
    integer :: istep, i, j

    ! 分配内存
    allocate(conk(Nx, Ny), dfdcon(Nx, Ny), dfdconk(Nx, Ny))
    ! 初始化浓度数组
    call micro_ch_pre(con, Nx, Ny, c0)
    call write_vtk_file("./VTK/dump_", Nx, Ny, dx, dy, con, 0)
    ! 初始化傅里叶变换的系数
    call prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)
    ! 对浓度数组进行傅里叶变换
    call fft_2d(input=con, output=conk, Nx=Nx, Ny=Ny)
    ! 开始进行时间积分
    do istep = 1, nstep
        ! 计算自由能导数
        do j = 1, Ny
            do i = 1, Nx
                dfdcon(i, j) = free_energ_ch_v1(con, i, j)
            end do
        end do
        ! 进行傅里叶变换
        call fft_2d(input=dfdcon, output=dfdconk, Nx=Nx, Ny=Ny)
        ! 进行时间积分
        do j = 1, Ny
            do i = 1, Nx
                numer = dt*mobility*k2(i, j)*dfdconk(i, j)
                denom = 1.0_wp + dt*coefA*mobility*grad_coef*k4(i, j)
                conk(i, j) = (conk(i, j) - numer) / denom
            end do
        end do
        ! 调用逆FFT将虚数变成实数
        call ifft_2d(input=conk, output=con, Nx=Nx, Ny=Ny, iNxNy=iNxNy)
        ! 数值误差校正
        where (con > 0.99999_WP)
            con = 0.99999_WP
        else where (con < 0.00001_WP)
            con = 0.00001_WP
        end where
        ! 输出vtk文件
        if (mod(istep, nprint) == 0) then
            write(stdout, *) "Done step: ", istep
            ! 写出vtk文件
            call write_vtk_file("./VTK/dump_", Nx, Ny, dx, dy, con, istep)
        end if
    end do

end program


program hello
    use case_10_V1, only: write_vtk_file, prepare_fft, fft_2d, ifft_2d, &
                           green_tensor_2d, solve_elasticity_v1, &
                           micro_ch_pre, FeCr_chem_potent_v1, dislo_strain, &
                           wp, cwp, stdout, free_energ_ch
    implicit none
    ! 模拟区域大小定义参数
    integer, parameter :: Nx = 128
    integer, parameter :: Ny = 128
    real(wp), parameter :: iNxNy = 1.0_wp / (Nx*Ny)
    real(wp), parameter :: dx = 1.0_wp
    real(wp), parameter :: dy = 1.0_wp
    ! 时间积分参数
    integer, parameter :: nstep = 10000
    integer, parameter :: nprint = 50
    real(wp), parameter :: dt = 1.0e-2_wp
    real(wp), parameter :: coefA = 2.0_wp
    ! 材料相关常数
    real(wp), parameter :: c0 = 0.2_wp
    real(wp), parameter :: mobility = 1.0_wp
    real(wp), parameter :: grad_coef = 0.5_wp
    real(wp), parameter :: tempr = 535.0_wp
    real(wp), parameter :: RT = 8.314462_wp*tempr
    ! 材料的弹性常数
    real(wp), parameter :: cm11 = 233.10e3_wp
    real(wp), parameter :: cm12 = 135.44e3_wp
    real(wp), parameter :: cm44 = 178.30e3_wp
    !--------------------------------------
    real(wp), parameter :: cp11 = 350.00e3_wp
    real(wp), parameter :: cp12 =  67.80e3_wp
    real(wp), parameter :: cp44 = 100.80e3_wp
    ! 特征应变
    real(wp), parameter :: ei0 = 0.006_wp
    ! 施加的应变
    real(wp), dimension(3), parameter :: ea = [0.0_wp, 0.0_wp, 0.0_wp]
    ! 仿真用到的数组
    real(wp), allocatable :: s11(:, :), s12(:, :), s22(:, :)
    real(wp), allocatable :: e11(:, :), e12(:, :), e22(:, :)
    real(wp), allocatable :: ed11(:, :), ed12(:, :), ed22(:, :)
    real(wp), allocatable :: cr(:, :), dfdcr(:, :), tmatx(:, :, :, :, :, :)
    real(wp), allocatable :: delsdc(:, :)
    real(wp), allocatable :: kx(:), ky(:), k2(:, :), k4(:, :)
    complex(cwp), allocatable :: crk(:, :), dfdcrk(:, :), delsdck(:, :)
    ! 仿真用到的局部变量和循环变量
    integer :: istep, i, j
    real(wp) :: denom
    complex(cwp) :: numer

    ! 分配内存
    allocate(s11(Nx, Ny), s12(Nx, Ny), s22(Nx, Ny), source=0.0_wp)
    allocate(e11(Nx, Ny), e12(Nx, Ny), e22(Nx, Ny), source=0.0_wp)
    allocate(ed11(Nx, Ny), ed12(Nx, Ny), ed22(Nx, Ny), source=0.0_wp)
    allocate(dfdcr(Nx, Ny), delsdc(Nx, Ny))
    allocate(crk(Nx, Ny), dfdcrk(Nx, Ny), delsdck(Nx, Ny))

    ! 施加位错周围的位移场
    call dislo_strain(ed11, ed22, ed12, Nx, Ny, idislo=2)
    !write(stdout, *) ed12(34:39, 64)
    ! 初始化各个节点的物质浓度
    call micro_ch_pre(Nx, Ny, c0, cr)
    call write_vtk_file("./VTK/dump_", Nx, Ny, dx, dy, cr, 0)

    ! 计算傅里叶变换的系数数组
    call prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)

    ! 计算green tensor
    call green_tensor_2d(tmatx, Nx, Ny, kx, ky, cm11, cm12, cm44, cp11, cp12, cp44)

    ! 开始演化
    do istep = 1, nstep
        ! 计算弹性势能的导数
        call solve_elasticity_v1(Nx, Ny, cm11, cm12, cm44, cp11, cp12, cp44, ed11, ed12, ed22, &
                                 ei0, ea, cr, s11, s12, s22, e11, e12, e22, delsdc, tmatx)
        ! 乘以系数变换
        delsdc(:, :) = delsdc(:, :) / RT
        ! 计算化学势能导数
        do j = 1, Ny
            do i = 1, Nx
                dfdcr(i, j) = FeCr_chem_potent_v1(cr, i, j, tempr)
            end do
        end do
        ! 进行傅里叶变换
        call fft_2d(input=cr, output=crk, Nx=Nx, Ny=Ny)
        call fft_2d(input=dfdcr, output=dfdcrk, Nx=Nx, Ny=Ny)
        call fft_2d(input=delsdc, output=delsdck, Nx=Nx, Ny=Ny)
        ! 进行时间积分
        do j = 1, Ny
            do i = 1, Nx
                numer = dt*mobility*k2(i, j)*(dfdcrk(i, j) + delsdck(i, j))
                denom = 1.0_wp + dt*coefA*mobility*grad_coef*k4(i, j)
                crk(i, j) = (crk(i, j) - numer) / denom
            end do
        end do
        ! 通过逆傅里叶变换将crk变为cr
        call ifft_2d(input=crk, output=cr, Nx=Nx, Ny=Ny, iNxNy=iNxNy)
        ! 数值误差校正
        where(cr > 0.9999_wp)
            cr = 0.9999_wp
        else where(cr < 0.00001_wp)
            cr = 0.00001_wp
        end where
        ! 输出求解信息
        if(mod(istep, nprint) == 0) then
            write(stdout, *) "Done step: ", istep
            call write_vtk_file("./VTK/dump_", Nx, Ny, dx, dy, cr, istep)
        end if
    end do

    print *, "-------------------------------------------All done!-------------------------------------------"

end program


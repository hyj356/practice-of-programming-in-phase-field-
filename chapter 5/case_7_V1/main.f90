
program hello
    use case_7_V1, only: wp, cwp, stdout, pi, coefA, &
        write_vtk_file, prepare_fft, fft_2d, ifft_2d, &
        init_grain_micro, free_energ_fd_ca_v1
    implicit none
    integer, parameter :: Nx = 64
    integer, parameter :: Ny = 64
    real(wp), parameter :: iNxNy = 1.0_wp / (Nx*Ny)
    real(wp), parameter :: dx = 0.5_wp
    real(wp), parameter :: dy = 0.5_wp
    real(wp), parameter :: dt = 0.005
    integer, parameter :: nstep = 5000
    integer, parameter :: nprint = 50
    integer, parameter :: iflag = 2     !< iflag为1时, 模型含有2个晶粒, iflag=2, 模型含有25个晶粒
    ! 材料相关的参数
    real(wp), parameter :: mobil = 5.0_wp
    real(wp), parameter :: grcoef = 0.1_wp
    ! 仿真的数据数组
    real(wp), allocatable, target :: etas(:, :, :)
    complex(cwp), allocatable :: etak(:, :)
    real(wp), allocatable :: eta2(:, :)
    real(wp), pointer :: eta(:, :)
    real(wp), allocatable :: dfdeta(:, :)
    complex(cwp), allocatable :: dfdetak(:, :)
    real(wp), allocatable :: kx(:), ky(:)
    real(wp), allocatable :: k2(:, :), k4(:, :)
    integer, allocatable :: glist(:)
    integer :: ngrain
    integer :: i, j, istep, igrain
    real(wp) :: grain_sum
    complex(cwp) :: numer       !< 欧拉半隐式时间积分公式中的分子
    real(wp) :: denom    !< 欧拉半隐式时间积分公式中的分母

    ! 分配内存
    allocate(dfdeta(Nx, Ny), dfdetak(Nx, Ny), eta2(Nx, Ny), etak(Nx, Ny))
    ! 初始化各个晶粒的序参数
    call init_grain_micro(etas, ngrain, glist, Nx, Ny, dx, dy, iflag)
    ! 准备二维傅里叶变换的对应系数
    call prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)
    ! 开始迭代
    do istep = 1, nstep
        ! 遍历所有晶粒
        do igrain = 1, ngrain
            ! 如果第igrain个晶粒还没有湮灭
            if(glist(igrain) == 1) then
                eta => etas(:, :, igrain)
                ! 计算每个节点的自由能导数
                do j = 1, Ny
                    do i = 1, Nx
                        dfdeta(i, j) = free_energ_fd_ca_v1(i, j, ngrain, etas, igrain)
                    end do
                end do
                ! 对第igrain个晶粒的序参数eta进行傅里叶变换
                call fft_2d(input=eta, output=etak, Nx=Nx, Ny=Ny)
                ! 对自由能导数进行傅里叶变换
                call fft_2d(input=dfdeta, output=dfdetak, Nx=Nx, Ny=Ny)
                ! 逐一对每个节点计算系数, 进行时间积分
                do j = 1, Ny
                    do i = 1, Nx
                        numer = dt*mobil*dfdetak(i, j)
                        denom = 1.0_wp + dt*coefA*mobil*grcoef*k2(i, j)
                        etak(i, j) = (etak(i, j) - numer) / denom
                    end do
                end do
                ! 将etak转到实数域中
                call ifft_2d(input=etak, output=eta, Nx=Nx, Ny=Ny, iNxNy=iNxNy)
                ! 数值误差校正
                where (eta > 0.99999_WP)
                    eta = 0.99999_WP
                else where (eta < 0.00001_WP)
                    eta = 0.00001_WP
                end where
                ! 计算第igrain个晶粒还有多少晶粒体积分数
                grain_sum = sum(eta) * iNxNy
                if (grain_sum < 0.001_wp) glist(igrain) = 0
            end if
        end do
        if(mod(istep, nprint) == 0) then
            write(stdout, *) "Done step: ", istep
            eta2(:, :) = sum(etas(:, :, :)*etas(:, :, :), dim=3)
            ! 写出vtk文件
            call write_vtk_file("./VTK_polycrystal/dump_", Nx, Ny, dx, dy, eta2, istep)
        end if
    end do

end program


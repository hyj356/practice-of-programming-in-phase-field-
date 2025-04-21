program hello
    use case_8_V1,only: &!以下为模组中的函数和子程序
        write_vtk_file, prepare_fft, fft_2d, ifft_2d, &
        init_FeCuNiMn, Fe_Cu_Mn_Ni_free_energy, &
        !以下为模组中定义的变量和常量
        wp, cwp, stdout, pi
    implicit none
    ! 模拟区域相关常量
    integer, parameter :: Nx = 128
    integer, parameter :: Ny = 128
    real(wp), parameter :: dx = 0.5_wp
    real(wp), parameter :: dy = 0.5_wp
    real(wp), parameter :: iNxNy = 1.0_wp / (Nx*Ny)
    ! 时间积分相关常量

    integer, parameter :: nstep = 5000
    integer, parameter :: nprint = 50
    real(wp), parameter :: dt = 1.0e-2_wp
    real(wp), parameter :: coefA = 1.0_WP
    ! 模拟的合金材料FeCuMnNi的相关常数
    real(wp), parameter :: Cu0 = 0.15_wp    !< Cu元素浓度
    real(wp), parameter :: Mn0 = 0.01_wp    !< Mn元素浓度
    real(wp), parameter :: Ni0 = 0.01_wp    !< Ni元素浓度
    real(wp), parameter :: grcoef_Cu = 0.68884_wp
    real(wp), parameter :: grcoef_Mn = 0.68884_wp
    real(wp), parameter :: grcoef_Ni = 0.68884_wp
    real(wp), parameter :: grcoef_or = 0.13729_wp
    ! 计算迁移率会用到的扩散系数
    real(wp), parameter :: gconst = 8.314472_wp
    real(wp), parameter :: tempr = 823.0_wp
    real(wp), parameter :: RT = gconst * tempr
    ! Cu元素
    real(wp), parameter :: QACu = 2.44e5_wp
    real(wp), parameter :: QGCu = 2.80e5_wp
    real(wp), parameter :: D0ACu = 4.7e-5_wp
    real(wp), parameter :: D0GCu = 4.3e-5_wp
    real(wp), parameter :: DCuA = (D0ACu*exp(-QACu/RT))
    real(wp), parameter :: DCuG = (D0GCu*exp(-QGCu/RT))/DCuA
    ! Ni元素
    real(wp), parameter :: QANi = 2.56e5_wp
    real(wp), parameter :: QGNi = 2.73e5_wp
    real(wp), parameter :: D0ANi = 1.4e-4_wp
    real(wp), parameter :: D0GNi = 1.08e-5_wp
    real(wp), parameter :: DNiA = (D0ANi*exp(-QANi/RT))/DCuA
    real(wp), parameter :: DNiG = (D0GNi*exp(-QGNi/RT))/DCuA
    ! Mn元素
    real(wp), parameter :: QAMn = 2.63e5_wp
    real(wp), parameter :: QGMn = 2.64e5_wp
    real(wp), parameter :: D0AMn = 1.49e-4_wp
    real(wp), parameter :: D0GMn = 2.78e-5_wp
    real(wp), parameter :: DMnA = (D0AMn*exp(-QAMn/RT))/DCuA
    real(wp), parameter :: DMnG = (D0GMn*exp(-QGMn/RT))/DCuA
    ! 迁移率常数
    real(wp) :: mcoef_Cu, mcoef_Mn, mcoef_Ni
    real(wp), parameter :: mcoef_orp = 0.1_wp
    ! 仿真的数据数组
    real(wp), allocatable :: Cu(:, :), Mn(:, :), Ni(:, :), orp(:, :)
    complex(cwp), allocatable :: Cuk(:, :), Mnk(:, :), Nik(:, :), orpk(:, :)
    real(wp), allocatable :: kx(:), ky(:), k2(:, :), k4(:, :)
    real(wp), allocatable :: dGdCu(:, :), dGdMn(:, :), dGdNi(:, :), dGdor(:, :)
    complex(cwp), allocatable :: dGdCuk(:, :), dGdMnk(:, :), dGdNik(:, :), dGdork(:, :)
    ! 循环变量
    integer :: istep, i, j

    ! 分配内存
    allocate(dGdCu(Nx, Ny), dGdMn(Nx, Ny), dGdNi(Nx, Ny), dGdor(Nx, Ny))
    allocate(Cuk(Nx, Ny), Mnk(Nx, Ny), Nik(Nx, Ny), orpk(Nx, Ny))
    allocate(dGdCuk(Nx, Ny), dGdMnk(Nx, Ny), dGdNik(Nx, Ny), dGdork(Nx, Ny))
    ! 初始化各个元素的浓度
    call init_FeCuNiMn(Cu, Mn, Ni, orp, Cu0, Mn0, Ni0, Nx, Ny)
    call write_vtk_file("./VTK_orp/dump_", Nx, Ny, dx, dy, orp, 0)
    ! 计算fft变换对应的系数
    call prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)
    ! 开始演化
    do istep = 1, nstep
        ! 进行正向傅里叶变换
        call fft_2d(input=Cu, output=Cuk, Nx=Nx, Ny=Ny)     !< Cu元素
        call fft_2d(input=Mn, output=Mnk, Nx=Nx, Ny=Ny)     !< Mn元素
        call fft_2d(input=Ni, output=Nik, Nx=Nx, Ny=Ny)     !< Ni元素
        call fft_2d(input=orp, output=orpk, Nx=Nx, Ny=Ny)   !< 序参数
        ! 计算自由能导数
        call Fe_Cu_Mn_Ni_free_energy(Nx, Ny, tempr, Cu, Mn, Ni, orp, dGdCu, dGdMn, dGdNi, dGdor)
        ! 对自由能导数进行正向傅里叶变换
        call fft_2d(input=dGdCu, output=dGdCuk, Nx=Nx, Ny=Ny)
        call fft_2d(input=dGdMn, output=dGdMnk, Nx=Nx, Ny=Ny)
        call fft_2d(input=dGdNi, output=dGdNik, Nx=Nx, Ny=Ny)
        call fft_2d(input=dGdor, output=dGdork, Nx=Nx, Ny=Ny)
        ! 对每个节点进行时间积分
        do j = 1, Ny
            do i = 1, Nx
                ! 计算不同节点的迁移率常数
                mcoef_Cu = Cu0*(1.0_wp-Cu0)*(1.0_wp*(1.0_wp-orp(i,j))+DCuG*orp(i,j))  !< Cu元素
                mcoef_Ni = Ni0*(1.0_wp-Ni0)*(DNiA*(1.0_wp-orp(i,j))+DNiG*orp(i,j))    !< Ni元素
                mcoef_Mn = Mn0*(1.0_wp-Mn0)*(DMnA*(1.0_wp-orp(i,j))+DMnG*orp(i,j))    !< Mn元素
                ! 使用半隐式欧拉法进行时间积分
                Cuk(i,j)=(Cuk(i,j)-dt*mcoef_Cu*dgdCuk(i,j)*k2(i,j))/(1.0_wp+dt*grcoef_Cu*k4(i,j)*mcoef_Cu)
                Nik(i,j)=(Nik(i,j)-dt*mcoef_Ni*dgdNik(i,j)*k2(i,j))/(1.0_wp+dt*grcoef_Ni*k4(i,j)*mcoef_Ni)
                Mnk(i,j)=(Mnk(i,j)-dt*mcoef_Mn*dgdMnk(i,j)*k2(i,j))/(1.0_wp+dt*grcoef_Mn*k4(i,j)*mcoef_Mn)
                orpk(i,j) = (orpk(i,j) - dt*dgdork(i,j)*mcoef_orp) / &
                            (1.0_wp + dt*grcoef_or*k2(i,j)*mcoef_orp)
            end do
        end do
        ! 进行逆傅里叶变换
        call ifft_2d(input=Cuk, output=Cu, Nx=Nx, Ny=Ny, iNxNy=iNxNy)     !< Cu元素
        call ifft_2d(input=Mnk, output=Mn, Nx=Nx, Ny=Ny, iNxNy=iNxNy)     !< Mn元素
        call ifft_2d(input=Nik, output=Ni, Nx=Nx, Ny=Ny, iNxNy=iNxNy)     !< Ni元素
        call ifft_2d(input=orpk, output=orp, Nx=Nx, Ny=Ny, iNxNy=iNxNy)   !< 序参数
        ! 输出相关信息
        if(mod(istep, nprint) == 0) then
            call write_vtk_file("./VTK_orp/dump_", Nx, Ny, dx, dy, orp, istep)
            write(stdout, *) "Done step: ", istep
        end if
    end do

    write(stdout, *) "---------------------------------All done!---------------------------------"

end program


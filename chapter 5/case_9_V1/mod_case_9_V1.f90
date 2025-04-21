module case_9_V1
    use, intrinsic:: iso_c_binding, only: c_double, c_int, c_double_complex
    use, intrinsic:: iso_fortran_env, only:real32, real64, output_unit
    use FFTW, only: fftw_plan_dft_2d, fftw_execute_dft, fftw_destroy_plan, &
                  FFTW_FORWARD, FFTW_BACKWARD, &
                  c_ptr, FFTW_ESTIMATE
    implicit none
    private
    public :: write_vtk_file, prepare_fft, fft_2d, ifft_2d, &
                green_tensor_2d, solve_elasticity_v1, &
                micro_ch_pre, free_energ_ch

    integer, parameter, public :: wp = c_double
    integer, parameter, public :: cwp = c_double_complex
    integer, parameter, public :: ci = c_int
    integer, parameter, public :: stdout = output_unit
    real(WP), parameter, public :: pi = 4.0_wp * atan(1.0_WP)
    integer, parameter :: niter = 10
    real(wp), parameter :: tolerence = 0.001_wp

contains

    pure function free_energ_ch(con, i, j) result(ret)
        !! 计算化学能函数f(c)的导数, 在书上的案例中, f(c) = A*c^2*(1-c)^2
        !! 那么, f'(c) = 2*A*[c*(1-c)^2 - c^2*(1-c)] = 2*A*c*(1-c)
        real(WP), intent(in) :: con(:, :) !< 浓度数组, 形状为con(Nx, Ny)
        integer, intent(in) :: i, j       !< 指明计算几行几列的节点
        real(WP) :: ret, A

        A = 1.0_wp
        ret = 2.0_WP * A * con(i, j) * (1.0_WP - con(i, j)) * (1.0_WP - 2*con(i, j))

    end function free_energ_ch

    subroutine micro_ch_pre(Nx, Ny, c0, con)
        !! 此函数根据用户的输入来初始化初始的浓度数组矩阵
        integer, intent(in) :: Nx   !< x方向上的节点个数
        integer, intent(in) :: Ny   !< y方向上的节点个数
        real(WP), intent(in) :: c0  !< 平均浓度
        real(WP), intent(inout), allocatable :: con(:, :) !< 需要初始化的浓度数组, 大小为con(Nx, Ny)
        real(WP) :: noise   !< 根据原文要求, 每个节点距离c0的噪音为0.02
        real(WP) :: rng
        integer :: i, j

        noise = 0.02_WP
        if(.not. allocated(con)) allocate(con(Nx, Ny))
        call random_seed()
        ! 开始遍历初始化
        do j = 1, Nx
          do i = 1, Ny
            call random_number(rng)
            con(j, i) = c0 + noise * (0.5_WP - rng)
          end do
        end do
    end subroutine micro_ch_pre

    subroutine solve_elasticity_v1(Nx, Ny, cm11, cm12, cm44, cp11, cp12, cp44, &
                                    ed11, ed12, ed22, ei0, ea, con, s11, s12,  &
                                    s22, e11, e12, e22, delsdc, tmatx)
        ! 函数的输入变量
        integer, intent(in) :: Nx, Ny
        real(wp),intent(in) :: tmatx(:, :, :, :, :, :)
        real(wp), intent(in) :: cm11, cm12, cm44
        real(wp), intent(in) :: cp11, cp12, cp44
        real(wp), intent(in) :: ei0, ea(3)
        real(wp), intent(in) :: con(:, :)       !< 形状为(Nx, Ny)
        real(wp), allocatable, intent(inout) :: delsdc(:, :)    !< 形状为(Nx, Ny)
        real(wp), intent(inout) :: s11(:, :), s12(:, :), s22(:, :)    !< 形状为(Nx, Ny)
        real(wp), intent(inout) :: e11(:, :), e12(:, :), e22(:, :)    !< 形状为(Nx, Ny)
        real(wp), intent(in) :: ed11(:, :), ed12(:, :), ed22(:, :)     !< 形状为(Nx, Ny)
        ! 函数的局部变量
        real(wp), allocatable :: ei11(:, :), ei22(:, :), ei33(:, :), ei12(:, :)
        real(wp), allocatable :: c11(:, :), c12(:, :), c44(:, :)
        real(wp), allocatable :: sum_stres(:, :)
        complex(cwp), allocatable :: s11k(:, :), s12k(:, :), s22k(:, :)
        complex(cwp), allocatable :: e11k(:, :), e12k(:, :), e22k(:, :)
        complex(cwp), allocatable :: smatx(:, :, :, :), ematx(:, :, :, :)
        integer :: i, j, iter
        integer :: kk, ll, ii, jj
        real(wp) :: normF, old_norm, conver
        real(wp) :: et11, et22, et12

        ! 分配内存
        if(.not. allocated(delsdc)) allocate(delsdc(Nx, Ny))
        allocate(ei11(Nx, Ny), ei22(Nx, Ny), ei33(Nx, Ny))
        allocate(ei12(Nx, Ny), source = 0.0_wp)
        allocate(c11(Nx, Ny), c12(Nx, Ny), c44(Nx, Ny), sum_stres(Nx, Ny))
        allocate(s11k(Nx, Ny), s12k(Nx, Ny), s22k(Nx, Ny))
        allocate(e11k(Nx, Ny), e12k(Nx, Ny), e22k(Nx, Ny))
        allocate(smatx(Nx, Ny, 2, 2), ematx(Nx, Ny, 2, 2))
        ! 进行初始化
        do j = 1, Ny
            do i = 1, Nx
                ! 计算特征应变, 其中ei12已经在一开始的分配内存时被初始化为0
                ei11(i,j) = ei0*con(i,j)
                ei22(i,j) = ei0*con(i,j)
                ei33(i,j) = ei0*con(i,j)
                ! 使用Vergards law计算有效弹性常数
                c11(i,j) = con(i,j)*cp11 +(1.0_wp - con(i,j))*cm11
                c12(i,j) = con(i,j)*cp12 +(1.0_wp - con(i,j))*cm12
                c44(i,j) = con(i,j)*cp44 +(1.0_wp - con(i,j))*cm44
            end do
        end do
        ! 开始迭代
        do iter = 1, niter
            ! 将应变转换到傅里叶空间中
            call fft_2d(input=e11, output=e11k, Nx=Nx, Ny=Ny)
            call fft_2d(input=e22, output=e22k, Nx=Nx, Ny=Ny)
            call fft_2d(input=e12, output=e12k, Nx=Nx, Ny=Ny)
            ! 将应力转换到傅里叶空间中
            call fft_2d(input=s11, output=s11k, Nx=Nx, Ny=Ny)
            call fft_2d(input=s22, output=s22k, Nx=Nx, Ny=Ny)
            call fft_2d(input=s12, output=s12k, Nx=Nx, Ny=Ny)
            ! 组装应力与应变张量
            do j = 1, Ny
                do i = 1, Nx
                    ! 组装应力张量
                    smatx(i,j,1,1) = s11k(i,j); smatx(i,j,1,2) = s12k(i,j)
                    smatx(i,j,2,1) = s12k(i,j); smatx(i,j,2,2) = s22k(i,j)
                    ! 组装应变张量
                    ematx(i,j,1,1) = e11k(i,j); ematx(i,j,1,2) = e12k(i,j)
                    ematx(i,j,2,1) = e12k(i,j); ematx(i,j,2,2) = e22k(i,j)
                end do
            end do
            ! 计算green tensor
            do j = 1, Ny
                do i = 1, Nx
                    ! 格林操作符
                    do kk = 1, 2
                        do ll = 1, 2
                            do ii = 1, 2
                                do jj = 1, 2
                                    ematx(i,j,ii,jj) = ematx(i,j,ii,jj)-tmatx( i,j,ii,jj,kk,ll)*smatx(i,j,kk,ll)
                                end do
                            end do
                        end do
                    end do
                    ! 格林操作符计算完成
                end do
            end do
            !---------------------------------------
            do j = 1, Ny
                do i = 1, Nx
                    e11k(i,j)=ematx(i,j,1,1)
                    e22k(i,j)=ematx(i,j,2,2)
                    e12k(i,j)=ematx(i,j,1,2)
                end do
            end do
            ! 通过逆傅里叶变换将结果转换到实数空间
            call ifft_2d(input=e11k, output=e11, Nx=Nx, Ny=Ny, iNxNy = 1.0_wp/(Nx*Ny))
            call ifft_2d(input=e22k, output=e22, Nx=Nx, Ny=Ny, iNxNy = 1.0_wp/(Nx*Ny))
            call ifft_2d(input=e12k, output=e12, Nx=Nx, Ny=Ny, iNxNy = 1.0_wp/(Nx*Ny))
            ! 计算应力
            do j = 1, Ny
                do i = 1, Nx
                    s11(i,j) = c11(i,j)*(ea(1)+e11(i,j)-ei11(i,j)-ed11(i,j))+ &
                               c12(i,j)*(ea(2)+e22(i,j)-ei22(i,j)-ed22(i,j))
                    s22(i,j) = c11(i,j)*(ea(2)+e22(i,j)-ei22(i,j)-ed22(i,j))+ &
                               c12(i,j)*(ea(1)+e11(i,j)-ei11(i,j)-ed11(i,j))
                    s12(i,j) = 2.0*c44(i,j)*(ea(3)+e12(i,j)-ei12(i,j)-ed12(i,j))
                end do
            end do
            ! 检查是否收敛
            sum_stres(:, :) = s11(:, :) + s12(:, :) + s22(:, :)
            normF = norm2(sum_stres)
            if (iter /= 1) conver = abs((normF - old_norm) / old_norm)
            if (conver < tolerence) exit
            old_norm = normF
        end do

        ! 计算应变能
        do j = 1, Ny
            do i = 1, Nx
                et11 = ea(1) + e11(i,j) - ei11(i,j) - ed11(i,j)
                et22 = ea(2) + e22(i,j) - ei22(i,j) - ed22(i,j)
                et12 = ea(3) + e12(i,j) - ei12(i,j) - ed12(i,j)
                delsdc(i,j) = 0.5*(et11*((cp12-cm12)*et22+(cp11-cm11)*et11      &
                              -c12(i,j)*ei0 - c11(i,j)*ei0)-ei0*(c12(i,j)       &
                              *et22+c11(i,j)*et11) + ((cp11-cm11)*et22+(cp12    &
                              -cm12)*et11-c12(i,j)*ei0-c11(i,j)*ei0)*et22-ei0   &
                              *(c11(i,j)*et22+c12(i,j)*et11)+2.0_wp*(cp44-cm44) &
                              *et12**2-4.0_wp*ei0*c44(i,j)*et12)
            end do
        end do
        ! 释放内存
        deallocate(ei11, ei22, ei33, ei12)
        deallocate(c11, c12, c44, sum_stres)
        deallocate(s11k, s12k, s22k)
        deallocate(e11k, e12k, e22k)
        deallocate(smatx, ematx)

    end subroutine solve_elasticity_v1

    pure subroutine green_tensor_2d(tmatx, Nx, Ny, kx, ky, cm11, cm12, cm44, cp11, cp12, cp44)
        ! 函数的输入变量
        real(wp),allocatable, intent(inout) :: tmatx(:, :, :, :, :, :)
        integer, intent(in) :: Nx, Ny
        real(wp), intent(in) :: kx(:), ky(:)
        real(wp), intent(in) :: cm11, cm12, cm44
        real(wp), intent(in) :: cp11, cp12, cp44
        ! 函数的局部变量
        real(wp) :: c11, c12, c44
        real(wp) :: chi, d0
        real(wp) :: rr
        real(wp), allocatable :: omega11(:, :), omega12(:, :), omega22(:, :)
        real(wp) :: gmatx(2, 2), dvect(2)
        integer :: i, j, kk, ll, ii, jj

        ! 分配内存
        if(.not. allocated(tmatx)) allocate(tmatx(Nx, Ny, 2, 2, 2, 2))
        allocate(omega11(Nx, Ny), omega12(Nx, Ny), omega22(Nx, Ny))
        ! 计算复合材料的弹性常数
        c11 = 0.5_wp*(cm11+cp11); c12 = 0.5_wp*(cm12+cp12); c44 = 0.5_wp*(cm44+cp44)
        chi = (c11-c12-2.0_wp*c44)/c44
        ! 开始循环计算
        do j = 1, Ny
            do i = 1, Nx
                rr = kx(i)*kx(i)+ky(j)*ky(j)
                d0 = c11*rr**3+chi*(c11+c12)*rr*(kx(i)**2 * ky(j)**2)
                if (rr < 1.0e-8_wp) d0 = 1.0_wp
                omega11(i, j) = (c44*rr**2 + (c11-c44)*rr*ky(j)**2) / (c44*d0);
                omega22(i, j) = (c44*rr**2 +( c11-c44)*rr*kx(i)**2) / (c44*d0);
                omega12(i, j) = -(c12+c44)*kx(i)*ky(j)*rr/(c44*d0)
            end do
        end do
        ! 计算green tensor
        do j = 1, Ny
            do i = 1, Nx
                gmatx(1,1) = omega11(i,j); gmatx(1,2) = omega12(i,j)
                gmatx(2,1) = omega12(i,j); gmatx(2,2) = omega22(i,j)
                dvect(1) = kx(i);          dvect(2)   = ky(j)
                ! 格林操作符
                do kk = 1, 2
                    do ll = 1, 2
                        do ii = 1, 2
                            do jj = 1, 2
                                tmatx(i,j,kk,ll,ii,jj) = 0.25_wp* (gmatx(ll,ii)*dvect(jj)*dvect(kk) &
                              + gmatx(kk,ii)*dvect(jj)*dvect(ll) + gmatx(ll,jj)*dvect(ii)*dvect(kk) &
                              + gmatx(kk,jj)*dvect(ii)*dvect(ll))
                            end do
                        end do
                    end do
                end do
                ! 格林操作符计算完成
            end do
        end do
        ! 释放局部内存
        deallocate(omega11, omega12, omega22)
    end subroutine

    subroutine write_vtk_file(filename, Nx, Ny, dx, dy, data, step)
        !! 将数组data中的数据写入到VTK文件中
        character(len=*), intent(in) :: filename  !< VTK文件的名字
        integer, intent(in) :: Nx   !< x方向上的节点个数
        integer, intent(in) :: Ny   !< y方向上的节点个数
        integer, intent(in) :: step !< 当前仿真步数
        real(WP), intent(in) :: dx  !< 网格节点在x方向上的间距
        real(WP), intent(in) :: dy  !< 网格节点在y方向上的间距
        real(WP), intent(in) :: data(:, :)  !< 需要写入VTK文件的数据
        character(len=10) :: cstep  !< 字符串形式的step
        integer :: fileid
        integer :: npoint
        integer :: i, j

        write(cstep, "(I0)") step
        open(newunit=fileid, file=filename//trim(adjustl(cstep))//'.vtk', &
        action='write')
        npoint = Nx * Ny
        ! 写入头部注释文字
        write(fileid, "(A)") "# vtk DataFile Version 2.0"
        write(fileid, "(A, I0)") "Timestep: ", step
        write(fileid, "(A)") "ASCII"
        write(fileid, "(A)") "DATASET STRUCTURED_GRID"
        ! 写入节点个数
        write(fileid, "(A, I5, I5, I5)") "DIMENSIONS ", Nx, Ny, 1
        ! 写入节点类型
        write(fileid, "(A, I0, A)") 'POINTS ', npoint, ' float'
        ! 写入节点坐标
        do i = 1, Nx
            do j = 1, Ny
                write(fileid, "(*(F14.6, 1x))") (i-1)*dx, (j-1)*dy, 0.0_WP
            end do
        end do
        ! 写入节点数据
        write(fileid, "(A, I0)") 'POINT_DATA ', npoint
        write(fileid, "(A)") 'SCALARS PHI float 1'
        write(fileid, "(A)") 'LOOKUP_TABLE default'
        ! 这里由于VTK文件是行优先读写, 只能按照行优先的顺序写出文件
        do j = 1, Ny
            do i = 1, Nx
                write(fileid, '(F14.6)') data(j, i)
            end do
        end do
        close(fileid)
    end subroutine write_vtk_file

    pure subroutine prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)
        !! 计算在对应的等间距节点下的傅里叶系数
        real(wp), allocatable, intent(out) :: kx(:), ky(:)
        real(wp), allocatable, intent(out) :: k2(:, :), k4(:, :)
        integer, intent(in) :: Nx, Ny
        real(wp), intent(in) :: dx, dy
        integer :: Nx21, Ny21
        integer :: Nx2, Ny2
        integer :: i, j
        real(wp) :: delkx, delky
        real(wp) :: factor

        ! 分配内存
        allocate(kx(Nx+1), ky(Ny+1), k2(Nx, Ny), k4(Nx, Ny))
        Nx21 = Nx / 2 + 1; Ny21 = Ny / 2 + 1
        Nx2 = Nx + 2;      Ny2  = Ny + 2
        delkx = (2*pi)/(Nx*dx); delky = (2*pi)/(Ny*dy)
        ! 计算kx
        do i = 1, Nx21
            factor = (i-1)*delkx
            kx(i) = factor
            kx(Nx2-i) = -factor
        end do
        ! 计算ky
        do i = 1, Ny21
            factor = (i-1)*delky
            ky(i) = factor
            ky(Ny2-i) = -factor
        end do
        ! 计算K2
        do j = 1, Ny
            do i = 1, Nx
                k2(i, j) = kx(i)*kx(i) + ky(j)*ky(j)
            end do
        end do
        ! 计算K4
        k4 = k2 * k2

    end subroutine

    subroutine fft_2d(input, output, Nx, Ny)
        !! 对输入的二维实数数组进行快速傅里叶变换
        real(wp), intent(inout) :: input(:, :)
        complex(cwp), intent(inout) :: output(:, :)
        integer, intent(in) :: Nx, Ny
        complex(cwp), allocatable :: c_input(:, :)
        type(c_ptr) :: plan

        ! 分配内存
        allocate(c_input(Nx, Ny))
        ! 进行初始化
        c_input = cmplx(input)
        ! 构建正向二维傅里叶变换的计划, 注意C与Fortran交换的时候数组维度是相反的
        plan = fftw_plan_dft_2d(Ny, Nx, c_input, output, FFTW_FORWARD, FFTW_ESTIMATE)
        ! 执行傅里叶变换
        call fftw_execute_dft(plan, c_input, output)
        ! 释放内存
        call fftw_destroy_plan(plan)
        deallocate(c_input)

    end subroutine fft_2d

    subroutine ifft_2d(input, output, Nx, Ny, iNxNy)
        !! 对输入的二维实数数组进行快速傅里叶变换
        complex(cwp), intent(inout) :: input(:, :)
        real(wp), intent(inout) :: output(:, :)
        integer, intent(in) :: Nx, Ny
        real(wp), intent(in) :: iNxNy
        complex(cwp), allocatable :: c_output(:, :)
        type(c_ptr) :: plan

        ! 分配内存
        allocate(c_output(Nx, Ny))
        ! 构建正向二维傅里叶变换的计划, 注意C与Fortran交换的时候数组维度是相反的
        plan = fftw_plan_dft_2d(Ny, Nx, input, c_output, FFTW_BACKWARD, FFTW_ESTIMATE)
        ! 执行傅里叶变换
        call fftw_execute_dft(plan, input, c_output)
        ! 将虚数转为实数
        output = real(c_output) * iNxNy
        ! 释放内存
        call fftw_destroy_plan(plan)
        deallocate(c_output)

    end subroutine ifft_2d
end module case_9_V1

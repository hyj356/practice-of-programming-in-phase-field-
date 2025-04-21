module case_6_V1
  use iso_c_binding, only: c_double, c_int, c_double_complex
  use iso_fortran_env, only:real32, real64, output_unit
  use FFTW, only: fftw_plan_dft_2d, fftw_execute_dft, fftw_destroy_plan, &
                  FFTW_FORWARD, FFTW_BACKWARD, &
                  c_ptr, FFTW_ESTIMATE
  implicit none
  private
  public :: fft_2d, ifft_2d, micro_ch_pre, write_vtk_file, &
            prepare_fft, free_energ_ch_v1
  integer, parameter, public :: wp = c_double
  integer, parameter, public :: cwp = c_double_complex
  integer, parameter, public :: ci = c_int
  integer, parameter, public :: stdout = output_unit
  real(WP), parameter, public :: pi = 4.0_wp * atan(1.0_WP)
  real(wp), parameter, public :: coefA = 1.0_wp

contains
    pure function free_energ_ch_v1(con, i, j) result(ret)
        !! 计算化学能函数f(c)的导数, 在书上的案例中, f(c) = A*c^2*(1-c)^2
        !! 那么, f'(c) = 2*A*[c*(1-c)^2 - c^2*(1-c)] = 2*A*c*(1-c)
        real(WP), intent(in) :: con(:, :) !< 浓度数组, 形状为con(Nx, Ny)
        integer, intent(in) :: i, j       !< 指明计算几行几列的节点
        real(WP) :: ret
        ret = 2.0_WP * coefA * con(i, j) * (1.0_WP - con(i, j)) * (1.0_WP - 2*con(i, j))

    end function free_energ_ch_v1

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

    subroutine micro_ch_pre(con, Nx, Ny, c0)
        !! 初始化合金的浓度数组
        real(WP), allocatable, intent(out) :: con(:, :)
        integer, intent(in) :: Nx, Ny
        real(wp), intent(in) :: c0
        integer :: i, j
        real(wp), allocatable :: rng(:, :)

        ! 分配内存
        if (.not. allocated(con)) allocate(con(Nx, Ny))
        allocate(rng(Nx, Ny))
        ! 初始化随机数种子
        call random_seed()
        call random_number(rng)
        rng = rng - 0.5_wp
        ! 开始初始化浓度, 噪音为±0.02
        do j = 1, Ny
            do i = 1, Nx
                con(i, j) = c0 + 0.02_wp * rng(i, j)
            end do
        end do
    end subroutine micro_ch_pre

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
end module

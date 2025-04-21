module case1
  use iso_fortran_env, only: WP => real64, stdout => output_unit
  implicit none
  private
  public :: WP, stdout, write_con, micro_ch_pre, free_energ_ch_v1, &
            write_vtk_file, calculate_energy, free_energ_ch_v2
  
  real(WP), parameter, public :: dt = 1.0e-2_WP  !< 单步步长
  real(WP), parameter, public :: mobility = 1.0_WP    !< 对应书上公式(4.16)中的M
  real(WP), parameter, public :: grad_coeff = 0.5_WP  !< 对应书上公式(4.17)中的κ
  real(WP), parameter, public :: dx = 1.0_WP !< x方向上相邻节点的距离为1
  real(WP), parameter, public :: dy = 1.0_WP !< y方向上相邻节点的距离为1
  real(WP), parameter, private :: A = 1.0_WP

contains 

  subroutine write_con(filename, con)
    character(len=*), intent(in) :: filename
    real(WP), intent(in) :: con(:, :)
    integer :: Ny
    integer :: fileid
    integer :: j

    Ny = size(con, dim=2)
    open(newunit=fileid, file=filename, action='write')
    do j = 1, Ny
        write(fileid, *) con(:, j)
    end do
    close(fileid)
  end subroutine write_con

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

  pure function calculate_energy(Nx, Ny, con) result(energy)
    real(WP), intent(in) :: con(:, :) !< 浓度数组, 形状为con(Nx, Ny)
    integer, intent(in) :: Nx   !< x方向上的节点个数
    integer, intent(in) :: Ny   !< y方向上的节点个数
    real(WP) :: energy
    integer :: i, j
    integer :: ip, jp

    energy = 0.0_WP
    do j = 1, Ny - 1
      jp = j + 1
      do i = 1, Nx - 1
        ip = i + 1
        energy = energy + con(i, j)*con(i, j) * (1.0_WP - con(i, j))*(1.0_WP - con(i, j)) + &
                 0.5_WP * grad_coeff * ( (con(ip, j) - con(i, j))**2 + (con(i, jp) - con(i, j))**2 )
      end do
    end do

  end function calculate_energy

  pure function free_energ_ch_v1(con, i, j) result(ret)
    !! 计算化学能函数f(c)的导数, 在书上的案例中, f(c) = A*c^2*(1-c)^2
    !! 那么, f'(c) = 2*A*[c*(1-c)^2 - c^2*(1-c)] = 2*A*c*(1-c)
    real(WP), intent(in) :: con(:, :) !< 浓度数组, 形状为con(Nx, Ny)
    integer, intent(in) :: i, j       !< 指明计算几行几列的节点
    real(WP) :: ret
    ret = 2.0_WP * A * con(i, j) * (1.0_WP - con(i, j)) * (1.0_WP - 2*con(i, j))
  
  end function free_energ_ch_v1

  pure elemental function free_energ_ch_v2(con) result(ret)
  !! 计算化学能函数f(c)的导数, 在书上的案例中, f(c) = A*c^2*(1-c)^2
  !! 那么, f'(c) = 2*A*[c*(1-c)^2 - c^2*(1-c)] = 2*A*c*(1-c)*(1-2*c)
    real(WP), intent(in) :: con !< 浓度数组, 形状为con(Nx, Ny)
    real(WP) :: ret
    ret = 2.0_WP * A * con * (1.0_WP - con) * (1.0_WP - 2*con)

  end function free_energ_ch_v2

  subroutine write_vtk_file(filename, Nx, Ny, data, step)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: Nx   !< x方向上的节点个数
    integer, intent(in) :: Ny   !< y方向上的节点个数
    integer, intent(in) :: step !< 当前仿真步数
    real(WP), intent(in) :: data(:, :)  !< 需要写入VTK文件的数据
    character(len=10) :: cstep  !< 字符串形式的step
    integer :: fileid
    integer :: npoint
    integer :: i, j
    
    write(cstep, "(I0)") step
    open(newunit=fileid, file=filename//trim(adjustl(cstep))//'.vtk', action='write')
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
    write(fileid, "(A)") 'SCALARS CON float 1'
    write(fileid, "(A)") 'LOOKUP_TABLE default'
    do j = 1, Ny
      do i = 1, Nx
        write(fileid, '(F14.6)') data(i, j)
      end do
    end do
    close(fileid)
  end subroutine write_vtk_file

end module
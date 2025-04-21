module case_3_V2
  use iso_fortran_env, only: real32, real64, output_unit
  implicit none

  private
  public :: write_eta, write_vtk_file, micro_sint_pre, free_energ_sint_v1, &
            laplacian_2, laplacian_4, free_energ_sint_v2
  
  integer, parameter, public :: WP = real32
  integer, parameter, public :: stdout = output_unit
  real(WP), parameter, private :: A = 16.0_WP
  real(WP), parameter, private :: B = 1.0_WP
  real(WP), parameter, public :: coefm = 5.0_WP !< 对应κ_ρ
  real(WP), parameter, public :: coefk = 2.0_WP
  real(WP), parameter, public :: coefl = 5.0_WP
  real(WP), parameter, public :: Dvol = 0.04_WP
  real(WP), parameter, public :: Dvap = 0.002_WP
  real(WP), parameter, public :: Dsur = 16.0_WP
  real(WP), parameter, public :: Dgrb = 1.6_WP

contains
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
    write(fileid, "(A)") 'SCALARS ETA float 1'
    write(fileid, "(A)") 'LOOKUP_TABLE default'
    do j = 1, Ny
      do i = 1, Nx
        write(fileid, '(F14.6)') data(i, j)
      end do
    end do
    close(fileid)
  end subroutine write_vtk_file

  subroutine write_eta(filename, eta)
    !! 将数组eta中的数据写入到文件中, 以便于python进行可视化
    character(len=*), intent(in) :: filename
    real(WP), intent(in) :: eta(:, :)
    integer :: Ny
    integer :: fileid
    integer :: j

    Ny = size(eta, dim=2)
    open(newunit=fileid, file=filename, action='write')
    do j = 1, Ny
        write(fileid, *) eta(:, j)
    end do
    close(fileid)
  end subroutine write_eta

  subroutine micro_sint_pre(Nx, Ny, npart, iflag, etas, con)
    !! 根据参数iflag, 用不同的方式初始化固态烧结的序参数和密度数组
    integer, intent(in) :: Nx     !< x方向上的节点个数
    integer, intent(in) :: Ny     !< y方向上的节点个数
    integer, intent(in) :: npart  !< 模型中有几个颗粒
    integer, intent(in) :: iflag  !< 按哪种方式初始化模型
    real(WP), allocatable, intent(out) :: etas(:, :, :)  !< 各晶粒的序参数的数组
    real(WP), allocatable, intent(out) :: con(:, :)      !< 浓度数组
    real(WP), dimension(9) :: xc, yc
    real(WP) :: Rx, R, R1, R2
    real(WP) :: x1, y1, y2
    real(WP) :: xx1, xx2 
    integer :: i, j, ipart !< 循环变量
    if(iflag == 1) then
      allocate(etas(Nx, Ny, npart), con(Nx, Ny), source=0.0_WP)
    end if
    ! 如果需要初始化的晶粒不是2个, 那么就初始化npart个晶粒
    if(npart /= 2) then
      R = 10.0_WP
      xc(1) = 29.0_WP; yc(1) = 50.0_WP
      xc(2) = 50.0_WP; yc(2) = 50.0_WP
      xc(3) = 71.0_WP; yc(3) = 50.0_WP
      xc(4) = 50.0_WP; yc(4) = 29.0_WP
      xc(5) = 50.0_WP; yc(5) = 71.0_WP
      xc(6) = 39.0_WP; yc(6) = 39.0_WP
      xc(7) = 61.0_WP; yc(7) = 39.0_WP
      xc(8) = 39.0_WP; yc(8) = 61.0_WP
      xc(9) = 61.0_WP; yc(9) = 61.0_WP
      do ipart = 1, npart
        Rx = merge(0.5_WP*R, R, ipart > 5)
        do j = 1, Ny
          do i = 1, Nx
            xx1 = sqrt( (i-xc(ipart))**2 + (j-yc(ipart))**2 )
            if (xx1 <= Rx) then
              con(i, j) = 0.9999999_WP
              etas(i, j, ipart) = 0.9999999_WP
            end if
          end do
        end do
      end do
    end if
    ! 如果npart等于2, 说明需要初始化两个晶粒
    if (npart == 2) then
      R1 = 20.0_WP; R2 = 0.5_WP * R1
      x1 = Nx/2; y1 = 40.0_WP; y2 = 70.0_WP
      do j = 1, Ny
        do i = 1, Nx
          xx1 = sqrt( (i-x1)**2 + (j-y1)**2 )
          xx2 = sqrt( (i-x1)**2 + (j-y2)**2 )
          if (xx1 <= R1) then
            con(i, j) = 0.9999999_WP 
            etas(i, j, 1) = 0.9999999_WP
          end if
          if(xx2 <= R2) then
            con(i, j) = 0.9999999_WP 
            etas(i, j, 1) = 0.0_WP
            etas(i, j, 2) = 0.9999999_WP
          end if
        end do
      end do
    end if
  end subroutine micro_sint_pre
  
  pure function free_energ_sint_v1(i, j, con, eta, etas, npart, iflag) result(dfunc)
    !! 此函数根据iflag等于1还是等于2, 
    !! 用于计算自由能F对浓度con和序参数eta的导数值 
    integer, intent(in) :: i, j, npart, iflag
    real(WP), intent(in) :: con(:, :), eta(:, :), etas(:, :, :)
    real(WP), dimension(2) :: dfunc !< 第一个等于dfdcon, 第二个等于dfdeta
    real(WP) :: sum2, sum3
    integer :: ipart

    dfunc(1:2) = 0.0_WP
    ! 计算dfdcon
    if (iflag == 1) then
      sum2 = 0.0_WP; sum3 = 0.0_WP
      do ipart = 1, npart
        sum2 = sum2 + etas(i, j, ipart)**2
        sum3 = sum3 + etas(i, j, ipart)**3
      end do
      dfunc(iflag) = B*(2*con(i, j) + 4*sum3-6*sum2) - &
                   2*A*con(i, j)*con(i, j)*(1.0_WP - con(i, j)) + &
                   2*A*con(i, j)*(1.0_WP - con(i, j))*(1.0_WP - con(i, j))
    end if
    ! 计算dfdeta
    if (iflag == 2) then
      sum2 = 0.0_WP
      do ipart = 1, npart
        sum2 = sum2 + etas(i, j, ipart)**2
      end do
      dfunc(iflag) = B*(-12*eta(i, j)**2 * (2.0_WP - con(i, j)) + &
                         12*eta(i, j)*(1.0_WP - con(i, j))      + &
                         12*eta(i, j)*sum2)
    end if

  end function free_energ_sint_v1

  pure subroutine free_energ_sint_v2(Nx, Ny, con, eta, etas, dfunc, iflag)
    !! 此函数根据iflag等于1还是等于2, 
    !! 用于计算自由能F对浓度con和序参数eta的导数值 
    integer, intent(in) :: Nx, Ny, iflag
    real(WP), intent(in) :: con(:, :), eta(:, :), etas(:, :, :)
    real(WP), allocatable, intent(inout) :: dfunc(:, :, :) !< 第一个等于dfdcon, 第二个等于dfdeta
    real(WP), allocatable :: sum2(:, :), sum3(:, :)

    if (.not. allocated(dfunc)) allocate(dfunc(Nx, Ny, 2))
    
    ! 计算dfdcon的数组形式
    if (iflag == 1) then
      sum2 = sum(etas(:, :, :)*etas(:, :, :), dim=3)
      sum3 = sum(etas(:, :, :)*etas(:, :, :)*etas(:, :, :), dim=3)
      dfunc(:, :, iflag) = B*(2*con(:, :) + 4*sum3(:, :)-6*sum2(:, :)) - &
                   2*A*con(:, :)*con(:, :)*(1.0_WP - con(:, :)) + &
                   2*A*con(:, :)*(1.0_WP - con(:, :))*(1.0_WP - con(:, :))
      deallocate(sum2, sum3)
    end if
    ! 计算dfdeta的数组形式
    if (iflag == 2) then
      sum2 = sum(etas(:, :, :)*etas(:, :, :), dim=3)
      dfunc(:, :, iflag) = B*(-12*eta(:, :)**2 * (2.0_WP - con(:, :)) + &
                         12*eta(:, :)*(1.0_WP - con(:, :))      + &
                         12*eta(:, :)*sum2(:, :))
      deallocate(sum2)
    end if

  end subroutine free_energ_sint_v2

  pure function laplacian_2(data, i, j, Nx, Ny, idxdy) result(ret)
    !! 对二维数组的第(i, j)个元素进行拉普拉斯算子操作
    real(WP), intent(in) :: data(:, :)  !< 对应的二维数据数组
    real(WP), intent(in) :: idxdy       !< 其值应该等于dx*dy的倒数
    integer, intent(in) :: i, j         !< 下标
    integer, intent(in) :: Nx, Ny       !< 对应二维数组的第一个和第二个维度大小
    real(WP) :: ret
    integer :: ip, im, jp, jm

    ip = i + 1; im = i - 1
    jp = j + 1; jm = j - 1
    ! 对下标数组进行修正
    if (ip == Nx + 1) ip = 1
    if (im == 0) im = Nx
    if (jp == Ny + 1) jp = 1
    if (jm == 0) jm = Ny

    ret = (data(ip, j) + data(im, j) + data(i, jp) + data(i, jm) - &
    4.0_WP*data(i, j)) * idxdy

  end function laplacian_2

  pure function laplacian_4(data, i, j, Nx, Ny, idxdy) result(ret)
    !! 对二维数组的第(i, j)个元素进行两次拉普拉斯算子操作
    real(WP), intent(in) :: data(:, :)  !< 对应的二维数据数组
    real(WP), intent(in) :: idxdy       !< 其值应该等于dx*dy的倒数
    integer, intent(in) :: i, j         !< 下标
    integer, intent(in) :: Nx, Ny       !< 对应二维数组的第一个和第二个维度大小
    real(WP) :: ret
    integer :: ip, im, jp, jm
    real(WP) :: lap_ip, lap_im, lap_jp, lap_jm, lap_ce

    ip = i + 1; im = i - 1
    jp = j + 1; jm = j - 1
    ! 对下标数组进行修正
    if (ip == Nx + 1) ip = 1
    if (im == 0) im = Nx
    if (jp == Ny + 1) jp = 1
    if (jm == 0) jm = Ny
    ! 计算5个节点的拉普拉斯算子值
    lap_ip = laplacian_2(data, ip, j, Nx, Ny, idxdy)
    lap_im = laplacian_2(data, im, j, Nx, Ny, idxdy)
    lap_jp = laplacian_2(data, i, jp, Nx, Ny, idxdy)
    lap_jm = laplacian_2(data, i, jm, Nx, Ny, idxdy)
    lap_ce = laplacian_2(data, i, j, Nx, Ny, idxdy)
    ! 计算最终的结果
    ret = (lap_ip + lap_im + lap_jp + lap_jm - 4*lap_ce) * idxdy

  end function laplacian_4
end module case_3_V2
module case_5_V2
  use iso_fortran_env, only: real32, real64, output_unit
  implicit none

  private
  public :: write_eta, write_vtk_file, laplacian_2, micro_poly_cell, free_energy_v1, &
            laplace_2, gradient_mat, free_energy_v2, gradient_matdx, gradient_matdy

  integer, parameter, public :: WP = real32
  integer, parameter, public :: stdout = output_unit
  real(WP), parameter, public :: lamda = 7.0_WP
  real(WP), parameter, public :: mu = 40.0_WP
  real(WP), parameter, public :: kisa = 1.5e3_WP
  real(WP), parameter, private :: velc = 0.2_WP
  
contains

  subroutine micro_poly_cell(Nx, Ny, ncell, R, phis, vac, nccel, ccell)
    !! 初始化细胞的序参数和自推进值数组
    integer, intent(in) :: Nx, Ny !< x和y方向上的节点数量
    real(WP), intent(in) :: R     !< cell的半径
    integer, intent(inout) :: ncell   !< 输入的是期望的cell数量, 输出的是实际的cell数量
    integer, intent(out) :: nccel !< 软细胞(soft cell)的数量
    integer, intent(out), allocatable :: ccell(:) !< soft cell的细胞数量
    real(WP), intent(out), allocatable :: phis(:, :, :)  !< 细胞的序参数数组
    real(WP), intent(out), allocatable :: vac(:)  !< 细胞的自推进值数组
    ! 子程序中的局部变量
    integer :: iflag
    integer :: iter
    integer :: i, j
    real(WP) :: R2, Rsq
    real(WP) :: xmin, ymin, xmax, ymax
    real(WP) :: xdist
    real(WP) :: xnc, ync
    real(WP) :: ix
    real(WP) :: irand, jrand
    real(WP), allocatable :: xc(:), yc(:)

    ! 分配内存并初始化
    if (.not. allocated(phis)) allocate(phis(Nx, Ny, ncell), source=0.0_WP)
    ! 计算细胞的半径
    ! 如果期望产生80个cell
    if (ncell == 80) then
      R2 = 2.0_WP*R; Rsq = R*R; ncell = 0
      xmin = 0.0_WP; ymin = 0.0_WP
      xmax = real(Nx, kind=WP); ymax = real(Ny, kind=WP)
      allocate(xc(80), yc(80), source=0.0_WP)
      call random_seed()
      do iter = 1, 500000
        call random_number(irand)
        call random_number(jrand)
        xnc = Nx*irand; ync = Ny*jrand
        iflag = 1
        ! 确保生成的圆会整个落在仿真区域内
        if((xnc - R) < xmin .or. (xnc + R) > xmax) iflag = 0
        if((ync - R) < ymin .or. (ync + R) > ymax) iflag = 0
        ! 通过上述判定之后, 接下来判断新的圆心是否会与旧的圆产生重叠
        if (iflag == 1) then
          do i = 1, ncell
            xdist = sqrt( (xc(i) - xnc)*(xc(i) - xnc) + (yc(i) - ync)*(yc(i) - ync) )
            if (xdist <= 1.6_WP*R) iflag = 0
          end do
        end if
        ! 通过上述两个判定之后, 确定新生成的点符合要求, 将其加入数组xc和yc中
        if (iflag == 1) then
          ncell = ncell + 1
          xc(ncell) = xnc; yc(ncell) = ync
          ! 遍历所有网格节点, 如果距离新生成的点的距离小于截断半径, 将其序参数设为1
          do j = 1, Ny
            do i = 1, Nx
              if( (i-xnc)*(i-xnc) + (j-ync)*(j-ync) < Rsq ) then
                phis(i, j, ncell) = 0.99999_WP
              end if
            end do
          end do
        end if
        ! 如果cell数量达到80个, 跳出循环
        if (ncell == 80) exit
      end do
      ! 分配内存
      if (.not. allocated(vac)) allocate(vac(ncell))
      ! 初始化细胞的自推进值
      do i = 1, ncell
        call random_number(ix)
        vac(i) = merge(-velc, velc, ix <= 0.5)
      end do
      ! 打印迭代结果
      write(stdout, '(A, I0)') "Iteration done: ", iter
      write(stdout, '(I0, A)') ncell, " of cell created."
      ! 初始化软细胞下标
      nccel = 5
      if (.not. allocated(ccell)) allocate(ccell(nccel))
      ! 这里我们选择第32, 11, 16, 21, 46个细胞为软细胞(soft cell)
      ccell(1) = 32
      ccell(2) = 11
      ccell(3) = 16
      ccell(4) = 21
      ccell(5) = 46
    end if

    ! 对于只有两个细胞的模拟
    if (ncell == 2) then
      ! 对于两个细胞的模拟, 不引入软细胞(soft cell)
      nccel = 0
      ! 分配内存
      if (.not. allocated(ccell)) allocate(ccell(1))
      if (.not. allocated(vac)) allocate(vac(ncell))
      vac(1) = 0.5_WP; vac(2) = -0.5_WP
      ccell(1) = 10000
      R2 = R*R
      allocate(xc(2), yc(2))
      ! 计算两个cell的中心坐标
      xc(1) = Nx/2 - 1.25*R; yc(1) = Ny / 2
      xc(2) = Nx/2 + 1.25*R; yc(2) = Ny / 2
      ! 遍历整个数组, 将位于细胞内的节点的序参数值改为1.0
      do j = 1, Ny
        do i = 1, Nx
          xdist = (xc(1) - i)*(xc(1) - i) + (yc(1) - j)*(yc(1) - j) 
          if ( xdist < R2 ) phis(i, j, 1) = 0.99999_WP
          xdist = (xc(2) - i)*(xc(2) - i) + (yc(2) - j)*(yc(2) - j) 
          if ( xdist < R2 ) phis(i, j, 2) = 0.99999_WP
        end do
      end do
      write(stdout, *) "Two cell generated in simulation region."
    end if

    ! 释放局部动态数组变量的内存
    if(allocated(xc)) deallocate(xc)
    if(allocated(yc)) deallocate(yc)

  end subroutine

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

  pure function free_energy_v1(i, j, icell, ncell, gamma_, kappa, phi, phis) result(dfdphi)
    integer, intent(in) :: i, j, icell, ncell
    real(WP), intent(in) :: gamma_, kappa
    real(WP), intent(in) :: phi(:, :), phis(:, :, :)
    ! 函数内部的局部变量
    real(WP) :: sum_phi
    integer :: jcell
    ! 该函数计算的结果变量
    real(WP) :: dfdphi

    sum_phi = 0.0_WP
    do jcell = 1, ncell
      if (jcell /= icell) then
        sum_phi  = sum_phi + phi(i, j)*phis(i, j, jcell)*phis(i, j, jcell)
      end if
    end do
    dfdphi = gamma_ * phi(i, j) * ( 1.0_WP - phi(i, j) ) * ( 1.0_WP - 2.0_WP*phi(i, j) ) + &
             2.0_WP * kappa * phi(i, j) * sum_phi

  end function free_energy_v1

  pure subroutine free_energy_v2(Nx, Ny, dfdphi, icell, ncell, gamma_, kappa, phi, phis)
    integer, intent(in) :: icell, ncell
    integer, intent(in) :: Nx, Ny
    real(WP), intent(in) :: gamma_, kappa
    real(WP), intent(in) :: phi(:, :), phis(:, :, :)
    real(WP), allocatable, intent(inout) :: dfdphi(:, :)
    ! 函数内部的局部变量
    real(WP), allocatable :: sum_phi(:, :)
    integer :: jcell

    ! 为该函数计算的结果变量分配内存
    if(.not. allocated(dfdphi)) allocate(dfdphi(Nx, Ny))
    allocate(sum_phi(Nx, Ny), source=0.0)
    sum_phi = 0.0_WP
    do jcell = 1, ncell
      if (jcell /= icell) then
        sum_phi(:, :)  = sum_phi(:, :) + phi(:, :)*phis(:, :, jcell)*phis(:, :, jcell)
      end if
    end do
    dfdphi(:, :) = gamma_ * phi(:, :) * ( 1.0_WP - phi(:, :) ) * ( 1.0_WP - 2.0_WP*phi(:, :) ) + &
            2.0_WP * kappa * phi(:, :) * sum_phi(:, :)
    deallocate(sum_phi)
  end subroutine free_energy_v2

  pure subroutine gradient_mat(Nx, Ny, dx, dy, matx, matdx, matdy)
    !! 此子程序用于计算输入二维矩阵matx在x和y方向上的导数矩阵
    integer, intent(in) :: Nx, Ny
    real(WP), intent(in) :: dx, dy
    real(WP), intent(in) :: matx(:, :)
    real(WP), allocatable, intent(inout) :: matdx(:, :)
    real(WP), allocatable, intent(inout) :: matdy(:, :)

    ! 分配内存
    if (.not. allocated(matdx)) allocate(matdx(Nx, Ny))
    if (.not. allocated(matdy)) allocate(matdy(Nx, Ny))

    ! 计算x方向上的梯度矩阵
    matdx(:, 2:Ny-1) = 0.5_WP*(matx(:, 3:Ny) - matx(:, 1:Ny-2))

    ! 计算y方向上的梯度矩阵
    matdy(2:Nx-1, :) = 0.5_WP*(matx(3:Nx, :) - matx(1:Nx-2, :))

    ! 根据周期性边界条件进行修正x方向上的梯度矩阵
    matdx(:, 1) = 0.5_WP*(matx(:, 2) - matx(:, Ny))
    matdx(:, Ny) = 0.5_WP*(matx(:, 1) - matx(:, Ny-1))

    ! 根据周期性边界条件进行修正y方向上的梯度矩阵
    matdy(1, :) = 0.5_WP*(matx(2, :) - matx(Nx, :))
    matdy(Nx, :) = 0.5_WP*(matx(1, :) - matx(Nx-1, :))

    ! 考虑网格大小, 加入系数修正
    matdx = 2.0_WP * matdx / dx
    matdy = 2.0_WP * matdy / dy
  end subroutine gradient_mat

  pure subroutine gradient_matdx(Nx, Ny, dx, dy, matx, matdx)
    !! 此子程序用于计算输入二维矩阵matx在x方向上的导数矩阵
    integer, intent(in) :: Nx, Ny
    real(WP), intent(in) :: dx, dy
    real(WP), intent(in) :: matx(:, :)
    real(WP), allocatable, intent(inout) :: matdx(:, :)

    ! 分配内存
    if (.not. allocated(matdx)) allocate(matdx(Nx, Ny))
    ! 计算x方向上的梯度矩阵
    matdx(:, 2:Ny-1) = 0.5_WP*(matx(:, 3:Ny) - matx(:, 1:Ny-2))
    ! 根据周期性边界条件进行修正x方向上的梯度矩阵
    matdx(:, 1) = 0.5_WP*(matx(:, 2) - matx(:, Ny))
    matdx(:, Ny) = 0.5_WP*(matx(:, 1) - matx(:, Ny-1))
    ! 考虑网格大小, 加入系数修正
    matdx = 2.0_WP * matdx / dx
  end subroutine gradient_matdx

  pure subroutine gradient_matdy(Nx, Ny, dx, dy, matx, matdy)
    !! 此子程序用于计算输入二维矩阵matx在x方向上的导数矩阵
    integer, intent(in) :: Nx, Ny
    real(WP), intent(in) :: dx, dy
    real(WP), intent(in) :: matx(:, :)
    real(WP), allocatable, intent(inout) :: matdy(:, :)

    ! 计算y方向上的梯度矩阵
    matdy(2:Nx-1, :) = 0.5_WP*(matx(3:Nx, :) - matx(1:Nx-2, :))
    ! 根据周期性边界条件进行修正y方向上的梯度矩阵
    matdy(1, :) = 0.5_WP*(matx(2, :) - matx(Nx, :))
    matdy(Nx, :) = 0.5_WP*(matx(1, :) - matx(Nx-1, :))
    ! 考虑网格大小, 加入系数修正
    matdy = 2.0_WP * matdy / dy

  end subroutine gradient_matdy

  pure subroutine laplace_2(Nx, Ny, idxdy, data, lap_data)
  !! 对二维数组data的所有元素进行拉普拉斯算子操作
  integer, intent(in) :: Nx, Ny       !< 对应二维数组的第一个和第二个维度大小
  real(WP), intent(in) :: idxdy       !< 其值应该等于dx*dy的倒数
  real(WP), intent(in) :: data(:, :)  !< (Nx, Ny), 对应的二维数据数组
  real(WP), allocatable, intent(inout) :: lap_data(:, :) !< (Nx, Ny), 进行Laplace算子之后的二维数据数组

    if(.not. allocated(lap_data)) allocate(lap_data(Nx, Ny))
    ! 计算中间格点的拉普拉斯算子值
    lap_data(2:Nx-1, 2:Ny-1) = ( data(1:Nx-2, 2:Ny-1) + data(3:Nx, 2:Ny-1) + &
                                data(2:Nx-1, 1:Ny-2) + data(2:Nx-1, 3:Ny) - &
                                4.0_WP * data(2:Nx-1, 2:Ny-1)&
                                )
    ! 考虑周期性边界条件, 对x = 1, 即第1列的情况进行修正
    lap_data(2:Nx-1, 1) = ( data(1:Nx-2, 1) + data(3:Nx, 1) + &
                            data(2:Nx-1, 2) + data(2:Nx-1, Ny) - &
                            4.0_WP * data(2:Nx-1, 1) &
                          )     
    ! 考虑周期性边界条件, 对x = Ny, 即第Ny列的情况进行修正
    lap_data(2:Nx-1, Ny) = ( data(1:Nx-2, Ny) + data(3:Nx, Ny) + &
                            data(2:Nx-1, 1) + data(2:Nx-1, Ny-1) - &
                            4.0_WP * data(2:Nx-1, Ny) &
                            )
    ! 考虑周期性边界条件, 对y = 1, 即第1行的情况进行修正 
    lap_data(1, 2:Ny-1) = ( data(1, 1:Ny-2) + data(1, 3:Ny) + &
                            data(2, 2:Ny-1) + data(Nx, 2:Ny-1) - &
                            4.0_WP * data(1, 2:Ny-1) &
                          ) 
    ! 考虑周期性边界条件, 对y = Nx, 即第Nx行的情况进行修正 
    lap_data(Nx, 2:Ny-1) = ( data(Nx, 1:Ny-2) + data(Nx, 3:Ny) + &
                            data(1, 2:Ny-1) + data(Nx-1, 2:Ny-1) - &
                            4.0_WP * data(Nx, 2:Ny-1) &
                          ) 
    ! 考虑周期性边界, 对上下左右共4个角点进行修正
    lap_data(1, 1) = ( data(1, 2) + data(2, 1) + data(Nx, 1) + data(1, Ny) - 4.0_WP*data(1,1))
    lap_data(1, Ny) = ( data(Nx, Ny) + data(2, Ny) + data(1, 1) + data(1, Ny-1) - 4.0_WP*data(1, Ny) )
    lap_data(Nx, Ny) = ( data(Nx-1, Ny) + data(1, Ny) + data(Nx, 1) + data(Nx, Ny-1) - 4.0_WP*data(Nx, Ny) )
    lap_data(Nx, 1) = ( data(Nx-1, 1) + data(1, 1) + data(Nx, 2) + data(Nx, Ny) - 4.0_WP*data(Nx, 1) )
    ! 最后除以(dx*dy)
    lap_data(:, :) = lap_data(:, :) * idxdy
  end subroutine laplace_2
end module case_5_V2
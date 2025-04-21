module case_4_V2
  use iso_fortran_env, only: real32, real64, output_unit
  implicit none

  private
  public :: gradient_mat, nucleus, laplacian_2, write_vtk_file, write_eta, &
  gradient_matdx, gradient_matdy, laplace_2
  ! 基础常用常数
  integer, parameter, public :: WP = real32
  integer, parameter, public :: stdout = output_unit
  real(WP), parameter, public :: pi = 4.0_WP*atan(1.0_WP)
  real(WP), parameter, public :: ipi = 1.0_WP/pi  !< pi的倒数
  ! 与仿真材料相关的常数
  real(WP), parameter, public :: tau = 0.0003_WP
  real(WP), parameter, public :: itau = 1.0_WP/tau  !< tau的倒数
  real(WP), parameter, public :: epsilonb = 0.01_WP
  real(WP), parameter, public :: mu = 1.0_WP
  real(WP), parameter, public :: kappa = 1.8_WP
  real(WP), parameter, public :: delta = 0.02_WP
  real(WP), parameter, public :: aniso = 4.0_WP !< 各向异性系数的模态数, 立方晶格为4, 六方晶格为6
  real(WP), parameter, public :: alpha = 0.9_WP
  real(WP), parameter, public :: gamma_ = 10.0_WP
  real(WP), parameter, public :: Teq = 1.0_WP
  real(WP), parameter, public :: theta0 = 0.2_WP
  

contains 

  pure subroutine nucleus(Nx, Ny, seed, phi, tempr)
    !! 初始化序参数数组phi, 给予一个初始微小的核
    integer, intent(in) :: Nx, Ny
    real(WP), intent(in) :: seed
    real(WP), intent(out), allocatable :: phi(:, :)
    real(WP), intent(out), allocatable :: tempr(:, :)
    integer :: i, j
    integer :: half_Nx, half_Ny

    if(.not. allocated(phi)) allocate(phi(Nx, Ny), source=0.0_WP)
    if(.not. allocated(tempr)) allocate(tempr(Nx, Ny), source=0.0_WP)

    half_Nx = Nx / 2
    half_Ny = Ny / 2

    do j = 1, Ny
      do i = 1, Nx
        if((i - half_Nx)*(i - half_Nx) + (j - half_Ny)*(j - half_Ny) < seed) then
          phi(i, j) = 1.0_WP
        end if
      end do
    end do

  end subroutine nucleus

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

  pure function laplacian_2(data, i, j, Nx, Ny, idxdy) result(ret)
    !! 对二维数组data的第(i, j)个元素进行拉普拉斯算子操作
    real(WP), intent(in) :: data(:, :)  !< (Nx, Ny), 对应的二维数据数组
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

end module case_4_V2
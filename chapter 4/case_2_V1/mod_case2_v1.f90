module case_2_v1
  use iso_fortran_env, only: WP => real32, stdout => output_unit
  implicit none
  public :: WP, stdout, write_vtk_file, init_grain_micro, &
            free_energ_fd_ca_v1, write_eta

  real(WP), public, parameter :: mobil = 5.0_WP
  real(WP), public, parameter :: grcoeff = 0.1_WP
  real(WP), public, parameter :: A = 1.0_WP
  real(WP), public, parameter :: B = 1.0_WP
  
contains
  subroutine write_vtk_file(filename, Nx, Ny, dx, dy, data, step)
    character(len=*), intent(in) :: filename
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

  subroutine init_grain_micro(etas, ngrain, glist, Nx, Ny, dx, dy, iflag)
    real(WP), allocatable, intent(out) :: etas(:, :, :)
    integer, intent(out) :: ngrain !< 晶粒个数
    integer, intent(inout), allocatable :: glist(:) !< 晶粒列表, 1表示晶粒存在
    integer, intent(in) :: Nx       !< x方向上的节点个数
    integer, intent(in) :: Ny       !< y方向上的节点个数
    real(WP), intent(in) :: dx      !< x方向上的节点间距
    real(WP), intent(in) :: dy      !< y方向上的节点间距
    integer, intent(in) :: iflag    !< 若为1, 则生成二维双晶模型, 若为2, 则生成二维多晶模型
    real(WP) :: x0, y0  !< 圆心坐标
    real(WP), allocatable :: vcord(:, :)  !< (2, nvpoin), voronoi多边形的顶点坐标
    real(WP) :: xv1, yv1, xv2, yv2, p1x, p1y, p2x, p2y
    real(WP) :: gx, gy, x1, x2, tx1, theta
    real(WP) :: twopi, epsilon_
    integer, allocatable :: vlnods(:, :)  !< (nvnode, nvelem), voronoi多边形的顶点下标
    integer, allocatable :: nnode2(:)     !< (nvelem), 记录每个voronoi多边形有几个顶点
    integer :: radius  !< 圆的半径
    integer :: xlength !< 节点距离圆心的距离
    integer :: i, j, jpoin
    integer :: fileid !< 读写inp文件的通道ID
    integer :: nvpoin !< 文件中一共有多少个voronoi多边形顶点
    integer :: nvnode !< 单个voronoi多边形最多有几个顶点
    integer :: nvelem !< 文件中一共有几个voronoi多边形
    integer :: ielem, igrain, inode, jnode, knode !< 循环变量
    integer :: mnode !< 循环中某个voronoi多边形的顶点数量

    ! 中心晶粒的半径为14.0
    radius = 14

    if(iflag == 1) then
      ngrain = 2
      x0 = Nx / 2.0_WP; y0 = Ny / 2.0_WP
      allocate(etas(Nx, Ny, ngrain))
      allocate(glist(ngrain), source=1)
      
      do j = 1, Ny
        do i = 1, Nx
          etas(i, j, 1) = 1.0_WP
          etas(i, j, 2) = 0.0_WP
          xlength = nint(sqrt((i - x0)**2 + (j - y0)**2))
          if(xlength <= radius)then
            etas(i, j, 1) = 0.0_WP    
            etas(i, j, 2) = 1.0_WP    ! eta为1表示其被第2个晶粒占据
          end if
        end do
      end do
    else if(iflag == 2) then
      ! 初始化固定变量
      twopi = 8.0_WP * atan(1.0_WP)
      epsilon_ = 1.0e-4_WP
      ! 打开文件
      open(newunit=fileid, file='./voronoi_out/grain_25.inp', action='read')
      ! 读取关键参数
      read(fileid, *) nvpoin, nvnode, nvelem, ngrain
      !write(stdout, *)nvpoin, nvnode, nvelem, ngrain
      ! 分配内存并进行初始化
      allocate(etas(Nx, Ny, ngrain), source=0.0_WP)
      allocate(glist(ngrain), source=1)
      ! 分配内存, 不进行初始化
      allocate(vcord(2, nvpoin))
      allocate(nnode2(nvelem))
      allocate(vlnods(nvnode + 1, nvelem)) ! 第nvnode + 1个元素表示voronoi单元属于哪个grain
      ! 读取voronoi多边形的顶点坐标
      do i = 1, nvpoin
        read(fileid, *) jpoin, vcord(1, i), vcord(2, i)
      end do
      ! 读取voronoi多边形的顶点坐标对应下标
      do i = 1, nvelem
        read(fileid, *) jpoin, vlnods(:, i)
      end do
      ! 统计每个voronoi多边形节点个数
      do i = 1, nvelem
        nnode2(i) = count(vlnods(:, i) /= 0) - 1
      end do
      ! 遍历所有etas的所有格点, 确定其属于哪个grain
      do j = 1, Ny
        do i = 1, Nx
          gx = i*dx; gy = j*dy
          do ielem = 1, nvelem
            igrain = vlnods(nvnode + 1, ielem)  ! 获取第ielem个voronoi单元是哪个晶粒
            theta = 0.0_WP
            mnode = nnode2(ielem) ! 获取第ielem个voronoi单元有多少个顶点
            do inode = 1, mnode
              knode = vlnods(inode, ielem)
              xv1 = vcord(1, knode); yv1 = vcord(2, knode)  ! 获取第knode个顶点的x, y坐标
              jnode = vlnods(inode + 1, ielem)
              if (inode == mnode) jnode = vlnods(1, ielem)  ! 对边界进行修正
              xv2 = vcord(1, jnode); yv2 = vcord(2, jnode)  ! 获取第jnode个顶点的x, y坐标
              ! 计算节点(i*dx, j*dy) 与第knode和第jnode之间形成的向量的分量
              p1x = xv1 - gx; p1y = yv1 - gy  
              p2x = xv2 - gx; p2y = yv2 - gy
              ! 计算上述向量的绝对长度
              x1 = sqrt(p1x*p1x + p1y*p1y)
              x2 = sqrt(p2x*p2x + p2y*p2y)
              ! 避免奇异性, 进行数值校正
              if(x1*x2 <= epsilon_) then
                theta = twopi
                tx1 = (p1x*p2x + p1y*p2y) / epsilon_
              else
                tx1 = (p1x*p2x + p1y*p2y) / (x1*x2)
              end if
              ! 数值误差校正
              if (tx1 >= 1.0_WP) tx1 = 0.9999999_WP
              theta = theta + acos(tx1)
            end do
            if (abs(theta - twopi) < epsilon_) then
              etas(i, j, igrain) = 1.0_WP
            end if
          end do
        end do
      end do
      !write(stdout, *) nnode2(1:5)
      !write(stdout, *) vcord(:, 1)
      !write(stdout, *) vlnods(:, 1)
      ! 关闭文件通道
      close(fileid)
    else 
      write(stdout, *) "Wrong value of iflag, the value of iflag must be 1 or 2."
      stop
    end if

  end subroutine init_grain_micro

  pure function free_energ_fd_ca_v1(i, j, ngrain, etas, igrain) result(dfdeta)
    integer, intent(in) :: i, j
    integer, intent(in) :: ngrain
    integer, intent(in) :: igrain
    real(WP),intent(in) :: etas(:, :, :)
    real(WP) :: dfdeta
    real(WP) :: sum
    integer :: k

    sum = 0.0
    do k = 1, ngrain
      if (k /= igrain) then
        sum = sum + etas(i, j, k) * etas(i, j, k)
      end if
    end do

    dfdeta = A*(2*B*etas(i, j, igrain)*sum + etas(i, j, igrain)**3 - etas(i, j, igrain))
    ! dfdeta = -A*etas(i, j, igrain) + B*etas(i, j, igrain)**3 + &
    !           2*etas(i, j, igrain)*sum
  end function free_energ_fd_ca_v1
end module case_2_v1
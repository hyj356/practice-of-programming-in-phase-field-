module case_8_V1
    use, intrinsic:: iso_c_binding, only: c_double, c_int, c_double_complex
    use, intrinsic:: iso_fortran_env, only:real32, real64, output_unit
    use FFTW, only: fftw_plan_dft_2d, fftw_execute_dft, fftw_destroy_plan, &
                  FFTW_FORWARD, FFTW_BACKWARD, &
                  c_ptr, FFTW_ESTIMATE
    implicit none
    private
    public :: write_vtk_file, prepare_fft, fft_2d, ifft_2d, init_FeCuNiMn, &
                Fe_Cu_Mn_Ni_free_energy

    integer, parameter, public :: wp = c_double
    integer, parameter, public :: cwp = c_double_complex
    integer, parameter, public :: ci = c_int
    integer, parameter, public :: stdout = output_unit
    real(WP), parameter, public :: pi = 4.0_wp * atan(1.0_WP)
    real(WP), parameter :: noise = 0.001_wp !< 期望的各个节点上浓度噪音
    real(WP), parameter :: R = 8.314472_wp  !< 气体常数
    real(WP), parameter :: coef = 10.0_wp   !< 多项式系数


contains

    subroutine init_FeCuNiMn(Cu, Mn, Ni, orp, Cu0, Mn0, Ni0, Nx, Ny)
        !! 初始化Fe，Cu，Ni和Mn四种元素的初始浓度
        real(wp), intent(out), allocatable :: Cu(:, :), Mn(:, :), Ni(:, :), orp(:, :)
        real(wp), intent(in) :: Cu0, Mn0, Ni0
        integer, intent(in) :: Nx, Ny
        real(wp), dimension(4) :: rng
        integer :: i, j

        ! 分配内存
        if(.not. allocated(Cu)) allocate(Cu(Nx, Ny))
        if(.not. allocated(Mn)) allocate(Mn(Nx, Ny))
        if(.not. allocated(Ni)) allocate(Ni(Nx, Ny))
        if(.not. allocated(orp)) allocate(orp(Nx, Ny))
        ! 初始化随机数种子
        call random_seed()
        ! 开始赋值
        do j = 1, Ny
            do i = 1, Nx
                call random_number(rng)
                Cu(i, j) = Cu0 + (0.5_wp - rng(1)) * noise
                Mn(i, j) = Mn0 + (0.5_wp - rng(2)) * noise
                Ni(i, j) = Ni0 + (0.5_wp - rng(3)) * noise
                orp(i, j) = 0.001_wp + noise * (0.5_wp - rng(4))
            end do
        end do
    end subroutine init_FeCuNiMn

    pure subroutine Fe_Cu_Mn_Ni_free_energy(Nx, Ny, tempr, Cu, Mn, Ni, orp, dGdCu, dGdMn, dGdNi, dGdor)
        ! 输入参数
        integer, intent(in) :: Nx, Ny
        real(wp), intent(in) :: tempr
        real(wp), intent(in) :: Cu(:, :)      !< Cu元素浓度在网格节点上的值
        real(wp), intent(in) :: Mn(:, :)      !< Mn元素浓度在网格节点上的值
        real(wp), intent(in) :: Ni(:, :)      !< Ni元素浓度在网格节点上的值
        real(wp), intent(in) :: orp(:, :)     !< 非保守序参数的值
        ! 输出参数
        real(wp), intent(inout) :: dGdCu(:, :)!< G对Cu元素浓度的导数
        real(wp), intent(inout) :: dGdMn(:, :)!< G对Mn元素浓度的导数
        real(wp), intent(inout) :: dGdNi(:, :)!< G对Ni元素浓度的导数
        real(wp), intent(inout) :: dGdor(:, :)!< G对非保守序参数的导数
        ! 函数内部的局部变量
        real(WP) :: RT
        real(wp) :: constw, conste
        real(wp) :: eta2, eta3, eta4
        real(wp) :: c02, c03, c04
        real(wp) :: c1, c2, c3, c4
        real(wp) :: funch, funcg
        real(wp) :: elaste
        real(wp) :: dgCua, delasCu, dgCug
        real(wp) :: dgMna, delasMn, dgMng
        real(wp) :: dgNia, delasNi, dgNig
        real(wp) :: gCuNia,gCuNig
        integer :: i, j

        ! 初始化计算过程中的常量的值
        RT = R * tempr
        constw = 5.0e3_wp / RT
        conste = 2.14e11_wp * 7.09e-6_wp / RT
        eta2 = 3.29e-2_wp; eta3 = 5.22e-4_wp; eta4 = 4.75e-4_wp;
        c02 = 0.15_wp  ! Cu的平均浓度
        c03 = 0.01_wp  ! Mn的平均浓度
        c04 = 0.01_wp  ! Ni的平均浓度

        ! 开始计算每个节点上的导数
        do j = 1, Ny
            do i = 1, Nx
                ! 计算第(i, j)个节点上4种元素的浓度值
                c2 = Cu(i, j); c3 = Mn(i, j); c4 = Ni(i, j); c1 = 1.0_wp - c2 - c3 - c4
                ! 如果该网格节点上, Fe元素浓度在0.05% - 99.95%之间
                if(c1 > 0.0005_wp .and. c1 < 0.9995_wp) then
                    funch = (3.0_wp - 2.0_wp*orp(i,j))*orp(i,j)**2
                    funcg = (1.0_wp - orp(i,j))*orp(i,j)
                    elaste = conste*((c4-c04)*eta4+(c3- c03)*eta3+(c2-c02)*eta2)**2
                    delasCu = 2.0_wp*conste*eta2*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)
                    delasMn = 2.0_wp*conste*eta3*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)
                    delasNi = 2.0_wp*conste*eta4*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)
                    ! 计算alpha相的相关系数
                    dgCua = 1.4613878411949395E-4_wp* &
                            (-201.3642400_wp*c1*c4-(201.364240_wp*(-2*c4-c3-c2+1.0_wp)-2016.04498_wp) &
                            *c4+10672.046_wp*c4+30000.0_wp*c3*c1+36076.894_wp &
                            *c1+6842.810456_wp*(log(c2) -log(c1))+(6252.0_wp-9865.0_wp*(c2-c3))*c3 &
                            -39865.0_wp*c2*c3+1740.949_wp*c3-36076.894_wp*c2+2984.135_wp)
                    ! 此变量第三行对应原文缺失一个运算符号, 根据上下文猜测用-代替
                    dgMna = 1.4613878411949395E-4_wp*(-201.364240_wp &
                             *c1*c4+(6276.0_wp*(c3-c4)-48642.59_wp) *c4-(201.364240_wp*(-2*c4-c3-c2+1.0_wp) &
                             -2016.04498_wp)*c4+6276.0_wp*c3*c4+30000.0_wp*c2*c1-1740.949_wp*c1 &
                             +6842.810456_wp*(log(c3)-log(c1)) - 20135.0_wp*c2*c3+1740.949_wp*c3+&
                             c2*(6252.0_wp-9865.0_wp*(c2-c3))-36076.894_wp*c2-33919.89130169578_wp)
                    dgNia = 1.4613878411949395E-4_wp*(6842.810456_wp &
                             *(log(c4)-log(c1))-402.7284800_wp*c1*c4-(201.364240_wp* &
                             (-2*c4-c3-c2+1.0_wp)-2016.04498_wp)*c4-6276.0_wp*c3*c4  &
                             +(201.3642400_wp*(-2*c4-c3-c2+1.0_wp)-2016.04498_wp)*   &
                             c1+c3*(6276.0_wp*(c3-c4)-48642.59_wp)- 30000.0_wp*c2*c3 &
                             +1740.949_wp*c3-25404.848_wp*c2 +5788.49600_wp)
                    ! 计算gamma相的相关系数
                    dgCug = 1.4613878411949395E-4_wp*(c1*(5672.8150_wp*(c4+c3+2*c2-1.0_wp)+42968.802_wp) - & ! 这里少一个运算符, 暂时用+代替
                             c2*(5672.8150_wp*(c4+c3+2*c2-1.0_wp)+42968.802_wp) &
                             +(1451.610348_wp*(-2*c4-c3-c2+1.0_wp)- 7419.147789_wp)*c1*c4-47841.3_wp*c1*c4 &
                             +(10672.046_wp-2868.3240_wp*(c2-c4))*c4-(-725.805174_wp*(-2*c4-c3-c2+1.0_wp)**2 &
                             +7419.147789_wp*(-2*c4-c3-c2+1.0_wp)- 9359.746009_wp)*c4+44972.976_wp*c2*c4 &
                             -26591.0_wp*c3*c1+11345.63_wp*c2*c1-c3*(-259.0_wp*(-c4-2*c3-c2+1.0_wp)-4581.105_wp)&
                             + 6842.810456_wp*(log(c2)-log(c1))+(-1969.5_wp*(c2-c3)**3-8131.0_wp*(c2-c3)+9927.1_wp)*c3 &
                             +c2*(-5908.5_wp*(c2-c3)**2-8131.0_wp)*c3+26850.0_wp*c2* c3+566.3008361308123_wp)
                    dgMng = 1.4613878411949395E-4_wp*(-c2*(5672.8150_wp*(c4+c3+2*c2-1.0_wp)+42968.802_wp)+ &
                             (1451.610348_wp*(-2*c4-c3-c2+1.0_wp)-7419.147789999999_wp)*c1*c4+ &
                             (6276.0_wp*(c3-c4)-49205.406_wp)*c4-(-725.805174_wp*(-2*c4-c3-c2+1.0_wp)**2 &
                             +7419.147789999999_wp*(-2*c4-c3-c2+1.0_wp)-9359.746009999999_wp)*c4+6276.0_wp*c3*c4 &
                             +47841.3_wp*c2*c4+(-259.0_wp*(-c4-2*c3-c2+1.0_wp)-4581.105_wp)* &
                             c1+518.0_wp*c3*c1-21177.185_wp*c2*c1-c3*(-259.0_wp*(-c4-2*c3-c2+1.0_wp)-4581.105_wp) &
                             +6842.810456_wp*(log(c3)-log(c1))+c2*(5908.5_wp*(c2-c3)**2+8131.0)*c3+26850.0_wp*c2*c3 &
                             +c2*(-1969.5_wp*(c2-c3)**3-8131.0_wp*(c2-c3)+9927.1_wp)-33766.35316921635_wp)
                    dgNig = 1.4613878411949395E-4_wp*(6842.810456_wp*(log(c4)-log(c1)) &
                             -c2*(5672.815000000001_wp*(c4+c3+2*c2-1.0_wp)+42968.802_wp)+(2903.220696_wp* &
                             (-2*c4-c3-c2+1.0_wp)-14838.29558_wp)*c1*c4-(-725.805174_wp*(-2*c4-c3-c2+1.0_wp)**2 &
                             +7419.147789999999_wp*(-2*c4-c3-c2+1.0_wp)-9359.746009999999_wp)*c4-6276.0_wp*c3*c4 &
                             +50709.624_wp*c2*c4+(-725.805174_wp*(-2*c4-c3-c2+1.0_wp)*2+7419.147789999999_wp* &
                             (-2*c4-c3-c2+1.0_wp)-9359.746009999999_wp)*c1+259.0_wp*c3*c1 &
                             -42168.485_wp*c2*c1+c3*(6276.0_wp*(c3-c4)-49205.406_wp)-c3*(-259.0_wp* &
                             (-c4-2*c3-c2+1.0_wp)-4581.105_wp)+c2*(10672.046_wp &
                             -2868.324000000001_wp*(c2-c4))+26850.0_wp*c2*c3+566.3008361308123_wp)
                    ! 计算alpha和gamma的自由能
                    gCuNia = 1.4613878411949395E-4_wp*(6842.810456_wp*(c4 &
                             *log(c4)+log(c1)*c1+c3*log(c3)+c2*log(c2))+(201.36424_wp*&
                             (-2*c4-c3-c2+1.0_wp)-2016.04498_wp)*c1*c4+c3*(6276.0_wp*(c3-c4)&
                             -48642.59_wp)*c4+10672.046_wp*c2*c4+5788.496_wp*c4+30000.0_wp*c2*c3*c1 &
                             -1740.949_wp*c3*c1+36076.894_wp*c2*c1+c2*(6252.0_wp &
                             -9865.0_wp*(c2-c3))*c3-33919.89130169578_wp*c3+2984.135_wp*c2)
                    gCuNig = 1.4613878411949395E-4_wp*(6842.810456_wp*(c4*log(c4)+log(c1)*c1 &
                             +c3*log(c3)+c2*log(c2))+c2*c1*(5672.815_wp*(c4+c3+2*c2-1.0_wp) &
                             +42968.802_wp)+(-725.805174_wp*(-2*c4-c3-c2+1.0_wp)**2 &
                             +7419.147789999999_wp*(-2*c4-c3-c2+1.0_wp) &
                             -9359.746009999999_wp)*c1*c4-47841.3_wp*c2*c1*c4+c3*(6276.0_wp*(c3-c4) &
                             -49205.406_wp)*c4+c2*(10672.046-2868.324000000001_wp &
                             *(c2-c4))*c4+c3*(-259.0_wp*(-c4-2*c3-c2+1.0_wp)-4581.105_wp)*c1-26850.0_wp*c2*c3*c1 &
                             -566.3008361308123_wp*c1+c2*(-1969.5_wp*(c2-c3)**3 &
                             -8131.0_wp*(c2-c3)+9927.1_wp)*c3-34332.65400534716_wp*c3)
                    ! 计算四种元素的自由能导数
                    dGdCu(i, j) = (1.0_wp - funch) * (dgCua + delasCu) + funcg * dgCug  !Cu
                    dGdMn(i, j) = (1.0_wp - funch) * (dgMna + delasMn) + funcg * dgMng  !Mn
                    dGdNi(i, j) = (1.0_wp - funch) * (dgNia + delasNi) + funcg * dgNig  !Ni
                    ! 计算自由能导数对序参数orp的导数
                    dgdor(i,j) = (gCuNia+elaste)*(2.0_wp*orp(i,j)**2-2.0_wp*(3.0_wp &
                                 -2.0_wp*orp(i,j))*orp(i,j))-2.0_wp*constw*(1.0_wp &
                                 -orp(i,j))*orp(i,j)**2-2.0_wp*gCuNig*orp(i,j)**2 &
                                 +2.0_wp*constw*(1.0_wp-orp(i,j))**2*orp(i,j)+2.0_wp*gCuNig &
                                 *(3.0_wp-2*orp(i,j))*orp(i,j)
                end if

                ! 如果Fe的浓度大于99.95%, 自由能的导数的计算就是另一套规则
                if (c1 > 0.9995_wp) then
                    funch = (3.0_wp-2.0_wp*orp(i,j))*orp(i,j)**2
                    funcg = (1.0_wp-orp(i,j))*orp(i,j)
                    elaste = conste*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)**2
                    ! 计算弹性常数相关导数
                    delasCu = 2.0_wp*conste*eta2*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)
                    delasMn = 2.0_wp*conste*eta3*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)
                    delasNi = 2.0_wp*conste*eta4*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)
                    ! 计算alpha和gamma的自由能
                    gCuNia = coef*((1.0_wp-c2-c3-c4)-0.9995_wp)**2
                    gCuNig = coef*((1.0_wp-c2-c3-c4)-0.9995_wp)**2
                    ! 计算Cu元素自由能对alpha和gamma相的导数
                    dgCua = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.9995_wp)
                    dgCug = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.9995_wp)
                    ! 计算Mn元素自由能对alpha和gamma相的导数
                    dgMna = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.9995_wp)
                    dgMng = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.9995_wp)
                    ! 计算Ni元素自由能对alpha和gamma相的导数
                    dgNia = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.9995_wp)
                    dgNig = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.9995_wp)
                    ! 计算3种元素自由能对各自元素浓度的导数
                    dgdCu(i, j) = (1.0_wp-funch)*(dgCua+delasCu)+funch*dgCug     ! Cu元素
                    dgdMn(i, j) = (1.0_wp-funch)*(dgMna+delasMn)+funch*dgMng     ! Mn元素
                    dgdNi(i, j) = (1.0_wp-funch)*(dgNia+delasNi)+funch*dgNig     ! Ni元素
                    ! 计算自由能导数对序参数orp的导数
                    dgdor(i, j) = (gCuNia+elaste)*(6.0_wp*orp(i,j)**2-6.0_wp*orp(i,j))+constw*(2.0_wp &
                                 *orp(i,j)-6.0_wp*orp(i,j)**2+4.0_wp*orp(i,j)**3)+gcunig &
                                 *(6.0_wp*orp(i,j)-6.0_wp*orp(i,j)**2)
                end if
                 ! 如果Fe的浓度小于0.05%, 自由能的导数的计算就是另一套规则
                if (c1 < 0.0005_wp) then
                    funch = (3.0_wp-2.0_wp*orp(i,j))*orp(i,j)**2
                    funcg = (1.0_wp-orp(i,j))*orp(i,j)
                    elaste = conste*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)**2
                    ! 计算弹性常数相关导数
                    delasCu = 2.0_wp*conste*eta2*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)    ! Cu元素
                    delasMn = 2.0_wp*conste*eta3*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)    ! Mn元素
                    delasNi = 2.0_wp*conste*eta4*((c4-c04)*eta4+(c3-c03)*eta3+(c2-c02)*eta2)    ! Ni元素
                    ! 计算alpha和gamma的自由能
                    gCuNia = coef*((1.0_wp-c2-c3-c4)-0.0005_wp)**2
                    gCuNig = coef*((1.0_wp-c2-c3-c4)-0.0005_wp)**2
                    ! 计算Cu元素自由能对alpha和gamma相的导数
                    dgCua = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.0005_wp)
                    dgCug = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.0005_wp)
                    ! 计算Mn元素自由能对alpha和gamma相的导数
                    dgMna = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.0005_wp)
                    dgMng = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.0005_wp)
                    ! 计算Ni元素自由能对alpha和gamma相的导数
                    dgNia = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.0005_wp)
                    dgNig = -2.0_wp*coef*((1.0_wp-c2-c3-c4)-0.0005_wp)
                    ! 计算3种元素自由能对各自元素浓度的导数
                    dgdCu(i, j) = (1.0_wp-funch)*(dgCua+delasCu)+funch*dgCug     ! Cu元素
                    dgdMn(i, j) = (1.0_wp-funch)*(dgMna+delasMn)+funch*dgMng     ! Mn元素
                    dgdNi(i, j) = (1.0_wp-funch)*(dgNia+delasNi)+funch*dgNig     ! Ni元素
                    ! 计算自由能导数对序参数orp的导数
                    dgdor(i, j) = (gCuNia+elaste)*(6.0_wp*orp(i,j)**2-6.0_wp*orp(i,j))+constw*(2.0_wp &
                                 *orp(i,j)-6.0_wp*orp(i,j)**2+4.0_wp*orp(i,j)**3)+gCuNig &
                                 *(6.0_wp*orp(i,j)-6.0_wp*orp(i,j)**2)
                end if
            end do
        end do

    end subroutine Fe_Cu_Mn_Ni_free_energy

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
end module

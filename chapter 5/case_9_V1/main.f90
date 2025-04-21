program hello
    use case_9_V1, only: &
        write_vtk_file, prepare_fft, fft_2d, ifft_2d, &
        green_tensor_2d, solve_elasticity_v1, &
        micro_ch_pre, free_energ_ch, &
        wp, cwp, stdout
    implicit none
    ! ģ�������С�������
    integer, parameter :: Nx = 256
    integer, parameter :: Ny = 256
    real(wp), parameter :: iNxNy = 1.0_wp / (Nx*Ny)
    real(wp), parameter :: dx = 1.0_wp
    real(wp), parameter :: dy = 1.0_wp
    ! ʱ����ֲ���
    integer, parameter :: nstep = 5000
    integer, parameter :: nprint = 25
    real(wp), parameter :: dt = 5.0e-2_wp
    real(wp), parameter :: coefA = 1.0_wp
    ! ������س���
    real(wp), parameter :: c0 = 0.4_wp
    real(wp), parameter :: mobility = 1.0_wp
    real(wp), parameter :: grad_coef = 0.5_wp
    ! ���ϵĵ��Գ���
    real(wp), parameter :: cm11 = 1400.0_wp
    real(wp), parameter :: cm12 = 600.0_wp
    real(wp), parameter :: cm44 = 400.0_wp
    !--------------------------------------
    real(wp), parameter :: cp11 = 2.0_wp*cm11
    real(wp), parameter :: cp12 = 2.0_wp*cm12
    real(wp), parameter :: cp44 = 2.0_wp*cm44
    ! ����Ӧ��
    real(wp), parameter :: ei0 = 0.01_wp
    ! ʩ�ӵ�Ӧ��
    real(wp), dimension(3), parameter :: ea = [0.0_wp, 0.01_wp, 0.0_wp]
    ! �����õ�������
    real(wp), allocatable :: s11(:, :), s12(:, :), s22(:, :)
    real(wp), allocatable :: e11(:, :), e12(:, :), e22(:, :)
    real(wp), allocatable :: ed11(:, :), ed12(:, :), ed22(:, :)
    real(wp), allocatable :: con(:, :), dfdcon(:, :), tmatx(:, :, :, :, :, :)
    real(wp), allocatable :: kx(:), ky(:), k2(:, :), k4(:, :)
    real(wp), allocatable :: delsdc(:, :)
    complex(cwp), allocatable :: conk(:, :), dfdconk(:, :), delsdck(:, :)
    ! �����õ��ľֲ�������ѭ������
    integer :: istep, i, j
    real(wp) :: denom
    complex(cwp) :: numer

    ! �����ڴ�
    allocate(s11(Nx, Ny), s12(Nx, Ny), s22(Nx, Ny), source=0.0_wp)
    allocate(e11(Nx, Ny), e12(Nx, Ny), e22(Nx, Ny), source=0.0_wp)
    allocate(ed11(Nx, Ny), ed12(Nx, Ny), ed22(Nx, Ny), source=0.0_wp)
    allocate(dfdcon(Nx, Ny), delsdc(Nx, Ny))
    allocate(conk(Nx, Ny), dfdconk(Nx, Ny), delsdck(Nx, Ny))
    ! ��ʼ��Ũ������
    call micro_ch_pre(Nx, Ny, c0, con)

    ! ���㸵��Ҷ�任��ϵ������
    call prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)

    ! ����green tensor
    call green_tensor_2d(tmatx, Nx, Ny, kx, ky, cm11, cm12, cm44, cp11, cp12, cp44)

    ! ��ʼ�ݻ�
    do istep = 1, nstep
        ! ���������ܶ�Ũ��c�ĵ���
        do j = 1, Ny
            do i = 1, Nx
                dfdcon(i, j) = free_energ_ch(con, i, j)
            end do
        end do
        ! ���㵯�����ܵĵ���
        call solve_elasticity_v1(Nx, Ny, cm11, cm12, cm44, cp11, cp12, cp44, ed11, ed12, ed22, &
                                 ei0, ea, con, s11, s12, s22, e11, e12, e22, delsdc, tmatx)
        ! ���и���Ҷ�任
        call fft_2d(input=con, output=conk, Nx=Nx, Ny=Ny)
        call fft_2d(input=dfdcon, output=dfdconk, Nx=Nx, Ny=Ny)
        call fft_2d(input=delsdc, output=delsdck, Nx=Nx, Ny=Ny)
        ! ����ʱ�����
        do j = 1, Ny
            do i = 1, Nx
                numer = dt*mobility*k2(i, j)*(dfdconk(i, j) + delsdck(i, j))
                denom = 1.0_wp + dt*coefA*mobility*grad_coef*k4(i, j)
                conk(i, j) = (conk(i, j) - numer) / denom
            end do
        end do
        ! ͨ���渵��Ҷ�任��conk��Ϊcon
        call ifft_2d(input=conk, output=con, Nx=Nx, Ny=Ny, iNxNy=iNxNy)
        ! ��ֵ���У��
        where(con > 0.9999_wp)
            con = 0.9999_wp
        else where(con < 0.00001_wp)
            con = 0.00001_wp
        end where
        ! ��������Ϣ
        if(mod(istep, nprint) == 0) then
            write(stdout, *) "Done step: ", istep
            call write_vtk_file("./VTK/dump_", Nx, Ny, dx, dy, con, istep)
        end if

    end do

    print *, "-------------------------------All done!-------------------------------"

end program


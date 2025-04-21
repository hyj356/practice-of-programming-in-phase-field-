
program hello
    use case_7_V1, only: wp, cwp, stdout, pi, coefA, &
        write_vtk_file, prepare_fft, fft_2d, ifft_2d, &
        init_grain_micro, free_energ_fd_ca_v1
    implicit none
    integer, parameter :: Nx = 64
    integer, parameter :: Ny = 64
    real(wp), parameter :: iNxNy = 1.0_wp / (Nx*Ny)
    real(wp), parameter :: dx = 0.5_wp
    real(wp), parameter :: dy = 0.5_wp
    real(wp), parameter :: dt = 0.005
    integer, parameter :: nstep = 5000
    integer, parameter :: nprint = 50
    integer, parameter :: iflag = 2     !< iflagΪ1ʱ, ģ�ͺ���2������, iflag=2, ģ�ͺ���25������
    ! ������صĲ���
    real(wp), parameter :: mobil = 5.0_wp
    real(wp), parameter :: grcoef = 0.1_wp
    ! �������������
    real(wp), allocatable, target :: etas(:, :, :)
    complex(cwp), allocatable :: etak(:, :)
    real(wp), allocatable :: eta2(:, :)
    real(wp), pointer :: eta(:, :)
    real(wp), allocatable :: dfdeta(:, :)
    complex(cwp), allocatable :: dfdetak(:, :)
    real(wp), allocatable :: kx(:), ky(:)
    real(wp), allocatable :: k2(:, :), k4(:, :)
    integer, allocatable :: glist(:)
    integer :: ngrain
    integer :: i, j, istep, igrain
    real(wp) :: grain_sum
    complex(cwp) :: numer       !< ŷ������ʽʱ����ֹ�ʽ�еķ���
    real(wp) :: denom    !< ŷ������ʽʱ����ֹ�ʽ�еķ�ĸ

    ! �����ڴ�
    allocate(dfdeta(Nx, Ny), dfdetak(Nx, Ny), eta2(Nx, Ny), etak(Nx, Ny))
    ! ��ʼ�����������������
    call init_grain_micro(etas, ngrain, glist, Nx, Ny, dx, dy, iflag)
    ! ׼����ά����Ҷ�任�Ķ�Ӧϵ��
    call prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)
    ! ��ʼ����
    do istep = 1, nstep
        ! �������о���
        do igrain = 1, ngrain
            ! �����igrain��������û������
            if(glist(igrain) == 1) then
                eta => etas(:, :, igrain)
                ! ����ÿ���ڵ�������ܵ���
                do j = 1, Ny
                    do i = 1, Nx
                        dfdeta(i, j) = free_energ_fd_ca_v1(i, j, ngrain, etas, igrain)
                    end do
                end do
                ! �Ե�igrain�������������eta���и���Ҷ�任
                call fft_2d(input=eta, output=etak, Nx=Nx, Ny=Ny)
                ! �������ܵ������и���Ҷ�任
                call fft_2d(input=dfdeta, output=dfdetak, Nx=Nx, Ny=Ny)
                ! ��һ��ÿ���ڵ����ϵ��, ����ʱ�����
                do j = 1, Ny
                    do i = 1, Nx
                        numer = dt*mobil*dfdetak(i, j)
                        denom = 1.0_wp + dt*coefA*mobil*grcoef*k2(i, j)
                        etak(i, j) = (etak(i, j) - numer) / denom
                    end do
                end do
                ! ��etakת��ʵ������
                call ifft_2d(input=etak, output=eta, Nx=Nx, Ny=Ny, iNxNy=iNxNy)
                ! ��ֵ���У��
                where (eta > 0.99999_WP)
                    eta = 0.99999_WP
                else where (eta < 0.00001_WP)
                    eta = 0.00001_WP
                end where
                ! �����igrain���������ж��پ����������
                grain_sum = sum(eta) * iNxNy
                if (grain_sum < 0.001_wp) glist(igrain) = 0
            end if
        end do
        if(mod(istep, nprint) == 0) then
            write(stdout, *) "Done step: ", istep
            eta2(:, :) = sum(etas(:, :, :)*etas(:, :, :), dim=3)
            ! д��vtk�ļ�
            call write_vtk_file("./VTK_polycrystal/dump_", Nx, Ny, dx, dy, eta2, istep)
        end if
    end do

end program


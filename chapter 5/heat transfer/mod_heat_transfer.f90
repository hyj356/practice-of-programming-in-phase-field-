module heat_transfer
  use iso_c_binding, only: c_double, c_int
  use iso_fortran_env, only:real32, real64, output_unit
  use FFTW, only:  &
                      ! FFTW���еĺ���
                      fftw_plan_dft_1d, fftw_execute_dft, fftw_destroy_plan, &
                      fftw_cleanup, &
                      ! ��C���Խ����ĳ���
                      cwp => c_double_complex, c_ptr, &
                      ! FFTW���еĳ���
                      FFTW_ESTIMATE, FFTW_FORWARD, FFTW_BACKWARD
  implicit none
  private
  public :: initial_temp, write_temp, fft_1d, ifft_1d
  integer, parameter, public :: wp = c_double
  integer, parameter, public :: ci = c_int
  integer, parameter, public :: stdout = output_unit
  real(WP), parameter, public :: pi = 4.0_wp * atan(1.0_WP)

contains
  pure subroutine initial_temp(temp, x, kx, k2, Nx, dx)
    ! ��ʼ����Ӧ�ĸ���Ҷϵ��
    real(WP), intent(out), allocatable :: temp(:), x(:), kx(:), k2(:)
    integer, intent(in) :: Nx
    real(WP), intent(in) ::dx
    real(WP) :: fk1, delkx
    integer :: i, Nx21, Nx2

    allocate(temp(Nx), x(Nx), kx(Nx+1), k2(Nx))
    ! ��ʼ�������ڵ���¶�
    do i = 1, Nx
      x(i) = i * dx
      temp(i) = 0.0_WP
      if(i >= 44 .and. i <= 84) temp(i) = 1.0_WP
    end do
    ! �����Ӧ�ڵ�ĸ���Ҷϵ��
    Nx21 = Nx / 2 + 1; Nx2 = 2 + Nx; delkx = (2.0*pi)/(Nx*dx);
    do i = 1, Nx21
      fk1 = (i-1) * delkx
      kx(i) = fk1
      kx(Nx2 - i) = -fk1
    end do

    do i = 1, Nx
      k2(i) = kx(i) * kx(i)
    end do

  end subroutine

  subroutine write_temp(filename, temp)
    !! ���¶�д�����ļ���
    character(len=*), intent(in) :: filename
    real(WP), intent(in) :: temp(:)
    integer :: i, N, fileid

    open(newunit=fileid, file=filename, action='write')
    N = size(temp)
    do i = 1, N
      write(fileid, *) temp(i)
    end do
    close(fileid)
  end subroutine

  subroutine fft_1d(input, output, Nx)
    ! ����һά�ĸ���Ҷ���ٱ任
    real(WP), intent(inout) :: input(:)          ! �����ʵ������
    complex(cwp), intent(inout) :: output(:)  ! ����ĸ�������
    complex(cwp), allocatable :: c_input(:)
    integer(c_int), intent(in) :: Nx
    type(c_ptr) :: plan

    ! �����ڴ沢����ת��Ϊ��������
    allocate(c_input(Nx))
    c_input = cmplx(input, kind=cwp)
    ! ������Ҷ���ٱ任, �ó��������ж�ʹ�ú����㷨����FFT
    plan = fftw_plan_dft_1d(Nx, c_input, output, sign=FFTW_FORWARD, flags=FFTW_ESTIMATE)
    ! ����FFTW����п��ٸ���Ҷ�任
    call fftw_execute_dft(plan, c_input, output)
    ! �������ƻ����ڴ�
    call fftw_destroy_plan(plan)
    call fftw_cleanup()
    deallocate(c_input)

  end subroutine

  subroutine ifft_1d(input, output, Nx)
    ! ����һά���渵��Ҷ���ٱ任
    complex(cwp), intent(inout) :: input(:)          ! ����ĸ�������
    real(wp), intent(inout) :: output(:)  ! �����ʵ������
    integer(c_int), intent(in) :: Nx
    type(c_ptr) :: plan
    complex(cwp), allocatable :: c_output(:)  ! ����ĸ�������

    allocate(c_output(Nx))
    ! ������Ҷ���ٱ任, �ó��������ж�ʹ�ú����㷨����FFT
    plan = fftw_plan_dft_1d(Nx, input, c_output, sign=FFTW_BACKWARD, flags=FFTW_ESTIMATE)
    ! ����FFTW����п��ٸ���Ҷ�任
    call fftw_execute_dft(plan, input, c_output)
    ! ������תΪʵ��
    output = real(c_output, kind=wp) / real(Nx, kind=wp)
    ! �������ƻ����ڴ�
    call fftw_destroy_plan(plan)
    call fftw_cleanup()
    deallocate(c_output)
  end subroutine
end module

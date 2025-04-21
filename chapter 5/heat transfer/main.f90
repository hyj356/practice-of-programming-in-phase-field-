program hello
    use iso_c_binding, only: cwp => c_double_complex, c_ptr
    use heat_transfer, only: initial_temp, fft_1d, ifft_1d, write_temp,  &
                              wp, stdout, ci
    implicit none
    integer, parameter :: Nx = 128
    real(WP), parameter :: dx = 1.0_WP
    real(WP), parameter :: dt = 0.2_WP
    integer, parameter :: nstep = 600
    integer, parameter :: nprint = 200
    real(WP), dimension(:), allocatable :: temp, x, kx, k2
    complex(cwp), allocatable :: temp_fft(:)
    integer :: istep, i
    type(c_ptr) ::plan

    ! 分配内存
    allocate(temp(Nx), temp_fft(Nx))
!    temp = [0.3_wp, 2.0_wp, 2.0_wp, 7.0_wp, 99.0_wp]
!    call fft_1d(input=temp, output=temp_fft, Nx=Nx)
!    do i = 1, Nx
!      write(stdout, *) temp(i), temp_fft(i)
!    end do
!    write(stdout, *) '----------------------------All done!----------------------------'
!    call ifft_1d(input=temp_fft, output=temp, Nx=Nx)
!    do i = 1, Nx
!      write(stdout, *) temp_fft(i), temp(i)
!    end do
    ! 初始化求解数组
    call initial_temp(temp, x, kx, k2, Nx, dx)
    call write_temp("./initial_temp.txt", temp)
    ! 将温度场转换到傅里叶空间中
    call fft_1d(input=temp, output=temp_fft, Nx=Nx)
    ! 开始迭代求解
    do istep = 1, nstep
      ! 进行时间积分
      do i = 1, Nx
        temp_fft(i) = temp_fft(i) - dt * k2(i) * temp_fft(i)
      end do
    end do
    ! 计算完成之后, 通过傅里叶逆变换求得温度场
    call ifft_1d(input=temp_fft, output=temp, Nx=Nx)
    ! 将结果存储到txt文件中
    call write_temp("./final_temp.txt", temp)

    write(stdout, *) '----------------------------All done!----------------------------'
end program


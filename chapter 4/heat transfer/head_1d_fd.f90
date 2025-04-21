program main
  use chapter4, only: WP, stdout, mu, dt, dx, Nx, initial_temp, write_temperature
  implicit none
  real(wp), allocatable :: temp(:)
  integer :: i
  real(WP) :: idx2  !< 其值等于1/(dx*dx)

  idx2 = 1.0_WP / (dx * dx)
  ! 初始化温度数组
  temp = initial_temp(Nx, 64-20, 64+20, 1.0_WP)
  ! 将初始化之后的数组写出
  call write_temperature("initial_temp.txt", temp)
  ! 开始循环仿真
  do i = 1, 200
    temp(2:Nx-1) = temp(2:Nx-1) + mu * dt * idx2 * &
                    (temp(1:Nx-2) - 2.0_WP * temp(2:Nx-1) + temp(3:Nx))
  end do

  ! 将迭代之后的数组写出
  call write_temperature("timestep_200.txt", temp)

  do i = 1, 200
    temp(2:Nx-1) = temp(2:Nx-1) + mu * dt * idx2 * &
                    (temp(1:Nx-2) - 2.0_WP * temp(2:Nx-1) + temp(3:Nx))
  end do

  ! 将迭代之后的数组写出
  call write_temperature("timestep_400.txt", temp)

  do i = 1, 200
    temp(2:Nx-1) = temp(2:Nx-1) + mu * dt * idx2 * &
                    (temp(1:Nx-2) - 2.0_WP * temp(2:Nx-1) + temp(3:Nx))
  end do

  ! 将迭代之后的数组写出
  call write_temperature("timestep_600.txt", temp)
  
end program main
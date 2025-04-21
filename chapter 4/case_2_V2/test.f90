program main
  !gfortran .\mod_case2_v2.f90 .\test.f90 -o test -O2; ./test.exe
  use case_2_v2, only: eyes, toepliz, WP, stdout, kron
  implicit none
  real(WP), allocatable :: E(:, :), T(:, :), temp(:, :)
  real(WP) :: A(2, 2), B(2, 2)
  integer :: i

  E = eyes(4)
  T = toepliz([1.0_WP, 2.0_WP, 3.0_WP])
  A = reshape([1, 2, 3, 4], [2, 2])
  B = reshape([0, 5, 6, 7], [2, 2])
  temp = kron(A, B)
  write(stdout, *) "eyes(4) = "
  do i = 1, 4
    write(stdout, *) E(i, :)
  end do
  write(stdout, *) "toepliz([1, 2, 3]) = "
  do i = 1, 3
    write(stdout, *) T(i, :)
  end do
  write(stdout, *) "A = [1 2; 3 4], B = [0 5; 6 7]; kron(A, B) = "

  do i = 1, 4
    write(stdout, *) temp(i, :)
  end do
  
end program main

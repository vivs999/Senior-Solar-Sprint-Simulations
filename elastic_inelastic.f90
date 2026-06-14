program momentum_carts
    implicit none

    ! Declare variables
    real(8) :: m1, m2
    real(8) :: v1_i, v2_i, v1_f, v2_f
    real(8) :: p1_i, p2_i, p1_f, p2_f
    real(8) :: p_sys_i, p_sys_f, ratio

    ! Input section
    print *, "Enter mass of cart 1 (kg):"
    read(*,*) m1

    print *, "Enter initial velocity of cart 1 (m/s):"
    read(*,*) v1_i

    print *, "Enter final velocity of cart 1 (m/s):"
    read(*,*) v1_f

    print *, "Enter mass of cart 2 (kg):"
    read(*,*) m2

    print *, "Enter initial velocity of cart 2 (m/s):"
    read(*,*) v2_i

    print *, "Enter final velocity of cart 2 (m/s):"
    read(*,*) v2_f

    ! Momentum calculations
    p1_i = m1 * v1_i
    p2_i = m2 * v2_i
    p_sys_i = p1_i + p2_i

    p1_f = m1 * v1_f
    p2_f = m2 * v2_f
    p_sys_f = p1_f + p2_f

    ! Output results
    print *
    write(*,'(A,1X,F10.3,1X,A)') "Initial momentum of cart 1:", p1_i, "kg·m/s"
    write(*,'(A,1X,F10.3,1X,A)') "Initial momentum of cart 2:", p2_i, "kg·m/s"
    write(*,'(A,1X,F10.3,1X,A)') "Initial system momentum:", p_sys_i, "kg·m/s"

    print *
    write(*,'(A,1X,F10.3,1X,A)') "Final momentum of cart 1:", p1_f, "kg·m/s"
    write(*,'(A,1X,F10.3,1X,A)') "Final momentum of cart 2:", p2_f, "kg·m/s"
    write(*,'(A,1X,F10.3,1X,A)') "Final system momentum:", p_sys_f, "kg·m/s"

    print *
    if (p_sys_i /= 0.0d0) then
        ratio = p_sys_f / p_sys_i
        write(*,'(A,1X,F10.3)') "Ratio (p_sys_final / p_sys_initial):", ratio
    else
        print *, "Ratio undefined (initial system momentum is zero)."
    end if

end program momentum_carts
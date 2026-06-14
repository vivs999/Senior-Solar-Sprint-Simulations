PROGRAM SPEED_RATIO_MONTECARLO
    USE omp_lib
    IMPLICIT NONE

    INTEGER, PARAMETER :: num_threads = 100
    INTEGER, PARAMETER :: num_iterations = 100000
    REAL, PARAMETER :: MAX_SIMULATION_TIME = 1000.0
    REAL, PARAMETER :: MIN_SPEED_THRESHOLD = 0.001
    INTEGER, PARAMETER :: STAGNATION_CHECKS = 500

    INTEGER :: i, bearing_count, zero_accel_count, failed_tests, valid_tests, stagnation_counter
    REAL :: g, wheel_radius, wheel_diameter, mass, gear_ratio, gear_efficiency
    REAL :: distance_remaining, time_step, total_time_sum, total_time_avg
    REAL :: final_speed_sum, final_speed_avg
    REAL :: motor_power_nominal, power_margin, stall_torque, free_speed_rpm
    REAL :: voltage_input, current_input, motor_efficiency, last_position

    REAL :: motor_rpm_max, motor_torque_stall
    REAL :: current_rpm, current_torque, current_speed, position, total_time
    REAL :: axle_normal_force, F_bearing_total, F_rolling, F_net, acceleration
    REAL :: F_drive, wheel_circumference, mu_bearing, mu_rolling
    REAL :: power_variation, angular_velocity, current_power, motor_power_max
    REAL :: F_air, rho, Cd, A
    REAL :: contact_chance, F_touch
    REAL, PARAMETER :: mu_wire = 0.05   ! PTFE teflon tape on nylon fishing line
    REAL, PARAMETER :: N_wire = 0.17    ! ~5% of car weight, guide geometry normal force
    LOGICAL :: is_in_contact

    REAL :: max_speed, min_speed, max_time, min_time

    ! Pitsco Motor 280
    stall_torque = 0.0127       ! 12.7 mN·m stall torque
    free_speed_rpm = 5200.0     ! no-load speed at 3V
    current_input = 1.55        ! current at max efficiency at 3V
    motor_efficiency = 0.60     ! RE-280 mechanical efficiency at max efficiency point
    voltage_input = 3.0         ! 2x AA batteries
    motor_power_nominal = voltage_input * current_input * motor_efficiency
    power_margin = 0.15 * motor_power_nominal

    ! Physical constants
    g = 9.81
    wheel_diameter = 0.0392
    wheel_radius = wheel_diameter / 2.0
    bearing_count = 4
    distance_remaining = 20.0
    time_step = 0.001
    mass = 0.297
    gear_ratio = 1.5
    gear_efficiency = 0.98

    ! Aerodynamics: rectangular box bluff body, 25mm x 100mm
    rho = 1.225
    Cd = 1.05
    A = 0.0025         ! box frontal area: 25mm x 100mm = 2500 mm²

    failed_tests = 0
    valid_tests = 0
    total_time_sum = 0.0
    final_speed_sum = 0.0
    max_speed = -1.0E30
    min_speed = 1.0E30
    max_time = -1.0E30
    min_time = 1.0E30

    CALL OMP_SET_NUM_THREADS(num_threads)

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& PRIVATE(i, mu_bearing, mu_rolling, motor_rpm_max, motor_torque_stall, &
    !$OMP& current_rpm, current_speed, position, total_time, wheel_circumference, &
    !$OMP& current_torque, F_drive, axle_normal_force, F_bearing_total, F_rolling, &
    !$OMP& F_net, acceleration, zero_accel_count, power_variation, current_power, F_touch, &
    !$OMP& angular_velocity, motor_power_max, last_position, stagnation_counter, F_air, &
    !$OMP& contact_chance, is_in_contact, current_input, motor_efficiency) &
    !$OMP& SHARED(motor_power_nominal, power_margin, g, wheel_diameter, wheel_radius, &
    !$OMP& bearing_count, distance_remaining, time_step, mass, gear_ratio, rho, Cd, A, &
    !$OMP& gear_efficiency, max_speed, min_speed, max_time, min_time) &
    !$OMP& REDUCTION(+:total_time_sum, final_speed_sum, valid_tests, failed_tests)
    DO i = 1, num_iterations
        CALL INIT_RANDOM_SEED()

        ! ZrO2 full ceramic ABEC3 bearings: 0.0005-0.002
        CALL RANDOM_NUMBER(mu_bearing)
        mu_bearing = 0.0005 + mu_bearing * 0.0015

        ! Thin latex band on ping pong ball, low deformation: 0.010-0.020
        CALL RANDOM_NUMBER(mu_rolling)
        mu_rolling = 0.004 + mu_rolling * 0.004  ! thin 0.5in latex on paper: 0.004-0.008

        ! Pitsco 280 no-load speed at 3V: 5000-5400 RPM
        CALL RANDOM_NUMBER(motor_rpm_max)
        motor_rpm_max = 5000.0 + motor_rpm_max * 400.0

        CALL RANDOM_NUMBER(motor_torque_stall)
        motor_torque_stall = 0.012 + motor_torque_stall * 0.002

        CALL RANDOM_NUMBER(power_variation)
        motor_power_max = motor_power_nominal - power_margin + 2.0 * power_margin * power_variation

        current_rpm = 0.0
        current_speed = 0.0
        position = 0.0
        total_time = 0.0
        zero_accel_count = 0
        last_position = 0.0
        stagnation_counter = 0
        wheel_circumference = 3.14159 * wheel_diameter

        simulation_loop: DO WHILE (position < distance_remaining)
            current_torque = motor_torque_stall * (1.0 - current_rpm / motor_rpm_max)

            angular_velocity = (2.0 * 3.14159 * current_rpm) / 60.0

            IF (angular_velocity > 0.0) THEN
                current_power = current_torque * angular_velocity
                IF (current_power > motor_power_max) THEN
                    current_torque = motor_power_max / angular_velocity
                END IF
            END IF

            current_torque = current_torque * gear_ratio * gear_efficiency

            ! Guide wire: teflon tape on fishing line, 80% contact chance
            CALL RANDOM_NUMBER(contact_chance)
            is_in_contact = (contact_chance <= 0.3)

            F_drive = current_torque / wheel_radius
            axle_normal_force = (mass * g) / bearing_count

            IF (is_in_contact) THEN
                F_touch = mu_wire * N_wire
            ELSE
                F_touch = 0.0
            END IF

            F_air = 0.5 * rho * Cd * A * current_speed**2
            F_bearing_total = mu_bearing * axle_normal_force * bearing_count
            F_rolling = mu_rolling * mass * g
            F_net = F_drive - (F_bearing_total + F_rolling + F_air + F_touch)

            acceleration = F_net / mass
            current_speed = current_speed + acceleration * time_step

            IF (current_speed < 0.0) THEN
                current_speed = 0.0
                zero_accel_count = zero_accel_count + 1
            ELSE IF (acceleration <= 0.0 .AND. current_speed == 0.0) THEN
                zero_accel_count = zero_accel_count + 1
            ELSE
                zero_accel_count = 0
            END IF

            position = position + current_speed * time_step
            total_time = total_time + time_step

            IF (total_time > MAX_SIMULATION_TIME) THEN
                failed_tests = failed_tests + 1
                EXIT simulation_loop
            END IF

            IF (ABS(position - last_position) < MIN_SPEED_THRESHOLD * time_step) THEN
                stagnation_counter = stagnation_counter + 1
                IF (stagnation_counter > STAGNATION_CHECKS) THEN
                    failed_tests = failed_tests + 1
                    EXIT simulation_loop
                END IF
            ELSE
                stagnation_counter = 0
                last_position = position
            END IF

            IF (zero_accel_count >= 100) THEN
                failed_tests = failed_tests + 1
                EXIT simulation_loop
            END IF

            current_rpm = ((current_speed * 60.0) / wheel_circumference) * gear_ratio
            current_rpm = MIN(current_rpm, motor_rpm_max)

        END DO simulation_loop

        IF (position >= distance_remaining) THEN
            total_time_sum = total_time_sum + total_time
            final_speed_sum = final_speed_sum + current_speed
            valid_tests = valid_tests + 1
            !$OMP CRITICAL
            IF (current_speed > max_speed) max_speed = current_speed
            IF (current_speed < min_speed) min_speed = current_speed
            IF (total_time > max_time) max_time = total_time
            IF (total_time < min_time) min_time = total_time
            !$OMP END CRITICAL
        END IF

    END DO
    !$OMP END PARALLEL DO

    IF (valid_tests > 0) THEN
        total_time_avg = total_time_sum / valid_tests
        final_speed_avg = final_speed_sum / valid_tests
    ELSE
        total_time_avg = 0.0
        final_speed_avg = 0.0
        max_speed = 0.0; min_speed = 0.0; max_time = 0.0; min_time = 0.0
    END IF

    PRINT *, '=== Final Results ==='
    WRITE(*, '(A,F6.4)') ' Gear ratio:      ', gear_ratio
    WRITE(*, '(A,F8.3,A)') ' Average Time:    ', total_time_avg, ' seconds'
    WRITE(*, '(A,F8.3,A)') ' Average Speed:   ', final_speed_avg, ' m/s'
    WRITE(*, '(A,F8.3,A)') ' Maximum Speed:   ', max_speed, ' m/s'
    WRITE(*, '(A,F8.3,A)') ' Minimum Speed:   ', min_speed, ' m/s'
    WRITE(*, '(A,F8.3,A)') ' Maximum Time:    ', max_time, ' seconds'
    WRITE(*, '(A,F8.3,A)') ' Minimum Time:    ', min_time, ' seconds'
    PRINT *,     'Failed Tests:    ', failed_tests
    PRINT *,     'Successful Tests:', valid_tests
    PRINT *,     'Total Iterations:', num_iterations

CONTAINS

    SUBROUTINE INIT_RANDOM_SEED()
        USE OMP_LIB
        IMPLICIT NONE
        INTEGER :: i, n, clock, tid
        INTEGER, ALLOCATABLE :: seed(:)

        tid = OMP_GET_THREAD_NUM()
        CALL RANDOM_SEED(SIZE = n)
        ALLOCATE(seed(n))
        CALL SYSTEM_CLOCK(COUNT=clock)
        seed = clock + 37 * (tid + 1)
        DO i = 1, n
            seed(i) = seed(i) + 37 * i
        END DO
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)
    END SUBROUTINE INIT_RANDOM_SEED

END PROGRAM SPEED_RATIO_MONTECARLO

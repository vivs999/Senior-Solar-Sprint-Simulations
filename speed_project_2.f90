PROGRAM SPEED_RATIO_MONTECARLO
    USE omp_lib
    IMPLICIT NONE

    ! Constants
    INTEGER, PARAMETER :: num_threads = 100
    INTEGER, PARAMETER :: num_iterations = 100000
    REAL, PARAMETER :: MAX_SIMULATION_TIME = 1000.0
    REAL, PARAMETER :: MIN_SPEED_THRESHOLD = 0.001
    INTEGER, PARAMETER :: STAGNATION_CHECKS = 500

    ! Simulation variables
    INTEGER :: i, bearing_count, zero_accel_count, failed_tests, valid_tests, stagnation_counter
    REAL :: g, wheel_radius, wheel_diameter, mass, gear_ratio
    REAL :: distance_remaining, time_step, total_time_sum, total_time_avg
    REAL :: final_speed_sum, final_speed_avg
    REAL :: motor_power_nominal, power_margin, stall_torque, free_speed_rpm
    REAL :: voltage_input, current_input, motor_efficiency, last_position

    ! Additional variables
    REAL :: motor_rpm_max, motor_torque_stall
    REAL :: current_rpm, current_torque, current_speed, position, total_time
    REAL :: axle_normal_force, F_bearing_total, F_rolling, F_net, acceleration
    REAL :: F_drive, wheel_circumference, mu_bearing, mu_rolling
    REAL :: power_variation, angular_velocity, current_power, motor_power_max
    REAL :: F_air, rho, Cd, A
    REAL :: contact_chance
    REAL :: F_touch
    LOGICAL :: is_in_contact
    REAL :: max_speed, min_speed
    REAL :: max_time, min_time

    ! Motor parameters
    stall_torque = 0.0127      ! 12.7 mN·m
    free_speed_rpm = 5200.0    ! No-load speed at 3V
    current_input = 1.55  ! Amperage (actual current draw)
    motor_efficiency = 0.8  ! Efficiency (80%)
    voltage_input = 3.0        ! Volts
    motor_power_nominal = voltage_input * current_input * motor_efficiency  ! Realistic mechanical power output
    power_margin = 0.15 * motor_power_nominal  ! 15% margin

    ! Initialize constants
    g = 9.81
    wheel_diameter = 0.0392     ! 39.2 mm diameter
    wheel_radius = wheel_diameter / 2.0
    bearing_count = 4
    distance_remaining = 20.0 ! meters
    time_step = 0.001
    mass = 0.57               ! kg
    gear_ratio = 50.0/22.0    ! gear ratio

    failed_tests = 0
    valid_tests = 0
    total_time_sum = 0.0
    final_speed_sum = 0.0
    max_speed = -1.0E30
    min_speed = 1.0E30
    max_time = -1.0E30
    min_time = 1.0E30

    rho = 1.225            ! Air density (kg/m^3)
    Cd = 0.09               ! Drag coefficient (typical value)
    A = 0.010               ! Frontal area (m^2)

    CALL OMP_SET_NUM_THREADS(num_threads)

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& PRIVATE(i, mu_bearing, mu_rolling, motor_rpm_max, motor_torque_stall, &
    !$OMP& current_rpm, current_speed, position, total_time, wheel_circumference, &
    !$OMP& current_torque, F_drive, axle_normal_force, F_bearing_total, F_rolling, &
    !$OMP& F_net, acceleration, zero_accel_count, power_variation, current_power, F_touch, &
    !$OMP& angular_velocity, motor_power_max, last_position, stagnation_counter, F_air, contact_chance, is_in_contact, &
    !$OMP& current_input, motor_efficiency) &
    !$OMP& SHARED(motor_power_nominal, power_margin, g, wheel_diameter, wheel_radius, &
    !$OMP& bearing_count, distance_remaining, time_step, mass, gear_ratio, rho, Cd, A, max_speed, min_speed, max_time, min_time) &
    !$OMP& REDUCTION(+:total_time_sum, final_speed_sum, valid_tests, failed_tests)
    DO i = 1, num_iterations
        CALL INIT_RANDOM_SEED()
        
        ! Randomize parameters
        CALL RANDOM_NUMBER(mu_bearing)
        mu_bearing = 0.05 + mu_bearing * 0.02          ! 0.05-0.07
        
        CALL RANDOM_NUMBER(mu_rolling)
        mu_rolling = 0.0075 + mu_rolling * 0.005          ! 0.01-0.02
        
        CALL RANDOM_NUMBER(motor_rpm_max)
        motor_rpm_max = 10500.0 + motor_rpm_max * 1000.0  ! 10500-11500 RPM
        
        CALL RANDOM_NUMBER(motor_torque_stall)
        motor_torque_stall = 0.012 + motor_torque_stall * 0.002  ! 0.012-0.014 N·m
        
        CALL RANDOM_NUMBER(power_variation)
        motor_power_max = motor_power_nominal - power_margin + 2.0 * power_margin * power_variation

        ! Reset simulation variables
        current_rpm = 0.0
        current_speed = 0.0
        position = 0.0
        total_time = 0.0
        zero_accel_count = 0
        last_position = 0.0
        stagnation_counter = 0
        wheel_circumference = 3.14159 * wheel_diameter

        simulation_loop: DO WHILE (position < distance_remaining)
            ! Calculate current torque
            current_torque = motor_torque_stall * (1.0 - current_rpm / motor_rpm_max)
            
            ! Convert RPM to angular velocity (rad/s)
            angular_velocity = (2.0 * 3.14159 * current_rpm) / 60.0
            
            ! Power limit enforcement
            if (angular_velocity > 0.0) then
                current_power = current_torque * angular_velocity
                if (current_power > motor_power_max) then
                    current_torque = motor_power_max / angular_velocity
                end if
            end if
            ! Calculate forces
            current_torque = current_torque * gear_ratio
            
            CALL RANDOM_NUMBER(contact_chance)
            is_in_contact = (contact_chance <= 0.1)

            F_drive = current_torque / wheel_radius
            axle_normal_force = (mass * g) / bearing_count
            
            IF (is_in_contact) THEN
                F_touch = 0.1
            ELSE
                F_touch = 0.0
            END IF
            
            ! Air resistance: F_air = 0.5 * rho * Cd * A * v^2
            F_air = 0.5 * rho * Cd * A * current_speed**2
            
            F_bearing_total = mu_bearing * axle_normal_force * bearing_count
            F_rolling = mu_rolling * mass * g

            F_net = F_drive - (F_bearing_total + F_rolling + F_air + F_touch)

            

            ! Update physics
            acceleration = F_net / mass
            current_speed = current_speed + acceleration * time_step

            ! Handle negative speed and stagnation
            if (current_speed < 0.0) then
                current_speed = 0.0
                zero_accel_count = zero_accel_count + 1
            else if (acceleration <= 0.0 .AND. current_speed == 0.0) then
                zero_accel_count = zero_accel_count + 1
            else
                zero_accel_count = 0
            end if

            ! Update position and time
            position = position + current_speed * time_step
            total_time = total_time + time_step

            ! Infinite loop safeguards
            if (total_time > MAX_SIMULATION_TIME) then
                failed_tests = failed_tests + 1
                exit simulation_loop
            end if
            
            if (ABS(position - last_position) < MIN_SPEED_THRESHOLD * time_step) then
                stagnation_counter = stagnation_counter + 1
                if (stagnation_counter > STAGNATION_CHECKS) then
                    failed_tests = failed_tests + 1
                    exit simulation_loop
                end if
            else
                stagnation_counter = 0
                last_position = position
            end if

            ! Failure condition check
            if (zero_accel_count >= 100) then
                failed_tests = failed_tests + 1
                exit simulation_loop
            end if

            ! Update RPM
            current_rpm = ((current_speed * 60.0) / wheel_circumference) * gear_ratio
            current_rpm = MIN(current_rpm, motor_rpm_max)

        END DO simulation_loop

        ! Record successful tests
        if (position >= distance_remaining) then
            total_time_sum = total_time_sum + total_time
            final_speed_sum = final_speed_sum + current_speed
            valid_tests = valid_tests + 1
            !$OMP CRITICAL
            if (current_speed > max_speed) max_speed = current_speed
            if (current_speed < min_speed) min_speed = current_speed
            if (total_time > max_time) max_time = total_time
            if (total_time < min_time) min_time = total_time
            !$OMP END CRITICAL
        end if


    END DO
    !$OMP END PARALLEL DO

    ! Calculate averages
    if (valid_tests > 0) then
        total_time_avg = total_time_sum / valid_tests
        final_speed_avg = final_speed_sum / valid_tests
    else
        total_time_avg = 0.0
        final_speed_avg = 0.0
    end if

    if (valid_tests == 0) then
        max_speed = 0.0
        min_speed = 0.0
        max_time = 0.0
        min_time = 0.0
    end if

    ! Output results
    PRINT *, '=== Final Results ==='
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

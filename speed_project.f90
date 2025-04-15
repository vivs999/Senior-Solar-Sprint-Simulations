PROGRAM SPEED_RATIO_MONTECARLO
    USE omp_lib
    IMPLICIT NONE

    ! Declare Constants and variables
    INTEGER, PARAMETER :: num_threads = 100
    INTEGER, PARAMETER :: num_iterations = 1000
    INTEGER :: i
    REAL :: g, wheel_radius, wheel_diameter, mass, motor_rpm_max, motor_torque_stall, motor_gear, wheel_gear
    REAL :: gear_ratio
    REAL :: current_rpm, current_torque, current_speed, position, time_step
    REAL :: distance_remaining, total_time, axle_normal_force, F_bearing_total
    REAL :: F_rolling, F_net, acceleration, F_drive, wheel_circumference
    REAL :: mu_bearing, mu_rolling
    REAL :: total_time_sum, total_time_avg, final_speed_sum, final_speed_avg
    INTEGER :: bearing_count, zero_accel_count, failed_tests, valid_tests

    ! Motor specifications
    REAL :: motor_power_nominal, power_margin, power_variation
    REAL :: stall_torque, free_speed_rpm, current_efficiency, voltage_input
    REAL :: angular_velocity, current_power, motor_power_max

    stall_torque = 0.0127  ! 12.7 mN·m stall torque
    free_speed_rpm = 5200.0  ! No-load speed at 3V
    current_efficiency = 1.55  ! Amperage at max efficiency at 3V
    voltage_input = 3.0  ! Input voltage
    motor_power_nominal = voltage_input * current_efficiency  ! 4.65 W
    power_margin = 0.15 * motor_power_nominal  ! 15% margin

    ! Initialize constants
    g = 9.81
    wheel_diameter = 0.0392
    wheel_radius = wheel_diameter / 2.0
    bearing_count = 4
    distance_remaining = 4000.0
    time_step = 0.001
    mass = 1.2

    motor_gear = 1.0
    wheel_gear = 50.0
    gear_ratio = wheel_gear / motor_gear

    failed_tests = 0
    valid_tests = 0

    total_time_sum = 0.0
    final_speed_sum = 0.0

    CALL RANDOM_SEED()
    CALL OMP_SET_NUM_THREADS(num_threads)

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& PRIVATE(i, mu_bearing, mu_rolling, motor_rpm_max, motor_torque_stall, current_rpm, current_speed, position, total_time,  wheel_circumference, current_torque, F_drive, axle_normal_force, F_bearing_total, F_rolling, F_net, acceleration, zero_accel_count, power_variation, current_power, angular_velocity, motor_power_max) &
    !$OMP& SHARED(motor_power_nominal, power_margin, g, wheel_diameter, wheel_radius, bearing_count, distance_remaining, time_step, mass, motor_gear, wheel_gear, gear_ratio) &
    !$OMP& REDUCTION(+:total_time_sum, final_speed_sum, valid_tests, failed_tests)
    DO i = 1, num_iterations

        CALL RANDOM_NUMBER(mu_bearing)
        mu_bearing = 0.05 + mu_bearing * 0.02

        CALL RANDOM_NUMBER(mu_rolling)
        mu_rolling = 0.01 + mu_rolling * 0.01

        CALL RANDOM_NUMBER(motor_rpm_max)
        motor_rpm_max = 10500.0 + motor_rpm_max * 1000.0

        CALL RANDOM_NUMBER(motor_torque_stall)
        motor_torque_stall = 0.012 + motor_torque_stall * 0.002

        CALL RANDOM_NUMBER(power_variation)
        motor_power_max = motor_power_nominal - power_margin + 2.0 * power_margin * power_variation

        current_rpm = 0.0
        current_speed = 0.0
        position = 0.0
        total_time = 0.0
        zero_accel_count = 0

        wheel_circumference = 3.14159 * wheel_diameter

        DO WHILE (position < distance_remaining)

            current_torque = motor_torque_stall * (1.0 - (current_rpm / motor_rpm_max)**2) * gear_ratio

            angular_velocity = (2.0 * 3.14159 * current_rpm) / 60.0
            current_power = current_torque * angular_velocity
            IF (current_power > motor_power_max .AND. angular_velocity > 0.0) THEN
                current_torque = motor_power_max / angular_velocity
            END IF

            F_drive = current_torque / wheel_radius

            axle_normal_force = (mass * g) / bearing_count
            F_bearing_total = mu_bearing * axle_normal_force * bearing_count
            F_rolling = mu_rolling * mass * g
            F_net = F_drive - (F_bearing_total + F_rolling)

            IF (F_net > 0.0) THEN
                acceleration = F_net / mass
                zero_accel_count = 0
            ELSE
                acceleration = 0.0
                zero_accel_count = zero_accel_count + 1
            END IF

            IF (zero_accel_count >= 10) THEN
                failed_tests = failed_tests + 1
                EXIT
            END IF

            current_speed = current_speed + acceleration * time_step
            position = position + current_speed * time_step

            current_rpm = ((current_speed * 60.0) / wheel_circumference) * gear_ratio
            IF (current_rpm >= motor_rpm_max) THEN
                current_rpm = motor_rpm_max
            END IF

            total_time = total_time + time_step

            IF (position >= distance_remaining) EXIT

        END DO

        IF (zero_accel_count < 10) THEN
            total_time_sum = total_time_sum + total_time
            final_speed_sum = final_speed_sum + current_speed
            valid_tests = valid_tests + 1
        END IF

    END DO
    !$OMP END PARALLEL DO

    total_time_avg = total_time_sum / valid_tests
    final_speed_avg = final_speed_sum / valid_tests

    PRINT *, 'Monte Carlo Simulation Results:'
    PRINT *, 'Average Total Time to Cross the Track: ', total_time_avg, ' seconds'
    PRINT *, 'Average Final Speed: ', final_speed_avg, ' m/s'
    PRINT *, 'Failed tests: ', failed_tests
    PRINT *, 'Number of iterations: ', num_iterations

END PROGRAM SPEED_RATIO_MONTECARLO

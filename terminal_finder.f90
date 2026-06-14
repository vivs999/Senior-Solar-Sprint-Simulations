PROGRAM TERMINAL_FINDER
    IMPLICIT NONE

    ! Fixed mean parameters (same as current friction_analysis / speed_project_2)
    REAL, PARAMETER :: g               = 9.81
    REAL, PARAMETER :: mass            = 0.297
    REAL, PARAMETER :: wheel_diameter  = 0.0392
    REAL, PARAMETER :: gear_ratio      = 1.5
    REAL, PARAMETER :: gear_efficiency = 0.98
    REAL, PARAMETER :: bearing_count   = 4.0
    REAL, PARAMETER :: distance        = 500.0   ! extended to find terminal v
    REAL, PARAMETER :: time_step       = 0.0001   ! finer step for precision

    ! Motor (Pitsco 280 mean values)
    REAL, PARAMETER :: motor_rpm_max     = 5200.0
    REAL, PARAMETER :: motor_torque_stall= 0.013
    REAL, PARAMETER :: voltage_input     = 3.0
    REAL, PARAMETER :: current_input     = 1.55
    REAL, PARAMETER :: motor_efficiency  = 0.60

    ! Friction (mean of ranges)
    REAL, PARAMETER :: mu_bearing        = 0.00125
    REAL, PARAMETER :: mu_rolling        = 0.006
    REAL, PARAMETER :: mu_wire           = 0.05
    REAL, PARAMETER :: N_wire            = 0.17
    REAL, PARAMETER :: wire_contact_prob = 0.3

    ! Aero
    REAL, PARAMETER :: rho = 1.225
    REAL, PARAMETER :: Cd  = 1.05
    REAL, PARAMETER :: A   = 0.0025

    ! Detection threshold: acceleration considered zero below this
    REAL, PARAMETER :: ACCEL_THRESHOLD  = 0.005   ! m/s²
    REAL, PARAMETER :: SPEED_TOL        = 0.0001  ! m/s change per step to confirm flatline

    REAL :: wheel_radius, wheel_circumference, motor_power_max
    REAL :: axle_normal_force, F_wire_avg

    REAL :: current_rpm, current_speed, position, total_time
    REAL :: current_torque, angular_velocity, current_power
    REAL :: F_drive, F_bearing, F_rolling, F_air, F_net, acceleration

    REAL :: terminal_speed, terminal_pos, terminal_time
    LOGICAL :: flatline_found

    REAL :: prev_speed
    INTEGER :: flatline_steps

    wheel_radius        = wheel_diameter / 2.0
    wheel_circumference = 3.14159 * wheel_diameter
    motor_power_max     = voltage_input * current_input * motor_efficiency
    axle_normal_force   = (mass * g) / bearing_count
    F_wire_avg          = mu_wire * N_wire * wire_contact_prob

    current_rpm   = 0.0
    current_speed = 0.0
    position      = 0.0
    total_time    = 0.0
    prev_speed    = 0.0
    flatline_steps = 0
    flatline_found = .FALSE.
    terminal_speed = 0.0
    terminal_pos   = 0.0
    terminal_time  = 0.0

    PRINT *, '=== Terminal Velocity Finder ==='
    PRINT *, ''
    WRITE(*, '(A8, A10, A12, A12, A12, A12)') &
        'Time(s)', 'Pos(m)', 'Speed(m/s)', 'Accel(m/s2)', 'F_drive(N)', 'F_drag(N)'
    WRITE(*, '(A)') REPEAT('-', 72)

    DO WHILE (position < distance)
        current_torque = motor_torque_stall * (1.0 - current_rpm / motor_rpm_max)

        angular_velocity = (2.0 * 3.14159 * current_rpm) / 60.0
        IF (angular_velocity > 0.0) THEN
            current_power = current_torque * angular_velocity
            IF (current_power > motor_power_max) THEN
                current_torque = motor_power_max / angular_velocity
            END IF
        END IF

        current_torque = current_torque * gear_ratio * gear_efficiency
        F_drive   = current_torque / wheel_radius

        F_bearing = mu_bearing * axle_normal_force * bearing_count
        F_rolling = mu_rolling * mass * g
        F_air     = 0.5 * rho * Cd * A * current_speed**2
        F_net     = F_drive - (F_bearing + F_rolling + F_air + F_wire_avg)

        acceleration = F_net / mass

        ! Print every 500 steps (0.05 s intervals)
        IF (MOD(NINT(total_time / time_step), 500) == 0) THEN
            WRITE(*, '(F8.3, F10.3, F12.4, F12.5, F12.5, F12.5)') &
                total_time, position, current_speed, acceleration, &
                F_drive, (F_bearing + F_rolling + F_air + F_wire_avg)
        END IF

        ! Detect flatline: acceleration near zero for 100 consecutive steps
        IF (.NOT. flatline_found) THEN
            IF (ABS(acceleration) < ACCEL_THRESHOLD .AND. &
                ABS(current_speed - prev_speed) < SPEED_TOL) THEN
                flatline_steps = flatline_steps + 1
                IF (flatline_steps >= 100) THEN
                    flatline_found = .TRUE.
                    terminal_speed = current_speed
                    terminal_pos   = position
                    terminal_time  = total_time
                END IF
            ELSE
                flatline_steps = 0
            END IF
        END IF

        prev_speed    = current_speed
        acceleration  = F_net / mass
        IF (current_speed + acceleration * time_step < 0.0) current_speed = 0.0
        current_speed = current_speed + acceleration * time_step
        position      = position + current_speed * time_step
        total_time    = total_time + time_step

        current_rpm = ((current_speed * 60.0) / wheel_circumference) * gear_ratio
        current_rpm = MIN(current_rpm, motor_rpm_max)

        IF (flatline_found) EXIT

        IF (total_time > 600.0) THEN
            PRINT *, 'ERROR: simulation timeout'
            STOP
        END IF
    END DO

    WRITE(*, '(A)') REPEAT('-', 72)
    PRINT *, ''
    PRINT *, '=== Terminal Velocity Results ==='
    IF (flatline_found) THEN
        WRITE(*, '(A,F7.4,A)') ' Terminal velocity reached at: ', terminal_pos,  ' m'
        WRITE(*, '(A,F7.3,A)') ' Time at terminal velocity:    ', terminal_time, ' s'
        WRITE(*, '(A,F7.4,A)') ' Terminal speed:               ', terminal_speed, ' m/s'
        WRITE(*, '(A,F7.4,A)') ' Distance remaining after:     ', distance - terminal_pos, ' m'
    ELSE
        PRINT *, ' Terminal velocity NOT reached within 20 m run.'
        WRITE(*, '(A,F7.4,A)') ' Final speed at finish:        ', current_speed, ' m/s'
    END IF
    PRINT *, ''
    WRITE(*, '(A,F7.3,A)') ' Total finish time:            ', total_time, ' s'
    WRITE(*, '(A,F7.4,A)') ' Final speed:                  ', current_speed, ' m/s'

END PROGRAM TERMINAL_FINDER

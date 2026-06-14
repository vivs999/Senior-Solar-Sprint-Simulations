PROGRAM FRICTION_ANALYSIS
    IMPLICIT NONE

    ! Fixed mean parameters (no Monte Carlo randomisation)
    REAL, PARAMETER :: g             = 9.81
    REAL, PARAMETER :: mass          = 0.297
    REAL, PARAMETER :: wheel_diameter= 0.0392
    REAL, PARAMETER :: gear_ratio    = 1.5
    REAL, PARAMETER :: gear_efficiency = 0.98
    REAL, PARAMETER :: bearing_count = 4.0
    REAL, PARAMETER :: distance      = 20.0
    REAL, PARAMETER :: time_step     = 0.001
    INTEGER, PARAMETER :: PRINT_EVERY = 100       ! print every 100 timesteps

    ! Motor (Pitsco 280 mean values)
    REAL, PARAMETER :: motor_rpm_max    = 5200.0
    REAL, PARAMETER :: motor_torque_stall = 0.013
    REAL, PARAMETER :: voltage_input   = 3.0
    REAL, PARAMETER :: current_input   = 1.55
    REAL, PARAMETER :: motor_efficiency= 0.60
    REAL, PARAMETER :: power_margin_pct= 0.15

    ! Friction parameters (mean of ranges)
    REAL, PARAMETER :: mu_bearing = 0.00125   ! ceramic midpoint 0.0005-0.002
    REAL, PARAMETER :: mu_rolling = 0.006     ! latex on paper midpoint 0.004-0.008
    REAL, PARAMETER :: mu_wire    = 0.05
    REAL, PARAMETER :: N_wire     = 0.17
    REAL, PARAMETER :: wire_contact_prob = 0.3  ! expected value: 30% contact

    ! Aero
    REAL, PARAMETER :: rho = 1.225
    REAL, PARAMETER :: Cd  = 1.05          ! rectangular box bluff body (25mm x 100mm)
    REAL, PARAMETER :: A   = 0.0025        ! box frontal area: 25mm x 100mm = 2500 mm²

    ! Derived constants
    REAL :: wheel_radius, wheel_circumference, motor_power_max
    REAL :: axle_normal_force

    ! State variables
    REAL :: current_rpm, current_speed, position, total_time
    REAL :: current_torque, angular_velocity, current_power
    REAL :: F_drive, F_bearing, F_rolling, F_air, F_wire_avg, F_net, acceleration

    ! Accumulators for averages
    REAL :: sum_F_bearing, sum_F_rolling, sum_F_air, sum_F_wire
    REAL :: sum_F_drive, sum_F_net
    INTEGER :: step_count

    wheel_radius       = wheel_diameter / 2.0
    wheel_circumference= 3.14159 * wheel_diameter
    motor_power_max    = voltage_input * current_input * motor_efficiency
    axle_normal_force  = (mass * g) / bearing_count

    ! Wire friction is stochastic — use expected value for deterministic run
    F_wire_avg = mu_wire * N_wire * wire_contact_prob

    ! Initialise state
    current_rpm   = 0.0
    current_speed = 0.0
    position      = 0.0
    total_time    = 0.0
    sum_F_bearing = 0.0
    sum_F_rolling = 0.0
    sum_F_air     = 0.0
    sum_F_wire    = 0.0
    sum_F_drive   = 0.0
    sum_F_net     = 0.0
    step_count    = 0

    WRITE(*, '(A5, A8, A10, A10, A10, A10, A10, A10, A8)') &
        'Step', 'Time(s)', 'Pos(m)', 'Speed', 'F_drive', &
        'F_bear', 'F_roll', 'F_air', 'F_wire'
    WRITE(*, '(A)') REPEAT('-', 85)

    DO WHILE (position < distance)
        ! Motor torque curve (linear)
        current_torque = motor_torque_stall * (1.0 - current_rpm / motor_rpm_max)

        ! Power cap
        angular_velocity = (2.0 * 3.14159 * current_rpm) / 60.0
        IF (angular_velocity > 0.0) THEN
            current_power = current_torque * angular_velocity
            IF (current_power > motor_power_max) THEN
                current_torque = motor_power_max / angular_velocity
            END IF
        END IF

        current_torque = current_torque * gear_ratio * gear_efficiency
        F_drive   = current_torque / wheel_radius

        ! Friction forces
        F_bearing = mu_bearing * axle_normal_force * bearing_count
        F_rolling = mu_rolling * mass * g
        F_air     = 0.5 * rho * Cd * A * current_speed**2

        F_net = F_drive - (F_bearing + F_rolling + F_air + F_wire_avg)

        ! Update accumulators
        sum_F_bearing = sum_F_bearing + F_bearing
        sum_F_rolling = sum_F_rolling + F_rolling
        sum_F_air     = sum_F_air     + F_air
        sum_F_wire    = sum_F_wire    + F_wire_avg
        sum_F_drive   = sum_F_drive   + F_drive
        sum_F_net     = sum_F_net     + F_net
        step_count    = step_count    + 1

        ! Print every PRINT_EVERY steps
        IF (MOD(step_count, PRINT_EVERY) == 0) THEN
            WRITE(*, '(I5, F8.3, F10.3, F10.4, F10.5, F10.6, F10.5, F10.5, F10.6)') &
                step_count, total_time, position, current_speed, &
                F_drive, F_bearing, F_rolling, F_air, F_wire_avg
        END IF

        ! Integrate
        acceleration  = F_net / mass
        IF (current_speed + acceleration * time_step < 0.0) current_speed = 0.0
        current_speed = current_speed + acceleration * time_step
        position      = position + current_speed * time_step
        total_time    = total_time + time_step

        current_rpm = ((current_speed * 60.0) / wheel_circumference) * gear_ratio
        current_rpm = MIN(current_rpm, motor_rpm_max)

        IF (total_time > 60.0) THEN
            PRINT *, 'ERROR: simulation timeout'
            STOP
        END IF
    END DO

    ! Summary
    WRITE(*, '(A)') REPEAT('-', 85)
    PRINT *, ''
    PRINT *, '=== Run Summary ==='
    WRITE(*, '(A,F7.3,A)') ' Finish time:          ', total_time, ' s'
    WRITE(*, '(A,F7.4,A)') ' Final speed:          ', current_speed, ' m/s'
    PRINT *, ''
    PRINT *, '=== Average Forces Over Full Run ==='
    WRITE(*, '(A,F10.6,A,F6.2,A)') ' F_drive   (motor):  ', sum_F_drive   / step_count, ' N  (', &
        100.0 * (sum_F_drive / step_count) / (sum_F_drive / step_count), '%)'
    WRITE(*, '(A,F10.6,A,F6.2,A)') ' F_bearing:          ', sum_F_bearing / step_count, ' N'
    WRITE(*, '(A,F10.6,A,F6.2,A)') ' F_rolling:          ', sum_F_rolling / step_count, ' N'
    WRITE(*, '(A,F10.6,A,F6.2,A)') ' F_air:              ', sum_F_air     / step_count, ' N'
    WRITE(*, '(A,F10.6,A,F6.2,A)') ' F_wire (expected):  ', sum_F_wire    / step_count, ' N'

    PRINT *, ''
    PRINT *, '=== Friction Share of Total Drag ==='
    BLOCK
        REAL :: total_drag, fb, fr, fa, fw
        fb = sum_F_bearing / step_count
        fr = sum_F_rolling / step_count
        fa = sum_F_air     / step_count
        fw = sum_F_wire    / step_count
        total_drag = fb + fr + fa + fw
        WRITE(*, '(A,F6.2,A)') ' Bearing:    ', 100.0 * fb / total_drag, '%'
        WRITE(*, '(A,F6.2,A)') ' Rolling:    ', 100.0 * fr / total_drag, '%'
        WRITE(*, '(A,F6.2,A)') ' Air drag:   ', 100.0 * fa / total_drag, '%'
        WRITE(*, '(A,F6.2,A)') ' Guide wire: ', 100.0 * fw / total_drag, '%'
        WRITE(*, '(A,F8.6,A)') ' Total drag: ', total_drag, ' N (avg)'
    END BLOCK

END PROGRAM FRICTION_ANALYSIS

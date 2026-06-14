PROGRAM SENSITIVITY
    IMPLICIT NONE

    REAL, PARAMETER :: g               = 9.81
    REAL, PARAMETER :: wheel_diameter  = 0.0392
    REAL, PARAMETER :: distance        = 20.0
    REAL, PARAMETER :: time_step       = 0.001
    REAL, PARAMETER :: PERTURB         = 0.10   ! 10% perturbation

    ! Baseline parameters
    REAL, PARAMETER :: BASE_mass            = 0.210
    REAL, PARAMETER :: BASE_gear_ratio      = 1.5
    REAL, PARAMETER :: BASE_gear_efficiency = 0.98
    REAL, PARAMETER :: BASE_bearing_count   = 4.0
    REAL, PARAMETER :: BASE_mu_bearing      = 0.00125
    REAL, PARAMETER :: BASE_mu_rolling      = 0.006
    REAL, PARAMETER :: BASE_mu_wire         = 0.05
    REAL, PARAMETER :: BASE_N_wire          = 0.17
    REAL, PARAMETER :: BASE_wire_prob       = 0.3
    REAL, PARAMETER :: BASE_motor_rpm_max   = 5200.0
    REAL, PARAMETER :: BASE_motor_torque    = 0.013
    REAL, PARAMETER :: BASE_rho             = 1.225
    REAL, PARAMETER :: BASE_Cd              = 1.05
    REAL, PARAMETER :: BASE_A               = 0.0025

    INTEGER, PARAMETER :: N_PARAMS = 9
    CHARACTER(LEN=20) :: param_names(N_PARAMS)
    REAL :: base_time, hi_time, lo_time, sens(N_PARAMS), abs_sens(N_PARAMS)
    INTEGER :: order(N_PARAMS), i, j, tmp_idx

    param_names(1) = 'mass'
    param_names(2) = 'Cd (drag coeff)'
    param_names(3) = 'frontal area A'
    param_names(4) = 'motor torque'
    param_names(5) = 'motor rpm max'
    param_names(6) = 'gear efficiency'
    param_names(7) = 'mu_rolling'
    param_names(8) = 'mu_bearing'
    param_names(9) = 'guide wire frict'

    base_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                        BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, &
                        BASE_N_wire, BASE_wire_prob, BASE_motor_rpm_max, &
                        BASE_motor_torque, BASE_Cd, BASE_A)

    ! 1. mass
    hi_time = run_sim(BASE_mass*(1+PERTURB), BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    lo_time = run_sim(BASE_mass*(1-PERTURB), BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    sens(1) = (hi_time - lo_time) / (2.0 * PERTURB * base_time) * 100.0

    ! 2. Cd
    hi_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd*(1+PERTURB), BASE_A)
    lo_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd*(1-PERTURB), BASE_A)
    sens(2) = (hi_time - lo_time) / (2.0 * PERTURB * base_time) * 100.0

    ! 3. frontal area A
    hi_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A*(1+PERTURB))
    lo_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A*(1-PERTURB))
    sens(3) = (hi_time - lo_time) / (2.0 * PERTURB * base_time) * 100.0

    ! 4. motor torque
    hi_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque*(1+PERTURB), BASE_Cd, BASE_A)
    lo_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque*(1-PERTURB), BASE_Cd, BASE_A)
    sens(4) = (hi_time - lo_time) / (2.0 * PERTURB * base_time) * 100.0

    ! 5. motor rpm max
    hi_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max*(1+PERTURB), BASE_motor_torque, BASE_Cd, BASE_A)
    lo_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max*(1-PERTURB), BASE_motor_torque, BASE_Cd, BASE_A)
    sens(5) = (hi_time - lo_time) / (2.0 * PERTURB * base_time) * 100.0

    ! 6. gear efficiency
    hi_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency*(1+PERTURB), &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    lo_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency*(1-PERTURB), &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    sens(6) = (hi_time - lo_time) / (2.0 * PERTURB * base_time) * 100.0

    ! 7. mu_rolling
    hi_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling*(1+PERTURB), BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    lo_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling*(1-PERTURB), BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    sens(7) = (hi_time - lo_time) / (2.0 * PERTURB * base_time) * 100.0

    ! 8. mu_bearing
    hi_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing*(1+PERTURB), BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    lo_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing*(1-PERTURB), BASE_mu_rolling, BASE_mu_wire, BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    sens(8) = (hi_time - lo_time) / (2.0 * PERTURB * base_time) * 100.0

    ! 9. guide wire friction (mu_wire * N_wire * prob)
    hi_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire*(1+PERTURB), BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    lo_time = run_sim(BASE_mass, BASE_gear_ratio, BASE_gear_efficiency, &
                      BASE_mu_bearing, BASE_mu_rolling, BASE_mu_wire*(1-PERTURB), BASE_N_wire, &
                      BASE_wire_prob, BASE_motor_rpm_max, BASE_motor_torque, BASE_Cd, BASE_A)
    sens(9) = (hi_time - lo_time) / (2.0 * PERTURB * base_time) * 100.0

    ! Sort by absolute sensitivity descending (bubble sort)
    DO i = 1, N_PARAMS
        order(i) = i
        abs_sens(i) = ABS(sens(i))
    END DO
    DO i = 1, N_PARAMS-1
        DO j = i+1, N_PARAMS
            IF (abs_sens(order(j)) > abs_sens(order(i))) THEN
                tmp_idx  = order(i)
                order(i) = order(j)
                order(j) = tmp_idx
            END IF
        END DO
    END DO

    WRITE(*, '(A)') ''
    WRITE(*, '(A,F6.3,A)') ' Baseline finish time: ', base_time, ' s  (mass=210g, ratio=1.5)'
    WRITE(*, '(A)') ''
    WRITE(*, '(A)') ' === Sensitivity: % change in finish time per 10% param improvement ==='
    WRITE(*, '(A)') '  Rank  Parameter             Sensitivity   Direction'
    WRITE(*, '(A)') '  ----  --------------------  -----------   ---------'
    DO i = 1, N_PARAMS
        j = order(i)
        WRITE(*, '(A,I2,A,A20,F8.2,A,A)') '   ', i, '   ', param_names(j), &
            -sens(j), '%', &
            MERGE('reduce to improve    ', 'increase to improve  ', sens(j) > 0.0)
    END DO
    WRITE(*, '(A)') ''
    WRITE(*, '(A)') '  (negative sensitivity = increasing param hurts time; shown as gain from reducing)'

CONTAINS

    REAL FUNCTION run_sim(mass, gear_ratio, gear_eff, mu_bear, mu_roll, &
                          mu_wire, N_wire, wire_prob, rpm_max, torque_stall, Cd, A)
        REAL, INTENT(IN) :: mass, gear_ratio, gear_eff, mu_bear, mu_roll
        REAL, INTENT(IN) :: mu_wire, N_wire, wire_prob, rpm_max, torque_stall, Cd, A

        REAL, PARAMETER :: g_loc          = 9.81
        REAL, PARAMETER :: wheel_diam     = 0.0392
        REAL, PARAMETER :: dist           = 20.0
        REAL, PARAMETER :: dt             = 0.001
        REAL, PARAMETER :: voltage        = 3.0
        REAL, PARAMETER :: current        = 1.55
        REAL, PARAMETER :: motor_eff      = 0.60
        REAL, PARAMETER :: bear_count     = 4.0

        REAL :: wheel_r, wheel_circ, power_max
        REAL :: spd, pos, t, rpm, torque, omega, power
        REAL :: F_drive, F_bear, F_roll, F_air, F_wire, F_net, accel

        wheel_r    = wheel_diam / 2.0
        wheel_circ = 3.14159265 * wheel_diam
        power_max  = voltage * current * motor_eff

        spd = 0.0; pos = 0.0; t = 0.0; rpm = 0.0

        DO WHILE (pos < dist)
            torque = torque_stall * (1.0 - rpm / rpm_max)
            omega  = (2.0 * 3.14159265 * rpm) / 60.0
            IF (omega > 0.0) THEN
                power = torque * omega
                IF (power > power_max) torque = power_max / omega
            END IF
            torque  = torque * gear_ratio * gear_eff
            F_drive = torque / wheel_r
            F_bear  = mu_bear * (mass * g_loc / bear_count) * bear_count
            F_roll  = mu_roll * mass * g_loc
            F_air   = 0.5 * 1.225 * Cd * A * spd**2
            F_wire  = mu_wire * N_wire * wire_prob
            F_net   = F_drive - (F_bear + F_roll + F_air + F_wire)
            accel   = F_net / mass
            IF (spd + accel * dt < 0.0) spd = 0.0
            spd = spd + accel * dt
            pos = pos + spd * dt
            t   = t   + dt
            rpm = MIN(((spd * 60.0) / wheel_circ) * gear_ratio, rpm_max)
            IF (t > 60.0) EXIT
        END DO
        run_sim = t
    END FUNCTION run_sim

END PROGRAM SENSITIVITY

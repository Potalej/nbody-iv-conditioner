!! TEST 1: CONDITIONERS FOR FIRST INTEGRALS
!!
!! here we test every conditioner separately to guarantee that
!! all works.
!!
SUBROUTINE test1 (ier)
    USE conditioners_mod
    USE utils_mod
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: ier

    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)

    REAL(pf) :: tol = 1e-10

    INTEGER :: N = 100
    REAL(pf), ALLOCATABLE :: ms(:), qs(:,:), ps(:,:)
    REAL(pf) :: G = 1.0_pf, soft

    REAL(pf) :: P_desired(3)
    REAL(pf) :: J_desired(3)
    REAL(pf) :: E_desired

    REAL(pf) :: com_final(3)
    REAL(pf) :: P_final(3)
    REAL(pf) :: J_final(3)
    REAL(pf) :: E_final

    REAL(pf) :: com_error
    REAL(pf) :: P_error
    REAL(pf) :: J_error
    REAL(pf) :: E_error
    REAL(pf) :: error

    ALLOCATE(ms(N))
    ALLOCATE(qs(3,N))
    ALLOCATE(ps(3,N))

    ! generating state vectors
    ms = 1.0_pf / N
    CALL RANDOM_NUMBER(qs)
    qs = 2.0_pf * qs - 1.0_pf
    CALL RANDOM_NUMBER(ps)

print *, ''
print *, '----------------------------------------------'
print *, '# TEST 1: separated conditioners'
print *, '----------------------------------------------'

    print *, '# CENTER OF MASS'
    print *, '[COM] original:   ', center_of_mass(ms, qs)
    CALL comori(ms,qs)
    com_final = center_of_mass(ms, qs)
    print *, '[COM] conditioned:', com_final

    print *, '----------------------------------------------'

    print *, '# TOTAL LINEAR MOMENTUM'
    print *, '[TLM] original:   ', total_linear_momentum(ps)
    P_desired(1) = 0.0_pf
    P_desired(2) = 0.0_pf
    P_desired(3) = 0.0_pf
    CALL contotlinmom(P_desired, ps, ms)
    P_final = total_linear_momentum(ps)
    print *, '[TLM] conditioned:', P_final

    print *, '----------------------------------------------'

    print *, '# TOTAL ANGULAR MOMENTUM'
    print *, '[TAM] original:   ', total_angular_momentum(qs, ps)
    J_desired(1) = 0.0_pf
    J_desired(2) = 0.0_pf
    J_desired(3) = 0.0_pf
    CALL contotangmom(J_desired, ms, qs, ps, ier)
    J_final = total_angular_momentum(qs, ps)
    print *, '[TAM] conditioned:', J_final

    print *, '----------------------------------------------'

    print *, '# TOTAL ENERGY (no soft)'
    soft = 0.5_pf
    
    print *, '[ENE] original:   ', total_energy(ms, qs, ps, G, soft)
    E_desired = -0.25_pf
    CALL contotene(E_desired, ms, qs, ps, G, soft)
    E_final = total_energy(ms, qs, ps, G, soft)
    print *, '[ENE] conditioned:', E_final

    print *, '----------------------------------------------'

    ! measuring the errors
    com_error = NORM2(com_final)
    P_error = NORM2(P_final - P_desired)
    J_error = NORM2(J_final - J_desired)
    E_error = ABS(E_final - E_desired)

    print *, '[COM] error: ', com_error
    print *, '[TLM] error: ', P_error
    print *, '[TAM] error: ', J_error
    print *, '[ENE] error: ', E_error

    error = NORM2((/com_error, P_error, J_error, E_error/))

    IF (error < tol) THEN
        ier = 0
    ELSE
        ier = 1
    ENDIF
END SUBROUTINE
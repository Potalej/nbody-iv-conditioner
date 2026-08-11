!! # CONDITIONERS
!!
!! this module provides routines to condition initial values for the
!! gravitational N-body problem given desired first integrals and starting
!! in equilibrium if desired.
!!
!! all the routines are separated and public, but the four principal are:
!! - 

MODULE conditioners_mod
    USE utils_mod
    IMPLICIT NONE
    INTEGER, PARAMETER :: pf  = SELECTED_REAL_KIND(15, 307)
    PUBLIC
    PRIVATE pf
CONTAINS

!> @brief conditioner for the total angular momentum
!!
!! Given an desired J=(Jx,Jy,Jz), this subroutine apply a rotation over the
!! ps so the new total angular momentum is J.
!!
!! @author oap
!! @date 2026-08-04
!!
!! @param[in]        J    desired angular momentum, size: 3
!! @param[in]       ms    masses (N)
!! @param[in]       qs    positions (3,N)
!! @param[in,out]   ps    momenta (3,N)
!! @param[out]     ier    error flag. if ier=0, success. if ier=1, singular
!!                        inertia tensor.
!! @calledby
SUBROUTINE contotangmom (J, ms, qs, ps, ier)
! CONditioner for the TOTal ANGular MOMentum
    REAL(pf), INTENT(IN) :: J(3), ms(:), qs(:,:)
    REAL(pf), INTENT(INOUT) :: ps(:,:)
    INTEGER,  INTENT(OUT) :: ier

    REAL(pf) :: J_now(3), inertia_tensor(3,3)
    REAL(pf) :: rotation(3)
    INTEGER  :: p

    J_now = total_angular_momentum(qs, ps)
    inertia_tensor = general_inertia_tensor(ms, qs)

    ! solving the system I w = J - J_now
    rotation = solve_linear_system_3(inertia_tensor, J - J_now)
    
    DO p = 1, SIZE(ms)
        ps(:,p) = ps(:,p) + ms(p) * cross_product(qs(:,p), rotation)
    END DO

    ! success
    ier = 0
END SUBROUTINE

!> @brief conditioner for the total energy
!!
!! Conditioner for the total energy given an desired energy E.
!!
!! @author oap
!! @date 2026-08-04
!!
!! @param[in]        E    desired total energy
!! @param[in]       ms    masses (N)
!! @param[in,out]   qs    positions (3,N)
!! @param[in,out]   ps    momenta (3,N)
!! @param[in]        G    gravitational constant
!! @param[in]     soft    potential softening
!!
!! @calledby
SUBROUTINE contotene (E, ms, qs, ps, G, soft)
! CONditioner for the TOTal ENErgy
    REAL(pf), INTENT(IN)    :: E
    REAL(pf), INTENT(IN)    :: ms(:)
    REAL(pf), INTENT(INOUT) :: qs(:,:), ps(:,:)
    REAL(pf), INTENT(IN)    :: G, soft

    REAL(pf) :: kinect, potential
    REAL(pf) :: scale_factor

    kinect = kinect_energy(ms, ps)
    potential = potential_energy(ms, qs, G, soft)

    ! evaluates the scale factor to achieve the energy
    IF (E > 0) THEN
        scale_factor = SQRT((-potential + E)/kinect)
    ELSE
        scale_factor = SQRT(-potential/kinect)
    ENDIF

    ! apply over the momenta
    ps = scale_factor * ps

    ! if the desired energy isnt positive, we need to scale the positions too
    IF (E < 0) THEN
        ! if soft==0, the potential energy is homogeneous so we can apply a homothety
        IF (soft == 0) THEN
            scale_factor = 1.0_pf / (E/potential + 1.0_pf)
            qs = scale_factor * qs
        ! if soft!=0, we need to apply an iterative method to achieve the desired energy
        ELSE
            CALL conpotenesof(E, ms, qs, ps, G, soft)
        ENDIF
    ENDIF
END SUBROUTINE

!> @brief conditioner for the softened potential energy
!!
!! Conditioner for the softened potential energy using the Newton's method,
!! necessary due to the loss of homogeneity.
!! This subroutine does not guarantee (virial) equilibrium, it only provides
!! that V = E - T.
!!
!! @author oap
!! @date 2026-08-04
!!
!! @param[in]                  E    desired total energy
!! @param[in]                 ms    masses (N)
!! @param[in,out]             qs    positions (3,N)
!! @param[in,out]             ps    momenta (3,N)
!! @param[in]                  G    gravitational constant
!! @param[in]               soft    potential softening
!! @param[in] ini_guess[Optional]   initial guess for the Newton iterator.
!!
!! @calledby contotene
SUBROUTINE conpotenesof (E, ms, qs, ps, G, soft, ini_guess)
! CONditioner for the POTential ENErgy SOFtened
    REAL(pf), INTENT(IN)    :: E, ms(:)
    REAL(pf), INTENT(INOUT) :: qs(:,:)
    REAL(pf), INTENT(IN)    :: ps(:,:), G, soft
    REAL(pf), INTENT(IN), OPTIONAL :: ini_guess

    INTEGER  :: N
    REAL(pf) :: pe ! potential energy
    REAL(pf) :: ke ! kinect energy
    INTEGER  :: i, p1, p2
    REAL(pf) :: qab(3) ! qab = qb - qa
    REAL(pf), ALLOCATABLE :: dist2(:) ! mutual distances squared
    REAL(pf), ALLOCATABLE :: ms2(:)   ! product of every two different masses

    REAL(pf) :: newton_cnst
    REAL(pf) :: newton_error
    REAL(pf) :: newton_approx
    INTEGER  :: newton_counter ! counter of newton iterations
    REAL(pf) :: pot_factor, pot_der_factor ! potential energy and its derivative wrt factor
    REAL(pf) :: den ! denominator
    REAL(pf), ALLOCATABLE :: qs_factor(:,:)

    N = SIZE(ms)
    ALLOCATE(dist2(INT(N*(N-1)/2)))
    ALLOCATE(ms2(INT(N*(N-1)/2)))
    ALLOCATE(qs_factor(3,N))

    ! vectors to facilitate the evaluation of the potential and its derivative
    ! and the kinect energy
    i = 1
    pe = 0.0_pf
    ke = DOT_PRODUCT(ps(:,1),ps(:,1))/ms(1)
    DO p1 = 2, N
        ke = ke + DOT_PRODUCT(ps(:,p1),ps(:,p1))/ms(p1)
        DO p2 = 1, p1 - 1
            qab = qs(:,p1) - qs(:,p2)
            dist2(i) = DOT_PRODUCT(qab, qab)
            ms2(i) = ms(p1) * ms(p2)
            pe = pe - G * ms2(i) / SQRT(dist2(i) + soft*soft)
            i = i + 1
        END DO
    END DO
    ke = 0.5_pf * ke

    !! application of the Newton method
    newton_cnst = - soft * (E - ke)

    ! initial guess
    newton_error = 1.0_pf
    IF (PRESENT(ini_guess)) THEN
        newton_approx = ini_guess / soft
    ELSE
        newton_approx = 1.0_pf / (soft * (E/pe + 1.0_pf))
    ENDIF
    ! iterations
    newton_counter = 1
    DO WHILE (newton_error > 1E-15 .AND. newton_counter < 50)
        ! potential energy and its derivative wrt the scale factor
        pot_factor = 0.0_pf
        pot_der_factor = 0.0_pf

        DO i = 1, SIZE(dist2)
            den = SQRT(newton_approx*newton_approx*dist2(i) + 1.0_pf)
            pot_factor     = pot_factor - G * ms2(i)/den
            pot_der_factor = pot_der_factor + G * ms2(i)/(den**3)
        END DO
        pot_der_factor = pot_der_factor * newton_approx

        ! newton iteration
        newton_approx = newton_approx - (pot_factor + newton_cnst)/pot_der_factor

        ! error
        qs_factor = newton_approx * soft * qs
        newton_error = ABS(E - total_energy(ms, qs_factor, ps, G, soft))

        newton_counter = newton_counter + 1
    END DO

    ! update the positions
    qs = qs_factor

    ! free memory
    DEALLOCATE(dist2, ms2, qs_factor)
END SUBROUTINE

!> @brief conditioner for the softened potential energy (with equilibrium)
!!
!! Conditioner for the softened potential energy using the Newton's method,
!! necessary due to the loss of homogeneity.
!! This subroutine guarantee (virial) equilibrium.
!!
!! @author oap
!! @date 2026-08-04
!!
!! @param[in]                  E    desired total energy
!! @param[in]                 ms    masses (N)
!! @param[in,out]             qs    positions (3,N)
!! @param[in]                  G    gravitational constant
!! @param[in]               soft    potential softening
!! @param[in] ini_guess[Optional]   initial guess for the Newton iterator.
!!
!! @calledby
SUBROUTINE conpotenesofequ (E, ms, qs, G, soft, ini_guess)
! CONditioner for the POTential ENErgy SOFtened in EQUilibria
    REAL(pf), INTENT(IN) :: E, ms(:)
    REAL(pf), INTENT(INOUT) :: qs(:,:)
    REAL(pf), INTENT(IN) :: G, soft
    REAL(pf), INTENT(IN), OPTIONAL :: ini_guess

    INTEGER :: N
    REAL(pf), ALLOCATABLE :: forces(:,:), forces_der(:,:), qs_factor(:,:)
    REAL(pf) :: pot_factor, pot_der_factor
    REAL(pf) :: newton_error, newton_approx
    INTEGER  :: newton_counter
    INTEGER  :: p1, p2
    REAL(pf) :: qab(3), dist2, den, f, f_der, Fab(3)

    N = SIZE(ms)
    ALLOCATE(forces(3,N))
    ALLOCATE(forces_der(3,N))
    ALLOCATE(qs_factor(3,N))

    ! Newton error
    newton_error = 1.0_pf
    newton_approx = 1.0_pf
    IF (PRESENT(ini_guess)) newton_approx = ini_guess
    newton_counter = 0
    ! iterations
    DO WHILE (newton_error > 1E-15 .AND. newton_counter < 10)
        forces = 0.0_pf
        forces_der = 0.0_pf
        pot_factor = 0.0_pf
        pot_der_factor = 0.0_pf
        DO p1 = 2, N
            DO p2 = 1, p1 - 1
                qab = qs(:,p2) - qs(:,p1)
                dist2 = DOT_PRODUCT(qab, qab)
                den = SQRT(dist2 + (soft/newton_approx)**2)
                f = G * ms(p1) * ms(p2)

                ! forces
                Fab = f * qab / (den**3)
                forces(:,p1) = forces(:,p1) + Fab
                forces(:,p2) = forces(:,p2) - Fab

                ! forces derivative wrt the parameter
                Fab = Fab / (den * den)
                forces_der(:,p1) = forces_der(:,p1) + Fab
                forces_der(:,p2) = forces_der(:,p2) - Fab

                ! potential
                pot_factor = pot_factor - f / den

                ! potential derivative wrt the parameter
                pot_der_factor = pot_der_factor - f / (den**3)
            END DO
        END DO

        f = 2.0_pf * (E - pot_factor/newton_approx)
        f_der = 2.0_pf * pot_factor / (newton_approx**2)
        f_der = f_der - 2.0_pf * pot_der_factor * soft * soft/(newton_approx**4)
        DO p1 = 1, N
            ! forces
            f = f + DOT_PRODUCT(forces(:,p1), qs(:,p1)) / newton_approx

            ! forces derivative wrt the parameter
            f_der = f_der - DOT_PRODUCT(forces(:,p1), qs(:,p1)) / (newton_approx**2)
            f_der = f_der + DOT_PRODUCT(forces_der(:,p1), qs(:,p1)) * 3.0_pf * soft * soft/(newton_approx**4)
        END DO

        ! Newton method
        newton_approx = newton_approx - f/f_der

        ! Newton error
        qs_factor = newton_approx * qs
        forces = 0.0_pf
        pot_factor = 0.0_pf
        DO p1 = 2, N
            DO p2 = 1, p1 - 1
                qab = qs_factor(:,p2) - qs_factor(:,p1)
                dist2 = DOT_PRODUCT(qab, qab)
                den = SQRT(dist2 + soft**2)
                f = G * ms(p1) * ms(p2)
                
                Fab = f * qab / (den**3)
                forces(:,p1) = forces(:,p1) + Fab
                forces(:,p2) = forces(:,p2) - Fab

                pot_factor = pot_factor - f / den
            END DO
        END DO
        newton_error = 2.0_pf * (E - pot_factor)
        DO p1 = 1, N
            newton_error = newton_error + DOT_PRODUCT(forces(:,p1), qs_factor(:,p1))
        END DO
        newton_error = ABS(newton_error) / (3.0_pf * N)
        
        newton_counter = newton_counter + 1
    END DO

    ! update the positions
    qs = qs_factor

    ! free memory
    DEALLOCATE(forces, forces_der, qs_factor)
END SUBROUTINE

!> @brief center of masses to origin
!!
!! This subroutine moves the center of masses of the system to the origin.
!!
!! @author oap
!! @date 2026-08-04
!!
!! @param[in]                 ms    masses (N)
!! @param[in,out]             qs    positions (3,N)
!!
!! @calledby
SUBROUTINE comori (ms, qs)
! COM to ORIgin
    REAL(pf), INTENT(IN) :: ms(:)
    REAL(pf), INTENT(INOUT) :: qs(:,:)
    REAL(pf) :: com(3)
    INTEGER  :: p

    com = center_of_mass(ms, qs)
    DO p = 1, SIZE(ms)
        qs(:,p) = qs(:,p) - com
    END DO
END SUBROUTINE

!> @brief conditioner for the total linear momentum
!!
!! Given an desired P=(Px,Py,Pz), this subroutine changes the
!! velocities proportionally to some weights vector
!!
!! @author oap
!! @date 2026-08-04
!!
!! @param[in]                  P    desired total linear momentum (3)
!! @param[in,out]             ps    linear momenta (3,N)
!! @param[in]                 ws    weight vector (N)
!!
!! @calledby
SUBROUTINE contotlinmom (P, ps, ws)
! CONditioner for the TOTal LINear MOMentum
    REAL(pf), INTENT(IN) :: P(3)
    REAL(pf), INTENT(INOUT) :: ps(:,:)
    REAL(pf), INTENT(IN) :: ws(:)
    REAL(pf) :: pcm(3)
    INTEGER  :: a

    pcm = (total_linear_momentum(ps) - P)/SUM(ws)
    
    DO a = 1, SIZE(ws)
        ps(:,a) = ps(:,a) - ws(a) * pcm
    END DO
END SUBROUTINE

!=====================================================================!
!=====================================================================!
!> SUBROUTINES WITH CONDITIONERS

!> @brief conditioner for initial values (iterative)
!!
!! Given desireds first integrals total linear momentum P, total angular
!! momentum J and total energy E, this subroutine apply the following
!! process:
!! 1. Moves the COM to the origin.
!! 2. Applies the P conditioner
!! 3. Applies the J conditioner
!! 4. Applies the E conditioner
!! 5. IF ($||(P, J, E) - (P_des, J_des, E_des)||_\infty > error_limit$)
!!    Go back to 2.
!!
!! @author oap
!! @date 2026-08-05
!!
!! @param[in]                    ms    masses (N)
!! @param[in,out]                qs    positions (3,N)
!! @param[in,out]                ps    linear momentum (3,N)
!! @param[in]                     G    gravitational constant
!! @param[in]                  soft    potential softening
!! @param[in]                   ier    error flag. if ier=0, success. else, error
!! @param[in]                     E    desired total energy
!! @param[in]        J_par[Optional]   desired total angular momentum (3)
!! @param[in]        P_par[Optional]   desired total linear momentum (3)
!! @param[in] nitermax_par[Optional]   max number of iterations
!!
!! @calledby
SUBROUTINE coninivalite (ms, qs, ps, G, soft, ier, E, J_par, P_par, nitermax_par)
! CONditioner for INItial VALues ITErative
    REAL(pf), INTENT(IN) :: ms(:)
    REAL(pf), INTENT(INOUT) :: qs(:,:), ps(:,:)
    REAL(pf), INTENT(IN) :: G, soft, E
    INTEGER,  INTENT(INOUT) :: ier
    REAL(pf), INTENT(IN), OPTIONAL :: J_par(3), P_par(3)
    INTEGER,  INTENT(IN), OPTIONAL :: nitermax_par
    
    REAL(pf) :: error_limit = 1E-8

    REAL(pf) :: J(3), P(3)
    INTEGER  :: nitermax
    REAL(pf) :: error_E, error_P, error_J, error_all
    INTEGER  :: step

    nitermax = 10
    IF (PRESENT(nitermax_par)) nitermax = nitermax_par

    J = 0.0_pf
    P = 0.0_pf
    IF (PRESENT(J_par)) J = J_par
    IF (PRESENT(P_par)) P = P_par

    ! move the com to origin
    CALL comori(ms, qs)

    ! conditioners
    CALL contotlinmom(P, ps, MS)
    CALL contotangmom(J, ms, qs, ps, ier)
    IF (ier > 0) RETURN
    CALL contotene(E, ms, qs, ps, G, soft)

    ! errors
    error_E = ABS(total_energy(ms, qs, ps, G, soft) - E)
    error_P = MAXVAL(ABS(total_linear_momentum(ps) - P))
    error_J = MAXVAL(ABS(total_angular_momentum(qs, ps) - J))
    error_all = MAXVAL((/error_E, error_P, error_J/))

    ! starts the iteration if needs
    step = 1
    IF (error_all >= error_limit) THEN
        DO WHILE (error_all >= error_limit .AND. step < nitermax)
            step = step + 1

            ! conditioning
            CALL contotlinmom(P, ps, ms)
            CALL contotangmom(J, ms, qs, ps, ier)
            IF (ier > 0) RETURN
            CALL contotene(E, ms, qs, ps, G, soft)

            ! errors
            error_E = ABS(total_energy(ms, qs, ps, G, soft) - E)
            error_P = MAXVAL(ABS(total_linear_momentum(ps) - P))
            error_J = MAXVAL(ABS(total_angular_momentum(qs, ps) - J))
            error_all = MAXVAL((/error_E, error_P, error_J/))
        END DO
    END IF

    ! success
    ier = 0
END SUBROUTINE

!> @brief conditioner for initial values (direct)
!!
!! Given desireds first integrals total linear momentum P, total angular
!! momentum J and total energy E, this subroutine condition the given
!! state vectors directly.
!!
!! IMPORTANT: this subroutine only works for the classic potential.
!!
!! @author oap
!! @date 2026-08-10
!!
!! @param[in]                 ms    masses vector (N)
!! @param[in,out]             ps    linear momenta (3,N)
!! @param[in,out]             qs    positions (3,N)
!! @param[in]                  G    gravitational constant
!! @param[out]               ier    error flag. ier=0 implies successful,
!!                                  ier=1/2 implies singular inertia tensor
!!                                  and ier=4 implies a zero coefficient
!! @param[in]                  E    desired energy
!! @param[in]              J_par    desired angular momentum (3)
!! @param[in]              P_par    desired linear momentum (3)
!! @calledby
SUBROUTINE coninivaldir (ms, qs, ps, G, ier, E, J_par, P_par)
! CONditioner for INItial VALues DIRect
    REAL(pf), INTENT(IN) :: ms(:)
    REAL(pf), INTENT(INOUT) :: qs(:,:), ps(:,:)
    REAL(pf), INTENT(IN) :: G, E
    REAL(pf), INTENT(IN), OPTIONAL :: J_par(3), P_par(3)
    INTEGER,  INTENT(OUT) :: ier

    REAL(pf) :: J(3), P(3)
    REAL(pf) :: energy, linear(3), angular(3), Minv, pot
    REAL(pf) :: inertia_tensor(3,3) ! W_T
    REAL(pf) :: rot(3), rot_(3)     ! omega, omega_tilde
    REAL(pf) :: alpha ! positions homothety factor
    REAL(pf) :: sigma_tilde
    REAL(pf) :: S1, S2, beta ! velocities factors
    INTEGER  :: a

    J = 0.0_pf
    P = 0.0_pf
    IF (PRESENT(J_par)) J = J_par
    IF (PRESENT(P_par)) P = P_par

    ! move the com to origin
    CALL comori(ms, qs)

    ! evaluate first integrals and other values
    energy  = total_energy(ms, qs, ps, G, 0.0_pf)
    linear  = total_linear_momentum(ps)
    angular = total_angular_momentum(qs, ps)
    Minv    = 1.0_pf / SUM(ms)
    pot     = energy - kinect_energy(ms, ps)

    ! to evaluate the beta factor
    inertia_tensor = general_inertia_tensor(ms, qs)
    rot  = solve_linear_system_3(inertia_tensor, angular)
    rot_ = solve_linear_system_3(inertia_tensor, J)

    ! for nonnegative energy, we just need to change the velocities
    ! for negative energy, its necessary an increment to beta exists
    alpha = 1.0_pf
    IF (E < 0) alpha = alpha + E / pot

    ! if is there angular (but not the linear) is necessary to add the
    ! rotation element
    IF (NORM2(J) > 0) THEN
        sigma_tilde = - DOT_PRODUCT(J, rot_)
        alpha = - pot / sigma_tilde

    ! if is there no angular momentum but have linear, is necessary to
    ! add the linear contribution to beta exists
    ELSE IF (NORM2(P) > 0) THEN
        alpha = alpha - 0.5_pf * Minv * NORM2(P)**2 / pot
    ENDIF

    ! evaluating beta
    S1 = (energy - pot) - 0.5_pf * Minv * NORM2(linear)**2 + 0.5_pf * DOT_PRODUCT(angular, rot)
    S2 = (0.5_pf * Minv * NORM2(P)**2 - 0.5_pf * alpha * alpha * DOT_PRODUCT(J, rot_))
    beta = SQRT((E - alpha * pot - S2)/S1)

    IF (alpha == 0.0_pf .OR. beta == 0.0_pf) THEN
        ier = 4
        RETURN
    ENDIF

    ! saving
    rot = solve_linear_system_3(inertia_tensor, angular - J*alpha/beta)

    DO a = 1, SIZE(ms)
        ps(:,a) = ps(:,a) - ms(a) * Minv * (linear - P / beta)
        ps(:,a) = ps(:,a) - ms(a) * cross_product(qs(:,a), rot)
        ps(:,a) = beta * ps(:,a)
    END DO
    qs = qs / alpha
    
    ier = 0
END SUBROUTINE

!> @brief conditioner for initial values (Aarseth conditioner)
!!
!! This conditioner provides the Henon units via Aarseths method. It
!! conditionates the c.o.m. to the origin, the momenta integrals to 
!! zero, the total energy to -1/4 and satisfies the virial ratios be
!! the unit.
!!
!! IMPORTANT: this subroutine only works for the classic potential.
!!            for softened potentials, use the modified version.
!!
!! @author oap
!! @date 2026-08-10
!!
!! @param[in]                 ms    masses vector (N)
!! @param[in,out]             ps    linear momenta (3,N)
!! @param[in,out]             qs    positions (3,N)
!! @param[in]                  G    gravitational constant
!! @param[out]               ier    error flag
!!
!! @calledby
SUBROUTINE coninivalaar (ms, qs, ps, G, ier)
! CONditioner for INItial VALues AARseth
    REAL(pf), INTENT(INOUT) :: ms(:), qs(:,:), ps(:,:)
    REAL(pf), INTENT(IN) :: G
    INTEGER,  INTENT(INOUT) :: ier
    
    REAL(pf) :: pot, kin ! potential and kinect energies
    REAL(pf) :: Qv   ! virial radius
    REAL(pf) :: beta ! scale factor
    
    ! normalize the masses to 1/N
    ms = 1.0_pf / SIZE(ms)

    ! conditionate the first integrals to zero
    CALL coninivaldir(ms, qs, ps, G, ier, 0.0_pf)
    IF (ier .NE. 0) THEN
        RETURN
    ENDIF

    ! Aarseth method
    pot  = potential_energy(ms, qs, G, 0.0_pf)
    kin  = kinect_energy(ms, ps)
    Qv   = SQRT(0.5_pf * ABS(pot) / kin)
    beta = 0.5_pf * pot / (-0.25_pf)

    ! conditionate
    ps = ps * Qv / SQRT(beta)
    qs = qs * beta
END SUBROUTINE

!> @brief conditioner for initial values (modified Aarseth conditioner)
!!
!! This conditioner provides the Henon units via a modified Aarseths 
!! method. It conditionates the c.o.m. to the origin, the momenta integrals 
!! to zero, the total energy to -1/4 and satisfies the virial ratios be
!! the unit.
!! If the softening parameter is zero, this subroutine provides the same
!! result of coninivalaar.
!!
!! @author oap
!! @date 2026-08-11
!!
!! @param[in]                 ms    masses vector (N)
!! @param[in,out]             ps    linear momenta (3,N)
!! @param[in,out]             qs    positions (3,N)
!! @param[in]                  G    gravitational constant
!! @param[in]               soft    softening parameter
!! @param[out]               ier    error flag
!!
!! @calledby
SUBROUTINE coninivalaarmod (ms, qs, ps, G, soft, ier)
! CONditioner for INItial VALues AARseth MODified
    REAL(pf), INTENT(INOUT) :: ms(:), qs(:,:), ps(:,:)
    REAL(pf), INTENT(IN)  :: G, soft
    INTEGER, INTENT(INOUT) :: ier
    
    REAL(pf) :: pot, kin ! potential and kinect energies
    REAL(pf) :: Qv   ! virial radius
    REAL(pf) :: beta ! scale factor
    
    ! normalize the masses to 1/N
    ms = 1.0_pf / SIZE(ms)

    ! conditionate the first integrals to zero
    IF (soft == 0.0_pf) THEN
        CALL coninivaldir(ms, qs, ps, G, ier, 0.0_pf)
    ELSE
        CALL coninivalite(ms, qs, ps, G, soft, ier, 0.0_pf)
    ENDIF

    ! Aarseth method
    pot  = potential_energy(ms, qs, G, soft)
    kin  = kinect_energy(ms, ps)
    Qv   = SQRT(0.5_pf * ABS(pot) / kin)
    beta = 0.5_pf * pot / (-0.25_pf)

    ! conditionate the velocities
    ps = ps * Qv / SQRT(beta)

    ! if no softening, just apply the same as the original method
    IF (soft == 0.0_pf) THEN
        qs = qs * beta

    ! if using softening, it needs to be iterative and uses beta/soft as initial guess
    ELSE
        CALL conpotenesofequ(-0.25_pf, ms, qs, G, soft, beta)
        pot = potential_energy(ms, qs, G, soft)

        ! verify if its possible to reconditionate
        IF (4.0_pf * ABS(pot) < 1.0_pf) THEN
            ier = 4
            RETURN
        ENDIF

        ! reconditionate the linear momentum
        beta = SQRT(4.0_pf * ABS(pot) - 1.0_pf)
        ps = ps * beta
    ENDIF
    
    ! success
    ier = 0
END SUBROUTINE

END MODULE
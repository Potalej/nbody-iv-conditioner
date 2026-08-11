MODULE utils_mod
    IMPLICIT NONE
    PUBLIC
    PRIVATE pf
    INTEGER, PARAMETER  :: pf  = SELECTED_REAL_KIND(15, 307)
CONTAINS

FUNCTION center_of_mass (ms, qs) RESULT(com)
    REAL(pf), INTENT(IN) :: ms(:), qs(:,:)
    REAL(pf) :: com(3), mtot
    INTEGER  :: p

    com = 0.0_pf
    mtot = 0.0_pf
    DO p = 1, SIZE(ms)
        com = com + ms(p) * qs(:,p)
        mtot = mtot + ms(p)
    END DO
    com = com / mtot
END FUNCTION

FUNCTION total_linear_momentum (ps) RESULT(Ptot)
    REAL(pf), INTENT(IN) :: ps(:,:)
    REAL(pf) :: Ptot(3)
    INTEGER  :: p

    Ptot(1) = SUM(ps(1,:))
    Ptot(2) = SUM(ps(2,:))
    Ptot(3) = SUM(ps(3,:))
END FUNCTION

FUNCTION kinect_energy (ms, ps) RESULT (T)
    REAL(pf), INTENT(IN) :: ms(:), ps(:,:)
    REAL(pf) :: T
    INTEGER  :: p

    T = 0.0_pf
    DO p = 1, SIZE(ms)
        T = T + (ps(1,p)*ps(1,p) + ps(2,p)*ps(2,p) + ps(3,p)*ps(3,p))/ms(p)
    END DO
    T = 0.5 * T
END FUNCTION

FUNCTION potential_energy (ms, qs, G, soft) RESULT (V)
    REAL(pf), INTENT(IN) :: ms(:), qs(:,:), G, soft
    
    REAL(pf) :: V
    REAL(pf) :: m1, x1, y1, z1
    REAL(pf) :: dx, dy, dz
    REAL(pf) :: invdist
    INTEGER  :: p1, p2

    V = 0.0_pf
    DO p1 = 2, SIZE(ms)
        m1 = ms(p1)
        x1 = qs(1,p1)
        y1 = qs(2,p1)
        z1 = qs(3,p1)
        DO p2 = 1, p1 - 1
            dx = qs(1,p2) - x1
            dy = qs(2,p2) - y1
            dz = qs(3,p2) - z1
            
            invdist = SQRT(dx*dx + dy*dy + dz*dz + soft*soft)
            invdist = 1.0_pf / invdist

            V = V - m1 * ms(p2) * invdist
        END DO
    END DO
    V = G * V
END FUNCTION

FUNCTION total_energy (ms, qs, ps, G, soft) RESULT(E)
    REAL(pf), INTENT(IN) :: ms(:), qs(:,:), ps(:,:), G, soft
    
    REAL(pf) :: E
    REAL(pf) :: m1, x1, y1, z1
    REAL(pf) :: dx, dy, dz
    REAL(pf) :: invdist
    INTEGER  :: p1, p2

    E = 0.5_pf * DOT_PRODUCT(ps(:,1), ps(:,1)) / ms(1)

    DO p1 = 2, SIZE(ms)
        E = E + 0.5_pf * DOT_PRODUCT(ps(:,p1), ps(:,p1)) / ms(p1)

        m1 = ms(p1)
        x1 = qs(1,p1)
        y1 = qs(2,p1)
        z1 = qs(3,p1)
        DO p2 = 1, p1 - 1
            dx = qs(1,p2) - x1
            dy = qs(2,p2) - y1
            dz = qs(3,p2) - z1
            
            invdist = SQRT(dx*dx + dy*dy + dz*dz + soft*soft)
            invdist = 1.0_pf / invdist

            E = E - G * m1 * ms(p2) * invdist
        END DO
    END DO
END FUNCTION

FUNCTION total_angular_momentum (qs, ps) RESULT(J)
! total angular momentum
    REAL(pf), INTENT(IN) :: qs(:,:), ps(:,:)
    
    REAL(pf) :: J(3)
    INTEGER  :: p

    J = 0.0_pf
    DO p = 1, SIZE(qs,2)
        J = J + cross_product(qs(:,p), ps(:,p))
    END DO
END FUNCTION

FUNCTION individual_inertia_tensor (m, q) RESULT(I)
! general inertia tensor
    REAL(pf), INTENT(IN) :: m, q(3)
    
    REAL(pf) :: I(3,3)
    INTEGER  :: a, b

    DO a = 2, 3
        DO b = 1, a-1
            I(a,b) = m * q(a) * q(b)
            I(b,a) = I(a,b)
        END DO
    END DO

    I(1,1) = - m * (q(2)**2 + q(3)**2)
    I(2,2) = - m * (q(1)**2 + q(3)**2)
    I(3,3) = - m * (q(1)**2 + q(2)**2)
END FUNCTION

FUNCTION general_inertia_tensor (ms, qs) RESULT(I)
! general inertia tensor
    REAL(pf), INTENT(IN) :: ms(:), qs(:,:)

    REAL(pf) :: I(3,3)
    INTEGER  :: a

    I = 0.0_pf
    DO a = 1, SIZE(ms)
        I = I + individual_inertia_tensor(ms(a), qs(:,a))
    END DO
END FUNCTION

FUNCTION solve_linear_system_3 (A_orig, b) RESULT(sol)
! solve linear system in R3
    REAL(pf), INTENT(IN)  :: A_orig(3,3), b(3)
    REAL(pf) :: tmp(4), sol(3), A(3,4)
    INTEGER  :: pivot

    ! matrix 3 x 4 [A | b]
    A = RESHAPE((/ A_orig, b /), [3,4])

    !> first line
    pivot = 1
    IF (ABS(A(pivot, 1)) < ABS(A(2, 1))) pivot = 2
    IF (ABS(A(pivot, 1)) < ABS(A(3, 1))) pivot = 3

    ! if pivot != 1
    IF (pivot .NE. 1) THEN
        tmp       = A(pivot,:)
        A(pivot,:) = A(1,:)
        A(1,:)    = tmp
    ENDIF

    ! elimination
    tmp = A(1,:) / A(1,1)
    A(3,:) = A(3,:) - tmp * A(3,1)
    A(2,:) = A(2,:) - tmp * A(2,1)
    A(1,:) = tmp

    !> second line
    ! pivot if needed
    IF (ABS(A(2, 2)) < ABS(A(3, 2))) THEN
    tmp = A(3,:)
    A(3,:) = A(2,:)
    A(2,:) = tmp
    ENDIF

    ! elimination
    A(3,:) = A(3,:) - A(2,:) * A(3,2) / A(2,2)
    A(2,:) = A(2,:) / A(2,2)

    !> third line
    A(3,:) = A(3,:) / A(3,3)

    !> solving
    sol(3) = A(3,4)
    sol(2) = A(2,4) - A(2,3) * sol(3)
    sol(1) = A(1,4) - A(1,3) * sol(3) - A(1,2) * sol(2)
END FUNCTION

FUNCTION cross_product (u, v)
! cross product
  REAL(pf), INTENT(IN) :: u(3), v(3)
  REAL(pf) :: cross_product(3)

  cross_product(1) =  u(2)*v(3)-v(2)*u(3)
  cross_product(2) = -u(1)*v(3)+v(1)*u(3)
  cross_product(3) =  u(1)*v(2)-v(1)*u(2)
END FUNCTION

END MODULE
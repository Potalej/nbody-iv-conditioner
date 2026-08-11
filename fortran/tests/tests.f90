PROGRAM tests
    IMPLICIT NONE
    INTEGER, PARAMETER  :: pf  = SELECTED_REAL_KIND(15, 307)
    
    INTEGER :: ier(3)
    INTEGER :: t

    ier = 0

! testing
    CALL test1(ier(1)) ! separated methods
    CALL test2(ier(2)) ! iterative (both classical and softened)
    CALL test3(ier(3)) ! direct

! now organizing the results
    print *, ''
    print *, '----------------------------------------------'
    print *, '> RESULTS'
    print *, 'ier =0: success'
    print *, 'ier!=0: error'
    print *, '----------------------------------------------'

    print *, 'test number |         ier'
    DO t = 1, SIZE(ier)
        PRINT *, t, '|', ier(t)
    END DO
END PROGRAM
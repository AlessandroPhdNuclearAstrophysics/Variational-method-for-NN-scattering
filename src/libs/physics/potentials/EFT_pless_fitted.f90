SUBROUTINE EFT_PLESS_PW_FITTED(ILB, R, L, S, J, TZ, V)
  IMPLICIT NONE
  INTEGER, PARAMETER :: NFIT = 4
  DOUBLE PRECISION, PARAMETER :: HTC = 197.32697D0
  DOUBLE PRECISION, PARAMETER :: PI= 4.D0*DATAN(1.D0)
  DOUBLE PRECISION, PARAMETER :: C2001(NFIT) = &
                    [                         &
                      -2.213333D0,            & ! Alessandro, grid, AV18 1 MeV CM
                      -2.144536513263479D0,   & ! Alessandro, PETSc, AV18 1 MeV CM
                      -1.11708432681862D1,    & ! Ylenia, MINUIT, Granada
                      -1.153170815331322D1    & ! Ylenia, MINUIT, AV18 1 MeV CM
                    ]
  DOUBLE PRECISION, PARAMETER :: RV00(NFIT) = &
                    [                        &
                      0.7700D0,              & ! Alessandro, grid, AV18 1 MeV CM
                      0.761781945701797D0,   & ! Alessandro, PETSc, AV18 1 MeV CM
                      0.23724000000000D1,    & ! Ylenia, MINUIT, Granada
                      0.27603D1              & ! Ylenia, MINUIT, AV18 1 MeV CM
                    ]

  INTEGER, INTENT(IN) :: L, S, J, TZ, ILB
  DOUBLE PRECISION, INTENT(IN) :: R
  DOUBLE PRECISION, INTENT(OUT):: V(2,2)

  INTEGER :: T

  T = MOD(MOD(L+S, 2) + 1, 2)
  IF (ABS(TZ) > T) THEN
    WRITE(*,*) "EFT_PLESS_PW_FITTED: Invalid TZ value"
    RETURN
  ENDIF

  IF (ABS(L-S)> J .OR. L+S < J) THEN
    WRITE(*,*) "EFT_PLESS_PW_FITTED: Invalid L, S, J combination"
    RETURN
  ENDIF

  V = 0.D0
  IF ( T.EQ.0 .AND. S.EQ.0 ) THEN
    V(1,1) = C2001(ILB)*(-4*R**2 + 6*RV00(ILB)**2)/(DEXP(R**2/RV00(ILB)**2)*PI**1.5*RV00(ILB)**7)
  ENDIF
  V = V*HTC

  RETURN
END SUBROUTINE EFT_PLESS_PW_FITTED
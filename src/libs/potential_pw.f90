!> \file potential_pw.f90
!! \brief Interfaccia per il calcolo dei potenziali nucleari a onde parziali.
!!
!! Questo modulo fornisce la subroutine POT_PW per il calcolo della matrice
!! del potenziale nucleone-nucleone in diversi modelli (AV14, AV18, EFT_PLESS, ...).
!!
!! \author Alessandro Grassi
!! \date 2025
MODULE POTENTIALS
  IMPLICIT NONE
  !> Espone la subroutine per il calcolo del potenziale a onde parziali.
  PUBLIC :: POT_PW
CONTAINS

  !> Calcola la matrice del potenziale nucleare a onde parziali per diversi modelli.
  !>
  !> In base al valore di IPOT seleziona il modello di potenziale desiderato e
  !> richiama la relativa subroutine di calcolo.
  !>
  !> \param[in] IPOT Identificatore del potenziale (14=AV14, 18=AV18, 19=EFT_PLESS, 21=EFT_PLESS_FITTED)
  !> \param[in] ILB Indice del canale (usato da alcuni potenziali)
  !> \param[in] LEMP Tipo di interazione elettromagnetica (0=nessuna, >0=presente)
  !> \param[in] L Momento angolare orbitale
  !> \param[in] S Momento di spin totale
  !> \param[in] J Momento angolare totale
  !> \param[in] T1Z Proiezione isospin nucleone 1 (+1=protone, -1=neutrone)
  !> \param[in] T2Z Proiezione isospin nucleone 2 (+1=protone, -1=neutrone)
  !> \param[in] R Distanza interparticellare (fm)
  !> \param[out] VPW Matrice 2x2 del potenziale in base accoppiata (canali misti)
  SUBROUTINE POT_PW(IPOT, ILB, LEMP, L, S, J, T1Z, T2Z, R, VPW)
    USE AV18
    USE AV14
    USE EFT_PLESS
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: IPOT, LEMP, L, S, J, T1Z, T2Z
    DOUBLE PRECISION, INTENT(IN) :: R
    DOUBLE PRECISION, INTENT(OUT):: VPW(2,2)
    INTEGER, INTENT(IN) :: ILB
  !
    INTEGER :: T
  !
    T = MOD( MOD(L+S,2) + 1 , 2)
  !
    SELECT CASE (IPOT)
    CASE (14)
      CALL AV14PW(LEMP, L, S, J, T1Z, T2Z, R, VPW)
      RETURN
    CASE (18)
      CALL AV18PW(ILB, L, S, J, T, T1Z, T2Z, R, VPW, LEMP)
      RETURN
    CASE (19)
      ! CALL EFT_PLESS_PW(LEMP, ILB, L, S, J, T1Z, T2Z, R, VPW)
      STOP "POT_PW: EFT_PLESS_PW not implemented"
      RETURN
    CASE (21)
      CALL EFT_PLESS_PW_FITTED(ILB, R, L, S, J, T1Z + T2Z, VPW)
    CASE DEFAULT
      STOP "ERROR, POTENTIAL NOT FOUND"
    END SELECT
    RETURN

  END SUBROUTINE POT_PW
END MODULE POTENTIALS
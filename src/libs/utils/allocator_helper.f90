!> \file allocator_helper.f90
!! \brief Utilities for (re)allocating allocatable arrays of various ranks and types.
!!
!! This module provides a set of helper subroutines to allocate or reallocate
!! 1D, 2D, 3D, and 4D arrays of DOUBLE PRECISION and INTEGER types.
!! If the array is not allocated or its size/shape does not match the requested one,
!! it will be (re)allocated accordingly.
!!
!! \author Alessandro
!! \date 2025
MODULE REALLOCATE_UTILS
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: REALLOCATE_1D_1_INT
  PUBLIC :: REALLOCATE_1D_1, REALLOCATE_1D_2, REALLOCATE_1D_3, REALLOCATE_1D_4, REALLOCATE_1D_5
  PUBLIC :: REALLOCATE_1D_6, REALLOCATE_1D_7, REALLOCATE_1D_8, REALLOCATE_1D_9, REALLOCATE_1D_10
  PUBLIC :: REALLOCATE_2D_1, REALLOCATE_2D_2, REALLOCATE_2D_3, REALLOCATE_2D_4, REALLOCATE_2D_5
  PUBLIC :: REALLOCATE_2D_6, REALLOCATE_2D_7, REALLOCATE_2D_8, REALLOCATE_2D_9, REALLOCATE_2D_10
  PUBLIC :: REALLOCATE_3D_1
  PUBLIC :: REALLOCATE_4D_1

CONTAINS

! -------- 1D --------
!> \brief Allocate or reallocate a 1D DOUBLE PRECISION array.
!! \param[inout] A1 Array to (re)allocate
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_1(A1, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:)
  INTEGER, INTENT(IN) :: N
  IF (.NOT. ALLOCATED(A1)) THEN
    ALLOCATE(A1(N))
  ELSEIF (SIZE(A1) .NE. N) THEN
    DEALLOCATE(A1)
    ALLOCATE(A1(N))
  END IF
END SUBROUTINE

!> \brief Allocate or reallocate a 1D INTEGER array.
!! \param[inout] A1 Array to (re)allocate
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_1_INT(A1, N)
  INTEGER, ALLOCATABLE, INTENT(INOUT) :: A1(:)
  INTEGER, INTENT(IN) :: N
  IF (.NOT. ALLOCATED(A1)) THEN
    ALLOCATE(A1(N))
  ELSEIF (SIZE(A1) .NE. N) THEN
    DEALLOCATE(A1)
    ALLOCATE(A1(N))
  END IF
END SUBROUTINE

!> \brief Allocate or reallocate two 1D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array to (re)allocate
!! \param[inout] A2 Second array to (re)allocate
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_2(A1, A2, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:), A2(:)
  INTEGER, INTENT(IN) :: N
  CALL REALLOCATE_1D_1(A1, N)
  CALL REALLOCATE_1D_1(A2, N)
END SUBROUTINE

!> \brief Allocate or reallocate three 1D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_3(A1, A2, A3, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:), A2(:), A3(:)
  INTEGER, INTENT(IN) :: N
  CALL REALLOCATE_1D_1(A1, N)
  CALL REALLOCATE_1D_1(A2, N)
  CALL REALLOCATE_1D_1(A3, N)
END SUBROUTINE

!> \brief Allocate or reallocate four 1D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_4(A1, A2, A3, A4, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:), A2(:), A3(:), A4(:)
  INTEGER, INTENT(IN) :: N
  CALL REALLOCATE_1D_1(A1, N)
  CALL REALLOCATE_1D_1(A2, N)
  CALL REALLOCATE_1D_1(A3, N)
  CALL REALLOCATE_1D_1(A4, N)
END SUBROUTINE

!> \brief Allocate or reallocate five 1D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_5(A1, A2, A3, A4, A5, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:), A2(:), A3(:), A4(:), A5(:)
  INTEGER, INTENT(IN) :: N
  CALL REALLOCATE_1D_1(A1, N)
  CALL REALLOCATE_1D_1(A2, N)
  CALL REALLOCATE_1D_1(A3, N)
  CALL REALLOCATE_1D_1(A4, N)
  CALL REALLOCATE_1D_1(A5, N)
END SUBROUTINE

!> \brief Allocate or reallocate six 1D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_6(A1, A2, A3, A4, A5, A6, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:), A2(:), A3(:), A4(:), A5(:), A6(:)
  INTEGER, INTENT(IN) :: N
  CALL REALLOCATE_1D_1(A1, N)
  CALL REALLOCATE_1D_1(A2, N)
  CALL REALLOCATE_1D_1(A3, N)
  CALL REALLOCATE_1D_1(A4, N)
  CALL REALLOCATE_1D_1(A5, N)
  CALL REALLOCATE_1D_1(A6, N)
END SUBROUTINE

!> \brief Allocate or reallocate seven 1D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[inout] A7 Seventh array
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_7(A1, A2, A3, A4, A5, A6, A7, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:), A2(:), A3(:), A4(:), A5(:), A6(:), A7(:)
  INTEGER, INTENT(IN) :: N
  CALL REALLOCATE_1D_1(A1, N)
  CALL REALLOCATE_1D_1(A2, N)
  CALL REALLOCATE_1D_1(A3, N)
  CALL REALLOCATE_1D_1(A4, N)
  CALL REALLOCATE_1D_1(A5, N)
  CALL REALLOCATE_1D_1(A6, N)
  CALL REALLOCATE_1D_1(A7, N)
END SUBROUTINE

!> \brief Allocate or reallocate eight 1D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[inout] A7 Seventh array
!! \param[inout] A8 Eighth array
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_8(A1, A2, A3, A4, A5, A6, A7, A8, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:), A2(:), A3(:), A4(:), A5(:), A6(:), A7(:), A8(:)
  INTEGER, INTENT(IN) :: N
  CALL REALLOCATE_1D_1(A1, N)
  CALL REALLOCATE_1D_1(A2, N)
  CALL REALLOCATE_1D_1(A3, N)
  CALL REALLOCATE_1D_1(A4, N)
  CALL REALLOCATE_1D_1(A5, N)
  CALL REALLOCATE_1D_1(A6, N)
  CALL REALLOCATE_1D_1(A7, N)
  CALL REALLOCATE_1D_1(A8, N)
END SUBROUTINE

!> \brief Allocate or reallocate nine 1D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[inout] A7 Seventh array
!! \param[inout] A8 Eighth array
!! \param[inout] A9 Ninth array
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_9(A1, A2, A3, A4, A5, A6, A7, A8, A9, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:), A2(:), A3(:), A4(:), A5(:), A6(:), A7(:), A8(:), A9(:)
  INTEGER, INTENT(IN) :: N
  CALL REALLOCATE_1D_1(A1, N)
  CALL REALLOCATE_1D_1(A2, N)
  CALL REALLOCATE_1D_1(A3, N)
  CALL REALLOCATE_1D_1(A4, N)
  CALL REALLOCATE_1D_1(A5, N)
  CALL REALLOCATE_1D_1(A6, N)
  CALL REALLOCATE_1D_1(A7, N)
  CALL REALLOCATE_1D_1(A8, N)
  CALL REALLOCATE_1D_1(A9, N)
END SUBROUTINE

!> \brief Allocate or reallocate ten 1D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[inout] A7 Seventh array
!! \param[inout] A8 Eighth array
!! \param[inout] A9 Ninth array
!! \param[inout] A10 Tenth array
!! \param[in] N Desired size
SUBROUTINE REALLOCATE_1D_10(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, N)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:), A2(:), A3(:), A4(:), A5(:), A6(:), A7(:), A8(:), A9(:), A10(:)
  INTEGER, INTENT(IN) :: N
  CALL REALLOCATE_1D_1(A1, N)
  CALL REALLOCATE_1D_1(A2, N)
  CALL REALLOCATE_1D_1(A3, N)
  CALL REALLOCATE_1D_1(A4, N)
  CALL REALLOCATE_1D_1(A5, N)
  CALL REALLOCATE_1D_1(A6, N)
  CALL REALLOCATE_1D_1(A7, N)
  CALL REALLOCATE_1D_1(A8, N)
  CALL REALLOCATE_1D_1(A9, N)
  CALL REALLOCATE_1D_1(A10, N)
END SUBROUTINE

! -------- 2D --------
!> \brief Allocate or reallocate a 2D DOUBLE PRECISION array.
!! \param[inout] A1 Array to (re)allocate
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_1(A1, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  IF (.NOT. ALLOCATED(A1) .OR. ANY(SHAPE(A1) /= [N1, N2])) THEN
    IF (ALLOCATED(A1)) DEALLOCATE(A1)
    ALLOCATE(A1(N1, N2))
  END IF
END SUBROUTINE

!> \brief Allocate or reallocate two 2D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_2(A1, A2, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:), A2(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  CALL REALLOCATE_2D_1(A1, N1, N2)
  CALL REALLOCATE_2D_1(A2, N1, N2)
END SUBROUTINE

!> \brief Allocate or reallocate three 2D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_3(A1, A2, A3, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:), A2(:,:), A3(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  CALL REALLOCATE_2D_1(A1, N1, N2)
  CALL REALLOCATE_2D_1(A2, N1, N2)
  CALL REALLOCATE_2D_1(A3, N1, N2)
END SUBROUTINE

!> \brief Allocate or reallocate four 2D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_4(A1, A2, A3, A4, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:), A2(:,:), A3(:,:), A4(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  CALL REALLOCATE_2D_1(A1, N1, N2)
  CALL REALLOCATE_2D_1(A2, N1, N2)
  CALL REALLOCATE_2D_1(A3, N1, N2)
  CALL REALLOCATE_2D_1(A4, N1, N2)
END SUBROUTINE

!> \brief Allocate or reallocate five 2D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_5(A1, A2, A3, A4, A5, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:), A2(:,:), A3(:,:), A4(:,:), A5(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  CALL REALLOCATE_2D_1(A1, N1, N2)
  CALL REALLOCATE_2D_1(A2, N1, N2)
  CALL REALLOCATE_2D_1(A3, N1, N2)
  CALL REALLOCATE_2D_1(A4, N1, N2)
  CALL REALLOCATE_2D_1(A5, N1, N2)
END SUBROUTINE

!> \brief Allocate or reallocate six 2D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_6(A1, A2, A3, A4, A5, A6, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:), A2(:,:), A3(:,:), A4(:,:), A5(:,:), A6(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  CALL REALLOCATE_2D_1(A1, N1, N2)
  CALL REALLOCATE_2D_1(A2, N1, N2)
  CALL REALLOCATE_2D_1(A3, N1, N2)
  CALL REALLOCATE_2D_1(A4, N1, N2)
  CALL REALLOCATE_2D_1(A5, N1, N2)
  CALL REALLOCATE_2D_1(A6, N1, N2)
END SUBROUTINE

!> \brief Allocate or reallocate seven 2D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[inout] A7 Seventh array
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_7(A1, A2, A3, A4, A5, A6, A7, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:), A2(:,:), A3(:,:), A4(:,:), A5(:,:), A6(:,:), A7(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  CALL REALLOCATE_2D_1(A1, N1, N2)
  CALL REALLOCATE_2D_1(A2, N1, N2)
  CALL REALLOCATE_2D_1(A3, N1, N2)
  CALL REALLOCATE_2D_1(A4, N1, N2)
  CALL REALLOCATE_2D_1(A5, N1, N2)
  CALL REALLOCATE_2D_1(A6, N1, N2)
  CALL REALLOCATE_2D_1(A7, N1, N2)
END SUBROUTINE

!> \brief Allocate or reallocate eight 2D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[inout] A7 Seventh array
!! \param[inout] A8 Eighth array
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_8(A1, A2, A3, A4, A5, A6, A7, A8, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:), A2(:,:), A3(:,:), A4(:,:), A5(:,:), A6(:,:), A7(:,:), A8(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  CALL REALLOCATE_2D_1(A1, N1, N2)
  CALL REALLOCATE_2D_1(A2, N1, N2)
  CALL REALLOCATE_2D_1(A3, N1, N2)
  CALL REALLOCATE_2D_1(A4, N1, N2)
  CALL REALLOCATE_2D_1(A5, N1, N2)
  CALL REALLOCATE_2D_1(A6, N1, N2)
  CALL REALLOCATE_2D_1(A7, N1, N2)
  CALL REALLOCATE_2D_1(A8, N1, N2)
END SUBROUTINE

!> \brief Allocate or reallocate nine 2D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[inout] A7 Seventh array
!! \param[inout] A8 Eighth array
!! \param[inout] A9 Ninth array
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_9(A1, A2, A3, A4, A5, A6, A7, A8, A9, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:), A2(:,:), A3(:,:), A4(:,:), A5(:,:), A6(:,:), A7(:,:), A8(:,:), A9(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  CALL REALLOCATE_2D_1(A1, N1, N2)
  CALL REALLOCATE_2D_1(A2, N1, N2)
  CALL REALLOCATE_2D_1(A3, N1, N2)
  CALL REALLOCATE_2D_1(A4, N1, N2)
  CALL REALLOCATE_2D_1(A5, N1, N2)
  CALL REALLOCATE_2D_1(A6, N1, N2)
  CALL REALLOCATE_2D_1(A7, N1, N2)
  CALL REALLOCATE_2D_1(A8, N1, N2)
  CALL REALLOCATE_2D_1(A9, N1, N2)
END SUBROUTINE

!> \brief Allocate or reallocate ten 2D DOUBLE PRECISION arrays.
!! \param[inout] A1 First array
!! \param[inout] A2 Second array
!! \param[inout] A3 Third array
!! \param[inout] A4 Fourth array
!! \param[inout] A5 Fifth array
!! \param[inout] A6 Sixth array
!! \param[inout] A7 Seventh array
!! \param[inout] A8 Eighth array
!! \param[inout] A9 Ninth array
!! \param[inout] A10 Tenth array
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
SUBROUTINE REALLOCATE_2D_10(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, N1, N2)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:), A2(:,:), A3(:,:), A4(:,:), A5(:,:), &
                                                  A6(:,:), A7(:,:), A8(:,:), A9(:,:), A10(:,:)
  INTEGER, INTENT(IN) :: N1, N2
  CALL REALLOCATE_2D_1(A1, N1, N2)
  CALL REALLOCATE_2D_1(A2, N1, N2)
  CALL REALLOCATE_2D_1(A3, N1, N2)
  CALL REALLOCATE_2D_1(A4, N1, N2)
  CALL REALLOCATE_2D_1(A5, N1, N2)
  CALL REALLOCATE_2D_1(A6, N1, N2)
  CALL REALLOCATE_2D_1(A7, N1, N2)
  CALL REALLOCATE_2D_1(A8, N1, N2)
  CALL REALLOCATE_2D_1(A9, N1, N2)
  CALL REALLOCATE_2D_1(A10, N1, N2)
END SUBROUTINE

! -------- 3D --------
!> \brief Allocate or reallocate a 3D DOUBLE PRECISION array.
!! \param[inout] A1 Array to (re)allocate
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
!! \param[in] N3 Third dimension
SUBROUTINE REALLOCATE_3D_1(A1, N1, N2, N3)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:,:)
  INTEGER, INTENT(IN) :: N1, N2, N3
  IF (.NOT. ALLOCATED(A1) .OR. ANY(SHAPE(A1) /= [N1, N2, N3])) THEN
    IF (ALLOCATED(A1)) DEALLOCATE(A1)
    ALLOCATE(A1(N1, N2, N3))
  END IF
END SUBROUTINE

! -------- 4D --------
!> \brief Allocate or reallocate a 4D DOUBLE PRECISION array.
!! \param[inout] A1 Array to (re)allocate
!! \param[in] N1 First dimension
!! \param[in] N2 Second dimension
!! \param[in] N3 Third dimension
!! \param[in] N4 Fourth dimension
SUBROUTINE REALLOCATE_4D_1(A1, N1, N2, N3, N4)
  DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: A1(:,:,:,:)
  INTEGER, INTENT(IN) :: N1, N2, N3, N4
  IF (.NOT. ALLOCATED(A1) .OR. ANY(SHAPE(A1) /= [N1, N2, N3, N4])) THEN
    IF (ALLOCATED(A1)) DEALLOCATE(A1)
    ALLOCATE(A1(N1, N2, N3, N4))
  END IF
END SUBROUTINE

END MODULE REALLOCATE_UTILS
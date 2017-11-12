      SUBROUTINE WEIGHT (T, W)
C----------------------------------------------------------------------
C
C     Computes the weights at 0 <= t <= 1 for the spline function and
C     its 1st derivative defined by:
C
C     x  (t) = x1*W11(t) + x2*W21(t) + dx1*W31(t) + dx2*W41(t)
C     x' (t) = x1*W12(t) + x2*W22(t) + dx1*W32(t) + dx2*W42(t)
C     x''(t) = x1*W13(t) + x2*W23(t) + dx1*W33(t) + dx2*W43(t)
C
C                           2                  2
C     W (t) = (1 + 2t) (t-1),   W (t)=  t (t-1)
C      11                        31
C                       2                     2
C     W (t) = (3 - 2t) t    ,   W (t)= (t-1) t
C      21                        41
C
C
C     For t=0 W (:,1)= [ 1, 0, 0, 0],
C             W (:,2)= [ 0, 0, 1, 0],
C             W (:,3)= [-6, 6,-4,-2]
C
C     For t=1 W (:,1)= [ 0, 1, 0, 0],
C             W (:,2)= [ 0, 0, 0, 1],
C             W (:,3)= [ 6,-6, 2, 4]
C
C     x(t), x'(t) and x''(t) are continuous over an element
C     truncation error is O(h^4)
C
C----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      REAL(KIND=RK) T
      REAL(KIND=RK) W (4, 2)
C
      W (1, 1) = (T - 1.0D0) ** 2 * (1.0D0 + 2.0D0 * T)
      W (2, 1) = T ** 2           * (3.0D0 - 2.0D0 * T)
      W (3, 1) = (T - 1.0D0) ** 2 *  T
      W (4, 1) = T ** 2           * (T - 1.0D0)

      W (1, 2) = 6.0D0 * T * (T - 1.0D0)
      W (2, 2) =-6.0D0 * T * (T - 1.0D0)
      W (3, 2) = (3.0D0 * T - 1.0D0) * (T - 1.0D0)
      W (4, 2) = (3.0D0 * T - 2.0D0) *  T
C
      RETURN
      END



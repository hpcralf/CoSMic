!     From Leonard J. Moss of SLAC:
!
! Here's a hybrid QuickSort I wrote a number of years ago.  It's
! based on suggestions in Knuth, Volume 3, and performs much better
! than a pure QuickSort on short or partially ordered input arrays.  
!
Module quicksort_nr

  Implicit None

Contains

  Subroutine SORTRX(N,Data,INDEX)

    !===================================================================
    !
    !     SORTRX -- SORT, Real input, indeX output
    !
    !
    !     Input:  N     INTEGER
    !             DATA  REAL
    !
    !     Output: INDEX INTEGER (DIMENSION N)
    !
    ! This routine performs an in-memory sort of the first N elements of
    ! array DATA, returning into array INDEX the indices of elements of
    ! DATA arranged in ascending order.  Thus,
    !
    !    DATA(INDEX(1)) will be the smallest number in array DATA;
    !    DATA(INDEX(N)) will be the largest number in DATA.
    !
    ! The original data is not physically rearranged.  The original order
    ! of equal input values is not necessarily preserved.
    !
    !===================================================================
    !
    ! SORTRX uses a hybrid QuickSort algorithm, based on several
    ! suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
    ! "pivot key" [my term] for dividing each subsequence is chosen to be
    ! the median of the first, last, and middle values of the subsequence;
    ! and the QuickSort is cut off when a subsequence has 9 or fewer
    ! elements, and a straight insertion sort of the entire array is done
    ! at the end.  The result is comparable to a pure insertion sort for
    ! very short arrays, and very fast for very large arrays (of order 12
    ! micro-sec/element on the 3081K for arrays of 10K elements).  It is
    ! also not subject to the poor performance of the pure QuickSort on
    ! partially ordered data.
    !
    ! Created:  15 Jul 1986  Len Moss
    !
    !===================================================================

    Integer   N,Index(N)
    Integer      Data(N)

    Integer   LSTK(31),RSTK(31),ISTK
    Integer   L,R,I,J,P,INDEXP,INDEXT
    Real      DATAP

    !     QuickSort Cutoff
    !
    !     Quit QuickSort-ing when a subsequence contains M or fewer
    !     elements and finish off at end with straight insertion sort.
    !     According to Knuth, V.3, the optimum value of M is around 9.

    Integer   M
    Parameter (M=9)

    !===================================================================
    !
    !     Make initial guess for INDEX

    Do I=1,N
       Index(I)=I
    End Do

    !     If array is short, skip QuickSort and go directly to
    !     the straight insertion sort.

    If (N.Le.M) Goto 900

    !===================================================================
    !
    !     QuickSort
    !
    !     The "Qn:"s correspond roughly to steps in Algorithm Q,
    !     Knuth, V.3, PP.116-117, modified to select the median
    !     of the first, last, and middle elements as the "pivot
    !     key" (in Knuth's notation, "K").  Also modified to leave
    !     data in place and produce an INDEX array.  To simplify
    !     comments, let DATA[I]=DATA(INDEX(I)).

    ! Q1: Initialize
    ISTK=0
    L=1
    R=N

200 Continue

    ! Q2: Sort the subsequence DATA[L]..DATA[R].
    !
    !     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
    !     r > R, and L <= m <= R.  (First time through, there is no
    !     DATA for l < L or r > R.)

    I=L
    J=R

    ! Q2.5: Select pivot key
    !
    !     Let the pivot, P, be the midpoint of this subsequence,
    !     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
    !     so the corresponding DATA values are in increasing order.
    !     The pivot key, DATAP, is then DATA[P].

    P=(L+R)/2
    INDEXP=Index(P)
    DATAP=Data(INDEXP)

    If (Data(Index(L)) .Gt. DATAP) Then
       Index(P)=Index(L)
       Index(L)=INDEXP
       INDEXP=Index(P)
       DATAP=Data(INDEXP)
    Endif

    If (DATAP .Gt. Data(Index(R))) Then
       If (Data(Index(L)) .Gt. Data(Index(R))) Then
          Index(P)=Index(L)
          Index(L)=Index(R)
       Else
          Index(P)=Index(R)
       Endif
       Index(R)=INDEXP
       INDEXP=Index(P)
       DATAP=Data(INDEXP)
    Endif

    !     Now we swap values between the right and left sides and/or
    !     move DATAP until all smaller values are on the left and all
    !     larger values are on the right.  Neither the left or right
    !     side will be internally ordered yet; however, DATAP will be
    !     in its final position.

300 Continue

    ! Q3: Search for datum on left >= DATAP
    !
    !     At this point, DATA[L] <= DATAP.  We can therefore start scanning
    !     up from L, looking for a value >= DATAP (this scan is guaranteed
    !     to terminate since we initially placed DATAP near the middle of
    !     the subsequence).

    I=I+1
    If (Data(Index(I)).Lt.DATAP) Goto 300

400 Continue

    ! Q4: Search for datum on right <= DATAP
    !
    !     At this point, DATA[R] >= DATAP.  We can therefore start scanning
    !     down from R, looking for a value <= DATAP (this scan is guaranteed
    !     to terminate since we initially placed DATAP near the middle of
    !     the subsequence).

    J=J-1
    If (Data(Index(J)).Gt.DATAP) Goto 400

    ! Q5: Have the two scans collided?

    If (I.Lt.J) Then

       ! Q6: No, interchange DATA[I] <--> DATA[J] and continue

       INDEXT=Index(I)
       Index(I)=Index(J)
       Index(J)=INDEXT
       Goto 300
    Else

       ! Q7: Yes, select next subsequence to sort
       !
       !     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
       !     for all L <= l < I and J < r <= R.  If both subsequences are
       !     more than M elements long, push the longer one on the stack and
       !     go back to QuickSort the shorter; if only one is more than M
       !     elements long, go back and QuickSort it; otherwise, pop a
       !     subsequence off the stack and QuickSort it.

       If (R-J .Ge. I-L .And. I-L .Gt. M) Then
          ISTK=ISTK+1
          LSTK(ISTK)=J+1
          RSTK(ISTK)=R
          R=I-1
       Else If (I-L .Gt. R-J .And. R-J .Gt. M) Then
          ISTK=ISTK+1
          LSTK(ISTK)=L
          RSTK(ISTK)=I-1
          L=J+1
       Else If (R-J .Gt. M) Then
          L=J+1
       Else If (I-L .Gt. M) Then
          R=I-1
       Else
          ! Q8: Pop the stack, or terminate QuickSort if empty
          If (ISTK.Lt.1) Goto 900
          L=LSTK(ISTK)
          R=RSTK(ISTK)
          ISTK=ISTK-1
       Endif
       Goto 200
    Endif

900 Continue

    !===================================================================
    !
    ! Q9: Straight Insertion sort

    Do I=2,N
       If (Data(Index(I-1)) .Gt. Data(Index(I))) Then
          INDEXP=Index(I)
          DATAP=Data(INDEXP)
          P=I-1
920       Continue
          Index(P+1) = Index(P)
          P=P-1
          If (P.Gt.0) Then
             If (Data(Index(P)).Gt.DATAP) Goto 920
          Endif
          Index(P+1) = INDEXP
       Endif
    End Do

    !===================================================================
    !
    !     All done

  End Subroutine SORTRX

  Subroutine SORTRX_REAL(N,Data,INDEX)

    !===================================================================
    !
    !     SORTRX -- SORT, Real input, indeX output
    !
    !
    !     Input:  N     INTEGER
    !             DATA  REAL
    !
    !     Output: INDEX INTEGER (DIMENSION N)
    !
    ! This routine performs an in-memory sort of the first N elements of
    ! array DATA, returning into array INDEX the indices of elements of
    ! DATA arranged in ascending order.  Thus,
    !
    !    DATA(INDEX(1)) will be the smallest number in array DATA;
    !    DATA(INDEX(N)) will be the largest number in DATA.
    !
    ! The original data is not physically rearranged.  The original order
    ! of equal input values is not necessarily preserved.
    !
    !===================================================================
    !
    ! SORTRX uses a hybrid QuickSort algorithm, based on several
    ! suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
    ! "pivot key" [my term] for dividing each subsequence is chosen to be
    ! the median of the first, last, and middle values of the subsequence;
    ! and the QuickSort is cut off when a subsequence has 9 or fewer
    ! elements, and a straight insertion sort of the entire array is done
    ! at the end.  The result is comparable to a pure insertion sort for
    ! very short arrays, and very fast for very large arrays (of order 12
    ! micro-sec/element on the 3081K for arrays of 10K elements).  It is
    ! also not subject to the poor performance of the pure QuickSort on
    ! partially ordered data.
    !
    ! Created:  15 Jul 1986  Len Moss
    !
    !===================================================================

    Integer   N,Index(N)
    Real      Data(N)

    Integer   LSTK(31),RSTK(31),ISTK
    Integer   L,R,I,J,P,INDEXP,INDEXT
    Real      DATAP

    !     QuickSort Cutoff
    !
    !     Quit QuickSort-ing when a subsequence contains M or fewer
    !     elements and finish off at end with straight insertion sort.
    !     According to Knuth, V.3, the optimum value of M is around 9.

    Integer   M
    Parameter (M=9)

    !===================================================================
    !
    !     Make initial guess for INDEX

    Do I=1,N
       Index(I)=I
    End Do

    !     If array is short, skip QuickSort and go directly to
    !     the straight insertion sort.

    If (N.Le.M) Goto 900

    !===================================================================
    !
    !     QuickSort
    !
    !     The "Qn:"s correspond roughly to steps in Algorithm Q,
    !     Knuth, V.3, PP.116-117, modified to select the median
    !     of the first, last, and middle elements as the "pivot
    !     key" (in Knuth's notation, "K").  Also modified to leave
    !     data in place and produce an INDEX array.  To simplify
    !     comments, let DATA[I]=DATA(INDEX(I)).

    ! Q1: Initialize
    ISTK=0
    L=1
    R=N

200 Continue

    ! Q2: Sort the subsequence DATA[L]..DATA[R].
    !
    !     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
    !     r > R, and L <= m <= R.  (First time through, there is no
    !     DATA for l < L or r > R.)

    I=L
    J=R

    ! Q2.5: Select pivot key
    !
    !     Let the pivot, P, be the midpoint of this subsequence,
    !     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
    !     so the corresponding DATA values are in increasing order.
    !     The pivot key, DATAP, is then DATA[P].

    P=(L+R)/2
    INDEXP=Index(P)
    DATAP=Data(INDEXP)

    If (Data(Index(L)) .Gt. DATAP) Then
       Index(P)=Index(L)
       Index(L)=INDEXP
       INDEXP=Index(P)
       DATAP=Data(INDEXP)
    Endif

    If (DATAP .Gt. Data(Index(R))) Then
       If (Data(Index(L)) .Gt. Data(Index(R))) Then
          Index(P)=Index(L)
          Index(L)=Index(R)
       Else
          Index(P)=Index(R)
       Endif
       Index(R)=INDEXP
       INDEXP=Index(P)
       DATAP=Data(INDEXP)
    Endif

    !     Now we swap values between the right and left sides and/or
    !     move DATAP until all smaller values are on the left and all
    !     larger values are on the right.  Neither the left or right
    !     side will be internally ordered yet; however, DATAP will be
    !     in its final position.

300 Continue

    ! Q3: Search for datum on left >= DATAP
    !
    !     At this point, DATA[L] <= DATAP.  We can therefore start scanning
    !     up from L, looking for a value >= DATAP (this scan is guaranteed
    !     to terminate since we initially placed DATAP near the middle of
    !     the subsequence).

    I=I+1
    If (Data(Index(I)).Lt.DATAP) Goto 300

400 Continue

    ! Q4: Search for datum on right <= DATAP
    !
    !     At this point, DATA[R] >= DATAP.  We can therefore start scanning
    !     down from R, looking for a value <= DATAP (this scan is guaranteed
    !     to terminate since we initially placed DATAP near the middle of
    !     the subsequence).

    J=J-1
    If (Data(Index(J)).Gt.DATAP) Goto 400

    ! Q5: Have the two scans collided?

    If (I.Lt.J) Then

       ! Q6: No, interchange DATA[I] <--> DATA[J] and continue

       INDEXT=Index(I)
       Index(I)=Index(J)
       Index(J)=INDEXT
       Goto 300
    Else

       ! Q7: Yes, select next subsequence to sort
       !
       !     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
       !     for all L <= l < I and J < r <= R.  If both subsequences are
       !     more than M elements long, push the longer one on the stack and
       !     go back to QuickSort the shorter; if only one is more than M
       !     elements long, go back and QuickSort it; otherwise, pop a
       !     subsequence off the stack and QuickSort it.

       If (R-J .Ge. I-L .And. I-L .Gt. M) Then
          ISTK=ISTK+1
          LSTK(ISTK)=J+1
          RSTK(ISTK)=R
          R=I-1
       Else If (I-L .Gt. R-J .And. R-J .Gt. M) Then
          ISTK=ISTK+1
          LSTK(ISTK)=L
          RSTK(ISTK)=I-1
          L=J+1
       Else If (R-J .Gt. M) Then
          L=J+1
       Else If (I-L .Gt. M) Then
          R=I-1
       Else
          ! Q8: Pop the stack, or terminate QuickSort if empty
          If (ISTK.Lt.1) Goto 900
          L=LSTK(ISTK)
          R=RSTK(ISTK)
          ISTK=ISTK-1
       Endif
       Goto 200
    Endif

900 Continue

    !===================================================================
    !
    ! Q9: Straight Insertion sort

    Do I=2,N
       If (Data(Index(I-1)) .Gt. Data(Index(I))) Then
          INDEXP=Index(I)
          DATAP=Data(INDEXP)
          P=I-1
920       Continue
          Index(P+1) = Index(P)
          P=P-1
          If (P.Gt.0) Then
             If (Data(Index(P)).Gt.DATAP) Goto 920
          Endif
          Index(P+1) = INDEXP
       Endif
    End Do

    !===================================================================
    !
    !     All done

  End Subroutine SORTRX_REAL
  
End Module quicksort_nr


C *******************************************************************
C COPYRIGHT (c) 1998 Council for the Central Laboratory
*                    of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
C SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
C
C Please note that for an ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date 16 Jan 1998
C MC60F/FD corrected 7 Feb. 2002
C MC60G/GD amended 7 Feb. 2002 to return mean frontal matrix size
C 20/2/02 Cosmetic changes applied to reduce single/double differences
C 17/7/04 MC60C/CD. Changed PAIR(2,NSUP/2) to PAIR(2,*) (problem
C         previusly if NSUP=1).

C 17th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC60AD(N,LIRN,IRN,ICPTR,ICNTL,IW,INFO)
      INTEGER N
      INTEGER LIRN
      INTEGER IRN(LIRN)
      INTEGER ICPTR(N+1)
      INTEGER ICNTL(2)
      INTEGER IW(N)
      INTEGER INFO(4)
      INTEGER CKP1,I,I1,I2,II,IOUT,IPOS,IREP,J,KZ,LP,NDIAG,NEWTAU
      LP = ICNTL(2)
      DO 5 J = 1,4
         INFO(J) = 0
    5 CONTINUE
      IF (N.LT.1) THEN
          INFO(1) = -1
          IF (LP.GT.0) WRITE (LP,'(/,A,I3/,A,I6)')
     *         ' MC60A/AD error: INFO(1) =', INFO(1), ' N =',N
          RETURN
      END IF
      IF (LIRN.LT.ICPTR(N+1)-1) THEN
          INFO(1) = -1
          IF (LP.GT.0) WRITE (LP,'(/,A,I3/,A)')
     *         ' MC60A/AD error:  INFO(1) =', INFO(1),
     *         ' LIRN is less than ICPTR(N+1)-1'
          RETURN
      END IF
      DO 10 I = 1,N
        IW(I) = 0
   10 CONTINUE
      IOUT = 0
      IREP = 0
      KZ = 0
      IF (ICNTL(1).EQ.1) THEN
         I1 = ICPTR(1)
         ICPTR(1) = 1
         DO 16 J = 1,N
            DO 15 II = I1,ICPTR(J+1) - 1
               I = IRN(II)
               IF (I.GT.N .OR. I.LT.J) THEN
                  IOUT = IOUT + 1
               ELSE IF (IW(I).EQ.J) THEN
                  IREP = IREP + 1
               ELSE
                  KZ = KZ + 1
                  IRN(KZ) = I
                  IW(I) = J
               END IF
   15       CONTINUE
            I1 = ICPTR(J+1)
            ICPTR(J+1) = KZ + 1
   16    CONTINUE
         IF (IOUT.GT.0) THEN
             INFO(1) = 1
             IF (LP.GT.0) WRITE (LP,'(/,A,I6,A)')
     *         ' MC60A/AD warning:',IOUT,' out-of-range entries ignored'
         END IF
         IF (IREP.GT.0) THEN
             INFO(1) = 1
             IF (LP.GT.0) WRITE (LP,'(/,A,I6,A)')
     *         ' MC60A/AD warning:',IREP,' duplicated entries ignored'
         END IF
         INFO(2) = IOUT
         INFO(3) = IREP
      ELSE
         I1 = ICPTR(1)
         DO 26 J = 1,N
            DO 25 II = I1,ICPTR(J+1) - 1
               I = IRN(II)
               IF (I.GT.N .OR. I.LT.J) THEN
                  IOUT = IOUT + 1
               ELSE IF (IW(I).EQ.J) THEN
                  IREP = IREP + 1
               ELSE
                  KZ = KZ + 1
                  IW(I) = J
               END IF
   25       CONTINUE
            I1 = ICPTR(J+1)
   26    CONTINUE
         IF (IOUT.GT.0 .OR. IREP.GT.0) THEN
            INFO(1) = -3
            IF (LP.GT.0) THEN
              WRITE (LP,'(/,A,I3)')
     *            ' MC60A/AD error:  INFO(1) =', INFO(1)
            IF (IOUT.GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',IOUT,' out-of-range entries'
            IF (IREP.GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',IREP,' duplicated entries'
            END IF
            INFO(2) = IOUT
            INFO(3) = IREP
            RETURN
         END IF
      END IF
      DO 30 J = 1,N
        IW(J) = 0
   30 CONTINUE
      NDIAG = 0
      DO 40 J = 1,N
        I1 = ICPTR(J)
        I2 = ICPTR(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 35 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1
          ELSE
            NDIAG = NDIAG + 1
          END IF
   35   CONTINUE
   40 CONTINUE
      NEWTAU = 2*KZ - NDIAG
      INFO(4) = NEWTAU
      IF (NEWTAU.GT.LIRN) THEN
          INFO(1) = -2
          IF (LP.GT.0) WRITE (LP,'(/,A)')
     *         ' MC60A/AD error: LIRN is too small'
          RETURN
      END IF
      I1 = KZ + 1
      CKP1 = NEWTAU + 1
      DO 60 J = N,1,-1
        I2 = I1 - 1
        I1 = ICPTR(J)
        IPOS = CKP1
          DO 50 II = I2,I1,-1
            IPOS = IPOS - 1
            IRN(IPOS) = IRN(II)
   50     CONTINUE
        ICPTR(J) = IPOS
        CKP1 = CKP1 - IW(J)
        IW(J) = I2 - I1 + 1
   60 CONTINUE
      DO 80 J = N,1,-1
        I1 = ICPTR(J)
        I2 = ICPTR(J) + IW(J) - 1
        IF(I1.LE.I2) THEN
          DO 70 II = I1,I2
            I = IRN(II)
            IF(I.EQ.J) GO TO 70
            ICPTR(I) = ICPTR(I) - 1
            IRN(ICPTR(I)) = J
   70     CONTINUE
        END IF
   80 CONTINUE
      ICPTR(N+1) = NEWTAU + 1
      END
      SUBROUTINE MC60BD(N,LIRN,IRN,ICPTR,NSUP,SVAR,VARS,IW)
      INTEGER N
      INTEGER LIRN
      INTEGER IRN(LIRN)
      INTEGER ICPTR(N+1)
      INTEGER NSUP
      INTEGER SVAR(N)
      INTEGER VARS(N)
      INTEGER IW(2*N)
      INTEGER FLAG
      EXTERNAL MC60OD,MC60PD
      FLAG = N+1
      CALL MC60OD(N,N,LIRN,IRN,ICPTR,SVAR,NSUP,IW,VARS,IW(FLAG))
      CALL MC60PD(N,NSUP,LIRN,IRN,ICPTR,SVAR,VARS,IW,IW(FLAG))
      END
      SUBROUTINE MC60CD(N,NSUP,LIRN,IRN,ICPTR,VARS,JCNTL,
     +                  PERMSV,WEIGHT,PAIR,INFO,IW,W)
      INTEGER N,NSUP,LIRN,IRN(LIRN),ICPTR(NSUP+1),VARS(NSUP),JCNTL(2)
      INTEGER PERMSV(NSUP),PAIR(2,*),INFO(4),IW(3*NSUP+1)
      DOUBLE PRECISION WEIGHT(2),W(NSUP)
      INTEGER DEGREE,I,IL,HINFO(6),J,K,LIST,LSTNUM,LWIDTH,LWDTH1,
     *        LZNUM,MAXPSV,MINPSV,NLVL,NLVL1,NODES,NSTOP,NSTRT,NVARS,XLS
      DOUBLE PRECISION NWGHT(2)
      EXTERNAL MC60HD,MC60JD,MC60LD
      XLS = NSUP
      LIST = 2*NSUP + 1
      NWGHT(1) = WEIGHT(1)
      NWGHT(2) = WEIGHT(2)
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      INFO(4) = 0
      NVARS = 0
      LSTNUM = 0
      LZNUM = NSUP+1
      IF(JCNTL(2).NE.2) THEN
         DO 5 I = 1,NSUP
            PERMSV(I) = 1
    5    CONTINUE
      END IF
      DO 6 I = 1,NSUP
         K = ICPTR(I)
         DEGREE = ICPTR(I+1) - K
         IF (DEGREE.LE.1) THEN
           IF (DEGREE.EQ.0) THEN
             LZNUM = LZNUM - 1
             PERMSV(I) = -LZNUM
           ELSE IF (IRN(K).EQ.I) THEN
              LSTNUM = LSTNUM + 1
              PERMSV(I) = -LSTNUM
            END IF
         END IF
    6 CONTINUE
      DO 30 I = 1, NSUP
        IF (LSTNUM.GE.LZNUM-1) GO TO 35
        IF(JCNTL(2).EQ.0) THEN
          CALL MC60HD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +               IW(XLS+1),IW(LIST+1),HINFO)
          NSTRT = HINFO(1)
          PAIR(1,I) = HINFO(1)
          PAIR(2,I) = HINFO(2)
          NLVL = HINFO(3)
          LWIDTH = HINFO(4)
          NVARS = HINFO(5)
          NODES = HINFO(6)
        ELSE IF(JCNTL(2).EQ.1) THEN
          NSTOP = PAIR(1,I)
          CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          LWDTH1 = LWIDTH
          NLVL1 = NLVL
          NODES = IW(XLS+NLVL+1)-1
          NSTRT = NSTOP
          NSTOP = PAIR(2,I)
          CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          IF (NLVL1.GT.NLVL .OR.
     +             (NLVL1.EQ.NLVL .AND. LWDTH1.LT.LWIDTH)) THEN
             NSTRT = NSTOP
             NSTOP = PAIR(1,I)
             CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          END IF
        ELSE
          MAXPSV = 0
          DO 10 J = 1,NSUP
            IF (PERMSV(J).GT.0) THEN
               IF (PERMSV(J).GT.MAXPSV) THEN
                  IF (MAXPSV.EQ.0) THEN
                    MINPSV = PERMSV(J)
                    NSTRT = J
                  END IF
                  MAXPSV = PERMSV(J)
               END IF
               IF (PERMSV(J).LT.MINPSV) THEN
                  MINPSV = PERMSV(J)
                  NSTRT = J
               END IF
            END IF
   10     CONTINUE
          CALL MC60LD(NSTRT,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          NODES = IW(XLS+NLVL+1)-1
          IF(MAXPSV.NE.MINPSV) THEN
            NWGHT(2) = (WEIGHT(2)*(NLVL-1))/(MAXPSV-MINPSV)
          ELSE
            NWGHT(2) = WEIGHT(2)
          END IF
        END IF
        INFO(1) = INFO(1) + 1
        IF (NVARS.GT.INFO(2)) THEN
           INFO(2) = NVARS
           INFO(3) = NLVL
           INFO(4) = LWIDTH
        END IF
        IF(JCNTL(1).EQ.1) THEN
           DO 11 J = NODES,1,-1
             LSTNUM = LSTNUM + 1
             PERMSV(IW(J)) = -LSTNUM
   11      CONTINUE
        ELSE
          IF(JCNTL(2).NE.2) THEN
            DO 15 IL = 1,NLVL
              DO 12 J = IW(XLS+IL), IW(XLS+IL+1) - 1
                PERMSV(IW(J)) = NLVL - IL
   12         CONTINUE
   15       CONTINUE
          END IF
          CALL MC60JD(NSUP,LIRN,NODES,NSTRT,LSTNUM,IRN,ICPTR,
     +         VARS,PERMSV,NWGHT,IW,IW(LIST+1),IW(XLS+1),W)
        END IF
   30 CONTINUE
   35 DO 40 I = 1,NSUP
          PERMSV(I) = -PERMSV(I)
   40 CONTINUE
      END
      SUBROUTINE MC60DD(N,NSUP,SVAR,VARS,PERMSV,PERM,POSSV)
      INTEGER N
      INTEGER NSUP
      INTEGER SVAR(N)
      INTEGER VARS(NSUP)
      INTEGER PERMSV(NSUP)
      INTEGER PERM(N)
      INTEGER POSSV(NSUP)
      INTEGER I,IS,L
      DO 10 IS = 1,NSUP
         PERM(PERMSV(IS)) = IS
   10 CONTINUE
      L = 1
      DO 20 I = 1,NSUP
         L = L + VARS(PERM(I))
         POSSV(PERM(I)) = L
   20 CONTINUE
      DO 30 I = 1,N
         IS = SVAR(I)
         L = POSSV(IS)-1
         POSSV(IS) = L
         PERM(I) = L
   30 CONTINUE
      END
      SUBROUTINE MC60ED(N,NSUP,LIRN,IRN,ICPTR,SVAR,VARS,PERMSV,PERM,IW)
      INTEGER N,NSUP,LIRN
      INTEGER IRN(LIRN),ICPTR(NSUP+1),PERM(N),PERMSV(NSUP),IW(NSUP),
     *        SVAR(N),VARS(NSUP)
      INTEGER I,IS,JS,K,L,M
      DO 10 IS = 1,NSUP
        JS = PERMSV(IS)
        IW(JS) = IS
        PERM(IS) = 0
   10 CONTINUE
      K = 0
      DO 30 L = 1,NSUP
        IS = IW(L)
        DO 20 M = ICPTR(IS), ICPTR(IS+1)-1
          JS = IRN(M)
          IF(PERM(JS).EQ.0) THEN
            K = K + 1
            PERM(JS) = K
            PERMSV(K) = JS
          END IF
   20   CONTINUE
        IF(PERM(IS).EQ.0) THEN
            K = K + 1
            PERM(IS) = K
            PERMSV(K) = IS
          END IF
   30 CONTINUE
      IF (N.EQ.NSUP) THEN
         DO 40 I = 1,N
            PERM(I) = PERMSV(I)
   40    CONTINUE
         RETURN
      END IF
      L = 1
      DO 45 IS = 1,NSUP
         JS = PERMSV(IS)
         L = L + VARS(JS)
         IW(JS) = L
   45 CONTINUE
      DO 50 I = 1,N
         IS = SVAR(I)
         L = IW(IS) - 1
         IW(IS) = L
         PERM(L) = I
   50 CONTINUE
      END
      SUBROUTINE MC60FD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,RINFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER LIRN,N,NSUP
      INTEGER IRN(LIRN),IW(2*NSUP+1),PERMSV(NSUP),ICPTR(NSUP+1),
     *        VARS(NSUP)
      DOUBLE PRECISION RINFO(4)
      INTEGER I,IMIN,J,JSTOP,JSTRT,K,NACTIV,NBR,NV
      INTRINSIC ABS,DBLE,MAX,MIN,SQRT
      DO 10 I = 1,4
         RINFO(I) = ZERO
  10  CONTINUE
      NACTIV = 0
      IW(NSUP+1) = 0
      DO 30 I = 1,NSUP
         J = PERMSV(I)
         IW(J) = I
   30 CONTINUE
      DO 80 I = 1,NSUP
         J = ABS(IW(I))
         NV = VARS(J)
         IW(NSUP+I+1) = IW(NSUP+I) + NV
         JSTRT = ICPTR(J)
         JSTOP = ICPTR(J+1) - 1
         IMIN = I + 1
         DO 50 K = JSTRT,JSTOP
            NBR = IRN(K)
            IMIN = MIN(IMIN,PERMSV(NBR))
            IF (IW(NBR).GT.0) THEN
               NACTIV = NACTIV +  VARS(NBR)
               IW(NBR) = -IW(NBR)
            END IF
   50    CONTINUE
         RINFO(3) = MAX(RINFO(3),DBLE(IW(NSUP+I+1)-IW(NSUP+IMIN)))
         IF (IW(J).GT.0) THEN
           IW(J) = -IW(J)
           RINFO(2) = MAX(RINFO(2),DBLE(NACTIV+1))
           RINFO(1) = RINFO(1) + NV*DBLE(NACTIV+1)
           RINFO(4) = RINFO(4) + NV*DBLE(NACTIV+1)**2
         ELSE
           RINFO(2) = MAX(RINFO(2),DBLE(NACTIV))
           DO 70 J = 1, NV
             RINFO(1) = RINFO(1) + DBLE(NACTIV)
             RINFO(4) = RINFO(4) + DBLE(NACTIV)**2
             NACTIV = NACTIV - 1
   70      CONTINUE
         END IF
   80 CONTINUE
      RINFO(4) = SQRT(RINFO(4)/DBLE(N))
      END
      SUBROUTINE MC60GD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,RINFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER LIRN,N,NSUP
      DOUBLE PRECISION RINFO(4)
      INTEGER ICPTR(NSUP+1),IRN(LIRN),IW(NSUP),PERMSV(NSUP),VARS(NSUP)
      INTEGER I,IEQ,J,JSTOP,JSTRT,K,KFRNT,LFRNT,MFR,NV
      INTRINSIC DBLE,MAX,SQRT
      DO 10 I = 1,4
         RINFO(I) = ZERO
   10 CONTINUE
      DO 20 I = 1,NSUP
         IW(I) = 0
   20 CONTINUE
      KFRNT = 0
      LFRNT = 0
      DO 40 IEQ = 1,NSUP
         I = PERMSV(IEQ)
         IW(I) = IEQ
         JSTRT = ICPTR(I)
         JSTOP = ICPTR(I+1) - 1
         DO 30 J = JSTRT,JSTOP
            MFR = IRN(J)
            IW(MFR) = IEQ
   30    CONTINUE
   40 CONTINUE
      DO 90 IEQ = 1,NSUP
         I = PERMSV(IEQ)
         JSTRT = ICPTR(I)
         JSTOP = ICPTR(I+1) - 1
         IF (JSTRT.GT.JSTOP) GO TO 90
         NV = VARS(I)
         KFRNT = KFRNT + NV
         IF (IW(I).GE.0) THEN
            LFRNT = LFRNT + NV
            IW(I) = -IW(I)
         END IF
         DO 60 J = JSTRT,JSTOP
            MFR = IRN(J)
            IF (IW(MFR).GE.0) THEN
               NV = VARS(MFR)
               LFRNT = LFRNT + NV
               IW(MFR) = -IW(MFR)
            END IF
   60    CONTINUE
         RINFO(1) = MAX(RINFO(1),DBLE(KFRNT))
         RINFO(2) = MAX(RINFO(2),DBLE(LFRNT))
         IF (-IW(I).EQ.IEQ) THEN
            IW(I) = 0
            NV = VARS(I)
            DO 65 K = 1,NV
               RINFO(3) = RINFO(3) + DBLE(KFRNT)**2
               RINFO(4) = RINFO(4) + DBLE(KFRNT)*DBLE(LFRNT)
               LFRNT = LFRNT - 1
               KFRNT = KFRNT - 1
   65       CONTINUE
         END IF
         DO 80 J = JSTRT,JSTOP
            MFR = IRN(J)
            IF (-IW(MFR).EQ.IEQ) THEN
               NV = VARS(MFR)
               DO 70 K = 1,NV
                  RINFO(3) = RINFO(3) + DBLE(KFRNT)**2
                  RINFO(4) = RINFO(4) + DBLE(KFRNT)*DBLE(LFRNT)
                  LFRNT = LFRNT - 1
                  KFRNT = KFRNT - 1
   70          CONTINUE
            END IF
   80    CONTINUE
   90 CONTINUE
      RINFO(3) = SQRT(RINFO(3)/DBLE(N))
      RINFO(4) = RINFO(4)/DBLE(N)
      END
      SUBROUTINE MC60HD(N,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,XLS,LIST,
     +                 INFO)
      INTEGER N,NSUP,LIRN
      INTEGER ICPTR(NSUP+1),IRN(LIRN),LIST(NSUP),LS(NSUP),
     +        MASK(NSUP),VARS(NSUP),XLS(NSUP+1),INFO(6)
      INTEGER DEGREE,I,J,LSIZE,LWIDTH,MAIN,MAXDEP,
     +        MINDEG,MINWID,NLSIZE,NLVL,NODE,NODES,NSTOP,NSTRT,NVARS
      EXTERNAL MC60LD
      MINDEG = N+1
      INFO(5) = 0
      DO 10 I = 1,NSUP
          IF (MASK(I).EQ.1) THEN
            INFO(5) = INFO(5) + VARS(I)
            DEGREE = ICPTR(I+1) - ICPTR(I)
            IF (DEGREE.LE.MINDEG) THEN
              IF (DEGREE.LT.MINDEG) THEN
                  NSTRT = I
                  MINDEG = DEGREE
              END IF
            END IF
          END IF
   10 CONTINUE
      CALL MC60LD(NSTRT,N,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,XLS,
     +            MAXDEP,LWIDTH,NVARS)
      NODES = XLS(MAXDEP+1) - 1
      NSTOP = 0
      DO 70 MAIN = 1, NODES
        INFO(4) = LWIDTH
        LSIZE = 0
        DO 30 I = XLS(MAXDEP),XLS(MAXDEP+1) - 1
          NODE = LS(I)
          LSIZE = LSIZE + 1
          LIST(LSIZE) = NODE
          XLS(NODE) = ICPTR(NODE+1) - ICPTR(NODE)
   30   CONTINUE
        DO 50 NLSIZE = 1,5
           MINDEG = N+1
           DO 41 I = NLSIZE,LSIZE
             IF(XLS(LIST(I)).LT.MINDEG)THEN
                J = I
                MINDEG = XLS(LIST(I))
             END IF
   41     CONTINUE
          IF(MINDEG.EQ.N+1) GO TO 55
          NODE = LIST(J)
          LIST(J) = LIST(NLSIZE)
          LIST(NLSIZE) = NODE
          DO 42 I = ICPTR(NODE), ICPTR(NODE+1)-1
             XLS(IRN(I)) = N+1
   42     CONTINUE
   50   CONTINUE
   55   NLSIZE = NLSIZE-1
        MINWID = N
        DO 60 I = 1,NLSIZE
          NODE = LIST(I)
          CALL MC60LD(NODE,MINWID,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,XLS,
     +               NLVL,LWIDTH,NVARS)
          IF (LWIDTH.LT.MINWID) THEN
            IF (NLVL.GT.MAXDEP) THEN
              NSTRT = NODE
              MAXDEP = NLVL
              GO TO 70
            ELSE
              NSTOP = NODE
              MINWID = LWIDTH
            END IF
          END IF
   60   CONTINUE
        GO TO 80
   70 CONTINUE
   80 IF (INFO(4) .LT. MINWID) THEN
         INFO(1) = NSTRT
         NSTRT = NSTOP
         NSTOP = INFO(1)
      END IF
      IF(NSTOP.NE.NODE) CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,
     +               MASK,LS,XLS,NLVL,LWIDTH,NVARS)
      INFO(1) = NSTRT
      INFO(2) = NSTOP
      INFO(3) = MAXDEP
      INFO(4) = LWIDTH
      INFO(5) = NVARS
      INFO(6) = NODES
      END
      SUBROUTINE MC60JD(NSUP,LIRN,NODES,NSTRT,LSTNUM,IRN,ICPTR,VARS,
     +                  STATUS,WEIGHT,NLIST,QUEUE,DEG,PRIOR)
      INTEGER LIRN,LSTNUM,NSUP,NODES,NSTRT
      INTEGER DEG(NSUP),ICPTR(NSUP+1),IRN(LIRN),NLIST(NSUP),
     +        QUEUE(0:NODES-1),STATUS(NSUP),VARS(NSUP)
      DOUBLE PRECISION WEIGHT(2),PRIOR(NSUP)
      INTEGER ADDRES,DEGREE,FATHER,I,ISTOP,ISTRT,J,JSTOP,JSTRT,J1,J2,
     +        K,L,NABOR,NBR,NEXT,NODE,NQ,QNODE,SON,THRESH
      PARAMETER (THRESH=100)
      DOUBLE PRECISION MAXPRT,PNODE,PRTY
      DO 10 I = 1,NODES
         NODE = NLIST(I)
         DEGREE = VARS(NODE)
         K = DEGREE
         VARS(NODE) = 0
         DO 7 J = ICPTR(NODE),ICPTR(NODE+1) - 1
            DEGREE = DEGREE + VARS(IRN(J))
    7    CONTINUE
         VARS(NODE) = K
         PRIOR(NODE) = -WEIGHT(1)*DEGREE-WEIGHT(2)*STATUS(NODE)
         STATUS(NODE) = 2
         DEG(NODE) = DEGREE
   10 CONTINUE
      NQ = 1
      QUEUE(NQ) = NSTRT
      QUEUE(0) = NSTRT
      STATUS(NSTRT) = 1
      PRIOR(NSTRT) = 1.0E30
      DO 70 L = 1, NODES
         IF(NQ.GT.THRESH)GO TO 100
         ADDRES = 1
         MAXPRT = PRIOR(QUEUE(1))
         DO 30 I = 2,NQ
            PRTY = PRIOR(QUEUE(I))
            IF (PRTY.GT.MAXPRT) THEN
               ADDRES = I
               MAXPRT = PRTY
            END IF
   30    CONTINUE
         NEXT = QUEUE(ADDRES)
         QUEUE(ADDRES) = QUEUE(NQ)
         NQ = NQ - 1
         ISTRT = ICPTR(NEXT)
         ISTOP = ICPTR(NEXT+1) - 1
         IF (STATUS(NEXT).EQ.1) THEN
            DO 40 I = ISTRT,ISTOP
               NBR = IRN(I)
               PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NEXT)
               DEG(NBR) = DEG(NBR) - VARS(NEXT)
               IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
               IF (STATUS(NBR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NBR
                  STATUS(NBR) = 1
                  PRIOR(NBR) = PRIOR(NBR)
               END IF
   40       CONTINUE
         END IF
         LSTNUM = LSTNUM + 1
         STATUS(NEXT) = -LSTNUM
         DO 60 I = ISTRT,ISTOP
            NBR = IRN(I)
            IF (STATUS(NBR).NE.1) GO TO 60
            PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NBR)
            STATUS(NBR) = -1
            DEG(NBR) = DEG(NBR) - VARS(NBR)
            IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
            JSTRT = ICPTR(NBR)
            JSTOP = ICPTR(NBR+1) - 1
            DO 50 J = JSTRT,JSTOP
               NABOR = IRN(J)
               IF (STATUS(NABOR).LT.0) GO TO 50
               PRIOR(NABOR) = PRIOR(NABOR) + WEIGHT(1)*VARS(NBR)
               DEG(NABOR) = DEG(NABOR) - VARS(NBR)
               IF(DEG(NABOR).EQ.0) PRIOR(NABOR) = 1.0E30
               IF (STATUS(NABOR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NABOR
                  STATUS(NABOR) = 1
               END IF
   50       CONTINUE
            STATUS(NBR) = 0
   60    CONTINUE
   70 CONTINUE
      RETURN
  100  DO 120 I = 1, NQ
         NBR = QUEUE(I)
         PNODE = PRIOR(NBR)
         J1 = I
         DO 116 K = 1, J1
            J2 = J1/2
            FATHER = QUEUE(J2)
            IF(PRIOR(FATHER).GE.PNODE) GO TO 118
            QUEUE(J1) = FATHER
            NLIST(FATHER) = J1
            J1 = J2
  116   CONTINUE
  118   QUEUE(J1) = NBR
        NLIST(NBR) = J1
  120 CONTINUE
      I = L
      DO 170 L =I, NODES
         NEXT = QUEUE(1)
         QNODE = QUEUE(NQ)
         PNODE = PRIOR(QNODE)
         NQ = NQ - 1
         J = 2
         J2 = 1
         IF(NQ.GT.1) QUEUE(NQ+1) = QUEUE(NQ)
         DO 125 I = 2, NQ
            IF(J.GT.NQ) GO TO 130
            IF( PRIOR(QUEUE(J)).LT.PRIOR(QUEUE(J+1)) ) J=J+1
            SON = QUEUE(J)
            IF(PNODE.GE.PRIOR(SON)) GO TO 130
            QUEUE(J2) = SON
            NLIST(SON) = J2
            J2 = J
            J = J*2
  125    CONTINUE
  130    QUEUE(J2) = QNODE
         NLIST(QNODE) = J2
         ISTRT = ICPTR(NEXT)
         ISTOP = ICPTR(NEXT+1) - 1
         IF (STATUS(NEXT).EQ.1) THEN
            DO 140 I = ISTRT,ISTOP
               NBR = IRN(I)
               IF (NBR.EQ.NEXT) GO TO 140
               PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NEXT)
               DEG(NBR) = DEG(NBR) - VARS(NEXT)
               IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
               IF (STATUS(NBR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NBR
                  STATUS(NBR) = 1
                  NLIST(NBR) = NQ
               END IF
               PNODE = PRIOR(NBR)
               J = NLIST(NBR)
               DO 133 K = 1, NQ
                  J2 = J/2
                  FATHER = QUEUE(J2)
                  IF(PRIOR(FATHER).GE.PNODE) GO TO 137
                  QUEUE(J) = FATHER
                  NLIST(FATHER) = J
                  J = J2
  133          CONTINUE
  137          QUEUE(J) = NBR
               NLIST(NBR) = J
  140       CONTINUE
         END IF
         LSTNUM = LSTNUM + 1
         STATUS(NEXT) = -LSTNUM
         DO 160 I = ISTRT,ISTOP
            NBR = IRN(I)
            IF (STATUS(NBR).NE.1) GO TO 160
            PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NBR)
            STATUS(NBR) = -1
            DEG(NBR) = DEG(NBR) - VARS(NBR)
            IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
            PNODE = PRIOR(NBR)
            J = NLIST(NBR)
            DO 142 K = 1, NQ
               J2 = J/2
               FATHER = QUEUE(J2)
               IF(PRIOR(FATHER).GE.PNODE) GO TO 144
               QUEUE(J) = FATHER
               NLIST(FATHER) = J
               J = J2
  142       CONTINUE
  144       QUEUE(J) = NBR
            NLIST(NBR) = J
            JSTRT = ICPTR(NBR)
            JSTOP = ICPTR(NBR+1) - 1
            DO 150 J = JSTRT,JSTOP
               NABOR = IRN(J)
               IF (STATUS(NABOR).LT.0) GO TO 150
               PRIOR(NABOR) = PRIOR(NABOR) + WEIGHT(1)*VARS(NBR)
               DEG(NABOR) = DEG(NABOR) - VARS(NBR)
               IF(DEG(NABOR).EQ.0) PRIOR(NABOR) = 1.0E30
               IF (STATUS(NABOR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NABOR
                  STATUS(NABOR) = 1
                  NLIST(NABOR) = NQ
               END IF
               PNODE = PRIOR(NABOR)
               J1 = NLIST(NABOR)
               J2 = J1/2
               FATHER = QUEUE(J2)
               IF(PRIOR(FATHER).GE.PNODE) GO TO 148
               QUEUE(J1) = FATHER
               NLIST(FATHER) = J1
               J1 = J2
               DO 146 K = 2, NQ
                  J2 = J1/2
                  FATHER = QUEUE(J2)
                  IF(PRIOR(FATHER).GE.PNODE) GO TO 148
                  QUEUE(J1) = FATHER
                  NLIST(FATHER) = J1
                  J1 = J2
  146          CONTINUE
  148          QUEUE(J1) = NABOR
               NLIST(NABOR) = J1
  150       CONTINUE
            STATUS(NBR) = 0
  160    CONTINUE
  170 CONTINUE
      END
      SUBROUTINE MC60LD(ROOT,MAXWID,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,
     +                 XLS,NLVL,LWIDTH,NVARS)
      INTEGER LWIDTH,MAXWID,NSUP,NLVL,LIRN,ROOT,NVARS
      INTEGER ICPTR(NSUP+1),IRN(LIRN),LS(NSUP),MASK(NSUP),
     *        VARS(NSUP),XLS(NSUP+1)
      INTEGER I,J,LBEGIN,LNBR,LVLEND,LW,NBR,NODE
      INTRINSIC MAX
      MASK(ROOT) = -MASK(ROOT)
      LS(1) = ROOT
      LVLEND = 0
      NVARS = 0
      LNBR = 1
      LWIDTH = VARS(ROOT)
      DO 35 NLVL = 1,NSUP
          LBEGIN = LVLEND + 1
          LVLEND = LNBR
          XLS(NLVL) = LBEGIN
          LW = 0
          DO 30 I = LBEGIN,LVLEND
              NODE = LS(I)
              DO 20 J = ICPTR(NODE), ICPTR(NODE+1) - 1
                  NBR = IRN(J)
                  IF (MASK(NBR).GT.0) THEN
                      LNBR = LNBR + 1
                      LS(LNBR) = NBR
                      MASK(NBR) = -MASK(NBR)
                      LW = LW + VARS(NBR)
                  END IF
   20         CONTINUE
   30     CONTINUE
          LWIDTH = MAX(LW, LWIDTH)
          NVARS = NVARS + LW
          IF (LNBR.EQ.LVLEND) GO TO 40
          IF (LWIDTH.GE.MAXWID) GO TO 40
   35 CONTINUE
   40 XLS(NLVL+1) = LVLEND + 1
      DO 50 I = 1,LNBR
          MASK(LS(I)) = ABS(MASK(LS(I)))
   50 CONTINUE
      END
      SUBROUTINE MC60OD(N,NC,LIRN,IRN,ICPTR,SVAR,NSUP,NEW,VARS,FLAG)
      INTEGER N
      INTEGER NC
      INTEGER LIRN
      INTEGER IRN(LIRN)
      INTEGER ICPTR(NC+1)
      INTEGER SVAR(N)
      INTEGER NSUP
      INTEGER NEW(N)
      INTEGER VARS(N)
      INTEGER FLAG(N)
      INTEGER I,IS,J,JS,K,K1,K2
      DO 10 I = 1,N
         SVAR(I) = 1
   10 CONTINUE
      VARS(1) = N
      FLAG(1) = 0
      NSUP = 1
      DO 40 J = 1,NC
         K1 = ICPTR(J)
         K2 = ICPTR(J+1) - 1
         DO 20 K = K1,K2
            IS = SVAR(IRN(K))
            VARS(IS) = VARS(IS) - 1
   20    CONTINUE
         DO 30 K = K1,K2
            I = IRN(K)
            IS = SVAR(I)
            IF (FLAG(IS).LT.J) THEN
               FLAG(IS) = J
               IF (VARS(IS).GT.0) THEN
                  NSUP = NSUP + 1
                  VARS(NSUP) = 1
                  FLAG(NSUP) = J
                  NEW(IS) = NSUP
                  SVAR(I) = NSUP
               ELSE
                  VARS(IS) = 1
                  NEW(IS) = IS
                  SVAR(I) = IS
               END IF
            ELSE
               JS = NEW(IS)
               VARS(JS) = VARS(JS) + 1
               SVAR(I) = JS
            END IF
   30    CONTINUE
   40 CONTINUE
      END
      SUBROUTINE MC60PD(N,NSUP,LIRN,IRN,ICPTR,SVAR,VARS,VAR,FLAG)
      INTEGER N
      INTEGER LIRN
      INTEGER IRN(LIRN)
      INTEGER NSUP
      INTEGER ICPTR(N+1)
      INTEGER SVAR(N)
      INTEGER VARS(NSUP)
      INTEGER VAR(NSUP)
      INTEGER FLAG(NSUP)
      INTEGER I,IS,J,JS,K,K1,L
      DO 10 IS = 1,NSUP
         FLAG(IS) = -1
         VARS(IS) = 1
   10 CONTINUE
      L = 1
      DO 20 I = 1,N
         IS = SVAR(I)
         JS = FLAG(IS)
         IF(JS.GT.0)THEN
            SVAR(I) = JS
            VARS(JS) = VARS(JS) + 1
         ELSE IF(JS.LT.0)THEN
            FLAG(IS) = L
            VAR(L) = I
            SVAR(I) = L
            L = L + 1
         END IF
   20 CONTINUE
      DO 30 IS = 1,NSUP
         FLAG(IS) = 0
   30 CONTINUE
      L = 1
      K1 = 1
      DO 60 JS = 1,NSUP
         J = VAR(JS)
         K1 = ICPTR(J)
         ICPTR(JS) = L
         DO 50 K = K1, ICPTR(J+1)-1
            IS = SVAR(IRN(K))
            IF(FLAG(IS).NE.JS)THEN
               FLAG(IS) = JS
               IRN(L) = IS
               L = L+1
            END IF
   50    CONTINUE
   60 CONTINUE
      ICPTR(JS) = L
      END

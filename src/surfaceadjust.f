      PROGRAM JACOBI
      IMPLICIT NONE
C---------------------------------------------------------------------
C
C     Jacobi Coordinates for colinear collisions, PES transformation
C     Reaction path definition
C
C   set 1  
C     
C     A ---- B - - - - - - C
C     ___xc___
C        __________X_______ 
C
C   set 2
C     A ---- B - - - - - - C
C     ___R__________
C            ______rc_______ 
C
C
C     rc= X - maxc/(ma+mb)
C     R= xc + mcrc/(mb+mc)  
C     Obs: Atomic units for axis
C          isotopic masses
C
C-------------------------------------------------------------------

      INTEGER I,J,K,NX,NY,KXORD, KYORD,NXG,NYG, LDZ
      PARAMETER (LDZ = 1000)
      REAL*8 POT(LDZ,LDZ), R(LDZ),rc(LDZ),X(LDZ),xc(LDZ)
      REAL*8 MA,MB,MC,Ri, Rf, stepx, Xn(LDZ), BCOEF(LDZ,LDZ),stepy
      REAL*8 XKNOT(LDZ), YKNOT(LDZ), POTN(LDZ,LDZ)
      REAL*8 FATE, FATL, DXTOOL, DYTOOL, rci, rcf
      CHARACTER*30 FIL, OUT, TYPE
      CHARACTER*3 COM
      CHARACTER*1 SPL
      FATL = +5.291772083D-1
      FATE = 1.0D+0

      WRITE(*,*) "enter PES file:"
      READ(*,*) FIL
      OPEN(10,STATUS='OLD',FILE=FIL)

      READ(10,*) COM, NX,NY,MA,MB,MC,OUT, type
c      IF (SPL.EQ.'Y' .OR.SPL.EQ.'y' ) THEN
c         WRITE(*,*) SPL
c         READ(10,*) COM, NXG, NYG
c         WRITE(*,*) "SPLINE PARAMETERS"
c         WRITE(*,*) NX, NXG, NY, NYG
c      ELSE IF (SPL.EQ.'N' .OR.SPL.EQ.'n') THEN
c         WRITE(*,*) "NO SPLINE REQUIRED"
         NXG=NX
         NYG=NY
c      ENDIF
      print*, type
      OPEN(20,FILE=OUT)
      KXORD=3.0D+0
      KYORD=3.0D+0
C---------------------------------------------------------------------
C     PES READING
C     NX is related with X and NY is related to xc

      DO I=1, NX
         DO J=1, NY
            READ(10,*) X(I),xc(J),POT(I,J)
         ENDDO
      ENDDO
C---------------------------------------------------------------------
C Grid
c      Ri=X(1) -  ((xc(1) * MB )/(MA+MB))
c      Rf=X(NX) -  ((xc(NY) * MB )/(MA+MB))
c      rci=X(1) - (ma*xc(1))/(ma+mb)
c      rcf=X(NX) - (ma*xc(NY))/(ma+mb)
c      Ri=xc(1) + (mc*rci)/(mb+mc)
c      Rf=xc(NY) + (mc*rcf)/(mb+mc)

      rci=X(1) - (ma*xc(1))/(ma+mb)
      rcf=X(NX) - (ma*xc(1))/(ma+mb)
c      rcf=X(NX) - (ma*xc(NY))/(ma+mb)
cc      rci=0.56026578331654875
cc      rcf=18.639218604832159 
cc      Ri=2.0196093024160793
cc      Rf=20.280132891658276
      Ri=xc(1) + (mc*rci)/(mb+mc)
      Rf=xc(1) + (mc*rcf)/(mb+mc)

      rci=rci+0.1
      rcf=rcf-0.1
      Ri=Ri+0.1
      Rf=Rf-0.1

      stepx= (Rf - Ri)/(NXG - 1)    
      stepy=( rcf- rci )/(NYG-1)

      IF(NX.LT.KXORD)KXORD = NX
      IF(NY.LT.KYORD)KYORD = NY
      
      WRITE(*,*) "TRANSFORMATION PARAMETERS"
      WRITE(*,*) MA, MB, MC
      WRITE(*,*) "GRID PARAMETERS"
c      WRITE(*,*) 'Ri','Rf', 'stepR', 'ri', 'rf', 'stepr'
      WRITE(*,*) Ri, Rf, stepx, rci,rcf, stepy
c      read(*,*) 
      CALL DBSNAK(NX, X, KXORD, XKNOT)
      CALL DBSNAK(NY, xc, KYORD, YKNOT)
      CALL DBS2IN(NX, X, NY, xc, POT, LDZ , KXORD, KYORD, 
     &        XKNOT, YKNOT, BCOEF)
CCCCCCCCCCCCCCCCCCCCCCCCCCC DEBUG
      WRITE(1,*) 'I','J','R','r','R`','r`'
        DO I=1, NYG
         rc(I)=rci + (I-1)*stepy
         DO J=1, NXG
c            write(*,*) I,J
            R(J)= Ri + (J-1)*stepx
            xc(I)=R(J) - (mc*rc(I))/(mb+mc)
            Xn(J)=rc(I) +  ((xc(I) * Ma )/(MA+MB))
cc            rc(I)=X(J) - (ma*xc(I))/(ma+mb)
cc            R(J)=xc(I) + (mc*rc(I))/(mb+mc)
            write(1,*) I,J, X(J),xc(I), R(J),rc(I)
         ENDDO
         write(1,*)
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCC DEBUG      

      DO I=1, NYG
         rc(I)=rci + (I-1)*stepy
         DO J=1, NXG
c            write(*,*) I,J
            R(J)= Ri + (J-1)*stepx
            xc(I)=R(J) - (mc*rc(I))/(mb+mc)
            Xn(J)=rc(I) +  ((xc(I) * Ma )/(MA+MB))
c
c     CALL DBS2GD (0, 0, 1, XVEC, 1, YVEC, KXORD, KYORD, 
c     &              XKNOT, YKNOT, NX, NY, BCOEF, XYDATA, LDZ)
            write(*,*) i,j
c
            CALL DBS2GD (0, 0, 1, Xn(J), 1, rc(I), KXORD, KYORD, 
     &              XKNOT, YKNOT, NX, NY, BCOEF, POTN(I,J), LDZ)
            WRITE(20,*) rc(I),R(J),POTN(I,J)
c            k=(I-1)*NXG + J
c            HG(k)=POTN(I,J)
c            print*, "k", k, HG(k),POTN(I,J)
         ENDDO
         IF ( TYPE.EQ."GNUPLOT") THEN
         WRITE(20,*)
         ENDIF
      ENDDO
c      read(*,*)
      CLOSE(20)
      WRITE(*,*) "SPLINE AND TRANSFORMATION DONE"
C---------------------------------------------------------------------

      STOP
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DERIVATIVE(Z,NX,NY,DER,AXIS,STEP,LDZ)
C
C     X - input - is a vector with the values of the x axis, 
C     Y - input - is a vector with the values of the y axis
C     Z - input - is a matrix with the values of the function
C     NX - input - is the number of points in the x axis
C     NY - input - is the number of points in the y axis
C     AXIS - input - is a char with the axis where the derivative must be calculated
C     DER - output - is a matrix with the derivative      
C
      IMPLICIT NONE
      REAL*8 Z(LDZ,LDZ), STEP
      REAL*8 DER(LDZ,LDZ)
      INTEGER I,J,K,NX,NY, L, LDZ
      CHARACTER*1 AXIS

      
      IF (AXIS .EQ. 'X' .OR. AXIS .EQ. 'x') THEN
      
         WRITE(*,*) "Hello there x"
         WRITE(*,*) AXIS, STEP
 
      DO J=1,NY
         DO I=2, NX-1
            DER(I,J)=(Z(I+1,J) - Z(I-1,J))/(2*STEP)
         ENDDO
      ENDDO
C-----------------------------------------------delete

c      DO J=1,NY
c         DO I=1, NX-1
c            WRITE(1,*) Y(J), X(I), DER(I,J)
c         ENDDO
c         WRITE(1,*)
c      ENDDO
C-------------------------------------------------

      ELSE IF (AXIS .EQ. 'Y' .OR. AXIS .EQ. 'y') THEN
         WRITE(*,*) "Hello there y"
         WRITE(*,*) AXIS, STEP
         DO I=1,NX
            DO J=2, NY-1
               DER(I,J)=(Z(I,J+1) - Z(I,J-1))/(2*STEP)
            ENDDO
         ENDDO

      ENDIF
      
      END

C--------------------------------------------------------------------
C--------------------------------------------------------------------

  

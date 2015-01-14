      subroutine readwholespl(Fil,IFLTH,X,Y,T,Re,Im,nx,ny,nf,xs,ys,ks
     &     ,nxr,nyr,nfr)
      IMPLICIT NONE
c------------------------------------------
c
c Purpose:
c This routine reads a entire wavepacket propagation into matrices, Re(y,x,t)
c and Im(y,x,t), and then generates a 3D spline representation for it
c
c
c     Goiania, 25th of september of 2014
c     Vinicius Vaz da Cruz
c-------------------------------------------
      INTEGER I,J,NX,NY,NF,K,ldz,xs,ys,ks,IFLTH
      INTEGER II,JJ,KK
      INTEGER NXR,NYR,NFR,NNX,NNY
      REAL*8 X(NX), Y(NY), T(NF), work
      REAL*8 checkRe, checkIm,DBS3VL,RMSD,RMSDRE,RMSDIM
      LOGICAL file_exists,file_open
      CHARACTER*14 wwork,wp,wo,pol
      CHARACTER*30 Fil, Filen
c--- 3D SPLINE PARAMETERS
      REAL*8 TKNOT(NFR+4), YKNOT(NYR+4), XKNOT(NXR+4)
      REAL*8 Re(NYR,NXR,NFR),Im(NYR,NXR,NFR)
      REAL*8 BCOEFRe(NYR,NXR,NFR),BCOEFIm(NYR,NXR,NFR)
C----ERROR ESTIMATION VARIABLES
      REAL*8 RECHK(NY,NX),IMCHK(NY,NX),XCHK(NX),YCHK(NY),TCHK
C---- AUXILIARY COUNTERS
      KK=0
      II=0
      JJ=0

      write(*,*) "Reading wavepacket files: "
      DO K = 1, NF, KS
c--------reads Psi(r,t) for a fixed R value.
         KK = KK +1
         CALL FILENAME(Filen,Fil,IFLTH,K)
         INQUIRE(FILE=Filen, EXIST=file_exists)
         IF (file_exists) THEN
            OPEN(unit=10,file=Filen)
         ELSE
            write(*,*) "could not open file ", Filen
            STOP
         ENDIF         

         READ(10,*) wwork, wp, T(KK), pol
         write(*,'(A20,A3,ES17.10,A5)')FILEN,',',T(KK),'fs'
         DO I=1,NX,1
             IF(MOD(I,XS).EQ.0)II = II +1
            DO J=1, NY,1
               IF(MOD(I,XS).EQ.0 .AND. MOD(J,YS).EQ.0 )THEN
                  JJ = JJ +1
                  READ(10,*) X(II), Y(JJ), Re(JJ,II,KK),Im(JJ,II,KK)
c                  write(*,'(3I4,4ES20.5)')II,JJ,KK,
c     &                 X(II), Y(JJ), Re(JJ,II,KK),Im(JJ,II,KK)
               ELSE
                  READ(10,*)
               ENDIF
            ENDDO
            JJ=0
         ENDDO
         II=0
         CLOSE(10)
c--------------------- end of K loop
      ENDDO

      write(*,*)
      write(*,*)'All data has been read!'
      write(*,*)'------------------------------------------------------'
      write(*,*)


c--------spline procedure
      write(*,*)
      write(*,*)'Generating 3D spline coefficients matrix'
      write(*,*)
   
      CALL DBSNAK(NYR, Y, 3, YKNOT)
      CALL DBSNAK(NXR, X, 3, XKNOT)
      CALL DBSNAK(NFR, T, 3, TKNOT)

      DO I =1, NYR,1
         write(1,*) Y(I)
      ENDDO
      write(1,*)
      DO I =1, NXR,1
         write(1,*) X(I)
      ENDDO
      write(1,*)
      DO I =1, NFR,1
         write(1,*) T(I)
      ENDDO

c      write(*,*) 'real part'
      CALL DBS3IN(NYR,Y,NXR,X,NFR,T,Re,NYR,NXR,
     &     3,3,3,YKNOT,XKNOT,TKNOT,BCOEFRe)
c      write(*,*) 'rimaginary part'
      CALL DBS3IN(NYR,Y,NXR,X,NFR,T,Im,NYR,NXR,
     &     3,3,3,YKNOT,XKNOT,TKNOT,BCOEFIm)

      write(*,*)
      write(*,*)'3D spline matrices calculated!'
      write(*,*)
      write(*,*)'estimating error in interpolation:'
      write(*,*)
      CALL ERRORSPL(FIL,IFLTH,Re,Im,BCOEFRE,BCOEFIM,Y,X,T
     &    ,NY,NX,NF,NYR,NXR,NFR,YKNOT,XKNOT,TKNOT)

      write(*,*)
      write(*,*)'Finished 3D spline section!'
      write(*,*)'------------------------------------------------------'
      
      END

      
      c----------- routine to perform the error analysis in the spline interpolation
      SUBROUTINE ERRORSPL(FIL,IFLTH,Re,Im,BCOEFRE,BCOEFIM,Y,X,T
     &    ,NY,NX,NF,NYR,NXR,NFR,YKNOT,XKNOT,TKNOT)
      IMPLICIT NONE
      INTEGER I,J,K,NX,NY,NF,NXR,NYR,NFR,IFLTH
      CHARACTER*30 Filen,Fil
      REAL*8 RMSD,RMSDRE,RMSDIM,ERRMAXRE,ERRMAXIM
      REAL*8 Re(NYR,NXR,NFR),Im(NYR,NXR,NFR),BCOEFRE(NYR,NXR,NFR)
      REAL*8 YKNOT(NYR),XKNOT(NXR),TKNOT(NFR)
      REAL*8 BCOEFIM(NYR,NXR,NFR),X(NXR),Y(NYR),T(NFR)
      REAL*8 XCHK(NX),YCHK(NY),TCHK,RECHK(NY,NX),IMCHK(NY,NX)



      write(*,*)'computing RMSE for included points'
      write(*,'(A38, A15,A3,A15)')'','Real part','','Imag. part'
      RMSDRE=RMSD(Re(1:NYR,1:NXR,1),BCOEFRE,Y,X,T(1),YKNOT,XKNOT,TKNOT,
     &     NYR,NXR,NYR,NXR,NFR,ERRMAXRE)
      RMSDIM=RMSD(Im(1:NYR,1:NXR,1),BCOEFIM,Y,X,T(1),YKNOT,XKNOT,TKNOT,
     &     NYR,NXR,NYR,NXR,NFR,ERRMAXIM)
      write(*,'(A35,I4,X2 2ES17.10)')'Root mean square errors for file',
     &     1,RMSDRE,RMSDIM
      write(*,'(A35,I4,X2 2ES17.10)')'maximum error for file',
     &     1,ERRMAXRE,ERRMAXIM
      
      RMSDRE=RMSD(Re(1:NYR,1:NXR,NFR-1),BCOEFRE,Y,X,T(NFR-1),YKNOT,XKNOT
     &     ,TKNOT,NYR,NXR,NYR,NXR,NFR,ERRMAXRE)
      RMSDIM=RMSD(Im(1:NYR,1:NXR,NFR-1),BCOEFIM,Y,X,T(NFR-1),YKNOT,XKNOT
     &     ,TKNOT,NYR,NXR,NYR,NXR,NFR,ERRMAXIM)
      write(*,'(A35,I4,X2 2ES17.10)')'Root mean square errors for file',
     &     NF-1,RMSDRE,RMSDIM
      write(*,'(A35,I4,X2 2ES17.10)')'maximum error for file',
     &     NF-1,ERRMAXRE,ERRMAXIM

C------------   
      write(*,*)
      write(*,*)'computing RMSE for all points'
      CALL READWHOLE(Fil,IFLTH,1,RECHK,IMCHK,XCHK,YCHK,NX,NY,TCHK)
      RMSDRE=RMSD(RECHK,BCOEFRE,YCHK,XCHK
     &     ,TCHK,YKNOT,XKNOT,TKNOT,NY,NX,NYR,NXR,NFR,ERRMAXRE)
      RMSDIM=RMSD(IMCHK,BCOEFIM,YCHK,XCHK
     &     ,TCHK,YKNOT,XKNOT,TKNOT,NY,NX,NYR,NXR,NFR,ERRMAXIM)
      write(*,'(A35,I4,X2 2ES17.10)')'Root mean square errors for file',
     &     1,RMSDRE,RMSDIM
      write(*,'(A35,I4,X2 2ES17.10)')'maximum error for file',
     &     1,ERRMAXRE,ERRMAXIM

      write(*,*)
      write(*,*)'computing RMSE for a file not included'
      CALL READWHOLE(Fil,IFLTH,NF-3,RECHK,IMCHK,XCHK,YCHK,NX,NY,TCHK)
      RMSDRE=RMSD(RECHK,BCOEFRE,YCHK,XCHK
     &     ,TCHK,YKNOT,XKNOT,TKNOT,NY,NX,NYR,NXR,NFR,ERRMAXRE)
      RMSDIM=RMSD(IMCHK,BCOEFIM,YCHK,XCHK
     &     ,TCHK,YKNOT,XKNOT,TKNOT,NY,NX,NYR,NXR,NFR,ERRMAXIM)
      write(*,'(A35,I4,X2 2ES17.10)')'Root mean square errors for file',
     &     NF-3,RMSDRE,RMSDIM
       write(*,'(A35,I4,X2 2ES17.10)')'maximum error for file',
     &     NF-3,ERRMAXRE,ERRMAXIM

c      RMSDIM=RMSD(Im(1:NYR,1:NXR,1),BCOEFIM,Y,X,T(1),YKNOT,XKNOT,TKNOT,
c     &     NYR,NXR,NYR,NXR,NFR)
c      write(*,'(A35,I4,X2 2ES17.10)')'Root mean square errors for file',
c     &     1,RMSDRE,RMSDIM

      END


c-------------- routine to print two matrices M1(J,I) and M2(I,J) into a file
      SUBROUTINE PRTWP(M1,M2,X,Y,NX,NY,FILEN)
      IMPLICIT NONE
      INTEGER I,J,NX,NY
      CHARACTER*30 FILEN
      REAL*8 M1(NY,NX),M2(NY,NX),X(NX),Y(NY)
      OPEN(unit=12,FILE=FILEN)

      DO I=1,NX,1
         DO J=1,NY,1
            WRITE(12,'(4ES17.10)')X(I),Y(J),M1(J,I),M2(J,I)
         ENDDO
      ENDDO
      CLOSE(12)
      END

c------------------ routine to compute the root mean square between
c------------------ a matrix M1, and the values generated by the
c------------------ spline coefficients matrix
      FUNCTION RMSD(M,MBCOEF,Y,X,T,YKNOT,XKNOT,TKNOT,NY,NX,
     &     NYR,NXR,NFR,ERRMAX)
      IMPLICIT NONE
      INTEGER I,J,NX,NY,NT,KX,KY,KZ
      INTEGER NXR,NYR,NFR
      REAL*8 VAR,RMSD,MSPL,M(NY,NX),ERRMAX
      REAL*8 X(NX),Y(NY),T,XKNOT(NXR),YKNOT(NYR),TKNOT(NFR)
      REAL*8 MBCOEF(NYR,NXR,NFR),DBS3VL,VAL
      
      KX=3.0D+0
      KY=3.0D+0
      KZ=3.0D+0
      NT=NY*NX
      VAR = 0.0D+0
      ERRMAX = 0.0

      DO I=3,NX-3,1
         DO J=3,NY-3,1
            MSPL=dbs3vl(Y(J),X(I),T,ky,kx,kz,
     &           yknot,xknot,tknot,nyr,nxr,nfr,mbcoef)
            VAL = (M(J,I) - MSPL)**2 
            IF(DSQRT(VAL).GT.ERRMAX) ERRMAX=DSQRT(VAL)
            VAR = VAR + VAL
         ENDDO
      ENDDO

      RMSD = DSQRT(VAR/(NT*1.0D+0))

      RETURN
      END

c------------------- routine to determine filename

      SUBROUTINE FILENAME(Filen,Fil,IFLTH,K)
      IMPLICIT NONE
      INTEGER K,IFLTH,LEN_TRIM
      CHARACTER*30 Filen,Fil
      CHARACTER*4 CHNUM

      WRITE(CHNUM,'(I4)')K
      IF(K.LT.10)THEN
         Filen = Fil(1:IFLTH)//'000'//CHNUM(4:4)//'.dat'
      ELSEIF(K.LT.100)THEN
         Filen = Fil(1:IFLTH)//'00'//CHNUM(3:4)//'.dat'
      ELSEIF(K.LT.1000)THEN
         Filen = Fil(1:IFLTH)//'0'//CHNUM(2:4)//'.dat'
      ELSEIF(K.LT.10000)THEN
         Filen = Fil(1:IFLTH)//CHNUM(1:4)//'.dat'
      ENDIF

      END

c------------------- routine to read whole wavepacket file

      SUBROUTINE READWHOLE(Fil,IFLTH,K,Re,Im,X,Y,NX,NY,T)
      IMPLICIT NONE
      INTEGER I,J,K,IFLTH,NX,NY
      REAL*8 Re(NY,NX), Im(NY,NX), X(NX),Y(NY),T
      CHARACTER*30 Fil, FILEN,wwork,wp,pol
      CHARACTER*4 CHNUM

      CALL FILENAME(FILEN,FIL,IFLTH,K)

      OPEN(UNIT=25,FILE=FILEN)
      READ(25,*) wwork, wp, T, pol
      DO I=1,NX,1
         DO J=1,NY,1
            READ(25,*) X(I), Y(J), Re(J,I),Im(J,I)
         ENDDO
      ENDDO
      CLOSE(25)
      END

      subroutine readpot(X,Fil,nx,ny,ldz)
      IMPLICIT NONE
      INTEGER I,J,nx, ny,K,ldz
      REAL*8 X(ldz,*), work
c      REAL*8 X(ny,nx), work
      LOGICAL file_exists
      CHARACTER*13 Fil, wwork, wp, wo
      INQUIRE(FILE=Fil, EXIST=file_exists)
      
c      print*, "inside fortran >", File

      IF (file_exists) THEN
c     debug
c       print*, "fortran routine ", Fil
       OPEN(unit=10,file=Fil)
      ELSE
      write(*,*) "could not open file ", Fil
      ENDIF

      READ(10,*) wwork, wp, wo
c     debug
c      print*, wwork, wp, wo

      DO I=1, nx
         DO J=1, ny
            READ(10,*) work, work, X(J,I)
c            write(1,*) X(J,I)
         ENDDO
      ENDDO

      CLOSE(10)
      END

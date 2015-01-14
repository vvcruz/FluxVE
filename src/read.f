      subroutine read(X,Y,Fil,nx,ny)
      IMPLICIT NONE
      INTEGER I,J,nx, ny,K
      REAL*8 X(ny,nx), Y(ny,nx), work
      LOGICAL file_exists
      CHARACTER*13 Fil, wwork, wp, wo, pol
      INQUIRE(FILE=Fil, EXIST=file_exists)
      
c      print*, "inside fortran >", File

      IF (file_exists) THEN
c     debug
c       print*, "fortran routine ", Fil
       OPEN(unit=10,file=Fil)
      ELSE
      write(*,*) "could not open file", Fil
      ENDIF

      READ(10,*) wwork, wp, wo, pol
c     debug
c      print*, wwork, wp, wo, pol

      DO I=1, nx
         DO J=1, ny
            READ(10,*) work, work, X(J,I),Y(J,I)
         ENDDO
      ENDDO

      CLOSE(10)
      END

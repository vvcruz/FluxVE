      subroutine readspl(Fil,X,Y,nx,ny,ldz)
      IMPLICIT NONE
      INTEGER I,J,nx, ny,K,ldz
      REAL*8 X(ldz,*), Y(ldz,*), work
      LOGICAL file_exists,file_open
      CHARACTER*14 Fil, wwork, wp, wo, pol
      INQUIRE(FILE=Fil, EXIST=file_exists)
c      INQUIRE(UNIT=10, OPENED=file_open)


c      IF (file_open) THEN
c     debugc
c         print*, "1st FILE IS OPEN ", Fil
c         CLOSE(10)
c      ELSE
c         write(*,*) "1st FILE IS CLOSED ", Fil
c      ENDIF

c      print*, "inside fortran >", File
      IF (file_exists) THEN
c     debug
c       print*, "fortran routine ", Fil
         OPEN(unit=10,file=Fil)
      ELSE
         write(*,*) "could not open file", Fil
      ENDIF

c      INQUIRE(UNIT=10, OPENED=file_open)
c      IF (file_open) THEN
c     debug
c         print*, "2nd FILE IS OPEN ", Fil
c      ELSE
c         write(*,*) "2nd FILE IS CLOSED ", Fil
c      ENDIF
      
      READ(10,*) wwork, wp, wo, pol
c     debug
c      print*, wwork, wp, wo, pol

      DO I=1, nx
         DO J=1, ny
            READ(10,*) work, work, X(J,I),Y(J,I)
         ENDDO
      ENDDO
      
      CLOSE(10)
      
c      INQUIRE(UNIT=10, OPENED=file_open)
c      IF (file_open) THEN
c     debug
c         print*, "3rd FILE IS OPEN ", Fil
c      ELSE
c         write(*,*) "3rd FILE IS CLOSED ", Fil
c      ENDIF
      
 
      END

      PROGRAM TEST
      IMPLICIT NONE
      INTEGER I,J,N,offset,FTELL

      OPEN(UNIT=10,FILE="testy.dat")

      DO I =1, 3
         READ(10,*)N
         WRITE(*,*)N
      ENDDO

      offset = FTELL(10)
      !WRITE(*,*) offset
      CLOSE(10)

      OPEN(UNIT=10,FILE="testy.dat")
      CALL FSEEK(10,offset,1)
      DO I =1, 3
         READ(10,*)N
         WRITE(*,*)N
      ENDDO

      STOP
      END

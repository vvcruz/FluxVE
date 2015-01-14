      program test
      real*8 var

      print*, 'type variable'
      read(*,*) var
      var=var/1D-30
 200  format(E15.6)
      write(*,200) var
      end
 

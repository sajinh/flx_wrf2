c-----------------------------------------------------------------------
c To be used, if the non-standard Fortran function erf does not exist on
c your machine
C
C     function sind
C     function cosd         
C
c-----------------------------------------------------------------------
      function sind(theta)
c
c      theta   angle(unit:degree)
c      sind    sin(theta)
c
      data pi /3.1415926536/

      sind = sin(theta*pi/180.0)
      return
      end
c
c----------------------------------------------------------------------
      function cosd(theta)
c
c      theta   angle(unit:degree)
c      cosd    cos(theta)
c

      data pi /3.1415926536/
      cosd = cos(theta*pi/180.0)
      return
      end


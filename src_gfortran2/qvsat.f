c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################

      FUNCTION F_QVSAT( p, t ) 

c     PURPOSE:
c
c     Calculate the saturation specific humidity using enhanced Teten's
c     formula.
c
c     AUTHOR: Yuhe Liu
c     01/08/1998
c
c     MODIFICATION HISTORY:
c
c     INPUT :
c       p        Pressure (Pascal)
c       t        Temperature (K)
c     OUTPUT:
c       f_qvsat  Saturation water vapor specific humidity (kg/kg).
c
c     Variable Declarations.
c
      implicit none
 
      real p         ! Pressure (Pascal)
      real t         ! Temperature (K)
      real f_qvsat   ! Saturation water vapor specific humidity (kg/kg)
      real f_esl,f_esi,fespt
      
      real rd        ! Gas constant for dry air  (m**2/(s**2*K))
      parameter( rd     = 287.0 )
      real rv        ! Gas constant for water vapor  (m**2/(s**2*K)).
      parameter( rv     = 461.0 )
      real rddrv    
      parameter( rddrv  = rd/rv )


C Change by A. Stohl to save computation time:
c      IF ( t.ge.273.15 ) THEN     ! for water
      if ( t.ge.253.15 ) then      ! modification Petra Seibert
                                   ! (supercooled water may be present)
        fespt=f_esl(p,t)
      else
        fespt=f_esi(p,t)
      endif

      f_qvsat = rddrv * fespt / (p-(1.0-rddrv)*fespt)

      RETURN
      END


      FUNCTION F_ESL( p, t )  

      implicit none
 
      real p         ! Pressure (Pascal)
      real t         ! Temperature (K)
      real f_esl     ! Saturation water vapor pressure over liquid water

      real f

c#######################################################################
c
c     Saturation specific humidity parameters used in enhanced Teten's
c     formula. (See A. Buck, JAM 1981)
c
c#######################################################################

      real satfwa, satfwb
      parameter ( satfwa = 1.0007 )
      parameter ( satfwb = 3.46e-8 )  ! for p in Pa

      real satewa, satewb, satewc
      parameter ( satewa = 611.21 )   ! es in Pa
      parameter ( satewb = 17.502 )
      parameter ( satewc = 32.18 )

      real satfia, satfib
      parameter ( satfia = 1.0003 )
      parameter ( satfib = 4.18e-8 )  ! for p in Pa

      real sateia, sateib, sateic
      parameter ( sateia = 611.15 )   ! es in Pa
      parameter ( sateib = 22.452 )
      parameter ( sateic = 0.6 )

      f = satfwa + satfwb * p
      f_esl = f * satewa * exp( satewb*(t-273.15)/(t-satewc) )

      RETURN
      END

      FUNCTION F_ESI( p, t )  

      implicit none
 
      real p         ! Pressure (Pascal)
      real t         ! Temperature (K)
      real f_esi     ! Saturation water vapor pressure over ice (Pa)

      real f

c#######################################################################
c
c     Saturation specific humidity parameters used in enhanced Teten's
c     formula. (See A. Buck, JAM 1981)
c
c#######################################################################
c
      real satfwa, satfwb
      parameter ( satfwa = 1.0007 )
      parameter ( satfwb = 3.46e-8 )  ! for p in Pa

      real satewa, satewb, satewc
      parameter ( satewa = 611.21 )   ! es in Pa
      parameter ( satewb = 17.502 )
      parameter ( satewc = 32.18 )

      real satfia, satfib
      parameter ( satfia = 1.0003 )
      parameter ( satfib = 4.18e-8 )  ! for p in Pa

      real sateia, sateib, sateic
      parameter ( sateia = 611.15 )   ! es in Pa
      parameter ( sateib = 22.452 )
      parameter ( sateic = 0.6 )

      f = satfia + satfib * p
      f_esi = f * sateia * exp( sateib*(t-273.15)/(t-sateic) )

      RETURN
      END

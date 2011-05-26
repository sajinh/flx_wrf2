      subroutine getfields(itime,nstop)
C                            i     o
********************************************************************************
*                                                                              *
*  This subroutine manages the 3 data fields to be kept in memory.             *
*  During the first time step of petterssen it has to be fulfilled that the    *
*  first data field must have |wftime|<itime, i.e. the absolute value of wftime*
*  must be smaller than the absolute value of the current time in [s].         *
*  The other 2 fields are the next in time after the first one.                *
*  Pointers (memind) are used, because otherwise one would have to resort the  *
*  wind fields, which costs a lot of computing time. Here only the pointers are*
*  resorted.                                                                   *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     29 April 1994                                                            *
*                                                                              *
*  Changes, Bernd C. Krueger, Feb. 2001:                                       *
*        Variables tth,qvh,tthn,qvhn (on eta coordinates) in common block.     *
*        Function of nstop extended.                                           *
*                                                                              *
*  Dec 2005, R. Easter -                                                       *
*          When "memtime(2) = itime = wftime(numbwf)", do not read a new file. *
*          This allows the ending date/time of the flexpart run to match       *
*          the date/time of the last met. file.                                *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* lwindinterval [s]    time difference between the two wind fields read in     *
* indj                 indicates the number of the wind field to be read in    *
* indmin               remembers the number of wind fields already treated     *
* memind(2)            pointer, on which place the wind fields are stored      *
* memtime(2) [s]       times of the wind fields, which are kept in memory      *
* itime [s]            current time since start date of trajectory calculation *
* nstop                > 0, if trajectory has to be terminated                 *
* nx,ny,nuvz,nwz       field dimensions in x,y and z direction                 *
* uu(0:nxmax,0:nymax,nuvzmax,2)   wind components in x-direction [m/s]         *
* vv(0:nxmax,0:nymax,nuvzmax,2)   wind components in y-direction [m/s]         *
* ww(0:nxmax,0:nymax,nwzmax,2)    wind components in z-direction [deltaeta/s]  *
* tt(0:nxmax,0:nymax,nuvzmax,2)   temperature [K]                              *
* ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                        *
*                                                                              *
* Constants:                                                                   *
* idiffmax             maximum allowable time difference between 2 wind fields *
*                                                                            *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer indj,indmin,itime,nstop,memaux
      save indmin

      real uuh(0:nxmax-1,0:nymax-1,nuvzmax)
      real vvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real pvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real wwh(0:nxmax-1,0:nymax-1,nwzmax)
      real uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)

      data indmin/1/

 
C Check, if wind fields are available for the current time step
***************************************************************

      nstop=0

      if ((ldirect*wftime(1).gt.ldirect*itime).or.
     +(ldirect*wftime(numbwf).lt.ldirect*itime)) then
        write(*,*) 'FLEXPART WARNING: NO WIND FIELDS ARE AVAILABLE.'
        write(*,*) 'A TRAJECTORY HAS TO BE TERMINATED.'
        nstop=4
        return
      endif


      if ((ldirect*memtime(1).le.ldirect*itime).and.
     +(ldirect*memtime(2).gt.ldirect*itime)) then

C The right wind fields are already in memory -> don't do anything
******************************************************************

        continue

c FLEXPART_WRF - following change allows the ending date/time 
c of the flexpart run to match that of the last met. file
      else if ( (ldirect*memtime(1).lt.ldirect*itime).and.
     +(memtime(2).eq.itime) .and. (wftime(numbwf).eq.itime) ) then

        continue

      else if ((ldirect*memtime(2).le.ldirect*itime).and.
     +(memtime(2).ne.999999999)) then
 

C Current time is after 2nd wind field
C -> Resort wind field pointers, so that current time is between 1st and 2nd
****************************************************************************

        memaux=memind(1)
        memind(1)=memind(2)
        memind(2)=memaux
        memtime(1)=memtime(2)


C Read a new wind field and store it on place memind(2)
*******************************************************

        do 30 indj=indmin,numbwf-1
           if (ldirect*wftime(indj+1).gt.ldirect*itime) then
              call readwind(indj+1,memind(2),uuh,vvh,wwh)
              call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn)
              call calcpar(memind(2),uuh,vvh,pvh)
              call calcpar_nests(memind(2),uuhn,vvhn,pvhn)
              call verttransform(memind(2),uuh,vvh,wwh,pvh)
              call verttransform_nests(memind(2),uuhn,vvhn,wwhn,pvhn)
              memtime(2)=wftime(indj+1)
              nstop = 1
              goto 40
           endif
 30     continue
 40     indmin=indj

      else

C No wind fields, which can be used, are currently in memory 
C -> read both wind fields
************************************************************

         do 50 indj=indmin,numbwf-1
            if ((ldirect*wftime(indj).le.ldirect*itime).and.
     +           (ldirect*wftime(indj+1).gt.ldirect*itime)) then
               memind(1)=1
               call readwind(indj,memind(1),uuh,vvh,wwh)
               call readwind_nests(indj,memind(1),uuhn,vvhn,wwhn)
               call calcpar(memind(1),uuh,vvh,pvh)
               call calcpar_nests(memind(1),uuhn,vvhn,pvhn)
               call verttransform(memind(1),uuh,vvh,wwh,pvh)
               call verttransform_nests(memind(1),uuhn,vvhn,wwhn,pvhn)
               memtime(1)=wftime(indj)
               memind(2)=2
               call readwind(indj+1,memind(2),uuh,vvh,wwh)
               call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn)
               call calcpar(memind(2),uuh,vvh,pvh)
               call calcpar_nests(memind(2),uuhn,vvhn,pvhn)
               call verttransform(memind(2),uuh,vvh,wwh,pvh)
               call verttransform_nests(memind(2),uuhn,vvhn,wwhn,pvhn)
               memtime(2)=wftime(indj+1)
               nstop = 1
               goto 60
            endif
 50      continue
 60      indmin=indj

      endif

      lwindinterv=abs(memtime(2)-memtime(1))

      if (lwindinterv.gt.idiffmax) nstop=3

      end

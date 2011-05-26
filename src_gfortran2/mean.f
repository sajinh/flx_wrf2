      subroutine mean(x,xm,xs,number)        
********************************************************************************
*                                                                              *
*  This subroutine calculates mean and standard deviation of a given element.  *
*                                                                              *
*      AUTHOR: Andreas Stohl, 25 January 1994                                  *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* x(number)           field of input data                                      *
* xm                  mean                                                     *
* xs                  standard deviation                                       *
* number              number of elements of field x                            *
*                                                                              *
* Constants:                                                                   *
* eps                 tiny number                                              *
*                                                                              *
********************************************************************************

      integer number,i
      real x(number),xm,xs,xl,xq,eps,xaux 
      parameter(eps=1.0e-30)

      xl=0.
      xq=0.
      do 10 i=1,number
        xl=xl+x(i)
10      xq=xq+x(i)*x(i)

      xm=xl/float(number)

      xaux=xq-xl*xl/float(number)

      if (xaux.lt.eps) then
        xs=0.
      else
        xs=sqrt(xaux/float(number-1))
      endif

      end

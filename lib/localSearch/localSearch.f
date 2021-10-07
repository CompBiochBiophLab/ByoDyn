c
c     Project: ByoDyn
c
c     Copyright (C) 2008 Adrian L. Garcia-Lomana & Jordi Villa-Freixa
c
c     Author: Adrian L. Garcia-Lomana
c
c     Created: 2005-10-17 by Adrian L. Garcia-Lomana
c
c     This application is free software; you can redistribute it and/or
c     modify it under the terms of the GNU General Public
c     License as published by the Free Software Foundation; either
c     version 2 of the License, or (at your option) any later version.
c
c     This application is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     Library General Public License for more details.
c     
c     You should have received a copy of the GNU General Public
c     License along with this library; if not, write to the Free
c     Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 USA.
c 
c $Id: localSearch.f,v 4.3 2008/12/03 16:37:38 alglomana Exp $
c
c     this is a program to call properly to dn2gb
c     we need to import this routine from python as follows,
c     f2py -c -m localSearch localSearch.f -lport
c     the flag -lport corresponds to the library PORT, available at:
c     rsync -avz netlib.org::netlib/port .
c     once you have the file localSearch.so, you just import it from a python
c     session as a module
c
      subroutine principal(x, b, s)

      implicit none

      include 'auxVar.dc'

      integer i, iv(ddliv), liv, lv, ui(1), p, n
      double precision b(2, XXDIM), v(ddlv), x(XXDIM) 
      double precision s
      external calca

cf2py intent(callback) pcalcr
      
      double precision xx, vv
      common /xxshared/ xx(XXDIM)
      common /vvshared/ vv

c**********************************************************
c*
c* all this must dissapear
c* WARNING: unnecessary function evaluation !!!

      double precision rr
      external pcalcr
      common /rrshared/ rr(RRDIM)
     
      call pcalcr(x, p)
            
c**********************************************************

c
c     determining some parameters
c

      p = XXDIM
      n = RRDIM

      liv = ddliv
      lv = ddlv
      
c
c     default values for iv and v input components (n2gb and n2fb only). 
c     check for dn2gb
c
      
c
c     initializing the control variables iv and v
c

      call divset(1, iv, liv, lv, v)

c
c     setting some values for the control variables
c

      iv(17) = 1000
      iv(18) = 2000
      v(31) = s
c
c     calling the function
c
c    n is the number of experimental points
c    p is the number of parameters to estimate plus init. cond.
c    x is the actual parameters values
c    b is the maximum and minum value for the parameters          

c     
c     no output on screen
c
      iv(21) = 0
      
      call dn2fb(n, p, x, b, calca, iv, liv, lv, v, ui)
           
c     writing the parameters in a common place for the python
c

      do 10 i = 1, p
         
         xx(i) = x(i)

 10   continue
            
c
c     determining the final function value, written in a common place for the python
c

      vv = v(10)

c
c     ending the program
c

      return
      end


      subroutine calca(n, p, x, nf, r)

c
c     this routine computes the residual vector, r = r(x) by calling a python function
c

c
c     ty, lty and nf must be defined

      implicit none

      include 'auxVar.dc'
      
      integer n, p, nf, i
      double precision x(p), r(n)


      double precision rr                                              
      external pcalcr
      common /rrshared/ rr(RRDIM)  
     
      call pcalcr(x, p)

      do 20 i = 1, n

         r(i) = sqrt(rr(i))

 20   continue

      return
      end

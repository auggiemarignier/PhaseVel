      subroutine forward_modelling(clm,avylm,pred)
      
      include 'parameters.h'

      real*8, dimension(MXPATHS,MXLENY) :: avylm,clm,temp 
      real*8, dimension(MXPATHS) :: pred

c     need to figure out how ylamv is ordered
c     and make sure clm is saved in the same order
c------------------------------------------------
c l = |  0 |      1      |   2         |...
c m = |  0 |  0   1   -1 | 0 1 -1 2 -2 |...
c
c if m > 0, index=l*l +2m -1
c if m =< 0, index = l*l + 2|m|
c
c
c------------------------------------------------

      temp = clm*avylm
      pred = sum(temp,dim=2)

      end subroutine
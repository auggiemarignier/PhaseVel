
      parameter (MXLH=200)
      parameter (MXLENY=(MXLH+1)**2)

      dimension ylmh(MXLENY),ylmt(MXLENY)
      double precision d((2*MXLH+1)**2)
      dimension sar(MXLH+1),wk1(MXLH+1),wk2(MXLH+1),wk3(MXLH+1)

      dimension  dum3(3)
      character*80 getunx,id


      call chekcl('| :r:2:Input data file and output A file'
     1          //'|-lmax:r:1:lmax'
     1          //'|')


      lmaxh=inunx('-lmax',1,llmax)
      ldummy=lmaxh
          
      if(lmaxh.gt.MXLH) then
        write(6,*) 'phasea: increase dimensions'
        call exit(2)
      endif



      open(1,file=getunx(' ',1,ll),status='old')
    !   open(2,file=getunx(' ',2,ll),form='unformatted')
      open(2,file=getunx(' ',2,ll),form='formatted')
      
      npath=0

      lenyh=(lmaxh+1)**2

      write(2,*) lmaxh,ldummy


   10 read(1,*,end=99) xlat1,xlon1,xlat2,xlon2,dc,wt,mima,nsim

      write(*,*)'test', xlat1,xlon1,xlat2,xlon2,dc,wt,mima,nsim

      if(err.eq.0) err=1.
      call delaz(xlat1,xlon1,xlat2,xlon2,del,azep,azst)

      npath=npath+1

      call ylmav(xlat1,xlon1,xlat2,xlon2,lmaxh,ylmh,ylmt
     1   ,wk1,wk2,wk3,sar,d)

c minor arcs only
      if(mima.eq.1) then
         write(2,*) (ylmt(i),i=1,lenyh)
      else
        ch=360./(360.-del)
        ct=-del/(360.-del)
        write(2) xlat1,xlon1,xlat2,xlon2,mima,wt,nsim
     1     ,(ylmh(i)*ch+ylmt(i)*ct,i=1,lenyh)

      endif

      if(mod(npath,100).eq.0) write(6,'(i8,'' paths read'')') npath
      goto 10
   99 continue
      end
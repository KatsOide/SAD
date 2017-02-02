      subroutine initbl
      use maccbk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACVAR.inc'
      include 'inc/MACPHYS.inc'
      include 'inc/MACFILE.inc'
c
      integer*4 i,idummy,sethtb,hsrch,mtaloc,sethtbr
c     external doline
      external doprin, doexpn, doread, dolist, docod, dostop, dotwis
      external dooffl, doonfl,dorvrs
c      external ActLie,ActTra,ActPlt,ActGRA,ActLin
      external ActLie,ActTra,ActPlt,ActGRA
c
c     external dotemp
c
      do i=1,HTMAX
        pname(i)=NULSTR
      enddo
      msglvl=8
c
       idummy=sethtb('LINE    ',icDEF,icLINE)
       idummy=sethtb('line    ',icDEF,icLINE)
       idummy=sethtb('CELL    ',icDEF,icCELL)
       idummy=sethtb('cell    ',icDEF,icCELL)
       idummy=sethtbr('EXPAND  ',doexpn)
       idummy=sethtbr('expand  ',doexpn)
       idummy=sethtbr('print   ',doprin)
       idummy=sethtbr('PRINT   ',doprin)
       idummy=sethtbr('READ    ',doread)
       idummy=sethtbr('read    ',doread)
       idummy=sethtbr('list    ',dolist)
       idummy=sethtbr('LIST    ',dolist)
c       idummy=sethtbr('cod     ',docod)
c       idummy=sethtbr('COD     ',docod)
       idummy=sethtbr('reverse ',dorvrs)
       idummy=sethtbr('REVERSE ',dorvrs)
       idummy=sethtbr('stop    ',dostop)
       idummy=sethtbr('STOP    ',dostop)
c      idummy=sethtbr('temp    ',dotemp)
c      idummy=sethtbr('TEMP    ',dotemp)
c       idummy=sethtbr('twiss   ',dotwis)
c       idummy=sethtbr('TWISS   ',dotwis)
       idummy=sethtbr('on      ',doonfl)
       idummy=sethtbr('ON      ',doonfl)
       idummy=sethtbr('off     ',dooffl)
       idummy=sethtbr('OFF     ',dooffl)
c
       idummy=sethtb('USE     ',icVAR,VarPt)
c
       idummy=sethtb('LIE     ',icACT,mtaloc(4))
       idummy=sethtb('lie     ',icACT,idval(idummy))
       ilist(1,idval(idummy))=2
       call setfnp(ilist(1,idval(idummy)+1),ActLie)
       ilist(1,idval(idummy)+2)=hsrch('USE')
c
       idummy=sethtb('OPTION  ',icVAR,VarStr)
       idummy=sethtb('OPTIONS ',icVAR,VarLst+VarStr)
c
c       idummy=sethtb('LINEAR  ',icACT,mtaloc(5))
c       idummy=sethtb('linear  ',icACT,idval(idummy))
c       call setfnp(ilist(1,idval(idummy)+1),ActLin)
c       ilist(1,idval(idummy)+2)=hsrch('USE     ')
c       ilist(1,idval(idummy)+3)=hsrch('OPTION  ')
c       ilist(1,idval(idummy)+4)=hsrch('OPTIONS ')
c       ilist(1,idval(idummy))=4
c
       idummy=sethtb('TURNS   ',icVAR,VarInt)
       idummy=sethtb('NPART   ',icVAR,VarInt)
       idummy=sethtb('XEMI    ',icVAR,VarRl)
       idummy=sethtb('YEMI    ',icVAR,VarRl)
       idummy=sethtb('SIGS    ',icVAR,VarRl)
       idummy=sethtb('SIGE    ',icVAR,VarRl)
       idummy=sethtb('ELEC    ',icVAR,VarLog)
       idummy=sethtb('POSI    ',icVAR,VarLog)
       idummy=sethtb('PROTON  ',icVAR,VarLog)
       idummy=sethtb('RADI    ',icVAR,VarLog)
       idummy=sethtb('SYNC    ',icVAR,VarLog)
       idummy=sethtb('PTYPE   ',icVAR,VarInt)
       idummy=sethtb('SPAN    ',icVAR,VarRl)
       idummy=sethtb('CENTER  ',icVAR,VarRl)
       idummy=sethtb('AFTER   ',icVAR,VarInt)
       idummy=sethtb('AXIS    ',icVAR,VarLst+ VarInt)
       idummy=sethtb('NX      ',icVAR,VarLst+ VarRL)
       idummy=sethtb('NY      ',icVAR,VarLst+ VarRL)
       idummy=sethtb('NZ      ',icVAR,VarLst+ VarRL)
c
       idummy=sethtb('TRACK   ',icACT,mtaloc(21))
       idummy=sethtb('track   ',icACT,idval(idummy))
       ilist(1,idval(idummy))=20
       call setfnp(ilist(1,idval(idummy)+1),ActTra)
       ilist(1,idval(idummy)+2)=hsrch('USE     ')
       ilist(1,idval(idummy)+3)=hsrch('TURNS   ')
       ilist(1,idval(idummy)+4)=hsrch('NPART   ')
       ilist(1,idval(idummy)+5)=hsrch('XEMI    ')
       ilist(1,idval(idummy)+6)=hsrch('YEMI    ')
       ilist(1,idval(idummy)+7)=hsrch('SIGS    ')
       ilist(1,idval(idummy)+8)=hsrch('SIGE    ')
       ilist(1,idval(idummy)+9)=hsrch('ELEC    ')
       ilist(1,idval(idummy)+10)=hsrch('POSI    ')
       ilist(1,idval(idummy)+11)=hsrch('PROTON  ')
       ilist(1,idval(idummy)+12)=hsrch('RADI    ')
       ilist(1,idval(idummy)+13)=hsrch('SPAN    ')
       ilist(1,idval(idummy)+14)=hsrch('CENTER  ')
       ilist(1,idval(idummy)+15)=hsrch('AFTER   ')
       ilist(1,idval(idummy)+16)=hsrch('AXIS    ')
       ilist(1,idval(idummy)+17)=hsrch('SYNC    ')
       ilist(1,idval(idummy)+18)=hsrch('NX      ')
       ilist(1,idval(idummy)+19)=hsrch('NY      ')
       ilist(1,idval(idummy)+20)=hsrch('NZ      ')
c
       idummy=sethtb('FFS     ',icACT,mtaloc(21))
       idummy=sethtb('ffs     ',icACT,idval(idummy))
       ilist(1,idval(idummy))=20
       call setfnp(ilist(1,idval(idummy)+1),ActTra)
       ilist(1,idval(idummy)+2)=hsrch('USE     ')
       ilist(1,idval(idummy)+3)=hsrch('TURNS   ')
       call setdfl(idummy,hsrch('TURNS   '),-1.d0)
       ilist(1,idval(idummy)+4)=hsrch('NPART   ')
       call setdfl(idummy,hsrch('NPART   '),-3.d0)
       ilist(1,idval(idummy)+5)=hsrch('XEMI    ')
       ilist(1,idval(idummy)+6)=hsrch('YEMI    ')
       ilist(1,idval(idummy)+7)=hsrch('SIGS    ')
       ilist(1,idval(idummy)+8)=hsrch('SIGE    ')
       ilist(1,idval(idummy)+9)=hsrch('ELEC    ')
       ilist(1,idval(idummy)+10)=hsrch('POSI    ')
       ilist(1,idval(idummy)+11)=hsrch('PROTON  ')
       ilist(1,idval(idummy)+12)=hsrch('RADI    ')
       ilist(1,idval(idummy)+13)=hsrch('SPAN    ')
       ilist(1,idval(idummy)+14)=hsrch('CENTER  ')
       ilist(1,idval(idummy)+15)=hsrch('AFTER   ')
       ilist(1,idval(idummy)+16)=hsrch('AXIS    ')
       ilist(1,idval(idummy)+17)=hsrch('SYNC    ')
       ilist(1,idval(idummy)+18)=hsrch('NX      ')
       ilist(1,idval(idummy)+19)=hsrch('NY      ')
       ilist(1,idval(idummy)+20)=hsrch('NZ      ')
c
       idummy=sethtb('QUICK   ',icACT,mtaloc(21))
       idummy=sethtb('quick   ',icACT,idval(idummy))
       ilist(1,idval(idummy))=20
       call setfnp(ilist(1,idval(idummy)+1),ActTra)
       ilist(1,idval(idummy)+2)=hsrch('USE     ')
       ilist(1,idval(idummy)+3)=hsrch('TURNS   ')
       call setdfl(idummy,hsrch('TURNS   '),-1.d0)
       ilist(1,idval(idummy)+4)=hsrch('NPART   ')
       call setdfl(idummy,hsrch('NPART   '),-3.d0)
       ilist(1,idval(idummy)+5)=hsrch('XEMI    ')
       ilist(1,idval(idummy)+6)=hsrch('YEMI    ')
       ilist(1,idval(idummy)+7)=hsrch('SIGS    ')
       ilist(1,idval(idummy)+8)=hsrch('SIGE    ')
       ilist(1,idval(idummy)+9)=hsrch('ELEC    ')
       ilist(1,idval(idummy)+10)=hsrch('POSI    ')
       ilist(1,idval(idummy)+11)=hsrch('PROTON  ')
       ilist(1,idval(idummy)+12)=hsrch('RADI    ')
       ilist(1,idval(idummy)+13)=hsrch('SPAN    ')
       ilist(1,idval(idummy)+14)=hsrch('CENTER  ')
       ilist(1,idval(idummy)+15)=hsrch('AFTER   ')
       ilist(1,idval(idummy)+16)=hsrch('AXIS    ')
       ilist(1,idval(idummy)+17)=hsrch('SYNC    ')
       ilist(1,idval(idummy)+18)=hsrch('NX      ')
       ilist(1,idval(idummy)+19)=hsrch('NY      ')
       ilist(1,idval(idummy)+20)=hsrch('NZ      ')
c
       idummy=sethtb('EMIT    ',icACT,mtaloc(21))
       idummy=sethtb('emit    ',icACT,idval(idummy))
       ilist(1,idval(idummy))=20
       call setfnp(ilist(1,idval(idummy)+1),ActTra)
       ilist(1,idval(idummy)+2)=hsrch('USE     ')
       ilist(1,idval(idummy)+3)=hsrch('TURNS   ')
       call setdfl(idummy,hsrch('TURNS   '),0.d0 )
       ilist(1,idval(idummy)+4)=hsrch('NPART   ')
       ilist(1,idval(idummy)+5)=hsrch('XEMI    ')
       ilist(1,idval(idummy)+6)=hsrch('YEMI    ')
       ilist(1,idval(idummy)+7)=hsrch('SIGS    ')
       ilist(1,idval(idummy)+8)=hsrch('SIGE    ')
       ilist(1,idval(idummy)+9)=hsrch('ELEC    ')
       ilist(1,idval(idummy)+10)=hsrch('POSI    ')
       ilist(1,idval(idummy)+11)=hsrch('PROTON  ')
       ilist(1,idval(idummy)+12)=hsrch('RADI    ')
       ilist(1,idval(idummy)+13)=hsrch('SPAN    ')
       ilist(1,idval(idummy)+14)=hsrch('CENTER  ')
       ilist(1,idval(idummy)+15)=hsrch('AFTER   ')
       ilist(1,idval(idummy)+16)=hsrch('AXIS    ')
       ilist(1,idval(idummy)+17)=hsrch('SYNC    ')
       ilist(1,idval(idummy)+18)=hsrch('NX      ')
       ilist(1,idval(idummy)+19)=hsrch('NY      ')
       ilist(1,idval(idummy)+20)=hsrch('NZ      ')
c
       call initb1
       return
       end

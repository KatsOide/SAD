      subroutine initbl
      use maccbk
      use mackw
      use macphys
      use macvar
      use macfile
      implicit none
c
      integer*8 ktcaloc
      integer*4 idummy,idummy1,hsrch,sethtbr
c     external doline
      external doprin, doexpn, doread, dolist, docod, dostop, dotwis
      external dooffl, doonfl,dorvrs
c      external ActLie,ActTra,ActPlt,ActGRA,ActLin
      external ActLie,ActTra,ActPlt,ActGRA
c
c     external dotemp
c
      pname=NULSTR
      lpname=0
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
       idummy=sethtb8('LIE     ',icACT,ktcaloc(4))
       idummy1=sethtb8('lie     ',icACT,idval(idummy))
       ilist(1,idval(idummy)-1)=2
       call setfnp(klist(idval(idummy)+1),ActLie)
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
       idummy=sethtb8('TRACK   ',icACT,ktcaloc(21))
       idummy1=sethtb8('track   ',icACT,idval(idummy))
       ilist(1,idval(idummy)-1)=20
       call setfnp(klist(idval(idummy)+1),ActTra)
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
       idummy=sethtb8('FFS     ',icACT,ktcaloc(21))
       idummy1=sethtb8('ffs     ',icACT,idval(idummy))
       ilist(1,idval(idummy)-1)=20
       call setfnp(klist(idval(idummy)+1),ActTra)
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
       call setdfl(idummy,hsrch('TURNS   '),-1.d0)
c
       idummy=sethtb8('QUICK   ',icACT,ktcaloc(21))
       idummy1=sethtb8('quick   ',icACT,idval(idummy))
       ilist(1,idval(idummy)-1)=20
       call setfnp(klist(idval(idummy)+1),ActTra)
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
       call setdfl(idummy,hsrch('TURNS   '),-1.d0)
       call setdfl(idummy,hsrch('NPART   '),-3.d0)
c
       idummy=sethtb8('EMIT    ',icACT,ktcaloc(21))
       idummy1=sethtb8('emit    ',icACT,idval(idummy))
       ilist(1,idval(idummy)-1)=20
       call setfnp(klist(idval(idummy)+1),ActTra)
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
       call setdfl(idummy,hsrch('TURNS   '),0.d0 )
c
       call initb1
       return
       end

#$header: /sad/cvsroot/oldsad/obj/makefile,v 1.31 1995/11/20 02:16:32 oide exp $
# makefile for sad1 in ./obj by n. yamamoto feb.9 '93
LIB=
OBJ=
LOPT=
OPT=
SLIB=
FOPT=
#
SPATH=../../src
INSTALLDIR=/SAD
#LIB="*.o"
VPATH=../../src

.suffixes: .o .f .c

.f.o:
	(cd $(SPATH);$(MAKE) $@ TARGETOS=$(TARGETOS))

.c.o:
	(cd $(SPATH);$(MAKE) $@ TARGETOS=$(TARGETOS))

.o:
	(cd $(SPATH);$(MAKE) $@ TARGETOS=$(TARGETOS))


OBJECT=$(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) \
       $(OBJAUTO) $(OBJF) $(OBJDEC) $(OBJDEPEND) $(OBJCA) \
       $(OBJRC) $(OBJTK) \
       $(OBJCKA)

OBJCA= bipipe_.o tfefun2.o  csetenv_.o \
       tpause.o  tftmpnam.o  pointer_.o \
       itregexp_.o cxlistfonts.o MAIN.o

OBJTK= cxlistfonts.o tftcl.o \
       tfpyarg.o tfefun3.o csocket.o itfgetclock_.o

OBJRC= tfefun1.o \
    tfsort.o    tfmodule.o  tfloor.o      tfsolvemember.o  itfmaloc.o \
    tfwrite.o   itfaloc.o   tfmemcheck.o  tfcanvasclip.o  tfdot.o \
    tfconvstr.o tfpart.o    tfreplace.o   tftake.o  tfeval1.o \
    tfearray.o  tfeexpr.o   tfsetlist.o   tfeeval.o  tfmap.o  tftable.o \
    itfdepth.o  tfdset.o    tfsameq.o     itfpmat.o tmatch.o \
    gamma.o     tfbessel.o  spkick.o      tsolqu.o tsolque.o tfbeamline.o\
	netSemaphore.o tfshared.o

OBJ0= qcav.o     JNLPRM.o    msolv1.o    prkick.o \
      prkick1.o  errset.o    errsav.o    twsdrw.o    tdrwdt.o    second.o \
      msortn.o   pamoeb.o    tdplt.o     bb.o \
      ttdr.o     undulator.o phsrot.o    matrix.o    itfork.o    fpipe_.o \
      autofg.o   tffswake.o  tfvectorize.o \
      itfmessage.o tfxlistfonts.o        tffindroot.o tfreadbuf.o spkick.o 

OBJ1=ActGRA.o    filbuf.o    pwmatq.o    ActLie.o    \
     msort.o     pwrite.o    ActTra.o    flmgr.o     msort2.o \
     pwrlat.o    tprofo.o    BEDGE.o     frand.o     mstack.o    tput.o \
     BMAG.o      mstack1.o   qbthin.o    tputbs.o    CM24.o      gaus3.o \
     mstack2.o   temp.o      tpwght.o    DM42.o      gaussn.o    mstat.o \
     qchg.o      tqente.o    FMUL2.o     geodrw.o    mstat2.o \
     qcoord.o    termes.o    tqfrie.o    GrCmd.o     mstatp.o  \
     GrLine.o    mstor1.o    tfadjst.o   tqlfre.o    INTGRL.o \
     gettok.o    mstore.o    IPAK.o      mwght.o     IgetGL.o    getvad.o \
     naiseki.o   tfbeam.o    LgetGL.o    nlfit.o     qdrift.o    \
     Lrdnum.o    MADD2.o     openP.o     hsrch.o     \
     MINV2.o     over.o      MMUL2.o     tfdset.o \
     p0bym.o     trackb.o    MMUL4.o     pack.o      qmdiag.o

OBJ2=packpi.o    tfdbun.o    tracke.o    NewGRF.o    packpim.o \
     trad.o      QUADRU.o    padd.o      R0CAL.o \
     palgn.o     tradp.o     R1INV.o     inifil.o    palgn1.o     \
     R2INV.o     palocx.o    initb1.o    pbset.o     SEXTU.o     initbl.o \
     pbump.o     tfgeti.o    SYMPS.o     iscan.o     pbump1.o \
     qthin.o     tfgetr.o    trotg.o     SYMTR2.o    pbumps.o    TCAL.o \
     pcbak.o     tscale.o    TEDGE.o     pclr.o \
     radlos.o    TMATRS.o    pclri.o     rdexpr.o    tsdrad.o \
     TROT.o      pclrp.o     rdkwdl.o    TWTRANS.o   lread.o \
     pcrmat.o    rdterm.o    UNITM2.o    lrflct.o    pcset.o      \
     tfltra.o    UNITM4.o    lstchk.o    pecorg.o    rslvin.o \
     tserad.o    VMUL4.o     lstelp.o    pecorr.o    rvrsln.o \
     tsetg.o     ZEROM2.o    ltype.o     pecorr1.o   savcod.o  \
     ZEROM4.o    matadd.o    ZEROV4.o    matclr.o \
     permes.o    tfrej.o     matcpy.o    petcod.o    setdfl.o    tbrad.o \
     tsmear.o    actPlt.o    matmax.o    petdif.o    setfnp.o

OBJ3=append.o    matmin.o    petdifa.o   sethtb.o    matmul.o    petfil.o \
     atof.o      matrst.o    petget.o    \
     matsqu.o    petmcut.o   tfsize.o    bint.o \
     matsub.o    petput.o    matswp.o \
     petsub.o    skipch.o    tspect.o    car.o       mbmp.o \
     pettime.o   mbmpf.o     pfcoup.o    cons.o \
     mbstr.o     pgauss.o    sols33.o    convrg.o    mbufw.o \
     pgaussj.o   spldrw.o    corfree.o   mc2to4.o    pgclr.o     spline.o \
     tstrad.o    corinit.o   mc4to2.o    phdrw.o     sprexl.o     \
     mcchar.o    phdrwa.o    sprlin.o   \
     mcchar1.o   pinner.o    mccor.o \
     tfzap.o     ttcav.o     mcepst.o    pkill.o \
     ttcave.o    plssvd.o    talign.o    mcfix.o \
     pltelm.o    tgfun.o     mckick.o    pltrng.o    cssetp.o    mclear.o \
     pmicad.o    tgrot.o     dAssgn.o    mclin.o     pmov.o      \
     data2.o     mcmon.o     pmovi.o     defflg.o    mcmon1.o \
     polpar.o

OBJ4=defglb.o    mcnrmc.o    poplst.o    mcrcod.o    ppair.o \
     tinitr.o    doelem.o    mcrcod1.o   pqcell.o    doexpn.o \
     mcrmat.o    prSad.o     tvcorr.o    doflag.o    mcrmda.o \
     prelem.o    title.o     doline.o    mcstr.o     prelm0.o \
     tcav.o      titrpt.o    dolist.o    mcstr1.o    prexln.o \
     tcave.o     tltrm.o     twbuf.o     mdpmax.o    prkwdv.o \
     metaer.o    prline.o    tchge.o     tluma.o     twelm.o     doprin.o \
     tmap.o      doread.o    mfdel.o     prnGlb.o  \
     tmast.o     dorvrs.o    mfdel1.o    prnflg.o   \
     twinil.o    dostop.o    mfdir.o     psolvg.o    tconv.o    synradcl.o\
     twinit.o    dotemp.o    mfnst.o     pstati.o    tconvm.o    drndsr.o \
     mhogal.o    pstati1.o   tcoord.o    drwkwd.o    mhogan.o    pstati2.o \
     tcoorde.o   drwln.o     pstati3.o \
     wexp.o      drwxbx.o    mkplst.o    psub.o      tday.o      tmuld.o \
     wfbin.o     monact.o    psvbk.o     tdcmd.o     wfres.o \
     eigs33.o    monel.o     psvd.o      wiord.o     monqu.o \
     ptimes.o    wjdraw.o    elname.o    ptol.o   \
     errmsg.o    mrecal.o    ptrace.o    tdinit.o    mrecal1.o \
     ptrim.o     toplvl.o    pundo.o     tdjin.o  \
     wsolvi.o    mrmb.o      push.o      tdnrm.o   \
     wsolvs.o    fexec.o     mrqcof.o    putghd.o    tdout.o     wtune.o \
     expln.o     mrqcov.o    putgrc.o    tphplt.o     \
     fdate1.o    mrqmin.o    putsti.o    xsin.o      datetime.o \
     msolb.o     putstr.o    tdrife.o    asinh.o \
     filaux.o    pvert.o     tprmpt.o    tfgetm.o    tfojit.o    tedrawf.o

OBJ5=trcoda.o    nalign.o    nistack.o    mcrmatb.o  ndelw.o     pwrtstr.o \
     pvbump.o    nfgetm.o    preadstr.o   pasex.o    msolvg.o    preadmon.o \
     mbufg.o     pwrtmon.o   petune.o     elnameK.o  yylex.o  \
     bison.simple.o \
     petunes.o   petune1.o   pvbump2.o    pmeas.o    pmbdata.o   pvbump3.o \
     pmbdrw.o    pmbump.o    pvbump1.o    preabuf.o  pmbdata1.o  \
     initdainterface.o       initda.o     tfepicsconstatcb.o

OBJAUTO=         abbrev.o    cputix.o     csinit.o     \
     doACT.o     eval.o      getbuf.o     getwrd.o   icsmrk.o \
     ielm.o      itfgeto.o   \
     pfalloc.o   pgflag.o    pgmast.o     pgrmat.o \
     pgsolvcond.o            qcell.o      qdbend.o     qdcell.o \
     qddrif.o    qdmdia.o    qdquad.o     qdthin.o   qdtwis.o \
     qgettr.o    qins.o      qmult.o      qquad.o    qquads.o   \
     qtent.o     qtwiss.o     tapert.o \
     tbal.o      tbdecoup.o  tbende.o     tbendi.o   tbfrie.o \
     tceigen.o    tcftr.o    tcod.o \
     tcorr.o     tdet.o      tdfun.o      tdgeo.o    tdlat.o \
     teigen.o    temap.o     temat.o      temit.o    temitf.o \
     temits.o    terror.o    tfaprt.o     tfattr.o   tfcalc.o     tfchgv.o \
     tfchro.o    tfcoup.o    tfdapert.o   tfdisp.o   tfemit.o \
     tffamsetup.o tffile.o \
     tffs.o      tffsa.o      tffscalc.o tffsfreefix.o \
     tffsmatch.o tfgeo.o     tfgetv.o     tfif.o     tfinit.o     tfinitn.o \
     tfkwrd.o    tflag.o     tflogi.o     tfltr1.o   \
     tfmat.o     tfoptics.o  tfprint.o   \
     tfsave.o    tfsetparam.o             tfsetv.o   tfshow.o    tftmat.o \
     tftrack.o   tftrad.o     tftrak.o    tftrb.o \
     tftyp1.o    tftype.o    tfvars.o     tfwake.o   tgauss.o \
     tgetfv.o    thess.o     tiniampl.o   tins.o     tinse.o \
     tintrb.o    tlum.o      tmov.o       tmulbs.o   tmulte.o \
     tmulti.o    tmulta.o    tdrift.o    tftwiss.o    trackd.o \
     tmultr.o    tnorm.o     tpara.o      tqr.o \
     tqrad.o     tquade.o    tquads.o     tquase.o   track.o \
     tracka.o    trade.o     tsconv.o     tsgeo.o \
     tshow.o     tsol.o      tsole.o      tsolvg.o \
     tsolvm.o    tspac.o     tspini.o     tsteee.o   tsteer.o     tsvdm.o \
     tthine.o    ttinit.o    ttstat.o     tcsvdm.o \
     tturne.o    txcalc.o    wjfit.o      wtunem.o     \
     tturn.o     tbend.o     tquad.o      tfsocket.o  tflifetrack.o

OBJF=tfeval.o    tfetok.o    itfcopy.o tfefun.o   italoc.o \
     tfmap.o     tfprinta.o

OBJCKA = CkaDebugArgList.o CkaEllipticRd.o

all:$(OBJECT) lib exe

test:
	echo $(FOO)

exe:$(SADLIB) MAIN.o
	($(FC) MAIN.o $(SLIB) $(OPT) $(FOPT) $(EXEOPT) $(LIB) -o $(OUT)  $(LOPT))

source:$(OBJECT:.o=.f)

$(SADLIB):$(OBJECT)
	ar rv $(SADLIB) $(?)

lib:$(SADLIB)

install:update exe
	cp ./$(OUT) $(INSTALLDIR)/$(TARGETOS)/

update:
	(cd ../;gmake update)


all_object:obj0 obj1 obj2 obj3 obj4 obj5 obj6 objdec objauto objdepend objca

obj0: $(OBJ0)

obj1: $(OBJ1)

obj2: $(OBJ2)

obj3: $(OBJ3)

obj4: $(OBJ4)

obj5: $(OBJ5)

obj6: $(OBJ6)

objdec: $(OBJDEC)

objdepend: $(OBJDEPEND)

objauto: $(OBJAUTO)

objca: $(OBJCA)

objcka: $(OBJCKA)

clean:
	-rm *.o

very_clean:clean
	-rm sadlib.a

debug:
	(cd ..; $(MAKE) debug)


include ../../src/Makefile.depend

# #$Log: Makefile,v $
# #Revision 1.106  2005/11/29 02:10:49  ckallen
# #merging of branch ckallen2 with main trunk
# #
# #Revision 1.106  
# #2005/11/30 ckallen
# # Merging of ckallen2 branch with main SAD branch
# #   Debugging of space charge additions to SAD
# #   Implementation of general matrix functions in SADScript
# #   Implementation of general Twiss parameter utilities in SADScript
# #
# #Revision 1.105  2005/11/07 05:03:18  amorita
# #2005/11/07:
# #  Backport from amorita's subversion repository
# #   Fix proto-type of yyerror()
# #	src/yylex.f		382 - 383
# #
# #   Fix API mismatch of undularor() invoke
# #    Remove unused `iturn' argument
# #    Rename argument name `work' to `pz'
# #	src/undulator.f		383 - 384
# #
# #   Fix API mismatch of mkmachinef() invoke
# #	src/initdainterface.f	384 - 385
# #
# #   Fix API mismatch of termes() invoke
# #	src/tfvars.f		386 - 387
# #
# #   Redirect API call prmmap@malloc.f/prmmap.f(obsolete) to pr_mem_map@pfalloc.f
# #	src/ActLie.f		387 - 388
# #	src/ActTra.f		387 - 388
# #	src/doACT.f		387 - 388
# #	src/dolist.f		387 - 388
# #
# #   Fix API mismatch of doelem() invoke
# #    Remove `dummy' argument
# #    Fix argument type(Integer -> Integer*4)
# #	src/doelem.f		388 - 389
# #
# #   Fix API mismatch of pasex() invoke
# #    mccor() API is changed to obtain pos argument for pasex
# #	src/mccor.f		389 - 390
# #	src/tffsa.f		389 - 390
# #
# #   Fix invalid use of itfmessageexp(enbugged 1.94 - 1.95 on SAD CVS repository)
# #	src/tfeval1.f		390 - 392
# #
# #   Add missing form argument for tfconvstr()
# #	src/mckick.f		392 - 393
# #	src/petune.f		392 - 393
# #	src/preadmon.f		392 - 393
# #	src/preadstr.f		392 - 393
# #	src/pstati.f		392 - 393
# #	src/tffile.f		392 - 393
# #
# #   Remove unused & broken code: msolv()
# #	src/msolv.f		393 - 394(Removed)
# #	obj/Makefile		based on mk/sad.obj.mk(393 - 394)
# #
# #   Fix API mismatch of tthin/tthine/ttstat() invoke
# #	src/tquad.f		394 - 395
# #	src/tquads.f		394 - 395
# #	src/tstrad.f		394 - 395
# #	src/tracka.f		394 - 395
# #	src/tturne.f		394 - 395
# #
# #   Fix API mismatch of qtwiss/qgettr() invoke[FFS:LTRACK command]
# #    enbugged by `dp' argument remove at src/qgettr.f(1.2 - 1.3: 1997/06/05)
# #    enbugged by `dp' argument remove at src/qtwiss.f(1.30 - 1.31: 1997/06/05)
# #    CAUTION: This fix is temporary fix against access violation.
# #             LTRACK command is still broken.
# #	src/tfltr1.f		395 - 397
# #
# #   Fix API mismatch of tftmat() invoke
# #    enbugged since src/qins.f 1.1
# #	src/qins.f		397 - 398
# #
# #   Fix API mismatch of pftmat() invoke
# #    enbugged by `dp' argument remove at src/pasex.f(1.1 - 1.2: 1997/06/05)
# #    CAUTION: This fix is temporary fix against access violation.
# #	src/pasex.f		398 - 399
# #
# #   Remove unlinked source codes
# #	src/tfearray1.f		(Removed)
# #	src/itfaloc1.f		(Removed)
# #	src/malloc.f		(Removed)
# #	src/nfread.f		(Removed)
# #
# #Revision 1.104  2005/09/25 03:25:25  oide
# #V1.0.8.21a
# #
# #Revision 1.103.5  2005/03/30 09:18:05  ckallen
# #SADScript space charge implementation version 1.0
# #
# #Revision 1.103  2003/01/13 02:20:02  oide
# #V1.0.8.17b
# #
# #Revision 1.102  2002/11/04 03:32:00  oide
# #V1.0.8.16a
# #
# #Revision 1.101  2000/09/23 07:39:04  oide
# #V1.0.8.7a
# #
# #Revision 1.100  2000/09/04 23:35:12  oide
# #V1.0.8.6.2b
# #
# #Revision 1.99  2000/08/19 02:48:03  oide
# #V1.0.8.6a
# #
# #Revision 1.98  2000/06/02 01:27:52  yamamoto
# #adding NetSem functions
# #
# #Revision 1.97  1999/11/20 16:03:32  oide
# #V1.0.8.2.5b3
# #
# #Revision 1.96  1999/10/06 02:24:22  oide
# #V1.0.8.1.9a
# #
# #Revision 1.95  1999/10/01 10:14:20  oide
# #V1.0.8.1.7b
# #
# #Revision 1.94  1999/09/30 06:58:08  oide
# #V1.0.8.1.7a
# #
# Revision 1.93  1999/09/29  04:43:58  oide
# V1.0.8.1.6a3
#
# #Revision 1.92  1999/09/28 07:12:42  oide
# #V1.0.8.1.6a1
# #
# #Revision 1.91  1999/09/28 02:12:15  oide
# #V1.0.8.1.6a
# #
# #Revision 1.90  1999/09/05 13:03:09  oide
# #V1.0.8.1.0a8
# #
# #Revision 1.89  1999/09/05 02:40:15  oide
# #V1.0.8.1.0a8
# #
# #Revision 1.88  1999/09/01 22:32:27  oide
# #V1.0.8.1.0a
# #
# #Revision 1.86  1999/08/21 06:32:38  oide
# #V1.0.8.0.1a1
# #
# #Revision 1.85  1999/05/28 15:07:12  oide
# #V1.0.7.12.9b
# #
# #Revision 1.84  1999/05/13 05:28:27  oide
# #V1.0.7.12.4a5
# #
# #Revision 1.83  1999/05/12 02:38:24  oide
# #V1.0.7.12.14a2
# #
# #Revision 1.82  1999/05/11 08:04:05  yamamoto
# #Add oracle library and gcc lib for linking
# #
# #Revision 1.81  1999/05/10 03:01:45  oide
# #V1.0.7.12.5a
# #
# #Revision 1.80  1999/04/30 03:07:44  yamamoto
# #now compatible with bison
# #
# #Revision 1.79  1999/04/30 02:13:21  oide
# #V1.0.7.12.1b
# #
# Revision 1.78  1998/10/10  03:11:20  akasaka
# cumulative changes
#
# #Revision 1.77  1998/09/29 10:53:38  oide
# #V1.0.7.8.6a3
# #
# #Revision 1.76  1998/09/04 14:14:51  oide
# #V1.0.7.7.15a
# #
# #Revision 1.75  1998/06/14 12:18:54  oide
# #V1.0.7.5.9b
# #
# #Revision 1.74  1998/06/12 12:52:52  oide
# #V1.0.7.5.7b
# #
# #Revision 1.73  1998/05/09 03:35:26  oide
# #V1.0.7.3a10
# #
# #Revision 1.72  1998/03/31 03:54:44  oide
# #V1.0.7.1b
# #
# #Revision 1.71  1998/03/05 00:23:30  oide
# #V1.0.6.3.2a
# #
# #Revision 1.70  1998/02/20 06:51:06  akasaka
# #restore omitted objdec
# #
# #Revision 1.69  1998/02/20 06:42:08  oide
# #V1.0.6.2b1
# #
# #Revision 1.68  1998/02/20 05:25:11  oide
# #V1.0.6.2b
# #
# #Revision 1.67  1998/02/17 03:31:39  oide
# #V1.0.6.1.3b
# #
# #Revision 1.66  1998/01/06 01:55:53  oide
# #V1.0.6a
# #
# #Revision 1.65  1997/11/30 09:21:23  oide
# #V1.0.5.9b
# #
# #Revision 1.64  1997/11/05 02:07:58  akasaka
# #*** empty log message ***
# #
# #Revision 1.63  1997/10/13 22:46:40  oide
# #V1.0.5.3b
# #
# #Revision 1.62  1997/10/07 01:33:21  oide
# #V1.0.5a
# #
# #Revision 1.61  1997/10/02 08:03:12  oide
# #V1.0.5a
# #
# #Revision 1.60  1997/09/30 09:57:46  akasaka
# #XListFonts support
# #
# #Revision 1.59  1997/09/25 11:18:04  oide
# #V1.0.4b
# #
# #Revision 1.58  1997/09/04 06:18:49  oide
# #V1.0.3.10.2b
# #
# #Revision 1.57  1997/07/15 09:28:47  oide
# #V1.0.3.9b
# #
# #Revision 1.56  1997/07/08 01:00:10  oide
# #V1.0.3.8b
# #
# #Revision 1.55  1997/06/26 07:28:26  oide
# #V1.0.3.7.2b
# #
# #Revision 1.54  1997/06/05 09:42:26  oide
# #V1.0.3b
# #
# #Revision 1.53  1997/05/14 02:38:51  oide
# #V1.0.2.7b
# #
# #Revision 1.52  1997/05/11 06:42:48  oide
# #V1.0.2.6.1b
# #
# #Revision 1.51  1997/05/10 09:49:22  oide
# #V1.0.2.6b
# #
# #Revision 1.50  1997/04/21 06:57:09  oide
# #V1.0.2.2b
# #
# #Revision 1.49  1997/04/03 13:18:10  yamamoto
# # Python->SAD interface
# #
# #Revision 1.48  1997/04/03 06:48:39  yamamoto
# #Python Interface
# #
# #Revision 1.47  1997/03/05 08:34:10  oide
# #oide
# #
# #Revision 1.46  1997/02/15 08:55:56  yamamoto
# #Improved CA access
# #
# #Revision 1.41  1996/07/28 09:17:25  yamamoto
# # .inc files dependece.
# #
# #Revision 1.40  1996/07/26 10:37:32  yamamoto
# #longer ( 32bytes length) symbol support
# #
# #Revision 1.39  1996/05/29 04:37:29  oide
# #oide
# #
# #Revision 1.38  1996/05/16 01:20:33  oide
# #oide
# #
# #Revision 1.37  1996/05/09 01:36:30  oide
# #oide
# #
# #Revision 1.36  1996/04/22 05:33:57  oide
# #oide
# #
# #Revision 1.35  1996/02/10 13:55:38  oide
# #oide
# #
# #Revision 1.34  1996/01/26 12:16:03  yamamoto
# #use standard pipe instead of named pipe(fifo)
# #
# #Revision 1.33  1996/01/24 02:46:24  oide
# #oide
# #
# #Revision 1.32  1996/01/24 01:10:49  yamamoto
# # working correctoly on HP-UX
# #
# #Revision 1.31  1995/11/20 02:16:32  oide
# #oide
# #
# #Revision 1.30  1995/11/18 10:25:53  oide
# #oide
# #
# #Revision 1.29  1995/11/16 10:01:34  oide
# #oide
# #
# #Revision 1.28  1995/11/16 01:59:12  oide
# #oide
# #
# #Revision 1.27  1995/11/15 06:28:35  oide
# #oide
# #
# #Revision 1.26  1995/11/05 14:59:46  oide
# #oide
# #
# #Revision 1.25  1995/11/02 22:30:58  oide
# #oide
# #
# Revision 1.24  1995/10/20  11:10:39  oide
# oide
#
# Revision 1.23  1995/10/18  07:21:43  oide
# oide
#
# #Revision 1.22  1995/10/18 05:56:41  oide
# #oide
# #
# #Revision 1.21  1995/10/10 02:29:27  oide
# #oide 10/10/1995
# #
# #Revision 1.20  1995/10/08 18:07:02  oide
# #oide 10/9/1995
# #
# #Revision 1.19  1995/10/06 21:57:21  oide
# #oide 10/7/1995
# #
# #Revision 1.18  1995/10/04 22:05:15  oide
# #oide 10/5/1995
# #
# #Revision 1.17  1995/10/04 20:08:35  oide
# #oide 10/5/1995
# #
# #Revision 1.16  1995/10/04 13:41:27  yamamoto
# #keep update
# #
# #Revision 1.15  1995/09/30 16:35:43  yamamoto
# #automatic check in by Make
# #
# #Revision 1.14  1995/09/21 10:06:41  oide
# #oide 9/21 19:01
# #
# #Revision 1.13  1995/09/21 06:00:46  yamamoto
# #Noboru's HP Work
# #
# #Revision 1.12  1995/09/21 05:43:20  yamamoto
# #Noboru HP-Work
# #
# #Revision 1.10  1995/09/19 14:35:44  yamamoto
# #makefile for HP-UX/OSF1
# # #Revision 1.9 1995/09/17 14:52:19 oide #oide
# 9/17 23:52 # #Revision 1.8 1995/09/17 09:09:14 oide #oide 9/17 18:10
# # #Revision 1.7 1995/09/16 09:53:30 oide #oide 9/16 18:53 #
# #Revision 1.6 1995/09/16 01:32:10 oide #oide 9/16 10:31 # #Revision
# 1.5 1995/09/15 08:16:42 yamamoto #working Alpha version 1.0 #
# #Revision 1.4 1995/09/14 23:49:39 oide #oide 9/15 8:49 # #Revision
# 1.3 1995/09/13 15:20:37 yamamoto #tweak Makefiles # ## changed by
# Y.Tange(DEC) Jan.11 '95 FC=f77 CC=cc OUT=sad1.exe SPATH=../../src
# OPATH=../obj/$(TARGETOS) SADLIB=sadlib.a VPATH=../src ../obj
# ../src/inc ../obj/$(TARGETOS) MAKE=gmake


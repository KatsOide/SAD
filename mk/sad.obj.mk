# New Makefile -- Common Object Part
#
all_object:	objmod obj0 obj1 obj2 obj3 obj4 obj5 obj6 \
		objauto objf objrc objfunc \
		objutil objosdep objsim objglue MAIN.o

OBJ_LIBSAD=$(OBJMOD) $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) \
	$(OBJAUTO) $(OBJF) $(OBJRC) \
	$(OBJFUNC) \
	$(OBJUTIL) $(OBJOSDEP)

OBJMOD= tfstk.o		tfstrbuf.o	tfreadbuf.o	toplvl.o	tintrb.o \
	tftwiss.o	tftake.o	tfsolvemember.o	tfetok.o	tfbeamline.o\
	temit.o		tfloor.o	 tffs.o tdjin.o 

OBJRC= 	tfefun1.o	tfsort.o	tfmodule.o	itfmaloc.o \
	tfwrite.o	itfaloc.o	tfmemcheck.o	tfdot.o		tfematrix.o \
	tfconvstr.o	tfpart.o	tfreplace.o	tftake.o	tfeval1.o \
	tfearray.o	tfeexpr.o	tfsetlist.o	tfeeval.o	tfmap.o	tftable.o \
	itfdepth.o	tfdset.o	tfsameq.o	itfpmat.o	tmatch.o \
	gamma.o		tfbessel.o	 spkick.o	tsolqu.o	tsolque.o \
	tfshared.o	tfinitn.o	tffsa.o photons.o

OBJ0= 	italoc.o qcav.o	JNLPRM.o msolv1.o \
	msortn.o tdplt.o bbstrhl.o bbtool.o \
	undulator.o phsrot.o matrix.o autofg.o tffswake.o tfvectorize.o \
	itfmessage.o  tffindroot.o spkick.o \
	autos_.o      gdtoa/g_dfmt.o gdtoa/gdtoa.o gdtoa/g__fmt.o \
	gdtoa/gmisc.o gdtoa/dmisc.o  gdtoa/misc.o

OBJ1=ActGRA.o    filbuf.o    pwmatq.o    ActLie.o    \
     pwrite.o    ActTra.o    flmgr.o   pwrlat.o   frand.o    tput.o \
     tputbs.o    CM24.o        DM42.o      gaussn.o \
     qchg.o      tqente.o    FMUL2.o     geodrw.o  \
     qcoord.o    termes.o    tqfrie.o  tfadjst.o   tqlfre.o  INTGRL.o \
     gettok.o    IPAK.o  IgetGL.o getvad.o \
     tfbeam.o    LgetGL.o         qdrift.o    \
     Lrdnum.o    openP.o     hsrch.o     \
     MINV2.o     over.o      MMUL2.o     tfdset.o \
     trackb.o    pack.o      qmdiag.o

OBJ2=packpi.o    tfdbun.o    tracke.o    NewGRF.o padd.o       \
     R1INV.o     inifil.o	R2INV.o     initb1.o        initbl.o \
     tfgeti.o    SYMPS.o  tfgetr.o   SYMTR2.o  \
     rdexpr.o    tsdrad.o rdkwdl.o    lread.o \
     rdterm.o    UNITM2.o lrflct.o    \
     tfltra.o    UNITM4.o lstchk.o    rslvin.o rvrsln.o \
     ZEROM2.o	ZEROM4.o ZEROV4.o \
     tfrej.o     setdfl.o    tsmear.o    actPlt.o

OBJ3=sethtb.o    atof.o      bint.o      skipch.o    tspect.o    pfcoup.o \
     sols33.o    pgaussj.o   spline.o     sprexl.o    tlinit.o \
     sprlin.o    pinner.o  tfzap.o     ttcav.o      \
     ttcave.o    talign.o    tgfun.o      \
     pmicad.o    tgrot.o     dAssgn.o     defflg.o

OBJ4=defglb.o     \
     tinitr.o    doelem.o    doexpn.o \
     prSad.o     tvcorr.o    doflag.o     \
     prelem.o    title.o     doline.o    prelm0.o \
     tcav.o      dolist.o    prexln.o \
     tcave.o     tltrm.o     twbuf.o         prkwdv.o \
     prline.o    tchge.o     tluma.o     doprin.o \
     tmap.o      doread.o prnGlb.o tmast.o dorvrs.o prnflg.o \
     twinil.o    dostop.o        tconv.o     synradcl.o \
     twinit.o    dotemp.o    tconvm.o    drndsr.o \
     tcoord.o    drwkwd.o   tcoorde.o \
     mkplst.o    psub.o      tday.o      tmuld.o \
     wfbin.o     tdcmd.o     wfres.o eigs33.o       wiord.o     \
     ptimes.o    elname.o errmsg.o ptrace.o    tdinit.o    mrecal1.o \
     push.o      mrqcof.o    wtune.o \
     expln.o     mrqcov.o    tphplt.o     \
     mrqmin.o    datetime.o tdrife.o \
     filaux.o     tfgetm.o    tfojit.o tedrawf.o

OBJ5=trcoda.o nalign.o ndelw.o nfgetm.o msolvg.o yylex_.o calc.o \
     initdainterface.o initda.o tfepicsconstatcb.o

OBJAUTO=         abbrev.o    cputix.o     csinit.o     \
     doACT.o     eval1_.o    getbuf.o     getwrd.o    \
     ielm.o      itfgeto.o   itfgetbuf_.o \
     pgflag.o    pgmast.o     pgrmat.o \
     pgsolvcond.o qcell.o qdbend.o qdcell.o \
     qddrif.o    qdmdia.o    qdquad.o     qdthin.o   qdtwis.o \
     qgettr.o    qins.o      qmult.o       \
     qtent.o     qtwiss.o     tapert.o \
     tbdecoup.o  tbende.o     tbendi.o   tbfrie.o \
     tceigen.o    tcftr.o    tcod.o \
     tcorr.o     tdet.o      tdfun.o      tdgeo.o    tdlat.o \
     teigen.o    temap.o     temat.o      temitf.o \
     temits.o    terror.o    tfaprt.o     tfattr.o   tfcalc.o     tfchgv.o \
     tfchro.o    tfcoup.o    tfdapert.o   tfdisp.o   tfemit.o \
     tffamsetup.o tffile.o      tffscalc.o tffsfreefix.o \
     tffsmatch.o tfgeo.o     tfgetv.o     tfif.o     tfinit.o     \
     tfkwrd.o    tflag.o     tflogi.o     tfltr1.o   \
     tfmat.o     tfoptics.o  tfprint.o   \
     tfsave.o    tfsetparam.o tfsetv.o   tfshow.o    tftmat.o \
     tftrack.o   tftrak.o    tftrb.o \
     tftype.o    tfvars.o    tgauss_.o \
     tgetfv.o    tins.o     tinse.o \
     tlum.o      tmov.o      tmulte.o \
     tmulti.o    tmulta.o    tdrift.o     trackd.o \
     tmultr.o    tnorm.o     \
     tquade.o    track.o \
     tracka.o    trade.o     tsconv.o     tsgeo.o \
     tshow.o     tsol.o      tsole.o      tsolvg.o \
     tsolvm.o    tspac.o     tspini.o     tsteee.o   tsteer.o     tsvdm.o \
     tthine.o    ttinit.o    ttstat.o     tcsvdm.o \
     tturne.o    txcalc.o    wjfit.o      wtunem.o     \
     tturn.o     tbend.o     tquad.o      tflifetrack.o \
     tspch_.o     ft.o        psn.o       spch.o \
     csrtest.o csrtrack.o csroy.o txwake.o tbbbrem.o

OBJF=tfeval.o    itfcopy.o tfefun.o  tfmap.o     tfprinta.o  tfestk.o

OBJUTIL=utils.o

OBJOBS = mbmp.o temp.o pstati3.o palgn.o pvbump.o pwrtmon.o pwrtstr.o \
	mfdir.o mfdel.o mfdel1.o mrecal.o mstor1.o putsti.o ttdr.o \
	prkick.o prkick1.o pstati.o pstai1.o mwght.o corfree.o \
	phdrw.o phdrwa.o pkill.o ppair.o pqcell.o mdpmax.o monqu.o \
	wjdraw.o ptol.o  msolb.o pvert.o preadstr.o preadmon.o \
	tdlat.o mstack.o mstack1.o mstack2.o msort.o mstore.o pstati2.o \
	mrmb.o mhogal.o mstat.o mstatp.o pmovi.o mhogan.o monel.o \
	mfnst.o mstat2.o metaer.o trotg.o nlfit.o gaus3.o corinit.o \
	pvbump2.o pmeas.o pmbdata.o pvbump3.o ptrim.o pundo.o \
	pmbdrw.o pmbump.o pvbump1.o preabuf.o pmbdata1.o 

# from sad.builtin.mk
OBJFUNC=tfDefFuncs_.o $(_SAD_FUNC_OBJS)
OPT_COPT+=$(_SAD_FUNC_COPT)
OPT_INCDIR+=$(_SAD_FUNC_IOPT)
OPT_RLIBDIR+=$(_SAD_FUNC_ROPT)
OPT_LIBS+=$(_SAD_FUNC_LIBS)

# for SIM objects
ifndef OBJDYNL
OBJDYNL=sim/dynl-dummy.o
endif

OBJSIM_DYNL=sim/dynl.o $(OBJDYNL)

OBJSIM= sim/unix_pointer_.o sim/unix_memory_.o sim/unix_memory8_.o \
	sim/sad_api.o sim/sad_functbl.o sim/sad_signal.o \
	sim/sad_xlib.o sim/sad_tcltk.o \
	$(OBJ_FORTRAN) \
	$(OBJSIM_DYNL) \
	$(FRAMEWORK_OBJS) \
	$(SYS_OBJS)

# Define individual object part target after OBJ* variable definition
objmod: $(OBJMOD)

objrc: $(OBJRC)

obj0: $(OBJ0)

obj1: $(OBJ1)

obj2: $(OBJ2)

obj3: $(OBJ3)

obj4: $(OBJ4)

obj5: $(OBJ5)

obj6: $(OBJ6)

objauto: $(OBJAUTO)

objf: $(OBJF)

objutil: $(OBJUTIL)

objfunc: $(OBJFUNC)

objsim: $(OBJSIM)

objosdep: $(OBJOSDEP)

objglue: $(OBJGLUE)

# End of File

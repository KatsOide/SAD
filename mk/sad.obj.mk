# New Makefile -- Common Object Part
#
all_object:	objmod obj0 obj1 obj2 obj3 obj4 obj5 obj6 \
		objauto objf objrc objfunc \
		objutil objosdep objsim objglue MAIN.o

OBJ_LIBSAD=$(OBJMOD) $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) \
	$(OBJAUTO) $(OBJF) $(OBJRC) \
	$(OBJFUNC) \
	$(OBJUTIL) $(OBJOSDEP)

OBJMOD= tfstk.o	tfstrbuf.o	tfreadbuf.o	toplvl.o	tintrb.o \
	tftwiss.o	tftake.o	tfsolvemember.o \
	tfetok.o	tfbeamline.o	temit.o	tbrad.o	tfloor.o

OBJRC= 	tfefun1.o	tfsort.o	tfmodule.o	itfmaloc.o \
	tfwrite.o	itfaloc.o	tfmemcheck.o	tfdot.o		tfematrix.o \
	tfconvstr.o	tfpart.o	tfreplace.o	tftake.o	tfeval1.o \
	tfearray.o	tfeexpr.o	tfsetlist.o	tfeeval.o	tfmap.o	tftable.o \
	itfdepth.o	tfdset.o	tfsameq.o	itfpmat.o	tmatch.o \
	gamma.o		tfbessel.o	 spkick.o	tsolqu.o	tsolque.o \
	tfshared.o	tfinitn.o	tffsa.o

OBJ0= 	italoc.o qcav.o	JNLPRM.o	msolv1.o	prkick.o \
	prkick1.o	twsdrw.o	tdrwdt.o	msortn.o   tdplt.o     bbstrhl.o bbtool.o \
	ttdr.o	 undulator.o phsrot.o    matrix.o \
	autofg.o   tffswake.o  tfvectorize.o \
	itfmessage.o           tffindroot.o spkick.o \
	autos_.o      gdtoa/g_dfmt.o gdtoa/gdtoa.o gdtoa/g__fmt.o \
	gdtoa/gmisc.o gdtoa/dmisc.o  gdtoa/misc.o

OBJ1=ActGRA.o    filbuf.o    pwmatq.o    ActLie.o    \
     msort.o     pwrite.o    ActTra.o    flmgr.o \
     pwrlat.o    BEDGE.o     frand.o     mstack.o    tput.o \
     BMAG.o      mstack1.o   tputbs.o    CM24.o      gaus3.o \
     mstack2.o   DM42.o      gaussn.o    mstat.o \
     qchg.o      tqente.o    FMUL2.o     geodrw.o    mstat2.o \
     qcoord.o    termes.o    tqfrie.o    mstatp.o  \
     mstor1.o    tfadjst.o   tqlfre.o    INTGRL.o \
     gettok.o    mstore.o    IPAK.o      mwght.o     IgetGL.o    getvad.o \
     tfbeam.o    LgetGL.o    nlfit.o     qdrift.o    \
     Lrdnum.o    openP.o     hsrch.o     \
     MINV2.o     over.o      MMUL2.o     tfdset.o \
     trackb.o    MMUL4.o     pack.o      qmdiag.o

OBJ2=packpi.o    tfdbun.o    tracke.o    NewGRF.o \
     trad.o      QUADRU.o    padd.o      R0CAL.o \
     R1INV.o     inifil.o	R2INV.o     initb1.o    SEXTU.o     initbl.o \
     tfgeti.o    SYMPS.o     \
     qthin.o     tfgetr.o    trotg.o     SYMTR2.o    TCAL.o \
     TEDGE.o     TMATRS.o    rdexpr.o    tsdrad.o \
     TROT.o      rdkwdl.o    TWTRANS.o   lread.o \
     rdterm.o    UNITM2.o    lrflct.o    \
     tfltra.o    UNITM4.o    lstchk.o    rslvin.o \
     tserad.o    VMUL4.o     rvrsln.o \
     ZEROM2.o	ZEROM4.o	ZEROV4.o \
     tfrej.o     setdfl.o    tsmear.o    actPlt.o

OBJ3=sethtb.o    atof.o      bint.o      skipch.o    tspect.o    mbmp.o       pfcoup.o \
     sols33.o    pgaussj.o   corfree.o   spline.o \
     tstrad.o    corinit.o   phdrw.o     sprexl.o     \
     phdrwa.o    sprlin.o    pinner.o     \
     tfzap.o     ttcav.o     pkill.o \
     ttcave.o    talign.o    tgfun.o      \
     pmicad.o    tgrot.o     dAssgn.o    pmovi.o     defflg.o

OBJ4=defglb.o    ppair.o \
     tinitr.o    doelem.o    pqcell.o    doexpn.o \
     prSad.o     tvcorr.o    doflag.o     \
     prelem.o    title.o     doline.o    prelm0.o \
     tcav.o      dolist.o    prexln.o \
     tcave.o     tltrm.o     twbuf.o     mdpmax.o    prkwdv.o \
     metaer.o    prline.o    tchge.o     tluma.o     twelm.o     doprin.o \
     tmap.o      doread.o    mfdel.o     prnGlb.o  \
     tmast.o     dorvrs.o    mfdel1.o    prnflg.o   \
     twinil.o    dostop.o    mfdir.o     tconv.o     synradcl.o \
     twinit.o    dotemp.o    temp.o      mfnst.o     pstati.o    tconvm.o    drndsr.o \
     mhogal.o    pstati1.o   tcoord.o    drwkwd.o    mhogan.o    pstati2.o \
     tcoorde.o   pstati3.o \
     mkplst.o    psub.o      tday.o      tmuld.o \
     wfbin.o     tdcmd.o     wfres.o \
     eigs33.o    monel.o     wiord.o     monqu.o \
     ptimes.o    wjdraw.o    elname.o    ptol.o   \
     errmsg.o    mrecal.o    ptrace.o    tdinit.o    mrecal1.o \
     ptrim.o     pundo.o     tdjin.o  \
     mrmb.o      push.o      \
     mrqcof.o    wtune.o \
     expln.o     mrqcov.o    tphplt.o     \
     fdate1.o    mrqmin.o    putsti.o    datetime.o \
     msolb.o     tdrife.o \
     filaux.o    pvert.o     tprmpt.o    tfgetm.o    tfojit.o    tedrawf.o

OBJ5=trcoda.o    nalign.o    ndelw.o    pwrtstr.o \
     palgn.o	pvbump.o    nfgetm.o    preadstr.o   msolvg.o    preadmon.o \
     pwrtmon.o   yylex_.o    calc.o \
     pvbump2.o    pmeas.o    pmbdata.o   pvbump3.o \
     pmbdrw.o    pmbump.o    pvbump1.o    preabuf.o  pmbdata1.o  \
     initdainterface.o       initda.o     tfepicsconstatcb.o

OBJAUTO=         abbrev.o    cputix.o     csinit.o     \
     doACT.o     eval1_.o    getbuf.o     getwrd.o    \
     ielm.o      itfgeto.o   \
     pgflag.o    pgmast.o     pgrmat.o \
     pgsolvcond.o            qcell.o      qdbend.o     qdcell.o \
     qddrif.o    qdmdia.o    qdquad.o     qdthin.o   qdtwis.o \
     qgettr.o    qins.o      qmult.o      qquad.o \
     qtent.o     qtwiss.o     tapert.o \
     tbal.o      tbdecoup.o  tbende.o     tbendi.o   tbfrie.o \
     tceigen.o    tcftr.o    tcod.o \
     tcorr.o     tdet.o      tdfun.o      tdgeo.o    tdlat.o \
     teigen.o    temap.o     temat.o      temitf.o \
     temits.o    terror.o    tfaprt.o     tfattr.o   tfcalc.o     tfchgv.o \
     tfchro.o    tfcoup.o    tfdapert.o   tfdisp.o   tfemit.o \
     tffamsetup.o tffile.o	tffs.o      tffscalc.o tffsfreefix.o \
     tffsmatch.o tfgeo.o     tfgetv.o     tfif.o     tfinit.o     \
     tfkwrd.o    tflag.o     tflogi.o     tfltr1.o   \
     tfmat.o     tfoptics.o  tfprint.o   \
     tfsave.o    tfsetparam.o             tfsetv.o   tfshow.o    tftmat.o \
     tftrack.o   tftrak.o    tftrb.o \
     tftype.o    tfvars.o     tgauss_.o \
     tgetfv.o    thess.o     tins.o     tinse.o \
     tlum.o      tmov.o       tmulbs.o   tmulte.o \
     tmulti.o    tmulta.o    tdrift.o     trackd.o \
     tmultr.o    tnorm.o     tqr.o \
     tqrad.o     tquade.o    tquads.o     tquase.o   track.o \
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

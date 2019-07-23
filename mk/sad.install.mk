# New Makefile -- Insatll Part
#

ifneq ($(SADDIR),$(SAD_ROOT))

LIST_BENCHMARK=bench.sad bench2.sad \
	bench-IO.sad bench-PRNG.sad bench-Tkinter.sad

LIST_EXAMPLES=BuildInfo.sad design_example.sad examples.sad \
	amida.sad ColorTable.sad SymbolFonts.sad

LIST_DEMOS=demo-Execve.sad demo-Font.sad demo-KBFrame1.sad \
	demo-LabelFrame.sad demo-Pipe.sad demo-SigAction.sad \
	demo-Wait4.sad

LIST_PACKAGES=init.n init.local.n.example \
	CEPICSRecord.n CaSad.n CaSad2.n CaSad2c.n \
	Canvas.n Canvas3D.n CanvasCM.n CANVASDRAW.n \
	CheckEntryList.n ChooseAxis.n Class.n Couple.n \
	CursorEntry.n CSR.n D.n EPICSComponents1.n \
	element.n FileDialog.n FileSystem.n Fit.n \
	GeometryPlot.n HelpMessages.n \
	KBFFileSelectionPanel.n KBFrame.n KBMainFrame.n LDSolve.n \
	Library.n LifeTrack.n ListPlot.n Matrix.n MessageName.n Misc.n \
	NISTACK.n Optics.n OpticsPlot.n \
	OptionDialog.n PkgConfig.n Plots.n ProcessStatus.n rational.n \
	SADInspect.n SADInspect.sad SADTerminal.n Speak.n TabFrame.n \
	TkinterCore.n Tkinter.n TkPhoto.n TopDrawer.n \
	beamline.n bump.n corutility.n correction.n \
	dapert.n emit.n functions.n help.n optimize.n tableform.n \
	CreatePVname.n ChangeDotToUS.n Oracle.n PyFormat.n \
	Optionmenu.py kekbdb.py orsad.py rdbtool.py sad.py

LIST_KBFRAME=KBFrame.Medium.res KBFrame.res KBFrame.Medium.Xft.res KBFrame.Xft.res \
	a.xbm b.xbm c.xbm d.xbm e.xbm 16ton.gif blstipple.xbm \
	KEKBlogo2.gif KEKBlogo4.gif KEKBlogo4d.gif bb.xbm \
	cursors.xbm inc.xbm dec.xbm \
	dirup.xbm folder.gif KBF_Recent.gif \
	dnarr.gif rtarr.gif uparrow.xbm downarrow.xbm downarrowmask.xbm \
	print.gif print2.gif printer2.gif \
	connect.xbm disconnect.xbm toioc.xbm \
	constat1.xbm constat2.xbm constat3.xbm constat4.xbm constat5.xbm \
	medm.gif medm.xbm sad.gif sad.xbm tcltk.gif text.gif 

DOC_HTML=index.html index-ja.html

DOC_KBFrame_HTML_JA=index-ja.html \
	hajime-ja.html hello-ja.html subframe-ja.html about-ja.html \
	menu-ja.html progbar-ja.html statline-ja.html cmparr-ja.html \
	dialog-ja.html messagebox-ja.html filedialog-ja.html \
	emmeas-ja.html smplepics-ja.html cmptepics-ja.html \
	hardcopy-ja.html manual.html optionmenu-ja.html \
	camonmanual.html epicsdbmanual.html cursorentrymanual.html \
	kblogrd-ja.html

DOC_KBFrame_GIF=KEKBlogo4.gif KEKBlogo4d.gif \
	hajime1.gif hello1.gif hello2.gif sub1.gif sub2.gif \
	about1.gif about2.gif menu1.gif progbar.gif statline1.gif \
	cmparr1.gif dialog1.gif dialog2.gif messagebox1.gif messagebox2.gif \
	fileopen.gif filesave.gif \
	emmeas1.gif emmeas2.gif emmeas3.gif \
	smplepics1.gif smplepics2.gif cmptepics1.gif \
	pagesetup.gif print.gif \
	optionmenu1.gif optionmenu2.gif \
	symbolfont.gif

DOC_KBFrame_PDF=5d043.pdf

DOC_Design_Example_PNG=cell.png suppressor.png ring.png chromaticity.png

install:	install-exe install-pkg install-doc install-example

install-exe:	install-exec install-std-modules

install-exec:	$(SAD_ARCH_DIR)/bin/$(SADEXE) \
		$(SAD_SHARE_ROOT)/Packages/init.local.n \
		$(SAD_ROOT)/bin/gs $(SAD_ROOT)/bin/HostArch

$(SAD_ARCH_DIR)/bin/$(SADEXE):	$(SADEXE)
	@mkdir -p $(SAD_ARCH_DIR)/bin
	if [ -f $@ ]; then mv -f $@ $@.old; fi
	$(INSTALL_PROG) $< $@

$(SAD_SHARE_ROOT)/Packages/init.local.n:	$(SRCDIR)/init.local.n.in
	@mkdir -p $(SAD_SHARE_ROOT)/Packages
	if [ -f $@ -a -f $@.default ]; then \
	    if cmp -s $@ $@.default; then \
		rm $@; \
	    fi \
	fi
	sed -e "s,xSADROOT,$(SAD_ROOT),g" \
	    -e "s,xSAD_SHARE_ROOT,$(SAD_SHARE_ROOT),g" \
	    -e "s,xScriptDirectory,examples,g" \
	    $(SRCDIR)/init.local.n.in > $@.default && chmod 644 $@.default
	echo $@*
	if [ ! -f $@ ]; then \
	    cp $@.default $@; \
	fi

$(SAD_ROOT)/bin/gs:	$(SRCDIR)/gs.in
	@mkdir -p $(SAD_ROOT)/bin
	sed -e "s,xSADROOT,$(SAD_ROOT),g" \
	    -e "s,xSAD_SHARE_ROOT,$(SAD_SHARE_ROOT),g" \
	    -e "s,xSAD_ARCH_ROOT,$(SAD_ARCH_ROOT),g" \
	    -e "s,xSAD_EXE_ROOT,$(SAD_ARCH_ROOT),g" \
	    -e "s,xSADEXE,bin/$(SADEXE),g" \
	    -e "s,xSAD_ARCH_OR_OBJ,arch,g" \
	    -e "s,xSAD_PACKAGES,share/Packages,g" \
	    $< > $@ && chmod 755 $@

$(SAD_ROOT)/bin/HostArch:	$(BINDIR)/HostArch
	@mkdir -p $(SAD_ROOT)/bin
	$(INSTALL_PROG) $< $@

install-doc: install-sad-images \
		$(SAD_ROOT)/Documents/SADHelp.html \
		$(SAD_ROOT)/Documents/example/design_example.html
	@mkdir -p $(SAD_ROOT)/Documents/KBFrame
	for i in $(DOC_KBFrame_HTML_JA) $(DOC_KBFrame_GIF) \
		$(DOC_KBFrame_PDF); do \
	    $(INSTALL_DATA) $(SADDIR)/doc/KBFrame/$$i $(SAD_ROOT)/Documents/KBFrame/$$i; \
	done
	for i in $(DOC_HTML); do \
	    $(INSTALL_DATA) $(SADDIR)/doc/$$i $(SAD_ROOT)/Documents/$$i; \
	done

install-sad-images:
	@mkdir -p $(SAD_ROOT)/Documents/example
	for i in sad.png; do \
	    $(INSTALL_DATA) $(SADDIR)/doc/images/$$i $(SAD_ROOT)/Documents/$$i; \
	done
	for i in $(DOC_Design_Example_PNG); do \
	    $(INSTALL_DATA) $(SADDIR)/doc/images/design-$$i \
		$(SAD_ROOT)/Documents/example/$$i; \
	done

$(SAD_ROOT)/Documents/SADHelp.html: \
	$(SADDIR)/Documents/SADHelp.html
	@mkdir -p $(SAD_ROOT)/Documents
	$(INSTALL_DATA) $< $@

$(SAD_ROOT)/Documents/example/design_example.html: \
	$(SADDIR)/Documents/example/design_example.html
	@mkdir -p $(SAD_ROOT)/Documents/example
	$(INSTALL_DATA) $< $@

install-example:
	@mkdir -p $(SAD_ROOT)/examples/Benchmark $(SAD_ROOT)/examples/Demos
	@printf '%s' "Installing benchmark script files:"
	@for i in $(LIST_BENCHMARK); do \
	    printf ' %s' $${i}; \
	    $(INSTALL_DATA) $(SADDIR)/script/$$i $(SAD_ROOT)/examples/Benchmark/$$i; \
	done
	@printf '... done\n'
	@printf '%s' "Installing demo script files:"
	@for i in $(LIST_DEMOS); do \
	    printf ' %s' $${i}; \
	    $(INSTALL_DATA) $(SADDIR)/script/$$i $(SAD_ROOT)/examples/Demos/$$i; \
	done
	@printf '... done\n'
	@printf '%s' "Installing example script files:"
	@for i in $(LIST_EXAMPLES); do \
	    printf ' %s' $${i}; \
	    $(INSTALL_DATA) $(SADDIR)/script/$$i $(SAD_ROOT)/examples/$$i; \
	done
	@printf '... done\n'

install-pkg:
	@mkdir -p $(SAD_SHARE_ROOT)/Packages
	@printf '%s' "Installing Packages script files:"
	@for i in $(LIST_PACKAGES); do \
	    printf ' %s' $${i}; \
	    $(INSTALL_DATA) $(SADDIR)/Packages/$$i $(SAD_SHARE_ROOT)/Packages/$$i; \
	done
	@printf '... done\n'
	@printf '%s' "Installing KBFrame resouce files:"
	@mkdir -p $(SAD_SHARE_ROOT)/KBFrame
	@for i in $(LIST_KBFRAME); do \
	    printf ' %s' $${i}; \
	    $(INSTALL_DATA) $(SADDIR)/KBFrame/$$i $(SAD_SHARE_ROOT)/KBFrame/$$i; \
	done
	@printf '... done\n'

install-std-modules:
	@for mod in /dev/null $(STD_MODULES); do \
	    if [ -x $(SADDIR)/extensions/Standard/$${mod} ]; then \
		(cd $(SADDIR)/extensions/Standard/$${mod}; \
		 unset OS_NAME MAKEFLAGS; \
		 env SADSRCDIR=$(SADDIR) $(MAKE) install); \
	    fi; \
	done

else # $(SADDIR) = $(SAD_ROOT)

install:
	$(warning Now SADDIR = SAD_ROOT, insatll target has noeffect!)

install-exe:
	$(warning Now SADDIR = SAD_ROOT, insatll-exe target has noeffect!)

install-pkg:
	$(warning Now SADDIR = SAD_ROOT, insatll-pkg target has noeffect!)

install-doc:
	$(warning Now SADDIR = SAD_ROOT, insatll-doc target has noeffect!)

install-example:
	$(warning Now SADDIR = SAD_ROOT, insatll-example target has noeffect!)

endif

# End of File

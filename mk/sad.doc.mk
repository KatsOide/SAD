# New Makefile -- Document Part
#

doc:	$(SADDIR)/Documents/SADHelp.html

$(SADDIR)/Documents/SADHelp.html: \
	$(SADDIR)/script/design_example.sad.result \
	$(SADEXE) $(PKGDIR)/init.local.n $(BINDIR)/gs
	mkdir -p $(SADDIR)/Documents/example
	$(BINDIR)/gs $(SADDIR)/script/help2HTML.sad

$(SADDIR)/Documents/example/design_example.html:
	$(SADDIR)/script/design_example.sad.result \
	$(SADEXE) $(PKGDIR)/init.local.n $(BINDIR)/gs
	mkdir -p $(SADDIR)/Documents/example
	$(BINDIR)/gs $(SADDIR)/script/help2HTML.sad

$(SADDIR)/script/design_example.sad.result: \
	$(SADDIR)/script/design_example.sad \
	$(SADEXE) $(PKGDIR)/init.local.n $(BINDIR)/gs
	mkdir -p $(SADDIR)/Documents/example
	(cd $(SADDIR)/Documents/example ; \
	    $(BINDIR)/gs $< > $@ ; rm -f a fort.9)

# End of File

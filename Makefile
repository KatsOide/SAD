# makefile for SAD1
# rewritten by A. Morita Apr. 16, 2002
#
SHELL=/bin/sh

all:

ifndef YACC
YACC=yacc
endif

ifndef OS_NAME
VPATH=src
MAKE_ENV=\
	SADDIR=\`pwd\` \
	\`./bin/HostArch -config\`

all_object lib exe all:

clean distclean mostlyclean depend patch debug sharedclean \
	all_object lib exe doc all \
	tcltk-patch  tcltk-build  tcltk \
	libtai-patch libtai-build libtai \
	install-exe install-pkg install-doc install-example install:
	@$(SHELL) -c "$(MAKE) $@ $(MAKE_ENV)"

update_timezone:	./bin/gen_timezone_table
	@rm -f src/timezone_table.h
	@./bin/gen_timezone_table >src/timezone_table.h

update_calc_y:		src/calc.h src/calc.c

src/calc.h src/calc.c:	calc.y
	$(YACC) -d -o src/calc.c $<

.f.o: $(@:.o=.f)
	@$(SHELL) -c "$(MAKE) $@ $(MAKE_ENV)"

.c.o: $(@:.o=.f)
	@$(SHELL) -c "$(MAKE) $@ $(MAKE_ENV)"

else
SADMK=$(SADDIR)/mk
SRCDIR=$(SADDIR)/src
SADCONF=$(SADDIR)/config
OBJDIR=$(SADDIR)/obj/${MACH_ARCH}
VPATH=$(SRCDIR)
$(OBJDIR)/Makefile:	$(SADMK)/Makefile
	@mkdir -p $(OBJDIR)/inc $(OBJDIR)/sim $(OBJDIR)/gdtoa
	@cp $(SADMK)/Makefile $(OBJDIR)/Makefile

$(OBJDIR)/config.mk:
	@./bin/HostArch -config		>  $(OBJDIR)/config.mk
	@echo "SADDIR=$(SADDIR)"		>> $(OBJDIR)/config.mk
ifdef CONF_FILE
	@echo "CONF_FILE=$(CONF_FILE)"	>> $(OBJDIR)/config.mk
endif

clean distclean mostlyclean depend patch debug \
	all_object lib exe doc all \
	tcltk-patch tcltk-build tcltk \
	libtai-patch libtai-build libtai \
	install-exe install-pkg install-doc install-example install: \
	$(OBJDIR)/Makefile $(OBJDIR)/config.mk
	@(cd $(OBJDIR); $(MAKE) $@)

.f.o: $(@:.o=.f) $(OBJDIR)/Makefile $(OBJDIR)/config.mk
	@(cd $(OBJDIR); $(MAKE) $@)

.c.o: $(@:.o=.f) $(OBJDIR)/Makefile $(OBJDIR)/config.mk
	@(cd $(OBJDIR); $(MAKE) $@)
endif

# target for extensions
MODULE_DIRS=extensions devel
modules:
	@for dir in $(MODULE_DIRS); do \
	    if [ -d $${dir} -a -f $${dir}/Makefile ]; then \
		(cd $${dir}; $(MAKE) module); \
	    fi; \
	done

depend-modules:
	@for dir in $(MODULE_DIRS); do \
	    if [ -d $${dir} -a -f $${dir}/Makefile ]; then \
		(cd $${dir}; $(MAKE) depend); \
	    fi; \
	done

install-modules:
	@for dir in $(MODULE_DIRS); do \
	    if [ -d $${dir} -a -f $${dir}/Makefile ]; then \
		(cd $${dir}; $(MAKE) install); \
	    fi; \
	done

clean-modules:
	@for dir in $(MODULE_DIRS); do \
	    if [ -d $${dir} -a -f $${dir}/Makefile ]; then \
		(cd $${dir}; $(MAKE) clean); \
	    fi; \
	done

distclean-modules:
	@for dir in $(MODULE_DIRS); do \
	    if [ -d $${dir} -a -f $${dir}/Makefile ]; then \
		(cd $${dir}; $(MAKE) distclean); \
	    fi; \
	done

mostlyclean-modules:
	@for dir in $(MODULE_DIRS); do \
	    if [ -d $${dir} -a -f $${dir}/Makefile ]; then \
		(cd $${dir}; $(MAKE) mostlyclean); \
	    fi; \
	done

# End of File
#

dockerimage:Dockerfile
	docker build -t oldsadk64:build .

prune:
	docker container prune -f 
	docker image prune -f

runimage:
	xhost +
	docker run -itd --rm --net host -e DISPLAY=:0 --name=oldsadk64 oldsadk64:build bash	
	docker attach oldsadk64

runimagewithXauthority:
	docker run -it --net host --rm -e DISPLAY=:0 --name=oldsadk64 -v /run/user/1000/gdm/Xauthority:/root/.Xauthority -v /opt/SAD/examples:/opt/SAD/input oldsadk64:build 


saveimage:
	docker save oldsadk64:build -o oldsadk64-image.tar
	gzip oldsadk64-image.tar

loadimage:
	ungzip oldsadk64-image.tar.gz
	docker load -i oldsadk64-image.tar

exportContainer:
	docker export 

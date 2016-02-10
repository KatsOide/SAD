# New Makefile -- Post Rule Part
#
ifndef AFTERPORTMK
AFTERPORTMK=	YES

include $(SADMK)/sad.port.mk

AFTERPORTMK=	NO
endif
# End of File

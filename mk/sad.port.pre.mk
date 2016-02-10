# New Makefile -- Pre Rule Part
#
ifndef BEFOREPORTMK
BEFOREPORTMK=	YES

include $(SADMK)/sad.port.mk

BEFOREPORTMK=	NO
endif
# End of File

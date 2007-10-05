include make.defs

SRCDIR     := $(LIBSRC)
OBJDIR     := obj
LIBDIR     := lib
INCLUDEDIR := include
TESTDIR    := test

SRC := $(SRCDIR)/block_connectivity.cc
OBJ := $(patsubst %.cc,$(OBJDIR)/%.$(O),$(notdir $(SRC)))
LIB := $(LIBDIR)/$(LIBBLKC)

#--------------------------------------------------------------------
# Construction targets
#--------------------------------------------------------------------

all : $(LIB)

$(LIB) : $(OBJDIR) $(LIBDIR) $(INCLUDEDIR) $(OBJ)
	-@$(RM) $@
	$(AR) $@ $(OBJ) $(REPOSOBJ)
	$(INSTALL_DATA) $(SRCDIR)/block_connectivity.h $(INCLUDEDIR)/block_connectivity.h
	$(INSTALL_DATA) $(SRCDIR)/coord_transform.h $(INCLUDEDIR)/coord_transform.h

$(OBJDIR) :
	-$(MKDIR) $(OBJDIR)

$(LIBDIR) :
	-$(MKDIR) $(LIBDIR)

$(INCLUDEDIR) :
	-$(MKDIR) $(INCLUDEDIR)

depend : $(OBJDIR)
	-@$(RM) $(OBJDIR)/$(DEPBLKC)
	@for X in $(SRC); do \
           $(CXX) $(DEPFLAGS) $$X | $(SED) '/\.o/s/^/$(OBJDIR)\//' >> \
           $(OBJDIR)/$(DEPBLKC); done

install : $(LIB)
	$(INSTALL_DATA) $(SRCDIR)/block_connectivity.h $(INSTINCLUDEDIR)/block_connectivity.h
	$(INSTALL_DATA) $(SRCDIR)/coord_transform.h $(INSTINCLUDEDIR)/coord_transform.h
	$(INSTALL_DATA) $< $(INSTLIBDIR)/$(LIBBLKC)

#--------------------------------------------------------------------
# Testing
#--------------------------------------------------------------------

check : $(LIB)
	@echo ""
	@echo "----- testing the library -----"
	@cd test && $(MAKE)

#--------------------------------------------------------------------
# Output a distribution
#--------------------------------------------------------------------

dist : distclean
	cd ../ && tar -cvjf blkc.tar.bz2 block_connectivity

#--------------------------------------------------------------------
# Cleaning targets
#--------------------------------------------------------------------

uninstall :
	-$(RM) $(INSTINCLUDEDIR)/block_connectivity.h
	-$(RM) $(INSTINCLUDEDIR)/coord_transform.h
	-$(RM) $(INSTLIBDIR)/$(LIBBLKC)

clean :
	-$(RM) $(OBJ) $(REPOSOBJ)
	-$(RM) $(SRCDIR)/*~
	@cd test && $(MAKE) clean

distclean : clean
	-$(RM) $(LIB)
	-$(RM) $(INCLUDEDIR)/block_connectivity.h
	-$(RM) $(INCLUDEDIR)/coord_transform.h
	-$(RMDIR) $(OBJDIR) $(patsubst %/*.o,%,$(REPOSOBJ))
	-$(RM) make.defs
	-$(RM) config.log config.status

#--------------------------------------------------------------------
# Pattern rules
#--------------------------------------------------------------------

$(OBJDIR)/%.$(O) : $(SRCDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(DEFS) -o $@ -c $<

#--------------------------------------------------------------------
# Dependencies
#--------------------------------------------------------------------

-include $(OBJDIR)/$(DEPBLKC)

include ../make.defs

SRCDIR     := .
OBJDIR     := .
LIBDIR     := ../lib
INCLUDEDIR := ../include

SRC := $(SRCDIR)/testblkc.cc
OBJ := $(patsubst %.cc,$(OBJDIR)/%.$(O),$(notdir $(SRC)))

CXXFLAGS += -I$(INCLUDEDIR)
LDFLAGS += -L$(LIBDIR) -l$(subst lib,,$(basename $(LIBBLKC)))

#--------------------------------------------------------------------
# Construction targets
#--------------------------------------------------------------------

all : testblkc run

testblkc : $(LIBDIR)/$(LIBBLKC) $(OBJ)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)

run : testblkc
	@./testblkc

#--------------------------------------------------------------------
# Cleaning targets
#--------------------------------------------------------------------

clean :
	-$(RM) $(OBJ)
	-$(RM) $(SRCDIR)/*~
	-$(RM) testblkc

distclean : clean

#--------------------------------------------------------------------
# Pattern rules
#--------------------------------------------------------------------

$(OBJDIR)/%.$(O) : $(SRCDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(DEFS) -o $@ -c $<

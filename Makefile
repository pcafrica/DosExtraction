############################################################
## Verbosity mode check
############################################################
ifdef V
MUTE =
else
MUTE = @
endif

############################################################
## Compiler settings
############################################################
CXX = g++
CXXFLAGS = -g -std=c++11 -Wall --pedantic -fopenmp
LIBS = -lboost_iostreams -lboost_system -lboost_filesystem

LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS), 64)
   # do 64 bit stuff here, like set some CFLAGS
else
   # do 32 bit stuff here
endif

############################################################
## Source, header, object files and executables
############################################################
EXEC  = simulate_dos
VPATH = src/

SRCS = $(wildcard *.c++) $(wildcard $(VPATH)*.c++)
HDRS = $(wildcard *.h) $(wildcard $(VPATH)*.h)

OBJDIR = obj/
OBJS   = $(addprefix $(OBJDIR), $(notdir $(SRCS:.c++=.o)))

############################################################
## Documentation generation settings
############################################################
DOCDIR = doc

############################################################
## Targets definition
############################################################
.PHONY: default debug all all-debug profile astyle doc clean distclean

default: astyle all

debug: clean astyle all-debug

all: CXXFLAGS := -O3 $(CXXFLAGS)	# Target-specific variable.
all: $(EXEC)
	# DONE! ✓

all-debug: CXXFLAGS := -O0 $(CXXFLAGS)
all-debug: $(EXEC)
	# DONE! ✓

profile: $(EXEC)
	valgrind --tool=callgrind ./$^

$(EXEC): $(OBJS)
	# Linking and creating executable...
	$(MUTE)$(CXX) $(CXXFLAGS) $(LIBS) $(OBJS) -o $(EXEC)

astyle:
	# Formatting source codes...
ifeq ($(MUTE), )
	$(MUTE) astyle    -A10 -t4 -C -S -N -f -p -H -E $(addprefix $(shell pwd)/, $(SRCS) $(HDRS))
else
	$(MUTE) astyle -q -A10 -t4 -C -S -N -f -p -H -E $(addprefix $(shell pwd)/, $(SRCS) $(HDRS))
endif

obj/%.o: %.c++ %.h	# Check if header file has been modified.
	$(MUTE)mkdir -p $(OBJDIR)
	# Compiling $<...
	$(MUTE)$(CXX) $(CXXFLAGS) $(LIBS) $< -c -o $@

obj/$(EXEC).o: $(EXEC).c++	# Except for sources without a related header file.
	$(MUTE)mkdir -p $(OBJDIR)
	# Compiling $<...
	$(MUTE)$(CXX) $(CXXFLAGS) $(LIBS) $< -c -o $@

doc:
	$(MUTE)$(MAKE) -C $(DOCDIR) -s --no-print-directory > /dev/null

doc-clean:
	$(MUTE)$(MAKE) -C $(DOCDIR) clean -s --no-print-directory > /dev/null

clean:
	# Deleting object files...
	$(MUTE)$(RM) $(OBJS)
	# DONE! ✓

distclean: clean
	# Deleting executables...
	$(MUTE)$(RM) $(EXEC)
	# DONE! ✓
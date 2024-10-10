# Compiler
FC = gfortran
# Compiler flags
FFLAGS = -Wall -O2

#Libraries
LIBS = -llapack -lblas

# Directories
SRCDIR = .
SUBDIR = SUBROUTINES
MODDIR = MODULES

# Main program source file
MAIN = dagen.f95

# Subroutine files
SUBS = $(wildcard $(SUBDIR)/*.f95)

# Module files
MODS = $(wildcard $(MODDIR)/*.f95)

# Object files
#OBJS = $(patsubst %.f95,%.o,$(SUBS) $(MODS) $(MAIN))
OBJS = $(patsubst %.f95,%.o,$(MODS) $(SUBS) $(MAIN))

# Executable name
EXEC = DAGen

# Default target
all: $(EXEC)

# Rule to create the executable
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

# Rule to compile main program
$(SRCDIR)/%.o: $(SRCDIR)/%.f95
	$(FC) $(FFLAGS) -c $< -o $@

# Rule to compile subroutines
$(SUBDIR)/%.o: $(SUBDIR)/%.f95
	$(FC) $(FFLAGS) -c $< -o $@

# Rule to compile modules
$(MODDIR)/%.o: $(MODDIR)/%.f95
	$(FC) $(FFLAGS) -c $< -o $@

# Clean up object and module files
clean:
	rm -f $(OBJS) $(EXEC)

# Clean all object files and executable
veryclean: clean
	rm -f *.mod

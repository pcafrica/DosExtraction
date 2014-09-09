############################################################
## Verbosity mode check
############################################################
ifdef V
MUTE =
else
MUTE = @
endif

############################################################
## Documentation generation settings
############################################################
DOCGEN   = doxygen
DOXYFILE = config.dox

############################################################
## Targets definition
############################################################
.PHONY: default doc clean

default: doc

doc:
	# Creating documentation files...
	$(MUTE)$(DOCGEN) $(DOXYFILE) > /dev/null
	#     Compiling LaTeX...
	$(MUTE)$(MAKE) -C ./latex -s --no-print-directory > /dev/null
	# DONE! ✓

clean:
	# Deleting documentation files...
	$(MUTE)$(RM) -r ./*/
	# DONE! ✓
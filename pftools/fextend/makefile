# Makefile for fextend module

# ###################################################
# Definitions
# ###################################################
SHELL := bash
HOSTNAME := mec050570

include $(HOSTNAME).env

.PHONY: all clean

all: fnc_reader snc_reader

fnc_reader:
	$(F2PYC) -c --fcompiler=$(FCV) --opt=-O3 -I$(NCFI) pf_params.f90 fnc_reader.f90 -m fnc_reader -L$(NCFL) -lnetcdff -L$(NCL) -lnetcdf

snc_reader:
	$(F2PYC) -c --fcompiler=$(FCV) --opt=-O3 -I$(NCFI) pf_params.f90 snc_reader.f90 -m snc_reader -L$(NCFL) -lnetcdff -L$(NCL) -lnetcdf

clean:
	rm fnc_reader$(EXT_SUFFIX)
	rm snc_reader$(EXT_SUFFIX)

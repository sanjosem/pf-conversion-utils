# Makefile for fextend module

# ###################################################
# Definitions
# ###################################################
SHELL := bash
HOSTNAME := kootenay

include $(HOSTNAME).env

all: fnc_reader snc_reader
	
fnc_reader:
	$(F2PYC) -c --fcompiler=$(FCV) -I$(NCFI) pf_params.f90 fnc_reader.f90 -m fnc_reader -L$(NCFL) -lnetcdff -L$(NCL) -lnetcdf
	
snc_reader:
	$(F2PYC) -c --fcompiler=$(FCV) -I$(NCFI) pf_params.f90 snc_reader.f90 -m snc_reader -L$(NCFL) -lnetcdff -L$(NCL) -lnetcdf

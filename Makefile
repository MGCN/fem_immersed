
###############################################################################
############ IMMERSED_BOUNDARY Application Standard Makefile ##################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project 
# MODULE_DIR       - Location of the MOOSE modules directory
# FRAMEWORK_DIR    - Location of the MOOSE framework
#

###############################################################################
IB_DIR             ?= $(shell pwd)
#MOOSE_DIR          ?= $(MOOSE_DIR)/moose
MOOSE_DIR          ?= $(shell dirname $(MECH_DIR))/moose
#MOOSE_DIR        ?= $(shell dirname $(MooseImpact_DIR))/moose (this for the cluster)
#MOOSE_DIR         ?= /opt/krause_builds/moose (this for the mac environment)
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
SRC_DIR            ?= $(IB_DIR)/src


###############################################################################
# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk
# dep apps
APPLICATION_DIR    := $(IB_DIR)
APPLICATION_NAME   := immersed_boundary
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here

include $(UTOPIA_FE_DIR)/config/utopia_fe_config.makefile
#include $(UTOPIA_DIR)/config/utopia-config.makefile

ADDITIONAL_INCLUDES+=$(UTOPIA_FE_INCLUDES)
ADDITIONAL_LIBS+=$(UTOPIA_FE_LIBRARIES)

clean:
	@./auxiliaryCleaner.sh;



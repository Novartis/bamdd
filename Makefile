# makefile written using https://yuukidach.github.io/p/makefile-for-projects-with-subdirectories/ as template 

TARGET = website

OUTDIR = ./docs
SRCDIR = ./src
##DIR_OBJ = ./obj

OUTDIR_ABS=$(abspath $(OUTDIR))
PROJROOT_ABS=$(abspath .)

# needed to ensure that sub make calls as used by cmdstanr do not emit "Entering/Leaving" directory messages
# these additional messages interfer otherwise with cmdstanr which
# infers via make print-STANCFLAGS - for example - the value of these variables
# Todo: consider moving this bit into the configuration of R
MAKEFLAGS += --no-print-directory

# scripts and includes
INCS = $(SRCDIR)/setup.R
INCS += $(SRCDIR)/cmq_brm.R
INCS += $(SRCDIR)/_macros.qmd
INCS += $(SRCDIR)/install_dependencies.R

SRCS = $(wildcard [0-9][0-9]*.qmd $(foreach fd, $(SRCDIR), $(fd)/[0-9][0-9]*.qmd))
NODIR_SRC = $(notdir $(SRCS))
OBJS = $(SRCS:qmd=html)

# OBJS = $(patsubst %.qmd, %.md, $(SRCS))

# set default R and quarto cmd's
R_CMD ?= `R RHOME`/bin/R
QUARTO_CMD ?= quarto
# 

#INC_DIRS = -I./ $(addprefix -I, $(SUBDIR))
#LIBS = -largp
#LIB_DIRS = -L/usr/local/Cellar/argp-standalone/1.3/lib

# rule which allows to monitor the R code as it gets executed; useful for debugging
%.Rlog : %.qmd
	cd $(@D); $(R_CMD_PREFIX) $(R_CMD) -q -e "knitr::purl('$(<F)', '$(patsubst %.qmd,%.R,$(<F))')" -e "source('$(patsubst %.qmd,%.R,$(<F))', echo=TRUE)" | tee $(@F) 2>&1

# tell makefile how to turn a qmd files into things we want

%.html : %.qmd
	@echo running per-document html quarto render $(<F) --to html --profile $(QUARTO_PROFILE)
	cd $(@D); $(QUARTO_CMD_PREFIX) $(QUARTO_CMD) render $(<F) --to html --profile $(QUARTO_PROFILE) $(QUARTO_CMD_POSTFIX)
PHONY := $(TARGET)

# enforce that we first setup the web-site and then run the remainder
# potentially with parallelism. TODO: enforcing the order should
# actually be possible to do via dependencies.
website: ## Render the full Quarto website
	$(MAKE) $(OUTDIR)/$(QUARTO_PROFILE)/bookdown-website

PHONY += website-public
website-public: build/docs/public/bookdown-website ## Safely render the public-profile website (fully content-hidden filtered) into build/docs/public
	@echo "Public website rendered into build/docs/public - safe for external export."

PHONY += public-export
public-export: build/brms-examples-public.tar.gz ## Safely build the filtered public sources tarball (for updating the public GitHub repo)
	@echo "Public sources tarball ready: build/brms-examples-public.tar.gz - safe to publish to the public GitHub repo."

PHONY += help
help: ## Show this help message
	@grep -E '^[a-zA-Z_%.-]+:.*##' $(MAKEFILE_LIST) \
	  | awk 'BEGIN {FS = ":.*## "}; {printf "  \033[36m%-22s\033[0m %s\n", $$1, $$2}'

# $(OUTDIR)/$(QUARTO_PROFILE)/%.html : %.qmd
# 	@echo running per-document html quarto render $< --profile $(QUARTO_PROFILE)
# 	$(QUARTO_CMD) render $< --profile $(QUARTO_PROFILE)

# note that the md's are run in the directory itself to cast the
# Quarto rendering process to be a standalone run. The intended effect
# is to populate the _cache directory of the respective document which
# can then get reused when rendering the full web-site
%.md : %.qmd
	@echo running per-document md quarto render $(<F) --to md --output $(@F) --profile $(QUARTO_PROFILE)
	cd $(@D); $(QUARTO_CMD) render $(<F) --to md --output $(@F) --profile $(QUARTO_PROFILE)

%.pdf : %.qmd
	@echo running per-document pdf quarto render $< --to pdf
	$(QUARTO_CMD_PREFIX) $(QUARTO_CMD) render $< --to pdf --output $(@F) --profile $(QUARTO_PROFILE) $(QUARTO_CMD_POSTFIX)

$(OUTDIR)/$(QUARTO_PROFILE)/bookdown-website : index.qmd _quarto.yml $(OBJS)
	@echo Rendering website
	$(QUARTO_CMD) render --to html --profile $(QUARTO_PROFILE)
	touch $(OUTDIR)/$(QUARTO_PROFILE)/bookdown-website

# ensure that R markdowns honor dependency on setup file
$(OBJS) : $(SRCDIR)/setup.R

#  \{\.content-hidden when-profile="public"\}[\s\S]*?::://g' $< > $@
# :::

#$(DIR_OBJ)/%.o: %.c $(INCS)
#    mkdir -p $(@D)
#    $(CC) -o $@ $(CFLAGS) -c $< $(INC_DIRS)

PHONY += clean
clean: ## Remove all build outputs and caches
	rm -rf $(OUTDIR)/*
	rm -rf src/*.html
	rm -rf brms-cache
	rm -rf .brms-cache
	rm -rf _brms-cache
	rm -rf src/*_files
	rm -rf src/*_cache
	rm -rf src/_freeze
	rm -rf _freeze
	rm -rf src/.quarto
	rm -rf .quarto
	rm -rf build

PHONY += echoes
echoes: ## Print INC, SRC, and OBJ file lists
	@echo "INC files: $(INCS)"
	@echo "SRC files: $(SRCS)"
	@echo "OBJ files: $(OBJS)"

## Debug target that allows you to print a variable (usage: make print-VARNAME)
print-%: ## Print the value of a make variable (e.g. make print-SRCS)
	@echo $* = $($*)

.PHONY = $(PHONY)

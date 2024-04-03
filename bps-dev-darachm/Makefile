# This is a self-documenting Makefile, so everything after a double octothorpe
# will print out when you do that.
#
# - meaning, if you change stuff, document it with the double octothorpe/hash
#   signs plz

.PHONY: all docs help install

## Help for this makefile in `bps-dev`:
## 
help: Makefile
	@sed -n 's/^##//p' $<

## install		This just tells you what to install! Not installing it for you.
install:
	@echo "Go install `nextflow` and `singularity` to run the pipelines, and `python3`, and `tox` for making the docs."

## docs			This builds the sphinx docs, again using tox
docs:
	python3 -m tox -e build_docs


## call			Do the BPS thing
call: main.nf nextflow.config config/config.yaml
	nextflow $(word 1,$^) -c $(word 2,$^) \
		-resume -ansi-log false -with-dag reports/dag.html \
		-profile common,blocky \
		--slurm_queue low30 \
		-params-file $(word 3,$^)

## conf			Just dump config
conf: main.nf nextflow.config config/config.yaml
	nextflow $(word 1,$^) -c $(word 2,$^) \
		-resume -ansi-log false -preview -with-dag reports/dag.html \
		-profile common,blocky \
		--just_config \
		--slurm_queue low30 \
		-params-file $(word 3,$^)

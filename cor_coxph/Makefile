## all           : clean saved results and get results to populate report
all: clean cor_coxph report

## clean         : delete saved bootstrap results
clean:
	rm -rf output/$(TRIAL)/*

## cor_coxph     : create tables and figures
cor_coxph:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_boost))
	Rscript code/cor_coxph_moderna_boost.R

endif

report:
	bash _build_chapter.sh


collate:
	Rscript code/cor_coxph_collate_meta_VE.R
	Rscript -e "rmarkdown::render('collate_ensemble.Rmd', output_file='covpn_correlates_cor_coxph_janssen_collated.pdf')"

meta:
	Rscript code/cor_coxph_meta_VE_curves.R
	Rscript -e "rmarkdown::render('meta_cve.Rmd', output_file='covpn_correlates_cor_coxph_1meta.pdf')"


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_coxph

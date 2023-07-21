## all            : clear current results and make figures
all: clean figs

## figs           : clean data and make figures
figs: 
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_boost))
	Rscript code/cor_data_preprocess.R # Yiwen's
	Rscript code/cor_graphics_all.R

endif

## clean          : clean up data_clean and output folder
clean:
	rm -f data_clean/* && rm -rf output/$(TRIAL)/*

## report        : create a standalone report
report:

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: figs clean all

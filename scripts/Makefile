
# run this with make -j 36 all for example so that it uses a decent number of processors.

pcfg_files = $(wildcard ../data/test*/*.pcfg)
wcfg_files = $(patsubst %.pcfg,%.wcfg,$(wildcard ../data/test*/*.pcfg))
json_files = $(patsubst %.pcfg,%.json,$(wildcard ../data/test*/*.pcfg))


../diagrams/boxplot_density.pdf : $(json_files)
	mkdir -p ../diagrams
	python ../testpcfg/plot_diagram_density.py

../diagrams/boxplot_log_density.pdf : $(json_files)
	mkdir -p ../diagrams
	python ../testpcfg/plot_diagram_density.py

../diagrams/hbar_lpf.pdf : $(json_files)
	mkdir -p ../diagrams
	python ../testpcfg/plot_hbar_lpf.py

../diagrams/hbar_extra_rules.pdf : $(json_files)
	mkdir -p ../diagrams
	python ../testpcfg/plot_hbar_extra_rules.py

../diagrams/hbar_lpf_iso.pdf : $(json_files)
	mkdir -p ../diagrams
	python ../testpcfg/plot_hbar_lpf_iso.py

../diagrams/density_pf.pdf : $(json_files)
	mkdir -p ../diagrams
	python ../testpcfg/plot_density_pf.py

../diagrams/anchors.pdf : $(json_files)
	mkdir -p ../diagrams
	python ../testpcfg/plot_diagram_anchor.py

%.wcfg : %.pcfg
	python ../testpcfg/oracle_test.py  $<  $@

%.json :  %.pcfg %.wcfg
	python ../testpcfg/evaluate.py   --json $@  $^

diagrams:  ../diagrams/anchors.pdf ../diagrams/density_pf.pdf ../diagrams/hbar_lpf_iso.pdf ../diagrams/boxplot_density.pdf ../diagrams/hbar_lpf.pdf ../diagrams/hbar_extra_rules.pdf

json : $(json_files)

clean : 
	rm $(wcfg_files)

.PHONY: all clean
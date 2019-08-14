

pcfg_files = $(wildcard ../data/test*/*.pcfg)
wcfg_files = $(patsubst %.pcfg,%.wcfg,$(wildcard ../data/test*/*.pcfg))

%.wcfg : %.pcfg
	python ../testpcfg/oracle_test.py  $<  $@

all : $(wcfg_files)

clean : 
	rm $(wcfg_files)

.PHONY: all clean
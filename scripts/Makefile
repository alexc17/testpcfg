
# run this with make -j 36 all for example so that it uses a decent number of processors.

pcfg_files = $(wildcard ../data/test*/*.pcfg)
wcfg_files = $(patsubst %.pcfg,%.wcfg,$(wildcard ../data/test*/*.pcfg))
json_files = $(patsubst %.pcfg,%.json,$(wildcard ../data/test*/*.pcfg))

%.wcfg : %.pcfg
	python ../testpcfg/oracle_test.py  $<  $@

%.json :  %.pcfg %.wcfg
	python ../testpcfg/evaluate.py   --json $@  $^

all : $(json_files)

clean : 
	rm $(wcfg_files)

.PHONY: all clean
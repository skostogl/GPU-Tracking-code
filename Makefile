$(VERBOSE).SILENT:

all: tracker NAFF

tracker:
	$(MAKE) -C $@

NAFF:
	$(MAKE) -C $@

install: all
	echo "Putting modules in place..."
	cp tracker/tracker.so modules/
	cp NAFF/NAFF.so modules/

clean:
	$(MAKE) -C tracker clean
	$(MAKE) -C NAFF clean
	rm -f $(wildcard modules/*.so)

.PHONY: tracker NAFF clean


include Makefile.common

$(VERBOSE).SILENT:

all: tracker NAFF

tracker:
	$(MAKE) -C $@

NAFF:
	$(MAKE) -C $@

install: all
	echo "Installing modules in $(MODULES_FOLDER)..."
	cp tracker/tracker.so $(MODULES_FOLDER)
	cp NAFF/NAFF.so $(MODULES_FOLDER)

clean:
	$(MAKE) -C tracker clean
	$(MAKE) -C NAFF clean
	rm -f $(wildcard $(MODULES_FOLDER)/*.so)

.PHONY: tracker NAFF clean


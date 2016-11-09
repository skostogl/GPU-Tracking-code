#SUBDIRS := tracker
SUBDIRS := tracker NAFF

all: $(SUBDIRS)
	cp tracker/tracker.so modules/
	cp NAFF/NAFF.so modules/

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: all $(SUBDIRS)


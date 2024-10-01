##
# MDS Scope
#

SHELL = bash
ALPHA ?= 2
K ?= 4

CPPFLAGS += -Wall
CXXFLAGS += -O3 -DNDEBUG -std=gnu++20 -pthread -DALPHA=$(ALPHA) -DK=$(K)
LDFLAGS += -pthread

BUILDDIR = A${ALPHA}K$(K)

# If libxxhash available on system, use that, otherwise download and install
LIBXXHASH ?= $(shell pkg-config -cflags libxxhash >& /dev/null && echo "sys" || echo "build")
LIBXXHASH_MK = $(BUILDDIR)/libxxhash_$(LIBXXHASH).mk
LIBXXHASH_URL = https://github.com/Cyan4973/xxHash/archive/refs/tags/v0.8.2.tar.gz

PROGRAMS = traverse_comp mdss2dot comp2rankdot fms2mds optimize_rem_path_len	\
mykkeltveit_set champarnaud_set sketch_components syncmer_set frac_set			\
create_seed sketch_histo old_champarnaud_set opt_canon

EXECS = $(addprefix $(BUILDDIR)/, $(PROGRAMS))

all: $(EXECS)

$(BUILDDIR):
	mkdir -p $@

$(BUILDDIR)/%.o: %.cc $(BUILDDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

$(BUILDDIR)/common.ar: $(BUILDDIR)/backtrace.o $(BUILDDIR)/common.o $(BUILDDIR)/sequence.o
	$(AR) $(ARFLAGS) $@ $^

$(BUILDDIR)/%: $(BUILDDIR)/%.o $(BUILDDIR)/common.ar
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm $(BUILDDIR)/*.o $(BUILDDIR)/*.ar

cleanall:
	rm -rf $(BUILDDIR)


# Handle Libxxhash sub-makefile. If on system, use pkg config. Otherwise, download, build and set
$(BUILDDIR)/libxxhash_sys.mk: $(BUILDDIR)
	{ echo -n "CXXFLAGS += "; pkg-config -cflags libxxhash; } > $@
	{ echo -n "LDFLAGS += "; pkg-config --libs-only-L libxxhash; } >> $@
	{ echo -n "LDFLAGS += "; pkg-config --libs-only-L libxxhash | sed 's/-L/-Wl,-rpath,/g'; } >> $@
	{ echo -n "LDLIBS += "; pkg-config --libs-only-l libxxhash;  } >> $@

$(BUILDDIR)/libxxhash_build.mk: $(BUILDDIR)
	mkdir -p $</xxHash
	curl -L $(LIBXXHASH_URL) | tar -zx --strip-components=1 -C $(BUILDDIR)/xxHash && cd $(BUILDDIR)/xxHash && $(MAKE) && $(MAKE) install PREFIX=$$(pwd)/inst
	./pkg2make.sh $$(find $(BUILDDIR)/xxHash/inst -iname '*.pc') > $@

include $(LIBXXHASH_MK)

##
# MDS Scope
#

ALPHA ?= 2
K ?= 4

CPPFLAGS += -Wall
CXXFLAGS += -O3 -DNDEBUG -std=gnu++20 -DALPHA=$(ALPHA) -DK=$(K)

BUILDDIR = A${ALPHA}K$(K)

# If libxxhash available on system, use that, otherwise download and install
LIBXXHASH ?= $(shell pkg-config -exists libxxhash && echo -n "sys" || echo -n "build")
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
	cd $< && curl -L $(LIBXXHASH_URL) | tar -zx --strip-components=1 -C xxHash && cd xxHash && make && make install PREFIX=$$(pwd)
	{ echo -n "CXXFLAGS += "; PKG_CONFIG_PATH=$(BUILDDIR)/xxHash/lib/pkgconfig pkg-config -cflags libxxhash; } > $@
	{ echo -n "LDFLAGS += "; PKG_CONFIG_PATH=$(BUILDDIR)/xxHash/lib/pkgconfig pkg-config --libs-only-L libxxhash; } >> $@
	{ echo -n "LDFLAGS += "; PKG_CONFIG_PATH=$(BUILDDIR)/xxHash/lib/pkgconfig pkg-config --libs-only-L libxxhash | sed 's/-L/-Wl,-rpath,/g'; } >> $@
	{ echo -n "LDLIBS += "; PKG_CONFIG_PATH=$(BUILDDIR)/xxHash/lib/pkgconfig pkg-config --libs-only-l libxxhash;  } >> $@

include $(LIBXXHASH_MK)

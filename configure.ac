AC_INIT([traverse_components_gen], [0.0.1], [gmarcais@cs.cmu.edu])

AS_IF([test -f configure],
      [AC_MSG_FAILURE([Must be run in a different directory that './configure.' The name of the directory is used as the configuration name.])])
AC_ARG_VAR([CONFNAME])
CONFNAME=$(basename $(pwd))
echo $CONFNAME


AC_ARG_VAR([ALPHA], [Alphabet size])
AC_ARG_VAR([K], [K-mer length])

AS_IF([test "x$ALPHA" = "x" -o "x$K" = "x"],
      [AC_MSG_FAILURE([ALPHA and K variables are required])])

# Check for yaggo
AC_ARG_VAR([YAGGO], [Yaggo switch parser generator])
AS_IF([test "x$YAGGO" = "x"], [AC_PATH_PROG([YAGGO], [yaggo])])
AS_IF([test "x$YAGGO" = "x"], [AC_MSG_FAILURE([yaggo parser not found])])

AC_SUBST([OPTFLAGS])
AC_ARG_ENABLE([debug], [Enable debugging symbols and disable optimizations], [OPTFLAGS="-O0 -g"], [OPTFLAGS="-O3 -DNDEBUG"])

PKG_CHECK_MODULES([XXHASH], [libxxhash])

AC_LIB_LINKFLAGS_FROM_LIBS([XXHASH_LDADD], [XXHASH_LIBS])

AC_CONFIG_FILES([${CONFNAME}.config:tup.config.in])
AC_OUTPUT()

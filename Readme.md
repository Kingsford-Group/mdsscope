# What is it?

Traverse the component graph, maybe optimizing for a particular property.
Generalized to arbitrary sized alphabets.

# Compilation
For the first time, create configure script: `autoreconf -fi`

For each desired configuration, run `configure` and `tup variant`.
For example:

```Shell
(NAME=A2K4; mkdir configs/$NAME; (cd configs/$NAME; ../../configure ALPHA=2 K=4); tup variant configs/$NAME/$NAME.config)
```

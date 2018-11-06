# iid
A set of mathematical and statistical libraries with no external dependencies.

## Motivation

One of the awesome things about the R programming environment is the access to high quality
numerical implementations of statistical functions, that provide excellent performance, especially
for the tails of the distributions, required in many advanced bioinformatics (and other) applications.

While there are some excellent libraries in other languages (the Boost C++ libraries are excellent),
creating appropriate glue into some of the interpreted languages is messy and/or limits portability.
For a fully portable set of mathematical and statistical functions, we want minimal dependecies on
external packages which are not always available. The set of implementations in this library depend
on a minimal set of existing library functions (sqrt, log, and exp), and use only core language
features in order to maximise portability.

The primary goal of this library is coverage and precision, not runtime performance, though we have
avoided unnecessarily slow algorithms where possible.

## Supported languages

### Python

Pure python allows the same code to run under any python implementation, including Pypy, which at
the time of writing does not have a full numpy/scipy stack.


### Lua (work in progress)

Pure Lua enables access even in embedded environments, and sandboxed environments where  linking to
external compiled libraries is not possible.

### Javascript (planned)

A Javascript implementation allows both server-side (e.g. nodejs) and client-side access to high quality
numerical implementations.

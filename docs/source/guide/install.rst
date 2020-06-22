.. _install:

Installation Guide
==================

Python
------
The easiest way to install RCR into Python is with ``pip``:

.. code-block:: python
    
    python3 -m pip3 install rcr

The only requirement is ``pybind11``, but that should be taken care of with the pip installer. If not,
run ``python3 -m pip3 install pybind11``.

Alternatively, one could also download the C++ source code (see below), 
clone the `pybind11 Github repository <https://github.com/pybind/pybind11>`_ into the same directory,
and compile everything into a Python module manually, but installing via pip as above is *much* easier.

C++
---
Because the RCR Python library uses pybind11 to wrap the original C++ source code seamlessly into Python, 
all of the speed of C++ is available through the Python library. However, if the C++ source code is desired, it
can be found at the `github repository <https://github.com/nickk124/RCR>`_ in 
the ``src`` directory (``RCR_python.cpp`` is only used for wrapping the C++ code into the rcr Python library, so it can be ignored).
Documentation specific to the C++ codebase can be found within the directory ``docs/cpp_docs`` of this repository.
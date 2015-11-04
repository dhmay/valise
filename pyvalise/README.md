Valise
======

My personal toolkit.

## Environment setup (Mac)

### non-Python

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew install swig
    xcode-select --install

Ensure the valise directory is on PYTHONPATH

#### for mayavi 3D charts

install VTK

install mayavi

### Python

Note: before installing lxml, run this magic line:

    export C_INCLUDE_PATH=/Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/usr/include/libxml2:$C_INCLUDE_PATH

OK, then do this:

    easy_install Cython
    pip install distribute numpy ipython nose sympy pandas matplotlib 
    pip install scipy 
    pip install matplotlib_venn lxml pyteomics biopython pandas patsy 
    pip install statsmodels pysam xlsxwriter sphinx 
    pip install sphinxarg sphinxcontrib-programoutput
    pip install rpy2

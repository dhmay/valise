Valise
======

My personal toolkit.

## Environment setup

### non-Python

    ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)"
    brew install swig
    xcode-select --install

Ensure the valise directory is on PYTHONPATH

### Python

    easy_install Cython
    pip install distribute numpy ipython nose sympy pandas matplotlib 
    pip install scipy 
    pip install matplotlib_venn lxml pyteomics biopython pandas patsy 
    pip install statsmodels pysam xlsxwriter sphinx 
    pip install sphinxarg sphinxcontrib-programoutput
    pip install rpy2

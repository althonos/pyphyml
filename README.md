# φ🌲 PyPhyML

*[Cython](https://cython.org/) bindings and Python interface to [PhyML 3.0](https://github.com/stephaneguindon/phyml), a method for estimating maximum-likelihood phylogenies.*

<!-- [![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyphyml/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyphyml/actions) -->
<!-- [![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyphyml?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pyphyml/) -->
[![License](https://img.shields.io/badge/license-GPL--3%2e0--or--later-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
<!-- [![PyPI](https://img.shields.io/pypi/v/pyphyml.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyphyml) -->
<!-- [![Bioconda](https://img.shields.io/conda/vn/bioconda/pyphyml?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pyphyml) -->
<!-- [![AUR](https://img.shields.io/aur/version/python-pyphyml?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyphyml) -->
<!-- [![Wheel](https://img.shields.io/pypi/wheel/pyphyml.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyphyml/#files) -->
<!-- [![Python Versions](https://img.shields.io/pypi/pyversions/pyphyml.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyphyml/#files) -->
<!-- [![Python Implementations](https://img.shields.io/pypi/implementation/pyphyml.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyphyml/#files) -->
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyphyml/)
[![Mirror](https://img.shields.io/badge/mirror-LUMC-003EAA.svg?maxAge=2678400&style=flat-square)](https://git.lumc.nl/mflarralde/pyphyml/)
[![Issues](https://img.shields.io/github/issues/althonos/pyphyml.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyphyml/issues)
<!-- [![Docs](https://img.shields.io/readthedocs/pyphyml/latest?style=flat-square&maxAge=600)](https://pyphyml.readthedocs.io) -->
<!-- [![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyphyml/blob/main/CHANGELOG.md) -->
<!-- [![Downloads](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fpepy.tech%2Fprojects%2Fpyphyml&search=%5B0-9%5D%2B.%5B0-9%5D%2B(k%7CM)&style=flat-square&label=downloads&color=303f9f&cacheSeconds=86400)](https://pepy.tech/project/pyphyml) -->

## 🗺️ Overview

PhyML is a software package developed by [Stephane Guindon](https://github.com/stephaneguindon)
and [Olivier Gascuel](https://scholar.google.com/citations?user=vlLa_zgAAAAJ&hl=fr)[\[1\]](#ref1) 
that uses modern statistical approaches to analyse alignments of nucleotide or 
amino acid sequences in a phylogenetic framework. The main PhyML tool builds 
phylogenies under the maximum likelihood criterion, and implements a large number 
of substitution models coupled to efficient options to search the space of 
phylogenetic tree topologies.

PyPhyML is a Python module that proivides bindings to PhyML using [Cython](https://cython.org/).
It allows creating alignments, inferring a phylogenetic tree, and 
processing the results through a Python API, without the need for external I/O.


## 🔧 Installing

PyPhyML is available for all modern Python versions (3.8+).


## 💡 Example

#### Build a tree from a Biopython alignment

```python
import Bio.AlignIO
from pyphyml import Alignment, TreeBuilder

# load example data from PhyML
msa = Bio.AlignIO.read("vendor/phyml/examples/nucleic", "phylip")

# create an alignment object
ali = Alignment(names=[seq.id for seq in msa], sequences=msa.alignment)

# create a tree builder, perform tree inference
builder = TreeBuilder(seed=42)
result = builder.build(ali)

# get a string with the tree in Newick format
newick = result.tree.dumps()
```

#### Use PyFAMSA to align sequence and infer a tree on-the-fly

```python
import Bio.SeqIO
import pyfamsa
import pyphyml
from scoring_matrices import ScoringMatrix

# load 16S sequences, convert from DNA to RNA
sequences = []
for record in Bio.SeqIO.parse('tests/data/16S.fna', 'fasta'):
    seq = record.seq.replace("U", "T")
    sequences.append(pyfamsa.Sequence(record.id.encode(), bytes(seq)))

# align sequences with PyFAMSA using a nucleotide scoring matrix
aligner = pyfamsa.Aligner(scoring_matrix=ScoringMatrix.from_name("NUC.4.4"))
msa = aligner.align(sequences)

# convert PyFAMSA alignment to PyPhyML alignment
ali = pyphyml.Alignment(
    names=[seq.id.decode() for seq in msa], 
    sequences=[seq.sequence.decode() for seq in msa]
)

# create tree builder and perform tree inference
builder = pyphyml.TreeBuilder()
result = builder.build(ali)
```

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue tracker](https://github.com/althonos/pyphyml/issues)
if you need to report or ask something. If you are filing in on a bug,
please include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.


<!-- ### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyphyml/blob/main/CONTRIBUTING.md)
for more details. -->


<!-- ## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pyphyml/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format. -->


## ⚖️ License

This library is provided under the [GNU General Public License 3.0 or later](https://choosealicense.com/licenses/gpl-3.0/).
The PhyML code is distributed under the [GNU General Public License 3.0](https://choosealicense.com/licenses/gpl-3.0/) as well. PhyML is developed at CNRS, consider [supporting its development](https://fondation-cnrs.org/faire-un-don/).

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the PhyML authors. It was developed by [Martin Larralde](https://github.com/althonos/) 
during his PhD project at the [Leiden University Medical Center](https://www.lumc.nl/en/) 
in the [Zeller team](https://github.com/zellerlab).*


## 📚 References

- <a id="ref1">\[1\]</a> Guindon, S., & Gascuel, O. (2003). A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood. Systematic biology, 52(5), 696–704. [doi:10.1080/10635150390235520](https://doi.org/10.1080/10635150390235520) [PMID:14530136](https://pubmed.ncbi.nlm.nih.gov/14530136/)
- <a id="ref2">\[2\]</a> Guindon, S., Dufayard, J. F., Lefort, V., Anisimova, M., Hordijk, W., & Gascuel, O. (2010). New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic biology, 59(3), 307–321. [doi:10.1093/sysbio/syq010](https://doi.org/10.1093/sysbio/syq010) [PMID:20525638](https://pubmed.ncbi.nlm.nih.gov/20525638/).

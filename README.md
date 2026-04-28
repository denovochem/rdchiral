[![PyPI version](https://badge.fury.io/py/rdchiral-plus.svg)](https://badge.fury.io/py/rdchiral-plus)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/denovochem/rdchiral_plus/graphs/commit-activity)
[![License](https://img.shields.io/pypi/l/rdchiral-plus)](https://github.com/denovochem/rdchiral_plus/blob/main/LICENSE)
[![Run Tests](https://github.com/denovochem/rdchiral_plus/actions/workflows/ci.yml/badge.svg)](https://github.com/denovochem/rdchiral_plus/actions/workflows/ci.yml)
[![Build Docs](https://github.com/denovochem/rdchiral_plus/actions/workflows/docs.yml/badge.svg)](https://github.com/denovochem/rdchiral_plus/actions/workflows/docs.yml)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/denovochem/rdchiral_plus/blob/main/examples/example_notebook.ipynb)

# rdchiral_plus
Wrapper for RDKit's RunReactants to improve stereochemistry handling

This repository is a fork of [rdchiral](https://github.com/connorcoley/rdchiral). It has been modified for improved performance, and is statically typed wherever possible so that it can be optionally compiled with [mypyc](https://mypyc.readthedocs.io/en/latest/introduction.html) for faster execution while maintaining consistency with the upstream library. These modifications provide comparable speed to the fast C++ version ([rdchiral_cpp](https://gitlab.com/ljn917/rdchiral_cpp)), with all of the benefits of being written in Python. This library is pip installable cross platform, and can be used as a drop-in replacement for the original rdchiral library.

The interface (`rdchiralRun`, `rdchiralRunText`, `rdchiralReaction`, `rdchiralReactants`, `rdchiralExtract`, etc.) remains largely unchanged from the original library, so existing code should work with minimal or no modifications. While behavior is mostly consistent with the original library, this fork includes several important fixes and improvements:

## Template application

- **Conjugated system bond direction correction**: Fixes corrupted single-bond directions (ENDUPRIGHT, ENDDOWNRIGHT) in conjugated systems. See [here](https://github.com/connorcoley/rdchiral/pull/40)
- **Broader stereochemistry handling**: Stereochemistry for tetrahedral centers with lone pairs is now properly handled
- **One-pot reactions**: Templates are initialized with parentheses where needed so that templates defining multiple reactions on the same product are properly handled
- **Recursive template application**: Templates can be recursively applied, useful for symmetric reactions

## Template extraction

- **Configurable template extraction**: Template extraction now supports configurable radius and special group handling. See [here](https://github.com/connorcoley/rdchiral/commit/78bbafaba040678b957497e7f2638e935104e3d7)
- **Stereochemistry tracking**: Inversions of tetrahedral centers are now counted as a changed atom, and included in the extracted template
- **Spectator tracking**: Spectator molecules are included in extracted template dictionaries

## General

- **Automatic dependency installation**: RDKit is now automatically installed as a dependency

Consistency with the original library:
- rdchiralRun (1000 templates applied to 1000 reactants): 99.98% consistent
- rdchiralRunText (1000 templates applied to 100 reactants): 99.97% consistent
- rdchiralExtract (templates extracted from 50,016 mapped reactions): 94.99% consistent

See [here](docs/consistency.md) for details on how consistency is measured against the original library and full details of what changes you can expect compared to the original rdchiral library.

## Requirements

- RDKit (version >= 2019)
- Python (version >= 3.10)

## Installation

Install rdchiral_plus from PyPI:

```bash
pip install rdchiral-plus
```

Or install rdchiral_plus with pip directly from this repo:

```bash
pip install git+https://github.com/denovochem/rdchiral_plus.git
```

For mypyc compilation:

```bash
RDCHIRAL_USE_MYPYC=1 pip install "git+https://github.com/denovochem/rdchiral_plus.git"
```


## Basic usage
```python
from rdchiral import rdchiralRunText, rdchiralReaction, rdchiralReactants

# Run directly from SMARTS and SMILES strings
# This is slower than pre-initializing rdchiralReaction and rdchiralReactants when
# processing a large number of reactions
reaction_smarts = '[C:1][OH:2]>>[C:1][O:2][C]'
reactant_smiles = 'OCC(=O)OCCCO'
outcomes = rdchiralRunText(reaction_smarts, reactant_smiles)
print(outcomes)

# Pre-initialize then run
rxn = rdchiralReaction(reaction_smarts)
reactants = rdchiralReactants(reactant_smiles)
outcomes = rdchiralRun(rxn, reactants)
print(outcomes)

# Get list of atoms that changed
outcomes, mapped_outcomes = rdchiralRun(rxn, reactants, return_mapped=True)
print(outcomes, mapped_outcomes)
```

## Documentation
Full documentation is available [here](https://denovochem.github.io/rdchiral_plus/)

## Contributing

- Feature ideas and bug reports are welcome on the Issue Tracker.
- Fork the [source code](https://github.com/denovochem/rdchiral_plus) on GitHub, make changes and file a pull request.

## License

rdchiral_plus is licensed under the [MIT license](https://github.com/denovochem/rdchiral_plus/blob/main/LICENSE).

## References

- [Original rdchiral paper](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00286)
- [rdchiral repo](https://github.com/connorcoley/rdchiral)
- [rdchiral_cpp repo](https://gitlab.com/ljn917/rdchiral_cpp)

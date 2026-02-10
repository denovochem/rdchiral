import functools
from typing import List

from rdkit import Chem


@functools.lru_cache(maxsize=128)
def mol_from_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return mol


@functools.lru_cache(maxsize=128)
def get_mol_atoms(mol: Chem.Mol) -> List[Chem.Atom]:
    return list(mol.GetAtoms())


@functools.lru_cache(maxsize=128)
def get_mol_bonds(mol: Chem.Mol) -> List[Chem.Bond]:
    return list(mol.GetBonds())

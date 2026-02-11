import functools

from rdkit import Chem


@functools.lru_cache(maxsize=128)
def mol_from_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return mol

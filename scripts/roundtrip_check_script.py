"""Roundtrip consistency check script - run via run_roundtrip_check.py."""

import os
from pathlib import Path

from rdkit import Chem

from rdchiral.main import rdchiralRunText

# Try fork-specific API first, fall back to original API
try:
    from rdchiral.template_extractor import extract_from_reaction_smiles

    def extract(rxn):
        return extract_from_reaction_smiles(rxn)
except ImportError:
    from rdchiral.template_extractor import extract_from_reaction

    def extract(rxn):
        reactants_side = rxn.split(">>")[0]
        products_side = rxn.split(">>")[1]
        return extract_from_reaction(
            {
                "reactants": reactants_side,
                "products": products_side,
                "_id": 0,
            }
        )


_script_dir = Path(__file__).resolve().parent
_default_repo_root = _script_dir.parent
_env_root = Path(os.environ.get("RDCHIRAL_REPO_ROOT", _default_repo_root))

# Resolve data file path
if _env_root.name == "scripts" and (_env_root.parent / "rdchiral").exists():
    _data_root = _env_root
else:
    _data_root = _env_root / "scripts"

REACTIONS_PATH = _data_root / "uspto_50k_mapped_reactions.txt"


def main():
    with open(REACTIONS_PATH, "r") as f:
        rxns = f.read().splitlines()

    consistent = 0
    for rxn in rxns:
        try:
            product_string = rxn.split(">>")[1]
            if "." in product_string:
                continue  ## Skip products with multiple fragments
            product_mol = Chem.MolFromSmiles(product_string)
            [atom.SetAtomMapNum(0) for atom in product_mol.GetAtoms()]
            product_string_no_mapping = Chem.MolToSmiles(product_mol)

            reactant_string = rxn.split(">>")[0]
            reactant_mol = Chem.MolFromSmiles(reactant_string)
            [atom.SetAtomMapNum(0) for atom in reactant_mol.GetAtoms()]
            reactant_string_no_mapping = Chem.MolToSmiles(reactant_mol)

            out = extract(rxn)
            rxn_smarts = out.get("reaction_smarts", "")
            if not rxn_smarts:
                continue

            rdchiral_generated_reactants = rdchiralRunText(
                rxn_smarts, product_string_no_mapping, combine_enantiomers=False
            )

            for ele in rdchiral_generated_reactants:
                if ele in reactant_string_no_mapping:
                    consistent += 1
                    break

        except Exception:
            pass

    print(consistent)


if __name__ == "__main__":
    main()

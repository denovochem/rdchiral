from typing import Dict, List, Tuple

import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem import rdChemReactions, rdmolops
from rdkit.Chem.rdchem import BondDir, BondStereo, ChiralType

from rdchiral.bonds import (
    bond_dirs_by_mapnum,
    enumerate_possible_cistrans_defs,
    get_atoms_across_double_bonds,
)
from rdchiral.chiral import template_atom_could_have_been_tetra
from rdchiral.function_cache import get_mol_atoms, get_mol_bonds, mol_from_smiles
from rdchiral.utils import PLEVEL

BondDirOpposite = {
    AllChem.BondDir.ENDUPRIGHT: AllChem.BondDir.ENDDOWNRIGHT,
    AllChem.BondDir.ENDDOWNRIGHT: AllChem.BondDir.ENDUPRIGHT,
}


class rdchiralReaction(object):
    """Class to store everything that should be pre-computed for a reaction. This
    makes library application much faster, since we can pre-do a lot of work
    instead of doing it for every mol-template pair

    Attributes:
        reaction_smarts (str): reaction SMARTS string
        rxn (rdkit.Chem.rdChemReactions.ChemicalReaction): RDKit reaction object.
            Generated from `reaction_smarts` using `initialize_rxn_from_smarts`
        template_r: Reaction reactant template fragments
        template_p: Reaction product template fragments
        atoms_rt_map (dict): Dictionary mapping from atom map number to RDKit Atom for reactants
        atoms_pt_map (dict): Dictionary mapping from atom map number to RDKit Atom for products
        atoms_rt_idx_to_map (dict): Dictionary mapping from atom idx to RDKit Atom for reactants
        atoms_pt_idx_to_map (dict): Dictionary mapping from atom idx to RDKit Atom for products

    Args:
        reaction_smarts (str): Reaction SMARTS string
    """

    def __init__(self, reaction_smarts: str):
        # Keep smarts, useful for reporting
        self.reaction_smarts = reaction_smarts

        # Initialize - assigns stereochemistry and fills in missing rct map numbers
        self.rxn: rdChemReactions.ChemicalReaction = initialize_rxn_from_smarts(
            reaction_smarts
        )

        # Combine template fragments so we can play around with mapnums
        self.template_r, self.template_p = get_template_frags_from_rxn(self.rxn)

        # Define molAtomMapNumber->atom dictionary for template rct and prd
        self.atoms_rt_map: Dict[int, Chem.Atom] = {
            a.GetAtomMapNum(): a
            for a in get_mol_atoms(self.template_r)
            if a.GetAtomMapNum()
        }
        self.atoms_pt_map: Dict[int, Chem.Atom] = {
            a.GetAtomMapNum(): a
            for a in get_mol_atoms(self.template_p)
            if a.GetAtomMapNum()
        }

        # Back-up the mapping for the reaction
        self.atoms_rt_idx_to_map: Dict[int, int] = {
            a.GetIdx(): a.GetAtomMapNum() for a in get_mol_atoms(self.template_r)
        }
        self.atoms_pt_idx_to_map: Dict[int, int] = {
            a.GetIdx(): a.GetAtomMapNum() for a in get_mol_atoms(self.template_p)
        }

        # Check consistency (this should not be necessary...)
        if any(
            self.atoms_rt_map[i].GetAtomicNum() != self.atoms_pt_map[i].GetAtomicNum()
            for i in self.atoms_rt_map
            if i in self.atoms_pt_map
        ):
            raise ValueError("Atomic identity should not change in a reaction!")

        # Call template_atom_could_have_been_tetra to pre-assign value to atom
        for a in get_mol_atoms(self.template_r):
            template_atom_could_have_been_tetra(a)
        for a in get_mol_atoms(self.template_p):
            template_atom_could_have_been_tetra(a)

        # Pre-list chiral double bonds (for copying back into outcomes/matching)
        self.rt_bond_dirs_by_mapnum = bond_dirs_by_mapnum(self.template_r)
        self.pt_bond_dirs_by_mapnum = bond_dirs_by_mapnum(self.template_p)

        # Enumerate possible cis/trans...
        self.required_rt_bond_defs, self.required_bond_defs_coreatoms = (
            enumerate_possible_cistrans_defs(self.template_r)
        )

    def reset(self):
        """Reset atom map numbers for template fragment atoms"""
        for idx, mapnum in self.atoms_rt_idx_to_map.items():
            self.template_r.GetAtomWithIdx(idx).SetAtomMapNum(mapnum)
        for idx, mapnum in self.atoms_pt_idx_to_map.items():
            self.template_p.GetAtomWithIdx(idx).SetAtomMapNum(mapnum)


class rdchiralReactants(object):
    """Class to store everything that should be pre-computed for a reactant mol
    so that library application is faster

    Attributes:
        reactant_smiles (str): Reactant SMILES string
        reactants (rdkit.Chem.rdchem.Mol): RDKit Molecule create from `initialize_reactants_from_smiles`
        atoms_r (dict): Dictionary mapping from atom map number to atom in `reactants` Molecule
        idx_to_mapnum (callable): callable function that takes idx and returns atom map number
        reactants_achiral (rdkit.Chem.rdchem.Mol): achiral version of `reactants`
        bonds_by_mapnum (list): List of reactant bonds
            (int, int, rdkit.Chem.rdchem.Bond)
        bond_dirs_by_mapnum (dict): Dictionary mapping from atom map number tuples to BondDir
        atoms_across_double_bonds (list): List of cis/trans specifications from `get_atoms_across_double_bonds`

    Args:
        reactant_smiles (str): Reactant SMILES string
        custom_reactant_mapping (bool): Whether to use custom reactant mapping
    """

    def __init__(self, reactant_smiles: str, custom_reactant_mapping: bool = False):
        # Keep original smiles, useful for reporting
        self.reactant_smiles = reactant_smiles
        self.custom_mapping = custom_reactant_mapping

        # Initialize into RDKit mol
        self.reactants = initialize_reactants_from_smiles(
            reactant_smiles, custom_reactant_mapping
        )

        # Set mapnum->atom dictionary
        # all reactant atoms must be mapped after initialization, so this is safe

        self.atoms_r: Dict[int, Chem.Atom] = {}
        for atom in get_mol_atoms(self.reactants):
            self.atoms_r[atom.GetAtomMapNum()] = atom

        self.idx_to_mapnum = lambda idx: self.reactants.GetAtomWithIdx(
            idx
        ).GetAtomMapNum()

        # Create copy of molecule without chiral information, used with
        # RDKit's naive runReactants
        self.reactants_achiral = Chem.Mol(self.reactants)
        for a in get_mol_atoms(self.reactants_achiral):
            a.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
        for b in get_mol_bonds(self.reactants_achiral):
            b.SetStereo(BondStereo.STEREONONE)
            b.SetBondDir(BondDir.NONE)

        # Pre-list reactant bonds (for stitching broken products)
        self.bonds_by_mapnum: List[Tuple[int, int, Chem.Bond]] = [
            (b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum(), b)
            for b in get_mol_bonds(self.reactants)
        ]

        # Pre-list chiral double bonds (for copying back into outcomes/matching)
        self.bond_dirs_by_mapnum: Dict[Tuple[int, int], BondDir] = {}
        for i, j, b in self.bonds_by_mapnum:
            if b.GetBondDir() != BondDir.NONE:
                self.bond_dirs_by_mapnum[(i, j)] = b.GetBondDir()
                self.bond_dirs_by_mapnum[(j, i)] = BondDirOpposite[b.GetBondDir()]

        # Get atoms across double bonds defined by mapnum
        self.atoms_across_double_bonds = get_atoms_across_double_bonds(self.reactants)


def initialize_rxn_from_smarts(
    reaction_smarts: str,
) -> rdChemReactions.ChemicalReaction:
    """Initialize RDKit reaction object from SMARTS string

    Args:
        reaction_smarts (str): Reaction SMARTS string

    Returns:
        rdkit.Chem.rdChemReactions.ChemicalReaction: RDKit reaction object
    """
    # Initialize reaction
    rxn = rdChemReactions.ReactionFromSmarts(reaction_smarts)
    rxn.Initialize()
    if rxn.Validate()[1] != 0:
        raise ValueError("validation failed")
    if PLEVEL >= 2:
        print("Validated rxn without errors")

    # Figure out if there are unnecessary atom map numbers (that are not balanced)
    # e.g., leaving groups for retrosynthetic templates. This is because additional
    # atom map numbers in the input SMARTS template may conflict with the atom map
    # numbers of the molecules themselves
    prd_maps: List[int] = [
        a.GetAtomMapNum()
        for prd in rxn.GetProducts()
        for a in get_mol_atoms(prd)
        if a.GetAtomMapNum()
    ]

    unmapped = 700
    for rct in rxn.GetReactants():
        rct.UpdatePropertyCache(strict=False)
        Chem.AssignStereochemistry(rct)
        # Fill in atom map numbers
        for a in get_mol_atoms(rct):
            if not a.GetAtomMapNum() or a.GetAtomMapNum() not in prd_maps:
                a.SetAtomMapNum(unmapped)
                unmapped += 1
    if PLEVEL >= 2:
        print("Added {} map nums to unmapped reactants".format(unmapped - 700))
    if unmapped > 800:
        raise ValueError(
            "Why do you have so many unmapped atoms in the template reactants?"
        )

    return rxn


def initialize_reactants_from_smiles(
    reactant_smiles: str, custom_reactant_mapping: bool
) -> Chem.Mol:
    """Initialize RDKit molecule from SMILES string

    Args:
        reactant_smiles (str): Reactant SMILES string

    Returns:
        rdkit.Chem.rdchem.Mol: RDKit molecule
    """
    # Initialize reactants
    reactants = mol_from_smiles(reactant_smiles)
    Chem.AssignStereochemistry(reactants, flagPossibleStereoCenters=True)
    reactants.UpdatePropertyCache(strict=False)
    # To have the product atoms match reactant atoms, we
    # need to populate the map number field, since this field
    # gets copied over during the reaction via reactant_atom_idx.
    if not custom_reactant_mapping:
        for i, a in enumerate(get_mol_atoms(reactants)):
            a.SetAtomMapNum(i + 1)
    if PLEVEL >= 2:
        print(
            "Initialized reactants, assigned map numbers, stereochem, flagpossiblestereocenters"
        )
    return reactants


def get_template_frags_from_rxn(
    rxn: rdChemReactions.ChemicalReaction,
) -> Tuple[Chem.Mol, Chem.Mol]:
    """Get template fragments from RDKit reaction object

    Args:
        rxn (rdkit.Chem.rdChemReactions.ChemicalReaction): RDKit reaction object

    Returns:
        (rdkit.Chem.rdchem.Mol, rdkit.Chem.rdchem.Mol): tuple of fragment molecules
    """
    # Copy reaction template so we can play around with map numbers
    for i, rct in enumerate(rxn.GetReactants()):
        if i == 0:
            template_r = rct
        else:
            template_r = rdmolops.CombineMols(template_r, rct)
    for i, prd in enumerate(rxn.GetProducts()):
        if i == 0:
            template_p = prd
        else:
            template_p = rdmolops.CombineMols(template_p, prd)
    return template_r, template_p

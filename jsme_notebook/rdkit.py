from rdkit import Chem
from . import JSMENotebook
from typing import Optional

import sys

if sys.version_info.minor <= 7:
    # backport of functools.singledispatchmethod to python <= 3.7
    from singledispatchmethod import singledispatchmethod
else:
    from functools import singledispatchmethod

class JSMERDKit(JSMENotebook):
    # yes, I known what mer+rd spells.
    """
    Subclass of JSMEJSMENotebook that uses RDKit to perform chemical operations.
    Initialisation accepts a smiles, mol block, rdkit.Chem.Mol, or None.
    The (read-only) property mol gives a rdkit.Chem.Mol object.
    """

    @singledispatchmethod
    def __init__(self, smiles: Optional[str]=None, cdn_url: Optional[str]=None):
        """
        Single dispatch method which can accept a smiles string or a rdkit mol.
        :param smiles:
        """
        if smiles is None:
            pass
        elif smiles.count('\n'):  # molblock
            mol_block = smiles
            smiles=Chem.MolToSmiles(Chem.MolFromMolBlock(mol_block))
        super().__init__(smiles, cdn_url=cdn_url)

    @__init__.register
    def _(self, mol: Chem.Mol, cdn_url: Optional[str]=None):
        super().__init__(Chem.MolToSmiles(mol), cdn_url=cdn_url)

    @property
    def mol(self) -> Chem.Mol:
        return Chem.MolFromSmiles(self.smiles)

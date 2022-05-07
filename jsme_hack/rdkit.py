from rdkit import Chem
from . import JSMEHack
from typing import Optional

import sys

if sys.version_info.minor <= 7:
    # backport of functools.singledispatchmethod to python <= 3.7
    from singledispatchmethod import singledispatchmethod
else:
    from functools import singledispatchmethod

class JSMERdkit(JSMEHack):
    # yes, I known what mer+rd spells.
    @singledispatchmethod
    def __init__(self, smiles: Optional[str]=None):
        if smiles is None:
            pass
        elif smiles.count('\n'):  # molblock
            mol_block = smiles
            smiles=Chem.MolToSmiles(Chem.MolFromMolBlock(mol_block))
        super().__init__(smiles)

    @__init__.register
    def _(self, mol: Chem.Mol):
        super().__init__(Chem.MolToSmiles(mol))
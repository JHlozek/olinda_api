from abc import ABC
from typing import Any, List

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

class Featurizer(ABC):
    def featurize(self: "Featurizer", batch: Any) -> Any:
        """Featurize input batch.

        Args:
            batch (Any): batch of smiles

        Returns:
            Any: featurized outputs
        """
        pass

class MorganFeaturizer(Featurizer):
    def __init__(self: "MorganFeaturizer") -> None:
        self.name = "morganfeaturizer"
        
    def featurize(self: "MorganFeaturizer", batch: Any) -> Any:
        """Featurize input batch.

        Args:
            batch (Any): batch of smiles

        Returns:
            Any: featurized outputs
        """
        mols = [Chem.MolFromSmiles(smi) for smi in batch if smi is not None]
        fps = [np.array(AllChem.GetMorganFingerprintAsBitVect(
                mol, radius=3, nBits=1024))
                if mol is not None else None
                for mol in mols
        ]
        return fps

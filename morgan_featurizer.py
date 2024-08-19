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
        mols = [Chem.MolFromSmiles(smi) for smi in batch if smi is not np.nan]
        ecfps = self.ecfp_counts(mols)
        return ecfps
    
    def ecfp_counts(self: "MorganFeaturizer", mols: List) -> List:
        """Create ECFPs from batch of smiles.

        Args:
            mols (List): batch of molecules

        Returns:
            List: batch of ECFPs
        """
        fps = [
            AllChem.GetMorganFingerprint(
                mol, radius=3, useCounts=True, useFeatures=True
            ) if mol is not None else None
            for mol in mols
        ]
        nfp = np.zeros((len(fps), 1024), np.float32)
        for i, fp in enumerate(fps):
            if fp is not None:
                for idx, v in fp.GetNonzeroElements().items():
                    nidx = idx % 1024
                    nfp[i, nidx] += int(v)
        return nfp
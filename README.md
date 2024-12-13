# Molecular Descriptor Script

This script calculates various molecular descriptors and provides information about a given compound based on its structure. It supports multiple file formats, including SMILES strings, MOL, SDF, MOL2, and PDB files. The descriptors include molecular weight, LogP, TPSA, QED, hydrogen bond donors and acceptors, rotatable bonds, and alerts for PAINS, Brenk's, and NIH.

## Features:
- Calculates the following descriptors:
  - **Molecular Weight** (MW)
  - **LogP** (Octanol-Water Partition Coefficient)
  - **TPSA** (Topological Polar Surface Area)
  - **QED** (Quantitative Estimation of Drug-likeness)
  - **Number of Hydrogen Bond Donors** (HBD)
  - **Number of Hydrogen Bond Acceptors** (HBA)
  - **Number of Rotatable Bonds**
  - **Molar Refractivity** (MR)
  - **Heavy Atom Count**
  - **Fraction of sp3 Carbons**
  - **Lipinski's Rule Violations**
  - **Veber's Criteria Violations**
  - **PAINS Alerts**
  - **Brenk's Alerts**
  - **NIH Alerts**

## Requirements:
- RDKit
- Pandas
- argparse
- tabulate

You can install the required libraries by running:

```bash
pip install rdkit pandas tabulate
```

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

## Usage:

To run the script, you need to pass either a SMILES string or a molecular file (in formats like `.mol`, `.sdf`, `.mol2`, `.pdb`, or `.smiles`).

### Running the Script with a SMILES String:
To get the molecular descriptor details of a compound from its SMILES string, use the following command:

```bash
python3 MolecularDescriptor.py 'SMILES_STRING'
```

Example:
```bash
python3 MolecularDescriptor.py 'CC(C)CC(C(=O)O)N'
```

### Running the Script with a Molecular File:
You can also pass molecular files in formats like `.sdf`, `.mol`, `.mol2`, `.pdb`, or `.smiles`. For example:

```bash
python3 MolecularDescriptor.py 'compound.sdf'
```
```bash
python3 MolecularDescriptor.py 'compound.mol2'
```
```bash
python3 MolecularDescriptor.py 'compound.pdb'
```

### Writing Results to a CSV File:
If you want to write the results to a CSV file, use the `--write` flag:

```bash
python3 MolecularDescriptor.py 'SMILES_STRING' --write
```
```bash
python3 MolecularDescriptor.py 'compound.sdf' --write
```
```bash
python3 MolecularDescriptor.py 'compound.mol2' --write
```

This will save the compound details in a CSV file (`compound_description.csv`).

## Output:
The output will display a table with the following molecular descriptors for the input compound:

- **SMILES**: SMILES string of the compound
- **Molecular Weight**: Molecular weight of the compound
- **LogP**: LogP value (octanol-water partition coefficient)
- **TPSA**: Topological Polar Surface Area
- **QED**: Quantitative Estimation of Drug-likeness
- **Number of H bond donors**: Number of hydrogen bond donors
- **Number of H bond acceptors**: Number of hydrogen bond acceptors
- **Number of rotatable bonds**: Number of rotatable bonds
- **Molar refractivity**: Molar refractivity of the compound
- **Heavy atom count**: Number of heavy atoms
- **Fraction of sp3 carbons**: Fraction of sp3 hybridized carbon atoms
- **Lipinski's violations**: Number of violations of Lipinski's rule of five
- **Veber's criteria**: Number of violations of Veber's criteria
- **Pains alert**: Whether the compound matches PAINS alerts (1 = matches, 0 = does not match)
- **Brenk's alert**: Whether the compound matches Brenk's alerts (1 = matches, 0 = does not match)
- **NIH alert**: Whether the compound matches NIH alerts (1 = matches, 0 = does not match)

### Example Table Output:
```
+---------------------------+------------------+
| Descriptor                 | Value            |
+---------------------------+------------------+
| SMILES                     | CC(C)CC(C(=O)O)N |
| Molecular Weight           | 173.2            |
| LogP                       | 1.45             |
| TPSA                       | 48.2             |
| QED                        | 0.85             |
| Number of H bond donors    | 1                |
| Number of H bond acceptors | 2                |
| Number of rotatable bonds  | 4                |
| Molar refractivity         | 75.3             |
| Heavy atom count           | 10               |
| Fraction of sp3 carbons    | 0.9              |
| Lipinski's violations      | 0                |
| Veber's criteria           | 0                |
| Pains alert                | 0                |
| Brenk's alert              | 0                |
| NIH alert                  | 0                |
+---------------------------+------------------+
```

This table will be displayed in the terminal when you run the script. If the `--write` flag is used, it will also be saved to a CSV file.


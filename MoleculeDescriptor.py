from rdkit import Chem
from rdkit.Chem import AllChem,Descriptors,Lipinski
from rdkit.Chem.FilterCatalog import FilterCatalog,FilterCatalogParams
import pandas as pd
import argparse
from tabulate import tabulate

def read_molecule(m:str):
    if not m.endswith(('.sdf','mol','mol2','.pdb','.smi','.pdb')):
        mol=Chem.MolFromSmiles(m)
        if mol is not None:
            return mol
        else:
            raise ValueError(f'Invalid SMILES string: {m}')
    if m.endswith('.sdf') or m.endswith('.mol'):
        mol=Chem.MolFromMolFile(m)
    elif m.endswith('.mol2'):
        mol=Chem.MolFromMol2File(m)
    elif m.endswith('.smiles') or m.endswith('smi'):
        with open(m,'r') as f:
            smiles=f.readline().strip()
            mol=Chem.MolFromSmiles(smiles)
    elif m.endswith('.pdb'):
        mol=Chem.MolFromPDBFile(m)
    else:
        raise ValueError(f'Unsupported file format: {m}')
    
    return mol

def LipinskiViolations(mol):
    mw=Descriptors.MolWt(mol)
    hba=Descriptors.NOCount(mol)
    hbd=Descriptors.NHOHCount(mol)
    logp=Descriptors.MolLogP(mol)
    conditions=[(mw<=500),(hba<=10),(hbd<=5),(logp<=5)]
    violations=sum([1 for condition in conditions if not condition])
    return violations

def VeberCriteria(mol):
    tpsa=Descriptors.TPSA(mol)
    rot=Descriptors.NumRotatableBonds(mol)
    conditions=[(tpsa<=140),(rot<=10)]
    violations=sum([1 for condition in conditions if not condition])
    return violations

def Pains(mol):
    params_pains=FilterCatalogParams()
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    catalog_pains=FilterCatalog(params_pains)
    flag=catalog_pains.HasMatch(mol)
    if flag:
        return 1
    return 0

def Brenks(mol):
    params_unwanted=FilterCatalogParams()
    params_unwanted.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalog_unwanted=FilterCatalog(params_unwanted)
    flag=catalog_unwanted.HasMatch(mol)
    if flag:
        return 1
    return 0

def Nih(mol):
    params_nih=FilterCatalogParams()
    params_nih.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
    catalog_nih=FilterCatalog(params_nih)
    flag=catalog_nih.HasMatch(mol)
    if flag:
        return 1
    return 0

def all_descriptors(m):
    mol = read_molecule(m)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    qed = Descriptors.qed(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    mr = Descriptors.MolMR(mol)
    heavy_atoms = Descriptors.HeavyAtomCount(mol)
    sp3_fraction = Descriptors.FractionCSP3(mol)
    lipinski_violations = LipinskiViolations(mol)
    veber_violations = VeberCriteria(mol)
    pains=Pains(mol)
    brenk=Brenks(mol)
    nih=Nih(mol)
    return {
        "SMILES": Chem.MolToSmiles(mol),
        "Molecular Weight": mw,
        "LogP": logp,
        "TPSA": tpsa,
        "QED": qed,
        "Number of H bond donors": hbd,
        "Number of H bond acceptors": hba,
        "Number of rotatable bonds": rot_bonds,
        "Molar refractivity": mr,
        "Heavy atom count": heavy_atoms,
        "Fraction of sp3 carbons": sp3_fraction,
        "Lipinski's violations": lipinski_violations,
        "Veber's criteria": veber_violations,
        "Pains alert": pains,
        "Brenk's alert": brenk,
        "NIH alert": nih
    }
    
def get_compound_details(m,write=False):
    info=all_descriptors(m)
    df=pd.DataFrame(list(info.items()),columns=['Descriptor','Value'])
    table=tabulate(df,headers='keys',tablefmt='fancy_grid')
    print(table)
    if write:
        df.to_csv('compound_description.csv',index=False)

def main():
    parser=argparse.ArgumentParser(description='Get details of a compound')
    parser.add_argument('m',type=str,help='Molecule file (mol,sdf,mol2,smiles or pdb format) of smiles string is passed')
    parser.add_argument('--write',action='store_true',help='Boolean flag to write the details to a csv file')
    args=parser.parse_args()

    get_compound_details(args.m,args.write)

if __name__=='__main__':
    main()

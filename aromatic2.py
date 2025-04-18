import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from rdkit.Chem import Draw

st.set_page_config(page_title="Molecule Aromaticity Checker", layout="centered")

st.title("🧪 Molecule Aromaticity Checker")
st.markdown(
    "Draw a molecule below and find out whether it is **aromatic**, **anti-aromatic**, "
    "or **non-aromatic** based on Hückel's Rule."
)

# -------------------- Classification Logic -------------------- #
def classify_molecule(mol):
    
    ri = mol.GetRingInfo()
    if ri.NumRings() == 0:
        return "Non-aromatic"

    aromatic_rings     = 0
    anti_aromatic_rings = 0

    for ring in ri.AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]

        # 1) Planarity assumption: skip anti‐aromatic check for large rings
        is_planar = len(ring_atoms) <= 7

        # 2) Check full conjugation of bonds
        bonds = [
            b for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in ring and b.GetEndAtomIdx() in ring
        ]
        is_fully_conjugated = all(b.GetIsConjugated() for b in bonds)

        # 3) Check RDKit’s own aromatic flag
        if all(a.GetIsAromatic() for a in ring_atoms):
            aromatic_rings += 1
            continue

        # 4) Anti‐aromatic candidate requires planarity + full conjugation
        if not (is_planar and is_fully_conjugated):
            continue

        # 5) Count π‑electrons in the ring
        pi_electrons = 0
        for atom in ring_atoms:
            hyb = atom.GetHybridization().name
            # only sp2 or sp can contribute one p‑orbital
            if hyb in ("SP2", "SP"):
                pi_electrons += 1
            # lone‐pair donors (pyrrole, furan, thiophene style)
            elif atom.GetSymbol() in ("N","O","S") and atom.GetTotalNumHs()==1:
                pi_electrons += 2

        # 6) Apply 4n rule for anti‑aromaticity
        if pi_electrons % 4 == 0 and pi_electrons > 0:
            anti_aromatic_rings += 1

    if aromatic_rings > 0:
        return "Aromatic"
    if anti_aromatic_rings > 0:
        return "Anti-aromatic"
    return "Non-aromatic"


# -------------------- Molecule Parser -------------------- #
def get_sanitized_mol(smiles):
    try:
        m = Chem.MolFromSmiles(smiles)
        if m: Chem.SanitizeMol(m)
        return m
    except:
        return None

# -------------------- SMILES Fixer -------------------- #
def fix_smiles(smiles):
    fixes = {
        "N1C=CC=C1": "c1cc[nH]c1",  # pyrrole
        "C1=CC=CC=C1": "c1ccccc1",  # benzene
        # add more if you spot other kekulé quirks...
    }
    return fixes.get(smiles, smiles)


# -------------------- Streamlit UI -------------------- #
smiles = st_ketcher()
if smiles:
    smiles = fix_smiles(smiles)
    st.markdown("🔬 **SMILES Representation:**")
    st.code(smiles)

    mol = get_sanitized_mol(smiles)
    if mol is None:
        st.error("Invalid molecule. Please check the structure.")
    else:
        classification = classify_molecule(mol)
        st.image(Draw.MolToImage(mol, size=(300,300)), caption="Molecule Structure")

        if classification == "Aromatic":
            st.success("# This molecule is **aromatic** 🌟")
        elif classification == "Anti-aromatic":
            st.warning("# This molecule is **anti-aromatic** ⚠️")
        else:
            st.info("# This molecule is **non-aromatic** 🧊")
else:
    st.info("👈 Use the editor above to draw a molecule.")

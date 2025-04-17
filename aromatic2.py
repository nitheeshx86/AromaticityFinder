import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit import Chem

st.set_page_config(page_title="Molecule Aromaticity Checker", layout="centered")

st.title("🧪 Molecule Aromaticity Checker")
st.markdown("Draw a molecule below and find out whether it is **aromatic**, **anti-aromatic**, or **non-aromatic**.")

# Molecule drawing canvas
smiles = st_ketcher()
if smiles:
    st.markdown("### 🔬 SMILES Representation:")
    st.code(smiles)

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.error("Invalid molecule. Please check the structure.")
        else:
            aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() 
                              if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]

            num_aromatic_rings = len(aromatic_rings)

            if num_aromatic_rings > 0:
                st.success("This molecule is **aromatic** 🌟")
            else:
                # Check if it's cyclic
                ri = mol.GetRingInfo()
                if ri.NumRings() > 0:
                    st.warning("This molecule is **anti-aromatic** ⚠️ ")
                else:
                    st.info("This molecule is **non-aromatic** 🧊 ")
    except Exception as e:
        st.error(f"Something went wrong: {e}")
else:
    st.info("👈 Use the editor above to draw a molecule.")

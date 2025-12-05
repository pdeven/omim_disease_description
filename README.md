**OMIM–MedGen JSON Builder**

This tool converts MedGen/OMIM biomedical resources into a structured JSON file for downstream annotation or research.

---

**Input Files**

**1. MedGen_HPO_OMIM_Mapping.txt.gz**  
- Contains OMIM ↔ MedGen ↔ HPO relationships  
- First line begins with `#`  
- Delimiter may be **pipe**, **tab**, or **comma**

**2. MGDEF.csv.gz**  
- Contains concept definitions  
- Columns: `CUI`, `DEF`, `source`, `SUPPRESS`  
- Delimiter may be **tab** or **comma**

---

**Output Format**

Each output entry is a dictionary:

{
"_id": "144750",
"omim_id": "144750",
"omim_disease": "ENDOSTEAL HYPEROSTOSIS, AUTOSOMAL DOMINANT",
"medgen_concept_id": "C0432273",
"medgen_disease_info": "Definition text from MGDEF mapped by CUI"
}

---

**What the Script Does**

- Reads MGDEF and builds a `{CUI → DEF}` lookup  
- Reads the MedGen mapping file and extracts OMIM/HPO relationships  
- Ensures unique OMIM_CUI entries  
- Fetches definition text from MGDEF  
- Falls back to `"NA"` if no definition exists  

---

**Requirements**

- Python 3.8+  
- pandas  



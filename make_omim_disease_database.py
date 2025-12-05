"""
Build list-of-dicts JSON from:
 - MedGen_HPO_OMIM_Mapping.txt.gz  (header line starts with '#', pipe or tab delimited)
 - MGDEF.csv.gz                     (has columns CUI, DEF, source, SUPPRESS; delimiter may be tab or comma)

Output: JSON list-of-dicts with fields:
  _id, omim_id, omim_disease, medgen_concept_id, medgen_disease_info
"""

import gzip
import json
import os
from io import StringIO
from typing import Dict, List

import pandas as pd


def read_text_gz_all(path: str) -> List[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        return fh.read().splitlines()


def detect_delim_from_header(header_line: str) -> str:
    """Detect delimiter from a header line: prefer tab, then comma, then pipe."""
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    if "|" in header_line:
        return "|"
    # fallback - default to tab if uncertain
    return "\t"


def read_mgdef_map(mgdef_path: str) -> Dict[str, str]:
    """
    Read MGDEF.csv.gz and return mapping {CUI: DEF}.
    Detects delimiter from header line and expects columns named 'CUI' and 'DEF'.
    """
    if not os.path.exists(mgdef_path):
        raise FileNotFoundError(f"MGDEF file not found: {mgdef_path}")

    lines = read_text_gz_all(mgdef_path)
    # find first non-empty header line
    idx = 0
    while idx < len(lines) and lines[idx].strip() == "":
        idx += 1
    if idx >= len(lines):
        raise ValueError("MGDEF file appears empty")

    header_line = lines[idx].strip()
    delim = detect_delim_from_header(header_line)
    # Build a CSV text starting at header_line
    data_text = "\n".join(lines[idx:])  # include header and rest
    df = pd.read_csv(
        StringIO(data_text),
        sep=delim,
        header=0,
        dtype=str,
        keep_default_na=False,
        engine="python",
    )

    # Normalize column names case
    cols_lower = {c.lower(): c for c in df.columns}
    if "cui" not in cols_lower or "def" not in cols_lower:
        raise ValueError(f"MGDEF must contain columns 'CUI' and 'DEF' (found: {list(df.columns)})")

    col_cui = cols_lower["cui"]
    col_def = cols_lower["def"]

    # Strip and filter
    df[col_cui] = df[col_cui].astype(str).str.strip()
    df[col_def] = df[col_def].astype(str).str.strip()
    df = df[(df[col_cui] != "") & (df[col_def] != "")]

    # If duplicate CUI entries exist, keep the first occurrence
    df = df.drop_duplicates(subset=[col_cui], keep="first")

    cui_def_map = pd.Series(df[col_def].values, index=df[col_cui].values).to_dict()
    return cui_def_map


def read_medgen_mapping(mapping_path: str) -> pd.DataFrame:
    """
    Read MedGen_HPO_OMIM_Mapping.txt.gz using the header (which may start with '#').
    Detects delimiter from header and returns a DataFrame with header names normalized.
    """
    if not os.path.exists(mapping_path):
        raise FileNotFoundError(f"Mapping file not found: {mapping_path}")

    lines = read_text_gz_all(mapping_path)
    # find first non-empty line which should be header (can start with '#')
    idx = 0
    while idx < len(lines) and lines[idx].strip() == "":
        idx += 1
    if idx >= len(lines):
        raise ValueError("Mapping file appears empty")

    header_line = lines[idx].strip()
    if header_line.startswith("#"):
        header_clean = header_line.lstrip("#").strip()
    else:
        header_clean = header_line

    delim = detect_delim_from_header(header_clean)
    header_cols = [c.strip() for c in header_clean.split(delim)]

    # Data lines begin after header
    data_lines = lines[idx + 1 :]
    data_text = "\n".join(data_lines).strip()
    if data_text == "":
        # create empty DataFrame with header_cols
        df = pd.DataFrame(columns=header_cols)
        return df

    df = pd.read_csv(
        StringIO(data_text),
        sep=delim,
        header=None,
        names=header_cols,
        dtype=str,
        keep_default_na=False,
        engine="python",
    )

    # Strip whitespace in all cells
    df = df.fillna("").astype(str)
    df = df.applymap(lambda x: x.strip())
    return df


def build_list_of_dicts(mapping_df: pd.DataFrame, cui_def_map: Dict[str, str]) -> List[Dict]:
    """
    For each row in mapping_df build a dict:
      {
        "_id": MIM_number,
        "omim_id": MIM_number,
        "omim_disease": OMIM_name,
        "medgen_concept_id": HPO_CUI or OMIM_CUI,
        "medgen_disease_info": DEF (from MGDEF mapped by CUI)
      }
    Skips rows missing MIM_number or OMIM_CUI.
    """
    results = []

    # Make a case-insensitive col lookup
    colmap = {c.lower(): c for c in mapping_df.columns}

    def get(row, name):
        key = name.lower()
        return row[colmap[key]] if key in colmap else ""

    ss = set()
    for _, row in mapping_df.iterrows():
        # skip completely empty rows
        if all((str(x).strip() == "") for x in row):
            continue

        omim_cui = get(row, "OMIM_CUI")
        mim_number = get(row, "MIM_number")
        omim_name = get(row, "OMIM_name")
        hpo_cui = get(row, "HPO_CUI")

        # Skip rows missing MIM number or OMIM_CUI
        if not mim_number or not omim_cui:
            continue

        if omim_cui in ss:
            continue
        ss.add(omim_cui)
        # medgen_concept_id prefer HPO_CUI, else OMIM_CUI
        medgen_concept_id = omim_cui

        # Look up definition using OMIM_CUI (user said: map #OMIM_CUI to CUI)
        # The mapping key is the same string (e.g., "C0432273")
        medgen_disease_info = cui_def_map.get(omim_cui, "NA")
        # If not found by OMIM_CUI, try HPO_CUI as fallback
        if not medgen_disease_info and hpo_cui:
            medgen_disease_info = cui_def_map.get(hpo_cui, "NA")

        entry = {
            "_id": mim_number,
            "omim_id": mim_number,
            "omim_disease": omim_name,
            "medgen_concept_id": medgen_concept_id,
            "medgen_disease_info": medgen_disease_info,
        }
        results.append(entry)

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Build OMIM-MedGen JSON from mapping + MGDEF")
    parser.add_argument("--mapping", default="MedGen_HPO_OMIM_Mapping.txt.gz", help="Path to MedGen_HPO_OMIM_Mapping.txt.gz")
    parser.add_argument("--mgdef", default="MGDEF.csv.gz", help="Path to MGDEF.csv.gz")
    parser.add_argument("--output", default="omim_medgen_data_new2.json", help="Output JSON file")
    args = parser.parse_args()

    print(f"Loading MGDEF from: {args.mgdef}")
    cui_def_map = read_mgdef_map(args.mgdef)
    print(f"Loaded {len(cui_def_map)} definitions from MGDEF")

    print(f"Loading mapping from: {args.mapping}")
    mapping_df = read_medgen_mapping(args.mapping)
    print(f"Mapping rows read: {len(mapping_df)}")

    results = build_list_of_dicts(mapping_df, cui_def_map)
    print(f"Built {len(results)} output records")

    # Write JSON array
    with open(args.output, "w", encoding="utf-8") as outf:
        json.dump(results, outf, ensure_ascii=False, indent=2)
    print(f"Wrote output to {args.output}")


if __name__ == "__main__":
    main()

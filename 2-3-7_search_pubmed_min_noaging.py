import os
import time
import argparse
import pandas as pd
from Bio import Entrez
from dotenv import load_dotenv
from pathlib import Path
import re  # **NEW10/28** add regex support

# --- Setup ---
load_dotenv()  # loads .env if present
Entrez.email = os.getenv("EMAIL", "apriloctober17@hotmail.com")
Entrez.api_key = os.getenv("NCBI_API_KEY","48035e8fafb81cfbe03a235c56c14ecc3909")

# Simple rate limit to be polite to NCBI
SLEEP_SEC = 0.34  # ~3 calls/sec if you have API key; 0.5-1s if no key

def search_pubmed(query, retmax=300, sort="relevance"):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax, sort=sort)
    results = Entrez.read(handle)
    handle.close()
    pmids = results.get("IdList", [])
    total = int(results.get("Count", 0))
    return pmids, total

def _species_block(species: str) -> str:
    species = (species or "any").lower()
    if species == "mice":
        return "(Mice[MeSH Terms] NOT Humans[MeSH Terms])"
    if species == "rats":
        return "(Rats[MeSH Terms] NOT Humans[MeSH Terms])"
    return "(Animals[MeSH Terms] NOT Humans[MeSH Terms])"

def _classify_tier(article) -> str:
    """
    Classify record as 'clinical', 'animal', 'review', or ''.
    """
    try:
        med = article.get("MedlineCitation", {})
        art = med.get("Article", {})
        clinical_types = {
            "Clinical Trial",
            "Randomized Controlled Trial",
            "Controlled Clinical Trial",
            "Clinical Study",
            "Pragmatic Clinical Trial",
            "Multicenter Study",
        }
        review_types = {
            "Review", "Systematic Review", "Meta-Analysis",
            "Scoping Review", "Evidence Synthesis", "Umbrella Review",
        }

        ptypes = [str(pt) for pt in art.get("PublicationTypeList", [])]

        mesh = med.get("MeshHeadingList", [])
        mesh_terms = {str(mh.get("DescriptorName", "")) for mh in mesh}
        has_humans = "Humans" in mesh_terms
        has_animals = "Animals" in mesh_terms

        # First, classify reviews
        if any(pt.lower() in {rt.lower() for rt in review_types} or "review" in pt.lower() for pt in ptypes):
            return "review"

        if any(pt in clinical_types for pt in ptypes) and has_humans:
            return "clinical"
        if has_animals and not has_humans:
            return "animal"
        return ""
    except Exception:
        return ""

def fetch_details(pmids):
    if not pmids:
        return []
    joined = ",".join(pmids)
    handle = Entrez.efetch(db="pubmed", id=joined, rettype="xml", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    articles = []
    for article in records["PubmedArticle"]:
        med = article.get("MedlineCitation", {})
        art = med.get("Article", {})
        journal = art.get("Journal", {})
        jl = journal.get("JournalIssue", {})
        pub_date = jl.get("PubDate", {})
        year = pub_date.get("Year") or pub_date.get("MedlineDate") or ""

        title = art.get("ArticleTitle", "")
        abstract = ""
        if "Abstract" in art and "AbstractText" in art["Abstract"]:
            parts = art["Abstract"]["AbstractText"]
            abstract = " ".join(str(p) for p in parts) if isinstance(parts, list) else str(parts)

        authors = []
        for a in art.get("AuthorList", []):
            last = a.get("LastName") or ""
            fore = a.get("ForeName") or ""
            if last or fore:
                authors.append(f"{last} {fore}".strip())
        authors_str = "; ".join(authors)

        journal_title = journal.get("Title", "")
        pmid = med.get("PMID", "")
        doi = ""
        art_ids = article.get("PubmedData", {}).get("ArticleIdList", [])
        for aid in art_ids:
            if aid.attributes.get("IdType") == "doi":
                doi = str(aid)

        tier = _classify_tier(article)
        has_review_kw = 1 if re.search(r"\breview\b", f"{title} {abstract}", flags=re.IGNORECASE) else 0

        articles.append({
            "PMID": str(pmid),
            "Year": str(year),
            "Title": str(title),
            "Journal": str(journal_title),
            "Authors": authors_str,
            "DOI": doi,
            "Abstract": abstract,
            "Tier": tier,
            "HasReviewKW": has_review_kw
        })
        time.sleep(SLEEP_SEC)
    return articles

def main():
    parser = argparse.ArgumentParser(description="Search PubMed for GlyNAC → Genomic Instability")
    parser.add_argument("--retmax", type=int, default=300, help="Max number of results")
    parser.add_argument("--sort", type=str, default="relevance",
                        choices=["relevance", "pub+date", "most+recent"], help="Sort order")
    parser.add_argument("--out", type=str,
                        default=r"C:\Users\apartment\OneDrive\Echo\Aging Project\17_GlyNAC_outputdata_no_aging.xlsx",
                        help="Output Excel filename (.xlsx). Requires `pip install openpyxl`.")
    parser.add_argument("--query", type=str, default=None, help="Custom PubMed query")

    parser.add_argument("--filter_priority", choices=["clinical", "animal", "none"], default="none",
                        help="Prioritize clinical trials (humans) or animal studies (non-human).")
    parser.add_argument("--species", choices=["any", "mice", "rats"], default="any",
                        help="Species focus when --filter_priority animal; 'any' = all non-human animals.")

    parser.add_argument("--include_reviews", action="store_true",
                        help="**NEW10/30_2** Include reviews/meta-analyses (do NOT require methods keywords).")

    args = parser.parse_args()

    # **NEW10/30_2** Split your query into TOPIC and METHODS
    topic_clause = r"""
    (
    "GlyNAC"[tiab] '
    OR (glycine[tiab] AND ("N-acetylcysteine"[tiab] OR "N acetylcysteine"[tiab] OR acetylcystein*[tiab] OR NAC[tiab]))'
    )
    AND
    (
    "oral"[tiab] OR "oral administration"[tiab] OR capsule*[tiab] OR tablet*[tiab] OR softgel*[tiab]
    OR lozenge*[tiab] OR sublingual*[tiab] OR syrup[tiab] OR solution[tiab] OR supplement*[tiab]
    )
    AND
    (
    proteostasis[tiab] OR "protein homeostasis"[tiab] OR "protein quality control"[tiab] OR PQC[tiab]
    OR "unfolded protein response"[tiab] OR UPR[tiab]
    OR "endoplasmic reticulum stress"[tiab] OR "ER stress"[tiab]
    OR "heat shock response"[tiab] OR HSR[tiab]
    OR "molecular chaperone*"[tiab] OR chaperonin*[tiab]
    OR "ubiquitin-proteasome system"[tiab] OR proteasome[tiab]
    OR "ER-associated degradation"[tiab] OR ERAD[tiab]
    OR "protein aggregation"[tiab] OR "misfolded protein*"[tiab] OR amyloid[tiab]
    OR "proteotoxic stress"[tiab]
    )
    )
    
    


 

    """


    methods_clause = r"""
    AND
    (randomized[tiab] OR randomised[tiab] OR trial[tiab] OR placebo[tiab] OR "double blind"[tiab]
     OR cohort[tiab] OR "case-control"[tiab] OR "cross-sectional"[tiab]
     OR prospective[tiab] OR retrospective[tiab] OR "phase I"[tiab] OR "phase II"[tiab] OR "phase III"[tiab]
     OR "pilot study"[tiab] OR participants[tiab] OR patients[tiab]
     OR mouse[tiab] OR mice[tiab] OR murine[tiab] OR rat[tiab] OR zebrafish[tiab] OR drosophila[tiab] OR "C. elegans"[tiab]
     OR "in vivo"[tiab] OR "in vitro"[tiab]
     OR "western blot"[tiab] OR qPCR[tiab] OR "RT-qPCR"[tiab] OR "RNA-seq"[tiab] OR "single-cell"[tiab]
     OR "flow cytometry"[tiab] OR ELISA[tiab] OR immunohistochemistry[tiab]
     OR CRISPR[tiab] OR siRNA[tiab] OR knockout[tiab] OR knockdown[tiab] OR overexpress*[tiab]
     OR assay[tiab] OR microscopy[tiab]
     OR "we measured"[tiab] OR "we investigated"[tiab] OR "we collected"[tiab] OR "sample size"[tiab] OR "n="[tiab])
    """

    # **NEW10/30_2** Build the final query depending on --include_reviews
    if args.include_reviews:
        query = topic_clause  # reviews allowed; don't require methods
    else:
        query = topic_clause + methods_clause  # keep primary-study signals

    # **NEW10/30_2** Only append clinical filter when not including reviews
    clinical_filter = r"""
    AND ("Clinical Trial"[pt] OR "Randomized Controlled Trial"[pt] OR "Controlled Clinical Trial"[pt]
         OR "Observational Study"[pt] OR "Comparative Study"[pt] OR "Evaluation Study"[pt]
         OR "Validation Study"[pt] OR "Multicenter Study"[pt])
    """
    if args.filter_priority == "clinical" and not args.include_reviews:
        query += clinical_filter

    # Animal priority still works either way
    if args.filter_priority == "animal":
        query += " AND " + _species_block(args.species)

    print(f"Using query:\n{query}\n")

    # Search
    pmids, total = search_pubmed(query, retmax=args.retmax, sort=args.sort)
    print(f"PubMed total matches for this query: {total}")
    print(f"IDs returned (up to retmax): {len(pmids)}")

    # Fetch + DataFrame
    articles = fetch_details(pmids)
    tier_order = {"clinical": 0, "animal": 1, "review": 2, "": 3}
    df = pd.DataFrame(
        articles,
        columns=["PMID", "Year", "Title", "Journal", "Authors", "DOI", "Abstract", "Tier", "HasReviewKW"]
    )
    df["Year_num"] = pd.to_numeric(df["Year"].str[:4], errors="coerce")
    df["Tier_rank"] = df["Tier"].map(tier_order).fillna(3).astype(int)
    df = df.sort_values(["Tier_rank", "Year_num"], ascending=[True, False]).drop(columns=["Year_num", "Tier_rank"])

    # **NEW10/30_3** Compute counts by Tier (treat empty as 'blank')
    tier_counts = (
        df.assign(Tier=df["Tier"].replace({"": "blank"}))
          .groupby("Tier")
          .size()
          .reindex(["clinical", "animal", "review", "blank"], fill_value=0)
    )
    total_n = int(tier_counts.sum())

    # **NEW10/30_3** Print a breakdown to console
    print("\nTier breakdown:")
    for k, v in tier_counts.items():
        print(f"  {k}: {v}")
    print(f"  total: {total_n}\n")

    # **NEW10/30_3** Prepare a Summary DataFrame for Excel
    summary_df = tier_counts.reset_index()
    summary_df.columns = ["Tier", "Count"]
    summary_df.loc[len(summary_df.index)] = ["total", total_n]

    # Save
    out_path = Path(args.out)
    if out_path.suffix.lower() != ".xlsx":
        out_path = out_path.with_suffix(".xlsx")
        print(f"Output extension corrected to: {out_path}")
    # **NEW10/30_3** Write two sheets: Results + Summary
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Results")
        summary_df.to_excel(writer, index=False, sheet_name="Summary")
    print(f"Saved {len(df)} records to {out_path}")
    print("Also wrote a 'Summary' sheet with counts for clinical/animal/review/blank.")

def new_func():
    return '('  # **NEW10/28-XLSX**

if __name__ == "__main__":
    main()

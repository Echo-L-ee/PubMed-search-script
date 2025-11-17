import os
import time
import argparse
import pandas as pd
from Bio import Entrez
from dotenv import load_dotenv
from pathlib import Path
import re  # regex for text-based classification


# === Helpers for Species + PaperType + scoring ===

def _collect_text_safe_1111(row):
    """
    Safely collect multiple text fields into a single lowercase string
    for text-based pattern matching.
    """
    parts = []
    for k in ["Title", "Abstract", "PublicationTypeList", "MeSHHeadingList", "Journal"]:
        try:
            v = row.get(k, "")
        except AttributeError:
            v = ""
        if v is None:
            continue
        if isinstance(v, list):
            v = " ".join(map(str, v))
        parts.append(str(v))
    return " ".join(parts).lower()


def _detect_papertype_1111(row):
    """
    Decide PaperType = ClinicalTrial / Observational / Review / Other / Research.

    Priority we want:
      ClinicalTrial > Observational > Research > Review > Other
    """
    species = str(row.get("SpeciesContext", "") or "").strip().lower()
    pub_type_str = (row.get("PublicationTypeList") or "").strip()
    pub_types_lower = pub_type_str.lower()

    # --- NEW: MeSH-based human/animal flags ---
    mesh_raw = str(row.get("MeSHHeadingList", "") or "").lower()
    mesh_human_terms = [" humans ", " human "]
    mesh_animal_terms = [
        " animals ", " animals, laboratory ",
        " mice ", " mouse ",
        " rats ", " rat ",
        " zebrafish ", " drosophila ",
        " c. elegans ",
        " disease models, animal "
    ]
    mesh_has_humans = any(t in mesh_raw for t in mesh_human_terms)
    mesh_has_animals = any(t in mesh_raw for t in mesh_animal_terms)
    pure_human_mesh = mesh_has_humans and not mesh_has_animals

    # --- 0) PublicationTypeList-based with explicit priority ---
    if pub_type_str:
        is_vet_trial = any(v in pub_types_lower for v in VET_TRIAL_PUBTYPES_1111)
        is_clin_trial = any(pt in pub_types_lower for pt in CLINICAL_TRIAL_PUBTYPES_1111)
        is_obs       = any(pt in pub_types_lower for pt in OBSERVATIONAL_PUBTYPES_1111)
        is_review    = any(pt in pub_types_lower for pt in REVIEW_PUBTYPES_1111)
        is_other     = any(pt in pub_types_lower for pt in OTHER_PUBTYPES_1111)

        # 0a. Veterinary trials → treat as generic Research
        if is_vet_trial:
            return "Research"

        # 1. Human clinical trials
        if species == "human" and is_clin_trial:
            return "ClinicalTrial"

        # 2. Human observational designs – ONLY pure human MeSH (no animals)
        if species == "human" and is_obs and pure_human_mesh:
            return "Observational"

        # 3. Reviews (including meta-analyses)
        if is_review:
            return "Review"

        # 4. Other non-primary content
        if is_other:
            return "Other"
        # else: fall through to text-based rules

    # --- 1) Fallback: text-based classification ---
    text = _collect_text_safe_1111(row).lower()
    tx = f" {text} "

    review_terms_any = [
        " systematic review ", " meta-analysis ", " meta analysis ",
        " scoping review ", " umbrella review ", " rapid review ", " living review ",
        " narrative review ", " integrative review ", " realist review ",
        " overview of reviews ", " evidence review ", " evidence map ", " evidence mapping ",
        " pooled analysis ", " meta-regression ", " metaregression ",
        " network meta-analysis ", " network meta analysis ",
        " bibliometric analysis ", " scientometric analysis ", " citation analysis ",
        " mini-review ", " mini review ", " literature review ", " state of the art "
    ]

    other_proto_methods = [
        " study protocol ", " trial protocol ", " protocol ", " statistical analysis plan ", " sap ",
        " pre-registration ", " preregistration ", " registered report ", " design and rationale ",
        " guideline ", " practice guideline ", " reporting guideline ", " prisma ",
        " consensus ", " consensus statement ", " position statement ", " recommendation ", " checklist ",
        " standard operating procedure ", " sop ", " standardization ", " harmonization ",
        " methods paper ", " methodology ", " methodological ",
        " assay development ", " assay validation ", " benchmarking ", " benchmark ",
        " workflow ", " pipeline ", " toolkit ", " software ", " package ", " web server ", " webserver ",
        " editorial ", " commentary ", " perspective ", " viewpoint ", " opinion ",
        " letter to the editor ", " correspondence ", " white paper ", " concept paper ", " roadmap "
    ]
    human_outcome_cues = [
        " participants ", " patients ", " subjects ", " volunteers ",
        " enrolled ", " recruited ", " randomized ", " randomised ", " assigned ", " allocation ",
        " primary outcome ", " secondary outcome ", " endpoint ", " endpoints ", " outcomes ",
        " adverse event ", " adverse events ", " safety ", " tolerability ",
        " efficacy ", " effectiveness ", " follow-up ", " follow up ", " arm ", " arms "
    ]

    trial_design_terms = [
        " randomized ", " randomised ", " randomization ", " randomisation ", " allocation ",
        " double blind ", " single blind ", " triple blind ", " quadruple blind ",
        " placebo ", " sham ", " active comparator ", " standard of care ",
        " crossover ", " cross-over ", " parallel group ", " factorial ",
        " phase i ", " phase ii ", " phase iii ", " phase iv ",
        " interventional ", " controlled trial ", " control group ",
        " intention-to-treat ", " per-protocol ", " clinicaltrials.gov ", " nct", " isrctn ",
        " eudract ", " umin-ctr ", " chictr ", " ctri "
    ]
    participant_cues = [
        " participants ", " patients ", " subjects ", " volunteers ",
        " enrolled ", " recruited ", " assigned ", " allocation ",
        " baseline ", " endpoint ", " outcomes ", " follow-up ", " follow up ",
        " primary outcome ", " secondary outcome ", " adverse events ", " safety ", " efficacy "
    ]
    observational_terms = [
        " cohort ", " case-control ", " case control ", " cross-sectional ", " cross sectional ",
        " prospective ", " retrospective ", " registry ", " real-world ", " population-based ", " population based "
    ]

    # animal cues to block observational when mixed
    animal_text_cues = [
        " mouse ", " mice ", " murine ", " rat ", " rats ",
        " zebrafish ", " drosophila ", " c. elegans ",
        " xenograft ", " animal model ", " disease models, animal "
    ]

    # --- 1c. Human clinical trial (text) ---
    if (
        species == "human"
        and any(t in tx for t in trial_design_terms)
        and any(c in tx for c in participant_cues)
        and " clinical trial, veterinary " not in tx
    ):
        return "ClinicalTrial"

    # --- 1d. Human observational (text) – ONLY if pure human and no animal cues ---
    if (
        species == "human"
        and any(t in tx for t in observational_terms)
        and not any(t in tx for t in trial_design_terms)
        and pure_human_mesh
        and not any(a in tx for a in animal_text_cues)
    ):
        return "Observational"

    # --- 1e. Reviews by text ---
    is_review_text = (
        any(s in tx for s in review_terms_any)
        or re.search(r"\breview(?:s)?\b", text)
        or " heart failure reviews " in tx
        or (" review " in tx and " clinical trial, veterinary " not in tx)
    )
    if is_review_text:
        return "Review"

    # --- 1f. "Other" – protocols/guidelines/methods without clear human outcomes ---
    if any(s in tx for s in other_proto_methods) and not any(h in tx for h in human_outcome_cues):
        return "Other"

    # --- 1g. Veterinary trials (text) → Research ---
    if " clinical trial, veterinary " in tx:
        return "Research"

    # --- 1h. Case reports / case series: treat as generic Research ---
    if " case report " in tx or " case series " in tx:
        return "Research"

    # --- 1i. Default: generic Research ---
    return "Research"


# PublicationTypeList-based sets
CLINICAL_TRIAL_PUBTYPES_1111 = {
    "clinical trial",
    "clinical trial, phase i",
    "clinical trial, phase ii",
    "clinical trial, phase iii",
    "clinical trial, phase iv",
    "controlled clinical trial",
    "randomized controlled trial",
    "pragmatic clinical trial",
    "clinical study",
    "multicenter study",
}

REVIEW_PUBTYPES_1111 = {
    "review",
    "systematic review",
    "meta-analysis",
    "scoping review",
    "evidence synthesis",
    "umbrella review",
}

OTHER_PUBTYPES_1111 = {
    "editorial",
    "comment",
    "letter",
    "news",
    "biography",
    "interview",
    "congresses",
    "case reports",
}

VET_TRIAL_PUBTYPES_1111 = {
    "clinical trial, veterinary"
}

# explanation: NEW set for human observational designs
OBSERVATIONAL_PUBTYPES_1111 = {
    "observational study",
    "comparative study",
    "cohort studies",
    "case-control studies",
    "cross-sectional studies",
    "prospective studies",
    "retrospective studies",
    "longitudinal studies",
    "follow-up studies",
    "comparative study",
}


def _detect_species_1111(row):
    """
    SpeciesContext:
      - 'Human'         # explanation: human participants / human population data
      - 'Animal'        # explanation: in vivo non-human animals (mice, rats, etc.)
      - 'cell_human'    # explanation: in vitro work using clearly human-derived cells
      - 'cell_animal'   # explanation: in vitro work using clearly animal-derived cells
      - 'cell_offscope' # explanation: generic in vitro / cell culture, species unclear
      - 'Other'         # explanation: not clearly classifiable; mixed / background-only

    Overall logic:
      1) First, look at MeSHHeadingList to catch "pure" human or animal:
         - If MeSH has Humans and NO Animals and NO cell terms → 'Human'.
         - If MeSH has Animals and NO Humans and NO cell terms → 'Animal'.
         - Otherwise (mixed humans+animals, or in vitro tags, or missing MeSH) → go to text scoring.

      2) Text-based scoring (internal scores, NOT your final 10/5/4/3/2/1):
           - Human participants (patients, participants, etc.) → Human score up to 8
           - Animal in vivo (mouse/rat/strain/in vivo phrases) → Animal score 5
           - Clearly human cell lines → cell_human score 10
           - Clearly animal cell lines → cell_animal score 4
           - Generic "in vitro / cell culture" → cell_offscope score 3

      3) Final decision rules (using these internal scores):
           a) Let cell win ONLY if a cell_* score is ≥ both Human and Animal.
              (prevents one small "in vitro" from overriding a big mouse/human study)
           b) If both Human and Animal have evidence → choose 'Animal'
              (conservative: avoid mislabeling animal work as human)
           c) If only one of Human / Animal has evidence → choose that one.
           d) Otherwise, fall back to the highest non-zero label with a simple priority.

      4) If absolutely nothing fires, but we see generic "animal" words → 'Other';
         else → 'Other'.
    """

    # --- Collect combined text (title, abstract, pub types, MeSH, journal) ---
    text = _collect_text_safe_1111(row)                # explanation: helper already lowercases & concatenates fields
    tx = f" {text.lower()} "                           # explanation: pad spaces so ' mouse ' doesn't match 'house'

    # === 1) MeSH-based quick decision (pure human / pure animal) ===
    mesh_raw = str(row.get("MeSHHeadingList", "") or "")
    mesh_tx = f" {mesh_raw.lower()} "

    # explanation: only count explicit "Humans"/"Human" MeSH, not age terms like 'Adult'
    has_humans = " humans " in mesh_tx or " human " in mesh_tx

    animal_mesh_terms = [
        " animals ", " animals, laboratory ",
        " mice ", " mouse ",
        " rats ", " rat ",
        " rabbits ", " rabbit ",
        " swine ", " pigs ", " pig ",
        " dogs ", " dog ",
        " cats ", " cat ",
        " sheep ", " cattle ", " cows ",
        " non-human ", " nonhuman ",
        " primates ", " macaques ", " monkeys ",
        " disease models, animal "
    ]
    has_animals = any(t in mesh_tx for t in animal_mesh_terms)   # explanation: any non-human animal MeSH

    mesh_cell_terms = [
        " cells, cultured ",
        " cell line ", " cell line, tumor ",
        " cell lines ",
        " single-cell analysis ", " single cell analysis ",
        " in vitro "
    ]
    has_cell = any(t in mesh_tx for t in mesh_cell_terms)        # explanation: MeSH suggests in vitro / cell work

    # Pure human: Humans present, NO Animals, NO cell terms
    if has_humans and not has_animals and not has_cell:
        return "Human"                                           # explanation: strongly human-only MeSH → Human

    # Pure animal: Animals present, NO Humans, NO cell terms
    if has_animals and not has_humans and not has_cell:
        return "Animal"                                          # explanation: strongly animal-only MeSH → Animal

    # Anything else (mixed Humans+Animals, any cell terms, or missing MeSH) → go to detailed scoring


    # === 2) Text-based multi-hit scoring ===
    scores = {
        "Human": 0,
        "Animal": 0,
        "cell_human": 0,
        "cell_animal": 0,
        "cell_offscope": 0,
    }

    def bump(label, value):
        """Increase internal score for a label if 'value' is higher than current."""
        if scores[label] < value:
            scores[label] = value                             # explanation: keep the strongest evidence seen so far

    # 2a. Human participants / populations
    # explanation: strong human evidence = patients/participants/etc, not just the word "human"
    strong_human_cues = [
        " participants ", " patients ", " subjects ", " volunteers ",
        " men ", " women ",
        " outpatient ", " inpatient ", " clinic ", " hospital ",
        " cohort ", " registry ", " follow-up ", " follow up ",
        " nhanes ", " uk biobank ", " framingham ",
        " nurses' health study ", " health professionals follow-up study "
    ]

    # explanation: if we see strong cues → high score; if only 'human(s)' present → smaller human score
    has_strong_human = any(term in tx for term in strong_human_cues)
    has_generic_human = (" humans " in tx) or (" human " in tx)

    if has_strong_human:
        bump("Human", 8)                                       # explanation: real human participant / cohort signal
    elif has_generic_human:
        bump("Human", 5)                                       # explanation: softer human signal (e.g., "human gene")

    # IVF safeguard: IVF with participants -> Human
    if "in vitro fertilization" in tx or " ivf " in tx:
        if any(w in tx for w in [" participants ", " patients ", " women ", " men ", " couples ", " infants ", " neonates "]):
            bump("Human", 8)                                   # explanation: IVF almost always implies human clinical context

    # 2b. Animal in vivo (score 5)
    animal_terms = [
        " mouse ", " mice ", " murine ", " m. musculus ",
        " rat ", " rats ", " rodent ", " rodentia ",
        " c57bl/6 ", " c57bl6 ", " balb/c ", " balb c ",
        " sprague-dawley ", " sprague dawley ", " wistar ",
        " fischer 344 ", " f344 ", " lewis rat ", " long-evans ",
        " rabbit ", " rabbits ", " hamster ", " hamsters ", " gerbil ", " gerbils ",
        " guinea pig ", " guinea pigs ", " cavia porcellus ",
        " zebrafish ", " danio rerio ", " medaka ", " oryzias latipes ",
        " xenopus ", " xenopus laevis ", " xenopus tropicalis ",
        " drosophila ", " caenorhabditis elegans ", " c. elegans ",
        " porcine ", " swine ", " pig ", " pigs ", " piglet ", " piglets ",
        " ovine ", " sheep ", " lamb ", " lambs ",
        " caprine ", " goat ", " goats ",
        " bovine ", " cow ", " cows ", " cattle ", " heifer ", " heifers ",
        " equine ", " horse ", " horses ", " pony ", " ponies ",
        " canine ", " dog ", " dogs ", " beagle ", " beagles ",
        " feline ", " cat ", " cats ",
        " primate ", " primates ", " macaque ", " macaques ", " rhesus ",
        " cynomolgus ", " marmoset ", " baboon ",
        " chick ", " chicks ", " chicken ", " chickens ", " broiler ", " broilers ",
        " animal model ", " animal models ",
        " nonhuman primate ", " non-human primate ",
        " veterinary ", " preclinical ", " pre-clinical ",
        " xenograft ", " allograft ", " orthotopic ", " syngeneic ", " pdx ",
        " scid ", " nude mice ", " nude mouse ", " nsg ",
        " intraperitoneal ", " ip injection ", " subcutaneous ", " s.c. ",
        " intramuscular ", " i.m. ", " intravenous ", " i.v. ", " tail vein ",
        " oral gavage ", " gavage ", " in vivo "
    ]
    animal_strains = [" c57bl/6 ", " c57bl6 ", " balb/c ", " wistar ", " sprague-dawley ", " lewis rat "]

    if any(t in tx for t in animal_terms) or any(st in tx for st in animal_strains):
        bump("Animal", 5)                                      # explanation: explicit in vivo animal work

    unknown_animal_cues = [
        " animal ", " animals ", " mammal ", " mammals ",
        " rodent ", " rodents ", " vertebrate ", " vertebrates ",
        " in vivo study ", " in vivo experiment ", " in-vivo study ", " in-vivo experiment "
    ]
    unknown_animal_hit = any(term in tx for term in unknown_animal_cues)  # explanation: very generic animal wording


    # 2c. In vitro / cell-line mapping
    cellline_species = {
        " mda-mb-231 ": "human", " mda mb 231 ": "human",
        " mcf-7 ": "human", " mcf7 ": "human",
        " hela ": "human",
        " a549 ": "human",
        " hek293 ": "human", " hek-293 ": "human", " hek 293 ": "human",
        " huvec ": "human",
        " thp-1 ": "human",
        " jurkat ": "human",
        " hle b-3 ": "human", " hle-b3 ": "human",  
        " raw264.7 ": "animal", " raw 264.7 ": "animal",
        " 3t3 ": "animal", " nih-3t3 ": "animal",
        " c2c12 ": "animal",
        " cho ": "animal", " cho-k1 ": "animal",
        " llc ": "animal",
    }


    vitro_terms = [
        " cell line ", " cell lines ", " cell-line ", " cell based ", " cell-based ",
        " cell culture ", " cell cultures ", " cultured cells ", " cells, cultured ",
        " in vitro ", " ex vivo ",
        " organoid ", " organoids ", " spheroid ", " spheroids ",
        " patient-derived organoid ", " pdo ",
        " colony formation assay ", " clonogenic assay ",
        " wound healing assay ", " scratch assay ", " transwell assay ",
        " proliferation assay ", " mtt assay ", " cck-8 assay ", " cck8 assay ",
        " sirna", " shrna ", " knockdown ", " overexpression ",
        " crispr ", " genome editing ",
        " rt-qpcr", " rt qpcr", " qrt-pcr", " qrt pcr",
        " western blot", " immunoblot ",
        " single-cell ", " single cell "
    ]

    human_vitro_hit = any(k in tx for k, sp in cellline_species.items() if sp == "human") \
                      or any(t in tx for t in [" human cell ", " human cells ", " human-derived ", " human derived "])
    animal_vitro_hit = any(k in tx for k, sp in cellline_species.items() if sp == "animal") \
                       or any(t in tx for t in [" murine cell ", " mouse cell ", " rat cell ", " hamster cell "])
        # treat MeSH cell tags as generic in vitro signal even if text doesn't say "cell line"
    generic_vitro_hit = has_cell or any(t in tx for t in vitro_terms)


    # --- NEW: use MeSH species info to resolve generic in vitro hits ---
    if human_vitro_hit:
        bump("cell_human", 10)   # clearly human cell-line work
    if animal_vitro_hit:
        bump("cell_animal", 4)   # clearly animal cell-line work

    if generic_vitro_hit and not (human_vitro_hit or animal_vitro_hit):
        # If MeSH says Humans-only and not Animals → treat generic in vitro as human cell work
        if has_humans and not has_animals:
            bump("cell_human", 7)
        # If MeSH says Animals-only and not Humans → treat generic in vitro as animal cell work
        elif has_animals and not has_humans:
            bump("cell_animal", 3)
        # Otherwise we really don't know the species → keep as off-scope
        else:
            bump("cell_offscope", 3)


    # === 3) Choose best label with priority: cell (if dominant) > Animal > Human ===
    max_score = max(scores.values())
    if max_score > 0:
        # 3a. Let cell_* win only if a cell score is at least as large as Human and Animal.
        #     This avoids "cell_offscope" stealing papers that are mostly mouse/human.
        cell_labels = ["cell_human", "cell_animal", "cell_offscope"]          # explanation: all cell-related contexts
        cell_scores = {lbl: scores[lbl] for lbl in cell_labels}               # explanation: sub-dict for cell scores
        best_cell_label = max(cell_scores, key=lambda k: cell_scores[k])      # explanation: which cell_* is highest
        best_cell_score = cell_scores[best_cell_label]                        # explanation: value of that best cell score

        human_score = scores["Human"]
        animal_score = scores["Animal"]

        if best_cell_score > 0 and best_cell_score >= human_score and best_cell_score >= animal_score:
            # explanation: within cell_* labels, break ties as cell_human > cell_animal > cell_offscope
            cell_priority = {"cell_human": 3, "cell_animal": 2, "cell_offscope": 1}
            best_cell_label = max(
                cell_scores.items(),
                key=lambda kv: (kv[1], cell_priority[kv[0]])
            )[0]
            return best_cell_label                                       # explanation: clearly dominated by cell evidence

        # 3b. If both Human and Animal have evidence, compare strengths.     
        #     - If Human score is clearly higher (>= Animal + 2)             
        #       → treat as Human (e.g., human trial mentioning mice).        
        #     - Otherwise stay conservative and call it Animal.              
        if human_score > 0 and animal_score > 0:                             
            if human_score >= animal_score + 2:                              
                return "Human"                                               
            else:
                return "Animal"                                              
                                           # explanation: avoid calling mixed animal papers "Human"

        # 3c. If only one of Human / Animal has evidence, choose that one
        if animal_score > 0:
            return "Animal"
        if human_score > 0:
            return "Human"

        # 3d. Otherwise, fall back to any remaining non-zero label, with cell > Animal > Human
        priority = {
            "cell_human": 5,
            "cell_animal": 4,
            "cell_offscope": 3,
            "Animal": 2,
            "Human": 1,
        }
        best_label = max(scores.items(), key=lambda kv: (kv[1], priority.get(kv[0], 0)))[0]
        return best_label

    # === 4) Final fallback ===
    if unknown_animal_hit:
        return "Other"                                                   # explanation: generic animal wording but no clear category
    return "Other"                                                       # explanation: nothing matched → leave as 'Other'





def _species_score_1111(v):
    return {
        "Human": 10,
        "Animal": 5,
        "cell_human": 4,
        "cell_animal": 3,
        "cell_offscope": 2,
        "Other": 1,
    }.get(v, 1)


def _papertype_score_1111(v):
    return {
        "ClinicalTrial": 10,
        "Observational": 7,     
        "Research": 5,
        "Review": 3,
        "Other": 1
    }.get(v, 1)



def _apply_scoring_1111(df):
    """
    Add SpeciesContext, PaperType, SpeciesScore, PaperTypeScore, FinalScore
    to the DataFrame and move them to the end of the column order.
    """
    if "Title" not in df.columns:
        df["Title"] = ""
    if "Abstract" not in df.columns:
        df["Abstract"] = ""

    df["SpeciesContext"] = df.apply(_detect_species_1111, axis=1)
    df["PaperType"] = df.apply(_detect_papertype_1111, axis=1)
    df["SpeciesScore"] = df["SpeciesContext"].map(_species_score_1111)
    df["PaperTypeScore"] = df["PaperType"].map(_papertype_score_1111)
    df["FinalScore"] = df["SpeciesScore"] + df["PaperTypeScore"]

    keep = [c for c in df.columns if c not in ["SpeciesContext", "PaperType",
                                               "SpeciesScore", "PaperTypeScore", "FinalScore"]]
    return df[keep + ["SpeciesContext", "PaperType", "SpeciesScore", "PaperTypeScore", "FinalScore"]]


# === Entrez setup ===

load_dotenv()
Entrez.email = os.getenv("EMAIL", "apriloctober17@hotmail.com")
Entrez.api_key = os.getenv("NCBI_API_KEY", "48035e8fafb81cfbe03a235c56c14ecc3909")

SLEEP_SEC = 0.34  # polite rate limit


def search_pubmed(query, retmax=50, sort="relevance"):
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


def fetch_details(pmids):
    """
    Fetch PubMed details and include:
      - PublicationTypeList (joined string)
      - MeSHHeadingList (joined string)
    """
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

        # PublicationTypeList + MeSHHeadingList
        pubtypes_raw = art.get("PublicationTypeList", [])
        publication_types = [str(pt) for pt in pubtypes_raw]
        publication_type_str = "; ".join(publication_types)

        mesh_list = med.get("MeshHeadingList", [])
        mesh_terms = []
        for mh in mesh_list:
            dn = mh.get("DescriptorName")
            if dn is not None:
                mesh_terms.append(str(dn))
        mesh_str = "; ".join(mesh_terms)

        articles.append({
            "PMID": str(pmid),
            "Year": str(year),
            "Title": str(title),
            "Journal": str(journal_title),
            "Authors": authors_str,
            "DOI": doi,
            "Abstract": abstract,
            "PublicationTypeList": publication_type_str,
            "MeSHHeadingList": mesh_str,
        })
        time.sleep(SLEEP_SEC)
    return articles


def main():
    """
    Main entry point:
      1. Parse CLI arguments.
      2. Build the PubMed query (unless a custom one is provided).
      3. Fetch PubMed IDs and article details.
      4. Classify SpeciesContext and PaperType, add scores.
      5. Save results + summary to an Excel file.
    """
    parser = argparse.ArgumentParser(description="Search PubMed for 5_zinc")
    parser.add_argument("--retmax", type=int, default=50, help="Max number of results")
    parser.add_argument("--sort", type=str, default="relevance",
                        choices=["relevance", "pub+date", "most+recent"], help="Sort order")
    parser.add_argument("--out", type=str,
                        default=r"C:\Users\apartment\OneDrive\Echo\Aging Project\NEW\5_zinc.xlsx",
                        help="Output Excel filename (.xlsx). Requires `pip install openpyxl`.")
    parser.add_argument("--query", type=str, default=None, help="Custom PubMed query")

    parser.add_argument("--filter_priority", choices=["clinical", "animal", "none"], default="none",
                        help="Prioritize clinical trials (humans) or animal studies (non-human).")
    parser.add_argument("--species", choices=["any", "mice", "rats"], default="any",
                        help="Species focus when --filter_priority animal; 'any' = all non-human animals.")

    parser.add_argument("--include_reviews", action="store_true",
                        help="Include reviews/meta-analyses without requiring methods keywords.")

    args = parser.parse_args()

    # === 1) Build base query ===
    # If user provides a custom query, use it directly.           
    if args.query:                                                 
        query = args.query                                        
    else:                                                          
        topic_clause = r"""
        (
        "zinc supplement"[tiab]
        OR "zinc supplements"[tiab]
        OR "zinc supplementation"[tiab]
        OR "zinc-supplemented"[tiab]
        OR "zinc repletion"[tiab]
        OR "zinc-repleted"[tiab]
        OR "oral zinc"[tiab]
        OR ("Zinc"[MeSH Terms] AND "Dietary Supplements"[MeSH Terms])
        )
        AND
        (
        "oral"[tiab]
        OR "oral administration"[tiab]
        OR "by mouth"[tiab]
        OR "p.o."[tiab]
        OR "per os"[tiab]

        OR capsule*[tiab]
        OR tablet*[tiab]
        OR softgel*[tiab]
        OR pill*[tiab]
        OR lozenge*[tiab]
        OR troche*[tiab]
        OR "chewable"[tiab]
        OR "chewable tablet*"[tiab]
        OR gummy*[tiab]

        OR "oral powder"[tiab]
        OR powder*[tiab]
        OR granule*[tiab]
        OR sachet*[tiab]
        OR "oral solution"[tiab]
        OR "oral suspension"[tiab]
        OR syrup*[tiab]

        OR supplement*[tiab]
        OR "dietary supplement*"[tiab]()*
        )
        AND
        (
        senescence[tiab]
        OR "Cellular Senescence"[MeSH Terms]
        OR "Aging, Cellular"[MeSH Terms]
        OR "Aged"[MeSH Terms]
        OR "Aged, 80 and over"[MeSH Terms]
        OR gerontology[tiab]
        OR geroscience[tiab]
        OR progeria[tiab]
        OR "Progeria"[MeSH Terms]
        OR "age-related"[tiab]
        OR ageing[tiab]
        OR "age-associated"[tiab]
        OR "age-dependent"[tiab]
        )
        """
        query = topic_clause  # **NEW**: start from topic_clause

    # === 2) human-clinical filter ===
    clinical_filter = r"""
    AND (
    "Clinical Trial"[pt]
    OR "Clinical Trial, Phase I"[pt]
    OR "Clinical Trial, Phase II"[pt]
    OR "Clinical Trial, Phase III"[pt]
    OR "Clinical Trial, Phase IV"[pt]
    OR "Pragmatic Clinical Trial"[pt]
    OR "Randomized Controlled Trial"[pt]
    OR "Controlled Clinical Trial"[pt]
    OR "Clinical Study"[pt]
    OR "Pilot Projects"[pt]
    OR "Feasibility Studies"[pt]
    OR "Multicenter Study"[pt]
    OR "Observational Study"[pt]
    OR "Comparative Study"[pt]
    OR "Evaluation Study"[pt]
    OR "Validation Study"[pt]
    OR "Case-Control Studies"[pt]
    OR "Cohort Studies"[pt]
    OR "Cross-Sectional Studies"[pt]
    OR "Prospective Studies"[pt]
    OR "Retrospective Studies"[pt]
    OR "Cross-Over Studies"[pt]
    )
    AND Humans[MeSH Terms]
    AND NOT "Clinical Trial, Veterinary"[pt]
    AND NOT (Animals[MeSH Terms] NOT Humans[MeSH Terms])
    """

    # Add filters on top of the base query
    if args.filter_priority == "clinical":                         # **REVISED**: simpler condition
        query += clinical_filter

    if args.filter_priority == "animal":
        query += " AND " + _species_block(args.species)

    print(f"Using query:\n{query}\n")


    pmids, total = search_pubmed(query, retmax=args.retmax, sort=args.sort)
    print(f"PubMed total matches for this query: {total}")
    print(f"IDs returned (up to retmax): {len(pmids)}")

    articles = fetch_details(pmids)
    df = pd.DataFrame(
        articles,
        columns=[
            "PMID", "Year", "Title", "Journal", "Authors", "DOI",
            "Abstract", "PublicationTypeList", "MeSHHeadingList"
        ]
    )

    # Sort by year (newest first)
    df["Year_num"] = pd.to_numeric(df["Year"].str[:4], errors="coerce")
    df = df.sort_values(["Year_num"], ascending=[False]).drop(columns=["Year_num"])

    # Apply Species / PaperType / scoring
    df = _apply_scoring_1111(df)

    # === NEW PRIORITY LOGIC: clinical + methods inside Python ===
    def _is_human_clinical(row):
        """
        Treat as 'clinical-ish' if:
          - SpeciesContext is Human, AND
          - PaperType is ClinicalTrial or Observational.
        """
        species = str(row.get("SpeciesContext", "")).strip()
        ptype = str(row.get("PaperType", "")).strip()
        return int(species == "Human" and ptype in ("ClinicalTrial", "Observational"))

    def _has_methods_terms(row):
        """
        True if the combined text contains method / design / assay words.
        (We reuse a shorter subset of your old methods_clause list.)
        """
        text = _collect_text_safe_1111(row)  # already lowercased
        return int(any(
            kw in text for kw in [
                " randomized ", " randomised ", " controlled trial ", " placebo ",
                " double blind ", " single blind ", " open-label ",
                " cohort ", " case-control ", " case control ", " cross-sectional ", " cross sectional ",
                " prospective ", " retrospective ", " longitudinal ",
                " pilot study ", " feasibility study ",
                " western blot", " qpcr", " rt-qpcr", " pcr",
                " rna-seq", " rna sequencing ",
                " single-cell ", " single cell ",
                " we measured ", " we investigated ", " we collected ",
                " we performed ", " we conducted ", " we analyzed ",
            ]
        ))

    df["ClinicalPriority"] = df.apply(_is_human_clinical, axis=1)   # 1 = human clinical/observational, 0 = other
    df["MethodsPriority"] = df.apply(_has_methods_terms, axis=1)     # 1 = methods-y, 0 = no methods words

    # Combine into one ranking score:
    #   - Strong weight for human clinical
    #   - Smaller bonus for methods words
    #   - Keep your existing FinalScore as part of it
    df["PriorityScore"] = (
        3 * df["ClinicalPriority"]
        + 1 * df["MethodsPriority"]
        + df["FinalScore"]
    )

    # Sort by PriorityScore (highest first) and then by Year (newest first if Year is clean)
    df["Year_num"] = pd.to_numeric(df["Year"].str[:4], errors="coerce")
    df = df.sort_values(["PriorityScore", "Year_num"], ascending=[False, False])
    df = df.drop(columns=["Year_num"])


    # PaperType breakdown (no "blank" bucket)
    papertype_counts = (
        df["PaperType"]
        .value_counts()
        .reindex(["ClinicalTrial", "Observational", "Research", "Review", "Other"], fill_value=0)  
    )


    total_n = int(papertype_counts.sum())

    print("\nPaperType breakdown:")
    for k, v in papertype_counts.items():
        print(f"  {k}: {v}")
    print(f"  total: {total_n}\n")

    summary_df = papertype_counts.reset_index()
    summary_df.columns = ["PaperType", "Count"]
    summary_df.loc[len(summary_df.index)] = ["total", total_n]


    out_path = Path(args.out)
    if out_path.suffix.lower() != ".xlsx":
        out_path = out_path.with_suffix(".xlsx")
        print(f"Output extension corrected to: {out_path}")

    # Reorder columns for Excel:
    # - drop PublicationTypeList
    # - drop ClinicalPriority / MethodsPriority / PriorityScore (internal only)
    # - keep MeSHHeadingList
    # - move MeSHHeadingList to between Abstract and SpeciesContext
    cols = [
        c for c in df.columns
        if c not in ["PublicationTypeList", "ClinicalPriority", "MethodsPriority", "PriorityScore"]
    ]


    if "MeSHHeadingList" in cols and "Abstract" in cols:                  # explanation: only reorder if both exist
        cols.remove("MeSHHeadingList")                                    # explanation: take MeSHHeadingList out temporarily
        insert_pos = cols.index("Abstract") + 1                           # explanation: position right after Abstract
        cols.insert(insert_pos, "MeSHHeadingList")                        # explanation: put MeSHHeadingList back in that spot

    df_excel = df[cols]                                                   # explanation: use this ordered column list for Excel



    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        df_excel.to_excel(writer, index=False, sheet_name="Results")
        summary_df.to_excel(writer, index=False, sheet_name="Summary")

    print(f"Saved {len(df_excel)} records to {out_path}")
    print("Also wrote a 'Summary' sheet with counts for PaperType.")


def new_func():
    return '('


if __name__ == "__main__":
    main()














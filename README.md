# PubMed-search-script
Implement Python script to query PubMed

## GlyNAC ↔ Genomic Instability (PubMed Search)
Python script that queries PubMed for **GlyNAC** and **genomic instability** papers and saves a CSV (PMID, title, authors, journal, DOI, abstract, year).

## Install
pip install biopython pandas python-dotenv

## Setup
Create `.env`:
EMAIL=XXX@hotmail.com
NCBI_API_KEY=XXXXXX

## Run
python pubmed_search.py

## Output
CSV with: PMID, Year, Title, Journal, Authors, DOI, Abstract.retmax=200

## Notes
- Default query targets GlyNAC + genomic/genome/chromosomal instability and DNA damage terms.
- Add an API key to avoid NCBI rate limits.

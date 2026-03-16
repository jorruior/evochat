# Ensembl REST API client for comparative genomics queries
# Includes HGNC integration for gene symbol history

import time
import requests
from typing import Optional


BASE_URL = "https://rest.ensembl.org"

# Rate limiting: Ensembl allows 15 requests/second
_last_request_time = 0.0
_MIN_INTERVAL = 0.07


def _rate_limit():
    global _last_request_time
    now = time.time()
    elapsed = now - _last_request_time
    if elapsed < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - elapsed)
    _last_request_time = time.time()


class EnsemblClient:
    # Client for the Ensembl REST API focused on comparative genomics
    #
    # Usage:
    #   client = EnsemblClient()
    #   info = client.lookup_gene("BRCA1", species="human")
    #   orthologs = client.get_orthologs("ENSG00000012048")

    def __init__(self, server: str = BASE_URL):
        self.server = server.rstrip("/")
        self.session = requests.Session()
        self.session.headers.update({"Content-Type": "application/json"})

    def _get(self, endpoint: str, params: Optional[dict] = None) -> dict:
        ## GET request with rate limiting and retry on 429
        _rate_limit()
        url = f"{self.server}{endpoint}"
        resp = self.session.get(url, params=params)
        if resp.status_code == 429:
            retry_after = float(resp.headers.get("Retry-After", 1))
            time.sleep(retry_after)
            return self._get(endpoint, params)
        resp.raise_for_status()
        return resp.json()

    # ── Gene Lookup ──────────────────────────────────────────────────────

    def lookup_gene(self, gene: str, species: str = "human", expand: bool = False) -> dict:
        # Look up a gene by symbol or Ensembl ID
        # Returns dict with id, display_name, species, biotype, description, location, etc.
        if gene.startswith("ENS"):
            endpoint = f"/lookup/id/{gene}"
            params = {"expand": 1} if expand else {}
        else:
            endpoint = f"/lookup/symbol/{species}/{gene}"
            params = {"expand": 1} if expand else {}
        return self._get(endpoint, params)

    # ── Homology ─────────────────────────────────────────────────────────

    def get_homologs(
        self,
        gene_id: str,
        type: Optional[str] = None,
        target_species: Optional[str] = None,
        target_taxon: Optional[int] = None,
        format: str = "full",
        sequence: str = "none",
        source_species: str = "human",
    ) -> dict:
        # Retrieve homologs (orthologs/paralogs) for a gene
        # type: 'orthologues', 'paralogues', or None for all
        # format: 'full' includes species/id info, 'condensed' is minimal
        # source_species goes in the URL path per Ensembl API spec
        endpoint = f"/homology/id/{source_species}/{gene_id}"
        params = {"format": format, "sequence": sequence}
        if type:
            params["type"] = type
        if target_taxon:
            params["target_taxon"] = target_taxon

        # Handle multiple target species (comma, semicolon, or list)
        if target_species:
            species_list = []
            if isinstance(target_species, list):
                species_list = target_species
            elif "," in target_species or ";" in target_species:
                species_list = [s.strip() for s in target_species.replace(",", ";").split(";") if s.strip()]
            else:
                species_list = [target_species]

            ## Multiple species: make separate calls and merge results
            if len(species_list) > 1:
                all_homologies = []
                for sp in species_list:
                    sp_params = {**params, "target_species": sp}
                    data = self._get(endpoint, sp_params)
                    for entry in data.get("data", [{}]):
                        all_homologies.extend(entry.get("homologies", []))
                return {"data": [{"homologies": all_homologies}]}
            else:
                params["target_species"] = species_list[0]

        return self._get(endpoint, params)

    def get_orthologs(
        self,
        gene_id: str,
        target_species: Optional[str] = None,
        target_taxon: Optional[int] = None,
        subtype: Optional[str] = None,
        sequence: str = "none",
    ) -> list:
        # Get orthologs as a flat list
        # subtype: 'one2one', 'one2many', 'many2many', or None for all
        orth_type = f"orthologues_{subtype}" if subtype else "orthologues"
        data = self.get_homologs(
            gene_id,
            type=orth_type,
            target_species=target_species,
            target_taxon=target_taxon,
            sequence=sequence,
        )
        homologies = data.get("data", [{}])
        if homologies:
            return homologies[0].get("homologies", [])
        return []

    def get_paralogs(self, gene_id: str, sequence: str = "none") -> list:
        # Get paralogs as a flat list
        data = self.get_homologs(gene_id, type="paralogues", sequence=sequence)
        homologies = data.get("data", [{}])
        if homologies:
            return homologies[0].get("homologies", [])
        return []

    def count_homologs(self, gene_id: str) -> dict:
        # Count orthologs and paralogs by subtype
        # Omits 'type' param to get all homolog types in one call
        data = self.get_homologs(gene_id)
        homologies = data.get("data", [{}])
        if not homologies:
            return {}

        counts = {
            "orthologs_one2one": 0,
            "orthologs_one2many": 0,
            "orthologs_many2many": 0,
            "paralogs_within_species": 0,
            "paralogs_other": 0,
            "other": 0,
        }
        for h in homologies[0].get("homologies", []):
            htype = h.get("type", "")
            if htype == "ortholog_one2one":
                counts["orthologs_one2one"] += 1
            elif htype == "ortholog_one2many":
                counts["orthologs_one2many"] += 1
            elif htype == "ortholog_many2many":
                counts["orthologs_many2many"] += 1
            elif htype == "within_species_paralog":
                counts["paralogs_within_species"] += 1
            elif "paralog" in htype:
                counts["paralogs_other"] += 1
            else:
                counts["other"] += 1

        counts["total_orthologs"] = (
            counts["orthologs_one2one"]
            + counts["orthologs_one2many"]
            + counts["orthologs_many2many"]
        )
        counts["total_paralogs"] = (
            counts["paralogs_within_species"] + counts["paralogs_other"]
        )
        counts["total"] = counts["total_orthologs"] + counts["total_paralogs"] + counts["other"]
        return counts

    # ── Gene Trees ───────────────────────────────────────────────────────

    def get_gene_tree(
        self,
        gene_id: str,
        nh_format: str = "simple",
        aligned: bool = False,
        sequence: str = "none",
    ) -> dict:
        # Retrieve gene tree JSON for a gene
        # Endpoint requires species in URL path
        endpoint = f"/genetree/member/id/homo_sapiens/{gene_id}"
        params = {
            "nh_format": nh_format,
            "aligned": int(aligned),
            "sequence": sequence,
        }
        return self._get(endpoint, params)

    def get_gene_tree_newick(self, gene_id: str, nh_format: str = "simple") -> str:
        # Get Newick string for a gene tree
        # Uses text/x-nh content type to get raw Newick
        _rate_limit()
        url = f"{self.server}/genetree/member/id/homo_sapiens/{gene_id}"
        resp = self.session.get(
            url,
            params={"nh_format": nh_format},
            headers={"Content-Type": "text/x-nh"},
        )
        resp.raise_for_status()
        return resp.text

    # ── Sequences ────────────────────────────────────────────────────────

    def get_sequence(
        self,
        seq_id: str,
        type: str = "genomic",
        species: Optional[str] = None,
    ) -> dict:
        # Retrieve sequence for an Ensembl ID
        # type: 'genomic', 'cdna', 'cds', 'protein'
        # For cdna/cds/protein on gene IDs, resolve to transcript first (see chat.py)
        endpoint = f"/sequence/id/{seq_id}"
        params = {"type": type}
        if species:
            params["species"] = species
        return self._get(endpoint, params)

    def get_ortholog_sequences(
        self,
        gene_id: str,
        seq_type: str = "protein",
        target_species: Optional[str] = None,
        target_taxon: Optional[int] = None,
    ) -> list:
        # Get sequences of orthologs for a gene in FASTA-ready format
        orthologs = self.get_orthologs(
            gene_id,
            target_species=target_species,
            target_taxon=target_taxon,
            sequence=seq_type,
        )
        sequences = []
        for h in orthologs:
            target = h.get("target", {})
            entry = {
                "species": target.get("species", ""),
                "gene_id": target.get("id", ""),
                "protein_id": target.get("protein_id", ""),
                "type": h.get("type", ""),
            }
            seq = target.get("align_seq") or target.get("seq")
            if seq:
                entry["sequence"] = seq.replace("-", "")
            sequences.append(entry)
        return sequences

    # ── Cross-references ─────────────────────────────────────────────────

    def get_xrefs(
        self,
        gene_id: str,
        external_db: Optional[str] = None,
        all_levels: bool = False,
    ) -> list:
        # Get cross-references (UniProt, RefSeq, HGNC, EntrezGene, etc.)
        endpoint = f"/xrefs/id/{gene_id}"
        params = {"all_levels": int(all_levels)}
        if external_db:
            params["external_db"] = external_db
        return self._get(endpoint, params)

    def get_uniprot_ids(self, gene_id: str) -> list:
        ## Get UniProt Swiss-Prot IDs
        xrefs = self.get_xrefs(gene_id, external_db="Uniprot/SWISSPROT", all_levels=True)
        return list({x["primary_id"] for x in xrefs})

    def get_refseq_ids(self, gene_id: str) -> dict:
        ## Get RefSeq mRNA and peptide IDs
        mrna = self.get_xrefs(gene_id, external_db="RefSeq_mRNA", all_levels=True)
        pep = self.get_xrefs(gene_id, external_db="RefSeq_peptide", all_levels=True)
        return {
            "mRNA": list({x["primary_id"] for x in mrna}),
            "peptide": list({x["primary_id"] for x in pep}),
        }

    # ── Species Info ─────────────────────────────────────────────────────

    def get_species_info(self, species: str = "human") -> dict:
        return self._get(f"/info/genomes/{species}")

    def list_species(self) -> list:
        # List all species available in Ensembl
        data = self._get("/info/species")
        return data.get("species", [])

    def search_species(self, query: str) -> Optional[str]:
        # Search for a species by common name, scientific name, or alias
        # Returns Ensembl species name (e.g. 'amazona_vittata') or None
        query_lower = query.lower().strip().replace("-", " ")
        species_list = self.list_species()

        ## Exact match on name, display_name, common_name, aliases
        for sp in species_list:
            name = sp.get("name", "").lower()
            display = sp.get("display_name", "").lower()
            common = sp.get("common_name", "").lower()
            aliases = [a.lower() for a in sp.get("aliases", [])]

            if query_lower in (name, display, common):
                return sp["name"]
            if query_lower in aliases:
                return sp["name"]

        ## Substring match
        for sp in species_list:
            name = sp.get("name", "").lower()
            display = sp.get("display_name", "").lower()
            common = sp.get("common_name", "").lower()

            if query_lower in common or query_lower in display or query_lower in name:
                return sp["name"]
            query_words = set(query_lower.split())
            common_words = set(common.split())
            if query_words and query_words.issubset(common_words):
                return sp["name"]

        return None

    # ── Archive / ID History ─────────────────────────────────────────────

    def get_archive(self, gene_id: str) -> dict:
        # Get archive/history info for an Ensembl stable ID
        return self._get(f"/archive/id/{gene_id}")

    # ── HGNC Gene Name History ───────────────────────────────────────────

    def get_hgnc_info(self, symbol: str) -> dict:
        # Query HGNC REST API for gene symbol info (human only)
        # Returns prev_symbol, alias_symbol, name, hgnc_id, status, etc.
        _rate_limit()
        url = f"https://rest.genenames.org/fetch/symbol/{symbol}"
        resp = self.session.get(url, headers={"Accept": "application/json"})
        if resp.status_code != 200:
            ## Try by previous symbol if direct fetch fails
            resp = self.session.get(
                f"https://rest.genenames.org/search/prev_symbol/{symbol}",
                headers={"Accept": "application/json"},
            )
        resp.raise_for_status()
        data = resp.json()
        docs = data.get("response", {}).get("docs", [])
        if not docs:
            return {}
        return docs[0]

    def get_gene_history(self, gene: str, species: str = "human") -> dict:
        # Get comprehensive gene name history combining HGNC + Ensembl archive
        # Handles old/retired symbols by searching HGNC prev_symbol and alias_symbol
        result = {
            "current_symbol": None,
            "ensembl_id": None,
            "name": None,
            "previous_symbols": [],
            "aliases": [],
            "hgnc_id": None,
            "status": None,
            "ensembl_archive": None,
        }

        # Try resolving in Ensembl first
        ensembl_id = None
        symbol = gene

        if gene.startswith("ENS"):
            ensembl_id = gene
            try:
                info = self.lookup_gene(gene)
                symbol = info.get("display_name", gene)
            except Exception:
                pass
        else:
            try:
                info = self.lookup_gene(gene, species=species)
                ensembl_id = info.get("id")
                symbol = info.get("display_name", gene)
            except Exception:
                ## Symbol not found in Ensembl — will try HGNC below
                pass

        # Get HGNC data (human only) — also resolves old/alias symbols
        hgnc = {}
        if species in ("human", "homo_sapiens"):
            hgnc = self.get_hgnc_info(symbol)
            if not hgnc:
                ## Try searching by previous symbol
                _rate_limit()
                try:
                    resp = self.session.get(
                        f"https://rest.genenames.org/search/prev_symbol/{symbol}",
                        headers={"Accept": "application/json"},
                    )
                    if resp.ok:
                        docs = resp.json().get("response", {}).get("docs", [])
                        if docs:
                            hgnc = docs[0]
                except Exception:
                    pass
            if not hgnc:
                ## Try by alias
                _rate_limit()
                try:
                    resp = self.session.get(
                        f"https://rest.genenames.org/search/alias_symbol/{symbol}",
                        headers={"Accept": "application/json"},
                    )
                    if resp.ok:
                        docs = resp.json().get("response", {}).get("docs", [])
                        if docs:
                            hgnc = docs[0]
                except Exception:
                    pass

        # If HGNC found a current symbol and we still lack an Ensembl ID, resolve it
        if hgnc and not ensembl_id:
            current_sym = hgnc.get("symbol", symbol)
            try:
                info = self.lookup_gene(current_sym, species=species)
                ensembl_id = info.get("id")
                result["current_symbol"] = info.get("display_name", current_sym)
            except Exception:
                result["current_symbol"] = current_sym
        else:
            result["current_symbol"] = symbol

        result["ensembl_id"] = ensembl_id

        # Apply HGNC data
        if hgnc:
            result["name"] = hgnc.get("name")
            result["previous_symbols"] = hgnc.get("prev_symbol", [])
            result["aliases"] = hgnc.get("alias_symbol", [])
            result["hgnc_id"] = hgnc.get("hgnc_id")
            result["status"] = hgnc.get("status")
            if hgnc.get("symbol") and hgnc["symbol"] != gene:
                result["current_symbol"] = hgnc["symbol"]

        # Get Ensembl archive info
        if ensembl_id:
            try:
                archive = self.get_archive(ensembl_id)
                result["ensembl_archive"] = {
                    "id": archive.get("id"),
                    "version": archive.get("version"),
                    "release": archive.get("release"),
                    "assembly": archive.get("assembly"),
                    "is_current": archive.get("is_current"),
                }
            except Exception:
                pass

        return result

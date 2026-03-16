# Query parsers for evochat
# LLMParser: uses Ollama for natural language understanding
# KeywordParser: rule-based fallback when no LLM is available

import json
import re
from typing import Optional

_ollama = None


def _get_ollama():
    global _ollama
    if _ollama is None:
        try:
            import ollama
            _ollama = ollama
        except ImportError:
            raise ImportError(
                "ollama package not installed. Install with: pip install ollama\n"
                "You also need the Ollama server running: https://ollama.ai"
            )
    return _ollama


# ── System prompt for the LLM ────────────────────────────────────────────

SYSTEM_PROMPT = """You are a query parser for a gene evolution chatbot called evochat.
Your job is to extract structured intent from natural language queries about genes, 
orthologs, paralogs, sequences, gene trees, and cross-references.

You must return ONLY a valid JSON object (no markdown, no explanation) with these fields:

{
    "action": "<one of the actions below>",
    "gene": "<gene symbol or Ensembl ID, or null>",
    "species": "<source species, default 'human'>",
    "target_species": "<target species filter as Ensembl name (e.g. 'mus_musculus', 'pan_troglodytes'). For multiple species use comma-separated string: 'mus_musculus,pan_troglodytes'. Use null for all species>",
    "target_taxon": "<taxon name like 'mammals' or taxon ID, or null>",
    "seq_type": "<'protein' or 'cdna', default 'protein'>",
    "subtype": "<'one2one', 'one2many', 'many2many', or null>",
    "external_db": "<'UniProt', 'RefSeq', or null>"
}

Available actions:
- "lookup": Look up basic gene information
- "count_homologs": Count orthologs and paralogs
- "list_orthologs": List ortholog gene IDs and species
- "list_paralogs": List paralog gene IDs
- "get_ortholog_sequences": Get protein or CDS sequences of orthologs
- "get_gene_tree": Retrieve gene tree (Newick)
- "get_xrefs": Get cross-references (UniProt, RefSeq, etc.)
- "get_sequence": Get the sequence of a specific gene/transcript
- "list_species": List all species available in Ensembl
- "gene_history": Get previous gene names/symbols and aliases (e.g. "previous names of TUG1", "gene history of BRCA1", "aliases of TP53", "what was X called before")
- "help": User is asking for help or what the bot can do
- "unknown": Cannot determine intent

Taxon mapping (use NCBI taxon IDs):
- "mammals" or "mammalia" → 40674
- "primates" → 9443
- "rodents" or "rodentia" → 9989
- "vertebrates" or "vertebrata" → 7742
- "fish" or "actinopterygii" → 7898
- "birds" or "aves" → 8782
- "insects" or "insecta" → 50557
- "plants" or "viridiplantae" → 33090

Examples:
- "How many orthologs does BRCA1 have?" → {"action": "count_homologs", "gene": "BRCA1", "species": "human", ...}
- "Get mouse ortholog of TP53" → {"action": "list_orthologs", "gene": "TP53", "species": "human", "target_species": "mus_musculus", ...}
- "Show me the protein sequences of LZTR1 orthologs in mammals" → {"action": "get_ortholog_sequences", "gene": "LZTR1", "species": "human", "target_taxon": 40674, "seq_type": "protein", ...}
- "What's the gene tree for ENSG00000012048?" → {"action": "get_gene_tree", "gene": "ENSG00000012048", ...}
- "UniProt IDs for BRCA2" → {"action": "get_xrefs", "gene": "BRCA2", "external_db": "UniProt", ...}
- "TTN orthologs in chimp and mouse" → {"action": "list_orthologs", "gene": "TTN", "species": "human", "target_species": "pan_troglodytes,mus_musculus", ...}
- "Is TUG1 annotated in parrot?" → {"action": "list_orthologs", "gene": "TUG1", "species": "human", "target_species": "parrot", ...}
- "Does BRCA1 exist in tuna?" → {"action": "list_orthologs", "gene": "BRCA1", "species": "human", "target_species": "tuna", ...}
- "Previous names of TUG1" → {"action": "gene_history", "gene": "TUG1", "species": "human", ...}

IMPORTANT: When the user asks "is gene X in species Y" or "is X annotated in Y",
always use action "list_orthologs" with species="human" and target_species=Y.
Do NOT change the "species" field to Y — keep it as "human" and put Y in "target_species".

Species name mapping (always use Ensembl names in target_species):
- mouse → mus_musculus
- rat → rattus_norvegicus
- chimp/chimpanzee → pan_troglodytes
- dog → canis_lupus_familiaris
- chicken → gallus_gallus
- zebrafish → danio_rerio
- fly/drosophila → drosophila_melanogaster
- pig → sus_scrofa
- cow → bos_taurus
- macaque/rhesus → macaca_mulatta
- gorilla → gorilla_gorilla

For species NOT in this list (e.g. "puerto rico parrot", "naked mole rat", "platypus"),
put the common name as-is in target_species — it will be resolved automatically.
If the user asks whether a gene "has a homolog in" a species, use action "list_orthologs" 
with that species as target_species.

Return ONLY the JSON object."""


# ── Taxon name → ID mapping ──────────────────────────────────────────────

TAXON_MAP = {
    "mammals": 40674, "mammalia": 40674, "mammalian": 40674,
    "primates": 9443, "primate": 9443,
    "rodents": 9989, "rodentia": 9989, "rodent": 9989,
    "vertebrates": 7742, "vertebrata": 7742, "vertebrate": 7742,
    "fish": 7898, "fishes": 7898, "actinopterygii": 7898,
    "birds": 8782, "aves": 8782, "bird": 8782,
    "insects": 50557, "insecta": 50557, "insect": 50557,
    "plants": 33090, "viridiplantae": 33090, "plant": 33090,
    "metazoa": 33208, "animals": 33208,
    "fungi": 4751,
    "amphibians": 8292, "amphibia": 8292,
    "reptiles": 8504, "reptilia": 8504,
}


def resolve_taxon(value) -> Optional[int]:
    ## Resolve a taxon name or ID to NCBI taxon ID
    if value is None:
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        if value.isdigit():
            return int(value)
        return TAXON_MAP.get(value.lower().strip())
    return None


# ── LLM Parser ───────────────────────────────────────────────────────────

class LLMParser:
    # Parse natural language gene queries using a local Ollama LLM

    def __init__(self, model: str = "llama3.2"):
        self.model = model
        self._ollama = None

    def _client(self):
        if self._ollama is None:
            self._ollama = _get_ollama()
        return self._ollama

    def parse(self, query: str) -> dict:
        # Parse a natural language query into a structured intent dict
        try:
            response = self._client().chat(
                model=self.model,
                messages=[
                    {"role": "system", "content": SYSTEM_PROMPT},
                    {"role": "user", "content": query},
                ],
                options={"temperature": 0.0},
            )
            text = response["message"]["content"].strip()
            ## Extract JSON from response (handle markdown fences)
            json_match = re.search(r'\{[^{}]*\}', text, re.DOTALL)
            if json_match:
                parsed = json.loads(json_match.group())
            else:
                parsed = json.loads(text)

            ## Post-process: resolve taxon names to IDs
            if "target_taxon" in parsed and parsed["target_taxon"]:
                parsed["target_taxon"] = resolve_taxon(parsed["target_taxon"])

            ## Defaults
            parsed.setdefault("action", "unknown")
            parsed.setdefault("gene", None)
            parsed.setdefault("species", "human")
            parsed.setdefault("target_species", None)
            parsed.setdefault("target_taxon", None)
            parsed.setdefault("seq_type", "protein")
            parsed.setdefault("subtype", None)
            parsed.setdefault("external_db", None)

            return parsed

        except Exception as e:
            return {"action": "error", "error": str(e), "gene": None}


# ── Keyword Parser (fallback, no LLM needed) ─────────────────────────────

class KeywordParser:
    # Rule-based parser using regex and keyword matching

    # Gene patterns: Ensembl IDs or gene symbols (uppercase with digits/hyphens)
    GENE_PATTERN = re.compile(
        r'(ENS[A-Z]*G\d{11}(?:\.\d+)?)'
        r'|'
        r'\b([A-Z][A-Z0-9\-]{1,15})\b'
    )

    SPECIES_ALIASES = {
        "human": "homo_sapiens", "homo sapiens": "homo_sapiens", "hsa": "homo_sapiens",
        "mouse": "mus_musculus", "mus musculus": "mus_musculus", "mus_musculus": "mus_musculus", "mmu": "mus_musculus",
        "rat": "rattus_norvegicus", "rattus norvegicus": "rattus_norvegicus", "rattus_norvegicus": "rattus_norvegicus",
        "zebrafish": "danio_rerio", "danio rerio": "danio_rerio", "danio_rerio": "danio_rerio",
        "fly": "drosophila_melanogaster", "drosophila": "drosophila_melanogaster", "drosophila_melanogaster": "drosophila_melanogaster",
        "worm": "caenorhabditis_elegans", "c. elegans": "caenorhabditis_elegans", "caenorhabditis_elegans": "caenorhabditis_elegans",
        "chicken": "gallus_gallus", "gallus gallus": "gallus_gallus", "gallus_gallus": "gallus_gallus",
        "dog": "canis_lupus_familiaris", "canis familiaris": "canis_lupus_familiaris", "canis_lupus_familiaris": "canis_lupus_familiaris",
        "pig": "sus_scrofa", "sus scrofa": "sus_scrofa", "sus_scrofa": "sus_scrofa",
        "cow": "bos_taurus", "bos taurus": "bos_taurus", "bos_taurus": "bos_taurus",
        "chimp": "pan_troglodytes", "chimpanzee": "pan_troglodytes", "pan_troglodytes": "pan_troglodytes",
        "macaque": "macaca_mulatta", "rhesus": "macaca_mulatta", "macaca_mulatta": "macaca_mulatta",
        "gorilla": "gorilla_gorilla", "gorilla_gorilla": "gorilla_gorilla",
    }

    # Words that look like gene symbols but aren't
    STOP_WORDS = {
        "THE", "AND", "FOR", "GET", "HOW", "MANY", "SHOW", "LIST", "WHAT",
        "DOES", "HAVE", "HAS", "ARE", "ALL", "CAN", "YOU", "GIVE", "FIND",
        "ME", "OF", "IN", "TO", "IS", "IT", "DO", "TREE", "GENE", "IDS",
        "ID", "SEQUENCE", "SEQUENCES", "PROTEIN", "CDS", "CDNA", "DNA",
        "ORTHOLOG", "ORTHOLOGS", "PARALOG", "PARALOGS", "HOMOLOG", "HOMOLOGS",
        "UNIPROT", "REFSEQ", "XREF", "XREFS", "CROSS", "REFERENCES",
        "COUNT", "NUMBER", "RETRIEVE", "HELP", "ABOUT",
        "ORTHOLOGS", "PARALOGS", "HOMOLOGS", "ORTHOLOGOUS", "PARALOGOUS",
        "MAMMALS", "PRIMATES", "VERTEBRATES", "FISH", "BIRDS", "INSECTS",
        "PLANTS", "RODENTS", "FROM", "WITH", "THIS", "THAT",
    }

    def parse(self, query: str) -> dict:
        q = query.strip()
        ql = q.lower()

        result = {
            "action": "unknown",
            "gene": None,
            "species": "human",
            "target_species": None,
            "target_taxon": None,
            "seq_type": "protein",
            "subtype": None,
            "external_db": None,
        }

        ## Help
        if ql in ("help", "?", "what can you do", "commands"):
            result["action"] = "help"
            return result

        ## Extract gene
        for match in self.GENE_PATTERN.finditer(q):
            ensembl_id = match.group(1)
            symbol = match.group(2)
            gene = ensembl_id or symbol
            if gene and gene.upper() not in self.STOP_WORDS:
                result["gene"] = gene
                break

        ## Extract target species
        for alias, species in self.SPECIES_ALIASES.items():
            if alias in ql and alias != "human":
                result["target_species"] = species
                break

        ## Extract taxon
        for taxon_name, taxon_id in TAXON_MAP.items():
            if taxon_name in ql:
                result["target_taxon"] = taxon_id
                break

        ## Detect sequence type
        if "cds" in ql or "coding sequence" in ql:
            result["seq_type"] = "cds"
        elif "cdna" in ql or "nucleotide" in ql or "mrna" in ql:
            result["seq_type"] = "cdna"

        ## Detect subtype
        if "one2one" in ql or "one-to-one" in ql or "1:1" in ql or "1to1" in ql:
            result["subtype"] = "one2one"
        elif "one2many" in ql or "one-to-many" in ql or "1:many" in ql:
            result["subtype"] = "one2many"

        # Detect action
        if any(w in ql for w in ["list species", "all species", "available species", "which species", "supported species"]):
            result["action"] = "list_species"
        elif any(w in ql for w in ["previous name", "prev name", "gene history", "name history",
                                    "old name", "former name", "alias", "aliases",
                                    "previous symbol", "called before", "renamed"]):
            result["action"] = "gene_history"
        elif any(w in ql for w in ["how many", "count", "number of"]):
            result["action"] = "count_homologs"
        elif any(w in ql for w in ["uniprot"]):
            result["action"] = "get_xrefs"
            result["external_db"] = "UniProt"
        elif any(w in ql for w in ["refseq"]):
            result["action"] = "get_xrefs"
            result["external_db"] = "RefSeq"
        elif any(w in ql for w in ["cross-ref", "xref", "cross ref"]):
            result["action"] = "get_xrefs"
        elif any(w in ql for w in ["tree", "phylogen", "newick"]):
            result["action"] = "get_gene_tree"
        elif any(w in ql for w in ["sequence", "fasta", "protein seq", "cds seq"]):
            if any(w in ql for w in ["ortholog", "orthologue"]):
                result["action"] = "get_ortholog_sequences"
            else:
                result["action"] = "get_sequence"
        elif any(w in ql for w in ["paralog", "paralogue"]):
            result["action"] = "list_paralogs"
        elif any(w in ql for w in ["ortholog", "orthologue", "homolog", "homologue"]):
            result["action"] = "list_orthologs"
        elif any(w in ql for w in ["look up", "lookup", "info", "tell me about", "what is"]):
            if result["gene"]:
                result["action"] = "lookup"
        elif result["gene"]:
            result["action"] = "lookup"

        return result

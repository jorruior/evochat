# Main chat engine: ties together parser, Ensembl client, and formatter

from evochat.ensembl import EnsemblClient
from evochat.parser import LLMParser, KeywordParser, TAXON_MAP
from evochat import formatter


class EvoChat:
    # Gene evolution chatbot engine
    # Parses natural language → calls Ensembl REST API → formats output
    #
    # Usage:
    #   chat = EvoChat(use_llm=False)
    #   print(chat.ask("How many orthologs does BRCA1 have?"))

    def __init__(self, use_llm: bool = True, model: str = "llama3.2", verbose: bool = False):
        self.client = EnsemblClient()
        self.verbose = verbose

        if use_llm:
            try:
                self.parser = LLMParser(model=model)
                self.parser.parse("test")
                self.parser_type = "llm"
            except Exception as e:
                if verbose:
                    print(f"  [LLM unavailable: {e}]")
                    print(f"  [Falling back to keyword parser]")
                self.parser = KeywordParser()
                self.parser_type = "keyword"
        else:
            self.parser = KeywordParser()
            self.parser_type = "keyword"

    def ask(self, query: str) -> str:
        # Process a natural language query and return formatted response
        intent = self.parser.parse(query)

        # Post-process target_species: catch taxon names and resolve unknown species
        ts = intent.get("target_species")
        if ts and isinstance(ts, str):
            ts_lower = ts.lower().strip()
            if ts_lower in TAXON_MAP:
                intent["target_taxon"] = TAXON_MAP[ts_lower]
                intent["target_species"] = None
            else:
                aliases = KeywordParser.SPECIES_ALIASES
                if ts_lower in aliases:
                    intent["target_species"] = aliases[ts_lower]
                elif ts_lower.replace(" ", "_") not in {v for v in aliases.values()}:
                    ## Unknown species — search Ensembl
                    resolved = self.client.search_species(ts)
                    if resolved:
                        intent["target_species"] = resolved

        if self.verbose:
            import json
            print(f"  [Parser: {self.parser_type}]")
            print(f"  [Intent: {json.dumps(intent, indent=2)}]")

        action = intent.get("action", "unknown")
        gene = intent.get("gene")
        species = intent.get("species", "human")

        # Resolve source species if it's not standard
        if species and species != "human":
            species_lower = species.lower().strip()
            aliases = KeywordParser.SPECIES_ALIASES
            if species_lower in aliases:
                species = aliases[species_lower]
            elif species_lower in TAXON_MAP:
                ## Taxon name as source species → convert to target filter
                intent["target_taxon"] = TAXON_MAP[species_lower]
                intent["target_species"] = None
                species = "human"
                if action == "lookup":
                    action = "list_orthologs"
            else:
                resolved = self.client.search_species(species)
                if resolved:
                    # User asked about a gene in a non-human species
                    # Convert to "find ortholog in that species" query
                    intent["target_species"] = resolved
                    species = "human"
                    if action == "lookup":
                        action = "list_orthologs"
                else:
                    return f"  Species '{species}' not found in Ensembl. Try 'list all species' to see available species."

        # Dispatch to handler
        try:
            if action == "help":
                return formatter.HELP_TEXT

            if action == "error":
                return f"  Error parsing query: {intent.get('error', 'unknown')}"

            if action == "unknown":
                return (
                    "  I couldn't understand that query. Try something like:\n"
                    "    'How many orthologs does BRCA1 have?'\n"
                    "    'List orthologs of TP53 in mammals'\n"
                    "  Type 'help' for more examples."
                )

            if action == "list_species":
                return self._handle_list_species()

            if not gene:
                return "  Please specify a gene symbol or Ensembl ID."

            # Resolve gene symbol to Ensembl ID
            gene_id = self._resolve_gene(gene, species)
            gene_display = gene

            if action == "lookup":
                return self._handle_lookup(gene, species)
            elif action == "count_homologs":
                return self._handle_count(gene_display, gene_id)
            elif action == "list_orthologs":
                return self._handle_list_orthologs(gene_display, gene_id, intent)
            elif action == "list_paralogs":
                return self._handle_list_paralogs(gene_display, gene_id)
            elif action == "get_ortholog_sequences":
                return self._handle_ortholog_sequences(gene_display, gene_id, intent)
            elif action == "get_gene_tree":
                return self._handle_gene_tree(gene_display, gene_id)
            elif action == "get_xrefs":
                return self._handle_xrefs(gene_display, gene_id, intent)
            elif action == "get_sequence":
                return self._handle_sequence(gene_display, gene_id, intent)
            elif action == "gene_history":
                return self._handle_gene_history(gene_display, species)
            else:
                return f"  Action '{action}' is not yet implemented."

        except Exception as e:
            return f"  Error: {e}"

    # ── Gene resolution ──────────────────────────────────────────────────

    def _resolve_gene(self, gene: str, species: str = "human") -> str:
        # Resolve gene symbol to Ensembl ID
        # Falls back to HGNC when Ensembl lookup fails (old/alias symbols)
        if gene.startswith("ENS"):
            return gene
        try:
            info = self.client.lookup_gene(gene, species=species)
            return info.get("id", gene)
        except Exception:
            ## Symbol not found — try HGNC to resolve old/alias names
            if species in ("human", "homo_sapiens"):
                current_symbol = self._resolve_via_hgnc(gene)
                if current_symbol and current_symbol != gene:
                    try:
                        info = self.client.lookup_gene(current_symbol, species=species)
                        return info.get("id", current_symbol)
                    except Exception:
                        pass
            raise ValueError(
                f"Gene '{gene}' not found in Ensembl. "
                f"Try using the current symbol or an Ensembl ID."
            )

    def _resolve_via_hgnc(self, symbol: str) -> str:
        # Try to resolve a gene symbol via HGNC (previous symbols, aliases)
        import requests as req

        ## Try fetching by symbol directly
        try:
            resp = req.get(
                f"https://rest.genenames.org/fetch/symbol/{symbol}",
                headers={"Accept": "application/json"},
            )
            if resp.ok:
                docs = resp.json().get("response", {}).get("docs", [])
                if docs:
                    return docs[0].get("symbol", symbol)
        except Exception:
            pass

        ## Try by previous symbol
        try:
            resp = req.get(
                f"https://rest.genenames.org/search/prev_symbol/{symbol}",
                headers={"Accept": "application/json"},
            )
            if resp.ok:
                docs = resp.json().get("response", {}).get("docs", [])
                if docs:
                    return docs[0].get("symbol")
        except Exception:
            pass

        ## Try by alias
        try:
            resp = req.get(
                f"https://rest.genenames.org/search/alias_symbol/{symbol}",
                headers={"Accept": "application/json"},
            )
            if resp.ok:
                docs = resp.json().get("response", {}).get("docs", [])
                if docs:
                    return docs[0].get("symbol")
        except Exception:
            pass

        return None

    # ── Handlers ─────────────────────────────────────────────────────────

    def _handle_lookup(self, gene: str, species: str) -> str:
        try:
            info = self.client.lookup_gene(gene, species=species)
        except Exception:
            if species in ("human", "homo_sapiens"):
                current = self._resolve_via_hgnc(gene)
                if current and current != gene:
                    info = self.client.lookup_gene(current, species=species)
                    result = formatter.format_gene_info(info)
                    return f"  Note: '{gene}' was resolved to current symbol '{current}'\n\n{result}"
            raise
        return formatter.format_gene_info(info)

    def _handle_count(self, gene_display: str, gene_id: str) -> str:
        counts = self.client.count_homologs(gene_id)
        return formatter.format_homolog_counts(gene_display, counts)

    def _handle_list_orthologs(self, gene_display: str, gene_id: str, intent: dict) -> str:
        orthologs = self.client.get_orthologs(
            gene_id,
            target_species=intent.get("target_species"),
            target_taxon=intent.get("target_taxon"),
            subtype=intent.get("subtype"),
        )
        return formatter.format_ortholog_list(gene_display, orthologs)

    def _handle_list_paralogs(self, gene_display: str, gene_id: str) -> str:
        paralogs = self.client.get_paralogs(gene_id)
        return formatter.format_paralog_list(gene_display, paralogs)

    def _handle_ortholog_sequences(self, gene_display: str, gene_id: str, intent: dict) -> str:
        seq_type = intent.get("seq_type", "protein")
        sequences = self.client.get_ortholog_sequences(
            gene_id,
            seq_type=seq_type,
            target_species=intent.get("target_species"),
            target_taxon=intent.get("target_taxon"),
        )
        return formatter.format_sequences_fasta(gene_display, sequences)

    def _handle_gene_tree(self, gene_display: str, gene_id: str) -> str:
        newick = self.client.get_gene_tree_newick(gene_id)
        return formatter.format_gene_tree(newick, gene_display)

    def _handle_xrefs(self, gene_display: str, gene_id: str, intent: dict) -> str:
        ext_db = intent.get("external_db")
        if ext_db and "uniprot" in ext_db.lower():
            ids = self.client.get_uniprot_ids(gene_id)
            return formatter.format_uniprot(gene_display, ids)
        elif ext_db and "refseq" in ext_db.lower():
            ids = self.client.get_refseq_ids(gene_id)
            return formatter.format_refseq(gene_display, ids)
        else:
            xrefs = self.client.get_xrefs(gene_id, all_levels=True)
            return formatter.format_xrefs(gene_display, xrefs)

    def _handle_sequence(self, gene_display: str, gene_id: str, intent: dict) -> str:
        seq_type = intent.get("seq_type", "protein")

        # For gene IDs, resolve to transcript/translation ID first
        # The sequence/id endpoint on gene IDs only supports 'genomic'
        fetch_id = gene_id
        if gene_id.startswith("ENSG") or (gene_id.startswith("ENS") and "G" in gene_id[:7]):
            if seq_type in ("cdna", "cds", "protein"):
                info = self.client.lookup_gene(gene_id, expand=True)
                transcripts = info.get("Transcript", [])
                if transcripts:
                    canonical = next(
                        (t for t in transcripts if t.get("is_canonical")),
                        transcripts[0],
                    )
                    if seq_type == "protein":
                        fetch_id = canonical.get("Translation", {}).get("id", canonical["id"])
                        seq_type = "protein"
                    else:
                        fetch_id = canonical["id"]

        data = self.client.get_sequence(fetch_id, type=seq_type)
        seq = data.get("seq", "")
        mol = data.get("molecule", seq_type)
        unit = "aa" if seq_type == "protein" else "bp"
        header = f">{fetch_id} | {gene_display} | {mol}"
        wrapped = "\n".join(seq[i:i+80] for i in range(0, len(seq), 80))
        return f"  Sequence for {gene_display} ({seq_type}, {len(seq)} {unit}):\n\n{header}\n{wrapped}"

    def _handle_gene_history(self, gene_display: str, species: str) -> str:
        history = self.client.get_gene_history(gene_display, species=species)
        return formatter.format_gene_history(gene_display, history)

    def _handle_list_species(self) -> str:
        species_list = self.client.list_species()
        lines = []
        lines.append(f"  Species in Ensembl ({len(species_list)} total):")
        lines.append(f"")
        lines.append(f"  {'Ensembl Name':<40} {'Common Name':<30} {'Taxon ID'}")
        lines.append(f"  {'─'*40} {'─'*30} {'─'*10}")

        for sp in sorted(species_list, key=lambda x: x.get("name", "")):
            name = sp.get("name", "?")
            common = sp.get("common_name", sp.get("display_name", "?"))
            taxon = sp.get("taxon_id", "?")
            lines.append(f"  {name:<40} {common:<30} {taxon}")

        return "\n".join(lines)

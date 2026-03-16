# Tests for evochat keyword parser and formatter

from evochat.parser import KeywordParser, resolve_taxon
from evochat import formatter


class TestKeywordParser:

    def setup_method(self):
        self.parser = KeywordParser()

    def test_help(self):
        assert self.parser.parse("help")["action"] == "help"
        assert self.parser.parse("?")["action"] == "help"

    def test_count_homologs(self):
        result = self.parser.parse("How many orthologs does BRCA1 have?")
        assert result["action"] == "count_homologs"
        assert result["gene"] == "BRCA1"

    def test_count_with_ensembl_id(self):
        result = self.parser.parse("Count homologs of ENSG00000012048")
        assert result["action"] == "count_homologs"
        assert result["gene"] == "ENSG00000012048"

    def test_list_orthologs(self):
        result = self.parser.parse("List orthologs of TP53")
        assert result["action"] == "list_orthologs"
        assert result["gene"] == "TP53"

    def test_orthologs_in_mammals(self):
        result = self.parser.parse("Orthologs of LZTR1 in mammals")
        assert result["action"] == "list_orthologs"
        assert result["gene"] == "LZTR1"
        assert result["target_taxon"] == 40674

    def test_one2one_orthologs(self):
        result = self.parser.parse("1:1 orthologs of BRCA2")
        assert result["action"] == "list_orthologs"
        assert result["gene"] == "BRCA2"
        assert result["subtype"] == "one2one"

    def test_mouse_ortholog(self):
        result = self.parser.parse("Get mouse ortholog of TP53")
        assert result["action"] == "list_orthologs"
        assert result["gene"] == "TP53"
        assert result["target_species"] == "mus_musculus"

    def test_protein_sequences(self):
        result = self.parser.parse("Protein sequences of BRCA1 orthologs in primates")
        assert result["action"] == "get_ortholog_sequences"
        assert result["gene"] == "BRCA1"
        assert result["target_taxon"] == 9443
        assert result["seq_type"] == "protein"

    def test_cds_sequences(self):
        result = self.parser.parse("CDS sequences of TP53 orthologs")
        assert result["action"] == "get_ortholog_sequences"
        assert result["gene"] == "TP53"
        assert result["seq_type"] == "cds"

    def test_gene_tree(self):
        result = self.parser.parse("Gene tree of BRCA1")
        assert result["action"] == "get_gene_tree"
        assert result["gene"] == "BRCA1"

    def test_uniprot(self):
        result = self.parser.parse("UniProt IDs for LZTR1")
        assert result["action"] == "get_xrefs"
        assert result["gene"] == "LZTR1"
        assert result["external_db"] == "UniProt"

    def test_refseq(self):
        result = self.parser.parse("RefSeq IDs for TP53")
        assert result["action"] == "get_xrefs"
        assert result["gene"] == "TP53"
        assert result["external_db"] == "RefSeq"

    def test_paralogs(self):
        result = self.parser.parse("Paralogs of BRCA1")
        assert result["action"] == "list_paralogs"
        assert result["gene"] == "BRCA1"

    def test_lookup(self):
        result = self.parser.parse("Tell me about BRCA1")
        assert result["action"] == "lookup"
        assert result["gene"] == "BRCA1"

    def test_gene_history(self):
        result = self.parser.parse("Gene history of TUG1")
        assert result["action"] == "gene_history"
        assert result["gene"] == "TUG1"

    def test_previous_names(self):
        result = self.parser.parse("Previous names of LINC01405")
        assert result["action"] == "gene_history"
        assert result["gene"] == "LINC01405"

    def test_list_species(self):
        result = self.parser.parse("List all species")
        assert result["action"] == "list_species"

    def test_unknown_no_gene(self):
        result = self.parser.parse("What is the meaning of life?")
        assert result["action"] == "unknown"


class TestResolveTaxon:
    def test_name_to_id(self):
        assert resolve_taxon("mammals") == 40674
        assert resolve_taxon("primates") == 9443

    def test_int_passthrough(self):
        assert resolve_taxon(40674) == 40674

    def test_string_number(self):
        assert resolve_taxon("9443") == 9443

    def test_none(self):
        assert resolve_taxon(None) is None

    def test_unknown(self):
        assert resolve_taxon("aliens") is None


class TestFormatter:
    def test_format_gene_info(self):
        info = {
            "display_name": "BRCA1",
            "id": "ENSG00000012048",
            "species": "homo_sapiens",
            "description": "BRCA1 DNA repair associated [Source:HGNC]",
            "biotype": "protein_coding",
            "seq_region_name": "17",
            "start": 43044295,
            "end": 43170245,
            "strand": -1,
        }
        result = formatter.format_gene_info(info)
        assert "BRCA1" in result
        assert "ENSG00000012048" in result
        assert "Homo Sapiens" in result

    def test_format_homolog_counts(self):
        counts = {
            "orthologs_one2one": 43,
            "orthologs_one2many": 172,
            "orthologs_many2many": 0,
            "paralogs_within_species": 1,
            "paralogs_other": 0,
            "total_orthologs": 215,
            "total_paralogs": 1,
            "total": 216,
        }
        result = formatter.format_homolog_counts("BRCA1", counts)
        assert "43" in result
        assert "215" in result

    def test_format_empty_orthologs(self):
        result = formatter.format_ortholog_list("TEST", [])
        assert "No orthologs" in result

    def test_format_gene_history(self):
        history = {
            "current_symbol": "VHRT",
            "ensembl_id": "ENSG00000185847",
            "name": "ventricular heart development associated lncRNA",
            "previous_symbols": ["LINC01405"],
            "aliases": ["MASCC1"],
            "hgnc_id": "HGNC:50688",
            "status": "Approved",
        }
        result = formatter.format_gene_history("LINC01405", history)
        assert "VHRT" in result
        assert "renamed" in result
        assert "LINC01405" in result

    def test_help_text(self):
        assert "evochat" in formatter.HELP_TEXT
        assert "BRCA1" in formatter.HELP_TEXT

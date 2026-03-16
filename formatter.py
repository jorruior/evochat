# Response formatters for terminal output

from io import StringIO
from typing import Optional


def format_gene_info(info: dict) -> str:
    name = info.get("display_name", info.get("id", "?"))
    ensembl_id = info.get("id", "?")
    species = info.get("species", "?").replace("_", " ").title()
    desc = info.get("description", "").split(" [")[0]
    biotype = info.get("biotype", "?")
    chrom = info.get("seq_region_name", "?")
    start = info.get("start", "?")
    end = info.get("end", "?")
    strand = "+" if info.get("strand", 1) == 1 else "-"

    lines = [
        f"  Gene:        {name}",
        f"  Ensembl ID:  {ensembl_id}",
        f"  Species:     {species}",
        f"  Description: {desc}",
        f"  Biotype:     {biotype}",
        f"  Location:    chr{chrom}:{start}-{end} ({strand})",
    ]
    return "\n".join(lines)


def format_homolog_counts(gene: str, counts: dict) -> str:
    lines = [
        f"  Homolog counts for {gene}:",
        f"",
        f"  Orthologs:",
        f"    1:1 orthologs      {counts.get('orthologs_one2one', 0):>6}",
        f"    1:many orthologs   {counts.get('orthologs_one2many', 0):>6}",
        f"    many:many orthologs{counts.get('orthologs_many2many', 0):>6}",
        f"    ─────────────────────────",
        f"    Total orthologs    {counts.get('total_orthologs', 0):>6}",
        f"",
        f"  Paralogs:",
        f"    Within species     {counts.get('paralogs_within_species', 0):>6}",
        f"    Other              {counts.get('paralogs_other', 0):>6}",
        f"    ─────────────────────────",
        f"    Total paralogs     {counts.get('total_paralogs', 0):>6}",
        f"",
        f"  Total homologs       {counts.get('total', 0):>6}",
    ]
    return "\n".join(lines)


def format_ortholog_list(gene: str, orthologs: list, max_show: int = 30) -> str:
    if not orthologs:
        return f"  No orthologs found for {gene}."

    lines = [
        f"  Orthologs of {gene} ({len(orthologs)} total):",
        f"",
        f"  {'Species':<35} {'Gene ID':<22} {'Symbol':<15} {'Type'}",
        f"  {'─'*35} {'─'*22} {'─'*15} {'─'*18}",
    ]

    for h in orthologs[:max_show]:
        ## Handle both full and condensed format structures
        target = h.get("target", h)
        species = (
            target.get("species") or h.get("target_species") or "?"
        ).replace("_", " ").title()
        gid = target.get("id") or target.get("gene_id") or h.get("id") or "?"
        symbol = target.get("protein_id") or target.get("perc_id") or "?"
        htype = h.get("type", "?").replace("ortholog_", "")
        lines.append(f"  {species:<35} {gid:<22} {symbol:<15} {htype}")

    if len(orthologs) > max_show:
        lines.append(f"  ... and {len(orthologs) - max_show} more.")

    return "\n".join(lines)


def format_paralog_list(gene: str, paralogs: list) -> str:
    if not paralogs:
        return f"  No paralogs found for {gene}."

    lines = [
        f"  Paralogs of {gene} ({len(paralogs)} total):",
        f"",
        f"  {'Gene ID':<22} {'Type':<30}",
        f"  {'─'*22} {'─'*30}",
    ]
    for h in paralogs:
        target = h.get("target", {})
        gid = target.get("id", "?")
        htype = h.get("type", "?")
        lines.append(f"  {gid:<22} {htype}")

    return "\n".join(lines)


def format_sequences_fasta(gene: str, sequences: list) -> str:
    if not sequences:
        return f"  No sequences retrieved for {gene} orthologs."

    buf = StringIO()
    count = 0
    for s in sequences:
        seq = s.get("sequence")
        if not seq:
            continue
        species = s.get("species", "unknown").replace("_", " ")
        gid = s.get("gene_id", "?")
        pid = s.get("protein_id", "?")
        htype = s.get("type", "?")
        buf.write(f">{gid} | {pid} | {species} | {htype}\n")
        for i in range(0, len(seq), 80):
            buf.write(seq[i:i+80] + "\n")
        count += 1

    header = f"  {count} ortholog sequences for {gene}:\n\n"
    return header + buf.getvalue()


def format_gene_tree(newick: str, gene: str) -> str:
    lines = [f"  Gene tree for {gene} (Newick format):", f""]
    if len(newick) > 2000:
        lines.append(f"  {newick[:2000]}...")
        lines.append(f"")
        lines.append(f"  [Tree truncated — {len(newick)} characters total]")
        lines.append(f"  Tip: save to file and visualize with FigTree, iTOL, or ETE3.")
    else:
        lines.append(f"  {newick}")
    return "\n".join(lines)


def format_xrefs(gene: str, xrefs: list, db_filter: Optional[str] = None) -> str:
    if not xrefs:
        return f"  No cross-references found for {gene}."

    by_db = {}
    for x in xrefs:
        db = x.get("dbname", "unknown")
        by_db.setdefault(db, []).append(x)

    lines = [f"  Cross-references for {gene}:", f""]
    for db, entries in sorted(by_db.items()):
        ids = sorted({e.get("primary_id", "?") for e in entries})
        lines.append(f"  {db}:")
        for pid in ids:
            lines.append(f"    {pid}")
        lines.append(f"")

    return "\n".join(lines)


def format_uniprot(gene: str, ids: list) -> str:
    if not ids:
        return f"  No UniProt (Swiss-Prot) IDs found for {gene}."
    lines = [f"  UniProt Swiss-Prot IDs for {gene}:"]
    for uid in sorted(ids):
        lines.append(f"    {uid}  →  https://www.uniprot.org/uniprot/{uid}")
    return "\n".join(lines)


def format_refseq(gene: str, ids: dict) -> str:
    lines = [f"  RefSeq IDs for {gene}:"]
    if ids.get("mRNA"):
        lines.append(f"  mRNA:")
        for rid in sorted(ids["mRNA"]):
            lines.append(f"    {rid}")
    if ids.get("peptide"):
        lines.append(f"  Peptide:")
        for rid in sorted(ids["peptide"]):
            lines.append(f"    {rid}")
    if not ids.get("mRNA") and not ids.get("peptide"):
        lines.append(f"  No RefSeq IDs found.")
    return "\n".join(lines)


def format_gene_history(gene: str, history: dict) -> str:
    lines = [f"  Gene history for {gene}:", f""]

    current = history.get("current_symbol", "?")
    ensembl_id = history.get("ensembl_id", "?")
    name = history.get("name", "?")
    hgnc_id = history.get("hgnc_id", "?")
    status = history.get("status", "?")
    prev = history.get("previous_symbols", [])
    aliases = history.get("aliases", [])

    lines.append(f"  Current symbol:     {current}")
    if current and gene and current.upper() != gene.upper():
        lines.append(f"  ⚠ Note:             '{gene}' has been renamed to '{current}'")
    lines.append(f"  Ensembl ID:         {ensembl_id}")
    lines.append(f"  Full name:          {name}")
    lines.append(f"  HGNC ID:            {hgnc_id}")
    lines.append(f"  Status:             {status}")
    lines.append(f"")

    if prev:
        lines.append(f"  Previous symbols:   {', '.join(prev)}")
    else:
        lines.append(f"  Previous symbols:   (none)")

    if aliases:
        lines.append(f"  Aliases:            {', '.join(aliases)}")
    else:
        lines.append(f"  Aliases:            (none)")

    archive = history.get("ensembl_archive")
    if archive:
        lines.append(f"")
        lines.append(f"  Ensembl archive:")
        lines.append(f"    Latest version:   {archive.get('version', '?')}")
        lines.append(f"    Release:          {archive.get('release', '?')}")
        lines.append(f"    Assembly:         {archive.get('assembly', '?')}")
        lines.append(f"    Current:          {archive.get('is_current', '?')}")

    return "\n".join(lines)


HELP_TEXT = """
  ╔══════════════════════════════════════════════════════════════╗
  ║                    evochat — Gene Evolution Chat             ║
  ╠══════════════════════════════════════════════════════════════╣
  ║                                                              ║
  ║  Example queries:                                            ║
  ║                                                              ║
  ║  Gene info:                                                  ║
  ║    "Tell me about BRCA1"                                     ║
  ║    "Lookup ENSG00000012048"                                  ║
  ║                                                              ║
  ║  Homolog counts:                                             ║
  ║    "How many orthologs does TP53 have?"                      ║
  ║    "Count homologs of LZTR1"                                 ║
  ║                                                              ║
  ║  Ortholog lists:                                             ║
  ║    "List orthologs of BRCA2"                                 ║
  ║    "Get mouse ortholog of TP53"                              ║
  ║    "1:1 orthologs of LZTR1 in mammals"                      ║
  ║                                                              ║
  ║  Sequences:                                                  ║
  ║    "Protein sequences of BRCA1 orthologs in primates"        ║
  ║    "CDS sequence of ENSG00000012048"                         ║
  ║                                                              ║
  ║  Gene trees:                                                 ║
  ║    "Gene tree of TP53"                                       ║
  ║                                                              ║
  ║  Cross-references:                                           ║
  ║    "UniProt IDs for BRCA1"                                   ║
  ║    "RefSeq IDs for TP53"                                     ║
  ║                                                              ║
  ║  Gene history:                                               ║
  ║    "Gene history of LINC01405"                               ║
  ║    "Previous names of TUG1"                                  ║
  ║                                                              ║
  ║  Species:                                                    ║
  ║    "List all species"                                        ║
  ║    "Is TUG1 annotated in parrot?"                            ║
  ║                                                              ║
  ║  Commands:  help | quit | exit                               ║
  ╚══════════════════════════════════════════════════════════════╝
"""

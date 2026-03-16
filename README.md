# 🧬 evochat (v0.1.0)

An LLM-powered chatbot for gene evolution queries, backed by the [Ensembl REST API](https://rest.ensembl.org) and [HGNC](https://www.genenames.org).

Ask natural language questions about orthologs, paralogs, gene trees, sequences, cross-references, and gene name history from your terminal.

## Features

- **Natural language queries** — parsed by a local LLM (Ollama) or keyword-based fallback
- **Ortholog/paralog counts** — broken down by type (1:1, 1:many, etc.)
- **Ortholog listing** — with species, Ensembl IDs, and relationship type
- **Sequence retrieval** — protein or CDS sequences for orthologs, in FASTA format
- **Gene trees** — Newick-format phylogenetic trees
- **Cross-references** — UniProt, RefSeq, HGNC, and more
- **Taxon filtering** — filter by mammals, primates, vertebrates, etc.
- **Species search** — resolve common names ("human", "tuna") against Ensembl's species database
- **Gene name history** — previous symbols and aliases via HGNC (e.g. LINC01405 → VHRT)
- **Old symbol resolution** — retired gene names are automatically resolved to current symbols

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/jorruior/evochat.git
cd evochat
```

### 2. Create a conda environment (recommended)

```bash
conda create -n evochat python=3.11
conda activate evochat
```

### 3. Install evochat

```bash
pip install -e .
```

This installs `requests`, `rich`, and `ollama` Python packages automatically.

### 4. Install Ollama (optional, for LLM-powered parsing)

Ollama runs LLMs locally. Without it, evochat uses keyword-based parsing which works well for structured queries.

**On a machine with sudo access:**

```bash
curl -fsSL https://ollama.com/install.sh | sh
```

**On HPC clusters without sudo (e.g. SLURM):**

```bash
# Download the tarball
curl -L https://ollama.com/download/ollama-linux-amd64.tar.zst -o /tmp/ollama-linux-amd64.tar.zst

# Extract to a local directory (adjust path as needed)
mkdir -p /fast/your_group/your_user/ollama
tar --use-compress-program=unzstd -xf /tmp/ollama-linux-amd64.tar.zst -C /fast/your_group/your_user/ollama

# Add to PATH and set model storage location (add to ~/.bashrc)
export PATH="/fast/your_group/your_user/ollama/bin:$PATH"
export OLLAMA_MODELS="/fast/your_group/your_user/.ollama/models"
mkdir -p $OLLAMA_MODELS
```

### 5. Start Ollama and pull a model

```bash
# Start the server (redirect logs to avoid terminal flooding)
ollama serve > /dev/null 2>&1 &
sleep 5

# Pull a small, fast model
ollama pull llama3.2
```

## Quick Start

### Interactive chat

```bash
evochat
```

```
evochat> How many orthologs does BRCA1 have?

  Homolog counts for BRCA1:

  Orthologs:
    1:1 orthologs         170
    1:many orthologs       29
    many:many orthologs     0
    ─────────────────────────
    Total orthologs       199

evochat> Gene history of LINC01405

  Gene history for LINC01405:

  Current symbol:     VHRT
  ⚠ Note:             'LINC01405' has been renamed to 'VHRT'
  Ensembl ID:         ENSG00000185847
  Full name:          ventricular heart development associated lncRNA
  HGNC ID:            HGNC:50688
  Status:             Approved

  Previous symbols:   LINC01405
  Aliases:            MASCC1, RP1-46F2.2

evochat> Is TUG1 annotated in parrot?

  Orthologs of TUG1 (1 total):
  ...
```

### Single query mode

```bash
evochat -q "List orthologs of TP53 in mammals"
```

### Without LLM (keyword parsing only)

```bash
evochat --no-llm
```

### Python API

```python
from evochat import EvoChat, EnsemblClient

# High-level chat interface
chat = EvoChat(use_llm=False)
print(chat.ask("How many orthologs does BRCA1 have?"))

# Low-level Ensembl client
client = EnsemblClient()
counts = client.count_homologs("ENSG00000012048")
orthologs = client.get_orthologs("ENSG00000012048", target_taxon=40674)
seqs = client.get_ortholog_sequences("ENSG00000012048", seq_type="protein")
tree = client.get_gene_tree_newick("ENSG00000012048")
uniprot = client.get_uniprot_ids("ENSG00000012048")
history = client.get_gene_history("LINC01405")
```

## CLI Options

| Option | Description |
|---|---|
| `--no-llm` | Use keyword parser instead of LLM |
| `--model MODEL` | Ollama model name (default: `llama3.2`) |
| `-v, --verbose` | Show parsed intents for debugging |
| `-q, --query TEXT` | Run a single query and exit |
| `--version` | Show version |

## Supported Queries

| Query type | Example |
|---|---|
| Gene info | `"Tell me about BRCA1"` |
| Homolog counts | `"How many orthologs does TP53 have?"` |
| Ortholog list | `"List orthologs of LZTR1 in mammals"` |
| 1:1 orthologs | `"1:1 orthologs of BRCA2"` |
| Paralog list | `"Paralogs of TP53"` |
| Protein sequences | `"Protein sequences of BRCA1 orthologs in primates"` |
| CDS sequences | `"CDS sequence of ENSG00000012048"` |
| Gene tree | `"Gene tree of TP53"` |
| UniProt IDs | `"UniProt IDs for BRCA1"` |
| RefSeq IDs | `"RefSeq IDs for TP53"` |
| All cross-refs | `"Cross-references for LZTR1"` |
| Gene history | `"Gene history of LINC01405"` |
| Previous names | `"Previous names of TUG1"` |
| Species check | `"Is TUG1 annotated in mouse?"` |
| List species | `"List all species"` |

## Taxon Filters

| Name | NCBI Taxon ID |
|---|---|
| mammals / mammalia | 40674 |
| primates | 9443 |
| rodents | 9989 |
| vertebrates | 7742 |
| fish | 7898 |
| birds / aves | 8782 |
| insects | 50557 |
| plants | 33090 |
| amphibians | 8292 |
| reptiles | 8504 |

## Data Sources

- **[Ensembl REST API](https://rest.ensembl.org)** — homology, gene trees, sequences, cross-references, species info
- **[HGNC REST API](https://rest.genenames.org)** — gene symbol history, previous names, aliases (human genes only)

## Requirements

- Python ≥ 3.9
- `requests` — HTTP client
- `rich` — terminal formatting
- `ollama` — local LLM interface (optional)
- [Ollama server](https://ollama.ai) — for LLM-powered parsing (optional)

## License

MIT

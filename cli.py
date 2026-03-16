# Interactive command-line interface for evochat

import argparse
import sys
from importlib.metadata import version as _get_version

from rich.console import Console
from rich.panel import Panel
from rich.text import Text

from evochat.chat import EvoChat

try:
    __version__ = _get_version("evochat")
except Exception:
    __version__ = "0.1.0"


BANNER = r"""
    ┌─────────────────────────────────────────────────┐
    │   🧬  [bold green]evochat[/bold green] v{version}                            │
    │   Gene Evolution Chat — powered by Ensembl      │
    │                                                 │
    │   Type a question about any gene, or 'help'     │
    │   Type 'quit' or 'exit' to leave                │
    └─────────────────────────────────────────────────┘
"""


def main():
    parser = argparse.ArgumentParser(
        prog="evochat",
        description="Interactive chatbot for gene evolution queries via Ensembl",
    )
    parser.add_argument(
        "--no-llm",
        action="store_true",
        help="Disable LLM parsing, use keyword-based parser only",
    )
    parser.add_argument(
        "--model",
        default="llama3.2",
        help="Ollama model to use (default: llama3.2)",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show parsed intent for debugging",
    )
    parser.add_argument(
        "--query", "-q",
        type=str,
        default=None,
        help="Run a single query and exit (non-interactive mode)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"evochat {__version__}",
    )
    args = parser.parse_args()

    console = Console()

    # Single query mode
    if args.query:
        chat = EvoChat(use_llm=not args.no_llm, model=args.model, verbose=args.verbose)
        result = chat.ask(args.query)
        console.print(result)
        return

    # Interactive mode
    console.print(BANNER.format(version=__version__))

    chat = EvoChat(use_llm=not args.no_llm, model=args.model, verbose=args.verbose)
    console.print(f"  Parser: [cyan]{chat.parser_type}[/cyan]")
    console.print()

    while True:
        try:
            query = console.input("[bold green]evochat>[/bold green] ").strip()
        except (KeyboardInterrupt, EOFError):
            console.print("\n  Goodbye! 🧬")
            break

        if not query:
            continue

        if query.lower() in ("quit", "exit", "q"):
            console.print("  Goodbye! 🧬")
            break

        with console.status("[dim]Querying Ensembl...[/dim]", spinner="dots"):
            result = chat.ask(query)

        console.print()
        console.print(result)
        console.print()


if __name__ == "__main__":
    main()

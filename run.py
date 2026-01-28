#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import random
import shutil
import time
from logging.handlers import RotatingFileHandler
from types import SimpleNamespace
from typing import Any, ClassVar, Optional, Tuple, Union

try:
    from .constants import Constants
    from .modules.make_query_table import filter_query_annotation, table_builder
except ImportError:
    from constants import Constants
    from modules.make_query_table import filter_query_annotation, table_builder

try:
    from rustools import convert
except ImportError as exc:
    raise ImportError(
        "ERROR: rustools extension not found; run `make configure` or "
        "`maturin develop` to build it."
    ) from exc

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.13"


class PostogaLogger:
    """Configure a shared console/file logger for the CLI."""

    NAME = "postoga"
    _LOG_FILE_ATTR = "_postoga_log_file"
    LEVELS: ClassVar = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warn": logging.WARNING,
        "off": logging.CRITICAL + 1,
    }

    @staticmethod
    def _formatter() -> logging.Formatter:
        return logging.Formatter(
            "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

    @classmethod
    def _resolve_level(cls, level: Optional[str]) -> Tuple[int, bool]:
        normalized = (level or "info").lower()
        if normalized not in cls.LEVELS:
            normalized = "info"
        level_value = cls.LEVELS[normalized]
        disabled = normalized == "off"
        return level_value, disabled

    @classmethod
    def _apply_level(
        cls, logger: logging.Logger, level_value: int, disabled: bool
    ) -> None:
        logger.setLevel(level_value)
        logger.disabled = disabled
        for handler in logger.handlers:
            handler.setLevel(level_value)

    @classmethod
    def get_logger(
        cls,
        log_file: Optional[Union[str, os.PathLike]] = None,
        level: Optional[str] = None,
    ) -> logging.Logger:
        level_value, disabled = cls._resolve_level(level)
        logger = logging.getLogger(cls.NAME)
        if not logger.handlers:
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(cls._formatter())
            logger.addHandler(stream_handler)
            logger.propagate = False

        cls._apply_level(logger, level_value, disabled)

        if log_file:
            log_file = os.path.abspath(log_file)
            if not any(
                isinstance(handler, RotatingFileHandler)
                and getattr(handler, cls._LOG_FILE_ATTR, None) == log_file
                for handler in logger.handlers
            ):
                file_handler = RotatingFileHandler(
                    log_file,
                    maxBytes=5_000_000,
                    backupCount=2,
                )
                file_handler.setFormatter(cls._formatter())
                setattr(file_handler, cls._LOG_FILE_ATTR, log_file)
                file_handler.setLevel(level_value)
                logger.addHandler(file_handler)
            else:
                for handler in logger.handlers:
                    if (
                        isinstance(handler, RotatingFileHandler)
                        and getattr(handler, cls._LOG_FILE_ATTR, None) == log_file
                    ):
                        handler.setLevel(level_value)

        return logger


class TogaDir:
    """A class to represent a TOGA results directory."""

    def __init__(
        self, args: Optional[argparse.Namespace] = None, **kwargs: Any
    ) -> None:
        """
        Constructs all the necessary attributes for the TogaDir object.

        A TogaDir object should expect the following dir structure from
        the user:

            ├── inactivating_mutations.tsv
            ├── loss_summary.tsv
            ├── orthology_classification.tsv
            ├── orthology_scores.tsv
            |── query_genes.tsv
            ├── protein_aln.fa.gz [ only if --extract is set ]
            ├── nucleotide.fa.gz [ only if --extract is set ]
            ├── query_annotation.bed
            └── query_annotation.with_utrs.bed

        Parameters:
        ----------
        args : argparse.Namespace
            The argparse namespace object with all the arguments.

        Returns:
        ----------
        None
        """
        if args is not None and kwargs:
            raise ValueError(
                "ERROR: Pass either 'args' or kwargs to TogaDir, not both."
            )

        if args is None:
            if not kwargs:
                raise ValueError(
                    "ERROR: You must provide either 'args' or kwargs to TogaDir"
                )
            # Behave like argparse.Namespace
            args = SimpleNamespace(**kwargs)

        self.log_level = args.log_level
        self.logger = PostogaLogger.get_logger(level=self.log_level)
        self.args = args
        self.togadir = args.togadir

        self.hash = self.__c_hash__()
        self.timestamp = time.strftime("%Y%m%d_%H%M%S")

        if args.outdir is None:
            self.outdir: Union[str, os.PathLike] = os.path.abspath(
                os.path.join(
                    self.togadir,
                    f"{Constants.FileNames.POSTOGA_OUTPUT_DIRECTORY}_{self.hash}_{self.timestamp}",
                )
            )
        else:
            self.outdir = os.path.abspath(
                os.path.join(
                    args.outdir,
                    f"{Constants.FileNames.POSTOGA_OUTPUT_DIRECTORY}_{self.hash}_{self.timestamp}",
                )
            )

        self.logger.debug(
            "Initialized run with outdir=%s, togadir=%s", self.outdir, self.togadir
        )

        # INFO: setting filters
        self.by_orthology_class = args.orthology_class
        self.by_loss_status = args.loss_status
        self.by_orthology_score = args.orthology_score
        self.by_paralog_score = args.min_paralog_score

        # INFO: setting modes
        self.extract = args.extract
        self.only_table = args.only_table
        self.only_convert = args.only_convert
        self.convert_to = args.to
        self.isoforms = args.with_isoforms
        self.depure = args.depure
        self.filter_bed = any(
            [
                self.by_orthology_class,
                self.by_loss_status,
                self.by_orthology_score,
                self.by_paralog_score,
            ]
        )

        # INFO: .bed class attributes [ mandatory ]
        if (
            os.path.isfile(os.path.join(self.togadir, Constants.FileNames.BED_UTR))
            and args.bed_type == "utr"
        ):
            self.query_annotation = os.path.join(
                self.togadir, Constants.FileNames.BED_UTR
            )
        elif (
            os.path.isfile(os.path.join(self.togadir, Constants.FileNames.BED))
            and args.bed_type == "bed"
        ):
            self.query_annotation = os.path.join(self.togadir, Constants.FileNames.BED)
        else:
            self.logger.error("No .bed file found in %s", self.togadir)
            raise FileNotFoundError(f"ERROR: no .bed file found in {self.togadir}!")

        # INFO: class attributes to build toga.table.gz
        # [ only mandatory if not --with-isoforms ]
        self.loss_summary = os.path.join(self.togadir, Constants.FileNames.LOSS)
        self.orthology_classification = os.path.join(
            self.togadir, Constants.FileNames.ORTHOLOGY
        )
        self.orthology_scores = os.path.join(self.togadir, Constants.FileNames.SCORES)
        self.query_genes = os.path.join(self.togadir, Constants.FileNames.QUERY_GENES)

        # INFO: asserts existence of required files
        self.__check()

    def __check(self) -> None:
        """
        Checks if the required files are present in the toga directory.

        Raises:
        -------
        FileNotFoundError
            If any of the required files is not present.

        Returns:
        --------
        None
        """
        if not os.path.isdir(self.togadir):
            self.logger.error("TOGA directory %s is missing", self.togadir)
            raise FileNotFoundError(f"ERROR: {self.togadir} is not a directory!")

        # INFO: .bed file assertion
        if not os.path.isfile(self.query_annotation):
            self.logger.error("Missing query annotation %s", self.query_annotation)
            raise FileNotFoundError(f"ERROR: {self.query_annotation} is not a file!")
        if not os.path.isfile(self.loss_summary):
            self.logger.error("Missing loss summary %s", self.loss_summary)
            raise FileNotFoundError(f"ERROR: {self.loss_summary} is not a file!")
        if not os.path.isfile(self.orthology_classification):
            self.logger.error(
                "Missing orthology classification %s", self.orthology_classification
            )
            raise FileNotFoundError(
                f"ERROR: {self.orthology_classification} is not a file!"
            )
        if not os.path.isfile(self.orthology_scores):
            self.logger.error("Missing orthology scores %s", self.orthology_scores)
            raise FileNotFoundError(f"ERROR: {self.orthology_scores} is not a file!")
        if not os.path.isfile(self.query_genes):
            self.logger.error("Missing query genes %s", self.query_genes)
            raise FileNotFoundError(f"ERROR: {self.query_genes} is not a file!")

        # INFO: isoforms assertion
        if self.isoforms is not None:
            if not os.path.isfile(self.isoforms):
                self.logger.error("Custom isoform file %s not found", self.isoforms)
                raise FileNotFoundError(f"ERROR: {self.isoforms} is not a file!")

        # INFO: extract assertion
        if self.extract:
            if not os.path.isfile(
                os.path.join(self.togadir, Constants.FileNames.CODON)
            ):
                self.logger.error(
                    "Extract requested but codon alignment %s is missing",
                    os.path.join(self.togadir, Constants.FileNames.CODON),
                )
                raise FileNotFoundError(
                    f"""ERROR: {
                        os.path.join(self.togadir, Constants.FileNames.CODON)
                    } is not a file!"""
                )
            if not os.path.isfile(
                os.path.join(self.togadir, Constants.FileNames.PROTEIN)
            ):
                self.logger.error(
                    "Extract requested but protein alignment %s is missing",
                    os.path.join(self.togadir, Constants.FileNames.PROTEIN),
                )
                raise FileNotFoundError(
                    f"""ERROR: {
                        os.path.join(self.togadir, Constants.FileNames.PROTEIN)
                    } is not a file!"""
                )

        return

    def __c_hash__(self) -> str:
        """
        Returns a 5-digit hash of the current run

        Returns:
        ----------
        int
            The hash of the current run
        """
        CHARSET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        return "".join(random.choices(CHARSET, k=5))

    def __depure(self) -> None:
        """
        Remove any trace of other postoga runs/files

        Returns:
        ----------
        None
        """
        parent_dir = os.path.dirname(os.path.abspath(self.outdir))
        patterns = [
            f"{Constants.FileNames.POSTOGA_OUTPUT_DIRECTORY}_*",
            Constants.FileNames.POSTOGA_OUTPUT_DIRECTORY,
            Constants.FileNames.POSTOGA_OUTPUT_DIRECTORY.lower(),
            Constants.FileNames.LOG,
            Constants.FileNames.TOGA_TABLE,
            f"{Constants.FileNames.TOGA_TABLE}.gz",
            Constants.FileNames.FILTERED_GTF,
            Constants.FileNames.GTF,
            Constants.FileNames.FILTERED_GFF,
            Constants.FileNames.GFF,
            Constants.FileNames.FRAGMENTED_BED,
            Constants.FileNames.FILTERED_BED,
        ]

        self.logger.info("Depuring previous outputs in %s", parent_dir)
        for pattern in patterns:
            for candidate in glob.glob(os.path.join(parent_dir, pattern)):
                if not os.path.exists(candidate):
                    continue
                try:
                    if os.path.isdir(candidate):
                        shutil.rmtree(candidate)
                        self.logger.debug("Removed directory %s", candidate)
                    else:
                        os.remove(candidate)
                        self.logger.debug("Removed file %s", candidate)
                except OSError as exc:
                    self.logger.warning(
                        "Could not remove %s while depuring: %s", candidate, exc
                    )

        return

    def __write_isoforms(self) -> Union[str, os.PathLike]:
        """
        Write isoforms to a file.

        Args:
            isoforms: Isoforms to write

        Returns:
            None
        """
        target = os.path.join(self.outdir, "isoforms.tsv")
        self.logger.debug("Writing computed isoforms to %s", target)
        self.isoforms[["query_gene", "id"]].to_csv(
            target,
            sep="\t",
            index=False,
            header=False,
        )

        return target

    def __write_table(self) -> None:
        """
        Write toga.table to a file.

        Returns:
        ----------
        None
        """
        target = os.path.join(self.outdir, Constants.FileNames.TOGA_TABLE)
        self.logger.info("Writing %s with %d rows", target, len(self.table))
        self.table.to_csv(
            target,
            sep="\t",
            index=False,
            header=True,
            compression="gzip",
        )

    def __convert(self) -> None:
        """
        Convert toga.table to gtf/gff.

        Returns:
        ----------
        None
        """
        if self.convert_to == "gtf":
            self.logger.info("Converting BED to GTF")
            self.gtf = os.path.join(
                self.outdir,
                f"{os.path.splitext(os.path.basename(self.query_annotation))[0]}.gtf.gz",
            )
            self.gene_model = convert(self.query_annotation, self.gtf, self.isoforms)
        elif self.convert_to == "gff":
            self.logger.info("Converting BED to GFF")
            self.gff = os.path.join(
                self.outdir,
                f"{os.path.splitext(os.path.basename(self.query_annotation))[0]}.gff.gz",
            )
            self.gene_model = convert(self.query_annotation, self.gff, self.isoforms)
        else:
            self.logger.warning(
                "Conversion target %s is not supported", self.convert_to
            )

        return

    def run(self) -> None:
        """
        Run the postoga pipeline.

        Returns:
        ----------
        None

        Examples:
        ----------
        >>> TogaDir(args).run()
        """
        if self.depure:
            self.__depure()

        os.makedirs(self.outdir, exist_ok=True)
        self.logger = PostogaLogger.get_logger(
            os.path.join(self.outdir, Constants.FileNames.LOG), level=self.log_level
        )
        self.logger.info("postoga started with args: %s", vars(self.args))

        self.table, isoforms = table_builder(
            self.query_annotation,
            self.loss_summary,
            self.orthology_classification,
            self.orthology_scores,
            self.query_genes,
        )
        self.logger.info("Loaded table with %d projections", len(self.table))

        if not self.isoforms:
            self.isoforms = isoforms

            if not self.only_table:
                self.isoforms = self.__write_isoforms()
        else:
            self.logger.debug("Using user-supplied isoforms at %s", self.isoforms)

        if any(
            [
                self.by_orthology_class,
                self.by_loss_status,
                self.by_orthology_score,
                self.by_paralog_score,
            ]
        ):
            self.filtered_query_annotation, self.custom_table = filter_query_annotation(
                self.table,
                self.by_orthology_class,
                self.by_loss_status,
                self.by_orthology_score,
                self.by_paralog_score,
                self.query_annotation,
                self.outdir,
            )

            self.query_annotation = self.filtered_query_annotation
            self.table = self.custom_table
            self.logger.info(
                "(class=%s, rel=%s, score=%s, paralog=%s); %d rows remaining",
                self.by_orthology_class,
                self.by_loss_status,
                self.by_orthology_score,
                self.by_paralog_score,
                len(self.table),
            )
        else:
            self.logger.debug("No filters applied to query annotation")

        if self.only_table:
            self.logger.info("only-table flag detected; skipping conversion step")
            self.__write_table()
            return
        elif self.only_convert:
            self.logger.info("only-convert flag detected; skipping table rebuild")
            self.__convert()
            return
        else:
            self.__write_table()
            self.__convert()
            self.logger.info("postoga run finished successfully")

        return


def parse_args() -> argparse.Namespace:
    """Argument parser for postoga"""
    app = argparse.ArgumentParser(add_help=False)
    app.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )

    app.add_argument(
        "--togadir",
        "-td",
        help="Path to TOGA results directory",
        required=True,
        type=str,
    )
    app.add_argument(
        "-bl",
        "--by-loss-status",
        dest="loss_status",
        help="Include certain orthology classes (FI, I, PI, UL, M, PM, L, UL)",
        required=False,
        type=str,
    )
    app.add_argument(
        "-bc",
        "--by-orthology-class",
        dest="orthology_class",
        help="Include certain orthology relationships (o2o, o2m, m2m, m2m, o2z)",
        required=False,
        type=str,
    )
    app.add_argument(
        "-bs",
        "--by-orthology-score",
        dest="orthology_score",
        help="Preserve orthology scores greater or equal to a threshold (0.0 - 1.0)",
        required=False,
        type=float,
    )
    app.add_argument(
        "-to",
        "--to",
        help="Specify the conversion format for .bed file (gtf, gff3)",
        required=False,
        type=str,
        choices=["gtf", "gff", "bed"],
        default="gtf",
    )
    app.add_argument(
        "-tg",
        "--target",
        dest="bed_type",
        help="Specify the .bed input file to used by the program",
        type=str,
        choices=["bed", "utr"],
        default="utr",
    )
    app.add_argument(
        "-bp",
        "--by-paralog-score",
        dest="min_paralog_score",
        help="Preserve transcripts with paralog projection probabilities",
        required=False,
        type=float,
    )
    app.add_argument(
        "-w",
        "--with-isoforms",
        help="Path to a custom isoform table (default: None)",
        required=False,
        default=None,
        type=str,
    )
    app.add_argument(
        "-ext",
        "--extract",
        help="Flag or option to extract sequences from the filtered projections. "
        "Can be 'query', 'reference', or just set as a flag (default: False).",
        required=False,
        default=False,
        nargs="?",
        choices=["query", "reference"],
    )
    app.add_argument(
        "--outdir",
        "-o",
        help="Path to posTOGA output directory",
        required=False,
        type=str,
        default=None,
    )
    app.add_argument(
        "--only-table",
        "-ot",
        help="Only produce the toga.table file",
        required=False,
        action="store_true",
    )
    app.add_argument(
        "--only-convert",
        "-oc",
        help="Only convert the toga.table file to gtf/gff",
        required=False,
        action="store_true",
    )
    app.add_argument(
        "-L",
        "--level",
        dest="log_level",
        help="Logging verbosity (debug, info, warn, off)",
        required=False,
        default="info",
        type=lambda value: value.lower(),
        choices=["debug", "info", "warn", "off"],
    )
    app.add_argument(
        "--depure",
        "-d",
        help="Remove any trace of other postoga runs/files",
        required=False,
        action="store_true",
    )
    app.add_argument(
        "--version",
        "-v",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    args = app.parse_args()

    return args

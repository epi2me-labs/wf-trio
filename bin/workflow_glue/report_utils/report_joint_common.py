"""Common functions to joint report."""
import re

from dominate.tags import li, p, ul
from dominate.util import raw
from ezcharts.layout.snippets.table import DataTable
import pandas as pd


def load_rtg(rtg_file, sample):
    """Read in RTG summary file and output to a table."""
    percs_conc_m = None
    percs_conc_f = None
    percs_conc_both = None
    total_variant = None
    total_records = None
    perc_variant = None
    ploidy = None
    total_conflicts = None
    perc_conflict = None
    incomplete = 0

    with open(rtg_file, "r") as rtg_report:
        for line in rtg_report:
            # Parse the report by iterating and regex particular lines.
            if "Concordance" in line:
                percs = re.findall(r"\d+\.\d+\%", line)
                if percs:
                    percs_conc_m = percs[0]
                    percs_conc_f = percs[1]
                    percs_conc_both = percs[2]
            elif "expected call ploidy" in line:
                total_records = re.search(r"(\d+)\/(\d+)", line)[2]
                ploidy = re.search(r"\d+\.\d+\%", line)[0]
            elif "violation of Mendelian" in line:
                totals = re.search(r"(\d+)\/(\d+)", line)
                total_variant = totals[2]
                total_conflicts = totals[1]
                perc_conflict = re.search(r"\d+\.\d+\%", line)[0]
            elif "due to incomplete calls" in line:
                incomplete = re.search(r"\d+\.\d+\%", line)[0]
        pv = format(
            (100.0 * int(total_variant) / int(total_records)), '.2f')
        perc_variant = "{}%".format(str(pv))

    rtg_table = pd.DataFrame.from_dict(
        {
            'Family': sample,
            'Total records': [total_records],
            'Total variants': [total_variant],
            'Perc variant': [perc_variant],
            'Concordance M': [percs_conc_m],
            'Concordance F': [percs_conc_f],
            'Concordance M+F': [percs_conc_both],
            'Unexpected ploidy': [ploidy],
            'Incomplete calls': [incomplete],
            'Total conflicts': [total_conflicts],
            'Perc conflict': [perc_conflict]
        }, orient='columns')
    return rtg_table


def add_pedigree_section(ped_file, report):
    """Add pedigree TSV as a table to report for reference."""
    with report.add_section('Pedigree summary', 'Pedigree'):
        colnames = [
            "Family ID",
            "Individual ID",
            "Paternal ID",
            "Maternal ID",
            "Sex",
            "Phenotype"]
        p("""
This table shows the pedigree TSV file that was input to the workflow
""")
        ul(
            li("Sex - (1=male; 2=female; other=unknown)"),
            li("Phenotype - (-9=missing;0=missing;1=unaffected;2=affected)."),
        )
        raw("""Find the Pedigree file format definition
        <a href="https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format">here</a>.""")  # noqa: E501
        df = pd.read_csv(
            ped_file, sep='\t', names=colnames, header=None,
            dtype={
                "Family ID": "string",
                "Individual ID": "string",
                "Paternal ID": "string",
                "Maternal ID": "string",
                "Sex": "string",
                "Phenotype": "string"}
            )
        DataTable.from_pandas(df, use_index=False)


def add_rtg_mendelian(rtg_file, report, sample_name):
    """Add RTG Mendelian section to report."""
    with report.add_section('Mendelian conflict summary', 'Mendelian'):
        p("""
This table provides stats from the RTG Mendelian Report
""")
        raw("""Find more information about the RTG Mendelian tool
        <a href="https://realtimegenomics.github.io/rtg-tools/rtg_command_reference.html#mendelian">here</a>.""")  # noqa: E501
        p("""
These are the definitions of the columns you will find in the table.
""")
        ul(
            li("Family: Family ID."),
            li("Total records: Total number of records in the VCF."),
            li("Total variants: Total number of variants identified \
               in the multi-sample VCF that were variant in \
               at least 1 family member."),
            li("Percent variant: Percentage of the VCF that \
               is variant in at least 1 family member."),
            li("Concordance M: Percentage of variants that \
               are in concordance with the Male."),
            li("Concordance F: Percentage of variants that are \
               in concordance with the Female."),
            li("Concordance M+F:  Percentage of variants that are in \
               concordance with both the Male and Female."),
            li("Unexpected ploidy: Percentage of records that \
               did not conform to expected ploidy."),
            li("Incomplete calls: Percentage of records that \
               had indeterminate consistency status due to incomplete calls."),
            li("Total conflicts: Number of records that \
               contained a violation of Mendelian constraints."),
            li("Percentage conflict: Percentage of records that \
            contained a violation of Mendelian constraints.")
        )

        rtg_table = load_rtg(rtg_file, sample_name)
        DataTable.from_pandas(rtg_table, use_index=False)

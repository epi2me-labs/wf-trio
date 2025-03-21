"""Create joint SV workflow report."""
import os

from ezcharts.components.reports import labs

from .report_utils.report_joint_common \
    import add_pedigree_section, add_rtg_mendelian  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("jointSVReport")
    report = labs.LabsReport(
        "Trio structural variation report", "wf-trio",
        args.params, args.versions, args.wf_version)
    if args.ped_file and os.path.exists(args.ped_file):
        add_pedigree_section(args.ped_file, report)
    if args.rtg_mendelian and os.path.exists(args.rtg_mendelian):
        add_rtg_mendelian(args.rtg_mendelian, report, args.sample_name)
    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("multisample_sv_report")
    parser.add_argument("report", help="Multi sample SV report output file")
    parser.add_argument(
        "--sample_name", default='Sample',
        help="Sample name"
    )
    parser.add_argument(
        "--versions", default=None, required=True,
        help="Directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--wf_version", default='unknown',
        help="Version of the executed workflow")
    parser.add_argument(
        "--ped_file", required=False,
        help="Pedigree TSV file")
    parser.add_argument(
        "--rtg_mendelian", required=False,
        help="RTG mendelian report")
    return parser

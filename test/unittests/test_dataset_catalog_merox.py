import csv
import hashlib
import io
import json
import unittest
import zipfile
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

from click.testing import CliRunner
from claudio.dataset_catalog.xifdr import read_merox_csms
from claudio.dataset_catalog.cli import main


FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "merox_pxd010895"
)


class TestMeroXStandardize(unittest.TestCase):
    def arguments(self, output):
        return [
            "standardize",
            "--results", str(FIXTURE / "velos.zhrm"),
            "--fasta", str(FIXTURE / "search.fasta"),
            "--dataset-accession", "PXD010895",
            "--output", str(output),
        ]

    def test_audit_retains_source_rows_terminal_types_and_spectrum_references(self):
        fixture_manifest = json.loads(
            (FIXTURE / "manifest.json").read_bytes()
        )
        before = {}
        for name, metadata in fixture_manifest["files"].items():
            body = (FIXTURE / name).read_bytes()
            self.assertEqual(
                hashlib.sha256(body).hexdigest(),
                metadata["sha256"],
                name,
            )
            before[name] = body

        with zipfile.ZipFile(FIXTURE / "velos.zhrm") as archive:
            source_rows = [
                row
                for row in csv.reader(
                    io.StringIO(
                        archive.read("Result.csv").decode("utf-8-sig")
                    ),
                    delimiter=";",
                )
                if any(value.strip() for value in row)
            ]

        with TemporaryDirectory() as directory:
            output = Path(directory) / "audit"
            result = CliRunner().invoke(
                main, self.arguments(output) + ["--audit-only"]
            )
            self.assertEqual(
                result.exit_code, 0,
                result.output + repr(result.exception),
            )

            report = json.loads(
                (output / "import_report.json").read_bytes()
            )
            manifest = json.loads(
                (output / "manifest.json").read_bytes()
            )
            records = report["records"]

            self.assertEqual(report["format"], "merox_csm")
            self.assertEqual(report["parser"]["name"], "pyXLMS")
            self.assertEqual(report["parser"]["parsed_records"], 112)
            self.assertEqual(len(records), 112)
            self.assertEqual(
                [record["source_record_number"] for record in records],
                list(range(1, 113)),
            )
            self.assertEqual(
                len({record["source_record_id"] for record in records}),
                112,
            )
            self.assertEqual(
                [record["source_values"] for record in records],
                source_rows,
            )

            spectra = [record["spectrum"] for record in records]
            self.assertEqual(
                len({spectrum["uuid"] for spectrum in spectra}), 103
            )
            for spectrum in spectra:
                self.assertEqual(
                    spectrum["reference_kind"], "merox_embedded_member"
                )
                self.assertIsNone(spectrum["original_file"])

            kinds = Counter(
                endpoint["kind"]
                for record in records
                for endpoint in record["endpoints"].values()
            )
            self.assertEqual(dict(kinds), {
                "residue": 212,
                "protein_n_terminus": 10,
                "protein_c_terminus": 2,
            })
            self.assertEqual(
                records[5]["endpoints"]["a"]["protein_position"], 1
            )
            self.assertEqual(
                records[9]["endpoints"]["b"]["protein_position"], 585
            )

            multiple = sum(
                len([
                    part
                    for part in record["source_fields"][
                        "All linkage positions"
                    ].split("#")
                    if part
                ]) > 1
                for record in records
            )
            self.assertEqual(multiple, 104)

            for record in records:
                self.assertEqual(record["decoy_status"], "unknown")
                self.assertEqual(
                    record["fdr"]["psm"]["status"], "missing"
                )

            self.assertEqual(manifest["counts"]["source_records"], 112)
            self.assertEqual(manifest["counts"]["parser_errors"], 0)
            self.assertFalse(manifest["input_exported"])
            self.assertFalse(manifest["assessment_performed"])
            self.assertFalse(manifest["fdr_filter_applied"])
            self.assertFalse((output / "claudio_input.csv").exists())

        for name, body in before.items():
            self.assertEqual((FIXTURE / name).read_bytes(), body)

    def test_export_requires_merox_policies_to_be_resolved(self):
        with TemporaryDirectory() as directory:
            output = Path(directory) / "input"
            result = CliRunner().invoke(main, self.arguments(output))

            self.assertNotEqual(result.exit_code, 0)
            self.assertIn("--audit-only", result.output)
            self.assertFalse(output.exists())
    def test_missing_spectrum_keeps_its_source_rows_in_the_audit(self):
        with TemporaryDirectory() as directory:
            damaged = Path(directory) / "missing_spectrum.zhrm"

            with zipfile.ZipFile(FIXTURE / "velos.zhrm") as original:
                rows = [
                    row
                    for row in csv.reader(
                        io.StringIO(
                            original.read("Result.csv").decode("utf-8-sig")
                        ),
                        delimiter=";",
                    )
                    if any(value.strip() for value in row)
                ]
                missing_uuid = rows[0][23]
                missing_member = f"spectra/{missing_uuid}.txt"

                with zipfile.ZipFile(damaged, "w") as output:
                    for member in original.infolist():
                        if member.filename != missing_member:
                            output.writestr(
                                member, original.read(member.filename)
                            )

            report = read_merox_csms(
                damaged,
                FIXTURE / "search.fasta",
                dataset_accession="PXD010895",
            )

        records = report["records"]
        affected = [
            record for record in records
            if record["source_fields"]["Spectrum UUID"] == missing_uuid
        ]

        self.assertEqual(len(records), 112)
        self.assertEqual(
            [record["source_values"] for record in records], rows
        )
        self.assertTrue(affected)
        for record in affected:
            self.assertEqual(
                record["record_reason"], "merox_source_unresolved"
            )
            self.assertIsNone(record["spectrum"])

        self.assertEqual(
            report["parser"]["parsed_records"],
            112 - len(affected),
        )

    def test_unrelated_zip_is_not_read_as_a_merox_container(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "other.zip"
            output = root / "audit"

            with zipfile.ZipFile(archive, "w") as handle:
                handle.writestr("readme.txt", "Unrelated archive")

            arguments = self.arguments(output) + ["--audit-only"]
            arguments[arguments.index("--results") + 1] = str(archive)
            result = CliRunner().invoke(main, arguments)

            self.assertNotEqual(result.exit_code, 0)
            self.assertIn(
                "not a supported native MeroX", result.output
            )
            self.assertFalse(output.exists())

if __name__ == "__main__":
    unittest.main()
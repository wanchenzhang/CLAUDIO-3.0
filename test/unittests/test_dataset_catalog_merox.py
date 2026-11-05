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
class TestMeroXLocalization(unittest.TestCase):
    def read_bundle(self):
        from claudio.dataset_catalog.xifdr import read_merox_csms

        return read_merox_csms(
            FIXTURE / "velos.zhrm",
            FIXTURE / "search.fasta",
            dataset_accession="PXD010895",
        )

    def localization_reader(self):
        from claudio.dataset_catalog import xifdr

        reader = getattr(xifdr, "_merox_localization", None)
        self.assertTrue(
            callable(reader),
            "Missing paired MeroX localization reader",
        )
        return reader

    def test_native_candidates_keep_pairs_order_and_best_membership(self):
        bundle = self.read_bundle()
        total_candidates = 0

        for record in bundle["records"]:
            self.assertIn("localization", record)
            localization = record["localization"]
            fields = record["source_fields"]
            tokens = fields["All linkage positions"].split("#")
            candidates = localization["candidates"]

            self.assertEqual(
                localization["raw"],
                fields["All linkage positions"],
            )
            self.assertEqual(
                [item["raw"] for item in candidates],
                tokens,
            )
            self.assertEqual(
                [item["source_candidate_number"] for item in candidates],
                list(range(1, len(tokens) + 1)),
            )
            self.assertEqual(
                [item["source_fields"] for item in candidates],
                [token.split("~") for token in tokens],
            )

            expected_best = [
                number
                for number, token in enumerate(tokens, 1)
                if token.split("~")[1:3] == [
                    fields["best linkage position peptide 1"],
                    fields["best linkage position peptide 2"],
                ]
            ]
            self.assertEqual(
                localization["source_best_candidate_numbers"],
                expected_best,
            )
            self.assertTrue(expected_best)

            for candidate in candidates:
                self.assertEqual(
                    candidate["status"], "valid", candidate
                )
                self.assertEqual(
                    candidate["is_source_best"],
                    candidate["source_candidate_number"] in expected_best,
                )

                for site in ("a", "b"):
                    endpoint = candidate["endpoints"][site]
                    sequence = bundle["proteins"][
                        endpoint["source_fasta_id"]
                    ]["sequence"]
                    start = endpoint["peptide_start"]
                    end = endpoint["peptide_end"]

                    self.assertEqual(
                        sequence[start - 1:end],
                        endpoint["peptide"],
                    )
                    self.assertEqual(
                        endpoint["protein_position"],
                        start + endpoint["link_position"] - 1,
                    )

            self.assertEqual(record["record_status"], "not_evaluable")
            total_candidates += len(candidates)

        self.assertEqual(len(bundle["records"]), 112)
        self.assertEqual(total_candidates, 554)

        first = bundle["records"][0]["localization"]["candidates"]
        self.assertEqual(
            [
                (
                    item["endpoints"]["a"]["link_position"],
                    item["endpoints"]["b"]["link_position"],
                )
                for item in first
            ],
            [(4, 9), (7, 9), (10, 8)],
        )

        n_terminal = bundle["records"][5]["localization"]["candidates"]
        self.assertTrue(any(
            item["endpoints"]["a"]["kind"] == "protein_n_terminus"
            and item["endpoints"]["a"]["protein_position"] == 1
            for item in n_terminal
        ))

        c_terminal = bundle["records"][9]["localization"]["candidates"]
        self.assertTrue(any(
            item["endpoints"]["b"]["kind"] == "protein_c_terminus"
            and item["endpoints"]["b"]["protein_position"] == 585
            for item in c_terminal
        ))

    def test_invalid_candidate_is_retained_without_losing_other_pairs(self):
        reader = self.localization_reader()
        bundle = self.read_bundle()
        fields = dict(bundle["records"][0]["source_fields"])

        fields["All linkage positions"] = (
            "80~E4~K9~0~47#64~E9999~K9~0~47"
        )
        original = dict(fields)

        localization = reader(fields, bundle["proteins"])

        self.assertEqual(fields, original)
        self.assertEqual(localization["status"], "unresolved")
        self.assertEqual(len(localization["candidates"]), 2)
        self.assertEqual(
            localization["candidates"][0]["status"], "valid"
        )

        invalid = localization["candidates"][1]
        self.assertEqual(invalid["status"], "invalid")
        self.assertEqual(invalid["raw"], "64~E9999~K9~0~47")
        self.assertTrue(invalid["reason"])
        self.assertEqual(
            localization["source_best_candidate_numbers"], [1]
        )

    def test_missing_source_best_does_not_select_another_candidate(self):
        reader = self.localization_reader()
        bundle = self.read_bundle()
        fields = dict(bundle["records"][0]["source_fields"])

        fields["All linkage positions"] = "13~K10~D8~1~4"

        localization = reader(fields, bundle["proteins"])

        self.assertEqual(localization["status"], "unresolved")
        self.assertEqual(
            localization["source_best_candidate_numbers"], []
        )
        self.assertEqual(len(localization["candidates"]), 1)
        self.assertEqual(
            localization["candidates"][0]["status"], "valid"
        )
        self.assertFalse(
            localization["candidates"][0]["is_source_best"]
        )
    def test_empty_or_malformed_candidates_remain_unresolved(self):
        reader = self.localization_reader()
        bundle = self.read_bundle()
        original = bundle["records"][0]["source_fields"]

        cases = [
            ("", []),
            ("80~E4~K9", ["80~E4~K9"]),
            (
                "80~E4~K9~0~47##64~E7~K9~0~47",
                ["80~E4~K9~0~47", "", "64~E7~K9~0~47"],
            ),
        ]

        for raw, expected_tokens in cases:
            with self.subTest(raw=raw):
                fields = dict(original)
                fields["All linkage positions"] = raw
                before = dict(fields)

                result = reader(fields, bundle["proteins"])

                self.assertEqual(fields, before)
                self.assertEqual(result["raw"], raw)
                self.assertEqual(result["status"], "unresolved")
                self.assertTrue(result["reason"])
                self.assertEqual(
                    [item["raw"] for item in result["candidates"]],
                    expected_tokens,
                )
                if expected_tokens:
                    self.assertTrue(any(
                        item["status"] == "invalid"
                        for item in result["candidates"]
                    ))

    def test_duplicate_pairs_keep_distinct_source_entries(self):
        reader = self.localization_reader()
        bundle = self.read_bundle()
        fields = dict(bundle["records"][0]["source_fields"])
        tokens = ["80~E4~K9~0~47", "64~E4~K9~0~33"]
        fields["All linkage positions"] = "#".join(tokens)

        result = reader(fields, bundle["proteins"])
        candidates = result["candidates"]

        self.assertEqual(len(candidates), 2)
        self.assertEqual(
            [item["raw"] for item in candidates], tokens
        )
        self.assertEqual(
            [item["source_candidate_number"] for item in candidates],
            [1, 2],
        )
        self.assertEqual(
            result["source_best_candidate_numbers"], [1, 2]
        )
        self.assertTrue(all(
            item["status"] == "valid" and item["is_source_best"]
            for item in candidates
        ))
        self.assertEqual(
            candidates[0]["endpoints"],
            candidates[1]["endpoints"],
        )
class TestMeroXNativeConfidence(unittest.TestCase):
    def test_saved_filter_evidence_does_not_promote_records(self):
        from claudio.dataset_catalog.xifdr import read_merox_csms

        report = read_merox_csms(
            FIXTURE / "velos.zhrm",
            FIXTURE / "search.fasta",
            dataset_accession="PXD010895",
        )

        self.assertEqual(len(report["records"]), 112)

        for record in report["records"]:
            self.assertIn("native_confidence", record)
            evidence = record["native_confidence"]

            self.assertEqual(
                evidence["status"], "within_saved_cutoff"
            )
            self.assertEqual(
                evidence["fdr_level"], "native_unspecified"
            )
            self.assertEqual(evidence["cutoff"], 0.05)
            self.assertEqual(
                evidence["value"],
                record["native_local_fdr"]["value"],
            )
            self.assertFalse(evidence["application_verified"])

            self.assertEqual(record["decoy_status"], "unknown")
            self.assertEqual(record["fdr"]["psm"]["status"], "missing")
            self.assertEqual(record["record_status"], "not_evaluable")
            self.assertEqual(
                record["source_fields"]["is Selected in Table"],
                "false",
            )

    def test_unusable_settings_or_values_do_not_pass_the_saved_filter(self):
        from claudio.dataset_catalog import xifdr

        reader = getattr(xifdr, "_merox_confidence_evidence", None)
        self.assertTrue(
            callable(reader),
            "Missing native MeroX confidence evidence reader",
        )

        settings = {
            "FDRCUTOFF": "0.05",
            "APPLYFDR": "1",
            "PERFORMDECOY": "1",
        }
        valid = {"raw": "0.05", "value": 0.05, "status": "valid"}

        cases = [
            ({}, valid, "within_saved_cutoff"),
            (
                {},
                {"raw": "0.06", "value": 0.06, "status": "valid"},
                "above_saved_cutoff",
            ),
            ({"APPLYFDR": "0"}, valid, "filter_disabled"),
            ({"APPLYFDR": "unknown"}, valid, "unresolved"),
            ({"PERFORMDECOY": "0"}, valid, "unresolved"),
            ({"FDRCUTOFF": ""}, valid, "unresolved"),
            ({"FDRCUTOFF": "nan"}, valid, "unresolved"),
            ({"FDRCUTOFF": "1.5"}, valid, "unresolved"),
            (
                {},
                {"raw": "", "value": None, "status": "missing"},
                "unresolved",
            ),
        ]

        for overrides, local_fdr, expected in cases:
            with self.subTest(overrides=overrides, local_fdr=local_fdr):
                supplied = dict(settings)
                supplied.update(overrides)
                before_settings = dict(supplied)
                before_fdr = dict(local_fdr)

                evidence = reader(supplied, local_fdr)

                self.assertEqual(evidence["status"], expected)
                self.assertFalse(evidence["application_verified"])
                self.assertEqual(
                    evidence["fdr_level"], "native_unspecified"
                )
                self.assertEqual(supplied, before_settings)
                self.assertEqual(local_fdr, before_fdr)
if __name__ == "__main__":
    unittest.main()
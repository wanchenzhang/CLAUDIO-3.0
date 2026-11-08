import csv
import tempfile
import unittest
from pathlib import Path

from claudio.dataset_catalog import cli, maxlynx


FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures" / "maxlynx_pxd053456"
)


class TestMaxLynxNative57(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        with (FIXTURE / "DSBU_crosslinkMsms.txt").open(
            newline="", encoding="utf-8-sig"
        ) as handle:
            self.columns = next(csv.reader(handle, delimiter="\t"))

        self.fasta = self.root / "synthetic.fasta"
        self.fasta.write_text(
            ">protein_a\nAAKAAA\n>protein_b\nAKAAAA\n",
            encoding="utf-8",
        )
        self.parameters = self.root / "synthetic_parameters.xml"
        self.parameters.write_text(
            "<MaxQuantParams>"
            "<crossLinker>TEST_LINKER</crossLinker>"
            "<identifierParseRule>^&gt;([^ ]+)</identifierParseRule>"
            "<psmFdrCrosslink>0.01</psmFdrCrosslink>"
            "</MaxQuantParams>",
            encoding="utf-8",
        )
        self.results = self.root / "synthetic.txt"

    def paired_row(self):
        row = dict.fromkeys(self.columns, "")
        row.update({
            "Raw file": "synthetic_run",
            "Scan number": "100",
            "Charge": "3",
            "Sequence1": "AAKAAA",
            "Sequence2": "AKAAAA",
            "Modified sequence1": "_AAKAAA_",
            "Modified sequence2": "_AKAAAA_",
            "Length1": "6",
            "Length2": "6",
            "Proteins1": "protein_a",
            "Proteins2": "protein_b",
            "Decoy": "forward",
            "Score": "20",
            "Partial score 1": "9",
            "Partial score 2": "11",
            "Delta score": "4",
            "PEP": "0.0001",
            "Crosslink product type": "Inter-protein link",
            "Pair of linked residues": "K3-K2",
            "Protein index of Crosslink 1": "3",
            "Peptide index of Crosslink 1": "3",
            "AA name of Crosslink 1": "K",
            "Protein index of Crosslink 2": "2",
            "Peptide index of Crosslink 2": "2",
            "AA name of Crosslink 2": "K",
        })
        return row

    def write_rows(self, rows, columns=None):
        with self.results.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns or self.columns, delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)

    def read(self):
        try:
            return maxlynx.read_maxlynx_csms(
                self.results, self.fasta,
                dataset_accession="PXD000001",
                parameters_path=self.parameters,
                backend="native",
            )
        except ValueError as error:
            self.fail(f"Reviewed native column profile is unsupported: {error}")

    def test_detects_the_actual_native_table(self):
        try:
            detected = cli._detect_psm_format(
                FIXTURE / "DSBU_crosslinkMsms.txt"
            )
        except ValueError as error:
            self.fail(f"Native MaxLynx layout is not detected: {error}")
        self.assertEqual(detected, "maxlynx_csm")

    def test_reads_geometry_without_rewriting_native_fields(self):
        paired = self.paired_row()
        rows = [paired]
        for product in (
            "Linear", "Mono-linked", "Loop-linked", "Mono- and loop-linked"
        ):
            rows.append(dict(paired, **{
                "Crosslink product type": product,
                "Sequence2": "-",
            }))
        rows.append(dict(paired, Decoy="decoy"))
        self.write_rows(rows)
        bundle = self.read()
        records = bundle["records"]

        self.assertEqual(len(records), 6)
        self.assertEqual(bundle["source_columns"], self.columns)
        self.assertEqual(
            [record["source_fields"] for record in records], rows
        )
        self.assertEqual(
            [record["source_values"] for record in records],
            [[row[column] for column in self.columns] for row in rows],
        )
        self.assertEqual(records[0]["record_status"], "valid")
        self.assertEqual(records[0]["record_kind"], "paired_crosslink_candidate")
        self.assertEqual(records[0]["endpoints"]["a"]["peptide_start"], 1)
        self.assertEqual(records[0]["endpoints"]["a"]["protein_position"], 3)
        self.assertEqual(records[0]["endpoints"]["a"]["link_position"], 3)
        self.assertEqual(records[0]["endpoints"]["b"]["peptide_start"], 1)
        self.assertEqual(records[0]["endpoints"]["b"]["protein_position"], 2)
        self.assertEqual(records[0]["endpoints"]["b"]["link_position"], 2)
        self.assertEqual(
            [record["record_kind"] for record in records[1:5]],
            ["linear_peptide", "dead_end", "loop_link", "mono_and_loop_link"],
        )
        self.assertTrue(all(
            record["record_status"] == "not_evaluable"
            for record in records[1:]
        ))
        self.assertEqual(records[-1]["decoy_status"], "source_decoy")
        self.assertEqual(records[-1]["record_reason"], "decoy_record")
        self.assertEqual(records[0]["source_pep"]["raw"], "0.0001")
        self.assertEqual(
            bundle["fdr_policy"]["verification_status"], "configured_only"
        )
        self.assertFalse(bundle["fdr_policy"]["per_record_fdr_applied"])

    def test_rejects_a_site_that_differs_from_the_peptide(self):
        row = self.paired_row()
        row["AA name of Crosslink 1"] = "A"
        row["Pair of linked residues"] = "A3-K2"
        self.write_rows([row])
        record = self.read()["records"][0]
        self.assertEqual(record["record_status"], "not_evaluable")
        self.assertEqual(record["source_fields"], row)

    def test_preserves_unknown_product_types_without_promoting_them(self):
        row = self.paired_row()
        row["Crosslink product type"] = "Unreviewed product"
        self.write_rows([row])
        record = self.read()["records"][0]
        self.assertEqual(record["record_kind"], "unsupported_product_type")
        self.assertEqual(record["record_status"], "not_evaluable")
        self.assertEqual(record["source_fields"], row)

    def test_conflicting_product_headers_are_not_detected_as_maxlynx(self):
        columns = self.columns + ["Crosslink Product Type"]
        row = self.paired_row()
        row["Crosslink Product Type"] = "SingleProtein"
        self.write_rows([row], columns)
        with self.assertRaises(ValueError):
            cli._detect_psm_format(self.results)

    def test_changed_native_source_is_rejected_after_import(self):
        self.write_rows([self.paired_row()])
        bundle = self.read()
        with self.results.open("ab") as handle:
            handle.write(b"\n")
        with self.assertRaisesRegex(ValueError, "(?i)checksum"):
            maxlynx.verify_maxlynx_sources(bundle)


if __name__ == "__main__":
    unittest.main()

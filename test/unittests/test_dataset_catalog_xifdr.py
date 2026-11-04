import csv
import hashlib
import json
import unittest
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, Mock, call, patch

import requests
from click.testing import CliRunner
from pyXLMS import parser as pyxlms_parser

from claudio.dataset_catalog.xifdr import (
    read_xifdr_psms, read_search_fasta, normalize_xifdr_peptide, target_proteins,
    map_reference_proteins, resolve_uniprot_proteins, prepare_af_input,
    write_xifdr_outputs, assess_af_complexes, compare_protein_sequences, search_uniprot_candidates,
    project_xifdr_candidates, write_projection_outputs, write_claudio_input,
)
from claudio.dataset_catalog.cli import main as catalog_main
from claudio.dataset_catalog import xifdr


COLUMNS = [
    "SearchID", "PSMID", "run", "scan", "Protein1", "Description1", "Decoy1",
    "Protein2", "Description2", "Decoy2", "PepSeq1", "PepSeq2", "PepPos1", "PepPos2",
    "PeptideLength1", "PeptideLength2", "LinkPos1", "LinkPos2", "ProteinLinkPos1",
    "ProteinLinkPos2", "Crosslinker", "Score", "isDecoy", "isTT", "isTD", "isDD",
    "fdrGroup", "fdr", "", "PeptidePairFDR", "Protein1FDR", "Protein2FDR", "LinkFDR", "PPIFDR",
]
FASTA = ">customA\nMACKRTAK\n>customB\nGGKPEPK\n"
REFERENCES = ">sp|P12345|A OX=9031\nMACKRTAK\n>sp|Q99999|B OX=9031\nGGKPEPK\n"
CLAUDIO_COLUMNS = [
    "peptide1", "peptide2", "position1", "position2",
    "k_pos1", "k_pos2", "entry1", "entry2",
]


def psm(**changes):
    row = dict.fromkeys(COLUMNS, "")
    row.update(SearchID="search-1", PSMID="psm-1", run="replicate-1", scan="42",
               Protein1="customA", Protein2="customB", Decoy1="false", Decoy2="false",
               PepSeq1="ACcmK", PepSeq2="KPEP", PepPos1="2", PepPos2="3",
               PeptideLength1="3", PeptideLength2="4", LinkPos1="3", LinkPos2="1",
               ProteinLinkPos1="4", ProteinLinkPos2="3", Crosslinker="BS3", Score="12.5",
               isDecoy="false", isTT="true", isTD="false", isDD="false",
               fdrGroup="between > 4 ", fdr="0.02", PeptidePairFDR="0.01",
               Protein1FDR="0", Protein2FDR="0", LinkFDR="0.005", PPIFDR="0.125")
    row.update(changes)
    return row


def response(url, payload):
    value = requests.Response()
    value.status_code = 200
    value.url = url
    value.headers["X-UniProt-Release"] = "example-release"
    value._content = (json.dumps(payload) if isinstance(payload, dict) else payload).encode()
    return value


class TestXiFDRPSM(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / "search.fasta"
        self.fasta.write_text(FASTA)
        self.refs = self.root / "uniprot.fasta"
        self.refs.write_text(REFERENCES)

    def read(self, rows=None, columns=None):
        path = self.root / "PSM.csv"
        columns = COLUMNS if columns is None else columns
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(columns)
            for row in [psm()] if rows is None else rows:
                writer.writerow([row.get(k, "") for k in columns] if isinstance(row, dict) else row)
        return read_xifdr_psms(path, self.fasta, dataset_accession="PXD000001")

    def mapping(self, bundle):
        return map_reference_proteins(bundle["proteins"], read_search_fasta(self.refs), taxon_id=9031)

    def export_pipeline(self, bundle, mappings=None, **options):
        return write_claudio_input(
            bundle, self.mapping(bundle) if mappings is None else mappings,
            self.root / "pipeline", fdr_level="link", fdr_cutoff=0.01,
            allowed_crosslinkers={"BS3"}, allowed_residues={"K"}, **options,
        )






    def test_pipeline_input_leaves_tracking_ids_for_m1_to_assign(self):
        bundle = self.read()
        output = self.export_pipeline(bundle)
        with (output / "claudio_input.csv").open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            columns, rows = reader.fieldnames, list(reader)
        self.assertEqual(columns[:8], CLAUDIO_COLUMNS)
        self.assertFalse({"source_record_id", "original_pos_a", "original_pos_b", "candidate_id"} & set(columns))
        self.assertFalse({"seq_a", "seq_b", "position_status_a", "position_status_b"} & set(columns))
        self.assertEqual(rows[0]["xlms_source_record_id"], bundle["records"][0]["source_record_id"])
        self.assertEqual(rows[0]["peptide1"], "ACK")
        self.assertEqual(rows[0]["k_pos1"], "3")
        self.assertEqual(rows[0]["entry1"], "P12345")
        self.assertEqual(rows[0]["source_crosslinker"], "BS3")
        self.assertEqual(json.loads(rows[0]["source_fields_json"])["PepSeq1"], "ACcmK")
        self.assertFalse((output / "af_input.csv").exists())

    def test_pipeline_export_keeps_repeated_psms_and_experimental_runs(self):
        bundle = self.read([psm(), psm(), psm(run="replicate-2")])
        output = self.export_pipeline(bundle)
        with (output / "claudio_input.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 3)
        self.assertEqual([r["source_run"] for r in rows], ["replicate-1", "replicate-1", "replicate-2"])
        self.assertEqual(len({r["xlms_source_record_id"] for r in rows}), 3)

    def test_pipeline_csv_is_accepted_by_the_existing_m1_reader(self):
        from claudio.module01.src.io.read_in import read_inputfile
        from claudio.module01.src.main import main as m1_main

        bundle = self.read([psm(), psm(run="replicate-2")])
        output = self.export_pipeline(bundle)
        columns = ["pep_a", "pep_b", "pos_a", "pos_b", "res_pos_a", "res_pos_b", "unip_id_a", "unip_id_b"]
        default_headers = next(option.default for option in m1_main.params
                               if option.name == "projections").split(",")
        with (output / "claudio_input.csv").open(newline="", encoding="utf-8") as handle:
            self.assertEqual(csv.DictReader(handle).fieldnames[:8], default_headers)
        data = read_inputfile(str(output / "claudio_input.csv"), dict(zip(default_headers, columns)))
        self.assertEqual(len(data), 2)
        self.assertEqual(data["xlms_source_record_id"].tolist(), [row["source_record_id"] for row in bundle["records"]])
        self.assertEqual(data["source_record_id"].tolist(), ["record_000001", "record_000002"])
        self.assertEqual(data["source_run"].tolist(), ["replicate-1", "replicate-2"])
        self.assertEqual(data["pos_a"].tolist(), [4, 4])
        self.assertEqual(data["pos_b"].tolist(), [3, 3])
        self.assertEqual(data["pep_a"].tolist(), ["ACK", "ACK"])
        self.assertEqual(data["pep_b"].tolist(), ["KPEP", "KPEP"])
        self.assertEqual(data["res_pos_a"].tolist(), [3.0, 3.0])
        self.assertEqual(data["res_pos_b"].tolist(), [1.0, 1.0])
        self.assertEqual(data["unip_id_a"].tolist(), ["P12345", "P12345"])
        self.assertEqual(data["unip_id_b"].tolist(), ["Q99999", "Q99999"])

    def test_pipeline_scope_does_not_discard_unselected_records_from_the_audit(self):
        selected = psm(Protein2="customA", PepSeq2="ACK", PepPos2="2", PeptideLength2="3",
                       LinkPos2="3", ProteinLinkPos2="4")
        bundle = self.read([selected, psm()])
        output = self.export_pipeline(bundle, source_protein_ids={"customA"})
        audit = json.loads((output / "psm_audit.json").read_bytes())
        self.assertEqual(len(audit), 2)
        self.assertEqual(audit[0]["input_status"], "ready")
        self.assertIn("outside_selected_source_proteins", audit[1]["input_reasons"])
        manifest = json.loads((output / "manifest.json").read_bytes())
        self.assertEqual(manifest["counts"]["pipeline_input_rows"], 1)

    def test_pipeline_residue_filter_keeps_the_original_coordinates_and_reason(self):
        bundle = self.read([psm(), psm(LinkPos1="1", ProteinLinkPos1="2")])
        output = self.export_pipeline(bundle)
        audit = json.loads((output / "psm_audit.json").read_bytes())
        self.assertIn("residue_not_selected_a", audit[1]["input_reasons"])
        self.assertEqual(audit[1]["endpoints"]["a"]["protein_position"], 2)
        self.assertEqual(audit[1]["source_fields"]["LinkPos1"], "1")

    def test_pipeline_crosslinker_filter_does_not_relabel_unknown_chemistry(self):
        bundle = self.read([psm(), psm(Crosslinker="other chemistry")])
        output = self.export_pipeline(bundle)
        audit = json.loads((output / "psm_audit.json").read_bytes())
        self.assertIn("crosslinker_not_selected", audit[1]["input_reasons"])
        self.assertEqual(audit[1]["source_fields"]["Crosslinker"], "other chemistry")

    def test_pipeline_export_does_not_promote_unresolved_identities(self):
        bundle = self.read()
        mappings = self.mapping(bundle)
        mappings["customA"].update(status="ambiguous", accession=None)
        output = self.export_pipeline(bundle, mappings)
        manifest = json.loads((output / "manifest.json").read_bytes())
        self.assertEqual(manifest["counts"]["pipeline_input_rows"], 0)
        with (output / "claudio_input.csv").open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            self.assertEqual(reader.fieldnames[:8], CLAUDIO_COLUMNS)
            self.assertEqual(list(reader), [])
        audit = json.loads((output / "psm_audit.json").read_bytes())
        self.assertIn("protein_identity_unresolved_a", audit[0]["input_reasons"])

    def test_pipeline_export_still_applies_the_requested_fdr_level(self):
        bundle = self.read([psm(), psm(LinkFDR="0.02"), psm(LinkFDR="Infinity")])
        output = self.export_pipeline(bundle)
        audit = json.loads((output / "psm_audit.json").read_bytes())
        self.assertEqual([row["input_status"] for row in audit], ["ready", "not_evaluable", "not_evaluable"])
        self.assertIn("selected_fdr_above_cutoff", audit[1]["input_reasons"])
        self.assertIn("selected_fdr_non_finite", audit[2]["input_reasons"])

    def test_pipeline_export_records_policy_checksums_and_requires_a_new_directory(self):
        bundle = self.read()
        output = self.export_pipeline(bundle)
        manifest = json.loads((output / "manifest.json").read_bytes())
        self.assertEqual(manifest["kind"], "xifdr_claudio_input")
        self.assertEqual(manifest["pipeline_input"]["projections"],
                         ",".join(CLAUDIO_COLUMNS))
        self.assertTrue(manifest["pipeline_input"]["requires_source_preservation"])
        self.assertEqual(manifest["selection_policy"]["allowed_residues"], ["K"])
        self.assertFalse(manifest["assessment_performed"])
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)
        with self.assertRaises(FileExistsError):
            self.export_pipeline(bundle)

    def test_pipeline_export_rejects_misspelled_or_empty_selection_policies(self):
        bundle = self.read()
        for field, value in [("allowed_crosslinkers", "BS3"), ("allowed_crosslinkers", set()),
                             ("allowed_residues", "K"), ("allowed_residues", {"Lys"}),
                             ("source_protein_ids", {"custom_typo"}), ("source_protein_ids", set())]:
            options = {"allowed_crosslinkers": {"BS3"}, "allowed_residues": {"K"}, field: value}
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                write_claudio_input(bundle, self.mapping(bundle), self.root / "pipeline",
                                    fdr_level="link", fdr_cutoff=0.01, **options)
        self.assertFalse((self.root / "pipeline").exists())

    def test_reads_both_endpoints_and_preserves_original_fields(self):
        bundle = self.read()
        record = bundle["records"][0]
        self.assertEqual(record["record_status"], "valid")
        self.assertEqual(record["endpoints"]["a"]["peptide"], "ACK")
        self.assertEqual(record["endpoints"]["a"]["protein_position"], 4)
        self.assertEqual(record["endpoints"]["b"]["protein_position"], 3)
        self.assertEqual(record["source_fields"]["PepSeq1"], "ACcmK")
        self.assertEqual(record["source_fields"]["fdrGroup"], "between > 4 ")
        self.assertEqual(record["source_record_number"], 1)
        self.assertEqual(bundle["provenance"]["psm_sha256"], hashlib.sha256((self.root / "PSM.csv").read_bytes()).hexdigest())

    def test_duplicate_psm_ids_and_experimental_runs_keep_distinct_source_rows(self):
        bundle = self.read([psm(), psm(run="replicate-2"), psm()])
        self.assertEqual(len(bundle["records"]), 3)
        self.assertEqual(len({r["source_record_id"] for r in bundle["records"]}), 3)
        self.assertEqual([r["source_fields"]["run"] for r in bundle["records"]], ["replicate-1", "replicate-2", "replicate-1"])

    def test_source_ids_are_stable_when_input_files_are_relocated(self):
        original = self.read()
        copied = self.root / "copy.csv"
        copied.write_bytes((self.root / "PSM.csv").read_bytes())
        moved = read_xifdr_psms(copied, self.fasta, dataset_accession="PXD000001")
        self.assertEqual(original["records"][0]["source_record_id"], moved["records"][0]["source_record_id"])

    def test_identity_lookup_does_not_require_all_unused_search_database_proteins(self):
        self.fasta.write_text(FASTA + ">unused\nAAAAAA\n>DECOY:customA\nKATRKCAM\n")
        bundle = self.read()
        self.assertEqual(set(target_proteins(bundle)), {"customA", "customB"})

    def test_residue_and_peptide_positions_must_agree(self):
        for change, reason in [({"ProteinLinkPos1": "5"}, "coordinate_inconsistent"),
                               ({"PepPos2": "4", "ProteinLinkPos2": "4"}, "peptide_sequence_mismatch"),
                               ({"LinkPos1": "0"}, "invalid_position"),
                               ({"PeptideLength1": "4"}, "peptide_length_mismatch")]:
            with self.subTest(change=change):
                row = self.read([psm(**change)])["records"][0]
                self.assertEqual(row["record_status"], "not_evaluable")
                self.assertIn(reason, [e["reason"] for e in row["endpoints"].values()])

    def test_missing_and_ambiguous_protein_ids_remain_in_the_audit(self):
        for value in ["not_in_fasta", "customA;customB"]:
            with self.subTest(value=value):
                record = self.read([psm(Protein1=value)])["records"][0]
                self.assertEqual(record["record_status"], "not_evaluable")
                self.assertEqual(record["source_fields"]["Protein1"], value)

    def test_known_xi_modifications_preserve_peptide_residue_numbering(self):
        sequence, modifications = normalize_xifdr_peptide("ACcmMoxKbs3ohR")
        self.assertEqual(sequence, "ACMKR")
        self.assertEqual([(m["position"], m["token"]) for m in modifications], [(2, "cm"), (3, "ox"), (4, "bs3oh")])
        self.assertEqual(normalize_xifdr_peptide("AKbs3nh2R")[0], "AKR")

    def test_unknown_modifications_are_not_silently_stripped(self):
        for peptide in ["AK[+138]R", "AKunknownR", "ack"]:
            with self.subTest(peptide=peptide), self.assertRaises(ValueError):
                normalize_xifdr_peptide(peptide)
        record = self.read([psm(PepSeq1="ACunknownK")])["records"][0]
        self.assertEqual(record["endpoints"]["a"]["reason"], "unsupported_peptide_annotation")

    def test_decoys_are_preserved_but_never_enter_af_input(self):
        bundle = self.read([psm(), psm(PSMID="decoy-1", Protein1="DECOY:customA", Decoy1="true",
                                      isDecoy="true", isTT="false", isTD="true")])
        prepared = prepare_af_input(bundle, self.mapping(bundle), fdr_level="link", fdr_cutoff=0.01)
        self.assertEqual(len(prepared["audit"]), 2)
        self.assertEqual(len(prepared["rows"]), 1)
        self.assertIn("decoy_record", prepared["audit"][1]["input_reasons"])

    def test_inconsistent_decoy_flags_do_not_become_targets(self):
        bundle = self.read([psm(isDecoy="true")])
        self.assertEqual(bundle["records"][0]["decoy_status"], "unknown")
        prepared = prepare_af_input(bundle, self.mapping(bundle), fdr_level="link", fdr_cutoff=0.01)
        self.assertFalse(prepared["rows"])
        self.assertIn("invalid_decoy_flags", prepared["audit"][0]["input_reasons"])

    def test_fdr_levels_are_kept_separate_and_require_an_explicit_policy(self):
        bundle = self.read()
        mappings = self.mapping(bundle)
        accepted = prepare_af_input(bundle, mappings, fdr_level="link", fdr_cutoff=0.01)
        rejected = prepare_af_input(bundle, mappings, fdr_level="psm", fdr_cutoff=0.01)
        self.assertEqual(len(accepted["rows"]), 1)
        self.assertEqual(rejected["rows"], [])
        self.assertEqual(bundle["records"][0]["fdr"]["protein_pair"]["value"], 0.125)
        for level, cutoff in [("wrong", 0.01), ("link", float("nan")), ("link", True), ("link", 5)]:
            with self.subTest(level=level, cutoff=cutoff), self.assertRaises(ValueError):
                prepare_af_input(bundle, mappings, fdr_level=level, fdr_cutoff=cutoff)

    def test_infinity_is_preserved_as_text_and_is_not_a_numeric_fdr(self):
        bundle = self.read([psm(Protein1FDR="Infinity")])
        value = bundle["records"][0]["fdr"]["protein_a"]
        self.assertEqual(value, {"raw": "Infinity", "value": None, "status": "non_finite"})
        json.dumps(bundle, allow_nan=False)
        prepared = prepare_af_input(bundle, self.mapping(bundle), fdr_level="protein_a", fdr_cutoff=0.01)
        self.assertFalse(prepared["rows"])

    def test_blank_header_cells_are_allowed_but_duplicate_named_columns_are_not(self):
        bundle = self.read(columns=COLUMNS + [""])
        self.assertEqual(len(bundle["records"][0]["source_values"]), len(COLUMNS) + 1)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.read(columns=COLUMNS + ["Protein1"])

    def test_ragged_rows_are_recorded_without_assigning_shifted_values(self):
        values = [psm().get(k, "") for k in COLUMNS][:-1]
        record = self.read([values])["records"][0]
        self.assertEqual(record["record_reason"], "row_width_mismatch")
        self.assertEqual(record["source_values"], values)

    def test_links_table_is_not_misread_as_a_psm_table(self):
        with self.assertRaisesRegex(ValueError, "PSM"):
            self.read(columns=["Protein1", "Protein2", "fromSite", "ToSite", "fdr"])

    def test_fasta_ids_are_checked_and_uniprot_pipe_ids_have_an_alias(self):
        path = self.root / "pipe.fasta"
        path.write_text(REFERENCES)
        proteins = read_search_fasta(path)
        self.assertEqual(proteins["sp|P12345|A"]["uniprot_accession"], "P12345")
        path.write_text(">same\nACK\n>same\nACK\n")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            read_search_fasta(path)

    def test_reference_sequences_map_custom_ids_without_changing_coordinates(self):
        bundle = self.read()
        mappings = self.mapping(bundle)
        self.assertEqual(mappings["customA"]["accession"], "P12345")
        prepared = prepare_af_input(bundle, mappings, fdr_level="link", fdr_cutoff=0.01)
        row = prepared["rows"][0]
        self.assertEqual((row["unip_id_a"], row["pos_a"], row["seq_a"]), ("P12345", 4, "MACKRTAK"))
        self.assertEqual(row["source_protein_id_a"], "customA")
        self.assertEqual(row["position_status_a"], "valid")
        self.assertEqual(row["source_fdr_group"], "between > 4 ")
        self.assertNotIn("XL_type", row)

    def test_identical_sequences_with_multiple_accessions_remain_ambiguous(self):
        self.refs.write_text(REFERENCES + ">sp|P11111|C OX=9031\nMACKRTAK\n")
        bundle = self.read()
        mappings = self.mapping(bundle)
        self.assertEqual(mappings["customA"]["status"], "ambiguous")
        self.assertEqual({c["accession"] for c in mappings["customA"]["candidates"]}, {"P12345", "P11111"})
        self.assertFalse(prepare_af_input(bundle, mappings, fdr_level="link", fdr_cutoff=0.01)["rows"])

    def test_tags_and_species_mismatches_are_not_automatically_removed(self):
        self.refs.write_text(">sp|P12345|A OX=9606\nMACKRTAK\n")
        self.assertEqual(self.mapping(self.read())["customA"]["status"], "unmapped")
        self.refs.write_text(">sp|P12345|A OX=9031\nMACKRTAKWSHPQFEK\n")
        identity = self.mapping(self.read())["customA"]
        self.assertEqual(identity["status"], "unmapped")
        self.assertEqual(identity["sources"][0]["sha256"], hashlib.sha256(self.refs.read_bytes()).hexdigest())

    def test_a_mapped_label_is_not_enough_if_reference_sequence_is_wrong(self):
        bundle = self.read()
        mappings = self.mapping(bundle)
        mappings["customA"]["sequence"] = "WRONG"
        prepared = prepare_af_input(bundle, mappings, fdr_level="link", fdr_cutoff=0.01)
        self.assertFalse(prepared["rows"])
        self.assertIn("protein_identity_unresolved_a", prepared["audit"][0]["input_reasons"])

    def test_changed_identity_sources_are_rejected_before_publishing_outputs(self):
        bundle = self.read()
        mappings = self.mapping(bundle)
        prepared = prepare_af_input(bundle, mappings, fdr_level="link", fdr_cutoff=0.01)
        self.refs.write_text(">changed\nACK\n")
        with self.assertRaisesRegex(ValueError, "source changed"):
            write_xifdr_outputs(bundle, mappings, prepared, self.root / "prepared")
        self.assertFalse((self.root / "prepared").exists())

    def test_exports_keep_all_sources_and_only_send_eligible_rows_to_af(self):
        bundle = self.read([psm(), psm(ProteinLinkPos1="5")])
        mappings = self.mapping(bundle)
        prepared = prepare_af_input(bundle, mappings, fdr_level="link", fdr_cutoff=0.01)
        output = self.root / "prepared"
        write_xifdr_outputs(bundle, mappings, prepared, output)
        audit = json.loads((output / "psm_audit.json").read_bytes())
        self.assertEqual(len(audit), 2)
        self.assertTrue((output / "af_input.csv").is_file())
        manifest = json.loads((output / "manifest.json").read_bytes())
        self.assertEqual(manifest["counts"]["af_input_rows"], 1)
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)
        before = (output / "manifest.json").read_bytes()
        with self.assertRaises(FileExistsError):
            write_xifdr_outputs(bundle, mappings, prepared, output)
        self.assertEqual((output / "manifest.json").read_bytes(), before)

    @patch("claudio.dataset_catalog.xifdr._af_backend", side_effect=RuntimeError("AF assessment branch required"))
    def test_missing_af_backend_is_reported_before_creating_assessment_outputs(self, backend):
        bundle = self.read()
        prepared = prepare_af_input(bundle, self.mapping(bundle), fdr_level="link", fdr_cutoff=0.01)
        policy = {"minimum": 0, "maximum": 35, "boundary": "inclusive"}
        with self.assertRaisesRegex(RuntimeError, "AF assessment branch"):
            assess_af_complexes(prepared, self.root / "models", xl_residues="K:CB:", allowed_crosslinkers={"BS3"},
                plddt_cutoff=70, topolink_bin="topolink", euclidean_policy=policy, topological_policy=policy)
        self.assertFalse((self.root / "models").exists())

    @patch("claudio.dataset_catalog.xifdr._af_backend")
    def test_crosslinker_policy_is_checked_before_running_models(self, backend):
        bundle = self.read([psm(Crosslinker="EDC-DE")])
        prepared = prepare_af_input(bundle, self.mapping(bundle), fdr_level="link", fdr_cutoff=0.01)
        policy = {"minimum": 0, "maximum": 35, "boundary": "inclusive"}
        with self.assertRaisesRegex(ValueError, "crosslinker"):
            assess_af_complexes(prepared, self.root / "models", xl_residues="K:CB:", allowed_crosslinkers={"BS3"},
                plddt_cutoff=70, topolink_bin="topolink", euclidean_policy=policy, topological_policy=policy)
        backend.assert_not_called()

    @patch("claudio.dataset_catalog.xifdr._af_backend")
    def test_af_adapter_renames_distances_and_preserves_model_relationship_classification(self, backend):
        bundle = self.read()
        prepared = prepare_af_input(bundle, self.mapping(bundle), fdr_level="link", fdr_cutoff=0.01)
        source_id = prepared["rows"][0]["source_record_id"]

        def table(name):
            value = MagicMock(name=name)
            value.columns = {"source_record_id"}
            value.__getitem__.return_value = [source_id]
            value.__len__.return_value = 1
            return value

        data, candidates, calculated, classified, assessed = [table(n) for n in ("data", "candidates", "calculated", "classified", "assessed")]
        pandas = Mock()
        pandas.DataFrame.return_value = data
        retrieve, distances = Mock(return_value=candidates), Mock(return_value=calculated)
        classify, assess = Mock(return_value=classified), Mock(return_value=assessed)
        backend.return_value = pandas, retrieve, distances, classify, assess
        calculated.pop.side_effect = [[10.0], [12.0]]
        assessed.to_csv.side_effect = lambda path, **kwargs: Path(path).write_text("source_record_id\n" + source_id + "\n")
        policy = {"minimum": 0, "maximum": 35, "boundary": "inclusive"}
        result = assess_af_complexes(prepared, self.root / "models", xl_residues="K:CB:", allowed_crosslinkers={"BS3"},
            plddt_cutoff=70, topolink_bin="topolink", euclidean_policy=policy, topological_policy=policy, topolink_seed=23)
        self.assertIs(result, assessed)
        self.assertEqual(calculated.__setitem__.call_args_list, [call("eucl_dist", [10.0]), call("topo_dist", [12.0])])
        classify.assert_called_once_with(calculated)
        assess.assert_called_once_with(classified, euclidean_policy=policy, topological_policy=policy, plddt_cutoff=70)
        self.assertEqual(distances.call_args.kwargs["topolink_seed"], 23)
        self.assertTrue((self.root / "models/manifest.json").is_file())


class TestUniProtSequenceLookup(unittest.TestCase):
    @patch("claudio.dataset_catalog.xifdr.requests.get")
    def test_exact_sequence_lookup_keeps_active_references_and_raw_sources(self, get):
        sequence = "MACKRTAK"
        references = [
            {"database": "UniProtKB/Swiss-Prot", "id": "P12345", "active": True, "organism": {"taxonId": 9031}},
            {"database": "UniProtKB/Swiss-Prot", "id": "P11111", "active": False, "organism": {"taxonId": 9031}},
            {"database": "UniProtKB/TrEMBL", "id": "Q99999", "active": True, "organism": {"taxonId": 9606}},
        ]
        get.side_effect = [response("search", {"results": [{"uniParcId": "UPI0000000001", "sequence": {"value": sequence}}]}),
            response("entry", {"uniParcId": "UPI0000000001", "sequence": {"value": sequence}, "uniParcCrossReferences": references}),
            response("fasta", ">sp|P12345|A OX=9031\n" + sequence + "\n")]
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "search.fasta"
            fasta.write_text(">customA\n" + sequence + "\n")
            mappings = resolve_uniprot_proteins(read_search_fasta(fasta), root / "lookup", taxon_id=9031)
            self.assertEqual(mappings["customA"]["status"], "mapped")
            self.assertEqual(mappings["customA"]["accession"], "P12345")
            self.assertEqual(get.call_count, 3)
            for source in mappings["customA"]["sources"]:
                self.assertEqual(hashlib.sha256((root / "lookup" / source["path"]).read_bytes()).hexdigest(), source["sha256"])
        self.assertEqual(get.call_args_list[0].kwargs["params"]["query"], "checksum:" + hashlib.md5(sequence.encode()).hexdigest().upper())

    @patch("claudio.dataset_catalog.xifdr.requests.get")
    def test_sequence_mismatch_is_not_accepted_just_because_the_checksum_query_returned_it(self, get):
        get.return_value = response("search", {"results": [{"uniParcId": "UPI0000000001", "sequence": {"value": "WRONG"}}]})
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "search.fasta"
            fasta.write_text(">customA\nMACKRTAK\n")
            mappings = resolve_uniprot_proteins(read_search_fasta(fasta), root / "lookup")
            self.assertEqual(mappings["customA"]["status"], "unmapped")
            self.assertIsNone(mappings["customA"]["accession"])
        self.assertEqual(get.call_count, 1)

    @patch("claudio.dataset_catalog.xifdr.requests.get", side_effect=requests.exceptions.Timeout("lookup timed out"))
    def test_network_failure_does_not_exclude_or_fabricate_a_protein_identity(self, get):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "search.fasta"
            fasta.write_text(">customA\nMACKRTAK\n")
            mappings = resolve_uniprot_proteins(read_search_fasta(fasta), root / "lookup")
            self.assertEqual(mappings["customA"]["status"], "lookup_failed")
            self.assertIn("timed out", mappings["customA"]["reason"])
            self.assertIsNone(mappings["customA"]["accession"])

    @patch("claudio.dataset_catalog.xifdr.requests.get")
    def test_an_incomplete_search_is_not_reported_as_a_unique_identity(self, get):
        first = response("search", {"results": []})
        first.headers["Link"] = '<https://rest.uniprot.org/uniparc/search?cursor=next>; rel="next"'
        get.side_effect = [first, requests.exceptions.Timeout("second page timed out")]
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "search.fasta"
            fasta.write_text(">customA\nMACKRTAK\n")
            mappings = resolve_uniprot_proteins(read_search_fasta(fasta), root / "lookup")
            self.assertEqual(mappings["customA"]["status"], "lookup_failed")
            self.assertIsNone(mappings["customA"]["accession"])

    @patch("claudio.dataset_catalog.xifdr.requests.get")
    def test_current_uniprot_sequence_must_still_match_the_uniparc_reference(self, get):
        sequence = "MACKRTAK"
        get.side_effect = [response("search", {"results": [{"uniParcId": "UPI0000000001", "sequence": {"value": sequence}}]}),
            response("entry", {"uniParcId": "UPI0000000001", "sequence": {"value": sequence}, "uniParcCrossReferences": [
                {"database": "UniProtKB/Swiss-Prot", "id": "P12345", "active": True, "organism": {"taxonId": 9031}}]}),
            response("fasta", ">sp|P12345|A OX=9031\nMACKRTAKG\n")]
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "search.fasta"
            fasta.write_text(">customA\n" + sequence + "\n")
            mappings = resolve_uniprot_proteins(read_search_fasta(fasta), root / "lookup", taxon_id=9031)
            self.assertEqual(mappings["customA"]["status"], "unmapped")
            self.assertIsNone(mappings["customA"]["accession"])


class TestReferenceCoordinateComparison(unittest.TestCase):
    def test_full_source_inside_reference_has_an_explicit_offset(self):
        report = compare_protein_sequences("ACDEFGHIK", "MNPQACDEFGHIK")
        self.assertEqual(report["mapping_method"], "full_source_exact_subsequence")
        self.assertEqual(report["identical_position_map"]["1"], 5)
        self.assertEqual(report["identical_position_map"]["9"], 13)
        self.assertTrue(report["coordinate_mapping_usable"])

    def test_terminal_tag_positions_do_not_receive_a_reference_residue(self):
        report = compare_protein_sequences("ACDEFGHIKWSHPQFEK", "ACDEFGHIK")
        self.assertEqual(report["mapping_method"], "full_reference_exact_subsequence")
        self.assertEqual(report["identical_position_map"]["9"], 9)
        self.assertNotIn("10", report["identical_position_map"])
        self.assertEqual(report["source_unmapped_positions"], list(range(10, 18)))

    def test_substitutions_have_coordinates_but_are_not_identical_residues(self):
        report = compare_protein_sequences("ACDEFGHIK", "ACDEYGHIK")
        self.assertTrue(report["coordinate_mapping_usable"])
        self.assertEqual(report["aligned_position_map"]["5"], 5)
        self.assertNotIn("5", report["identical_position_map"])
        self.assertEqual(report["mismatches"], [{"source_position": 5, "reference_position": 5, "source_residue": "F", "reference_residue": "Y"}])

    def test_internal_reference_insertions_shift_only_the_following_positions(self):
        report = compare_protein_sequences("ACDEFGHIK", "ACDEWFGHIK")
        self.assertEqual(report["identical_position_map"]["4"], 4)
        self.assertEqual(report["identical_position_map"]["5"], 6)
        self.assertEqual(report["identical_position_map"]["9"], 10)
        self.assertEqual(report["reference_unmapped_positions"], [5])
        self.assertAlmostEqual(report["identity_percent"], 90)

    def test_repeated_subsequences_do_not_get_an_arbitrary_position_map(self):
        report = compare_protein_sequences("AAAA", "AAAAA")
        self.assertFalse(report["coordinate_mapping_usable"])
        self.assertEqual(report["identical_position_map"], {})
        self.assertIn("ambiguous", report["mapping_method"])

    def test_identity_denominator_includes_gap_columns(self):
        report = compare_protein_sequences("ACDEFGHIK", "MNPQACDEFGHIK")
        self.assertEqual(report["alignment_columns"], 13)
        self.assertEqual(report["identical_residues"], 9)
        self.assertAlmostEqual(report["identity_percent"], 100 * 9 / 13)
        self.assertEqual(report["source_coverage_percent"], 100)
        self.assertEqual(report["identity_scope"], "global_alignment_including_terminal_gaps")

    def test_invalid_or_pre_gapped_sequences_are_rejected(self):
        for source, target in [("", "ACK"), ("AC-K", "ACK"), ("ack", "ACK")]:
            with self.subTest(source=source), self.assertRaises(ValueError):
                compare_protein_sequences(source, target)


class TestGeneCandidateSearch(unittest.TestCase):
    def protein(self):
        return {"customA": {"id": "customA", "sequence": "ACDEFGHIK", "uniprot_accession": None}}

    @patch("claudio.dataset_catalog.xifdr.requests.get")
    def test_candidates_keep_taxonomy_gene_hint_and_alignment_evidence_without_mapping_identity(self, get):
        get.return_value = response("search", {"results": [{"primaryAccession": "P12345",
            "organism": {"taxonId": 9031}, "genes": [{"geneName": {"value": "FANCA"}}],
            "sequence": {"value": "MNPQACDEFGHIK"}, "entryType": "UniProtKB unreviewed (TrEMBL)"}]})
        with TemporaryDirectory() as tmp:
            result = search_uniprot_candidates(self.protein(), Path(tmp) / "candidates", taxon_id=9031, gene_hints={"customA": "FANCA"})
            report = result["customA"]
            self.assertEqual(report["status"], "candidates_for_review")
            self.assertEqual(report["gene_hint"], "FANCA")
            self.assertEqual(report["candidates"][0]["accession"], "P12345")
            self.assertEqual(report["candidates"][0]["comparison"]["identical_position_map"]["1"], 5)
            self.assertNotIn("accession", report)
            self.assertTrue((Path(tmp) / "candidates/candidate_report.json").is_file())
        self.assertEqual(get.call_args.kwargs["params"]["query"], "(gene_exact:FANCA) AND (organism_id:9031)")

    @patch("claudio.dataset_catalog.xifdr.requests.get")
    def test_different_species_or_missing_gene_hint_does_not_get_an_accepted_candidate(self, get):
        get.return_value = response("search", {"results": [{"primaryAccession": "P12345", "organism": {"taxonId": 9606}, "sequence": {"value": "ACDEFGHIK"}}]})
        with TemporaryDirectory() as tmp:
            result = search_uniprot_candidates(self.protein(), Path(tmp) / "wrong_species", taxon_id=9031, gene_hints={"customA": "FANCA"})
            self.assertEqual(result["customA"]["status"], "no_candidates")
            get.reset_mock()
            result = search_uniprot_candidates(self.protein(), Path(tmp) / "no_hint", taxon_id=9031, gene_hints={})
            self.assertEqual(result["customA"]["status"], "gene_hint_required")
            get.assert_not_called()

    @patch("claudio.dataset_catalog.xifdr.requests.get")
    def test_query_hints_are_checked_before_network_or_output_creation(self, get):
        with TemporaryDirectory() as tmp:
            target = Path(tmp) / "candidates"
            with self.assertRaises(ValueError):
                search_uniprot_candidates(self.protein(), target, taxon_id=9031, gene_hints={"customA": "FANCA OR *"})
            self.assertFalse(target.exists())
        get.assert_not_called()

    @patch("claudio.dataset_catalog.xifdr.requests.get", side_effect=requests.exceptions.Timeout("candidate lookup timed out"))
    def test_candidate_query_failure_is_kept_separate_from_no_hits(self, get):
        with TemporaryDirectory() as tmp:
            result = search_uniprot_candidates(self.protein(), Path(tmp) / "candidates", taxon_id=9031, gene_hints={"customA": "FANCA"})
            self.assertEqual(result["customA"]["status"], "lookup_failed")
            self.assertIn("timed out", result["customA"]["reason"])


class TestIdentityFailureReasons(unittest.TestCase):
    @patch("claudio.dataset_catalog.xifdr.requests.get")
    def test_absent_sequences_and_sequences_without_uniprot_references_have_different_reasons(self, get):
        proteins = {"customA": {"id": "customA", "sequence": "MACKRTAK", "uniprot_accession": None}}
        with TemporaryDirectory() as tmp:
            get.return_value = response("search", {"results": []})
            absent = resolve_uniprot_proteins(proteins, Path(tmp) / "absent")
            self.assertEqual(absent["customA"]["reason"], "sequence_not_found_in_uniparc")
            get.side_effect = [response("search", {"results": [{"uniParcId": "UPI0000000001", "sequence": {"value": "MACKRTAK"}}]}),
                response("entry", {"uniParcId": "UPI0000000001", "sequence": {"value": "MACKRTAK"},
                    "uniParcCrossReferences": [{"database": "RefSeq", "id": "XP_example", "active": True}]})]
            missing_reference = resolve_uniprot_proteins(proteins, Path(tmp) / "missing_reference")
            self.assertEqual(missing_reference["customA"]["reason"], "no_active_uniprot_reference")
            self.assertEqual(missing_reference["customA"]["lookup_details"]["exact_uniparc_entries"], 1)


class TestXiFDRCandidateProjection(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / "search.fasta"
        self.refs = self.root / "uniprot.fasta"
        self.refs.write_text(">sp|Q99999|B OX=9031\nGGKPEPK\n")

    def bundle(self, source="MACKRTAK", rows=None):
        self.fasta.write_text(">customA\n" + source + "\n>customB\nGGKPEPK\n")
        path = self.root / "PSM.csv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(COLUMNS)
            for row in [psm()] if rows is None else rows:
                writer.writerow([row.get(k, "") for k in COLUMNS])
        return read_xifdr_psms(path, self.fasta, dataset_accession="PXD000001")

    def exact(self, bundle):
        return map_reference_proteins(bundle["proteins"], read_search_fasta(self.refs), taxon_id=9031)

    def reports(self, *references):
        return {"customA": {"status": "candidates_for_review", "taxon_id": 9031, "gene_hint": "FANCA", "sources": [],
                "candidates": [{"accession": accession, "sequence": sequence, "taxon_id": 9031} for accession, sequence in references]}}

    def test_exact_reference_offsets_change_only_reference_coordinates(self):
        bundle = self.bundle()
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "GGGMACKRTAK")))
        row = result["records"][0]
        self.assertEqual(row["status"], "projected_candidates")
        candidate = row["endpoint_candidates"]["a"][0]
        self.assertEqual(candidate["source_protein_position"], 4)
        self.assertEqual(candidate["reference_protein_position"], 7)
        self.assertEqual(candidate["reference_peptide_start"], 5)
        self.assertEqual(bundle["records"][0]["endpoints"]["a"]["protein_position"], 4)
        self.assertEqual(row["pairs"][0]["identity_scope"], "candidate_hypothesis")
        self.assertFalse(result["protein_identity_acceptance_performed"])

    def test_ambiguity_outside_the_peptide_does_not_block_a_stable_projection(self):
        row = psm(PepSeq1="FGHIK", PepPos1="6", PeptideLength1="5", LinkPos1="5", ProteinLinkPos1="10")
        bundle = self.bundle("ACDSSFGHIK", [row])
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "ACDSFGHIK")))
        candidate = result["records"][0]["endpoint_candidates"]["a"][0]
        self.assertGreater(candidate["optimal_alignment_count"], 1)
        self.assertEqual(candidate["status"], "projected")
        self.assertEqual(candidate["reference_protein_position"], 9)
        self.assertEqual(candidate["distinct_projection_count"], 1)

    def test_a_peptide_spanning_the_missing_residue_is_not_projected(self):
        row = psm(PepSeq1="SSFGHIK", PepPos1="4", PeptideLength1="7", LinkPos1="7", ProteinLinkPos1="10")
        bundle = self.bundle("ACDSSFGHIK", [row])
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "ACDSFGHIK")))
        candidate = result["records"][0]["endpoint_candidates"]["a"][0]
        self.assertEqual(candidate["status"], "not_evaluable")
        self.assertEqual(candidate["compatible_alignment_count"], 0)
        self.assertEqual(result["records"][0]["pairs"], [])

    def test_repeated_peptides_with_multiple_reference_positions_remain_ambiguous(self):
        row = psm(PepSeq1="FGHIK", PepPos1="5", PeptideLength1="5", LinkPos1="5", ProteinLinkPos1="9")
        bundle = self.bundle("ACDEFGHIK", [row])
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "ACDEFGHIKACDEFGHIK")))
        candidate = result["records"][0]["endpoint_candidates"]["a"][0]
        self.assertEqual(candidate["reason"], "ambiguous_reference_position")
        self.assertGreater(candidate["distinct_projection_count"], 1)
        self.assertEqual(result["records"][0]["pairs"], [])

    def test_a_tag_peptide_is_retained_without_a_native_protein_position(self):
        row = psm(PepSeq1="WSHPQFEK", PepPos1="9", PeptideLength1="8", LinkPos1="8", ProteinLinkPos1="16")
        bundle = self.bundle("MACKRTAKWSHPQFEK", [row])
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "MACKRTAK")))
        candidate = result["records"][0]["endpoint_candidates"]["a"][0]
        self.assertEqual(candidate["status"], "not_evaluable")
        self.assertIsNone(candidate["reference_protein_position"])
        self.assertIn("peptide_overlaps_unmapped_source_residues", candidate["rejected_alignment_reasons"])

    def test_sequence_changes_inside_the_observed_peptide_are_not_accepted(self):
        bundle = self.bundle()
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "MACERTAK")))
        candidate = result["records"][0]["endpoint_candidates"]["a"][0]
        self.assertEqual(candidate["status"], "not_evaluable")
        self.assertIn("peptide_sequence_differs_from_reference", candidate["rejected_alignment_reasons"])

    def test_changes_outside_the_observed_peptide_keep_a_conditional_projection(self):
        bundle = self.bundle()
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "MACKRVAK")))
        self.assertEqual(result["records"][0]["endpoint_candidates"]["a"][0]["status"], "projected")
        self.assertEqual(result["records"][0]["pairs"][0]["identity_scope"], "candidate_hypothesis")

    def test_one_source_protein_cannot_become_two_different_reference_identities_in_a_pair(self):
        row = psm(Protein2="customA", PepSeq2="RTAK", PepPos2="5", PeptideLength2="4", LinkPos2="4", ProteinLinkPos2="8")
        bundle = self.bundle(rows=[row])
        reports = self.reports(("P12345", "MACKRTAK"), ("P11111", "MACKRTAK"))
        result = project_xifdr_candidates(bundle, self.exact(bundle), reports)
        pairs = result["records"][0]["pairs"]
        self.assertEqual(len(pairs), 2)
        self.assertTrue(all(p["unip_id_a"] == p["unip_id_b"] for p in pairs))
        self.assertEqual(result["records"][0]["discarded_cross_identity_pair_count"], 2)

    def test_decoys_invalid_rows_and_duplicate_sources_stay_in_the_projection_audit(self):
        bundle = self.bundle(rows=[psm(), psm(), psm(ProteinLinkPos1="5"),
            psm(Decoy1="true", isDecoy="true", isTT="false", isTD="true")])
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "MACKRTAK")))
        self.assertEqual([r["source_record_id"] for r in result["records"]], [r["source_record_id"] for r in bundle["records"]])
        self.assertEqual(len(result["records"]), 4)
        self.assertEqual(result["records"][2]["status"], "not_evaluable")
        self.assertEqual(result["records"][3]["status"], "excluded_decoy")
        ids = [p["projection_candidate_id"] for r in result["records"] for p in r["pairs"]]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertFalse(result["fdr_filter_applied"])

    def test_alignment_budget_is_reported_instead_of_choosing_the_first_alignment(self):
        row = psm(PepSeq1="FGHIK", PepPos1="6", PeptideLength1="5", LinkPos1="5", ProteinLinkPos1="10")
        bundle = self.bundle("ACDSSFGHIK", [row])
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "ACDSFGHIK")), max_optimal_alignments=1)
        candidate = result["records"][0]["endpoint_candidates"]["a"][0]
        self.assertEqual(candidate["reason"], "optimal_alignment_limit_exceeded")
        self.assertEqual(candidate["status"], "not_evaluable")
        self.assertEqual(result["records"][0]["pairs"], [])

    def test_projection_exports_preserve_source_rows_and_hashes(self):
        bundle = self.bundle()
        result = project_xifdr_candidates(bundle, self.exact(bundle), self.reports(("P12345", "GGGMACKRTAK")))
        output = self.root / "projection"
        write_projection_outputs(bundle, result, output)
        stored = json.loads((output / "projection_report.json").read_bytes())
        self.assertEqual(stored["records"], result["records"])
        self.assertTrue((output / "projection_pairs.csv").is_file())
        manifest = json.loads((output / "manifest.json").read_bytes())
        self.assertFalse(manifest["model_assessment_performed"])
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)
        with self.assertRaises(FileExistsError):
            write_projection_outputs(bundle, result, output)

    def test_retry_history_keeps_successes_and_does_not_hide_failed_attempts(self):
        bundle = self.bundle()
        success = self.reports(("P12345", "MACKRTAK"))
        failure = {"customA": {"status": "lookup_failed", "reason": "HTTP 500", "gene_hint": "FANCA", "taxon_id": 9031, "candidates": [], "sources": []}}
        result = project_xifdr_candidates(bundle, self.exact(bundle), [success, failure])
        self.assertEqual(result["records"][0]["status"], "projected_candidates")
        self.assertEqual(result["identity_evidence"]["candidate_report_history"][1]["customA"]["reason"], "HTTP 500")
        self.assertEqual(result["protein_lookup_states"]["customA"]["selected_candidate_report_index"], 0)
        recovered = project_xifdr_candidates(bundle, self.exact(bundle), [failure, success])
        self.assertEqual(recovered["protein_lookup_states"]["customA"]["selected_candidate_report_index"], 1)

class TestXiFDRPyXLMS(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.psms = self.root / "PSM.csv"
        self.fasta = self.root / "search.fasta"
        self.fasta.write_text(FASTA, encoding="utf-8")

    def write_psms(self, rows):
        # Older exports can omit CrosslinkerModMass while retaining charge.
        columns = COLUMNS + ["exp charge"]
        with self.psms.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(columns)
            for row in rows:
                writer.writerow([
                    row.get(column, "3" if column == "exp charge" else "")
                    for column in columns
                ])

    def read(self):
        return read_xifdr_psms(
            self.psms, self.fasta, dataset_accession="PXD000001", backend="pyxlms"
        )

    def reversed_psm(self, **changes):
        row = psm(
            Protein1="customB", Protein2="customA",
            PepSeq1="KPEP", PepSeq2="ACcmK", PepPos1="3", PepPos2="2",
            PeptideLength1="4", PeptideLength2="3", LinkPos1="1", LinkPos2="3",
            ProteinLinkPos1="3", ProteinLinkPos2="4",
        )
        row.update(changes)
        return row

    def test_reads_without_mass_and_keeps_source_modifications_and_fdr(self):
        self.write_psms([psm()])
        original = self.psms.read_bytes()
        with patch.object(pyxlms_parser, "read_xi", wraps=pyxlms_parser.read_xi) as reader:
            bundle = self.read()

        self.assertEqual(reader.call_count, 1)
        self.assertIs(reader.call_args.kwargs["parse_modifications"], False)
        self.assertEqual(bundle["parser"]["name"], "pyXLMS")
        self.assertEqual(bundle["parser"]["version"], version("pyXLMS"))
        self.assertFalse(bundle["parser"]["modification_masses_assessed"])
        record = bundle["records"][0]
        self.assertEqual(record["record_status"], "valid")
        self.assertEqual(record["endpoints"]["a"]["peptide"], "ACK")
        self.assertEqual(record["endpoints"]["a"]["modifications"], [
            {"position": 2, "residue": "C", "token": "cm"},
        ])
        self.assertEqual(record["source_fields"]["PepSeq1"], "ACcmK")
        self.assertNotIn("CrosslinkerModMass", record["source_fields"])
        self.assertEqual(record["fdr"]["psm"]["value"], 0.02)
        self.assertEqual(record["fdr"]["link"]["value"], 0.005)
        self.assertEqual(record["fdr"]["protein_pair"]["raw"], "0.125")
        self.assertEqual(self.psms.read_bytes(), original)

    def test_restores_the_original_endpoint_order(self):
        self.write_psms([self.reversed_psm()])
        bundle = self.read()
        record = bundle["records"][0]

        self.assertEqual(record["parser_orientation"], "swapped")
        self.assertEqual(record["endpoints"]["a"]["peptide"], "KPEP")
        self.assertEqual(record["endpoints"]["a"]["source_protein_id"], "customB")
        self.assertEqual(record["endpoints"]["a"]["protein_position"], 3)
        self.assertEqual(record["endpoints"]["b"]["source_peptide"], "ACcmK")
        self.assertEqual(record["endpoints"]["b"]["protein_position"], 4)
        self.assertEqual(record["source_fields"]["Protein1"], "customB")
        self.assertEqual(record["record_status"], "valid")

    def test_identical_peptides_are_matched_using_their_positions(self):
        self.fasta.write_text(">customA\nMACKACK\n", encoding="utf-8")
        self.write_psms([psm(
            Protein2="customA", PepSeq1="ACcmK", PepSeq2="ACcmK",
            PepPos1="2", PepPos2="5", PeptideLength2="3",
            LinkPos1="3", LinkPos2="1", ProteinLinkPos1="4", ProteinLinkPos2="5",
        )])
        record = self.read()["records"][0]

        self.assertEqual(record["record_status"], "valid")
        self.assertEqual(record["parser_orientation"], "swapped")
        self.assertEqual(record["endpoints"]["a"]["link_position"], 3)
        self.assertEqual(record["endpoints"]["b"]["link_position"], 1)
        self.assertEqual(record["endpoints"]["a"]["protein_position"], 4)
        self.assertEqual(record["endpoints"]["b"]["protein_position"], 5)

    def test_decoy_flags_follow_the_original_source_endpoints(self):
        self.write_psms([self.reversed_psm(
            Decoy1="true", isDecoy="true", isTT="false", isTD="true",
        )])
        bundle = self.read()
        record = bundle["records"][0]

        self.assertEqual(record["parser_orientation"], "swapped")
        self.assertEqual(record["decoy_status"], "target_decoy")
        self.assertEqual(record["source_fields"]["Decoy1"], "true")
        self.assertEqual(record["source_fields"]["Decoy2"], "false")
        self.assertEqual(record["endpoints"]["a"]["source_protein_id"], "customB")
        self.assertEqual(prepare_af_input(bundle, {}, fdr_level="link", fdr_cutoff=0.01)["rows"], [])

    def test_duplicate_psms_and_different_runs_keep_distinct_source_rows(self):
        row = self.reversed_psm()
        self.write_psms([row, row, self.reversed_psm(run="replicate-2")])
        bundle = self.read()
        records = bundle["records"]

        self.assertEqual(len(records), 3)
        self.assertEqual([r["source_record_number"] for r in records], [1, 2, 3])
        self.assertEqual(len({r["source_record_id"] for r in records}), 3)
        self.assertEqual([r["source_fields"]["run"] for r in records], [
            "replicate-1", "replicate-1", "replicate-2",
        ])
        self.assertTrue(all(r["record_status"] == "valid" for r in records))

    def test_invalid_and_ambiguous_rows_remain_in_the_audit(self):
        self.write_psms([
            psm(), psm(LinkPos1="0"),
            psm(Protein1="customA;customB", PepPos1="2;3", ProteinLinkPos1="4;5"),
        ])
        bundle = self.read()
        records = bundle["records"]

        self.assertEqual(len(records), 3)
        self.assertEqual(records[0]["record_status"], "valid")
        self.assertEqual(records[1]["record_status"], "not_evaluable")
        self.assertEqual(records[2]["record_status"], "not_evaluable")
        self.assertEqual(records[1]["source_fields"]["LinkPos1"], "0")
        self.assertEqual(records[2]["source_fields"]["Protein1"], "customA;customB")
        self.assertEqual([r["source_record_number"] for r in records], [1, 2, 3])

    def test_parser_failure_is_reported_without_accepting_native_results(self):
        self.write_psms([psm()])
        with patch.object(pyxlms_parser, "read_xi", side_effect=ValueError("unsupported export")):
            bundle = self.read()
        record = bundle["records"][0]

        self.assertEqual(record["record_status"], "not_evaluable")
        self.assertEqual(record["record_reason"], "pyxlms_parse_failed")
        self.assertEqual(record["source_fields"]["PSMID"], "psm-1")
        self.assertIn("unsupported export", record["parser_error"])

    def test_a_different_parsed_protein_cannot_pass_the_endpoint_check(self):
        self.write_psms([psm()])
        read_xi = pyxlms_parser.read_xi

        def change_protein(*args, **kwargs):
            result = read_xi(*args, **kwargs)
            first = result["crosslink-spectrum-matches"][0]
            fields = (
                "peptide", "peptide_crosslink_position", "proteins",
                "proteins_crosslink_positions", "proteins_peptide_positions", "decoy",
            )
            parsed = {f"{side}_{field}": first[f"{side}_{field}"]
                      for side in ("alpha", "beta") for field in fields}
            parsed["additional_information"] = first["additional_information"]
            parsed["alpha_proteins"] = ["another_protein"]
            return {"crosslink-spectrum-matches": [parsed], "crosslinks": None}

        with patch.object(pyxlms_parser, "read_xi", side_effect=change_protein):
            bundle = self.read()
        record = bundle["records"][0]

        self.assertEqual(record["record_status"], "not_evaluable")
        self.assertEqual(record["record_reason"], "pyxlms_endpoint_mismatch")
        self.assertEqual(record["source_fields"]["Protein1"], "customA")

    def test_decoy_prefix_is_preserved_as_part_of_the_fasta_identifier(self):
        self.fasta.write_text(FASTA + ">decoy:customA\nMACKRTAK\n", encoding="utf-8")
        self.write_psms([psm(
            Protein1="decoy:customA", Decoy1="true",
            isDecoy="true", isTT="false", isTD="true",
        )])
        record = self.read()["records"][0]

        self.assertEqual(record["record_status"], "valid")
        self.assertEqual(record["endpoints"]["a"]["source_fasta_id"], "decoy:customA")
        self.assertEqual(record["endpoints"]["a"]["source_protein_id"], "decoy:customA")
        self.assertEqual(record["decoy_status"], "target_decoy")

    def test_parser_result_order_does_not_change_source_associations(self):
        self.write_psms([psm(), self.reversed_psm(run="replicate-2")])
        read_xi = pyxlms_parser.read_xi

        def reverse_results(*args, **kwargs):
            result = read_xi(*args, **kwargs)
            return {"crosslink-spectrum-matches": result["crosslink-spectrum-matches"][::-1],
                    "crosslinks": None}

        with patch.object(pyxlms_parser, "read_xi", side_effect=reverse_results):
            records = self.read()["records"]

        self.assertTrue(all(record["record_status"] == "valid" for record in records))
        self.assertEqual([r["source_fields"]["run"] for r in records], ["replicate-1", "replicate-2"])
        self.assertEqual([r["endpoints"]["a"]["peptide"] for r in records], ["ACK", "KPEP"])
        self.assertTrue(all("__claudio_source_record_id" not in r["source_fields"] for r in records))

    def test_a_row_missing_from_parser_output_stays_unresolved(self):
        self.write_psms([psm(), psm(run="replicate-2")])
        read_xi = pyxlms_parser.read_xi

        def drop_first_result(*args, **kwargs):
            result = read_xi(*args, **kwargs)
            return {"crosslink-spectrum-matches": result["crosslink-spectrum-matches"][1:],
                    "crosslinks": None}

        with patch.object(pyxlms_parser, "read_xi", side_effect=drop_first_result):
            records = self.read()["records"]

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["record_status"], "not_evaluable")
        self.assertEqual(records[0]["record_reason"], "pyxlms_parse_failed")
        self.assertEqual(records[1]["record_status"], "valid")
        self.assertEqual(records[1]["source_fields"]["run"], "replicate-2")

class TestDatasetStandardizeCLI(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.psms = self.root / "PSM.csv"
        self.fasta = self.root / "search.fasta"
        self.references = self.root / "uniprot.fasta"
        self.mapping_path = self.root / "mapping.json"
        self.output = self.root / "standardized"
        self.fasta.write_text(FASTA, encoding="utf-8")
        self.references.write_text(REFERENCES, encoding="utf-8")
        self.prepare_inputs([psm(), psm(run="replicate-2"), psm(LinkPos1="0")])

    def prepare_inputs(self, rows):
        columns = COLUMNS + ["exp charge"]
        with self.psms.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(columns)
            for row in rows:
                writer.writerow([row.get(key, "3" if key == "exp charge" else "") for key in columns])
        bundle = read_xifdr_psms(self.psms, self.fasta, dataset_accession="PXD000001")
        mappings = map_reference_proteins(bundle["proteins"], read_search_fasta(self.references), taxon_id=9031)
        self.mapping_path.write_text(json.dumps(mappings), encoding="utf-8")

    def arguments(self, mapping=None, output=None):
        return [
            "standardize", "--results", str(self.psms), "--fasta", str(self.fasta),
            "--mapping", str(self.mapping_path if mapping is None else mapping),
            "--dataset-accession", "PXD000001",
            "--fdr-level", "link", "--fdr-cutoff", "0.01",
            "--crosslinker", "BS3", "--residue", "K",
            "--output", str(self.output if output is None else output),
        ]

    def run_command(self, arguments):
        result = CliRunner().invoke(catalog_main, arguments)
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        return result

    def test_exports_default_m1_columns_with_pyxlms_and_keeps_the_full_audit(self):
        result = self.run_command(self.arguments())

        self.assertIn("Source PSM records: 3", result.output)
        self.assertIn("Standardized input rows: 2", result.output)
        with (self.output / "claudio_input.csv").open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            self.assertEqual(reader.fieldnames[:8], CLAUDIO_COLUMNS)
            rows = list(reader)
        self.assertEqual([row["source_run"] for row in rows], ["replicate-1", "replicate-2"])
        self.assertEqual(len({row["xlms_source_record_id"] for row in rows}), 2)
        self.assertNotIn("source_record_id", reader.fieldnames)
        audit = json.loads((self.output / "psm_audit.json").read_text())
        self.assertEqual(len(audit), 3)
        self.assertEqual(audit[2]["record_status"], "not_evaluable")
        manifest = json.loads((self.output / "manifest.json").read_text())
        self.assertEqual(manifest["parser"]["name"], "pyXLMS")
        self.assertEqual(manifest["fdr_policy"]["level"], "link")
        self.assertFalse(manifest["assessment_performed"])
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((self.output / name).read_bytes()).hexdigest(), digest)

    def test_reads_identity_sources_from_a_previous_export(self):
        previous = self.root / "previous"
        self.run_command(self.arguments(output=previous))
        self.references.unlink()

        self.run_command(self.arguments(mapping=previous / "protein_mapping.json"))

        self.assertEqual((self.output / "claudio_input.csv").read_bytes(),
                         (previous / "claudio_input.csv").read_bytes())
        manifest = json.loads((self.output / "manifest.json").read_text())
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((self.output / name).read_bytes()).hexdigest(), digest)

    def test_existing_output_is_left_untouched(self):
        self.output.mkdir()
        marker = self.output / "keep.txt"
        marker.write_text("previous output", encoding="utf-8")

        result = CliRunner().invoke(catalog_main, self.arguments())

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("exist", result.output.lower())
        self.assertEqual(marker.read_text(), "previous output")
        self.assertEqual(list(self.output.iterdir()), [marker])

    def test_requires_the_user_to_specify_the_fdr_level_and_cutoff(self):
        arguments = self.arguments()
        for flag in ("--fdr-level", "--fdr-cutoff"):
            index = arguments.index(flag)
            del arguments[index:index + 2]

        result = CliRunner().invoke(catalog_main, arguments)

        self.assertEqual(result.exit_code, 2)
        self.assertIn("fdr", result.output.lower())
        self.assertFalse(self.output.exists())

    def test_reports_invalid_mapping_json_without_creating_output(self):
        self.mapping_path.write_text("[]", encoding="utf-8")

        result = CliRunner().invoke(catalog_main, self.arguments())

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("mapping", result.output.lower())
        self.assertFalse(self.output.exists())

    def test_rejects_an_identity_source_changed_after_mapping(self):
        self.references.write_text(REFERENCES + "\n", encoding="utf-8")

        result = CliRunner().invoke(catalog_main, self.arguments())

        self.assertNotEqual(result.exit_code, 0)
        self.assertRegex(result.output.lower(), r"checksum|source|changed")
        self.assertFalse(self.output.exists())

    def test_applies_the_selected_source_proteins_without_dropping_audit_rows(self):
        self.prepare_inputs([
            psm(),
            psm(Protein2="customA", PepSeq2="RTAK", PepPos2="5", PeptideLength2="4",
                LinkPos2="4", ProteinLinkPos2="8"),
        ])
        self.run_command(self.arguments() + ["--source-protein", "customA"])

        with (self.output / "claudio_input.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["entry1"], rows[0]["entry2"]), ("P12345", "P12345"))
        audit = json.loads((self.output / "psm_audit.json").read_text())
        self.assertEqual(len(audit), 2)
        self.assertIn("outside_selected_source_proteins", audit[0]["input_reasons"])

class TestPlinkPSM(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        fixture = Path(__file__).resolve().parents[1] / "fixtures" / "plink_pxd079671"
        manifest = json.loads((fixture / "manifest.json").read_text())
        for name, source in manifest["files"].items():
            raw = (fixture / name).read_bytes()
            self.assertEqual(hashlib.sha256(raw).hexdigest(), source["sha256"])
            (self.root / name).write_bytes(raw)
        self.psms = self.root / "spectra.csv"
        self.fasta = self.root / "search.fasta"
        with self.psms.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            self.columns = reader.fieldnames
            self.original_rows = list(reader)

    def read(self, path=None, **options):
        return xifdr.read_plink_psms(
            self.psms if path is None else path, self.fasta,
            dataset_accession="PXD079671", **options,
        )

    def toy_row(self, **changes):
        row = dict(self.original_rows[0])
        row.update(Peptide="ACK(3)-KPEP(3)", Proteins="customA(4)-customB(5)/",
                   Modifications="null", ReportedPSMFDR="0.005")
        row.update(changes)
        return row

    def write_rows(self, rows, columns=None):
        columns = self.columns + ["ReportedPSMFDR"] if columns is None else columns
        with self.psms.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    def test_real_psms_keep_their_titles_and_missing_confidence_information(self):
        with patch.object(pyxlms_parser, "parse_sites_from_plink", wraps=pyxlms_parser.parse_sites_from_plink) as parse_sites:
            bundle = self.read()

        self.assertEqual(parse_sites.call_count, 3)
        self.assertEqual(bundle["format"], "plink_psm")
        self.assertEqual(bundle["parser"]["name"], "pyXLMS")
        self.assertEqual(bundle["parser"]["version"], version("pyXLMS"))
        self.assertFalse(bundle["parser"]["modification_masses_assessed"])
        self.assertEqual(bundle["provenance"]["psm_sha256"], hashlib.sha256(self.psms.read_bytes()).hexdigest())
        self.assertEqual(bundle["source_columns"], self.columns)
        records = bundle["records"]
        self.assertEqual(len(records), 3)
        self.assertEqual([r["source_fields"] for r in records], self.original_rows)
        self.assertEqual([r["source_record_number"] for r in records], [1, 2, 3])
        for record in records:
            self.assertIsNone(record.get("source_scan"))
            self.assertEqual(record["decoy_status"], "unknown")
            self.assertTrue(all(value["status"] == "missing" and value["value"] is None
                                for value in record["fdr"].values()))
        prepared = prepare_af_input(bundle, {}, fdr_level="psm", fdr_cutoff=0.01)
        self.assertEqual(prepared["rows"], [])
        self.assertEqual(len(prepared["audit"]), 3)
        self.assertTrue(all("selected_fdr_missing" in row["input_reasons"] for row in prepared["audit"]))

    def test_peptide_and_position_matches_do_not_silently_rename_a_protein(self):
        record = self.read()["records"][0]
        first, second = record["endpoints"]["a"], record["endpoints"]["b"]

        self.assertEqual(first["source_protein_id"], "sp|P02787|TRFE_HUMAN")
        self.assertEqual((first["peptide"], first["link_position"], first["protein_position"]),
                         ("KDSGFQMNQLRGK", 1, 122))
        self.assertEqual(second["source_protein_id"], "BH2531|CThomas|Karst|")
        self.assertEqual((second["peptide"], second["peptide_start"], second["protein_position"]),
                         ("DEVSVQGHQK", 354, 354))
        self.assertIsNone(second["source_fasta_id"])
        self.assertEqual(second["status"], "not_evaluable")
        self.assertEqual(second["reason"], "protein_not_in_search_fasta")
        self.assertEqual(record["record_status"], "not_evaluable")

    def test_reported_modifications_are_retained_without_claiming_they_were_parsed(self):
        record = self.read()["records"][1]

        self.assertEqual(record["source_fields"]["Modifications"],
                         "Oxidation[M](5);Carbamidomethyl[C](15);Carbamidomethyl[C](16)")
        self.assertEqual(record["modification_status"], "not_parsed")
        self.assertIsNone(record["endpoints"]["a"]["modifications"])
        self.assertIsNone(record["endpoints"]["b"]["modifications"])
        self.assertEqual(record["endpoints"]["a"]["peptide"], "PVATMWESKGK")
        self.assertEqual(record["endpoints"]["b"]["peptide"], "CCTKPESER")

    def test_exports_validated_sites_with_explicit_fdr_and_identity_evidence(self):
        self.fasta.write_text(FASTA, encoding="utf-8")
        self.write_rows([self.toy_row()])
        bundle = self.read(decoy_prefix="REV_", fdr_columns={"psm": "ReportedPSMFDR"})
        record = bundle["records"][0]

        self.assertEqual(record["record_status"], "valid")
        self.assertEqual(record["decoy_status"], "target_target")
        self.assertEqual(record["fdr"]["psm"]["value"], 0.005)
        self.assertEqual(record["fdr"]["link"]["status"], "missing")
        references = self.root / "references.fasta"
        references.write_text(REFERENCES, encoding="utf-8")
        mappings = map_reference_proteins(bundle["proteins"], read_search_fasta(references), taxon_id=9031)
        output = write_claudio_input(
            bundle, mappings, self.root / "input", fdr_level="psm", fdr_cutoff=0.01,
            allowed_crosslinkers={"EDC-DE"}, allowed_residues={"K", "E"},
        )
        with (output / "claudio_input.csv").open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            self.assertEqual(reader.fieldnames[:8], CLAUDIO_COLUMNS)
            rows = list(reader)
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["entry1"], rows[0]["entry2"]), ("P12345", "Q99999"))
        self.assertEqual((rows[0]["position1"], rows[0]["position2"]), ("4", "5"))
        self.assertEqual((rows[0]["k_pos1"], rows[0]["k_pos2"]), ("3", "3"))
        self.assertEqual(rows[0]["source_crosslinker"], "EDC-DE")
        self.assertEqual(rows[0]["source_psm_id"], "1")
        self.assertEqual(rows[0]["source_spectrum_title"], self.original_rows[0]["Title"])
        self.assertEqual(rows[0]["source_scan"], "")
        manifest = json.loads((output / "manifest.json").read_text())
        self.assertEqual(manifest["kind"], "plink_claudio_input")
        self.assertEqual(manifest["parser"]["name"], "pyXLMS")
        self.assertFalse(manifest["assessment_performed"])
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)

    def test_duplicate_orders_keep_distinct_and_relocatable_source_records(self):
        self.fasta.write_text(FASTA, encoding="utf-8")
        row = self.toy_row()
        self.write_rows([row, row])
        first = self.read(decoy_prefix="REV_")
        relocated = self.root / "renamed.csv"
        relocated.write_bytes(self.psms.read_bytes())
        second = self.read(path=relocated, decoy_prefix="REV_")
        ids = [record["source_record_id"] for record in first["records"]]

        self.assertEqual(len(ids), 2)
        self.assertNotEqual(ids[0], ids[1])
        self.assertEqual(ids, [record["source_record_id"] for record in second["records"]])
        self.assertEqual([record["source_fields"]["Order"] for record in first["records"]], ["1", "1"])

    def test_protein_alternatives_do_not_create_unreported_pairings(self):
        self.fasta.write_text(FASTA + ">customC\nMACKRTAK\n>customD\nGGKPEPK\n", encoding="utf-8")
        self.write_rows([self.toy_row(Proteins="customA(4)-customB(5)/customC(4)-customD(5)/")])
        record = self.read(decoy_prefix="REV_")["records"][0]
        pairs = [(pair["a"]["source_protein_id"], pair["b"]["source_protein_id"])
                 for pair in record["protein_pairs"]]

        self.assertEqual(pairs, [("customA", "customB"), ("customC", "customD")])
        self.assertEqual(record["record_status"], "not_evaluable")
        self.assertIsNone(record["endpoints"]["a"]["source_fasta_id"])
        self.assertIsNone(record["endpoints"]["b"]["source_fasta_id"])

    def test_ragged_records_stay_in_the_audit_without_shifted_fields(self):
        with self.psms.open("a", encoding="utf-8", newline="") as handle:
            csv.writer(handle).writerow(["broken", "short"])
        records = self.read()["records"]

        self.assertEqual(len(records), 4)
        self.assertEqual(records[3]["record_reason"], "row_width_mismatch")
        self.assertEqual(records[3]["source_fields"], {})
        self.assertEqual(records[3]["source_values"], ["broken", "short"])
        self.assertEqual(records[0]["source_fields"], self.original_rows[0])

    def test_incorrect_protein_coordinates_are_not_realigned_automatically(self):
        self.fasta.write_text(FASTA, encoding="utf-8")
        self.write_rows([self.toy_row(Proteins="customA(5)-customB(5)/")])
        record = self.read(decoy_prefix="REV_")["records"][0]

        self.assertEqual(record["record_status"], "not_evaluable")
        self.assertEqual(record["endpoints"]["a"]["protein_position"], 5)
        self.assertEqual(record["endpoints"]["a"]["reason"], "peptide_sequence_mismatch")
        self.assertEqual(record["source_fields"]["Proteins"], "customA(5)-customB(5)/")

    def test_a_crosslink_level_table_is_not_imported_as_psms(self):
        columns = [column for column in self.columns if column != "Title"]
        self.write_rows([self.original_rows[0]], columns=columns)

        with self.assertRaisesRegex(ValueError, "PSM|spectra"):
            self.read()

class TestStandardizeFormatsCLI(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        fixture = Path(__file__).resolve().parents[1] / "fixtures" / "plink_pxd079671"
        self.psms = self.root / "renamed_results.csv"
        self.fasta = self.root / "search.fasta"
        self.psms.write_bytes((fixture / "spectra.csv").read_bytes())
        self.fasta.write_bytes((fixture / "search.fasta").read_bytes())
        with self.psms.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            self.columns = reader.fieldnames
            self.original_rows = list(reader)
        self.mapping = self.root / "mapping.json"
        self.mapping.write_text("{}", encoding="utf-8")
        self.output = self.root / "output"

    def arguments(self):
        return [
            "standardize", "--results", str(self.psms), "--fasta", str(self.fasta),
            "--dataset-accession", "PXD079671", "--output", str(self.output),
        ]

    def run_command(self, arguments):
        result = CliRunner().invoke(catalog_main, arguments)
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        return result

    def save_mapping(self):
        references = self.root / "references.fasta"
        references.write_text(REFERENCES, encoding="utf-8")
        mappings = map_reference_proteins(read_search_fasta(self.fasta), read_search_fasta(references), taxon_id=9031)
        self.mapping.write_text(json.dumps(mappings), encoding="utf-8")

    def prepare_plink_example(self):
        self.fasta.write_text(FASTA, encoding="utf-8")
        row = dict(self.original_rows[0])
        row.update(Peptide="ACK(3)-KPEP(3)", Proteins="customA(4)-customB(5)/",
                   Modifications="null", ReportedPSMFDR="0.005")
        with self.psms.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.columns + ["ReportedPSMFDR"])
            writer.writeheader()
            writer.writerow(row)
        self.save_mapping()

    def prepare_xifdr_example(self):
        self.fasta.write_text(FASTA, encoding="utf-8")
        columns = COLUMNS + ["exp charge"]
        with self.psms.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerow({**psm(), "exp charge": "3"})
        self.save_mapping()

    def plink_export_arguments(self):
        return self.arguments() + [
            "--mapping", str(self.mapping), "--fdr-level", "psm", "--fdr-cutoff", "0.01",
            "--fdr-column", "ReportedPSMFDR", "--decoy-prefix", "REV_",
            "--crosslinker", "EDC-DE", "--residue", "K", "--residue", "E",
        ]

    def test_real_plink_audit_requires_no_mapping_or_filter_policy(self):
        result = self.run_command(self.arguments() + ["--audit-only"])

        self.assertIn("Detected format: pLink PSM", result.output)
        self.assertIn("Source PSM records: 3", result.output)
        report = json.loads((self.output / "import_report.json").read_text())
        self.assertEqual(report["format"], "plink_psm")
        self.assertEqual([r["source_fields"] for r in report["records"]], self.original_rows)
        self.assertTrue(all(r["fdr"]["psm"]["status"] == "missing" for r in report["records"]))
        self.assertFalse((self.output / "claudio_input.csv").exists())
        manifest = json.loads((self.output / "manifest.json").read_text())
        self.assertEqual(manifest["kind"], "xlms_psm_import")
        self.assertEqual(manifest["format"], "plink_psm")
        self.assertEqual(manifest["counts"]["source_records"], 3)
        self.assertEqual(manifest["counts"]["unresolved_source_records"], 3)
        self.assertFalse(manifest["fdr_filter_applied"])
        self.assertFalse(manifest["identity_mapping_applied"])
        self.assertFalse(manifest["input_exported"])
        self.assertFalse(manifest["assessment_performed"])
        self.assertIsNone(manifest.get("fdr_policy"))
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((self.output / name).read_bytes()).hexdigest(), digest)

    def test_xifdr_audit_preserves_reported_fdr_without_selecting_a_level(self):
        self.prepare_xifdr_example()
        result = self.run_command(self.arguments() + ["--audit-only"])

        self.assertIn("Detected format: xiFDR PSM", result.output)
        report = json.loads((self.output / "import_report.json").read_text())
        self.assertEqual(report["format"], "xifdr_psm")
        self.assertEqual(report["parser"]["name"], "pyXLMS")
        self.assertEqual(len(report["records"]), 1)
        self.assertEqual(report["records"][0]["fdr"]["psm"]["value"], 0.02)
        self.assertEqual(report["records"][0]["fdr"]["link"]["value"], 0.005)
        manifest = json.loads((self.output / "manifest.json").read_text())
        self.assertFalse(manifest["fdr_filter_applied"])
        self.assertFalse(manifest["input_exported"])
        self.assertFalse((self.output / "claudio_input.csv").exists())

    def test_plink_export_uses_explicit_fdr_and_saved_identities(self):
        self.prepare_plink_example()
        result = self.run_command(self.plink_export_arguments())

        self.assertIn("Detected format: pLink PSM", result.output)
        self.assertIn("Standardized input rows: 1", result.output)
        with (self.output / "claudio_input.csv").open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            self.assertEqual(reader.fieldnames[:8], CLAUDIO_COLUMNS)
            rows = list(reader)
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["entry1"], rows[0]["entry2"]), ("P12345", "Q99999"))
        self.assertEqual((rows[0]["position1"], rows[0]["position2"]), ("4", "5"))
        self.assertEqual(rows[0]["source_spectrum_title"], self.original_rows[0]["Title"])
        manifest = json.loads((self.output / "manifest.json").read_text())
        self.assertEqual(manifest["kind"], "plink_claudio_input")
        self.assertEqual(manifest["parser"]["fdr_columns"], {"psm": "ReportedPSMFDR"})
        self.assertEqual(manifest["parser"]["decoy_prefix"], "REV_")
        self.assertEqual(manifest["fdr_policy"]["level"], "psm")
        self.assertFalse(manifest["assessment_performed"])

    def test_plink_export_requires_a_reported_fdr_column_and_decoy_prefix(self):
        self.prepare_plink_example()
        for flag in ("--fdr-column", "--decoy-prefix"):
            with self.subTest(missing=flag):
                arguments = self.plink_export_arguments()
                index = arguments.index(flag)
                del arguments[index:index + 2]
                result = CliRunner().invoke(catalog_main, arguments)
                self.assertEqual(result.exit_code, 2)
                self.assertIn(flag, result.output)
                self.assertFalse(self.output.exists())

    def test_audit_rejects_selection_options_instead_of_ignoring_them(self):
        for flag, value in (("--fdr-cutoff", "0.01"), ("--mapping", str(self.mapping)),
                            ("--source-protein", "bh2531")):
            with self.subTest(option=flag):
                result = CliRunner().invoke(catalog_main, self.arguments() + ["--audit-only", flag, value])
                self.assertEqual(result.exit_code, 2)
                self.assertIn("--audit-only", result.output)
                self.assertIn(flag, result.output)
                self.assertFalse(self.output.exists())

    def test_audit_keeps_existing_outputs_untouched(self):
        self.output.mkdir()
        marker = self.output / "keep.txt"
        marker.write_text("previous output", encoding="utf-8")

        result = CliRunner().invoke(catalog_main, self.arguments() + ["--audit-only"])

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("exist", result.output.lower())
        self.assertEqual(marker.read_text(), "previous output")
        self.assertEqual(list(self.output.iterdir()), [marker])

    def test_unknown_table_format_is_rejected_without_creating_output(self):
        self.psms.write_text("ID,Value\n1,2\n", encoding="utf-8")

        result = CliRunner().invoke(catalog_main, self.arguments() + ["--audit-only"])

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("format", result.output.lower())
        self.assertFalse(self.output.exists())

    def test_report_writer_checks_source_hashes_before_creating_output(self):
        bundle = xifdr.read_plink_psms(self.psms, self.fasta, dataset_accession="PXD079671")
        self.psms.write_bytes(self.psms.read_bytes() + b"\n")

        with self.assertRaisesRegex(ValueError, "source|checksum|changed"):
            xifdr.write_import_report(bundle, self.output)

        self.assertFalse(self.output.exists())

    def test_xifdr_export_rejects_plink_interpretation_options(self):
        self.prepare_xifdr_example()
        arguments = self.arguments() + [
            "--mapping", str(self.mapping), "--fdr-level", "link", "--fdr-cutoff", "0.01",
            "--crosslinker", "BS3", "--residue", "K",
        ]
        for flag, value in (("--decoy-prefix", "REV_"), ("--fdr-column", "LinkFDR")):
            with self.subTest(option=flag):
                result = CliRunner().invoke(catalog_main, arguments + [flag, value])
                self.assertEqual(result.exit_code, 2)
                self.assertIn("pLink", result.output)
                self.assertFalse(self.output.exists())

    def test_plink_audit_can_label_a_declared_fdr_field_without_filtering(self):
        self.prepare_plink_example()
        self.run_command(self.arguments() + [
            "--audit-only", "--decoy-prefix", "REV_",
            "--fdr-level", "psm", "--fdr-column", "ReportedPSMFDR",
        ])

        report = json.loads((self.output / "import_report.json").read_text())
        self.assertEqual(report["records"][0]["fdr"]["psm"]["value"], 0.005)
        self.assertEqual(report["parser"]["fdr_columns"], {"psm": "ReportedPSMFDR"})
        manifest = json.loads((self.output / "manifest.json").read_text())
        self.assertFalse(manifest["fdr_filter_applied"])
        self.assertIsNone(manifest.get("fdr_policy"))
        self.assertFalse((self.output / "claudio_input.csv").exists())


class TestPlinkNativeEvidence(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        fixture = Path(__file__).resolve().parents[1] / "fixtures" / "plink_pxd079671"
        self.psms = self.root / "spectra.csv"
        self.fasta = self.root / "search.fasta"
        self.native = self.root / "unfiltered.csv"
        self.parameters = self.root / "search.plink"
        for destination, relative in (
            (self.psms, "spectra.csv"), (self.fasta, "search.fasta"),
            (self.native, "native/unfiltered.csv"),
            (self.parameters, "native/parameters.plink"),
        ):
            destination.write_bytes((fixture / relative).read_bytes())
        checksums = json.loads((fixture / "native/source_checksums.json").read_text())
        for name, path in (("BH2531_2024.10.15.csv", self.native),
                           ("pLink_task_2024.10.15.11.45.41_G1_3A.plink", self.parameters)):
            expected = next(item for item in checksums if item["file_name"] == name)
            self.assertEqual(path.stat().st_size, expected["bytes"])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), expected["sha256"])
        with self.native.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            self.columns, self.native_rows = reader.fieldnames, list(reader)
        with self.psms.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            self.source_columns, self.source_rows = reader.fieldnames, list(reader)
        self.output = self.root / "output"

    def read(self, **options):
        return xifdr.read_plink_psms(
            self.psms, self.fasta, dataset_accession="PXD079671",
            unfiltered_path=self.native, parameters_path=self.parameters, **options,
        )

    def save_native(self, rows=None, columns=None):
        with self.native.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns or self.columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(self.native_rows if rows is None else rows)

    def change_first_match(self, **changes):
        rows = [dict(row) for row in self.native_rows]
        next(row for row in rows if row["Order"] == "1194").update(changes)
        self.save_native(rows)

    def arguments(self):
        return [
            "standardize", "--results", str(self.psms), "--fasta", str(self.fasta),
            "--dataset-accession", "PXD079671", "--plink-unfiltered", str(self.native),
            "--plink-parameters", str(self.parameters), "--output", str(self.output),
        ]

    def prepare_complete_example(self):
        self.fasta.write_text(FASTA, encoding="utf-8")
        source = dict(self.source_rows[0])
        source.update(Peptide="ACK(3)-KPEP(3)", Proteins="customA(4)-customB(5)/", Modifications="null")
        with self.psms.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.source_columns)
            writer.writeheader()
            writer.writerow(source)
        native = dict(self.native_rows[1193])
        native.update({key: source[key] for key in (
            "Title", "Charge", "Peptide", "Proteins", "FileID", "Modifications"
        )})
        native.update({"Target_Decoy": "2", "Q-value": "0.005", "Peptide_Type": "3"})
        self.save_native([native])
        references = self.root / "references.fasta"
        references.write_text(REFERENCES, encoding="utf-8")
        mappings = map_reference_proteins(read_search_fasta(self.fasta), read_search_fasta(references), taxon_id=9031)
        mapping = self.root / "mapping.json"
        mapping.write_text(json.dumps(mappings), encoding="utf-8")
        return mappings, mapping

    def test_real_matches_add_confidence_without_accepting_missing_sequences(self):
        before = xifdr.read_plink_psms(self.psms, self.fasta, dataset_accession="PXD079671")
        bundle = self.read()

        self.assertEqual([r["source_record_id"] for r in bundle["records"]],
                         [r["source_record_id"] for r in before["records"]])
        self.assertEqual([r["source_fields"] for r in bundle["records"]], self.source_rows)
        self.assertEqual([r["endpoints"] for r in bundle["records"]],
                         [r["endpoints"] for r in before["records"]])
        self.assertEqual(bundle["native_filter_policy"]["cutoff"], 0.05)
        self.assertEqual(bundle["native_filter_policy"]["level"], "psm")
        self.assertEqual(bundle["native_filter_policy"]["grouping_raw"], "EFT_INTRA_INTER_SEPARATE_FDR")
        self.assertEqual(bundle["native_filter_policy"]["software_version"], "2.3.11")
        self.assertEqual(bundle["native_filter_policy"]["crosslinkers"], ["EDC-DE"])
        self.assertEqual([r["native_evidence"]["native_order"] for r in bundle["records"]],
                         ["1194", "1205", "1206"])
        for record in bundle["records"]:
            self.assertEqual(record["native_evidence"]["status"], "verified")
            self.assertTrue(record["native_evidence"]["within_reported_fdr_cutoff"])
            self.assertEqual(record["decoy_status"], "target_target")
            self.assertEqual(record["fdr"]["psm"], {"raw": "0.000000", "value": 0.0, "status": "valid"})
            self.assertEqual(record["fdr"]["link"]["status"], "missing")
            self.assertEqual(record["record_status"], "not_evaluable")
            self.assertEqual(record["record_reason"], "source_endpoint_unresolved")
        prepared = prepare_af_input(bundle, {}, fdr_level="psm", fdr_cutoff=0.05)
        self.assertEqual(prepared["rows"], [])
        self.assertEqual(len(prepared["audit"]), 3)

    def test_native_table_order_does_not_determine_source_associations(self):
        first = self.read()
        self.save_native(list(reversed(self.native_rows)))
        second = self.read()
        self.assertEqual([r["source_record_id"] for r in first["records"]],
                         [r["source_record_id"] for r in second["records"]])
        self.assertEqual([r["native_evidence"]["native_order"] for r in second["records"]],
                         ["1194", "1205", "1206"])

    def test_the_same_title_with_other_sites_is_not_a_match(self):
        extra = dict(self.native_rows[1193])
        extra.update(Order="other", Proteins="differentA(122)-differentB(354)/")
        self.save_native([extra] + self.native_rows)
        record = self.read()["records"][0]
        self.assertEqual(record["native_evidence"]["native_order"], "1194")
        self.assertEqual(record["native_evidence"]["status"], "verified")

    def test_duplicate_native_matches_stay_ambiguous(self):
        self.save_native(self.native_rows + [self.native_rows[1193]])
        records = self.read()["records"]
        self.assertEqual(records[0]["native_evidence"]["reason"], "ambiguous_native_psm")
        self.assertEqual(len(records[0]["native_evidence"]["candidate_record_numbers"]), 2)
        self.assertEqual(records[0]["decoy_status"], "unknown")
        self.assertEqual(records[0]["fdr"]["psm"]["status"], "missing")
        self.assertEqual(records[1]["native_evidence"]["status"], "verified")

    def test_a_missing_native_record_does_not_remove_the_source(self):
        self.save_native([r for r in self.native_rows if r["Order"] != "1194"])
        records = self.read()["records"]
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0]["native_evidence"]["reason"], "native_psm_not_found")
        self.assertEqual(records[0]["source_fields"], self.source_rows[0])
        self.assertEqual(records[0]["decoy_status"], "unknown")

    def test_missing_source_identifiers_are_not_completed_from_native_rows(self):
        rows = [dict(row) for row in self.source_rows]
        rows[0]["Title"] = ""
        with self.psms.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.source_columns)
            writer.writeheader()
            writer.writerows(rows)
        self.change_first_match(Title="")
        record = self.read()["records"][0]
        self.assertEqual(record["source_fields"]["Title"], "")
        self.assertEqual(record["native_evidence"]["reason"], "missing_source_join_fields")
        self.assertEqual(record["decoy_status"], "unknown")

    def test_native_q_values_are_not_rescaled_or_replaced(self):
        for raw, status in (("100.000000", "out_of_range"), ("NaN", "non_finite"),
                            ("", "missing"), ("bad", "invalid")):
            with self.subTest(raw=raw):
                self.change_first_match(**{"Q-value": raw})
                record = self.read()["records"][0]
                self.assertEqual(record["native_evidence"]["status"], "not_evaluable")
                self.assertEqual(record["native_evidence"]["q_value"]["raw"], raw)
                self.assertEqual(record["native_evidence"]["q_value"]["status"], status)
                self.assertIsNone(record["native_evidence"]["within_reported_fdr_cutoff"])
                self.assertEqual(record["fdr"]["psm"]["raw"], raw)

    def test_decoy_codes_are_read_without_guessing_a_prefix(self):
        for code, expected in (("0", "decoy_decoy"), ("1", "target_decoy"),
                               ("2", "target_target"), ("7", "unknown")):
            with self.subTest(code=code):
                self.change_first_match(Target_Decoy=code)
                record = self.read()["records"][0]
                self.assertEqual(record["decoy_status"], expected)
                self.assertEqual(record["native_evidence"]["source_fields"]["Target_Decoy"], code)
                if code == "7":
                    self.assertEqual(record["native_evidence"]["reason"], "unknown_native_decoy_code")

    def test_cutoff_is_reported_separately_from_record_retention(self):
        self.change_first_match(**{"Q-value": "0.08"})
        bundle = self.read()
        self.assertEqual(len(bundle["records"]), 3)
        self.assertEqual(bundle["native_filter_policy"]["cutoff"], 0.05)
        record = bundle["records"][0]
        self.assertEqual(record["fdr"]["psm"]["value"], 0.08)
        self.assertFalse(record["native_evidence"]["within_reported_fdr_cutoff"])

    def test_isFilterIn_alone_is_not_used_as_fdr_acceptance(self):
        self.change_first_match(**{"Q-value": "0.08", "isFilterIn": "1"})
        record = self.read()["records"][0]
        self.assertEqual(record["native_evidence"]["source_fields"]["isFilterIn"], "1")
        self.assertFalse(record["native_evidence"]["within_reported_fdr_cutoff"])

    def test_a_modification_conflict_is_not_joined_by_title(self):
        self.change_first_match(Modifications="Oxidation[M](1)")
        record = self.read()["records"][0]
        self.assertEqual(record["native_evidence"]["reason"], "native_psm_not_found")
        self.assertEqual(record["decoy_status"], "unknown")

    def test_conflicting_native_peptide_or_protein_types_do_not_supply_confidence(self):
        for changes, reason in (({"Peptide_Type": "0"}, "native_record_not_crosslinked"),
                                ({"Protein_Type": "1"}, "native_protein_type_conflict")):
            with self.subTest(changes=changes):
                self.change_first_match(**changes)
                record = self.read()["records"][0]
                self.assertEqual(record["native_evidence"]["reason"], reason)
                self.assertEqual(record["decoy_status"], "unknown")
                self.assertEqual(record["fdr"]["psm"]["status"], "missing")

    def test_parameters_for_another_crosslinker_do_not_supply_confidence(self):
        self.parameters.write_text(self.parameters.read_text().replace("EDC-DE", "DSS"), encoding="utf-8")
        record = self.read()["records"][0]
        self.assertEqual(record["native_evidence"]["reason"], "crosslinker_conflict")
        self.assertEqual(record["decoy_status"], "unknown")
        self.assertEqual(record["fdr"]["psm"]["status"], "missing")

    def test_invalid_parameter_fdr_is_not_replaced_with_a_default(self):
        original = self.parameters.read_text()
        for value in ("NaN", "5", ""):
            with self.subTest(value=value):
                self.parameters.write_text(original.replace("FDR = 0.05", "FDR = " + value), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "FDR"):
                    self.read()

    def test_native_headers_must_be_complete_and_unique(self):
        for columns in ([c for c in self.columns if c != "Target_Decoy"], self.columns + ["Q-value"]):
            with self.subTest(columns=columns):
                self.save_native(columns=columns)
                with self.assertRaisesRegex(ValueError, "native|unfiltered|duplicate"):
                    self.read()

    def test_native_files_are_checked_again_before_writing_the_audit(self):
        for changed in (self.native, self.parameters):
            with self.subTest(path=changed.name):
                bundle = self.read()
                original = changed.read_bytes()
                changed.write_bytes(original + b"\n")
                with self.assertRaisesRegex(ValueError, "source|checksum|changed"):
                    xifdr.write_import_report(bundle, self.output)
                self.assertFalse(self.output.exists())
                changed.write_bytes(original)

    def test_native_files_must_be_paired_and_not_mixed_with_manual_policies(self):
        for options in (
            {"unfiltered_path": self.native}, {"parameters_path": self.parameters},
            {"unfiltered_path": self.native, "parameters_path": self.parameters, "decoy_prefix": "REV_"},
            {"unfiltered_path": self.native, "parameters_path": self.parameters, "fdr_columns": {"psm": "Score"}},
        ):
            with self.subTest(options=options):
                with self.assertRaisesRegex(ValueError, "unfiltered|parameters|native"):
                    xifdr.read_plink_psms(self.psms, self.fasta, dataset_accession="PXD079671", **options)

    def test_complete_sequences_can_be_exported_using_native_evidence(self):
        mappings, _ = self.prepare_complete_example()
        bundle = self.read()
        self.assertEqual(bundle["records"][0]["record_status"], "valid")
        output = write_claudio_input(
            bundle, mappings, self.output, fdr_level="psm", fdr_cutoff=0.01,
            allowed_crosslinkers={"EDC-DE"}, allowed_residues={"K", "E"},
        )
        with (output / "claudio_input.csv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["entry1"], rows[0]["entry2"]), ("P12345", "Q99999"))
        manifest = json.loads((output / "manifest.json").read_text())
        self.assertEqual(manifest["native_filter_policy"]["cutoff"], 0.05)
        self.assertEqual(manifest["fdr_policy"]["cutoff"], 0.01)

    def test_export_rejects_native_evidence_changed_after_import(self):
        mappings, _ = self.prepare_complete_example()
        bundle = self.read()
        self.native.write_bytes(self.native.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "native.*source.*changed"):
            write_claudio_input(
                bundle, mappings, self.output, fdr_level="psm", fdr_cutoff=0.01,
                allowed_crosslinkers={"EDC-DE"}, allowed_residues={"K", "E"},
            )
        self.assertFalse(self.output.exists())

    def test_cli_audit_reports_verified_confidence_but_exports_no_input(self):
        result = CliRunner().invoke(catalog_main, self.arguments() + ["--audit-only"])
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertIn("Native confidence verified: 3", result.output)
        report = json.loads((self.output / "import_report.json").read_text())
        self.assertEqual(report["native_filter_policy"]["cutoff"], 0.05)
        manifest = json.loads((self.output / "manifest.json").read_text())
        self.assertEqual(manifest["native_evidence_counts"], {"verified": 3})
        self.assertEqual(manifest["counts"]["unresolved_source_records"], 3)
        for flag in ("fdr_filter_applied", "identity_mapping_applied", "input_exported", "assessment_performed"):
            self.assertFalse(manifest[flag])
        self.assertFalse((self.output / "claudio_input.csv").exists())
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((self.output / name).read_bytes()).hexdigest(), digest)

    def test_cli_requires_both_native_files_and_rejects_manual_confidence_flags(self):
        for option in ("--plink-unfiltered", "--plink-parameters"):
            arguments = self.arguments() + ["--audit-only"]
            index = arguments.index(option)
            del arguments[index:index + 2]
            result = CliRunner().invoke(catalog_main, arguments)
            self.assertEqual(result.exit_code, 2)
            self.assertIn("together", result.output)
            self.assertFalse(self.output.exists())
        for extra in (["--decoy-prefix", "REV_"], ["--fdr-column", "Q-value", "--fdr-level", "psm"]):
            result = CliRunner().invoke(catalog_main, self.arguments() + ["--audit-only"] + extra)
            self.assertEqual(result.exit_code, 2)
            self.assertIn("native", result.output.lower())
            self.assertFalse(self.output.exists())

    def test_cli_native_export_uses_the_psm_level_and_requires_an_explicit_cutoff(self):
        _, mapping = self.prepare_complete_example()
        selection = ["--mapping", str(mapping), "--fdr-level", "psm", "--fdr-cutoff", "0.01",
                     "--crosslinker", "EDC-DE", "--residue", "K", "--residue", "E"]
        result = CliRunner().invoke(catalog_main, self.arguments() + selection)
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertIn("Standardized input rows: 1", result.output)
        other = self.root / "other"
        for level, cutoff in (("link", "0.01"), ("psm", None)):
            arguments = self.arguments() + selection
            arguments[arguments.index(str(self.output))] = str(other)
            arguments[arguments.index("--fdr-level") + 1] = level
            if cutoff is None:
                index = arguments.index("--fdr-cutoff")
                del arguments[index:index + 2]
            result = CliRunner().invoke(catalog_main, arguments)
            self.assertEqual(result.exit_code, 2)
            self.assertFalse(other.exists())

    def test_cli_rejects_native_plink_evidence_for_xifdr_inputs(self):
        with self.psms.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS + ["exp charge"])
            writer.writeheader()
            writer.writerow({**psm(), "exp charge": "3"})
        result = CliRunner().invoke(catalog_main, self.arguments() + ["--audit-only"])
        self.assertEqual(result.exit_code, 2)
        self.assertIn("pLink", result.output)
        self.assertFalse(self.output.exists())


class TestPlinkTitleMetadata(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        fixture = Path(__file__).resolve().parents[1] / "fixtures" / "plink_pxd014520"
        manifest = json.loads((fixture / "manifest.json").read_text())
        for name, source in manifest["files"].items():
            raw = (fixture / name).read_bytes()
            self.assertEqual(len(raw), source["bytes"])
            self.assertEqual(hashlib.sha256(raw).hexdigest(), source["sha256"])
            (self.root / name).write_bytes(raw)
        self.psms = self.root / "spectra.csv"
        self.fasta = self.root / "search.fasta"
        with self.psms.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            self.columns, self.raw_rows = next(reader), list(reader)
        self.first = dict(zip(
            self.columns,
            [self.raw_rows[0][0], ",".join(self.raw_rows[0][1:3])] + self.raw_rows[0][3:],
        ))

    def read(self, **options):
        return xifdr.read_plink_psms(
            self.psms, self.fasta, dataset_accession="PXD014520", **options,
        )

    def save_raw_rows(self, rows, columns=None):
        with self.psms.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(self.columns if columns is None else columns)
            writer.writerows(rows)

    def test_real_rows_keep_original_ids_and_cells_while_restoring_the_title(self):
        original = self.psms.read_bytes()
        bundle = self.read()
        records = bundle["records"]
        self.assertEqual(len(records), 257)
        self.assertEqual(bundle["source_columns"], self.columns)
        self.assertEqual([r["source_values"] for r in records], self.raw_rows)
        self.assertEqual([r["source_record_number"] for r in records], list(range(1, 258)))
        digest = hashlib.sha256(original).hexdigest()
        self.assertEqual(records[0]["source_record_id"], f"PXD014520:{digest}:000001")
        self.assertEqual(records[-1]["source_record_id"], f"PXD014520:{digest}:000257")
        self.assertEqual(records[0]["source_fields"], self.first)
        self.assertEqual(records[0]["source_fields"]["Charge"], "3")
        self.assertEqual(records[0]["source_fields"]["Peptide"], "SVLISLKQAPLVH(7)-LCDFNPKSSK(7)")
        self.assertEqual(records[0]["source_fields"]["Linker"], "BS3")
        self.assertEqual(records[0]["source_fields"]["Protein_Type"], "Intra-Protein")
        for record in records:
            interpretation = record["source_csv_interpretation"]
            self.assertEqual(interpretation["method"], "plink_unquoted_title_metadata")
            self.assertEqual(interpretation["original_field_count"], 22)
            self.assertEqual(interpretation["interpreted_field_count"], 21)
        self.assertEqual(self.psms.read_bytes(), original)

    def test_all_reported_endpoint_alternatives_match_the_search_sequences(self):
        records = self.read()["records"]
        pairs = [pair for record in records for pair in record["protein_pairs"]]
        self.assertEqual(len(pairs), 260)
        endpoints = [pair[side] for pair in pairs for side in ("a", "b")]
        self.assertEqual(len(endpoints), 520)
        self.assertTrue(all(endpoint["status"] == "valid" for endpoint in endpoints))
        self.assertEqual(records[0]["endpoints"]["a"]["source_fasta_id"], "sp|P09496|CLCA_HUMAN")
        self.assertEqual(records[0]["endpoints"]["a"]["protein_position"], 242)
        self.assertEqual(records[0]["endpoints"]["b"]["protein_position"], 223)
        for number in (166, 187, 217):
            record = records[number - 1]
            self.assertEqual(len(record["protein_pairs"]), 2)
            self.assertEqual(record["endpoints"]["b"]["reason"], "protein_identity_ambiguous")

    def test_format_recovery_does_not_invent_fdr_decoy_or_scan_values(self):
        records = self.read()["records"]
        self.assertTrue(all(record["decoy_status"] == "unknown" for record in records))
        self.assertTrue(all(record["fdr"]["psm"]["status"] == "missing" for record in records))
        self.assertTrue(all(record["source_scan"] is None for record in records))
        self.assertTrue(all(record["record_status"] == "not_evaluable" for record in records))

    def test_a_correctly_quoted_title_is_read_without_a_recovery_marker(self):
        self.save_raw_rows([[self.first[name] for name in self.columns]])
        record = self.read()["records"][0]
        self.assertEqual(record["source_fields"], self.first)
        self.assertEqual(len(record["source_values"]), 21)
        self.assertNotIn("source_csv_interpretation", record)
        self.assertEqual(record["endpoints"]["a"]["status"], "valid")

    def test_an_extra_cell_outside_the_title_stays_unresolved(self):
        row = list(self.raw_rows[0])
        row.insert(5, "unexplained cell")
        self.save_raw_rows([row])
        record = self.read()["records"][0]
        self.assertEqual(record["source_values"], row)
        self.assertEqual(record["source_fields"], {})
        self.assertEqual(record["record_reason"], "row_width_mismatch")
        self.assertNotIn("source_csv_interpretation", record)

    def test_arbitrary_extra_title_cells_are_not_merged(self):
        row = list(self.raw_rows[0])
        row[2] = "some unrelated annotation"
        self.save_raw_rows([row])
        record = self.read()["records"][0]
        self.assertEqual(record["record_reason"], "row_width_mismatch")
        self.assertEqual(record["source_fields"], {})

    def test_invalid_alignment_anchors_prevent_title_recovery(self):
        for column, value in (("Charge", "bad"), ("Score", "NaN"),
                              ("Peptide_Type", "unknown"), ("FileID", "0")):
            with self.subTest(column=column):
                row = list(self.raw_rows[0])
                row[self.columns.index(column) + 1] = value
                self.save_raw_rows([row])
                record = self.read()["records"][0]
                self.assertEqual(record["record_reason"], "row_width_mismatch")
                self.assertEqual(record["source_fields"], {})

    def test_a_valid_title_does_not_fix_inconsistent_protein_positions(self):
        row = list(self.raw_rows[0])
        index = self.columns.index("Proteins") + 1
        row[index] = row[index].replace("(242)", "(243)")
        self.save_raw_rows([row])
        record = self.read()["records"][0]
        self.assertIn("source_csv_interpretation", record)
        self.assertEqual(record["endpoints"]["a"]["protein_position"], 243)
        self.assertEqual(record["endpoints"]["a"]["reason"], "peptide_sequence_mismatch")
        self.assertEqual(record["record_reason"], "source_endpoint_unresolved")

    def test_additional_reported_fdr_columns_keep_their_original_values(self):
        self.save_raw_rows([self.raw_rows[0] + ["0.02"]], self.columns + ["ReportedPSMFDR"])
        record = self.read(fdr_columns={"psm": "ReportedPSMFDR"})["records"][0]
        self.assertEqual(record["source_fields"]["ReportedPSMFDR"], "0.02")
        self.assertEqual(record["fdr"]["psm"]["value"], 0.02)
        self.assertEqual(record["source_fields"]["Linker"], "BS3")

    def test_a_bad_row_does_not_prevent_other_psms_from_being_imported(self):
        self.save_raw_rows([self.raw_rows[0], ["broken", "short"], self.raw_rows[1]])
        records = self.read()["records"]
        self.assertEqual(len(records), 3)
        self.assertEqual(records[1]["source_values"], ["broken", "short"])
        self.assertEqual(records[1]["record_reason"], "row_width_mismatch")
        self.assertEqual(records[0]["endpoints"]["a"]["status"], "valid")
        self.assertEqual(records[2]["endpoints"]["a"]["status"], "valid")

    def test_cli_audit_keeps_all_sources_and_their_interpretation(self):
        output = self.root / "audit"
        result = CliRunner().invoke(catalog_main, [
            "standardize", "--results", str(self.psms), "--fasta", str(self.fasta),
            "--dataset-accession", "PXD014520", "--audit-only", "--output", str(output),
        ])
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertIn("Source PSM records: 257", result.output)
        self.assertIn("Records with parser errors: 0", result.output)
        report = json.loads((output / "import_report.json").read_text())
        self.assertEqual(len(report["records"]), 257)
        self.assertEqual(report["records"][0]["source_values"], self.raw_rows[0])
        self.assertEqual(report["records"][0]["source_fields"], self.first)
        manifest = json.loads((output / "manifest.json").read_text())
        self.assertFalse(manifest["fdr_filter_applied"])
        self.assertFalse(manifest["input_exported"])
        self.assertFalse((output / "claudio_input.csv").exists())


class TestAuthorFilteredPSMs(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fixture = Path(__file__).resolve().parents[1] / "fixtures"
        self.psms = self.root / "spectra.csv"
        self.fasta = self.root / "search.fasta"
        self.policy_path = self.root / "author_filter.json"
        self.note = self.root / "source_note.json"
        self.note.write_text('{"scope": "Reviewed source filtering for this test"}\n', encoding="utf-8")
        source = self.fixture / "plink_pxd079671"
        self.psms.write_bytes((source / "spectra.csv").read_bytes())
        self.fasta.write_bytes((source / "search.fasta").read_bytes())
        with self.psms.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            self.columns, self.source_rows = reader.fieldnames, list(reader)
        self.output = self.root / "output"
        self.write_policy()

    def write_policy(self, **changes):
        policy = {
            "schema_version": 1, "kind": "xlms_author_filter_evidence", "review_status": "reviewed",
            "dataset_accession": "PXD079671",
            "psm_sha256": hashlib.sha256(self.psms.read_bytes()).hexdigest(),
            "fasta_sha256": hashlib.sha256(self.fasta.read_bytes()).hexdigest(),
            "result_scope": "author_filtered_psm_table", "target_only": True,
            "fasta_role": "target_search_database", "filter_level": "psm",
            "reported_fdr_cutoff": 0.05, "per_record_fdr_available": False,
            "sources": [{"path": self.note.name, "sha256": hashlib.sha256(self.note.read_bytes()).hexdigest(),
                         "role": "reviewed_source_note"}],
        }
        policy.update(changes)
        self.policy_path.write_text(json.dumps(policy), encoding="utf-8")
        return policy

    def read(self, **options):
        return xifdr.read_plink_psms(
            self.psms, self.fasta, dataset_accession="PXD079671",
            author_filter_path=self.policy_path, **options,
        )

    def prepare_complete_example(self, count=1):
        self.fasta.write_text(FASTA, encoding="utf-8")
        row = dict(self.source_rows[0])
        row.update(Peptide="ACK(3)-KPEP(3)", Proteins="customA(4)-customB(5)/", Modifications="null")
        with self.psms.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.columns)
            writer.writeheader()
            writer.writerows([row] * count)
        self.write_policy()
        reference = self.root / "references.fasta"
        reference.write_text(REFERENCES, encoding="utf-8")
        mappings = map_reference_proteins(read_search_fasta(self.fasta), read_search_fasta(reference), taxon_id=9031)
        mapping = self.root / "mapping.json"
        mapping.write_text(json.dumps(mappings), encoding="utf-8")
        return mappings, mapping

    def export(self, bundle, mappings):
        return write_claudio_input(
            bundle, mappings, self.output,
            allowed_crosslinkers={"EDC-DE"}, allowed_residues={"K", "E"},
        )

    def arguments(self):
        return [
            "standardize", "--results", str(self.psms), "--fasta", str(self.fasta),
            "--dataset-accession", "PXD079671", "--author-filter-evidence", str(self.policy_path),
            "--output", str(self.output),
        ]

    def test_real_source_table_keeps_fdr_missing_and_retains_ambiguous_proteins(self):
        source = self.fixture / "plink_pxd014520"
        bundle = xifdr.read_plink_psms(
            source / "spectra.csv", source / "search.fasta", dataset_accession="PXD014520",
            author_filter_path=source / "author_filter/policy.json",
        )
        records = bundle["records"]
        self.assertEqual(len(records), 257)
        self.assertEqual(sum(r["record_status"] == "valid" for r in records), 254)
        self.assertTrue(all(r["source_filter_status"] == "reported_accepted" for r in records))
        self.assertTrue(all(r["decoy_status"] == "target_target" for r in records))
        self.assertTrue(all(r["fdr"]["psm"]["status"] == "missing" for r in records))
        self.assertEqual(bundle["author_filter_policy"]["reported_fdr_cutoff"], 0.05)
        for number in (166, 187, 217):
            self.assertEqual(records[number - 1]["record_status"], "not_evaluable")
            self.assertEqual(records[number - 1]["endpoints"]["b"]["reason"], "protein_identity_ambiguous")

    def test_export_records_the_author_cutoff_without_fabricating_per_psm_values(self):
        mappings, _ = self.prepare_complete_example()
        bundle = self.read()
        output = self.export(bundle, mappings)
        with (output / "claudio_input.csv").open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            self.assertEqual(reader.fieldnames[:8], CLAUDIO_COLUMNS)
            rows = list(reader)
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["entry1"], rows[0]["entry2"]), ("P12345", "Q99999"))
        self.assertEqual(rows[0]["selection_basis"], "author_filtered_psm_table")
        self.assertEqual(rows[0]["source_filter_level"], "psm")
        self.assertEqual(float(rows[0]["source_filter_cutoff"]), 0.05)
        for name in ("selected_fdr_level", "selected_fdr_value", "selected_fdr_cutoff"):
            self.assertEqual(rows[0][name], "")
        audit = json.loads((output / "psm_audit.json").read_text())
        self.assertIsNone(audit[0]["fdr"]["psm"]["value"])
        manifest = json.loads((output / "manifest.json").read_text())
        self.assertEqual(manifest["fdr_policy"]["kind"], "author_filtered_psm_table")
        self.assertFalse(manifest["fdr_policy"]["per_record_fdr_applied"])
        self.assertEqual(manifest["author_filter_policy"]["reported_fdr_cutoff"], 0.05)

    def test_unresolved_source_proteins_cannot_be_accepted_by_a_filter_declaration(self):
        bundle = self.read()
        self.assertTrue(all(r["record_status"] == "not_evaluable" for r in bundle["records"]))
        self.assertTrue(all(r["source_filter_status"] == "not_evaluable" for r in bundle["records"]))
        self.assertTrue(all(r["decoy_status"] == "unknown" for r in bundle["records"]))
        prepared = prepare_af_input(bundle, {})
        self.assertEqual(prepared["rows"], [])
        self.assertEqual(len(prepared["audit"]), 3)

    def test_exact_uniprot_sequence_identity_is_still_required(self):
        self.prepare_complete_example()
        bundle = self.read()
        prepared = prepare_af_input(bundle, {})
        self.assertEqual(prepared["rows"], [])
        self.assertIn("protein_identity_unresolved_a", prepared["audit"][0]["input_reasons"])

    def test_a_declared_filter_does_not_realign_incorrect_source_positions(self):
        self.prepare_complete_example()
        self.psms.write_text(self.psms.read_text().replace("customA(4)", "customA(5)"), encoding="utf-8")
        self.write_policy()
        record = self.read()["records"][0]
        self.assertEqual(record["record_status"], "not_evaluable")
        self.assertEqual(record["endpoints"]["a"]["reason"], "peptide_sequence_mismatch")

    def test_duplicate_psms_remain_distinct_and_are_not_collapsed_to_residue_pairs(self):
        mappings, _ = self.prepare_complete_example(count=2)
        bundle = self.read()
        self.assertNotEqual(bundle["records"][0]["source_record_id"], bundle["records"][1]["source_record_id"])
        prepared = prepare_af_input(bundle, mappings)
        self.assertEqual(len(prepared["rows"]), 2)
        self.assertEqual(len(prepared["audit"]), 2)

    def test_malformed_rows_are_not_accepted_by_a_table_level_declaration(self):
        mappings, _ = self.prepare_complete_example()
        with self.psms.open("a", encoding="utf-8", newline="") as handle:
            csv.writer(handle).writerow(["broken", "short"])
        self.write_policy()
        bundle = self.read()
        prepared = prepare_af_input(bundle, mappings)
        self.assertEqual(len(prepared["rows"]), 1)
        self.assertEqual(len(prepared["audit"]), 2)
        bad = prepared["audit"][1]
        self.assertEqual(bad["record_reason"], "row_width_mismatch")
        self.assertEqual(bad["source_filter_status"], "not_evaluable")
        self.assertIn("author_filter_acceptance_unresolved", bad["input_reasons"])

    def test_the_declaration_is_bound_to_dataset_and_source_checksums(self):
        for changes in ({"dataset_accession": "PXD000001"}, {"psm_sha256": "0" * 64},
                        {"fasta_sha256": "0" * 64}):
            with self.subTest(changes=changes):
                self.write_policy(**changes)
                with self.assertRaisesRegex(ValueError, "author.*(source|dataset)|checksum"):
                    self.read()

    def test_only_reviewed_target_only_psm_declarations_are_supported(self):
        for changes in ({"review_status": "pending"}, {"target_only": False},
                        {"target_only": "true"}, {"result_scope": "final_residue_pairs"},
                        {"fasta_role": "target_and_decoy_database"}, {"filter_level": "link"},
                        {"per_record_fdr_available": True}, {"kind": "other"}):
            with self.subTest(changes=changes):
                self.write_policy(**changes)
                with self.assertRaisesRegex(ValueError, "author|Author"):
                    self.read()

    def test_missing_or_invalid_author_cutoffs_are_not_replaced_with_defaults(self):
        for cutoff in (None, "0.05", True, float("nan"), -1, 5):
            with self.subTest(cutoff=cutoff):
                self.write_policy(reported_fdr_cutoff=cutoff)
                with self.assertRaises(ValueError):
                    self.read()

    def test_supporting_sources_must_exist_and_match_their_checksums(self):
        for sources in ([], [{"path": "missing.json", "sha256": "0" * 64, "role": "publication"}],
                        [{"path": self.note.name, "sha256": "0" * 64, "role": "publication"}],
                        [{"path": "../outside.json", "sha256": "0" * 64, "role": "publication"}]):
            with self.subTest(sources=sources):
                self.write_policy(sources=sources)
                with self.assertRaises((ValueError, OSError)):
                    self.read()

    def test_author_evidence_is_not_mixed_with_native_or_manual_fdr_policies(self):
        for options in ({"decoy_prefix": "REV_"}, {"fdr_columns": {"psm": "Score"}},
                        {"unfiltered_path": self.psms, "parameters_path": self.policy_path}):
            with self.subTest(options=options):
                with self.assertRaisesRegex(ValueError, "author|Author"):
                    self.read(**options)
        self.prepare_complete_example()
        bundle = self.read()
        with self.assertRaisesRegex(ValueError, "author|Author"):
            prepare_af_input(bundle, {}, fdr_level="psm", fdr_cutoff=0.01)

    def test_audit_snapshots_can_replay_the_same_author_declaration(self):
        bundle = self.read()
        output = xifdr.write_import_report(bundle, self.output)
        manifest = json.loads((output / "manifest.json").read_text())
        snapshot = output / manifest["author_filter_evidence_snapshot"]
        replay = xifdr.read_plink_psms(
            self.psms, self.fasta, dataset_accession="PXD079671", author_filter_path=snapshot,
        )
        self.assertEqual(replay["author_filter_policy"], bundle["author_filter_policy"])
        self.assertEqual(replay["records"], bundle["records"])
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)
        for flag in ("fdr_filter_applied", "identity_mapping_applied", "input_exported", "assessment_performed"):
            self.assertFalse(manifest[flag])

    def test_changed_sources_are_rejected_before_creating_outputs(self):
        for path in (self.policy_path, self.note, self.psms, self.fasta):
            with self.subTest(path=path.name):
                bundle = self.read()
                raw = path.read_bytes()
                path.write_bytes(raw + b"\n")
                with self.assertRaisesRegex(ValueError, "source|checksum|changed"):
                    xifdr.write_import_report(bundle, self.output)
                self.assertFalse(self.output.exists())
                path.write_bytes(raw)

    def test_export_snapshots_include_the_reviewed_evidence(self):
        mappings, _ = self.prepare_complete_example()
        output = self.export(self.read(), mappings)
        manifest = json.loads((output / "manifest.json").read_text())
        self.assertEqual((output / manifest["author_filter_evidence_snapshot"]).read_bytes(), self.policy_path.read_bytes())
        self.assertEqual((output / "author_filter_sources/source_note.json").read_bytes(), self.note.read_bytes())
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)

    def test_cli_exports_without_per_record_fdr_arguments(self):
        _, mapping = self.prepare_complete_example()
        result = CliRunner().invoke(catalog_main, self.arguments() + [
            "--mapping", str(mapping), "--crosslinker", "EDC-DE", "--residue", "K", "--residue", "E",
        ])
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertIn("Standardized input rows: 1", result.output)

    def test_cli_audit_reports_author_filtering_without_exporting_inputs(self):
        self.prepare_complete_example()
        result = CliRunner().invoke(catalog_main, self.arguments() + ["--audit-only"])
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertIn("Author-filtered PSM records: 1", result.output)
        report = json.loads((self.output / "import_report.json").read_text())
        self.assertEqual(report["records"][0]["fdr"]["psm"]["status"], "missing")
        self.assertFalse((self.output / "claudio_input.csv").exists())

    def test_cli_rejects_per_record_and_native_options_with_author_evidence(self):
        for extra in (["--fdr-level", "psm"], ["--fdr-cutoff", "0.01"],
                      ["--fdr-column", "Score"], ["--decoy-prefix", "REV_"],
                      ["--plink-unfiltered", str(self.psms), "--plink-parameters", str(self.policy_path)]):
            with self.subTest(extra=extra):
                result = CliRunner().invoke(catalog_main, self.arguments() + ["--audit-only"] + extra)
                self.assertEqual(result.exit_code, 2)
                self.assertIn("author", result.output.lower())
                self.assertFalse(self.output.exists())

    def test_cli_rejects_an_author_plink_declaration_for_xifdr_results(self):
        with self.psms.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS + ["exp charge"])
            writer.writeheader()
            writer.writerow({**psm(), "exp charge": "3"})
        result = CliRunner().invoke(catalog_main, self.arguments() + ["--audit-only"])
        self.assertEqual(result.exit_code, 2)
        self.assertIn("pLink", result.output)
        self.assertFalse(self.output.exists())

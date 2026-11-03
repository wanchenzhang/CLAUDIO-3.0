import csv
import hashlib
import json
import unittest
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, Mock, call, patch

import requests
from pyXLMS import parser as pyxlms_parser

from claudio.dataset_catalog.xifdr import (
    read_xifdr_psms, read_search_fasta, normalize_xifdr_peptide, target_proteins,
    map_reference_proteins, resolve_uniprot_proteins, prepare_af_input,
    write_xifdr_outputs, assess_af_complexes, compare_protein_sequences, search_uniprot_candidates,
    project_xifdr_candidates, write_projection_outputs, write_claudio_input,
)


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

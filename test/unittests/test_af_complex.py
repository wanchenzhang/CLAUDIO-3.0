import json
import unittest
from pathlib import Path

from claudio.utils.af_complex import read_complex_chains


FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "af_complex"
)


class TestReadComplexChains(unittest.TestCase):

    def test_preserves_both_chains_and_their_protein_identities(self):
        structure_path = FIXTURE_DIR / "AF-0000000208747029.cif"

        records = json.loads(
            (FIXTURE_DIR / "metadata.json").read_text(encoding="utf-8")
        )

        chains = read_complex_chains(structure_path, records)

        self.assertEqual(len(chains), 2)

        by_id = {chain["chain_id"]: chain for chain in chains}
        self.assertEqual(set(by_id), {"A", "B"})

        self.assertEqual(by_id["A"]["protein_id"], "A0A0P0Y832")
        self.assertEqual(by_id["A"]["residue_count"], 346)

        self.assertEqual(by_id["B"]["protein_id"], "P0C2Z6")
        self.assertEqual(by_id["B"]["residue_count"], 507)

        for record in records:
            chain = by_id[record["chainId"]]
            self.assertEqual(chain["sequence"], record["sequence"])
    def test_rejects_sequence_mismatch(self):
        structure_path = FIXTURE_DIR / "AF-0000000208747029.cif"

        records = json.loads(
            (FIXTURE_DIR / "metadata.json").read_text(encoding="utf-8")
        )

        for record in records:
            if record["chainId"] == "A":
                record["sequence"] = "A" * len(record["sequence"])

        with self.assertRaises(ValueError):
            read_complex_chains(structure_path, records)
    def test_rejects_missing_chain_metadata(self):
        structure_path = FIXTURE_DIR / "AF-0000000208747029.cif"

        records = json.loads(
            (FIXTURE_DIR / "metadata.json").read_text(encoding="utf-8")
        )

        # Keep only chain A's metadata, leaving chain B without matching metadata.
        records = [
            record for record in records
            if record["chainId"] == "A"
        ]

        with self.assertRaisesRegex(ValueError, r"chain B"):
            read_complex_chains(structure_path, records)


if __name__ == "__main__":
    unittest.main()

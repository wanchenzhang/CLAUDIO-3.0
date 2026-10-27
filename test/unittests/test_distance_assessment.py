"""Assess exact distances and strict lower bounds against an explicit interval."""
import unittest

import pandas as pd

from claudio.module04.src.algorithm.distance_assessment import (
    assess_distance_interval,
)


class TestDistanceAssessment(unittest.TestCase):
    def assess(self, exact=None, bound=None, boundary="inclusive"):
        return assess_distance_interval(
            exact, bound, minimum=5.0, maximum=35.0, boundary=boundary,
        )

    def test_classifies_exact_distances_inside_below_and_above_the_interval(self):
        cases = (
            (20.0, "satisfied", "exact_within_range"),
            (4.0, "violated", "exact_below_range"),
            (36.0, "violated", "exact_above_range"),
        )
        for distance, status, reason in cases:
            with self.subTest(distance=distance):
                result = self.assess(distance)
                self.assertEqual((result["status"], result["reason"]), (status, reason))

    def test_interval_endpoints_follow_the_explicit_boundary_policy(self):
        cases = (
            ("inclusive", 5.0, "satisfied", "exact_within_range"),
            ("inclusive", 35.0, "satisfied", "exact_within_range"),
            ("exclusive", 5.0, "violated", "exact_below_range"),
            ("exclusive", 35.0, "violated", "exact_above_range"),
        )
        for policy, distance, status, reason in cases:
            with self.subTest(policy=policy, distance=distance):
                result = self.assess(distance, boundary=policy)
                self.assertEqual((result["status"], result["reason"]), (status, reason))

    def test_strict_lower_bound_at_or_above_the_maximum_proves_a_violation(self):
        for policy in ("inclusive", "exclusive"):
            for bound in (35.0, 75.328):
                with self.subTest(policy=policy, bound=bound):
                    result = self.assess(bound=bound, boundary=policy)
                    self.assertEqual(result["status"], "violated")
                    self.assertEqual(result["reason"], "lower_bound_above_range")

    def test_a_lower_bound_below_the_maximum_cannot_prove_satisfaction(self):
        for bound in (0.0, 4.0, 5.0, 20.0, 34.999):
            with self.subTest(bound=bound):
                result = self.assess(bound=bound)
                self.assertEqual(result["status"], "indeterminate")
                self.assertEqual(result["reason"], "lower_bound_inconclusive")

    def test_absent_distances_are_not_evaluable(self):
        for missing in (None, float("nan"), pd.NA):
            with self.subTest(missing=missing):
                result = self.assess(missing, missing)
                self.assertEqual(result["status"], "not_evaluable")
                self.assertEqual(result["reason"], "distance_unavailable")

    def test_invalid_values_are_not_used_as_distances(self):
        for invalid in (-1, float("inf"), float("-inf"), True, "20", ">20"):
            for field in ("exact", "bound"):
                with self.subTest(value=invalid, field=field):
                    result = self.assess(**{field: invalid})
                    self.assertEqual(result["status"], "not_evaluable")
                    self.assertEqual(result["reason"], "invalid_distance_value")

    def test_simultaneous_exact_and_bounded_values_are_not_silently_combined(self):
        result = self.assess(20.0, 10.0)
        self.assertEqual(result["status"], "not_evaluable")
        self.assertEqual(result["reason"], "conflicting_distance_values")

    def test_records_thresholds_boundary_and_rule_version(self):
        result = self.assess(20.0, boundary="exclusive")
        self.assertEqual(result["minimum"], 5.0)
        self.assertEqual(result["maximum"], 35.0)
        self.assertEqual(result["boundary"], "exclusive")
        self.assertEqual(result["rule_version"], "distance-interval-v1")

    def test_rejects_invalid_interval_configuration(self):
        cases = (
            {"minimum": -1.0}, {"minimum": float("nan")},
            {"maximum": float("inf")}, {"minimum": True},
            {"minimum": "5"}, {"maximum": 5.0}, {"maximum": 4.0},
            {"boundary": "unknown"}, {"boundary": None},
        )
        for changes in cases:
            with self.subTest(changes=changes):
                options = {"minimum": 5.0, "maximum": 35.0, "boundary": "inclusive"}
                options.update(changes)
                with self.assertRaisesRegex(ValueError, "interval|boundary"):
                    assess_distance_interval(20.0, None, **options)

    def test_zero_distance_is_a_numeric_value_rather_than_missing(self):
        result = assess_distance_interval(
            0.0, None, minimum=0.0, maximum=35.0, boundary="inclusive",
        )
        self.assertEqual(result["status"], "satisfied")
        exclusive = assess_distance_interval(
            0.0, None, minimum=0.0, maximum=35.0, boundary="exclusive",
        )
        self.assertEqual(exclusive["status"], "violated")

    def test_boundary_policy_has_no_implicit_default(self):
        with self.assertRaises(TypeError):
            assess_distance_interval(20.0, None, minimum=5.0, maximum=35.0)


if __name__ == "__main__":
    unittest.main()

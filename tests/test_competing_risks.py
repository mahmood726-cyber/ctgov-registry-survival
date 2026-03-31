"""Tests for Aalen-Johansen cumulative incidence function estimator."""
import numpy as np
import pandas as pd
import pytest

from rse.cohort import build_cohort, classify_events
from rse.competing_risks import aalen_johansen, cif_at_times


class TestAalenJohansen:
    def test_cif_sum_leq_one(self, synthetic_cohort):
        """Sum of all CIFs must be <= 1.0 at every time point."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        aj = aalen_johansen(events, event_codes=[1, 2, 3])
        for t_idx in range(len(aj["times"])):
            total = sum(aj["cif"][code][t_idx] for code in [1, 2, 3])
            assert total <= 1.0 + 1e-10, f"CIF sum > 1 at t={aj['times'][t_idx]}: {total}"

    def test_cif_monotonically_increasing(self, synthetic_cohort):
        """Each CIF must be non-decreasing."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        aj = aalen_johansen(events, event_codes=[1, 2, 3])
        for code in [1, 2, 3]:
            cif = aj["cif"][code]
            for i in range(len(cif) - 1):
                assert cif[i] <= cif[i + 1] + 1e-10, \
                    f"CIF for event {code} decreases at index {i}"

    def test_known_two_event_scenario(self):
        """Hand-calculated Aalen-Johansen on a minimal dataset."""
        df = pd.DataFrame({
            "time_days": [1, 2, 2, 3, 3],
            "event_code": [1, 0, 1, 2, 1],
        })
        aj = aalen_johansen(df, event_codes=[1, 2])
        # Data: obs0(t=1,e=1), obs1(t=2,e=0), obs2(t=2,e=1), obs3(t=3,e=2), obs4(t=3,e=1)
        #
        # t=1: n=5, d1=1, d_total=1
        #   CIF1 += S(0)*d1/n = 1.0 * 1/5 = 0.2        -> CIF1=0.2, CIF2=0.0
        #   S(1) = 1.0 * (1 - 1/5) = 0.8
        #
        # t=2: n=4 (obs1..4 still at risk), d1=1, d_total=1, censored obs1 does NOT reduce d_total
        #   CIF1 += S(1)*d1/n = 0.8 * 1/4 = 0.2        -> CIF1=0.4, CIF2=0.0
        #   S(2) = 0.8 * (1 - 1/4) = 0.6
        #
        # t=3: n=2 (obs3,obs4), d1=1, d2=1, d_total=2
        #   CIF1 += S(2)*d1/n = 0.6 * 1/2 = 0.3        -> CIF1=0.7
        #   CIF2 += S(2)*d2/n = 0.6 * 1/2 = 0.3        -> CIF2=0.3
        #   S(3) = 0.6 * (1 - 2/2) = 0.0
        cif1_final = aj["cif"][1][-1]
        cif2_final = aj["cif"][2][-1]
        assert abs(cif1_final - 0.7) < 1e-10, f"CIF1 = {cif1_final}, expected 0.7"
        assert abs(cif2_final - 0.3) < 1e-10, f"CIF2 = {cif2_final}, expected 0.3"

    def test_cif_plus_survival_equals_one(self, synthetic_cohort):
        """At each time, sum(CIFs) + S(t) should approximately equal 1."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        aj = aalen_johansen(events, event_codes=[1, 2, 3])
        for t_idx in range(len(aj["times"])):
            total = sum(aj["cif"][code][t_idx] for code in [1, 2, 3])
            s = aj["survival"][t_idx]
            assert abs(total + s - 1.0) < 1e-6, \
                f"CIF + S != 1 at t={aj['times'][t_idx]}: {total} + {s} = {total + s}"

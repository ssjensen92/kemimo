import unittest
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]


class ThreePhaseEwtScalingTests(unittest.TestCase):
    def test_spin_and_nospin_templates_use_per_column_scaling(self):
        templates = (
            "kemimo_ode_include_H2.f90",
            "kemimo_ode_include_H2_nospin.f90",
        )
        for filename in templates:
            with self.subTest(filename=filename):
                source = (
                    ROOT_DIR / "f90templates" / filename).read_text()
                self.assertIn("ewt_fac(j) = 1d0", source)
                self.assertIn(
                    "ewt_fac(j) = abs(pdj(idx_surface_mask) * n(j) / dn_surface)",
                    source,
                )
                self.assertIn(
                    "ewt_scale = max(1d0, ewt_fac(i))", source)
                self.assertIn(
                    "abs(ycur(i))/ewt_scale", source)
                self.assertNotIn(
                    "ewt_fac(:) = 1.0d0/ewt_fac(:)", source)


if __name__ == "__main__":
    unittest.main()

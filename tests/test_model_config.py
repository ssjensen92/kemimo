import sys
import tempfile
import unittest
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src_py"
sys.path.insert(0, str(SRC_DIR))

from model_config import load_model_config


VALID_CONFIG = """\
&kemimo
  datadir = './data/'
  nlayers = 2
  layer_thickness = 4d0
  ignore_missing_h = 1
  ignore_missing_ebind = 1
  all_species_to_dust = 0
  limit_2body = 1
  respect_gasphase_limits = 0
  multiprocessing = 0
  do_swap = 0
  h2_spin = 1
  encounter_desorption = 0
  reaction_diffusion_competition = 0
  thin_ice_approximation = 1
/
"""


class ModelConfigTests(unittest.TestCase):
    def write_config(self, content):
        directory = tempfile.TemporaryDirectory()
        path = Path(directory.name) / "config.nml"
        path.write_text(content)
        self.addCleanup(directory.cleanup)
        return path

    def test_loads_numeric_booleans_and_fortran_float(self):
        config = load_model_config(self.write_config(VALID_CONFIG))
        self.assertEqual(config.datadir, "./data/")
        self.assertEqual(config.nlayers, 2)
        self.assertEqual(config.layer_thickness, 4.0)
        self.assertTrue(config.h2_spin)
        self.assertFalse(config.reaction_diffusion_competition)
        self.assertTrue(config.thin_ice_approximation)

    def test_thin_ice_approximation_defaults_are_enabled(self):
        paths = list((ROOT_DIR / "models").glob("**/config.nml"))
        for path in paths:
            with self.subTest(path=path):
                self.assertTrue(load_model_config(path).thin_ice_approximation)

    def test_three_phase_thin_ice_templates_are_selectable(self):
        database_source = (SRC_DIR / "database.py").read_text()
        templates = (
            "kemimo_ode_include_H2_thin_ice_fast.f90",
            "kemimo_ode_include_H2_nospin_thin_ice_fast.f90",
        )
        self.assertIn("if self.thinIceApproximation", database_source)
        for filename in templates:
            with self.subTest(filename=filename):
                self.assertIn(filename, database_source)

    def test_rejects_missing_variable(self):
        content = VALID_CONFIG.replace("  do_swap = 0\n", "")
        with self.assertRaisesRegex(ValueError, "do_swap"):
            load_model_config(self.write_config(content))

    def test_rejects_invalid_numeric_boolean(self):
        content = VALID_CONFIG.replace(
            "reaction_diffusion_competition = 0",
            "reaction_diffusion_competition = 2",
        )
        with self.assertRaisesRegex(ValueError, "must be 1 or 0"):
            load_model_config(self.write_config(content))

    def test_all_model_namelists_are_valid(self):
        paths = list((ROOT_DIR / "models").glob("**/config.nml"))
        self.assertGreaterEqual(len(paths), 2)
        for path in paths:
            with self.subTest(path=path):
                load_model_config(path)

    def test_all_fortran_templates_have_a_selection_path(self):
        database_source = (SRC_DIR / "database.py").read_text()
        templates = (ROOT_DIR / "f90templates").glob("*.f90")

        for template in templates:
            with self.subTest(template=template.name):
                self.assertIn(template.name, database_source)

    def test_fixed_h2_templates_and_option_are_removed(self):
        fixed_templates = (
            "kemimo_fixed_H2.f90",
            "kemimo_ode_fixed_H2.f90",
            "kemimo_ode_twophase_fixed_H2.f90",
            "kemimo_twophase_fixed_H2.f90",
        )
        for filename in fixed_templates:
            self.assertFalse((ROOT_DIR / "f90templates" / filename).exists())

        database_source = (SRC_DIR / "database.py").read_text()
        self.assertNotIn("self.include_H2", database_source)
        self.assertNotIn("include_H2,", database_source)
        self.assertNotIn("fixed_H2", database_source)


if __name__ == "__main__":
    unittest.main()

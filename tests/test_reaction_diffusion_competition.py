import inspect
import math
import re
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch


try:
    import numpy  # noqa: F401
except ModuleNotFoundError:
    sys.modules["numpy"] = types.SimpleNamespace(
        nan=math.nan,
        isnan=math.isnan,
    )

try:
    import scipy  # noqa: F401
except ModuleNotFoundError:
    scipy = types.ModuleType("scipy")
    scipy_interpolate = types.ModuleType("scipy.interpolate")
    scipy_interpolate.interp1d = None
    scipy.interpolate = scipy_interpolate
    sys.modules["scipy"] = scipy
    sys.modules["scipy.interpolate"] = scipy_interpolate

SRC_DIR = Path(__file__).resolve().parents[1] / "src_py"
sys.path.insert(0, str(SRC_DIR))

import database as database_module
from database import database
from reaction import reaction


class SurfaceSpecies:
    def __init__(self, name, mass, eice, ediff, enthalpy):
        self.name = name + "_surface"
        self.namebase = name
        self.dictname = name
        self.mass = mass
        self.Eice = eice
        self.Ediff = ediff
        self.dH = enthalpy
        self.isGas = False
        self.exploded = [name]
        self.layer = 1
        self.natoms = 1
        self.dof = 3


def evaluate_f90_expression(expression, **variables):
    expression = expression.replace("&\n", "")
    expression = re.sub(r"(?<=\d)d(?=[+-]?\d)", "e", expression)
    return eval(
        expression,
        {"__builtins__": {}, "exp": math.exp, "max": max},
        variables,
    )


class ReactionDiffusionCompetitionTests(unittest.TestCase):
    def setUp(self):
        self.reactant_i = SurfaceSpecies("A", 1.0, 400.0, 200.0, 0.0)
        self.reactant_j = SurfaceSpecies("B", 4.0, 800.0, 400.0, 0.0)
        self.product = SurfaceSpecies("C", 5.0, 1200.0, 600.0, 1.0)
        self.activation_energy = 2500.0
        barriers = [{
            "reactants": ["A", "B"],
            "products": ["C"],
            "Ea": self.activation_energy,
        }]
        self.reaction = reaction(
            [self.reactant_i, self.reactant_j],
            [self.product],
            barriers,
            barrierWidths=[],
            Bratios=[],
        )

    def expected_rates(self, temperature):
        kb = 1.38064852e-16
        hbar = 1.0545718e-27
        mp = 1.6726219e-24
        site_density = 1.5e15
        barrier_width = 1e-8
        reactants = [self.reactant_i, self.reactant_j]
        attempt_rates = [
            math.sqrt(
                2.0 * site_density * species.Eice * kb
                / math.pi**2 / (species.mass * mp)
            )
            for species in reactants
        ]
        reduced_mass = (
            self.reactant_i.mass * self.reactant_j.mass
            / (self.reactant_i.mass + self.reactant_j.mass)
        )
        tunnel_probability = math.exp(
            -2.0 * barrier_width / hbar
            * math.sqrt(
                2.0 * reduced_mass * mp * kb * self.activation_energy
            )
        )
        barrier_probability = max(
            tunnel_probability,
            math.exp(-self.activation_energy / temperature),
        )
        diffusion_rate = sum(
            attempt_rate * math.exp(-species.Ediff / temperature)
            for attempt_rate, species in zip(attempt_rates, reactants)
        )
        reaction_rate = max(attempt_rates) * barrier_probability
        return (
            reaction_rate / (reaction_rate + diffusion_rate) * diffusion_rate,
            barrier_probability * diffusion_rate,
        )

    def test_generated_rate_matches_reduced_competition_formula(self):
        self.reaction.rate(reactionDiffusionCompetition=True)
        for temperature in (20.0, 1000.0):
            with self.subTest(temperature=temperature):
                actual = evaluate_f90_expression(
                    self.reaction.krateF90,
                    invTd=1.0 / temperature,
                    indns=1.0,
                )
                expected, _ = self.expected_rates(temperature)
                self.assertAlmostEqual(actual / expected, 1.0, places=7)

    def test_competition_can_be_disabled(self):
        temperature = 20.0
        self.reaction.rate(reactionDiffusionCompetition=False)
        actual = evaluate_f90_expression(
            self.reaction.krateF90,
            invTd=1.0 / temperature,
            indns=1.0,
        )
        _, expected = self.expected_rates(temperature)
        self.assertAlmostEqual(actual / expected, 1.0, places=7)

    def test_database_forwards_configured_value(self):
        self.assertEqual(
            list(inspect.signature(database.__init__).parameters), ["self"])
        calls = []
        fake_reaction = types.SimpleNamespace(
            rate=lambda **kwargs: calls.append(kwargs))
        model = database.__new__(database)
        model.reactions = [fake_reaction]
        model.reactionDiffusionCompetition = False
        model.computeRates()
        self.assertEqual(
            calls, [{"reactionDiffusionCompetition": False}])


class SpinTemplateSelectionTests(unittest.TestCase):
    def selected_templates(self, h2_spin):
        model = database.__new__(database)
        model.nlayers = 2
        model.H2spin = h2_spin
        copies = []

        def record_copy(source, destination):
            copies.append((source, destination))
            if len(copies) == 4:
                raise RuntimeError("template selection captured")

        with patch.object(
                database_module, "copyfile", side_effect=record_copy):
            with self.assertRaisesRegex(
                    RuntimeError, "template selection captured"):
                model.preproc()
        return copies

    def test_spin_and_nonspin_templates_both_include_dynamic_h2(self):
        spin_templates = self.selected_templates(h2_spin=True)
        nonspin_templates = self.selected_templates(h2_spin=False)
        self.assertEqual(
            spin_templates[0][0],
            "./f90templates/kemimo_ode_include_H2.f90",
        )
        self.assertEqual(
            nonspin_templates[0][0],
            "./f90templates/kemimo_ode_include_H2_nospin.f90",
        )
        self.assertEqual(
            spin_templates[1][0],
            "./f90templates/kemimo_include_H2.f90",
        )
        self.assertEqual(nonspin_templates[1], spin_templates[1])
        self.assertEqual(
            spin_templates[3][0],
            "./f90templates/kemimo_rates.f90",
        )
        self.assertEqual(
            nonspin_templates[3][0],
            "./f90templates/kemimo_rates_nospin.f90",
        )

if __name__ == "__main__":
    unittest.main()

"""Software checks using sourced clay coordinates, NOT equilibrium validation."""
from argparse import Namespace
import json
from pathlib import Path
import tempfile
import unittest

import ase.io
import numpy as np
import openmm as mm
from openmm import unit

import test36_aa as aa


class ExplicitClayChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="test36-aa-check-")
        cls.input = Path(cls.temporary.name) / "input"
        aa.build(Namespace(source=aa.SOURCE, output=cls.input, nx=4, ny=4,
                           layers=2, waters_per_cell=12, spacing_nm=1.55, seed=3601))
        cls.meta = json.loads((cls.input / "provenance.json").read_text())
        cls.positions = np.load(cls.input / "initial_positions_nm.npy")

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_composition_geometry_and_mapping(self):
        self.assertEqual(self.meta["atoms"], 2444)
        self.assertAlmostEqual(self.meta["net_charge_e"], 0., places=9)
        self.assertAlmostEqual(self.meta["clay_charge_e"], -16., places=9)
        self.assertTrue(self.meta["mg_coordination_verified"])
        numbers = [a["atomic_number"] for a in self.meta["atom_table"]]
        self.assertEqual(numbers.count(14), 256)
        self.assertEqual(numbers.count(13), 112)
        self.assertEqual(numbers.count(12), 16)
        self.assertEqual(numbers.count(11), 8)
        self.assertEqual(numbers.count(20), 4)
        water = self.positions[1292:].reshape(-1, 3, 3)
        np.testing.assert_allclose(np.linalg.norm(water[:, 1:]-water[:, :1], axis=-1), .1, atol=1e-14)
        np.testing.assert_allclose(np.linalg.norm(water[:, 1]-water[:, 2], axis=-1), .1633, atol=1e-14)
        mapping = json.loads((self.input / "mapping.json").read_text())
        self.assertFalse(mapping["source_is_equilibrium"])
        self.assertEqual([s["indices"] for s in mapping["sites"]], [[i] for i in range(1292)])
        self.assertFalse({s["name"] for s in mapping["species"]} & {"OW", "HW"})
        self.assertEqual(len(mapping["species"]), 10)
        self.assertTrue((self.input / "ClayCode.LICENSE.txt").is_file())

    def test_imported_forcefield_and_exclusions(self):
        top, system = aa.make_system(self.input, 1.)
        aa.validate_system(self.input, system, self.meta)
        nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
        self.assertFalse(nb.getUseDispersionCorrection())
        self.assertTrue(nb.getExceptionsUsePeriodicBoundaryConditions())
        params = {r[0]: r for r in aa.sections(self.input / "system.top")["atomtypes"]}
        for i, atom in enumerate(self.meta["atom_table"]):
            q, sigma, epsilon = nb.getParticleParameters(i)
            expected = params[atom["type"]]
            self.assertAlmostEqual(sigma.value_in_unit(unit.nanometer), float(expected[5]))
            self.assertAlmostEqual(epsilon.value_in_unit(unit.kilojoule_per_mole), float(expected[6]))
        bonds = next(f for f in system.getForces() if isinstance(f, mm.HarmonicBondForce))
        self.assertTrue(bonds.usesPeriodicBoundaryConditions())
        self.assertEqual(bonds.getNumBonds(), 128)
        for i in range(bonds.getNumBonds()):
            a, b, d, k = bonds.getBondParameters(i)
            self.assertAlmostEqual(d.value_in_unit(unit.nanometer), .1)
            self.assertAlmostEqual(k.value_in_unit(unit.kilojoule_per_mole/unit.nanometer**2), 463532.8)
        self.assertEqual({a.element.symbol for a in top.topology.atoms()}, {"Si", "Al", "Mg", "O", "H", "Na", "Ca"})

    def test_periodic_single_point_and_short_md(self):
        _, system = aa.make_system(self.input, 1.)
        integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, .0005*unit.picosecond)
        context = mm.Context(system, integrator, mm.Platform.getPlatformByName("CPU"), {"Threads": "1"})
        context.setPositions(self.positions*unit.nanometer)
        initial = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        # Translation by an entire periodic cell must not change the energy.
        translated = self.positions + np.array(self.meta["box_nm"])
        context.setPositions(translated*unit.nanometer)
        moved = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        self.assertAlmostEqual(initial, moved, delta=.1)
        context.setVelocitiesToTemperature(50*unit.kelvin, 33)
        integrator.step(4)
        state = context.getState(getPositions=True, getEnergy=True)
        self.assertTrue(np.isfinite(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).all())
        self.assertTrue(np.isfinite(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)))
        del context, integrator

    def test_refuse_existing_output(self):
        with self.assertRaises(ValueError):
            aa.fresh(self.input)


if __name__ == "__main__":
    unittest.main()

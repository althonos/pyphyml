
import unittest

import Bio.AlignIO

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files

from pyphyml import Alignment, TreeBuilder, Tree, ModelPrototype



class TestTreeBuilder(unittest.TestCase):

    TEST_SEED = 0

    @classmethod
    def setUpClass(cls):
        with resource_files(__package__).joinpath("data", "nucleic").open() as f:
            msa = Bio.AlignIO.read(f, "phylip")
        cls.nucleic = Alignment(names=[seq.id for seq in msa], sequences=msa.alignment)

    def test_nucleic_default(self):
        builder = TreeBuilder(seed=self.TEST_SEED)
        result = builder.build(self.nucleic)

        with resource_files(__package__).joinpath("data", "nucleic.default.tree").open() as f:
            expected = f.read().strip()

        self.assertAlmostEqual(result.log_likelihood, -5593.07669, places=4)
        self.assertMultiLineEqual(result.tree.dumps().strip(), expected)

    def test_nucleic_alpha(self):
        model = ModelPrototype.from_name("HKY85", alpha=0.1)
        builder = TreeBuilder(model, seed=self.TEST_SEED)
        result = builder.build(self.nucleic)

        with resource_files(__package__).joinpath("data", "nucleic.alpha0_1.tree").open() as f:
            expected = f.read().strip()

        self.assertAlmostEqual(result.log_likelihood, -5479.85392, places=4)
        # self.assertMultiLineEqual(result.tree.dumps().strip(), expected)  # FIXME

    def test_nucleic_model_F84(self):
        model = ModelPrototype.from_name("F84")
        builder = TreeBuilder(model, seed=self.TEST_SEED)
        result = builder.build(self.nucleic)

        with resource_files(__package__).joinpath("data", "nucleic.modelF84.tree").open() as f:
            expected = f.read().strip()

        self.assertAlmostEqual(result.log_likelihood, -5593.07669, places=4)
        self.assertMultiLineEqual(result.tree.dumps().strip(), expected)  # FIXME



import unittest

import Bio.AlignIO
import Bio.Phylo.NewickIO

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files

from pyphyml import Alignment, TreeBuilder, Tree, ModelPrototype



class TestTreeBuilder(unittest.TestCase):

    TEST_SEED = 42

    @classmethod
    def setUpClass(cls):
        with resource_files(__package__).joinpath("data", "nucleic").open() as f:
            msa = Bio.AlignIO.read(f, "phylip")
        cls.nucleic = Alignment(names=[seq.id for seq in msa], sequences=msa.alignment)

    def assertTreeAlmostEqual(self, t1, t2):
        clades1 = sorted(t1.find_elements(terminal=True), key=lambda clade: clade.name)
        clades2 = sorted(t2.find_elements(terminal=True), key=lambda clade: clade.name)
        self.assertEqual(len(clades1), len(clades2))
        for c1, c2 in zip(clades2, clades2):
            self.assertAlmostEqual(c1.name, c2.name)
            self.assertAlmostEqual(c1.branch_length, c2.branch_length, places=4)

    def test_nucleic_default(self):
        builder = TreeBuilder(seed=self.TEST_SEED)
        result = builder.build(self.nucleic)

        with resource_files(__package__).joinpath("data", "nucleic.default.tree").open() as f:
            expected = f.read().strip()

        tree_exp = next(Bio.Phylo.NewickIO.Parser.from_string(str(expected)).parse())
        tree_actual = next(Bio.Phylo.NewickIO.Parser.from_string(result.tree.dumps()).parse())

        self.assertAlmostEqual(result.log_likelihood, -5593.07667604, places=4)
        self.assertTreeAlmostEqual(tree_actual, tree_exp)

    def test_nucleic_alpha(self):
        model = ModelPrototype.from_name("HKY85", alpha=0.1)
        builder = TreeBuilder(model, seed=self.TEST_SEED)
        result = builder.build(self.nucleic)

        with resource_files(__package__).joinpath("data", "nucleic.alpha0_1.tree").open() as f:
            expected = f.read().strip()

        tree_exp = next(Bio.Phylo.NewickIO.Parser.from_string(str(expected)).parse())
        tree_actual = next(Bio.Phylo.NewickIO.Parser.from_string(result.tree.dumps()).parse())

        self.assertAlmostEqual(result.log_likelihood, -5488.198436, places=4)
        self.assertTreeAlmostEqual(tree_actual, tree_exp)

    def test_nucleic_model_F84(self):
        model = ModelPrototype.from_name("F84")
        builder = TreeBuilder(model, seed=self.TEST_SEED)
        result = builder.build(self.nucleic)

        with resource_files(__package__).joinpath("data", "nucleic.modelF84.tree").open() as f:
            expected = f.read().strip()
        
        tree_exp = next(Bio.Phylo.NewickIO.Parser.from_string(str(expected)).parse())
        tree_actual = next(Bio.Phylo.NewickIO.Parser.from_string(result.tree.dumps()).parse())

        self.assertAlmostEqual(result.log_likelihood, -5593.076676046, places=4)
        self.assertTreeAlmostEqual(tree_actual, tree_exp)


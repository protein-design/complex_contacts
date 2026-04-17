
from .helper import *
from src.contacts import CalcSasa


@ddt
class TestCalcSasa(TestCase):

    
    @data(
        ['7S8L', None, None, (1096, 12)],
        ['7S8L', 'S', None, (1096, 12)],
        ['7S8L', None, True, (1096, 12)],
    )
    @unpack
    def test_freesasa(self, pdb_id, method, hetatm, expect):
        pdb_file = os.path.join(DATA_DIR, f'{pdb_id}.pdb')
        c = CalcSasa(pdb_file)
        result = c.freesasa(method, hetatm)
        df = c.collect_freesasa_result(result)
        assert df.shape == expect



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

    @data(
        [
            '1I8I_0_A-B', '1I8I_0_A', ['chain_id', 'res_no', 'res', 'aa', 'value'],
            [{'seq': 'DIELTQSPA', 'sig_res': 2, 'start': 0, 'end': 8, 'cutoff': 1, 'max_span': 10}, 
                {'seq': 'NWYQQKPGEPPKFLISEGNTLR', 'sig_res': 10, 'start': 33, 'end': 54, 'cutoff': 1, 'max_span': 10}, 
                {'seq': 'DYYCLQSFNVPLTFGCGTK', 'sig_res': 13, 'start': 84, 'end': 102, 'cutoff': 1, 'max_span': 10}]
        ],
    )
    @unpack
    def test_delta_sasa(self, pdb_id, unbounded_id, expect_delta, expect_motif):
        # bounded
        pdb_file = os.path.join(DATA_DIR, f'{pdb_id}.pdb')
        c = CalcSasa(pdb_file)
        result = c.freesasa()
        bounded_df = c.collect_freesasa_result(result)
        # unbounded
        unbounded_file = os.path.join(DATA_DIR, f'{unbounded_id}.pdb')
        c = CalcSasa(unbounded_file)
        result = c.freesasa()
        unbounded_df = c.collect_freesasa_result(result)
        
        # delta SASA
        delta = CalcSasa.delta_sasa(bounded_df, unbounded_df)
        assert list(delta) == expect_delta

        #parse seq
        motif = CalcSasa.parse_seq(delta, cutoff=1, max_span=10)
        assert motif == expect_motif
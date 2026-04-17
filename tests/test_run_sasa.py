
from .helper import *
from src.contacts import RunSasa


@ddt
class TestRunSasa(TestCase):

    
    @data(
        ['1I8I_0_A-B', ['1I8I_0_A', '1I8I_0_B'], 1706.71],
        ['wrong', ['1I8I_0_A', '1I8I_0_B'], None],
        ['1I8I_0_A-B', ['1I8I_0_A', 'wrong'], None],
    )
    @unpack
    def test_freesasa(self, bounded_pdb_id, unbounded_pdb_ids, expect_delta):
        pdb_file = os.path.join(DATA_DIR, f'{bounded_pdb_id}.pdb')
        unbounded_pdb_files = [os.path.join(DATA_DIR, f'{pdb_id}.pdb') for pdb_id in unbounded_pdb_ids]
        c = RunSasa()
        res = c.pair_freesasa(pdb_file, unbounded_pdb_files)
        assert res.get('delta_total_sasa') == expect_delta

    @data(
        [
            '1I8I_0_A-B', ['1I8I_0_A', '1I8I_0_B'], 1706.71,
            [{'seq': 'DIELTQSPASLSVATGEKVTIRCMTSTDIDDDMNWYQQKPGEPPKFLISEGNTLRPGVPSRFSSSGTGTDFVFTIENTLSEDVGDYYCLQSFNVPLTFGCGTK', 
                'sig_res': 25, 'start': 0, 'end': 102, 'cutoff': 1, 'max_span': 30}]
        ],
    )
    @unpack
    def test_run(self, bounded_pdb_id, unbounded_pdb_ids, expect_delta, expect_motif):
        pdb_file = os.path.join(DATA_DIR, f'{bounded_pdb_id}.pdb')
        unbounded_pdb_files = [os.path.join(DATA_DIR, f'{pdb_id}.pdb') for pdb_id in unbounded_pdb_ids]
        pool = [(pdb_file, unbounded_pdb_files),]
        pfile = os.path.join(RESULTS_DIR, f'{bounded_pdb_id}.pkl')
        res = RunSasa().run(pool, pfile)
        assert len(res) == 1
        assert res[0].get('delta_total_sasa') == expect_delta
        assert res[0].get('unbounded')[0].get('motif') == expect_motif
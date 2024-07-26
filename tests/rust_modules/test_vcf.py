import pytest
from pathlib import Path
from viola import rs
from viola.rust.vcf import Vcf
import threading

class TestRustVcf:

    @pytest.fixture
    def vcf(self):
        path = Path(__file__).parent / "data/test.manta.vcf"
        rs_vcf = rs.read_vcf(path, "manta")
        return Vcf(rs_vcf)
    
    def test_sv_count(self, vcf):
        assert vcf.sv_count == 6
    
    def test_contigs(self, vcf):
        assert vcf.contigs == ["chr1", "chr2", "chr11"]
    
    def test_ids(self, vcf):
        assert vcf.ids == ["test1", "test2", "test4_1", "test3", "test4_2", "test5"]
    
    def test_get_positions(self, vcf):
        positions = vcf.get_positions_table()
    
    def test_get_str_info(self, vcf):
        str_info = vcf.get_str_info_table()
    
    def test_get_str_format(self, vcf):
        str_format = vcf.get_str_format_table()
        


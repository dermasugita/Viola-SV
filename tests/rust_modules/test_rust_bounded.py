from viola import rs
from pathlib import Path



def test_rust_bounded() -> None:
    assert rs.rust_bounded("1") == "1"

def test_read_vcf() -> None:
    res = rs.read_vcf(Path(__file__).parent / "../../examples/demo_merge/test.merge.manta.vcf")
    print(res.get_by_id("MD1"))
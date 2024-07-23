from viola import rs


def test_rust_bounded() -> None:
    assert rs.rust_bounded("1") == "1"
import viola

def test_get_inslen_and_insseq_from_alt():
    result1 = viola.get_inslen_and_insseq_from_alt(']chr8:53306291]TTTTTTTTTGTTGT')
    assert result1 == (13, 'TTTTTTTTTGTTG')
    result2 = viola.get_inslen_and_insseq_from_alt('TCTCCCTCCCTCC[chr11:20592279[')
    assert result2 == (12, 'CTCCCTCCCTCC')
    result3 = viola.get_inslen_and_insseq_from_alt(']chr8:53306291]T')
    assert result3 == (0, '')
    result4 = viola.get_inslen_and_insseq_from_alt('T[chr11:20592279[')
    assert result4 == (0, '')
    result5 = viola.get_inslen_and_insseq_from_alt('<DEL>')
    assert result5 == (0, '')

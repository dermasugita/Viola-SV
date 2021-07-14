import viola 

# non-overlap examples
test_param1 = {
    'chr1h': 'chr1',
    'chr2h': 'chr1',
    'pos1h': 10,
    'pos2h': 20,
    'chr1w': 'chr1',
    'chr2w': 'chr1',
    'pos1w': 21,
    'pos2w': 31
}
test_param2 = {
    'chr1h': 'chr1',
    'chr2h': 'chr1',
    'pos1h': 20,
    'pos2h': 10,
    'chr1w': 'chr1',
    'chr2w': 'chr1',
    'pos1w': 21,
    'pos2w': 31
}
test_param3 = {
    'chr1h': 'chr1',
    'chr2h': 'chr1',
    'pos1h': 20,
    'pos2h': 10,
    'chr1w': 'chr1',
    'chr2w': 'chr1',
    'pos1w': 31,
    'pos2w': 21
}
test_param4 = {
    'chr1h': 'chr1',
    'chr2h': 'chr2',
    'pos1h': 10,
    'pos2h': 20,
    'chr1w': 'chr1',
    'chr2w': 'chr2',
    'pos1w': 11,
    'pos2w': 21
}

# overlap examples
test_param5 = {
    'chr1h': 'chr1',
    'chr2h': 'chr1',
    'pos1h': 10,
    'pos2h': 20,
    'chr1w': 'chr1',
    'chr2w': 'chr1',
    'pos1w': 11,
    'pos2w': 21
}
test_param6 = {
    'chr1h': 'chr1',
    'chr2h': 'chr1',
    'pos1h': 20,
    'pos2h': 10,
    'chr1w': 'chr1',
    'chr2w': 'chr1',
    'pos1w': 11,
    'pos2w': 21
}
test_param7 = {
    'chr1h': 'chr1',
    'chr2h': 'chr1',
    'pos1h': 20,
    'pos2h': 10,
    'chr1w': 'chr1',
    'chr2w': 'chr1',
    'pos1w': 21,
    'pos2w': 11
}
test_param8 = {
    'chr1h': 'chr1',
    'chr2h': 'chr1',
    'pos1h': 10,
    'pos2h': 20,
    'chr1w': 'chr1',
    'chr2w': 'chr1',
    'pos1w': 10,
    'pos2w': 20
}
test_param9 = {
    'chr1h': 'chr1',
    'chr2h': 'chr1',
    'pos1h': 10,
    'pos2h': 20,
    'chr1w': 'chr1',
    'chr2w': 'chr1',
    'pos1w': 9,
    'pos2w': 21
}
test_param10 = {
    'chr1h': 'chr1',
    'chr2h': 'chr1',
    'pos1h': 10,
    'pos2h': 20,
    'chr1w': 'chr1',
    'chr2w': 'chr1',
    'pos1w': 11,
    'pos2w': 19 
}

def test_nonoverlap():
    assert viola.Bedpe._nonoverlap('self', test_param1) == True
    assert viola.Bedpe._nonoverlap('self', test_param2) == True
    assert viola.Bedpe._nonoverlap('self', test_param3) == True
    assert viola.Bedpe._nonoverlap('self', test_param4) == False
    assert viola.Bedpe._nonoverlap('self', test_param5) == False
    assert viola.Bedpe._nonoverlap('self', test_param6) == False
    assert viola.Bedpe._nonoverlap('self', test_param7) == False
    assert viola.Bedpe._nonoverlap('self', test_param8) == False
    assert viola.Bedpe._nonoverlap('self', test_param9) == False
    assert viola.Bedpe._nonoverlap('self', test_param10) == False
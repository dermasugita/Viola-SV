import viola
import pandas as pd
from pandas.testing import assert_frame_equal
from io import StringIO
data1	=	"""chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	info_0	info_1
chr1	10	11	chr1	20	21	test1	60	+	-	A	1
chr1	150	151	chr1	300	301	test2	60	+	-	B	1
chr3	10	11	chr4	20	21	test3	60	-	-	F	4
"""

data2	=	"""chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	info_0	info_1
chr1	100	101	chr1	250	251	test1	60	+	-	A	2
chr2	180	181	chr2	3000	3001	test2	60	-	-	S	3
chr2	10	11	chr5	20	21	test3	60	+	-	L	6
"""
data3	=	"""chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	info_1	info_2
chr2	10	11	chr5	20	21	test1	60	+	-	7	T
chr3	10	11	chr4	20	21	test2	60	-	-	1	F
"""
data4	=	"""chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	info_1	info_2
chr2	15	16	chr2	20	21	test1	60	-	-	10	T
chr3	10	11	chr4	20	21	test2	60	-	-	1	F
"""
data5	=	"""chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	info_1	info_3
chr1	105	106	chr1	290	291	test1	60	+	-	2	Y
chr2	100	101	chr2	280	281	test2	60	-	-	2	Y
"""

data6	=	"""chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	info_1	info_3
"""


def test_concat1():
    bedpe1 = viola.read_bedpe(StringIO(data1), patient_name='patient1')
    bedpe2 = viola.read_bedpe(StringIO(data2), patient_name='patient2')
    multi_bedpe1 = viola.MultiBedpe([bedpe1, bedpe2], ['patient1', 'patient2'])
    bedpe3 = viola.read_bedpe(StringIO(data3), patient_name='patient3')
    concat_bedpe = viola.concat([multi_bedpe1, bedpe3])
    df_id_expected = pd.read_csv(StringIO("""global_id,patient_id,id
patient1_test1,0,test1
patient1_test2,0,test2
patient1_test3,0,test3
patient2_test1,1,test1
patient2_test2,1,test2
patient2_test3,1,test3
patient3_test1,2,test1
patient3_test2,2,test2
"""))
    df_patients_expected = pd.read_csv(StringIO("""id,patients
0,patient1
1,patient2
2,patient3
"""))
    df_positions_expected = pd.read_csv(StringIO("""id,chrom1,pos1,chrom2,pos2,strand1,strand2,ref,alt,qual,svtype
patient1_test1,chr1,11,chr1,21,+,-,N,<DEL>,60,DEL
patient1_test2,chr1,151,chr1,301,+,-,N,<DEL>,60,DEL
patient1_test3,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient2_test1,chr1,101,chr1,251,+,-,N,<DEL>,60,DEL
patient2_test2,chr2,181,chr2,3001,-,-,N,<INV>,60,INV
patient2_test3,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test1,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test2,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
"""))
    df_id = concat_bedpe.get_table('global_id')
    df_patients = concat_bedpe.get_table('patients')
    df_positions = concat_bedpe.get_table('positions')
    assert_frame_equal(df_id, df_id_expected)
    assert_frame_equal(df_patients, df_patients_expected)
    assert_frame_equal(df_positions, df_positions_expected)
    print(concat_bedpe)
    
def test_concat2():
    bedpe1 = viola.read_bedpe(StringIO(data1), patient_name='patient1')
    bedpe2 = viola.read_bedpe(StringIO(data2), patient_name='patient2')
    multi_bedpe1 = viola.MultiBedpe([bedpe1, bedpe2], ['patient1', 'patient2'])
    bedpe3 = viola.read_bedpe(StringIO(data3), patient_name='patient3')
    bedpe4 = viola.read_bedpe(StringIO(data4), patient_name='patient4')
    multi_bedpe2 = viola.MultiBedpe([bedpe3, bedpe4], ['patient3', 'patient4'])
    bedpe5 = viola.read_bedpe(StringIO(data5), patient_name='patient5')
    concat_bedpe = viola.concat([multi_bedpe1, multi_bedpe2, bedpe5])
    df_id_expected = pd.read_csv(StringIO("""global_id,patient_id,id
patient1_test1,0,test1
patient1_test2,0,test2
patient1_test3,0,test3
patient2_test1,1,test1
patient2_test2,1,test2
patient2_test3,1,test3
patient3_test1,2,test1
patient3_test2,2,test2
patient4_test1,3,test1
patient4_test2,3,test2
patient5_test1,4,test1
patient5_test2,4,test2
"""))
    df_patients_expected = pd.read_csv(StringIO("""id,patients
0,patient1
1,patient2
2,patient3
3,patient4
4,patient5
"""))
    df_positions_expected = pd.read_csv(StringIO("""id,chrom1,pos1,chrom2,pos2,strand1,strand2,ref,alt,qual,svtype
patient1_test1,chr1,11,chr1,21,+,-,N,<DEL>,60,DEL
patient1_test2,chr1,151,chr1,301,+,-,N,<DEL>,60,DEL
patient1_test3,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient2_test1,chr1,101,chr1,251,+,-,N,<DEL>,60,DEL
patient2_test2,chr2,181,chr2,3001,-,-,N,<INV>,60,INV
patient2_test3,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test1,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test2,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient4_test1,chr2,16,chr2,21,-,-,N,<INV>,60,INV
patient4_test2,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient5_test1,chr1,106,chr1,291,+,-,N,<DEL>,60,DEL
patient5_test2,chr2,101,chr2,281,-,-,N,<INV>,60,INV
"""))
    df_id =  concat_bedpe.get_table('global_id')
    df_patients = concat_bedpe.get_table('patients')
    df_positions = concat_bedpe.get_table('positions')
    assert_frame_equal(df_id, df_id_expected)
    assert_frame_equal(df_patients, df_patients_expected)
    assert_frame_equal(df_positions, df_positions_expected)

def test_concat3():
    bedpe1 = viola.read_bedpe(StringIO(data1), patient_name='patient1')
    bedpe2 = viola.read_bedpe(StringIO(data2), patient_name='patient2')
    multi_bedpe1 = viola.MultiBedpe([bedpe1, bedpe2], ['patient1', 'patient2'])
    bedpe3 = viola.read_bedpe(StringIO(data3), patient_name='patient3')
    bedpe4 = viola.read_bedpe(StringIO(data4), patient_name='patient4')
    multi_bedpe2 = viola.MultiBedpe([bedpe3, bedpe4], ['patient3', 'patient4'])
    bedpe6 = viola.read_bedpe(StringIO(data6), patient_name='patient6')
    concat_bedpe = viola.concat([multi_bedpe1, multi_bedpe2, bedpe6])
    df_id_expected = pd.read_csv(StringIO("""global_id,patient_id,id
patient1_test1,0,test1
patient1_test2,0,test2
patient1_test3,0,test3
patient2_test1,1,test1
patient2_test2,1,test2
patient2_test3,1,test3
patient3_test1,2,test1
patient3_test2,2,test2
patient4_test1,3,test1
patient4_test2,3,test2
"""))
    df_patients_expected = pd.read_csv(StringIO("""id,patients
0,patient1
1,patient2
2,patient3
3,patient4
4,patient6
"""))
    df_positions_expected = pd.read_csv(StringIO("""id,chrom1,pos1,chrom2,pos2,strand1,strand2,ref,alt,qual,svtype
patient1_test1,chr1,11,chr1,21,+,-,N,<DEL>,60,DEL
patient1_test2,chr1,151,chr1,301,+,-,N,<DEL>,60,DEL
patient1_test3,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient2_test1,chr1,101,chr1,251,+,-,N,<DEL>,60,DEL
patient2_test2,chr2,181,chr2,3001,-,-,N,<INV>,60,INV
patient2_test3,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test1,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test2,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient4_test1,chr2,16,chr2,21,-,-,N,<INV>,60,INV
patient4_test2,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
"""))
    df_id = concat_bedpe.get_table('global_id')
    df_patients = concat_bedpe.get_table('patients')
    df_positions = concat_bedpe.get_table('positions')
    assert_frame_equal(df_id, df_id_expected)
    assert_frame_equal(df_patients, df_patients_expected)
    assert_frame_equal(df_positions, df_positions_expected)

def test_concat4():
    bedpe1 = viola.read_bedpe(StringIO(data1), patient_name='patient1')
    bedpe2 = viola.read_bedpe(StringIO(data2), patient_name='patient2')
    multi_bedpe1 = viola.MultiBedpe([bedpe1, bedpe2], ['patient1', 'patient2'])
    bedpe3 = viola.read_bedpe(StringIO(data3), patient_name='patient3')
    bedpe4 = viola.read_bedpe(StringIO(data4), patient_name='patient4')
    multi_bedpe2 = viola.MultiBedpe([bedpe3, bedpe4], ['patient3', 'patient4'])
    bedpe6 = viola.read_bedpe(StringIO(data6), patient_name='patient6')
    concat_bedpe = viola.concat([bedpe6, multi_bedpe1, multi_bedpe2])
    df_id_expected = pd.read_csv(StringIO("""global_id,patient_id,id
patient1_test1,1,test1
patient1_test2,1,test2
patient1_test3,1,test3
patient2_test1,2,test1
patient2_test2,2,test2
patient2_test3,2,test3
patient3_test1,3,test1
patient3_test2,3,test2
patient4_test1,4,test1
patient4_test2,4,test2
"""))
    df_patients_expected = pd.read_csv(StringIO("""id,patients
0,patient6
1,patient1
2,patient2
3,patient3
4,patient4
"""))
    df_positions_expected = pd.read_csv(StringIO("""id,chrom1,pos1,chrom2,pos2,strand1,strand2,ref,alt,qual,svtype
patient1_test1,chr1,11,chr1,21,+,-,N,<DEL>,60,DEL
patient1_test2,chr1,151,chr1,301,+,-,N,<DEL>,60,DEL
patient1_test3,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient2_test1,chr1,101,chr1,251,+,-,N,<DEL>,60,DEL
patient2_test2,chr2,181,chr2,3001,-,-,N,<INV>,60,INV
patient2_test3,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test1,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test2,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient4_test1,chr2,16,chr2,21,-,-,N,<INV>,60,INV
patient4_test2,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
"""))
    df_id = concat_bedpe.get_table('global_id')
    df_patients = concat_bedpe.get_table('patients')
    df_positions = concat_bedpe.get_table('positions')
    assert_frame_equal(df_id, df_id_expected)
    assert_frame_equal(df_patients, df_patients_expected)
    assert_frame_equal(df_positions, df_positions_expected)
def test_concat5():
    bedpe1 = viola.read_bedpe(StringIO(data1), patient_name='patient1')
    bedpe2 = viola.read_bedpe(StringIO(data2), patient_name='patient2')
    multi_bedpe1 = viola.MultiBedpe([bedpe1, bedpe2], ['patient1', 'patient2'])
    bedpe3 = viola.read_bedpe(StringIO(data3), patient_name='patient3')
    bedpe4 = viola.read_bedpe(StringIO(data4), patient_name='patient4')
    multi_bedpe2 = viola.MultiBedpe([bedpe3, bedpe4], ['patient3', 'patient4'])
    bedpe6 = viola.read_bedpe(StringIO(data6), patient_name='patient6')
    concat_bedpe = viola.concat([multi_bedpe1, bedpe6, multi_bedpe2])
    df_id_expected = pd.read_csv(StringIO("""global_id,patient_id,id
patient1_test1,0,test1
patient1_test2,0,test2
patient1_test3,0,test3
patient2_test1,1,test1
patient2_test2,1,test2
patient2_test3,1,test3
patient3_test1,3,test1
patient3_test2,3,test2
patient4_test1,4,test1
patient4_test2,4,test2
"""))
    df_patients_expected = pd.read_csv(StringIO("""id,patients
0,patient1
1,patient2
2,patient6
3,patient3
4,patient4
"""))
    df_positions_expected = pd.read_csv(StringIO("""id,chrom1,pos1,chrom2,pos2,strand1,strand2,ref,alt,qual,svtype
patient1_test1,chr1,11,chr1,21,+,-,N,<DEL>,60,DEL
patient1_test2,chr1,151,chr1,301,+,-,N,<DEL>,60,DEL
patient1_test3,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient2_test1,chr1,101,chr1,251,+,-,N,<DEL>,60,DEL
patient2_test2,chr2,181,chr2,3001,-,-,N,<INV>,60,INV
patient2_test3,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test1,chr2,11,chr5,21,+,-,N,N[chr5:21[,60,BND
patient3_test2,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
patient4_test1,chr2,16,chr2,21,-,-,N,<INV>,60,INV
patient4_test2,chr3,11,chr4,21,-,-,N,[chr4:21[N,60,BND
"""))
    df_id = concat_bedpe.get_table('global_id')
    df_patients = concat_bedpe.get_table('patients')
    df_positions = concat_bedpe.get_table('positions')
    assert_frame_equal(df_id, df_id_expected)
    assert_frame_equal(df_patients, df_patients_expected)
    assert_frame_equal(df_positions, df_positions_expected)
import pandas as pd
import viola
from viola.core.bedpe import Bedpe
from viola.core.vcf import Vcf
from viola.core.cohort import MultiBedpe, MultiVcf
from collections import OrderedDict
from viola._exceptions import DuplicatedPatientIDError

def _detect_patient_name_duplication(ls_left, ls_right):
    ls_left = [str(i) for i in ls_left]
    ls_right = [str(i) for i in ls_right]
    if len(set(ls_left) | set(ls_right)) != len(set(ls_left)) + len(set(ls_right)):
        ls_duplicated_patients = list(set(ls_left) & set(ls_right))
        raise DuplicatedPatientIDError(', '.join(ls_duplicated_patients))

def concat_bedpe(ls_objs):
    dict_ls_df_info = dict() 
    loop_count = 0
    ls_patients = []
    ls_df_global_id = []
    ls_df_patients = []
    ls_df_positions = []
    max_patient_id = -1
    for obj in ls_objs:
        # checking type
        if not isinstance(obj, Bedpe) | isinstance(obj, MultiBedpe):
            raise TypeError('Input values should be Bedpe or MultiBedpe class.')
        if type(obj) is Bedpe:
            obj = MultiBedpe([obj], [obj.patient_name])
        if obj.sv_count == 0:
            df_patients_current = obj.get_table('patients')
            ls_patients_current = df_patients_current.loc[:, 'patients'].to_list()
            df_patients_current.loc[:, 'id'] = df_patients_current.loc[:, 'id'] + max_patient_id + 1
            ls_df_patients.append(df_patients_current)
            max_patient_id = df_patients_current.loc[:, 'id'].max()
            _detect_patient_name_duplication(ls_patients, ls_patients_current)
            ls_patients += ls_patients_current
            continue
        df_patients_current = obj.get_table('patients')
        ls_patients_current = df_patients_current.loc[:, 'patients'].to_list()
        df_patients_current.loc[:, 'id'] = df_patients_current.loc[:, 'id'] + max_patient_id + 1
        ls_df_patients.append(df_patients_current)

        _detect_patient_name_duplication(ls_patients, ls_patients_current)

        ls_patients += ls_patients_current

        df_global_id_current = obj.get_table('global_id')
        df_global_id_current.loc[:, 'patient_id'] = df_global_id_current.loc[:, 'patient_id'] + max_patient_id + 1
        ls_df_global_id.append(df_global_id_current)
        ls_df_positions.append(obj.get_table('positions'))

        max_patient_id = df_patients_current.loc[:, 'id'].max() # update max_patient_id

        if loop_count == 0:
            for key, value in obj._odict_df_info.items():
                dict_ls_df_info[key] = [value.copy()]
            loop_count = loop_count + 1
        else:
            for key, value in obj._odict_df_info.items():
                value = value.copy()
                if dict_ls_df_info.get(key) is None:
                    dict_ls_df_info[key] = [value]
                else:
                    dict_ls_df_info[key].append(value)

    df_concat_global_id = pd.concat(ls_df_global_id, ignore_index=True)
    df_concat_patients = pd.concat(ls_df_patients, ignore_index=True)
    df_concat_positions = pd.concat(ls_df_positions, ignore_index=True)
    odict_df_info = OrderedDict()
    for key, value in dict_ls_df_info.items():
        odict_df_info[key] = pd.concat(value)
    
    return MultiBedpe(direct_tables=[df_concat_global_id, df_concat_patients, df_concat_positions, odict_df_info])
        


def concat(ls_objs):
    """
    concat(ls_objs)
    Concatenate multiple Bedpe/MultiBedpe objects into single MultiBedpe object.

    Parameters
    --------------
    ls_objs: List[Bedpe or MultiBedpe]
        A list of Bedpe/MultiBedpe objects to be concatenated.
    
    Returns
    ---------
    MultiBedpe
        The MultiBedpe object including all patients in the input data.

    Examples
    ---------
    >>> import viola
    >>> from io import StringIO
    >>> data1	=	\"\"\"chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	info_0	info_1
    ... chr1	10	11	chr1	20	21	test1	60	+	-	A	1
    ... chr1	150	151	chr1	300	301	test2	60	+	-	B	1
    ... chr3	10	11	chr4	20	21	test3	60	-	-	F	4
    ... \"\"\"
    >>> data2	=	\"\"\"chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	info_0	info_1
    ... chr1	100	101	chr1	250	251	test1	60	+	-	A	2
    ... chr2	180	181	chr2	3000	3001	test2	60	-	-	S	3
    ... chr2	10	11	chr5	20	21	test3	60	+	-	L	6
    ... \"\"\"
    >>> data3	=	\"\"\"chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	info_1	info_2
    ... chr2	10	11	chr5	20	21	test1	60	+	-	7	T
    ... chr3	10	11	chr4	20	21	test2	60	-	-	1	F
    ... \"\"\"
    >>> bedpe1 = viola.read_bedpe(StringIO(data1), patient_name='patient1')
    >>> bedpe2 = viola.read_bedpe(StringIO(data2), patient_name='patient2')
    >>> multi_bedpe1 = viola.MultiBedpe([bedpe1, bedpe2], ['patient1', 'patient2'])
    >>> bedpe3 = viola.read_bedpe(StringIO(data3), patient_name='patient3')
    >>> concat_bedpe = viola.concat([multi_bedpe1, bedpe3])
    >>> print(concat_bedpe)
    INFO=svlen,svtype,cipos,ciend,info_0,info_1,info_2
    Documentation of Bedpe object ==> https://dermasugita.github.io/ViolaDocs/docs/html/reference/multi_bedpe.html
                   id       be1        be2 strand  qual svtype
    0  patient1_test1   chr1:11    chr1:21     +-    60    DEL
    1  patient1_test2  chr1:151   chr1:301     +-    60    DEL
    2  patient1_test3   chr3:11    chr4:21     --    60    BND
    3  patient2_test1  chr1:101   chr1:251     +-    60    DEL
    4  patient2_test2  chr2:181  chr2:3001     --    60    INV
    5  patient2_test3   chr2:11    chr5:21     +-    60    BND
    6  patient3_test1   chr2:11    chr5:21     +-    60    BND
    7  patient3_test2   chr3:11    chr4:21     --    60    BND
    """
    return concat_bedpe(ls_objs)
    
def select_Maxwell_2011_d220_N1_laminaCurv_date(in1):
    '''    
    '''
    if in1 == 'Date_1':
        out1 = 271
	out2 ='c:\\openaleapkg\\adel\\adel\\data\Maxwell_2011_d220_N1_laminaCurv_date1.RData'
    elif in1 == 'Date_2':
        out1 = 466
	out2 ='c:\\openaleapkg\\adel\\adel\\data\\Maxwell_2011_d220_N1_laminaCurv_date2.RData'
    elif in1 == 'Date_3':
        out1 = 637
	out2 ='c:\\openaleapkg\\adel\\adel\\data\\Maxwell_2011_d220_N1_laminaCurv_date3.RData'
    elif in1 == 'Date_4':
        out1 = 922
	out2 ='c:\\openaleapkg\\adel\\adel\\data\\Maxwell_2011_d220_N1_laminaCurv_date4.RData'
    elif in1 == 'Date_5':
        out1 = 1125
	out2 ='c:\\openaleapkg\\adel\\adel\\data\\Maxwell_2011_d220_N1_laminaCurv_date5.RData'
    elif in1 == 'Date_6':
        out1 = 1428
	out2 ='c:\\openaleapkg\\adel\\adel\\data\\Maxwell_2011_d220_N1_laminaCurv_date6.RData'
    else: 
        out1 = 0
	out2 ='ERREUR'
    return out1, out2,

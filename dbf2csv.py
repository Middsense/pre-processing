# Myles Stokowski
# 07.03.2020
# quick .dbf to .csv conversion

import sys
import csv
from dbfread import DBF

def dbf2csv(indbf, outcsv):
    table = DBF(indbf)

    with open(outcsv, 'w', newline = '') as f:
        writer = csv.writer(f)
        writer.writerow(table.field_names)  
        for record in table:
            writer.writerow(list(record.values()))
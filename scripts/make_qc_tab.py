#!/usr/bin/env python
import sys

# Parse through the output log from hisat2 in SE mode
#outfile_h = open('../logs/alnment_table.tex'

def make_row(row):
    '''replace spaces with single & for latex table'''
    final_row = ' & '.join(row.strip().split('\t')) + '\\\\'
    return final_row

def escape_underscore(row):
    row = row.replace('_', '\_')
    return row

def escape_percent(row):
    row = row.replace('%', '\%')
    return row


header = '''
\\begin{table}[H]
\\begin{tabular}{ l c c }
\hline
\\textbf{Metric} & \\textbf{Reads} & \\textbf{Bases} \\\\
   '''
footer = '''
\hline
\end{tabular}
\caption{Summary of read QC results.}
\label{table:qc}
\end{table}
    '''

result = []
files = []

with open(sys.argv[1]) as f:
    for row in f:
        if row.find('in=') >= 0:
            seqf = row[row.find('in=')+13:].split(',')[0]
            if len(seqf.split(' ')) > 1:
                seqf = seqf.split(' ')[0]
            if seqf not in files:
                files.append(seqf)
                seqf = escape_underscore(seqf)
                seqf = '\hline\n' + seqf + '& & \\\\\n\hline'
                result.append(seqf)
        if row.find('Result:') >= 0:
            row = make_row(row)
            row = escape_percent(row)
            result.append(row)
        if row.find('QTrimmed:') >= 0:
            row = make_row(row)
            row = escape_percent(row)
            result.append(row)
        if row.find('KTrimmed:') >= 0:
            row = make_row(row)
            row = escape_percent(row)
            result.append(row)
        if row.find('Total Removed') >= 0:
            row = make_row(row)
            row = escape_percent(row)
            result.append(row)
        if row.find('Results:') >= 0:
            row = make_row(row)
            row = escape_percent(row)
            result.append(row)
print(header)
for row in result:
    print(row)
print(footer)

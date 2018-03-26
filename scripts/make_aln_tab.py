#!/usr/bin/env python
import sys

# Parse through the output log from hisat2

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
\\begin{tabular}{ | l | r |}
   '''
footer = '''
\hline
\end{tabular}
\caption{Summary of read alignment results.}
\label{table:qc}
\end{table}
    '''

result = []

with open(sys.argv[1]) as f:
    for row in f:
        if row.find('bams') >= 0:
            seqf = escape_underscore(row.strip())
            seqf = '\hline\n' + seqf[5:] + '&\\\\\n\hline'
            result.append(seqf)
        else:
            row = escape_percent(row)
            row = '&' + row.strip() + '\\\\'
            result.append(row)
print(header)
for row in result:
    print(row)
print(footer)

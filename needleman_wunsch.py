'''
Needleman Wunsch Algorithm
This algorithm uses the following scores:
Match=+1
Mismatch=-1
Gap=-1

There are 3 steps:
1. Intialization
2. Matrix Filling
3. Backtracking

This code has been developed without importing numpy and is an easy to understand approach for biologists new to coding!
'''

#!/usr/bin/env python3

#NOTE: A gap is represented by space for this code. For example, a sequence with a gap in second position can be represented as A TTGC.
#NOTE: The scores set for this are match=+1,mismatch=-1 and gap=-1
import sys
seq1=sys.argv[1]
seq2=sys.argv[2]

def nw_match(row_no1,col_no1,matrix):
    diagonal_val=matrix[row_no1-1][col_no1-1]
    diagonal_val=diagonal_val+1 #giving match score for gap
    top_val=matrix[row_no1-1][col_no1]
    top_val=top_val-1 #giving gap score for top
    left_val=matrix[row_no1][col_no1-1]
    left_val=left_val-1 #giving gap score for left
    value1=max(diagonal_val,top_val,left_val)
    if diagonal_val >= top_val and diagonal_val >= left_val:
        value4="diagonal"
    elif top_val >=left_val:
        value4="top"
    elif left_val > top_val:
        value4="left"
    return value1,value4

def nw_mismatch(row_no,col_no,matrix): #not assigning individual scores as gap score = mismatch score = -1
    diagonal_val=matrix[row_no-1][col_no-1]
    top_val=matrix[row_no-1][col_no]
    left_val=matrix[row_no][col_no-1]
    max_val=max(diagonal_val,top_val,left_val)
    value2=max_val-1
    value2=int(value2)
    if diagonal_val >= top_val and diagonal_val >= left_val:
        value5="diagonal"
    elif top_val >= left_val:
        value5="top"
    elif left_val > top_val:
        value5="left"
    return value2,value5

def nw_gap(row_no,col_no,matrix): #not assigning individual score as gap score = mismatch score = -1
    diagonal_val = matrix[row_no - 1][col_no - 1]
    top_val = matrix[row_no - 1][col_no]
    left_val = matrix[row_no][col_no - 1]
    max_val = max(diagonal_val, top_val, left_val)
    value3 = max_val - 1
    value3=int(value3)
    if diagonal_val >= top_val and diagonal_val >= left_val:
        value6="diagonal"
    elif top_val >= left_val:
        value6="top"
    elif left_val > top_val:
        value6="left"
    return value3,value6

#INITIALIZATION
#extracting only the sequences from the fasta files as 2 variables sequence_1 and sequence_2
with open(seq1, 'r') as file:
    seq=[line.strip() for line in file]
    sequence_1=seq[1]

with open(seq2, 'r') as file:
    seq=[line.strip() for line in file]
    sequence_2=seq[1]
#making a matrix with sequence_1 as the first column and sequence_2 as the first row
rows=len(sequence_1)+1
columns=len(sequence_2)+1
matrix= [0] * rows
for x in range(rows):
    matrix[x] = [0] * columns
arrow_matrix=[""] * rows #to keep a track of arrows for traceback
for x in range(rows):
    arrow_matrix[x]=[""] *columns
sequence_1=" "+sequence_1
sequence_2=" "+sequence_2
match=1 #variable not used in the code, only for us to remember values while coding
mismatch=-1 #variable not used in the code, only for us to remember values while coding
gap=-1
#MATRIX FILLING

for z in range(columns):
    matrix[0][0] = 0
    matrix[0][z]=matrix[0][z-1] + gap
    arrow_matrix[0][z]="left"
for z in range(rows):
    matrix[0][0]=0
    matrix[z][0]=matrix[z-1][0] + gap
    arrow_matrix[z][0]="top"

y=1
z=1
for y in range(1,rows):
    nucl_seq1=sequence_1[y]
    for z in range(1,columns):
        nucl_seq2=sequence_2[z]
        if nucl_seq1 == nucl_seq2 and nucl_seq1 !=" " and nucl_seq2 !=" ": #match
            row_no1=y
            col_no1=z
            value1,value4= nw_match(row_no1,col_no1,matrix)
            matrix[y][z]=value1
            arrow_matrix[y][z]=value4
        elif nucl_seq1 != nucl_seq2 and nucl_seq1 !=" " and nucl_seq2 !=" ": #mismatch
            row_no=y
            col_no=z
            value2,value5=nw_mismatch(row_no,col_no,matrix)
            matrix[row_no][col_no]=value2
            arrow_matrix[row_no][col_no]=value5
        elif nucl_seq1 == nucl_seq2 and nucl_seq1 == " " or nucl_seq2 == " ": #gap
            row_no=y
            col_no=z
            value3,value6=nw_gap(row_no,col_no,matrix)
            matrix[row_no][col_no]=value3
            arrow_matrix[row_no][col_no]=value6

#BACK TRACKING
seq1_out=""
seq2_out=""
alignment=""
current_rowno=rows-1
current_nuclno_1=current_rowno
current_colno=columns-1
current_nuclno_2=current_colno
score=matrix[current_rowno][current_colno]
while current_rowno >0 and current_colno >0:
    arrow=arrow_matrix[current_rowno][current_colno]
    current_cell = matrix[current_rowno][current_colno]
    diagonal_cell = matrix[current_rowno-1][current_colno-1]
    top_cell = matrix[current_rowno-1][current_colno]
    left_cell = matrix[current_rowno][current_colno-1]
    if arrow == "diagonal" and sequence_1[current_nuclno_1] == sequence_2[current_nuclno_2]: #tracking diagonal arrow when bases are same
        seq1_out=sequence_1[current_nuclno_1]+seq1_out
        seq2_out=sequence_2[current_nuclno_2]+seq2_out
        alignment="|"+alignment
        current_cell=matrix[current_rowno-1][current_colno-1]
        current_rowno=current_rowno-1
        current_colno=current_colno-1
        current_nuclno_1 = current_nuclno_1 - 1
        current_nuclno_2 = current_nuclno_2 - 1
    elif arrow == "diagonal" and sequence_1[current_nuclno_1] != sequence_2[current_nuclno_2]: #tracking diagonal arrow when bases are different
        seq1_out = sequence_1[current_nuclno_1] + seq1_out
        seq2_out = sequence_2[current_nuclno_2] + seq2_out
        alignment = "*" + alignment
        current_cell = matrix[current_rowno - 1][current_colno - 1]
        current_rowno = current_rowno - 1
        current_colno = current_colno - 1
        current_nuclno_1 = current_nuclno_1 - 1
        current_nuclno_2 = current_nuclno_2 - 1
    elif arrow == "top": #tracking top arrow
        seq1_out = sequence_1[current_nuclno_1] + seq1_out
        seq2_out="-"+seq2_out
        alignment=" "+alignment
        current_cell=matrix[current_rowno-1][current_colno]
        current_rowno=current_rowno-1
        current_nuclno_1 = current_nuclno_1 - 1
    elif arrow == "left": #tracking left arrow
        seq1_out="-"+seq1_out
        seq2_out = sequence_2[current_nuclno_2] + seq2_out
        alignment=" "+alignment
        current_cell=matrix[current_rowno][current_colno-1]
        current_colno=current_colno-1
        current_nuclno_2 = current_nuclno_2 - 1

print(seq1_out)
print(alignment)
print(seq2_out)
print("Alignment score:%d" %score)

'''
Smith Watermann Algorithm - Local Alignment
Easy to understand for biologists new to coding!
This code has been developed without using numpy!
'''


#!/usr/bin/env python3

#NOTE: A gap is represented by space for this code. For example, a sequence with a gap in second position can be represented as A TTGC.
#NOTE: The scores set for this are match=+1,mismatch=-1 and gap=-1
import sys
seq1=sys.argv[1]
seq2=sys.argv[2]

def sw_match(row_no,col_no,matrix):
    diagonal_val = matrix[row_no - 1][col_no - 1]
    diagonal_val = diagonal_val + 1  # adding match points
    top_val = matrix[row_no - 1][col_no]
    top_val = top_val - 1  # subtracting gap points
    left_val = matrix[row_no][col_no - 1]
    left_val = left_val - 1  # subtracting gap points
    max_val = max(diagonal_val, top_val, left_val)
    if max_val > 0:
        value1 = max_val
    else:
        value1 = 0
    if diagonal_val >= top_val and diagonal_val >= left_val:
        value4 = "diagonal"
    elif top_val >= left_val:
        value4 = "top"
    elif left_val > top_val:
        value4 = "left"
    return value1, value4

def sw_mismatch(row_no,col_no,matrix):
    diagonal_val=matrix[row_no-1][col_no-1]
    diagonal_val=diagonal_val-1 #subtracting mismatch points
    top_val=matrix[row_no-1][col_no]
    top_val=top_val-1 #subtracting gap points
    left_val=matrix[row_no][col_no-1]
    left_val=left_val-1 #subtracting gap points
    max_val=max(diagonal_val,top_val,left_val)
    if max_val > 0:
        value2=max_val
    else:
        value2=0
    if diagonal_val >= top_val and diagonal_val >= left_val:
        value5="diagonal"
    elif top_val >= left_val:
        value5="top"
    elif left_val > top_val:
        value5="left"
    return value2,value5

def sw_gap(row_no,col_no,matrix):
    diagonal_val = matrix[row_no - 1][col_no - 1]
    diagonal_val=diagonal_val-1 #subtracting gap points
    top_val = matrix[row_no - 1][col_no]
    top_val=top_val-1 #subtracting gap points
    left_val = matrix[row_no][col_no - 1]
    left_val=left_val-1 #subtracting gap points
    max_val = max(diagonal_val, top_val, left_val)
    if max_val > 0:
        value3 = max_val
    else:
        value3 = 0
    if diagonal_val >= top_val and diagonal_val >= left_val:
        value6 = "diagonal"
    elif top_val >= left_val:
        value6 = "top"
    elif left_val > top_val:
        value6 = "left"
    return value3, value6


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
match=3 #variable not used in the code, only for us to remember values while coding
mismatch=-3 #variable not used in the code, only for us to remember values while coding
gap=-2
#MATRIX FILLING
for z in range(columns):
    matrix[0][0] = 0
    flag=matrix[0][z-1] + gap
    if flag > 0:
        matrix[0][z]=matrix[0][z-1] + gap
    else:
        matrix[0][z] = 0
    arrow_matrix[0][z]="left"
for z in range(rows):
    matrix[0][0]=0
    flag=matrix[z-1][0] + gap
    if flag > 0:
        matrix[z][0]=matrix[z-1][0] + gap
    else:
        matrix[z][0] = 0
    arrow_matrix[z][0]="top"
y=1
z=1
for y in range(1,rows):
    nucl_seq1=sequence_1[y]
    for z in range(1,columns):
        nucl_seq2=sequence_2[z]
        if nucl_seq1 == nucl_seq2 and nucl_seq1 != " " and nucl_seq2 != " ": #match
            row_no=y
            col_no=z
            value1, value4=sw_match(row_no, col_no, matrix)
            matrix[y][z]=value1
            arrow_matrix[y][z]=value4
        elif nucl_seq1 != nucl_seq2 and nucl_seq1 !=" " and nucl_seq2 != " ": #mismatch
            row_no = y
            col_no = z
            value2, value5 = sw_mismatch(row_no, col_no, matrix)
            matrix[y][z] = value2
            arrow_matrix[y][z] = value5
        elif nucl_seq1 != nucl_seq2 and nucl_seq1 == " " or nucl_seq2 == " ": #gap
            row_no = y
            col_no = z
            value3, value6 = sw_gap(row_no, col_no, matrix)
            matrix[y][z] = value3
            arrow_matrix[y][z] = value6

#BACKTRACKING
seq1_out=""
seq2_out=""
alignment=""
highest=0
#finding score
for x in range(rows):
    for y in range(columns):
        if matrix[x][y] >= highest:
            highest=matrix[x][y]
            score_row=x
            score_col=y
current_rowno=score_row
current_colno=score_col
current_nuclno_1=score_row
current_nuclno_2=score_col
while current_rowno > 0 and current_colno > 0:
    arrow = arrow_matrix[current_rowno][current_colno]
    current_cell = matrix[current_rowno][current_colno]
    diagonal_cell = matrix[current_rowno - 1][current_colno - 1]
    top_cell = matrix[current_rowno - 1][current_colno]
    left_cell = matrix[current_rowno][current_colno - 1]
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
seq1_out=seq1_out[1:] #removing the first character as its value is 0
seq2_out=seq2_out[1:] #removing the first character as its value is 0
alignment=alignment[1:]
print(seq1_out)
print(alignment)
print(seq2_out)
print("Alignment score:%d"%highest)

# Needleman-Wunsch
# This program takes in a single FASTA file to search and prompts the user for a query sequence.
# This program returns a list of full matches to the command line, along with their start index.


# Function to handle scoring character comparisons
def score(char1, char2):
    if char1 == char2:
        return 2
    else:
        return -1


# Function to fill in any given grid space
def calculate_score(table, row, column, seq1, seq2):
    if row == 0:
        value = column*-1
        move = 6 # Move left
    elif column == 0:
        value = row*-1
        move = 4 # Move up
    else:
        char1 = seq1[row-1]
        char2 = seq2[column-1]
        a = table[row-1][column-1] + score(char1, char2)
        b = table[row-1][column] - 1
        c = table[row][column - 1] - 1
        if a > b and a > c:
            value = a
            # Diagonal move only
            move = 0
        elif a > b and a == c:
            value = a
            # Diagonal and left move
            move = 1
        elif a > c and a == b:
            value = a
            # Diagonal and vertical move
            move = 2
        elif a == b and a == c:
            value = a
            # Diagonal and left and vertical move
            move = 3
        elif b > a and b > c:
            value = b
            # Vertical move only
            move = 4
        elif b > a and b == c:
            value = b
            # Left and vertical move
            move = 5
        else:
            value = c
            # Vertical move only
            move = 6
    return value, move


# The Needleman-Wunsch algorithm
def needleman_wunsch(seq1, seq2):
    score_table = []
    move_table = []
    for i in range(len(seq1)+1):
        score_table.append([])
        move_table.append([])
        for j in range(len(seq2)+1):
            scores = calculate_score(score_table, i, j, seq1, seq2)
            score_table[i].append(scores[0])
            move_table[i].append(scores[1])
            # Loops through table filling in values
    return score_table[len(seq1)][len(seq2)], move_table

def get_alighnment(table, seq1, seq2):
    alighnment_seq1 = ""
    alighnment_seq2 = ""
    i, j = len(seq1), len(seq2)
    while i != 0  or j != 0:
        x = table[i][j]
        if x == 0:
            alighnment_seq1 += seq1[i-1]
            alighnment_seq2 += seq2[j-1]
            i -= 1
            j -= 1
        elif x == 1:
            alighnment_seq1 += seq1[i-1]
            alighnment_seq2 += seq2[j-1]
            i -= 1
            j -= 1
        elif x == 2:
            alighnment_seq1 += seq1[i-1]
            alighnment_seq2 += seq2[j-1]
            i -= 1
            j -= 1
        elif x == 3:
            alighnment_seq1 += seq1[i-1]
            alighnment_seq2 += seq2[j-1]
            i -= 1
            j -= 1
        elif x == 4:
            alighnment_seq1 += seq1[i-1]
            alighnment_seq2 += '_'
            i -= 1
        elif x == 5:
            alighnment_seq1 += seq1[i-1]
            alighnment_seq2 += '_'
            i -= 1
        else: # Case 6
            alighnment_seq1 += '_'
            alighnment_seq2 += seq2[j-1]
            j -= 1
    return alighnment_seq1, alighnment_seq2



def main():
    """
    Prompts user for file path and reads in FASTA sequences.
    Calls Needleman-Wunsch.
    Prints the output alignments.
    """
    # Handel the file opening of the Fasta file to be searched
    file = open(input("Please enter a file path: "))
    # Create an empty list for all sequences
    sequences = []
    # Initialize temp variables for sequence concatenation, and loop control
    seq = ""
    next_seq = False
    # Read in and concatenate sequence lines from Fasta file
    for line in file:
        # If header detected and seq not trivial save seq to sequences list
        # Else append line to seq
        if line.find('>') != -1:
            if not next_seq:
                next_seq = True
                continue
            else:
                sequences.append(seq)
                seq = ""
        else:
            line = line.strip().lower()
            seq += line
    # Append seq to sequences list
    sequences.append(seq)
    # Call algorithm
    result = needleman_wunsch(sequences[0], sequences[1])
    max_score = result[0]
    alighnments = get_alighnment(result[1], sequences[0], sequences[1])
    print("Max score: " + str(max_score))
    alighned1 = alighnments[0][::-1]
    alighned2 = alighnments[1][::-1]
    compare_seq = ""
    for i in range(len(alighned1)):
        if alighned1[i] == alighned2[i]:
            compare_seq += '|'
        else:
            compare_seq += '*'
    print(alighned1)
    print(compare_seq)
    print(alighned2)


# Call main
main()

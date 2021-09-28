# Needleman-Wunsch
# This program takes in a single FASTA file to search and prompts the user for a query sequence.
# This program returns a list of full matches to the command line, along with their start index.


# Function to handle scoring character comparisons
def score(char1, char2):
    if char1 == '_' and char2 == '_':
        raise NameError("Matching Blanks!")
    elif char1 == char2:
        return 2
    else:
        return -1


# Function to fill in any given grid space
def calculate_score(table, row, column, seq1, seq2):
    if row == 0:
        return column*-1
    elif column == 0:
        return row*-1
    else:
        a = table[row-1][column-1]
        b = table[row-1][column]
        c = table[row][column - 1]
        char1 = seq1[row]
        char2 = seq2[column]
        if a > b and a > c:
            return a + score(char1, char2)
            # Diagonal move only
            # Put code to fill out path table here
        elif a > b and a == c:
            return a + score(char1, char2)
            # Diagonal and left move
            # Put code to fill out path table here
        elif a > c and a == b:
            return a + score(char1, char2)
            # Diagonal and vertical move
            # Put code to fill out path table here
        elif a == b and a == c:
            return a + score(char1, char2)
            # Diagonal and left and vertical move
            # Put code to fill out path table here
        elif b > a and b > c:
            return b + score(char1, char2)
            # Vertical move only
            # Put code to fill out path table here
        elif b > a and b == c:
            return b + score(char1, char2)
            # Left and vertical move
            # Put code to fill out path table here
        else:
            return c + score(char1, char2)
            # Vertical move only
            # Put code to fill out path table here


# The Needleman-Wunsch algorithm
def needleman_wunsch(seq1, seq2):
    table = []
    for i in range(len(seq1)):
        table.append([])
        for j in range(len(seq2)):
            table[i].append(calculate_score(table, i, j, seq1, seq2))
            # Loops through table filling in values
    return table[len(seq1)-1][len(seq2)-1]


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
    print(needleman_wunsch(sequences[0], sequences[1]))


# Call main
main()

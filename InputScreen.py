# Import the required libraries
from tkinter import ttk
from tkinter import *
import numpy as np
import DnaAlignment
import ProteinAlignment


def get_PAM_Matrix():
    try:
        data = open("PAM.txt", 'r')
        dimensions = data.readline()
        try:
            n = int(dimensions)
        except ValueError:
            print("Wrong file")
            exit()
        letters = data.readline()
        letters = letters.replace('\n', '')
        letters_arr = letters.split('\t')
        score_matrix = np.zeros((n, n))
        for i in range(0, n):
            arr = data.readline().split("\t")
            for j in range(0, n):
                score_matrix[i][j] = float(arr[j])
        return score_matrix
    except FileNotFoundError:
        print("There is no matrix file.")
        exit()


def draw_matrix(Protein_LF):
    matrix = get_PAM_Matrix()
    f = Frame(Protein_LF)
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix)):
            Label(f, text=matrix[i, j], borderwidth=2).grid(row=i, column=j, padx=1, pady=1)
    f.place(x=600, y=50)


def protein_button_click(seq1, seq2, gap, window):
    ProteinAlignment.get_argments_Protein(seq1, seq2, gap,window)

def protein_LF(window, Protein_LF):
    Label(Protein_LF, text="Hint:", foreground="#800080", font=('Georgia', 10, "bold")).place(x=25, y=50, height=25)

    Label(Protein_LF, text="Use nucleotide [ A ,C ,G, T ] only ", foreground="#800080", font=("Arial", 10)).place(x=100,y=50,height=25)
    seq1 = StringVar()
    Label(Protein_LF, text="First Sequence", font=("Aerial", 10)).place(x=25, y=100, height=25)
    TOP_SEQ = Entry(Protein_LF, textvariable=seq1).place(x=200, y=100, width=200, height=25)
    seq2 = StringVar()

    Label(Protein_LF, text="Second Sequence", font=("Aerial", 10)).place(x=25, y=150, height=25)
    BOTTOM_SEQ = Entry(Protein_LF, textvariable=seq2).place(x=200, y=150, width=200, height=25)

    Label(Protein_LF, text="Scoring:-", font=("Aerial", 10)).place(x=25, y=200, height=25)

    gap = IntVar(value=-8)
    Label(Protein_LF, text="gap", font=("Aerial", 10)).place(x=200, y=200, height=25)

    g = Spinbox(Protein_LF, from_=-10, to=0, textvariable=gap).place(x=250, y=200, height=25, width=50)

    Button(Protein_LF, text="Result", background="#800080", foreground="white", cursor="hand2",
           command=lambda: protein_button_click(seq1.get(), seq2.get(), gap.get(), window)).place(x=200, y=300,
                                                                                                  height=50, width=100)


def dna_button_click(seq1, seq2, match, mismatch, gap, window):
    DnaAlignment.get_argments_DNA(seq1, seq2, match, mismatch, gap,window)



def dna_LF(window, DNA_LF):
    Label(DNA_LF, text="Hint:", foreground="#800080", font=('Georgia', 10, "bold")).place(x=25, y=50, height=25)
    Label(DNA_LF, text="Use nucleotide [ A ,C ,G, T ] only ", foreground="#800080", font=("Arial", 10)).place(x=100,
                                                                                                              y=50,
                                                                                                              height=25)
    seq1 = StringVar()
    Label(DNA_LF, text="First Sequence", font=("Aerial", 10)).place(x=25, y=100, height=25)
    TOP_SEQ = Entry(DNA_LF, textvariable=seq1).place(x=200, y=100, width=200, height=25)
    seq2 = StringVar()
    Label(DNA_LF, text="Second Sequence", font=("Aerial", 10)).place(x=25, y=150, height=25)
    BOTTOM_SEQ = Entry(DNA_LF, textvariable=seq2).place(x=200, y=150, width=200, height=25)

    Label(DNA_LF, text="Scoring:-", font=("Aerial", 10)).place(x=25, y=200, height=25)
    match = IntVar(value=1)
    Label(DNA_LF, text="Match", font=("Aerial", 10)).place(x=200, y=200, height=25)
    m = Spinbox(DNA_LF, from_=1, to=15, font=('sans-serif', 14), textvariable=match).place(x=250, y=200, height=25,
                                                                                           width=50)
    mismatch = IntVar(value=-1)
    Label(DNA_LF, text="Mismatch", font=("Aerial", 10)).place(x=350, y=200, height=25)
    mis = Spinbox(DNA_LF, from_=-10, to=0, font=('sans-serif', 14), textvariable=mismatch).place(x=420, y=200,
                                                                                                 height=25, width=50)

    gap = IntVar(value=-2)
    Label(DNA_LF, text="Gap", font=("Aerial", 10)).place(x=515, y=200, height=25)
    g = Spinbox(DNA_LF, from_=-10, to=0, font=('sans-serif', 14), textvariable=gap).place(x=550, y=200, height=25,width=50)

    B1 = Button(DNA_LF, text="Result", background="#800080", foreground="white", cursor="hand2",
                command=lambda: dna_button_click(seq1.get(), seq2.get(), match.get(), mismatch.get(), gap.get(),window)).place(x=400, y=300, height=50, width=100)


# Define a function for switching the frames
def change_to_DNA(window, DNA_LF, Protein_LF):
    dna_LF(window, DNA_LF)
    DNA_LF.pack(fill='both', padx=50, pady=20)
    Protein_LF.pack_forget()


def change_to_PROTEIN(window, DNA_LF, Protein_LF):
    protein_LF(window, Protein_LF)
    draw_matrix(Protein_LF)
    Protein_LF.pack(fill='both', padx=50, pady=20)
    DNA_LF.pack_forget()


# Create fonts for making difference in the frame

def CREATE_WINDOW_FRAMES():
    window = Tk()
    # Set the size of the window
    window.geometry("1500x900")
    window.title("LOCAL ALIGNMENT APP")
    label = Label(window, text="LOCAL SEQUENCE ALIGNMEN", fg="#800080", font=("Georgia", 15, "bold"))
    label.pack(pady=30)
    frame = ttk.LabelFrame(window, height=100, width=1000, text="DNA OR PROTEIN")
    DNA_LF = ttk.LabelFrame(window, height=500, width=1000, text="DNA")
    Protein_LF = ttk.LabelFrame(window, height=500, width=1000, text="PROTEIN")
    return window, frame, DNA_LF, Protein_LF


def DNA_PROTEIN_LF():
    # Create an instance of tkinter frame or window
    window, frame, DNA_LF, Protein_LF = CREATE_WINDOW_FRAMES()
    var = IntVar()
    var.set(1)
    dna_LF(window, DNA_LF)
    DNA = ttk.Radiobutton(frame, text="DNA", value=1, variable=var,
                          command=lambda: change_to_DNA(window, DNA_LF, Protein_LF)).grid(row=1, column=1, padx=50)
    protein = ttk.Radiobutton(frame, text="PROTEIN", value=2, variable=var,
                              command=lambda: change_to_PROTEIN(window, DNA_LF, Protein_LF)).grid(row=1, column=2)
    frame.pack(fill='both', padx=50)
    DNA_LF.pack(fill='both', expand=1, padx=50, pady=20)
    window.mainloop()

# Import the required libraries
from tkinter import ttk
from tkinter import *
import InputScreen



def draw_matrix(matrix,align_locs,f,seq1,seq2):

    frame = ttk.LabelFrame(f, height=500, width=500, text="Matrix")

    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[0])):
            loc = (i, j)

            if (loc in align_locs):
                Label(frame, text=matrix[i, j], foreground="white", background="#800080", borderwidth=2).grid(row=i, column=j, padx=1, pady=1)
            else:
                Label(frame, text=matrix[i, j], background="white", borderwidth=5).grid(row=i,column=j,padx=1,pady=1)
    frame.pack(pady=10)

def back_window(window):
    window.destroy()
    InputScreen.DNA_PROTEIN_LF()


def output_information(matrix,align_locs,align1,align2,seq1,seq2,score):
    window = Tk()
    window.geometry()
    window.geometry("1500x900")
    f = ttk.LabelFrame(window, height=500, width=500, text="Alignment Output",)
    align1_label= Label(f, text="Alignment First sequence:- "+align1,foreground="#800080" ,font=('Georgia',10)).pack(pady=10,padx=50)
    align1_label= Label(f, text="Alignment Second Sequence:- " +align2,foreground="#800080" ,font=('Georgia',10)).pack(pady=10,padx=50)
    align1_label = Label(f, text="Score :- " + str(score), foreground="#800080", font=('Georgia',10)).pack(pady=10, padx=50)
    draw_matrix(matrix,align_locs,f,seq1,seq2)
    back = Button(f, text="Back", foreground="white",background="#800080", command=lambda :back_window(window)).pack()
    f.pack(fill="both", padx=50, pady=50)
    window.mainloop()

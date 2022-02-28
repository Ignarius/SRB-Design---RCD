################################################################################################
# If you can read python code, the code below is actual dogshit so don't judge it I know :P    #
# Created By: Ignarius                                                                         #
################################################################################################

from tkinter import *
from tkinter import ttk
from tkinter import messagebox
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.pyplot as plt
import pandas as pd

def solve_beam(maxMoment, baseBeam, ccBeamT, ccBeamS, f_cBeam, f_yBeam, depthBeamT, diameter):
    """Solves for the specifications of the beam

    Args:
        maxMoment (float): List of moments arranging from left, center, right
        baseBeam (float): Base of beam
        ccBeamT (float): Concrete cover from centroid to the bottom of beam
        ccBeamS (float): Concrete cover from centroid to the side of beam
        f_cBeam (float): Concrete Stress
        f_yBeam (float): Steel yield Stress
        depthBeamT (float): Depth of beam
        diameter (float): Diameter of steel bars

    Returns:
        float: Specifications of the beam
    """
    if f_cBeam <= 28:
        beta = 0.85
    elif f_cBeam < 55:
        beta = 0.85-(0.05*((f_cBeam-28)/7))
    else:
        beta = 0.65

    if ((2*np.abs(maxMoment)/(0.9*baseBeam*np.power(depthBeamT, 2))) / (0.85*f_cBeam)) > 1:
        messagebox.showerror("Moment is too big", "Increase Beam size!")
        return

    # Steel Ratios
    steel_ratio_max = 0.75*((0.85*beta*f_cBeam)/f_yBeam)*(600/(600+f_yBeam))
    steel_ratio_min1 = np.sqrt(f_cBeam)/(4*f_yBeam)
    steel_ratio_min2 = 1.4/f_yBeam
    steel_ratio = ((0.85*f_cBeam)/f_yBeam) * (1 - np.sqrt(1 - ((2*np.abs(maxMoment)/(0.9*baseBeam*np.power(depthBeamT, 2))) / (0.85*f_cBeam))))

    # Selecting the minimum value for steel ratio
    if steel_ratio_min1 > steel_ratio_min2:
        steel_ratio_min = steel_ratio_min1
    else:
        steel_ratio_min = steel_ratio_min2

    # Checking if steel ratio is less than min
    if steel_ratio <= steel_ratio_min:
        steel_ratio = steel_ratio_min

    if steel_ratio > steel_ratio_max:
        messagebox.showerror("ERROR", "BEAM IS DOUBLY REINFORCED!")
        return

    else:
        # Checking Spacing
        if diameter == "Specify Diameter":
            messagebox.showerror("Define Diameter", "Define a diameter!")
            return

        elif diameter == "Auto":
            for i in diameters:
                d_beam_T = i
                area_steel_T = steel_ratio*baseBeam*depthBeamT
                area_T = (np.pi/4)*np.power(i, 2)
                no_of_bars_T = int(area_steel_T / area_T) + 1
                if no_of_bars_T == 1:
                    messagebox.showerror("Error!", "Cannot Design Beam")
                    return
                else:
                    spacing_T = (baseBeam-2*ccBeamT-no_of_bars_T*i) / (no_of_bars_T-1)
                    if spacing_T > 25:
                        break
                    if spacing_T < 25 and diameter == 57.3:
                        messagebox.showerror("Error!", "Cannot Design Beam")
                        return

        else:
            d_beam_T = float(diameter)
            area_steel_T = steel_ratio*baseBeam*depthBeamT
            area_T = (np.pi/4)*np.power(d_beam_T, 2)
            no_of_bars_T = int(area_steel_T / area_T) + 1
            if no_of_bars_T == 1:
                no_of_bars_T = 2
            spacing_T = (baseBeam-2*ccBeamT-no_of_bars_T*d_beam_T) / (no_of_bars_T-1)
            if spacing_T < 25:
                messagebox.showerror("Check Diameter", "Steel Spacing is less than minimum!\nIncrease Steel diameter")
                return


    # Capacity of Beam
    steel_area_T = no_of_bars_T*(np.pi/4)*np.power(d_beam_T, 2)
    a_val = (steel_area_T*f_yBeam) / (0.85*f_cBeam*baseBeam)
    c_val = a_val / beta

    deflection_T = 0.003*((depthBeamT-c_val) / c_val)
    if deflection_T >= 0.005:
        red_factor = 0.9
    elif deflection_T < (f_yBeam / 200000):
        red_factor = 0.65
    else:
        red_factor = 0.65 + (deflection_T-(f_yBeam / 200000)) * (0.25/(0.005-(f_yBeam / 200000)))

    max_capacity = 0.85*f_cBeam*beta*c_val*baseBeam*(depthBeamT-((beta*c_val)/2))
    capacity = red_factor*max_capacity

    return d_beam_T, steel_ratio_max, steel_ratio_min1, steel_ratio_min2, steel_ratio, area_steel_T, no_of_bars_T, spacing_T, capacity, beta, steel_ratio_min


def solve_stirrups(max_shear, mid_shear, lengthBeam, baseBeam, depthBeamT, f_cBeam, f_ytBeam, lambda_concrete, stirrup_d):
    maximum_shear = max_shear * np.power(10, 3)
    middle_shear = mid_shear * np.power(10, 3)
    mid_length = lengthBeam / 2

    V_c = 0.75*0.17*lambda_concrete*np.sqrt(f_cBeam)*baseBeam*depthBeamT
    if (1/2)*V_c < middle_shear:
        dist_of_stirrup = mid_length
        centerline = "Yes"
    else:
        try:
            # From the center measurement
            cent_dist_of_stirrup = ((1/2)*V_c-middle_shear)*(mid_length / (maximum_shear-middle_shear))
            # From the left measurement
            dist_of_stirrup = mid_length - cent_dist_of_stirrup
            dist_of_stirrup = int(dist_of_stirrup / 5) * 5

        except OverflowError:
            dist_of_stirrup = 0
        centerline = "No"

    V_u = (((maximum_shear-middle_shear) / mid_length) * (mid_length-depthBeamT)) + middle_shear
    V_s = (V_u-V_c) / 0.75

    if V_s > 0.66*np.sqrt(f_cBeam)*baseBeam*depthBeamT:
        messagebox.showerror("Shear is too big", "Increase Beam Size!")
        return

    # S_max solution
    if V_s <= 0.33*np.sqrt(f_cBeam)*baseBeam*depthBeamT:
        S_max_1 = depthBeamT / 2
        S_max_2 = 600
    else:
        S_max_1 = depthBeamT / 4
        S_max_2 = 300

    if S_max_1 < S_max_2:
        S_max = int(S_max_1 / 5) * 5 # Change to 25
    else:
        S_max = int(S_max_2 / 5) * 5 # Change to 25

    # A_v solution
    A_v_min1 = 0.062*np.sqrt(f_cBeam)*((baseBeam*S_max) / f_ytBeam)
    A_v_min2 = 0.35*((baseBeam*S_max) / f_ytBeam)

    if A_v_min1 < A_v_min2:
        A_v_min = A_v_min2
    else:
        A_v_min = A_v_min1

    if stirrup_d == "Specify Diameter":
        messagebox.showerror("Define Diameter", "Define a diameter!")
        return

    elif stirrup_d == "Auto":
        for i in diameters:
            d_stirrup = i
            A_v = 2*(np.pi/4)*np.power(d_stirrup, 2)
            if A_v > A_v_min:
                break
    else:
        d_stirrup = float(stirrup_d)
        A_v = 2*(np.pi/4)*np.power(d_stirrup, 2)
        if A_v < A_v_min:
            messagebox.showerror("Stirrup Diameter too small", "Increase Stirrup Diameter!")
            return

    S_value = (A_v*f_ytBeam*depthBeamT) / (V_s)

    if V_s < 0:
        S_value = S_max
    else:
        if S_value > S_max:
            S_value = S_max
        else:
            S_value = int(S_value / 25) * 25

    return V_c, centerline, V_u, V_s, S_max, A_v_min, S_value, dist_of_stirrup, A_v, d_stirrup


def draw_beam(array_num):
    """Draws the cross sectional area of the beam

    Args:
        array_num (int): Beam number
    """
    for widget in wrapper3.winfo_children():
        widget.destroy()

    fig = Figure(figsize = (5, 10), dpi = 100)
    plot = [fig.add_subplot(131, aspect=1), fig.add_subplot(132, aspect=1), fig.add_subplot(133, aspect=1)]
    plot[0].axis("off")
    plot[1].axis("off")
    plot[2].axis("off")
    
    # Beam Exterior Drawing
    place = ['Left', 'Center', 'Right']
    for i in range(3):
        if place[i] == 'Center':
            # Center Drawing
            bot_spacing = (data_values[array_num][i][1] - (2*data_values[array_num][i][4])) / (data_values[array_num][i][14]-1)
            plot[i].plot([0, 0], [0, data_values[array_num][i][2]], 'k')
            plot[i].plot([data_values[array_num][i][1], data_values[array_num][i][1]], [0, data_values[array_num][i][2]], 'k')
            plot[i].plot([0, data_values[array_num][i][1]], [0, 0], 'k')
            plot[i].plot([0, data_values[array_num][i][1]], [data_values[array_num][i][2], data_values[array_num][i][2]], 'k')

            # Add Circles
            for ii in range(int(data_values[array_num][i][14])):
                circle = plt.Circle((data_values[array_num][i][4]+(ii*bot_spacing), data_values[array_num][i][2]-data_values[array_num][i][7]), data_values[array_num][i][9] / 2, color='k')
                plot[i].add_patch(circle)

            # Add Text
            plot[i].text(data_values[array_num][i][1]/2, data_values[array_num][i][2], place[i], style='oblique', ha='center', wrap=True, va='bottom')

        else:
            # Left and right Drawings
            bot_spacing = (data_values[array_num][i][1] - (2*data_values[array_num][i][4])) / (data_values[array_num][i][14]-1)
            plot[i].plot([0, 0], [0, data_values[array_num][i][2]], 'k')
            plot[i].plot([data_values[array_num][i][1], data_values[array_num][i][1]], [0, data_values[array_num][i][2]], 'k')
            plot[i].plot([0, data_values[array_num][i][1]], [0, 0], 'k')
            plot[i].plot([0, data_values[array_num][i][1]], [data_values[array_num][i][2], data_values[array_num][i][2]], 'k')

            # Add Circles
            for ii in range(int(data_values[array_num][i][14])):
                circle = plt.Circle((data_values[array_num][i][4]+(ii*bot_spacing), data_values[array_num][i][2]-data_values[array_num][i][3]), data_values[array_num][i][9] / 2, color='k')
                plot[i].add_patch(circle)

            # Add Text
            plot[i].text(data_values[array_num][i][1]/2, data_values[array_num][i][2], place[i], style='oblique', ha='center', wrap=True, va='bottom')

    canvas = FigureCanvasTkAgg(fig, master = wrapper3)
    canvas.draw()
    canvas.get_tk_widget().place(x=10, y=10)
    toolbar = NavigationToolbar2Tk(canvas, wrapper3)
    toolbar.update()
    canvas.get_tk_widget().pack()


def draw_stirrups(array_num):
    for widget in wrapper4.winfo_children():
        widget.destroy()

    length = data_values[array_num][3][6]
    height = data_values[array_num][1][2]
    cc = data_values[array_num][1][3]
    depth = data_values[array_num][1][7]
    dist_of_stirrup = data_values[array_num][3][7]
    S_value = data_values[array_num][3][5]
    V_c = data_values[array_num][3][0]
    V_u = data_values[array_num][3][1]

    no_of_stirrups = int((dist_of_stirrup - depth) / S_value) # Add 50

    fig = Figure(figsize = (5, 10), dpi = 100)
    plot = fig.add_subplot(111, aspect=1)
    plot.axis("off")

    plot.plot([0, length / 2], [0, 0], 'k')
    plot.plot([0, length / 2], [height, height], 'k')
    plot.plot([0, 0], [0, height], 'k')
    plot.plot([length / 2, length / 2], [0, (3/2)*height], '--k')
    plot.plot([0, length / 2], [cc, cc], 'k')
    plot.plot([0, length / 2], [height-cc, height-cc], 'k')

    plot.plot([50, 50], [cc, height-cc], 'k')
    if (50+depth) < (length/2):
        plot.plot([50+depth, 50+depth], [cc, height-cc], 'k')

    if V_u > V_c:
        plot.plot([dist_of_stirrup, dist_of_stirrup], [cc, height-cc], 'k')
        for ii in range(no_of_stirrups):
            plot.plot([depth+50+(S_value*(ii+1)), depth+50+(S_value*(ii+1))], [cc, height-cc], 'k')

    plot.text(length/2+50, height/2, "C. L.", style='oblique', ha='left', va='center', wrap=True)
    
    canvas = FigureCanvasTkAgg(fig, master = wrapper4)
    canvas.draw()
    canvas.get_tk_widget().place(x=10, y=10)
    toolbar = NavigationToolbar2Tk(canvas, wrapper4)
    toolbar.update()
    canvas.get_tk_widget().pack()


def update_solved(array_num):
    """Updates the labels used for the specification

    Args:
        array_num (int): Beam number
    """
    label_1.config(text=data_values[array_num][0][9])
    label_2.config(text=data_values[array_num][0][10])
    label_3.config(text=data_values[array_num][0][-2])
    label_4.config(text=data_values[array_num][0][-1])
    label_5.config(text=data_values[array_num][0][12])
    label_6.config(text=data_values[array_num][0][13])
    label_7.config(text=data_values[array_num][0][14])
    label_8.config(text=data_values[array_num][0][15])
    label_9.config(text=data_values[array_num][0][16])
    
    label_10.config(text=data_values[array_num][1][9])
    label_11.config(text=data_values[array_num][1][10])
    label_12.config(text=data_values[array_num][1][-2])
    label_13.config(text=data_values[array_num][1][-1])
    label_14.config(text=data_values[array_num][1][12])
    label_15.config(text=data_values[array_num][1][13])
    label_16.config(text=data_values[array_num][1][14])
    label_17.config(text=data_values[array_num][1][15])
    label_18.config(text=data_values[array_num][1][16])

    label_19.config(text=data_values[array_num][2][9])
    label_20.config(text=data_values[array_num][2][10])
    label_21.config(text=data_values[array_num][2][-2])
    label_22.config(text=data_values[array_num][2][-1])
    label_23.config(text=data_values[array_num][2][12])
    label_24.config(text=data_values[array_num][2][13])
    label_25.config(text=data_values[array_num][2][14])
    label_26.config(text=data_values[array_num][2][15])
    label_27.config(text=data_values[array_num][2][16])
    
    label_28.config(text=data_values[array_num][3][0])
    label_29.config(text=data_values[array_num][3][-1])
    label_30.config(text=data_values[array_num][3][1])
    label_31.config(text=data_values[array_num][3][2])
    label_32.config(text=data_values[array_num][3][3])
    label_33.config(text=data_values[array_num][3][4])
    label_34.config(text=data_values[array_num][3][8])
    label_35.config(text=data_values[array_num][3][5])
    label_36.config(text=data_values[array_num][3][9])

def design():
    """Gets the values in the Entry boxes and passes them to a function"""
    try:
        maxMoment_L = np.abs(float(Mu_L.get()) * np.power(10, 6))
        maxMoment_C = np.abs(float(Mu_C.get()) * np.power(10, 6))
        maxMoment_R = np.abs(float(Mu_R.get()) * np.power(10, 6))
        baseBeam = float(base.get())
        heightBeam = float(height.get())
        ccBeamT = float(ccTension.get())
        ccBeamS = float(ccSide.get())
        f_cBeam = float(f_c.get())
        f_yBeam = float(f_y.get())
        depthBeamT = heightBeam - ccBeamT
        diameter_L = diameterBoxT_L.get()
        diameter_C = diameterBoxT_C.get()
        diameter_R = diameterBoxT_R.get()
        lengthBeam = float(length.get())
        max_shear = np.abs(float(shear_max.get()))
        mid_shear = np.abs(float(shear_mid.get()))
        lambda_concrete = float(concrete_lambda.get())
        f_ytBeam = float(f_yt.get())
        stirrup_d = diameterBoxS.get()


    except ValueError:
        messagebox.showerror("Check Values", "Check Inputed Values!")
        return


    moments = [maxMoment_L, maxMoment_C, maxMoment_R]
    diameter = [diameter_L, diameter_C, diameter_R]
    beam_array = []

    if data_values:
        data_values.pop(-1)
    try:
        for i in range(3):
            # Beam parameters solution
            d_beam_T, steel_ratio_max, steel_ratio_min1, steel_ratio_min2, steel_ratio, area_steel_T, no_of_bars_T, spacing_T, capacity, beta, steel_ratio_min = solve_beam(moments[i], baseBeam, ccBeamT, ccBeamS, f_cBeam, f_yBeam, depthBeamT, diameter[i])

            # Append data values
            append_val = [moments[i], baseBeam, heightBeam, ccBeamT, ccBeamS, f_cBeam, f_yBeam,
                depthBeamT, beta, d_beam_T, steel_ratio_max, steel_ratio_min, steel_ratio,
                area_steel_T, no_of_bars_T, spacing_T, capacity/np.power(10, 6), steel_ratio_min1, steel_ratio_min2]

            ######################################################################################################################################
            # Index Value of append_val list in the first 3 of the second index i.e. [:][1,2,3]:                                                 #
            # moment = 0, baseBeam = 1, heightBeam = 2, ccBeamT = 3, ccBeamS = 4, f_cBeam = 5, f_yBeam = 6, depthBeamT = 7, beta = 8             #
            # d_beam_T = 9, steel_ratio_max = 10, steel_ratio_min = 11, steel_ratio = 12, area_steel_T = 13, no_of_bars_T = 14, spacing_T = 15   #
            # capacity/np.power(10, 6) = 16, steel_ratio_min1 = 17, steel_ratio_min2 = 18                                                        #
            ######################################################################################################################################

            append_val = [round(num, 12) for num in append_val]
            beam_array.append(append_val)
            data_values.append(beam_array)

        # Stirrup parameters solution plus append values
        V_c, centerline, V_u, V_s, S_max, A_v_min, S_value, dist_of_stirrup, A_v, d_stirrup = solve_stirrups(max_shear, mid_shear, lengthBeam, baseBeam, depthBeamT, f_cBeam, f_ytBeam, lambda_concrete, stirrup_d)
        append_val = [V_c, V_u, V_s, S_max, A_v_min, S_value, lengthBeam, dist_of_stirrup, A_v, d_stirrup]
        append_val = [round(num, 12) for num in append_val]
        append_val.append(centerline)
        beam_array.append(append_val)

    except TypeError:
        return

    # Update Values in wrapper 2
    update_solved(-1)

    # Draw the recent inputed values for the cross section
    draw_beam(-1)

    # Draw the recent inputed values for the stirrups
    draw_stirrups(-1)


mainWindow = Tk()
mainWindow.title("SRB Beam Designer (Metric Units)")
mainWindow.geometry("1180x760+0+0")
mainWindow.resizable(width=0, height=0)

# Wrappers
wrapper1 = LabelFrame(mainWindow, text="Beam Properties")
wrapper2 = LabelFrame(mainWindow, text="Design Specifications")
wrapper3 = LabelFrame(mainWindow, text="Cross Section Drawing")
wrapper4 = LabelFrame(mainWindow, text="Stirrup Design Drawing")

wrapper1.place(x=10, y=10, width=630, height=400)
wrapper2.place(x=10, y=420, width=1160, height=330)
wrapper3.place(x=650, y=10, width=520, height=200)
wrapper4.place(x=650, y=220, width=520, height=190)


# Input Area
Mu_L = StringVar()
Mu_C = StringVar()
Mu_R = StringVar()
base = StringVar()
height = StringVar()
ccTension= StringVar()
ccSide = StringVar()
f_c = StringVar()
f_y = StringVar()
length = StringVar()
shear_max = StringVar()
shear_mid = StringVar()
concrete_lambda = StringVar()
f_yt = StringVar()


# Global Variables
diameters = [9.5, 12.7, 15.9, 19.1, 22.2, 25.4, 28.7, 32.3, 35.8, 43, 57.3]
data_values = [] # [ Beam No. ][ Placement (Left, Center, Right, Stirrups) ][ Specifications ]
table_values = []


# Wrapper 1 Widgets
Label(wrapper1, text="Left Mu:").place(x=10, y=10)
Label(wrapper1, text="Center Mu:").place(x=10, y=40)
Label(wrapper1, text="Right Mu:").place(x=10, y=70)
Label(wrapper1, text="Base:").place(x=10, y=100)
Label(wrapper1, text="Height:").place(x=10, y=130)
Label(wrapper1, text="Bottom Bar cc:").place(x=10, y=160)
Label(wrapper1, text="Side cc:").place(x=10, y=190)
Label(wrapper1, text="f'c:").place(x=10, y=220)
Label(wrapper1, text="fy:").place(x=10, y=250)
Label(wrapper1, text="Left Diameter:").place(x=10, y=280)
Label(wrapper1, text="Center Diameter:").place(x=10, y=310)
Label(wrapper1, text="Right Diameter:").place(x=10, y=340)

Label(wrapper1, text="Length:").place(x=320, y=10)
Label(wrapper1, text="Max Shear:").place(x=320, y=40)
Label(wrapper1, text="Middle Shear:").place(x=320, y=70)
Label(wrapper1, text="Concrete Factor (λ):").place(x=320, y=100)
Label(wrapper1, text="fyt:").place(x=320, y=130)
Label(wrapper1, text="Stirrup Diameter:").place(x=320, y=160)


showBtn = Button(wrapper1, text="Calculate", command=design, width=15)
showBtn.place(x=340, y=200, width=100)


### Entry Boxes
Entry(wrapper1, textvariable=Mu_L).place(x=130, y=10)
Entry(wrapper1, textvariable=Mu_C).place(x=130, y=40)
Entry(wrapper1, textvariable=Mu_R).place(x=130, y=70)
Entry(wrapper1, textvariable=base).place(x=130, y=100)
Entry(wrapper1, textvariable=height).place(x=130, y=130)
Entry(wrapper1, textvariable=ccTension).place(x=130, y=160)
Entry(wrapper1, textvariable=ccSide).place(x=130, y=190)
Entry(wrapper1, textvariable=f_c).place(x=130, y=220)
Entry(wrapper1, textvariable=f_y).place(x=130, y=250)

Entry(wrapper1, textvariable=length).place(x=450, y=10)
Entry(wrapper1, textvariable=shear_max).place(x=450, y=40)
Entry(wrapper1, textvariable=shear_mid).place(x=450, y=70)
Entry(wrapper1, textvariable=concrete_lambda).place(x=450, y=100)
Entry(wrapper1, textvariable=f_yt).place(x=450, y=130)


diameterBoxT_L = ttk.Combobox(wrapper1, values=["Specify Diameter", "Auto", 9.5, 12.7, 15.9, 19.1, 22.2, 25.4, 28.7, 32.3, 35.8, 43, 57.3])
diameterBoxT_L.place(x=130, y=280, width=125)
diameterBoxT_L.current(0)

diameterBoxT_C = ttk.Combobox(wrapper1, values=["Specify Diameter", "Auto", 9.5, 12.7, 15.9, 19.1, 22.2, 25.4, 28.7, 32.3, 35.8, 43, 57.3])
diameterBoxT_C.place(x=130, y=310, width=125)
diameterBoxT_C.current(0)

diameterBoxT_R = ttk.Combobox(wrapper1, values=["Specify Diameter", "Auto", 9.5, 12.7, 15.9, 19.1, 22.2, 25.4, 28.7, 32.3, 35.8, 43, 57.3])
diameterBoxT_R.place(x=130, y=340, width=125)
diameterBoxT_R.current(0)

diameterBoxS = ttk.Combobox(wrapper1, values=["Specify Diameter", "Auto", 9.5, 12.7, 15.9, 19.1, 22.2, 25.4, 28.7, 32.3, 35.8, 43, 57.3])
diameterBoxS.place(x=450, y=160, width=125)
diameterBoxS.current(0)


Label(wrapper1, text="kN-m").place(x=260, y=10)
Label(wrapper1, text="kN-m").place(x=260, y=40)
Label(wrapper1, text="kN-m").place(x=260, y=70)
Label(wrapper1, text="mm").place(x=260, y=100)
Label(wrapper1, text="mm").place(x=260, y=130)
Label(wrapper1, text="mm").place(x=260, y=160)
Label(wrapper1, text="mm").place(x=260, y=190)
Label(wrapper1, text="MPa").place(x=260, y=220)
Label(wrapper1, text="MPa").place(x=260, y=250)
Label(wrapper1, text="mm").place(x=260, y=280)
Label(wrapper1, text="mm").place(x=260, y=310)
Label(wrapper1, text="mm").place(x=260, y=340)


Label(wrapper1, text="mm").place(x=580, y=10)
Label(wrapper1, text="kN").place(x=580, y=40)
Label(wrapper1, text="kN").place(x=580, y=70)
Label(wrapper1, text="MPa").place(x=580, y=130)
Label(wrapper1, text="mm").place(x=580, y=160)

# 2nd Window
Label(wrapper2, text="Diameter (mm):").place(x=10, y=10)
Label(wrapper2, text="Max Steel Ratio:").place(x=10, y=40)
Label(wrapper2, text="Min Steel Ratio - 1:").place(x=10, y=70)
Label(wrapper2, text="Min Steel Ratio - 2:").place(x=10, y=100)
Label(wrapper2, text="Steel Ratio:").place(x=10, y=130)
Label(wrapper2, text="Steel Area (sq. mm):").place(x=10, y=160)
Label(wrapper2, text="No. of Bars:").place(x=10, y=190)
Label(wrapper2, text="Bars Spacing (mm):").place(x=10, y=220)
Label(wrapper2, text="Beam Capacity (kN-m):").place(x=10, y=250)


label_1 = Label(wrapper2, text="")
label_2 = Label(wrapper2, text="")
label_3 = Label(wrapper2, text="")
label_4 = Label(wrapper2, text="")
label_5 = Label(wrapper2, text="")
label_6 = Label(wrapper2, text="")
label_7 = Label(wrapper2, text="")
label_8 = Label(wrapper2, text="")
label_9 = Label(wrapper2, text="")

label_1.place(x=150, y=10)
label_2.place(x=150, y=40)
label_3.place(x=150, y=70)
label_4.place(x=150, y=100)
label_5.place(x=150, y=130)
label_6.place(x=150, y=160)
label_7.place(x=150, y=190)
label_8.place(x=150, y=220)
label_9.place(x=150, y=250)


Label(wrapper2, text="Diameter (mm):").place(x=300, y=10)
Label(wrapper2, text="Max Steel Ratio:").place(x=300, y=40)
Label(wrapper2, text="Min Steel Ratio - 1:").place(x=300, y=70)
Label(wrapper2, text="Min Steel Ratio - 2:").place(x=300, y=100)
Label(wrapper2, text="Steel Ratio:").place(x=300, y=130)
Label(wrapper2, text="Steel Area (sq. mm):").place(x=300, y=160)
Label(wrapper2, text="No. of Bars:").place(x=300, y=190)
Label(wrapper2, text="Bars Spacing (mm):").place(x=300, y=220)
Label(wrapper2, text="Beam Capacity (kN-m):").place(x=300, y=250)


label_10 = Label(wrapper2, text="")
label_11 = Label(wrapper2, text="")
label_12 = Label(wrapper2, text="")
label_13 = Label(wrapper2, text="")
label_14 = Label(wrapper2, text="")
label_15 = Label(wrapper2, text="")
label_16 = Label(wrapper2, text="")
label_17 = Label(wrapper2, text="")
label_18 = Label(wrapper2, text="")

label_10.place(x=440, y=10)
label_11.place(x=440, y=40)
label_12.place(x=440, y=70)
label_13.place(x=440, y=100)
label_14.place(x=440, y=130)
label_15.place(x=440, y=160)
label_16.place(x=440, y=190)
label_17.place(x=440, y=220)
label_18.place(x=440, y=250)


Label(wrapper2, text="Diameter (mm):").place(x=590, y=10)
Label(wrapper2, text="Max Steel Ratio:").place(x=590, y=40)
Label(wrapper2, text="Min Steel Ratio - 1:").place(x=590, y=70)
Label(wrapper2, text="Min Steel Ratio - 2:").place(x=590, y=100)
Label(wrapper2, text="Steel Ratio:").place(x=590, y=130)
Label(wrapper2, text="Steel Area (sq. mm):").place(x=590, y=160)
Label(wrapper2, text="No. of Bars:").place(x=590, y=190)
Label(wrapper2, text="Bars Spacing (mm):").place(x=590, y=220)
Label(wrapper2, text="Beam Capacity (kN-m):").place(x=590, y=250)


label_19 = Label(wrapper2, text="")
label_20 = Label(wrapper2, text="")
label_21 = Label(wrapper2, text="")
label_22 = Label(wrapper2, text="")
label_23 = Label(wrapper2, text="")
label_24 = Label(wrapper2, text="")
label_25 = Label(wrapper2, text="")
label_26 = Label(wrapper2, text="")
label_27 = Label(wrapper2, text="")

label_19.place(x=730, y=10)
label_20.place(x=730, y=40)
label_21.place(x=730, y=70)
label_22.place(x=730, y=100)
label_23.place(x=730, y=130)
label_24.place(x=730, y=160)
label_25.place(x=730, y=190)
label_26.place(x=730, y=220)
label_27.place(x=730, y=250)


Label(wrapper2, text="V_c (N):").place(x=880, y=10)
Label(wrapper2, text="Reach centerline:").place(x=880, y=40)
Label(wrapper2, text="V_u (N):").place(x=880, y=70)
Label(wrapper2, text="V_s (N):").place(x=880, y=100)
Label(wrapper2, text="Max Spacing (mm):").place(x=880, y=130)
Label(wrapper2, text="Minimum A_v (sq. mm):").place(x=880, y=160)
Label(wrapper2, text="Used A_v:").place(x=880, y=190)
Label(wrapper2, text="Spacing Used (mm):").place(x=880, y=220)
Label(wrapper2, text="Diameter (mm):").place(x=880, y=250)

label_28 = Label(wrapper2, text="")
label_29 = Label(wrapper2, text="")
label_30 = Label(wrapper2, text="")
label_31 = Label(wrapper2, text="")
label_32 = Label(wrapper2, text="")
label_33 = Label(wrapper2, text="")
label_34 = Label(wrapper2, text="")
label_35 = Label(wrapper2, text="")
label_36 = Label(wrapper2, text="")

label_28.place(x=1020, y=10)
label_29.place(x=1020, y=40)
label_30.place(x=1020, y=70)
label_31.place(x=1020, y=100)
label_32.place(x=1020, y=130)
label_33.place(x=1020, y=160)
label_34.place(x=1020, y=190)
label_35.place(x=1020, y=220)
label_36.place(x=1020, y=250)

Label(wrapper2, text="---------Left Moment---------").place(x=80, y=280)
Label(wrapper2, text="---------Center Moment---------").place(x=360, y=280)
Label(wrapper2, text="---------Right Moment---------").place(x=650, y=280)
Label(wrapper2, text="---------For Stirrups---------").place(x=940, y=280)


# Treeview for Data

mainWindow.mainloop()

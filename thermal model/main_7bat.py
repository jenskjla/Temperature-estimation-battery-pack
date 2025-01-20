"""
==================================================================================================
Author: Yangxiao Xiang @ CityU , yaxiang@cityu.edu.hk
Note:
==================================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt 
import EKF_Estimation
import pandas as pd
import datetime

"""
Main Script - Handles the jacobians, Motion,Measurement Models,reading the data file and helper functions. 
"""

srcData = "dataset_pack7_NewParam"
Ta = 35
nRe1 = 5 * 1e-3 * 0.85
nRe2 = 5 * 1e-3 * 1
nRe3 = 5 * 1e-3 * 1.15
nRe4 = 5 * 1e-3 * 0.9
nRe5 = 5 * 1e-3 * 0.85
nRe6 = 5 * 1e-3 * 1.1
nRe7 = 5 * 1e-3 * 1.0
nRe = [nRe1, nRe2, nRe3, nRe4, nRe5, nRe6, nRe7 ]
df = pd.read_csv("./thermal model/Simulation_data/Q_values.csv")
starttime = 1
endtime = 25001
dataframeQ = df.iloc[starttime:endtime, :].reset_index()
Qv = dataframeQ[["Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"]].values
nCc = 6.7 * 1e1
nCs = 3.115 * 1e0
nRc = 1.83 * 1e0
nRu = 4.03 * 1e0
nRcc = 4 * 1e-1


def motion_model(control_input,x_k_1, delta_t, k):

    Cc = x_k_1[14][0]
    Cs = x_k_1[15][0]
    Rc = x_k_1[16][0]
    Ru = x_k_1[17][0]
    Rcc = x_k_1[18][0]

    state_matrix = np.array([[-1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [1/(Rc*Cs), -1*(1/Cs)*(1/Ru+1/Rc+1/Rcc), 0, 1/(Rcc*Cs), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 1/(Rcc*Cs), 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc+2/Rcc), 0, 1/(Rcc*Cs), 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 1/(Rcc*Cs), 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc+2/Rcc), 0, 1/(Rcc*Cs), 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 1/(Rcc*Cs), 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc+2/Rcc), 0, 1/(Rcc*Cs), 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 1/(Rcc*Cs), 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc+2/Rcc), 0, 1/(Rcc*Cs), 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(Rcc*Cs), 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc+2/Rcc), 0, 1/(Rcc*Cs)],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc)],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(Rcc*Cs), 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc+1/Rcc)],
                             ], dtype='float')


    B_matrix = np.array([[1/Cc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 1/(Cs*Ru), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #Change dimension of the matrix 14x14
                         [0, 0, 1/Cc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 1/(Cs*Ru), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1/Cc, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 1/(Cs*Ru), 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1/Cc, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1/(Cs*Ru), 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1/Cc, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(Cs*Ru), 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Cc, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(Cs*Ru), 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Cc, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(Cs*Ru)]])

    prev_state = x_k_1[:14, :].reshape([-1, 1])
    new_state = prev_state + (state_matrix.dot(prev_state) + B_matrix.dot(control_input)) * delta_t # Next State matrix
    new_state = np.concatenate((new_state, x_k_1[14:, :]), axis=0)
    return new_state


def jacobian_measurement_model(x_k_1):


    H = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]])

    M = np.identity(2)

    y_out = H.dot(x_k_1)

    return H ,M, y_out

def jacobian_motion_model(x_k_1, ISquare, delta_t, Qv, k):

    #  Motion Model Jacobian w.r.t Posterior
    Tc1 = x_k_1[0][0]
    Ts1 = x_k_1[1][0]
    Tc2 = x_k_1[2][0]
    Ts2 = x_k_1[3][0]
    Tc3 = x_k_1[4][0]
    Ts3 = x_k_1[5][0]
    Tc4 = x_k_1[6][0]
    Ts4 = x_k_1[7][0]
    Tc5 = x_k_1[8][0]
    Ts5 = x_k_1[9][0]
    Tc6 = x_k_1[10][0]
    Ts6 = x_k_1[11][0]
    Tc7 = x_k_1[12][0]
    Ts7 = x_k_1[13][0]
    Cc = x_k_1[14][0]
    Cs = x_k_1[15][0]
    Rc = x_k_1[16][0]
    Ru = x_k_1[17][0]
    Rcc = x_k_1[18][0]
    ISquare = ISquare[0]

    state_matrix1 = np.array([[1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [delta_t / (Rc * Cs), 1-delta_t * (1 / Cs) * (1 / Ru + 1 / Rc + 1 / Rcc), 0, delta_t / (Rcc * Cs), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, delta_t / (Rcc * Cs), delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Ru + 1 / Rc + 2 / Rcc), 0, delta_t / (Rcc * Cs), 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, delta_t / (Rcc * Cs), delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Ru + 1 / Rc + 2 / Rcc), 0, delta_t / (Rcc * Cs), 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, delta_t / (Rcc * Cs), delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Ru + 1 / Rc + 2 / Rcc), 0, delta_t / (Rcc * Cs), 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, delta_t / (Rcc * Cs), delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Ru + 1 / Rc + 2 / Rcc), 0, delta_t / (Rcc * Cs), 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, delta_t / (Rcc * Cs), delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Ru + 1 / Rc + 2 / Rcc), 0, delta_t / (Rcc * Cs)],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc)],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, delta_t / (Rcc * Cs), delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Ru + 1 / Rc + 1 / Rcc)],
                             ], dtype='float')


    Q1 = Qv[k, 0] #Change this to proper equation of Q
    Q2 = Qv[k, 1]
    Q3 = Qv[k, 2]
    Q4 = Qv[k, 3]
    Q5 = Qv[k, 4]
    Q6 = Qv[k, 5]
    Q7 = Qv[k, 6]

    A11 = Tc1 / (Rc * Cc * Cc) - Ts1 / (Rc * Cc * Cc) - Q1 / (Cc * Cc)
    A12 = 0
    A13 = Tc1 / (Rc * Rc * Cc) - Ts1 / (Rc * Rc * Cc)
    A14 = 0
    A15 = 0
    A21 = 0
    A22 = -Tc1 / (Rc * Cs * Cs) + Ts1 * (1 / Ru + 1 / Rc + 1 / Rcc) / (Cs * Cs) - Ta / (Cs * Cs * Ru) - Ts2 / (Rcc * Cs * Cs)
    A23 = -Tc1 / (Rc * Rc * Cs) + Ts1 / (Rc * Rc * Cs)
    A24 = Ts1 / (Cs * Ru * Ru) - Ta / (Cs * Ru * Ru)
    A25 = (Ts1 - Ts2) / (Cs * Rcc * Rcc)
    A31 = Tc2 / (Rc * Cc * Cc) - Ts2 / (Rc * Cc * Cc) - Q2 / (Cc * Cc)
    A32 = 0
    A33 = Tc2 / (Rc * Rc * Cc) - Ts2 / (Rc * Rc * Cc)
    A34 = 0
    A35 = 0
    A41 = 0
    A42 = -Ts1 / (Rcc * Cs * Cs) - Tc2 / (Rc * Cs * Cs) + Ts2 * (1 / Ru + 1 / Rc + 2 / Rcc) / (Cs * Cs) - Ts3 / (Rcc * Cs * Cs) - Ta / (Cs * Cs * Ru)
    A43 = (Ts2 - Tc2) / (Rc * Rc * Cs)
    A44 = (Ts2 - Ta) / (Ru * Ru * Cs)
    A45 = -Ts1 / (Rcc * Rcc * Cs) + 2 * Ts2 / (Rcc * Rcc * Cs) - Ts3 / (Rcc * Rcc * Cs)
    A51 = Tc3 / (Rc * Cc * Cc) - Ts3 / (Rc * Cc * Cc) - Q3 / (Cc * Cc)
    A52 = 0
    A53 = Tc3 / (Rc * Rc * Cc) - Ts3 / (Rc * Rc * Cc)
    A54 = 0
    A55 = 0

    A61 = 0
    A62 = -Ts2 / (Rcc * Cs * Cs) - Tc3 / (Rc * Cs * Cs) + Ts3 * (1 / Ru + 1 / Rc + 2 / Rcc) / (Cs * Cs) - Ts4 / (Rcc * Cs * Cs) - Ta / (Cs * Cs * Ru)
    A63 = (Ts3 - Tc3) / (Rc * Rc * Cs)
    A64 = (Ts3 - Ta) / (Ru * Ru * Cs)
    A65 = -Ts2 / (Rcc * Rcc * Cs) + 2 * Ts3 / (Rcc * Rcc * Cs) - Ts4 / (Rcc * Rcc * Cs)

    A71 = Tc4 / (Rc * Cc * Cc) - Ts4 / (Rc * Cc * Cc) - Q4 / (Cc * Cc)
    A72 = 0
    A73 = Tc4 / (Rc * Rc * Cc) - Ts4 / (Rc * Rc * Cc)
    A74 = 0
    A75 = 0

    A81 = 0
    A82 = -Ts3 / (Rcc * Cs * Cs) - Tc4 / (Rc * Cs * Cs) + Ts4 * (1 / Ru + 1 / Rc + 2 / Rcc) / (Cs * Cs) - Ts5 / (Rcc * Cs * Cs) - Ta / (Cs * Cs * Ru)
    A83 = (Ts4 - Tc4) / (Rc * Rc * Cs)
    A84 = (Ts4 - Ta) / (Ru * Ru * Cs)
    A85 = -Ts3 / (Rcc * Rcc * Cs) + 2 * Ts4 / (Rcc * Rcc * Cs) - Ts5 / (Rcc * Rcc * Cs)

    A91 = Tc5 / (Rc * Cc * Cc) - Ts5 / (Rc * Cc * Cc) - Q5 / (Cc * Cc)
    A92 = 0
    A93 = Tc5 / (Rc * Rc * Cc) - Ts5 / (Rc * Rc * Cc)
    A94 = 0
    A95 = 0

    A101 = 0
    A102 = -Ts4 / (Rcc * Cs * Cs) - Tc5 / (Rc * Cs * Cs) + Ts5 * (1 / Ru + 1 / Rc + 2 / Rcc) / (Cs * Cs) - Ts6 / (Rcc * Cs * Cs) - Ta / (Cs * Cs * Ru)
    A103 = (Ts5 - Tc5) / (Rc * Rc * Cs)
    A104 = (Ts5 - Ta) / (Ru * Ru * Cs)
    A105 = -Ts4 / (Rcc * Rcc * Cs) + 2 * Ts5 / (Rcc * Rcc * Cs) - Ts6 / (Rcc * Rcc * Cs)

    A111 = Tc6 / (Rc * Cc * Cc) - Ts6 / (Rc * Cc * Cc) - Q6 / (Cc * Cc)
    A112 = 0
    A113 = Tc6 / (Rc * Rc * Cc) - Ts6 / (Rc * Rc * Cc)
    A114 = 0
    A115 = 0

    A121 = 0
    A122 = -Ts5 / (Rcc * Cs * Cs) - Tc6 / (Rc * Cs * Cs) + Ts6 * (1 / Ru + 1 / Rc + 2 / Rcc) / (Cs * Cs) - Ts7 / (Rcc * Cs * Cs) - Ta / (Cs * Cs * Ru)
    A123 = (Ts6 - Tc6) / (Rc * Rc * Cs)
    A124 = (Ts6 - Ta) / (Ru * Ru * Cs)
    A125 = -Ts5 / (Rcc * Rcc * Cs) + 2 * Ts6 / (Rcc * Rcc * Cs) - Ts7 / (Rcc * Rcc * Cs)

    A131 = Tc7 / (Rc * Cc * Cc) - Ts7 / (Rc * Cc * Cc) - Q7 / (Cc * Cc)
    A132 = 0
    A133 = Tc7 / (Rc * Rc * Cc) - Ts7 / (Rc * Rc * Cc)
    A134 = 0
    A135 = 0

    A141 = 0
    A142 = -Tc7 / (Rc * Cs * Cs) + Ts7 * (1 / Ru + 1 / Rc + 1 / Rcc) / (Cs * Cs) - Ta / (Cs * Cs * Ru) - Ts6 / (Rcc * Cs * Cs)
    A143 = (Ts7 - Tc7) / (Rc * Rc * Cs)
    A144 = (Ts7 - Ta) / (Ru * Ru * Cs)
    A145 = (Ts7 - Ts6) / (Cs * Rcc * Rcc)

    UpperRight = np.array([[A11 * delta_t, A12 * delta_t, A13 * delta_t, A14 * delta_t, A15 * delta_t],
                         [A21 * delta_t, A22 * delta_t, A23 * delta_t, A24 * delta_t, A25 * delta_t],
                         [A31 * delta_t, A32 * delta_t, A33 * delta_t, A34 * delta_t, A35 * delta_t],
                         [A41 * delta_t, A42 * delta_t, A43 * delta_t, A44 * delta_t, A45 * delta_t],
                         [A51 * delta_t, A52 * delta_t, A53 * delta_t, A54 * delta_t, A55 * delta_t],
                         [A61 * delta_t, A62 * delta_t, A63 * delta_t, A64 * delta_t, A65 * delta_t],
                           [A71 * delta_t, A72 * delta_t, A73 * delta_t, A74 * delta_t, A75 * delta_t],
                           [A81 * delta_t, A82 * delta_t, A83 * delta_t, A84 * delta_t, A85 * delta_t],
                           [A91 * delta_t, A92 * delta_t, A93 * delta_t, A94 * delta_t, A95 * delta_t],
                           [A101 * delta_t, A102 * delta_t, A103 * delta_t, A104 * delta_t, A105 * delta_t],
                           [A111 * delta_t, A112 * delta_t, A113 * delta_t, A114 * delta_t, A115 * delta_t],
                           [A121 * delta_t, A122 * delta_t, A123 * delta_t, A124 * delta_t, A125 * delta_t],
                           [A131 * delta_t, A132 * delta_t, A133 * delta_t, A134 * delta_t, A135 * delta_t],
                           [A141 * delta_t, A142 * delta_t, A143 * delta_t, A144 * delta_t, A145 * delta_t],])
    LowerLeft = np.zeros((5, 14))
    LowerRight = np.identity(5)

    Upper = np.concatenate((state_matrix1, UpperRight), axis=1)
    Lower = np.concatenate((LowerLeft, LowerRight), axis=1)
    F = np.concatenate((Upper,Lower), axis=0)

    # Motion Model Jacobian w.r.t Posterior

    L = np.array([[1]])
    
    return F, L 


def data_initialization(QState, QParam, RState, PState, PParam):

    Cc_init = nCc * 1.1
    Cs_init = nCs * 0.9
    Rc_init = nRc * 2
    Ru_init = nRu * 0.7
    Rcc_init = nRcc * 1.5
    

    df = pd.read_csv("./thermal model/Simulation_data/" + srcData + ".csv")
    starttime = 1
    endtime = 25001
    dataframe = df.iloc[starttime:endtime, :].reset_index()
    
    data_t = dataframe["t"].values
    data_ib = dataframe["ib"].values
    States_True = dataframe[["Tc1", "Ts1", "Tc2", "Ts2", "Tc3", "Ts3", "Tc4", "Ts4", "Tc5", "Ts5", "Tc6", "Ts6", "Tc7", "Ts7"]]

    state_init = np.array([[Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Cc_init, Cs_init, Rc_init, Ru_init, Rcc_init]])
    ibSquare = data_ib.reshape(1, -1) * data_ib.reshape(1, -1)
    Control = []
    for i in range(len(data_t)):
        tmp = []
        for n in range(len(nRe)):
            Re = nRe[n]
            tmp.append([ibSquare[0, i]*Re])
            tmp.append([Ta])
        Control.append(np.array(tmp))

    Q = np.diag([QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QParam, QParam, QParam, QParam, QParam])  # Intializing Process Noise Co-Variance Matrix
    R = np.diag([RState, RState, RState, RState]) # Intializing Measurement Noise Co-Variance Matrix

    x_est = np.zeros([len(data_t), state_init.shape[1]])  # Initializing State Estimate Matrix
    x_est[0, :] = state_init
    P_est = np.zeros([len(data_t), state_init.shape[1], state_init.shape[1]])  # Initializing Co-Variance Matrix
    P_est[0] = np.diag([PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PParam, PParam, PParam, PParam, PParam])

    return Control, data_t, Q, R, x_est, P_est, States_True



if __name__ == '__main__':

    QState = 1e-8
    QParam = 1e-10 #Higher -> Trusts more in model
    RState = 1e-8
    PState = 1e-6
    PParam = 5e-6

    print("START")
    Breakflag = False
    Ap = ""
    ekf = EKF_Estimation.EKF() # EKF Object

    controls, t, Q, R, x_est, P_est, States_True= data_initialization(QState, QParam, RState, PState, PParam)

    x_k_1 = x_est[0,:].reshape(-1,1)
    P_k_1 = P_est[0]

    # Main Time loop for Prediction and Measurement

    for k in range(1,len(t)):

        delta_t = t[k] - t[k-1] # Time Period

        # Prediction Step
        x_k_1,P_k_1= ekf.prediction_step(delta_t,controls[k],x_k_1,P_k_1,Q,jacobian_motion_model,motion_model, Qv, k)

        x_k_1,P_k_1, Breakflag = ekf.measurement_update(States_True[["Ts1", "Ts2", "Ts5", "Ts7"]].values[k, :],P_k_1,x_k_1,Q,R,jacobian_measurement_model)
        if Breakflag:
            break

        # Estimated States

        x_est[k, :] = x_k_1.T
        P_est[k,:,:] = P_k_1

    # Plotting the estimated Values.
    Tc1_pred = x_est[:, 0]
    Ts1_pred = x_est[:, 1]
    Tc2_pred = x_est[:, 2]
    Ts2_pred = x_est[:, 3]
    Tc3_pred = x_est[:, 4]
    Ts3_pred = x_est[:, 5]
    Tc4_pred = x_est[:, 6]
    Ts4_pred = x_est[:, 7]
    Tc5_pred = x_est[:, 8]
    Ts5_pred = x_est[:, 9]
    Tc6_pred = x_est[:, 10]
    Ts6_pred = x_est[:, 11]
    Tc7_pred = x_est[:, 12]
    Ts7_pred = x_est[:, 13]
    Cc_pred = x_est[:, 14][-1]
    Cs_pred = x_est[:, 15][-1]
    Rc_pred = x_est[:, 16][-1]
    Ru_pred = x_est[:, 17][-1]
    Rcc_pred = x_est[:, 18][-1]
    Trace = P_est[-1].trace()

    Tc1_true = States_True["Tc1"].values
    Ts1_true = States_True["Ts1"].values
    Tc2_true = States_True["Tc2"].values
    Ts2_true = States_True["Ts2"].values
    Tc3_true = States_True["Tc3"].values
    Ts3_true = States_True["Ts3"].values
    Tc4_true = States_True["Tc4"].values
    Ts4_true = States_True["Ts4"].values
    Tc5_true = States_True["Tc5"].values
    Ts5_true = States_True["Ts5"].values
    Tc6_true = States_True["Tc6"].values
    Ts6_true = States_True["Ts6"].values
    Tc7_true = States_True["Tc7"].values
    Ts7_true = States_True["Ts7"].values
    Tc1Error = np.mean(np.abs(Tc1_pred - Tc1_true))
    Ts1Error = np.mean(np.abs(Ts1_pred - Ts1_true)) * 1e12
    Tc2Error = np.mean(np.abs(Tc2_pred - Tc2_true))
    Ts2Error = np.mean(np.abs(Ts2_pred - Ts2_true)) * 1e12
    Tc3Error = np.mean(np.abs(Tc3_pred - Tc3_true))
    Ts3Error = np.mean(np.abs(Ts3_pred - Ts3_true))
    Tc4Error = np.mean(np.abs(Tc4_pred - Tc4_true))
    Ts4Error = np.mean(np.abs(Ts4_pred - Ts4_true))
    Tc5Error = np.mean(np.abs(Tc5_pred - Tc5_true))
    Ts5Error = np.mean(np.abs(Ts5_pred - Ts5_true)) * 1e12
    Tc6Error = np.mean(np.abs(Tc6_pred - Tc6_true))
    Ts6Error = np.mean(np.abs(Ts6_pred - Ts6_true))
    Tc7Error = np.mean(np.abs(Tc7_pred - Tc7_true))
    Ts7Error = np.mean(np.abs(Ts7_pred - Ts7_true)) * 1e12
    CcError = Cc_pred / nCc * 100 - 100
    CsError = Cs_pred / nCs * 100 - 100
    RcError = Rc_pred / nRc * 100 - 100
    RuError = Ru_pred / nRu * 100 - 100
    RccError = Rcc_pred / nRcc * 100 - 100
    meanError = np.mean([np.abs(CcError), np.abs(CsError), np.abs(RcError), np.abs(RuError), np.abs(RccError)])

    print("Tc1Error: %.3f, Ts1Error: %.0f, Tc2Error: %.3f, Ts2Error: %.0f, Tc3Error: %.3f, Ts3Error: %.3f, Tc4Error: %.3f, Ts4Error: %.3f, Tc5Error: %.3f, Ts5Error: %.0f, Tc6Error: %.3f, Ts6Error: %.3f, Tc7Error: %.3f, Ts7Error: %.0f, Cc: %.1f, Cs: %.1f, Rc: %.1f, Ru: %.1f, Rcc: %.1f, MeanError: %.3f, Trace: %.3f"
          % (Tc1Error, Ts1Error, Tc2Error, Ts2Error, Tc3Error, Ts3Error, Tc4Error, Ts4Error, Tc5Error, Ts5Error, Tc6Error, Ts6Error, Tc7Error, Ts7Error, CcError, CsError, RcError, RuError, RccError, meanError, Trace))
    Ap += 'Tc1Error: %.3f, Ts1Error: %.0f, Tc2Error: %.3f, Ts2Error: %.0f, Tc3Error: %.3f, Ts3Error: %.3f, Tc4Error: %.3f, Ts4Error: %.3f, Tc5Error: %.3f, Ts5Error: %.3f, Tc6Error: %.3f, Ts6Error: %.0f, Tc7Error: %.3f, Ts7Error: %.0f, Cc: %.1f, Cs: %.1f, Rc: %.1f, Ru: %.1f, Rcc: %.1f, MeanError: %.3f, Trace: %.3f' % (
            Tc1Error, Ts1Error, Tc2Error, Ts2Error, Tc3Error, Ts3Error, Tc4Error, Ts4Error, Tc5Error, Ts5Error, Tc6Error, Ts6Error, Tc7Error, Ts7Error, CcError, CsError, RcError, RuError, RccError, meanError, Trace) + "\n"
    # plt.plot(t, Ts2_pred)  # Tc1
    # plt.plot(t, Ts2_true)
    # plt.plot(t, Tc2_pred)  # Tc1
    # plt.plot(t, Tc2_true)
    # plt.show()

    for generation in range(1, 10):
        import copy
        x_est_save = copy.deepcopy(x_est)
        P_est[0, :, :] = P_est[-1, :, :]
        x_k_1 = np.concatenate((np.array([Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta]), x_est[-1, 14:]))
        x_k_1 = x_k_1.reshape(-1, 1)
        P_k_1 = P_est[0]
        for k in range(1,len(t)):

            delta_t = t[k] - t[k-1] # Time Period

            # Prediction Step
            x_k_1,P_k_1 = ekf.prediction_step(delta_t,controls[k],x_k_1,P_k_1,Q,jacobian_motion_model,motion_model, Qv, k)

            x_k_1,P_k_1, Breakflag = ekf.measurement_update(States_True[["Ts1", "Ts2", "Ts5", "Ts7"]].values[k, :],P_k_1,x_k_1,Q,R,jacobian_measurement_model)
            if Breakflag:
                break

            # Estimated States

            x_est[k, :] = x_k_1.T
            P_est[k,:,:] = P_k_1
        if Breakflag:
            break
        # Plotting the estimated Values.
        Tc1_pred = x_est[:, 0]
        Ts1_pred = x_est[:, 1]
        Tc2_pred = x_est[:, 2]
        Ts2_pred = x_est[:, 3]
        Tc3_pred = x_est[:, 4]
        Ts3_pred = x_est[:, 5]
        Tc4_pred = x_est[:, 6]
        Ts4_pred = x_est[:, 7]
        Tc5_pred = x_est[:, 8]
        Ts5_pred = x_est[:, 9]
        Tc6_pred = x_est[:, 10]
        Ts6_pred = x_est[:, 11]
        Tc7_pred = x_est[:, 12]
        Ts7_pred = x_est[:, 13]
        Cc_pred = x_est[:, 14][-1]
        Cs_pred = x_est[:, 15][-1]
        Rc_pred = x_est[:, 16][-1]
        Ru_pred = x_est[:, 17][-1]
        Rcc_pred = x_est[:, 18][-1]
        Trace = P_est[-1].trace()

        Tc1_true = States_True["Tc1"].values
        Ts1_true = States_True["Ts1"].values
        Tc2_true = States_True["Tc2"].values
        Ts2_true = States_True["Ts2"].values
        Tc3_true = States_True["Tc3"].values
        Ts3_true = States_True["Ts3"].values
        Tc4_true = States_True["Tc4"].values
        Ts4_true = States_True["Ts4"].values
        Tc5_true = States_True["Tc5"].values
        Ts5_true = States_True["Ts5"].values
        Tc6_true = States_True["Tc6"].values
        Ts6_true = States_True["Ts6"].values
        Tc7_true = States_True["Tc7"].values
        Ts7_true = States_True["Ts7"].values
        Tc1Error = np.mean(np.abs(Tc1_pred - Tc1_true))
        Ts1Error = np.mean(np.abs(Ts1_pred - Ts1_true)) * 1e12
        Tc2Error = np.mean(np.abs(Tc2_pred - Tc2_true))
        Ts2Error = np.mean(np.abs(Ts2_pred - Ts2_true)) * 1e12
        Tc3Error = np.mean(np.abs(Tc3_pred - Tc3_true))
        Ts3Error = np.mean(np.abs(Ts3_pred - Ts3_true))
        Tc4Error = np.mean(np.abs(Tc4_pred - Tc4_true))
        Ts4Error = np.mean(np.abs(Ts4_pred - Ts4_true))
        Tc5Error = np.mean(np.abs(Tc5_pred - Tc5_true))
        Ts5Error = np.mean(np.abs(Ts5_pred - Ts5_true)) * 1e12
        Tc6Error = np.mean(np.abs(Tc6_pred - Tc6_true))
        Ts6Error = np.mean(np.abs(Ts6_pred - Ts6_true))
        Tc7Error = np.mean(np.abs(Tc7_pred - Tc7_true))
        Ts7Error = np.mean(np.abs(Ts7_pred - Ts7_true)) * 1e12
        CcError = Cc_pred / nCc * 100 - 100
        CsError = Cs_pred / nCs * 100 - 100
        RcError = Rc_pred / nRc * 100 - 100
        RuError = Ru_pred / nRu * 100 - 100
        RccError = Rcc_pred / nRcc * 100 - 100
        meanError = np.mean([np.abs(CcError), np.abs(CsError), np.abs(RcError), np.abs(RuError), np.abs(RccError)])

        print(
            "Iterations: %d, Tc1Error: %.3f, Ts1Error: %.0f, Tc2Error: %.3f, Ts2Error: %.0f, Tc3Error: %.3f, Ts3Error: %.3f, Tc4Error: %.3f, Ts4Error: %.3f, Tc5Error: %.3f, Ts5Error: %.0f, Tc6Error: %.3f, Ts6Error: %.3f, Tc7Error: %.3f, Ts7Error: %.0f, Cc: %.1f, Cs: %.1f, Rc: %.1f, Ru: %.1f, Rcc: %.1f, MeanError: %.3f, Trace: %.3f"
            % (generation, Tc1Error, Ts1Error, Tc2Error, Ts2Error, Tc3Error, Ts3Error, Tc4Error, Ts4Error, Tc5Error, Ts5Error, Tc6Error, Ts6Error, Tc7Error, Ts7Error, CcError, CsError, RcError, RuError, RccError, meanError, Trace))
        Ap += 'Iterations: %d, Tc1Error: %.3f, Ts1Error: %.0f, Tc2Error: %.3f, Ts2Error: %.0f, Tc3Error: %.3f, Ts3Error: %.3f, Tc4Error: %.3f, Ts4Error: %.3f, Tc5Error: %.3f, Ts5Error: %.3f, Tc6Error: %.3f, Ts6Error: %.0f, Tc7Error: %.3f, Ts7Error: %.0f, Cc: %.1f, Cs: %.1f, Rc: %.1f, Ru: %.1f, Rcc: %.1f, MeanError: %.3f, Trace: %.3f' % (
            generation, Tc1Error, Ts1Error, Tc2Error, Ts2Error, Tc3Error, Ts3Error, Tc4Error, Ts4Error, Tc5Error, Ts5Error, Tc6Error, Ts6Error, Tc7Error, Ts7Error, CcError, CsError, RcError, RuError, RccError, meanError, Trace) + "\n"

        # plt.plot(t, Ts2_pred)  # Tc1
        # plt.plot(t, Ts2_true)
        # plt.plot(t, Tc2_pred)  # Tc1
        # plt.plot(t, Tc2_true)
        # plt.show()
    NAME = './thermal model/outptKF/Pcontinue_QState' + "{:.0e}".format(QState)  \
          + "QParam" + "{:.0e}".format(QParam) \
            + "RState" + "{:.0e}".format(RState) \
            + "PState" + "{:.0e}".format(PState) \
            + "RParam" + "{:.0e}".format(PParam) + '.txt'
    current_datetime = datetime.datetime.now()
    with open(NAME, 'w') as f:  # 设置文件对象
        f.write(current_datetime.strftime('%Y-%m-%d %H:%M:%S')+'\n')
        f.write(Ap)








"""
Author - Yangxiao Xiang
Co-author - Jens Kjaersgaard Larrañaga
"""
import numpy as np 
import matplotlib.pyplot as plt 
import EKF_Estimation
import pandas as pd
import os
import sympy as sp

"""
Main Script - Handles the jacobians, Motion,Measurement Models,reading the data file and helper functions. 
"""

srcData = "targetWave_elec"
Ta = 21
nRe1 = 1 * 1e-1 * 0.85
nRe2 = 1 * 1e-1 * 1
nRe3 = 1 * 1e-1 * 1.15
nRe4 = 1 * 1e-1 * 0.9
nRe5 = 1 * 1e-1 * 0.85
nRe6 = 1 * 1e-1 * 1.1
nRe7 = 1 * 1e-1 * 1.
nRe = [nRe1, nRe2, nRe3, nRe4, nRe5, nRe6, nRe7]
project_dir = os.path.dirname(os.path.abspath(__file__))
df_Q = pd.read_csv(os.path.join(project_dir, "..", "electrical_model", "Simulation_data", "Qvalues_elec.csv"))
Qv = df_Q[["Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"]].values
nCc = 6.7 * 1e1
nCs = 3.115 * 1e0
nRc = 1.83 * 1e0
nRu = 4.03 * 1e0

def motion_model(control_input,x_k_1, delta_t, Qv, k):

    Cc = x_k_1[14][0]
    Cs = x_k_1[15][0]
    Rc = x_k_1[16][0]
    Ru = x_k_1[17][0]

    state_matrix = np.array([[-1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [1/(Rc*Cs), -1*(1/Cs)*(1/Ru+1/Rc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc), 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc), 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc), 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc), 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc), 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/(Rc*Cc), 1/(Rc*Cc)],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(Rc*Cs), -1/Cs*(1/Ru+1/Rc)],
                             ], dtype=np.float64)


    B_matrix = np.array([[1/Cc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 1/(Cs*Ru), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
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
    new_state = prev_state+(state_matrix.dot(prev_state) + B_matrix.dot(control_input)) * delta_t  # Next State matrix
    new_state = np.concatenate((new_state, x_k_1[14:, :]), axis=0)
    return new_state


def jacobian_measurement_model(x_k_1):


    H = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]])

    M = np.identity(2)

    y_out = H.dot(x_k_1)

    return H ,M, y_out

def jacobian_motion_model(x_k_1, ISquare, delta_t, Qv, k):

    # Define symbolic variables for differentiation
    Tc, Ts, Q, Ta_meas = sp.symbols('Tc Ts Q Ta')
    Cc, Cs, Rc, Rs = sp.symbols('Cc Cs Rc Rs')

    #Create the 
    UpperRight = []
    # State functions
    state_function_1 = -Tc/(Rc*Cc) + Ts/(Rc*Cc) + Q/Cc
    state_function_2 = Tc/(Rc*Cs) + (-Ts * (1/Rs + 1/Rc))/(Cs) + Ta_meas/(Cs * Rs)

    J1 = sp.Matrix([state_function_1]).jacobian([Cc, Cs, Rc, Rs])
    J2 = sp.Matrix([state_function_2]).jacobian([Cc, Cs, Rc, Rs])

    for i in range(1, 8):
        # Substitute values from x_k_1 for each parameter and state
        subs_dict = {
            Tc: x_k_1[(i-1)*2][0],
            Ts: x_k_1[(i-1)*2+1][0],
            Q: Qv[k, i-1],
            Ta_meas: Ta,
            Cc: x_k_1[14][0],
            Cs: x_k_1[15][0],
            Rc: x_k_1[16][0],
            Rs: x_k_1[17][0]
        }

        # Evaluate rows directly:
        row1 = [float((J1[0, j].subs(subs_dict))) * delta_t for j in range(4)]
        row2 = [float((J2[0, j].subs(subs_dict))) * delta_t for j in range(4)]

        UpperRight.append(row1)
        UpperRight.append(row2)
    UpperRight = np.array(UpperRight) 

    Cc = x_k_1[14][0]
    Cs = x_k_1[15][0]
    Rc = x_k_1[16][0]
    Rs = x_k_1[17][0]


    state_matrix1 = np.array([[1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [delta_t / (Rc * Cs), 1-delta_t * (1 / Cs) * (1 / Rs + 1 / Rc ), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Rs + 1 / Rc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Rs + 1 / Rc), 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Rs + 1 / Rc), 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Rs + 1 / Rc), 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc), 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Rs + 1 / Rc), 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1-delta_t / (Rc * Cc), delta_t / (Rc * Cc)],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, delta_t / (Rc * Cs), 1-delta_t / Cs * (1 / Rs + 1 / Rc )],
                             ], dtype=np.float64)
    LowerLeft = np.zeros((4, 14))
    LowerRight = np.identity(4)

    Upper = np.concatenate((state_matrix1, UpperRight), axis=1)
    Lower = np.concatenate((LowerLeft, LowerRight), axis=1)
    F = np.concatenate((Upper,Lower), axis=0)

    # Motion Model Jacobian w.r.t Posterior

    L = np.array([[1]])
    
    return F, L 


def data_initialization(QState, QParam, RState, PState, PParam):

    
    Cc_init = nCc *1.1
    Cs_init = nCs *0.9
    Rc_init = nRc *2
    Ru_init = nRu *0.7

    
    
    df = pd.read_csv(os.path.join(project_dir, "..", "electrical_model", "Simulation_data", srcData + ".csv"))
    starttime = 1
    endtime = 25001
    dataframe = df.iloc[starttime:endtime, :].reset_index()

    data_t = dataframe["t"].values
    States_True = dataframe[["Tc1", "Ts1", "Tc2", "Ts2", "Tc3", "Ts3", "Tc4", "Ts4", "Tc5", "Ts5", "Tc6", "Ts6", "Tc7", "Ts7"]]

    state_init = np.array([[Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Ta, Cc_init, Cs_init, Rc_init, Ru_init]])
    #ibSquare = data_ib.reshape(1, -1) * data_ib.reshape(1, -1)
    Control = []
    num = 7
    for i in range(len(data_t-1)):
        #Control.append(np.array([[ibSquare[0, i]], [Ta]]))
        tmp = []
        for n in range(num):
            #Qv = ibSquare[0][i] * nRe[n]
            Qvi = Qv[i][n]
            tmp.append([Qvi])
            tmp.append([Ta])
        Control.append(np.array(tmp))

        
    Q = np.diag([QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QState, QParam, QParam, QParam, QParam])  # Intializing Process Noise Co-Variance Matrix
    R = np.diag([RState, RState]) # Intializing Measurement Noise Co-Variance Matrix

    x_est = np.zeros([len(data_t), state_init.shape[1]])  # Initializing State Estimate Matrix
    x_est[0, :] = state_init
    P_est = np.zeros([len(data_t), state_init.shape[1], state_init.shape[1]])  # Initializing Co-Variance Matrix
    # P_est[0] = np.eye(state_init.shape[1])*6e-4
    P_est[0] = np.diag([PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PState, PParam, PParam, PParam, PParam])

    return Control, data_t, Q, R, x_est, P_est, States_True



if __name__ == '__main__':
    QStateList = [1e-8]
    QParamList = [1e-10]
    RState = 1e-8
    PStateList = [1e-6]
    PParamList = [5e-6]
    for QState in QStateList:
        for QParam in QParamList:
            for PState in PStateList:
                for PParam in PParamList:
                    if (QParam > QState) or (PParam < PState):
                        continue
                    print("START")
                    Breakflag = False
                    Ap = ""
                    ekf = EKF_Estimation.EKF() # EKF Object

                    controls, t, Q, R, x_est, P_est, States_True = data_initialization(QState, QParam, RState, PState, PParam)

                    x_k_1 = x_est[0,:].reshape(-1,1)
                    P_k_1 = P_est[0]

                    # Main Time loop for Prediction and Measurement

                    for k in range(1,len(t)):

                        delta_t = t[k] - t[k-1] # Time Period

                        # Prediction Step
                        x_k_1,P_k_1 = ekf.prediction_step(delta_t,controls[k],x_k_1,P_k_1,Q,jacobian_motion_model,motion_model, Qv, k)

                        x_k_1,P_k_1, Breakflag = ekf.measurement_update(States_True[["Ts1","Ts5"]].values[k, :],P_k_1,x_k_1,Q,R,jacobian_measurement_model)
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
                    CcError = ((Cc_pred-nCc) / nCc)*100
                    CsError = ((Cs_pred-nCs) / nCs)*100
                    RcError = ((Rc_pred-nRc) / nRc)*100
                    RuError = ((Ru_pred-nRu) / nRu) * 100
                    meanError = np.mean([np.abs(CcError), np.abs(CsError), np.abs(RcError), np.abs(RuError)])

                    print("Tc1Error: %.3f, Ts1Error: %.0f, Tc2Error: %.3f, Ts2Error: %.0f, Tc3Error: %.3f, Ts3Error: %.3f, Tc4Error: %.3f, Ts4Error: %.3f, Tc5Error: %.3f, Ts5Error: %.0f, Tc6Error: %.3f, Ts6Error: %.3f, Tc7Error: %.3f, Ts7Error: %.0f, CcError: %.1f%%, CsError: %.1f%%, RcError: %.1f%%, RuError: %.1f%%, MeanError: %.3f, Trace: %.3f"
                          % (Tc1Error, Ts1Error, Tc2Error, Ts2Error, Tc3Error, Ts3Error, Tc4Error, Ts4Error, Tc5Error, Ts5Error, Tc6Error, Ts6Error, Tc7Error, Ts7Error, CcError, CsError, RcError, RuError, meanError, Trace))
                    Ap += 'Tc1Error: %.3f, Ts1Error: %.0f, Tc2Error: %.3f, Ts2Error: %.0f, Tc3Error: %.3f, Ts3Error: %.3f, Tc4Error: %.3f, Ts4Error: %.3f, Tc5Error: %.3f, Ts5Error: %.3f, Tc6Error: %.3f, Ts6Error: %.0f, Tc7Error: %.3f, Ts7Error: %.0f, CcError: %.1f%%, CsError: %.1f%%, RcError: %.1f%%, RuError: %.1f%%, MeanError: %.3f, Trace: %.3f' % (
                            Tc1Error, Ts1Error, Tc2Error, Ts2Error, Tc3Error, Ts3Error, Tc4Error, Ts4Error, Tc5Error, Ts5Error, Tc6Error, Ts6Error, Tc7Error, Ts7Error, CcError, CsError, RcError, RuError, meanError, Trace) + "\n"

                    for generation in range(1, 20):
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

                            x_k_1,P_k_1, Breakflag = ekf.measurement_update(States_True[["Ts1", "Ts5"]].values[k, :],P_k_1,x_k_1,Q,R,jacobian_measurement_model)
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
                        meanError = np.mean([np.abs(CcError), np.abs(CsError), np.abs(RcError), np.abs(RuError)])

                        print(
                            "Iterations: %d, Tc1Error: %.3f, Ts1Error: %.0f, Tc2Error: %.3f, Ts2Error: %.0f, Tc3Error: %.3f, Ts3Error: %.3f, Tc4Error: %.3f, Ts4Error: %.3f, Tc5Error: %.3f, Ts5Error: %.0f, Tc6Error: %.3f, Ts6Error: %.3f, Tc7Error: %.3f, Ts7Error: %.0f, CcError: %.1f%%, CsError: %.1f%%, RcError: %.1f%%, RuError: %.1f%%, MeanError: %.3f, Trace: %.3f"
                            % (generation, Tc1Error, Ts1Error, Tc2Error, Ts2Error, Tc3Error, Ts3Error, Tc4Error, Ts4Error, Tc5Error, Ts5Error, Tc6Error, Ts6Error, Tc7Error, Ts7Error, CcError, CsError, RcError, RuError, meanError, Trace))
                        Ap += 'Iterations: %d, Tc1Error: %.3f, Ts1Error: %.0f, Tc2Error: %.3f, Ts2Error: %.0f, Tc3Error: %.3f, Ts3Error: %.3f, Tc4Error: %.3f, Ts4Error: %.3f, Tc5Error: %.3f, Ts5Error: %.3f, Tc6Error: %.3f, Ts6Error: %.0f, Tc7Error: %.3f, Ts7Error: %.0f, CcError: %.1f%%, CsError: %.1f%%, RcError: %.1f%%, RuError: %.1f%%, MeanError: %.3f, Trace: %.3f' % (
                            generation, Tc1Error, Ts1Error, Tc2Error, Ts2Error, Tc3Error, Ts3Error, Tc4Error, Ts4Error, Tc5Error, Ts5Error, Tc6Error, Ts6Error, Tc7Error, Ts7Error, CcError, CsError, RcError, RuError, meanError, Trace) + "\n"

                    k_vec = np.arange(k+1)
                    plt.plot(k_vec, Ts4_pred, c='g')  # Ts4 predicted
                    plt.plot(k_vec, Ts4_true, c='r')  # Ts4 true
                    plt.plot(k_vec, Tc4_pred, c='b')  # Tc4 predicted
                    plt.plot(k_vec, Tc4_true, c='y')  # Tc4 true
                    plt.xlabel('Samples $n$')
                    plt.ylabel(r'Temperature $T$ (°C)')
                    plt.title(r'Temperature $T$ simulation vs estimation')
                    plt.show()

                    # Plot error between Ts4_pred and Ts4_true (core temperature error)
                    plt.figure()
                    error_vec = Ts4_pred - Ts4_true  # Core temperature error
                    plt.plot(k_vec, error_vec, color='red')
                    plt.xlabel('Samples $n$')
                    plt.ylabel(r'Surface temperature error $\Delta Ts4$ (°C)')
                    plt.title(r'Surface temperature error $Ts4$')
                    plt.ylim([-1, 1])
                    plt.grid(True)
                    plt.show()

                    # Plot error between Tc4_pred and Tc4_true (core temperature error)
                    plt.figure()
                    error_vec = Tc4_pred - Tc4_true  # Core temperature error
                    plt.plot(k_vec, error_vec, color='red')
                    plt.xlabel('Samples $n$')
                    plt.ylabel(r'Core temperature error $\Delta Tc4$ (°C)')
                    plt.title(r'Core temperature error $Tc4$')
                    plt.ylim([-1, 1])
                    plt.grid(True)
                    plt.show()

                    out_dir = 'thermal model/outptKF/'
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir)
                    NAME = out_dir + 'EstimationPerformance.txt'
                    with open(NAME, 'w') as f:  
                        f.write(Ap)








"""
==================================================================================================
Author: Yangxiao Xiang @ CityU , yaxiang@cityu.edu.hk
Note:
==================================================================================================
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
np.random.seed(1234)
random.seed(1234)


class RLSIdentification:
    def __init__(self, num_params, theta, forgetting_factor=0.99, delta=1.0):
        theta = theta.reshape(-1, 1)
        self.num_params = num_params
        self.forgetting_factor = forgetting_factor
        self.theta = theta  # Parameter vector to be estimated
        self.P = np.eye(num_params) * delta  # Covariance matrix, initialized large

    def update(self, phi, y):
        # phi is the input, theta is the parameters to be identified
        phi = phi.reshape(-1, 1)
        y = np.array([[y]])

        # Gain computation
        P_phi = np.dot(self.P, phi)
        gain = P_phi / (self.forgetting_factor + np.dot(phi.T, P_phi))

        # Update parameter estimates
        error = y - np.dot(phi.T, self.theta)
        print(error)
        self.theta += gain * error

        # Update covariance matrix
        self.P = (self.P - np.dot(gain, phi.T @ self.P)) / self.forgetting_factor

        return self.theta


def dataProcess(dataframe, ib_integration):

    output_x = []
    output_y = []
    data_ib = dataframe["ib"].values
    data_vb = dataframe["vb1"].values

    for i in range(1, len(dataframe)-1):
        if data_ib[i] != data_ib[i+1]:
            continue
        elif data_ib[i] != data_ib[i-1]:
            continue
        else:
            output_x.append([data_vb[i], data_ib[i], ib_integration[i, 0], 1])
            output_y.append([data_vb[i+1]])
    output_x = np.array(output_x)
    output_y = np.array(output_y)

    return output_x, output_y


def main():
    # Simulate some data
    dt = 0.1
    df = pd.read_csv("./Simulation_data/targetWave.csv")
    starttime = 0
    endtime = -50
    df = df.iloc[starttime:endtime, :].reset_index()
    ib = df[["ib"]].values
    ib_integration = [0]
    for current in ib:
        tmp = ib_integration[-1] + current * dt
        ib_integration.append(tmp)
    ib_integration = np.array(ib_integration[1:])
    X, Y = dataProcess(df, ib_integration)

    # True parameters
    R_true = 5e-2 * 0.85
    R1_true = 1.3e0 * 1.2
    C1_true = 5e2 * 1.1
    C_true = 3e0 * 3600
    SOCinit_true = 0.81
    k = 0.824
    b = 3.347
    theta1_true = (R1_true * C1_true / dt - 1) / (R1_true * C1_true / dt)
    theta2_true = (1 * R1_true * C1_true * k / C_true + R_true + R1_true) / (R1_true * C1_true / dt)
    theta3_true = (k / C_true) / (R1_true * C1_true / dt)
    theta4_true = (k * SOCinit_true + b) / (R1_true * C1_true / dt)

    # Initial guess
    R_init = R_true * 1.2
    R1_init = R1_true * 0.8
    C1_init = C1_true * 0.9
    C_init = C_true * 0.95
    SOCinit_init = SOCinit_true * 1.1
    # R_init = R_true * 1.
    # R1_init = R1_true * 0.999999
    # C1_init = C1_true * 0.999999
    # C_init = C_true * 0.999999
    # SOCinit_init = SOCinit_true * 1.

    theta1 = (R1_init * C1_init / dt - 1) / (R1_init * C1_init / dt)
    theta2 = (1 * R1_init * C1_init * k / C_init + R_init + R1_init) / (R1_init * C1_init / dt)
    theta3 = (k / C_init) / (R1_init * C1_init / dt)
    theta4 = (k * SOCinit_init + b) / (R1_init * C1_init / dt)

    # Recursive Least Squares initialization
    rls = RLSIdentification(num_params=4, theta=np.array([theta1, theta2, theta3, theta4]), forgetting_factor=1, delta=1000000)

    Error1 = abs(rls.theta[0] / theta1_true) * 100 - 100
    Error2 = abs(rls.theta[1] / theta2_true) * 100 - 100
    Error3 = abs(rls.theta[2] / theta3_true) * 100 - 100
    Error4 = abs(rls.theta[3] / theta4_true) * 100 - 100
    print("Error1: %.3f, Error2: %.3f, Error3: %.3f, Error4: %.3f" % (Error1, Error2, Error3, Error4))

    # Identified parameters over time
    estimated_params = np.zeros((len(Y), 4))

    # Recursive identification process
    for t in range(0, len(Y)):
        # Feature vector
        input = X[t, :]

        # Update RLS
        rls_theta = rls.update(input, Y[t, 0])

        # evaluation
        est_theta1 = rls_theta[0]
        est_theta2 = rls_theta[1]
        est_theta3 = rls_theta[2]
        est_theta4 = rls_theta[3]

        Error1 = abs(est_theta1 / theta1_true) * 100 - 100
        Error2 = abs(est_theta2 / theta2_true) * 100 - 100
        Error3 = abs(est_theta3 / theta3_true) * 100 - 100
        Error4 = abs(est_theta4 / theta4_true) * 100 - 100
        fitting_error = np.mean(np.abs(Y - np.dot(X, rls_theta))) * 1e10
        print("Error1: %.3f, Error2: %.3f, Error3: %.3f, Error4: %.3f, fitting error: %.3f" % (Error1, Error2, Error3, Error4, fitting_error))

        # estimated_params[t, :] = rls_theta.ravel()

    # Plotting the results
    # plt.figure(figsize=(10, 6))
    #
    # plt.subplot(2, 1, 1)
    # plt.plot(time, V_b, label="True Voltage")
    # plt.title("Battery Voltage")
    # plt.xlabel("Time [s]")
    # plt.ylabel("Voltage [V]")
    #
    # plt.subplot(2, 1, 2)
    # plt.plot(time, estimated_params[:, 0], label="VOC")
    # plt.plot(time, estimated_params[:, 1], label="R")
    # plt.plot(time, estimated_params[:, 2], label="R1")
    # plt.plot(time, estimated_params[:, 3], label="C1")
    # plt.title("Estimated Parameters")
    # plt.xlabel("Time [s]")
    # plt.ylabel("Parameter Value")
    # plt.legend()
    #
    # plt.tight_layout()
    # plt.show()


if __name__ == "__main__":
    main()
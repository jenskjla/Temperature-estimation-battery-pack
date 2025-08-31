"""
Author - Kailash Nagarajan
Date - 3-Oct-2020
email-ID - kailashnagarajan@gmail.com 
"""

import numpy as np 

"""
Extended Kalman Filter - Linearizing Non-Linear Motion Models using Jacobians.
1. Predict the states using the Motion Model and Jacobian of the Motion Model.
2. Correct the states using the Kalman Gain, Sensor Values and Jacobian of the Measurement Model.
"""

class EKF:

	def __init__(self):
		pass

	def measurement_update(self,y_meas,P_k_1,x_k_1,Q,R,jacobian):
		"""
		The Measurement Update Function Updates the predicted states using the Sensor Values and Kalman Gain.
		
		Parameters : 
		P_k_1 -> Predicted Co-variance Matrix - Array - (18,18)
		x_k_1 -> Predicted States Matrix - Array - (18,1)
		Q - Process Noise Co-Variance - Array - (18,18)
		R - Measurement Noise Co-variance - Array - (2,2)
		jacobian - Function that computes the measurement jacobians (Both Noise and Last State) - Function

		Returns : 
		x_k - Corrected State - Array - (18,1)
		P_k - Corrected Covariance - Array - (18,18)

		"""

		# Assigning predicted state values to respective states.


		H_k,M_k,y_out = jacobian(x_k_1)

		"""
		H_k - Jacobian w.r.t last state - Array - (2,18).
		M_k - Jacobian w.r.t noise - Array - (2,2).
		"""
		
		# Calculation of Kalman Gain
		tmp = H_k.dot(P_k_1).dot(H_k.T) + R
		tmp = tmp.astype(np.float64)
		det = np.linalg.det(tmp)
		if det != 0:
			K_k = P_k_1.dot(H_k.T).dot(np.linalg.inv(tmp))
			# K_k = P_k_1.dot(H_k.T).dot(np.linalg.inv(H_k.dot(P_k_1).dot(H_k.T) + M_k.dot(R).dot(M_k.T)))


			# State and Co-Variance Correction Step

			x_k = x_k_1 + K_k.dot(y_meas.reshape(-1, 1) - y_out) # State Matrix

			P_k = (np.identity(len(K_k))-K_k.dot(H_k)).dot(P_k_1) # Co-Variance Matrix

			BreakFlag = False
		else:
			x_k = x_k_1
			P_k = P_k_1
			BreakFlag = True


		return x_k,P_k, BreakFlag
	

	def prediction_step(self,delta_t,ControlSet,x_k_1,P_k_1,Q,jacobian,motion_model, Qv, k):

		"""
		The Prediction Update Function, Predicts the next step using the previous state and control and 
		the Motion and Noise Model.
		
		Parameters : 

		delta_t -> Time Period - float val  
		x_k_1 -> Current State - Array - (18,1)
		P_k_1 -> Current Covariance - Array - (18,18)
		jacobian -> Function that computes the process jacobians (Both Noise and Last State) - Function
		motion_model -> Function that computes the output states for a given input states - Function 

		Returns : 
		x_k_1 : Predicted States - Array - (18,1)
		P_k_1 : Predicted states - Array - (18,18)

		"""
		# Jacobian of Motion Model w.r.t last state and Noise
		F, L = jacobian(x_k_1, ControlSet[0], delta_t, Qv, k)

		# Motion Model Returns the states [x,y,theta]
		x_k_1 = motion_model(ControlSet,x_k_1,delta_t, Qv, k)

		# Predicted Co-Variance
		P_k_1 = F.dot((P_k_1).dot(F.T)) + Q
		# P_k_1 = F.dot((P_k_1).dot(F.T)) + L.dot((Q).dot(L.T))

		return x_k_1,P_k_1 

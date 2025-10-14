# 🔋 Surface & Core Temperature Estimation of Large-Scale Li-Ion Battery Packs with Reduced Sensors

This repository contains the full code, models, and documentation developed for my Bachelor's Final Year Project at **City University of Hong Kong**, titled:  
**“Surface and Core Temperature Estimation of Large-Scale Lithium-Ion Battery Packs with Reduced Sensors.”**

The work focuses on building a **reduced-sensor temperature estimation framework** for lithium-ion battery packs, combining **electrical-thermal coupled modeling**, **Genetic Algorithm (GA)** parameter identification, and **state estimation algorithms** such as the **Extended Kalman Filter (EKF)**.

---

## 🚀 Motivation

With the increasing demand for **electric vehicles** and **energy storage systems**, precise temperature monitoring of lithium-ion batteries is critical for:

- 🔒 **Safety:** Preventing thermal runaway and degradation.  
- ⚙️ **Performance:** Maintaining consistent internal temperature distribution.  
- 💰 **Cost reduction:** Minimizing the number of physical sensors required.  

Traditional battery systems rely on multiple embedded temperature sensors, increasing complexity and cost.  
This project proposes a **model-based approach** capable of estimating **core and surface temperatures** using only a subset of sensors, maintaining high accuracy through advanced filtering methods.

---

## 🎯 Objectives

1. **Develop** a coupled **electrical–thermal model** of lithium-ion cells in MATLAB/Simulink.  
2. **Identify** optimal parameters (Rc, Cc, Cs, Ru, Rcc) using a **Genetic Algorithm**.  
3. **Design and implement** temperature estimation algorithms (RLS, KF, EKF) in Python.  
4. **Validate** results with both **simulation and experimental data** from a 7-cell pack.  
5. **Evaluate** the algorithms by comparing:
   - Mean Absolute Error (MAE)    
   - Sensor-reduction effectiveness  

---

## 🧩 Project Roadmap

| Phase | Description | Tools |
|:------|:-------------|:------|
| **1. Modeling** | Creation of electrical & lumped-thermal models | MATLAB / Simulink |
| **2. Parameter Identification** | GA optimization of thermal–electrical parameters | Python |
| **3. Data Acquisition** | Experimental current, voltage, and temperature measurements | CSV logs |
| **4. Temperature Estimation** | RLS and EKF implementations | Python |
| **5. Validation** | Comparison of estimated vs measured results | MATLAB / Python |
| **6. Reporting** | Documentation and visualization of results | LaTeX / Markdown |

---

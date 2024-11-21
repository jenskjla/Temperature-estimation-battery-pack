"""
==================================================================================================
Author: Yangxiao Xiang @ CityU , yaxiang@cityu.edu.hk
Note:
     * The code is developed based on the github repository
       https://github.com/ezstoltz/genetic-algorithm, where more useful information can be found therein.
==================================================================================================
"""

import numpy as np
import random
import operator
import pandas as pd
import time
import shutil
np.random.seed(1234)
random.seed(1234)

import os

class parameter_information():
    def __init__(self):
        self.parameter_names = ["R", "R1", "C1", "C", "SOCinit"]
        self.target_names = ["loss"]
        self.True_ParamValue = {"R": 5e-2 * 0.85, "R1": 1.3e0 * 1.2, "C1": 5e2 * 1.1, "C": 3e0, "SOCinit": 0.81}
        self.upper_limits = {"R": self.True_ParamValue["R"]*2,
                             "R1": self.True_ParamValue["R1"] * 2,
                             "C1": self.True_ParamValue["C1"] * 1.3,
                             "C": self.True_ParamValue["C"] * 1.3,
                             "SOCinit": self.True_ParamValue["SOCinit"] * 1.3}
        self.lower_limits = {"R": self.True_ParamValue["R"]*0.2,
                             "R1": self.True_ParamValue["R1"] * 0.2,
                             "C1": self.True_ParamValue["C1"] * 0.7,
                             "C": self.True_ParamValue["C"] * 0.7,
                             "SOCinit": self.True_ParamValue["SOCinit"] * 0.7}
        self.dataDir = "./outputGA"
        self.history_Population_pandas_form = pd.DataFrame(columns=self.parameter_names + self.target_names + ["fitness", "waiting_time"])

        df = pd.read_csv("./electrical_model/Simulation_data/targetWave.csv")
        starttime = 0
        endtime = -50
        df = df.iloc[starttime:endtime, :].reset_index()
        self.ib = df["ib"].values
        self.vb = df["vb1"].values
        self.dt = 0.1


def getFitness(loss):
    return - (loss)

def calculate_error(estimated_values, ParamNames, TargetParamValue):

    ground_truth = np.array([TargetParamValue[ParamName] for ParamName in ParamNames]).reshape(1, -1)
    ground_truth = np.tile(ground_truth, (len(estimated_values), 1))

    error = np.abs(estimated_values / ground_truth * 100 - 100)
    meanError = np.mean(error, axis=1)
    return meanError


class geneticAlgrithm(parameter_information):
    def __init__(self,):
        super(geneticAlgrithm, self).__init__()
        self.max_value = 1e6

    def initialPopulation(self,popSize):
        # randomly generate children from the whole domain
        population = []
        for i in range(0, popSize):
            new_pop = {key: None for key in self.target_names}
            for p in self.parameter_names:
                new_pop[p] = (self.upper_limits[p]-self.lower_limits[p])*random.random() + self.lower_limits[p]
            population.append(new_pop)
        return ([], population)

    def rankRoutes(self,population):
        fitnessResults = {}
        for i in range(0,len(population)):
            fitnessResults[i] = getFitness(population[i]["loss"])
        return sorted(fitnessResults.items(), key=operator.itemgetter(1), reverse=True)

    def selection(self,popRanked, eliteSize):
        selectionResults = []
        df = pd.DataFrame(np.array(popRanked), columns=["Index", "Fitness"])
        df['cum_sum'] = df.Fitness.cumsum()
        df['cum_perc'] = 100 * df.cum_sum / df.Fitness.sum()

        for i in range(0, eliteSize):
            selectionResults.append(popRanked[i][0])
        for i in range(0, len(popRanked) - eliteSize):
            pick = 100 * random.random()
            for i in range(0, len(popRanked)):
                if pick <= df.iat[i, 3]:
                    selectionResults.append(popRanked[i][0])
                    break
        return selectionResults

    def matingPool(self,population, selectionResults):
        matingpool = []
        for i in range(0, len(selectionResults)):
            index = selectionResults[i]
            matingpool.append(population[index])
        return matingpool

    def breed(self,parent1, parent2):
        child1 = {key: None for key in self.target_names}
        child2 = {key: None for key in self.target_names}
        for feature in self.parameter_names:
            while_jump_flag = False
            while not while_jump_flag:
                u = random.random()
                if u > 0.5:
                    gamma = (1 / (2 * (1 - u))) ** 0.5
                else:
                    gamma = (2 * u) ** 0.5
                new_data1 = 0.5 * ((1 + gamma) * parent1[feature] + (1 - gamma) * parent2[feature])
                new_data2 = 0.5 * ((1 - gamma) * parent1[feature] + (1 + gamma) * parent2[feature])
                if new_data1 >= self.lower_limits[feature] and new_data1 <= self.upper_limits[feature] and new_data2 >= self.lower_limits[feature] and new_data2 <= self.upper_limits[feature]:
                    while_jump_flag = True
            child1[feature] = new_data1
            child2[feature] = new_data2

        return child1, child2

    def breedPopulation(self,matingpool, eliteSize):
        children = []
        length = len(matingpool) - eliteSize
        pool = random.sample(matingpool, len(matingpool))

        for i in range(0, eliteSize):
            children.append(matingpool[i]) # add elite children

        children_another_group = []
        for i in range(0, length):
            child1, child2 = self.breed(pool[i], pool[len(matingpool) - i - 1])
            children_another_group.append(child1)
            children_another_group.append(child2)

        children += random.sample(children_another_group, length)
        return children

    def mutate(self, individual):
        for feature in self.parameter_names:
            p = random.random()
            if p > 0.75:
                individual[feature] = \
                    (self.upper_limits[feature] - individual[feature]) * random.random() + individual[feature]
            elif p < 0.25:
                individual[feature] = \
                    individual[feature] - (individual[feature] - self.lower_limits[feature]) * random.random()
            else:
                pass
        return individual

    def mutatePopulation(self,population, mutationRate, eliteSize):
        mutatedPop = []
        nochangedPop = []
        for ind in range(0, len(population)):
            if ind < eliteSize:  # reserve best children
                nochangedPop.append(population[ind])
            else:
                if (random.random() < mutationRate):
                    mutatedPop.append(self.mutate(population[ind]))
                else:
                    mutatedPop.append(population[ind])
        return (nochangedPop, mutatedPop)

    def getNextGeneration(self, currentGen, eliteSize, mutationRate):
        popRanked = self.rankRoutes(currentGen)
        selectionResults = self.selection(popRanked, eliteSize) # Roulette selection
        matingpool = self.matingPool(currentGen, selectionResults)
        children = self.breedPopulation(matingpool, eliteSize) #  simulated binary crossover ref: K. Deb and R. B. Agrawal, “Simulated binary crossover for continuous search space,” Complex Syst., vol. 9, no. 2, pp. 115–148, 1995.
        nextGeneration = self.mutatePopulation(children, mutationRate, eliteSize) # simplest uniform mutation
        return nextGeneration

    def soc_voc(self, soc):
        return 0.824 * soc + 3.347

    def calculate_loss(self, param):

        R = param["R"]
        R1 = param["R1"]
        C1 = param["C1"]
        C = param["C"]
        SOCinit = param["SOCinit"]

        vo_list = []
        SOC_list = [SOCinit]
        voc_list = [self.soc_voc(SOC_list[0])]
        i_integration = 0
        v1_integration = 0
        SOC_integration = 0
        v1_init = 0
        for n in range(0, len(self.ib)):
            v1 = (i_integration - v1_integration / R1) / C1 + v1_init
            v1_integration = v1_integration + v1 * self.dt
            i_integration = i_integration + self.ib[n] * self.dt
            SOC_integration = SOC_integration + self.ib[n] * self.dt / (C * 3600)
            SOC = SOC_integration + SOC_list[0]
            voc = self.soc_voc(SOC)
            SOC_list.append(SOC)
            voc_list.append(voc)
            vo_list.append(voc + v1 + R * self.ib[n])
        SOC_list = SOC_list[1:]
        voc_list = voc_list[1:]
        vo_list = np.array(vo_list)

        loss = np.sum(np.square(self.vb - vo_list))

        return loss

    def getPerformance(self, population):

        population_new = []
        for pop in population:
            pop[self.target_names[0]] = self.calculate_loss(pop)
            population_new.append(pop)
        return population_new

if __name__ == '__main__':

    popSize = 250
    eliteSize = round(popSize * 0.15) # save the top eliteSize children without breed and mutation
    mutationRate = 0.15
    generations = 300 # run "generations" iterations

    GA = geneticAlgrithm()
    if os.path.exists(GA.dataDir):
        shutil.rmtree(GA.dataDir)
    os.mkdir(GA.dataDir)

    time_start = time.time()

    training_log = []
    for generation_cnt in range(0, generations):
        if generation_cnt == 0:
            (noChangedPop, toBeTestedPop) = GA.initialPopulation(popSize)
        else:
            (noChangedPop, toBeTestedPop) = GA.getNextGeneration(pop, eliteSize, mutationRate)
        toBeTestedPop = GA.getPerformance(toBeTestedPop)
        pop = noChangedPop + toBeTestedPop

        population_pandas_form = pd.DataFrame(columns=GA.parameter_names + GA.target_names + ["fitness", "error", "waiting_time"])
        population_pandas_form = pd.concat([population_pandas_form, pd.DataFrame(pop)], ignore_index=True)

        population_pandas_form["fitness"] = getFitness(population_pandas_form["loss"].values)
        population_pandas_form["error"] = calculate_error(population_pandas_form[GA.parameter_names].values, GA.parameter_names, GA.True_ParamValue)
        time_end = time.time()
        population_pandas_form["waiting_time"] = (time_end-time_start)*np.ones((len(population_pandas_form),))
        population_pandas_form = population_pandas_form.sort_values(by="fitness", ascending=False, inplace=False)

        if generation_cnt % 50 == 0:
            population_pandas_form.to_csv(os.path.join(GA.dataDir, "%d.csv" % (generation_cnt)), index=False)
        if generation_cnt % 5 == 0:
            training_log.append({"generation": generation_cnt, "loss": population_pandas_form["loss"].values[0], "mean": population_pandas_form["error"].values[0], "waiting_time": population_pandas_form["waiting_time"].values[0]})
        print("Generation: %d, loss: %f, error: %f, time: %f" %(generation_cnt, population_pandas_form["loss"].values[0], population_pandas_form["error"].values[0], population_pandas_form["waiting_time"].values[0]))
    df = pd.DataFrame(training_log)
    df.to_csv("./electrical_model/outputGA/results_GA.csv", index=False)
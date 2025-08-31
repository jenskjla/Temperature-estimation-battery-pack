import numpy as np
import random
import operator
import pandas as pd
import os
import shutil
import time

# Reproducibility
np.random.seed(1234)
random.seed(1234)

class parameter_information:
    """
    Handles data loading and parameter bounds for multi-battery GA.
    """
    def __init__(self, num_batteries=7, data_file="targetWave.csv"):
        self.num_batteries = num_batteries
        # Base parameter names for each battery
        base_names = ["R", "R1", "C1", "C", "SOCinit"]
        # Dynamic parameter names: name_1, ..., name_N
        self.parameter_names = [f"{name}{i+1}" for i in range(num_batteries) for name in base_names]
        self.target_names = ["loss"]

        # Synthetic true values (for error calc)
        default_true = {"R": 5e-2, "R1": 1.3e0, "C1": 5e2, "C": 3e3, "SOCinit": 0.1}
        self.True_ParamValue = {}
        self.upper_limits = {}
        self.lower_limits = {}
        for i in range(num_batteries):
            for name, val in default_true.items():
                key = f"{name}{i+1}"
                self.True_ParamValue[key] = val
                if name in ["R", "R1"]:
                    self.upper_limits[key] = val * 2
                    self.lower_limits[key] = val * 0.2
                else:
                    self.upper_limits[key] = val * 1.3
                    self.lower_limits[key] = val * 0.7

        # Load data from project directory electrical_model/Simulation_data/
        data_path = os.path.join(os.path.dirname(__file__), "Simulation_data", data_file)
        df = pd.read_csv(data_path)
        self.ib = df['ib'].values
        self.vb = np.vstack([df[f"vo{i+1}"].values for i in range(num_batteries)]).T
        if 't' in df.columns:
            dt_arr = np.diff(df['t'].values)
            self.dt = float(np.median(dt_arr)) if len(dt_arr)>0 else 0.1
        else:
            self.dt = 0.1

        # Output directory
        self.dataDir = "electrical_model/outputGA"


def getFitness(losses):
    total = np.sum(losses)
    return -total

class geneticAlgorithm(parameter_information):
    def __init__(self, num_batteries=7, data_file="targetWave_elec.csv"):
        super().__init__(num_batteries, data_file)

    def initialPopulation(self, popSize):
        population = []
        for _ in range(popSize):
            individual = {key: self.lower_limits[key] + (self.upper_limits[key] - self.lower_limits[key]) * random.random()
                          for key in self.parameter_names}
            population.append(individual)
        return ([], population)

    def rankRoutes(self, population):
        fitnessResults = {i: getFitness(ind['losses']) for i, ind in enumerate(population)}
        return sorted(fitnessResults.items(), key=operator.itemgetter(1), reverse=True)

    def selection(self, popRanked, eliteSize):
        df = pd.DataFrame(popRanked, columns=['Index', 'Fitness'])
        df['cum_sum'] = df.Fitness.cumsum()
        df['cum_perc'] = 100 * df.cum_sum / df.Fitness.sum()

        selectionResults = [idx for idx, _ in popRanked[:eliteSize]]
        for _ in range(len(popRanked) - eliteSize):
            pick = 100 * random.random()
            for idx, row in df.iterrows():
                if pick <= row['cum_perc']:
                    selectionResults.append(int(row['Index']))
                    break
        return selectionResults

    def matingPool(self, population, selectionResults):
        return [population[i] for i in selectionResults]

    def breed(self, parent1, parent2):
        child1, child2 = {}, {}
        for feature in self.parameter_names:
            u = random.random()
            beta = (1 / (2 * (1 - u)))**0.5 if u > 0.5 else (2 * u)**0.5
            val1 = 0.5*((1+beta)*parent1[feature] + (1-beta)*parent2[feature])
            val2 = 0.5*((1-beta)*parent1[feature] + (1+beta)*parent2[feature])
            lo, hi = self.lower_limits[feature], self.upper_limits[feature]
            child1[feature] = min(max(val1, lo), hi)
            child2[feature] = min(max(val2, lo), hi)
        return child1, child2

    def breedPopulation(self, matingpool, eliteSize):
        children = matingpool[:eliteSize]
        pool = random.sample(matingpool, len(matingpool))
        for i in range(len(matingpool) - eliteSize):
            c1, c2 = self.breed(pool[i], pool[-i-1])
            children.append(c1)
            if len(children) < len(matingpool):
                children.append(c2)
        return children[:len(matingpool)]

    def mutate(self, individual, mutationRate):
        for feature in self.parameter_names:
            if random.random() < mutationRate:
                individual[feature] = self.lower_limits[feature] + (self.upper_limits[feature] - self.lower_limits[feature]) * random.random()
        return individual

    def mutatePopulation(self, population, mutationRate, eliteSize):
        new_pop = population[:eliteSize]
        for ind in population[eliteSize:]:
            new_pop.append(self.mutate(ind.copy(), mutationRate))
        return new_pop

    def getNextGeneration(self, pop, eliteSize, mutationRate):
        # Selection and breeding
        popRanked = self.rankRoutes(pop)
        selectionResults = self.selection(popRanked, eliteSize)
        matingpool = self.matingPool(pop, selectionResults)
        children = self.breedPopulation(matingpool, eliteSize)
        nextGen = self.mutatePopulation(children, mutationRate, eliteSize)
        # Split elites (no change) and rest for evaluation
        noChangePop = nextGen[:eliteSize]
        toBeTested = nextGen[eliteSize:]
        return noChangePop, toBeTested

    def soc_voc(self, soc):
        return 0.824 * soc + 3.347

    def calculate_losses(self, param):
        losses = []
        for i in range(self.num_batteries):
            R   = param[f"R{i+1}"]
            R1  = param[f"R1{i+1}"]
            C1  = param[f"C1{i+1}"]
            C   = param[f"C{i+1}"]
            SOCinit = param[f"SOCinit{i+1}"]

            ib = self.ib
            true_v = self.vb[:, i]
            sim_v = []

            v1_int = 0.0
            i_int  = 0.0
            soc_int= 0.0
            prev_v1= 0.0
            for n in range(len(ib)):
                v1 = (i_int - v1_int/R1)/C1 + prev_v1
                v1_int += v1 * self.dt
                i_int  += ib[n] * self.dt
                soc_int+= ib[n] * self.dt / (C * 3600)
                soc    = SOCinit + soc_int
                sim_v.append(self.soc_voc(soc) + v1 + R*ib[n])
                prev_v1 = v1

            loss_i = (np.sum((true_v - np.array(sim_v))**2))/len(self.ib)
            losses.append(loss_i)
        return losses

    def getPerformance(self, population):
        for ind in population:
            ind['losses'] = self.calculate_losses(ind)
            # also store total for convenience
            ind['total_loss'] = sum(ind['losses'])
        return population

if __name__ == '__main__':
    popSize = 250
    eliteSize = int(0.15 * popSize)
    mutationRate = 0.15
    generations = 100

    GA = geneticAlgorithm(num_batteries=7, data_file="targetWave.csv")
    if os.path.exists(GA.dataDir): shutil.rmtree(GA.dataDir)
    os.makedirs(GA.dataDir, exist_ok=True)

    # Initialize best-loss trackers
    best_losses = [np.inf] * GA.num_batteries

    start_time = time.time()
    noChangePop, toBeTested = GA.initialPopulation(popSize)
    for gen in range(generations):
        if gen > 0:
            noChangePop, toBeTested = GA.getNextGeneration(pop, eliteSize, mutationRate)
        toBeTested = GA.getPerformance(toBeTested)
        pop = noChangePop + toBeTested

        fitness = [getFitness(ind['losses']) for ind in pop]
        best_idx = int(np.argmax(fitness))
        best = pop[best_idx]
        # Print per-battery losses and parameters
        loss_str = ", ".join([f"B{i+1}:{l:.4e}" for i, l in enumerate(best['losses'])])
        print(f"Gen {gen}: Losses = [{loss_str}]")
        param_strs = []
        for i in range(GA.num_batteries):
            params_i = {name: best[name] for name in GA.parameter_names if name.endswith(str(i+1))}
            param_strs.append(f"B{i+1}:" + ";".join([f"{k}={v:.3e}" for k, v in params_i.items()]))
        

        # Check and update best per-battery files
        for i, loss_i in enumerate(best['losses']):
            if loss_i < best_losses[i]:
                best_losses[i] = loss_i
                # Save parameters for battery i+1
                params_i = {k: best[k] for k in GA.parameter_names if k.endswith(str(i+1))}
                df_i = pd.DataFrame([params_i])
                file_path = os.path.join(GA.dataDir, f"battery{i+1}_best.csv")
                df_i.to_csv(file_path, index=False)

    print("GA completed.")
    print("Final best losses per battery:")
    for i, l in enumerate(best_losses):
        print(f"Battery {i+1}: {l:.4e}")
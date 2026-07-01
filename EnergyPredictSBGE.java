import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

public class EnergyPredictSBGE {

    // GE chromosome parameters
    static final int CODON_LENGTH = 20; // number of codons per individual
    static final int CODON_MAX = 255; // codons are integers
    static final int MAX_WRAPS = 3; // max chromosome wraps

    // GE mapping parameters
    static final int MAPPING_MAX_DEPTH = 10; // max recursion depth during mapping

    // Evolutionary parameters
    static final int POP_SIZE = 100; // population size
    static final int MAX_GEN = 50; // generations per run
    static final int TOURNAMENT_SIZE = 5; // tournament pool size

    static final double MUTATION_PROB = 0.01;

    // Adaptive mutation thresholds
    static final double MUTATION_HIGH = 0.30; // explore, that means crossover = 0.80
    static final double MUTATION_LOW = 0.20; // exploit, that means crossover = 0.60

    // Structural similarity thresholds
    static final double SIM_THRESHOLD = 0.75; // similarity cutoff
    static final double CVG_FRACTION = 0.60; // converged

    static final double ALPHA = 0.4; // weight given to structural diversity

    // Time-series load configuration
    static final int N_loadS = 7; // number of load inputs
    static final int SLOTS_PER_DAY = 96; // 96 x 15-min intervals per day

    static final int SEED = 200;
    static final int NO_OF_RUNS = 10;
    static Random rng = new Random(SEED);

    // Success criterion
    static final double SUCCESS_THRESHOLD = 0.05; // MSE target

    static final String CSV_FILENAME = "Residential_Energy_Dataset_UK- 2014-2020.csv";

    // Train / test splits [G]
    static double[][] X_train, X_test;
    static double[] y_train, y_test;
    static final int n = 7;

    static int numForDataset = 1; // 0 = runs reads half of the dataset, 1 = runs 10k rows of the dataset, 2 going
                                  // up = runs the full dataset
    // Shared instance
    static final ExpressionEvaluator evaluator = new ExpressionEvaluator();

    // - Global tracking -
    static List<Integer> successGenerations = new ArrayList<>();
    static long startTime = 0;
    static long endTime = 0;

    public static void main(String[] args) {

        printBanner();

        // Build train/test dataset [G][K]
        buildDataset();

        // Initialise the BNF grammar with the correct number of load variables
        BNFGrammar.init(N_loadS, MAPPING_MAX_DEPTH, MAX_WRAPS);

        printConfig();

        tick();

        // Accumulators across all runs [F]
        List<Double> allBestTrain = new ArrayList<>();
        List<Double> allBestTest = new ArrayList<>();
        List<Double> allAvgMSE = new ArrayList<>();
        List<Long> allRunTimes = new ArrayList<>();
        List<String> allEquations = new ArrayList<>();
        // Keep the best individual object for final prediction demo
        GEIndividual overallBestInd = null;
        double overallBestMSE = Double.MAX_VALUE;

        // RUN LOOP
        for (int run = 0; run < NO_OF_RUNS; run++) {

            rng = new Random(SEED + run);

            System.out.printf("  RUN %2d / %d    seed=%d%n",
                    run, NO_OF_RUNS - 1, SEED + run);
            System.out.println(sep('-', 70));

            long runStart = System.currentTimeMillis();
            double bestMSERun = Double.MAX_VALUE;
            GEIndividual bestIndRun = null;
            double sumGenBestMSE = 0.0;
            int successGen = -1;

            // Initialise population
            List<GEIndividual> population = initPopulation();

            // Map all genotypes to phenotypes, evaluate fitness
            mapAndEvaluateAll(population);

            // GENERATIONAL LOOP
            for (int gen = 0; gen < MAX_GEN; gen++) {

                GEIndividual bestThisGen = getBest(population);
                double bestMSEGen = bestThisGen.fitness;

                // Track success
                if (successGen == -1 && bestMSEGen < SUCCESS_THRESHOLD)
                    successGen = gen;

                sumGenBestMSE += bestMSEGen;

                // Keep run-best
                if (bestMSEGen < bestMSERun) {
                    bestMSERun = bestMSEGen;
                    bestIndRun = bestThisGen.copy();
                }

                // STEP 3: Compute phenotype similarity
                int[] simCounts = computeSimilarityMatrix(population);
                int cvgCount = 0;
                for (int sc : simCounts)
                    if (sc > 0)
                        cvgCount++;
                double cvgFrac = (double) cvgCount / POP_SIZE;

                // Adaptive mutation rate
                double curMutProb = (cvgFrac >= CVG_FRACTION)
                        ? MUTATION_HIGH // explore
                        : MUTATION_LOW; // exploit

                int remaining = POP_SIZE - 1;

                int mutSlots = (int) Math.round(remaining * curMutProb);
                int crossSlots = remaining - mutSlots;
                // int mutSlots = remaining - crossSlots;

                List<GEIndividual> newPop = new ArrayList<>();

                // Elitism
                newPop.add(bestThisGen.copy());

                int crossUsed = 0;
                int mutUsed = 0;

                while (newPop.size() < POP_SIZE) {

                    if (crossUsed < crossSlots) {
                        // -- CROSSOVER -
                        GEIndividual pA = selectCombined(population, newPop);
                        GEIndividual pB = selectCombined(population, newPop);

                        GEIndividual[] children = crossover(pA, pB);
                        GEIndividual cA = children[0];
                        GEIndividual cB = children[1];

                        // Re-map phenotypes for modified children
                        BNFGrammar.mapGenotype(cA);
                        BNFGrammar.mapGenotype(cB);
                        evaluateFitness(cA);
                        evaluateFitness(cB);

                        if (newPop.size() < POP_SIZE) {
                            newPop.add(cA);
                            crossUsed++;
                        }
                        if (newPop.size() < POP_SIZE) {
                            newPop.add(cB);
                            crossUsed++;
                        }

                    } else {
                        // -- MUTATION-ONLY REPRODUCTION --
                        GEIndividual parent = selectCombined(population, newPop);
                        GEIndividual child = parent.copy();

                        mutateCodons(child, curMutProb);
                        BNFGrammar.mapGenotype(child);
                        evaluateFitness(child);

                        newPop.add(child);
                        mutUsed++;
                    }
                }

                population = newPop;
            }
            // end generational loop

            // Record run outcome
            successGenerations.add(successGen == -1 ? MAX_GEN + 1 : successGen);

            long runTime = System.currentTimeMillis() - runStart;
            double testMSE = (bestIndRun != null && bestIndRun.valid)
                    ? computeMSE(bestIndRun.phenotype, y_test, X_test)
                    : 1e9;
            double avgMSE = sumGenBestMSE / MAX_GEN;

            allBestTrain.add(bestMSERun);
            allBestTest.add(testMSE);
            allAvgMSE.add(avgMSE);
            allRunTimes.add(runTime);
            allEquations.add(bestIndRun != null ? bestIndRun.phenotype : "INVALID");

            if (bestIndRun != null && bestMSERun < overallBestMSE) {
                overallBestMSE = bestMSERun;
                overallBestInd = bestIndRun.copy();
            }

            System.out.println("\n  " + sep('-', 62));
            System.out.printf("  Run %2d Summary%n", run);
            System.out.println("  " + sep('-', 62));
            System.out.printf("  Best train MSE   : %.6f%n", bestMSERun);
            System.out.printf("  Best test  MSE   : %.6f%n", testMSE);
            System.out.printf("  Avg best-gen MSE : %.6f%n", avgMSE);
            System.out.printf("  Run time (ms)    : %d%n", runTime);
            System.out.println("  Best phenotype   : "
                    + (bestIndRun != null ? bestIndRun.phenotype : "INVALID"));
            System.out.println("  Best codons (first 20): "
                    + Arrays.toString(Arrays.copyOf(
                            bestIndRun != null ? bestIndRun.codons : new int[0], 20)));
            if (successGen >= 0)
                System.out.printf("  First success    : generation %d%n", successGen);
            else
                System.out.println("  First success    : not reached this run");
        }
        // end multi-run loop

        tock();

        // REPORT
        printReport(allBestTrain, allBestTest, allAvgMSE,
                allRunTimes, allEquations);

        // COMPUTATIONAL EFFORT
        computeEffort(successGenerations, POP_SIZE);

        // PREDICTION
        if (overallBestInd != null)
            printPrediction(overallBestInd);
    }

    // DATASET CONSTRUCTION
    static void buildDataset() {
        List<Double> electricLoad = new ArrayList<>();

        File csv = new File(CSV_FILENAME);
        if (csv.exists()) {
            System.out.println("DATASET Reading CSV: " + CSV_FILENAME);
            try (Scanner sc = new Scanner(csv)) {
                sc.nextLine(); // skip header line
                while (sc.hasNextLine()) {
                    String line = sc.nextLine().trim();
                    if (line.isEmpty())
                        continue;
                    String[] parts = line.split(",");
                    if (parts.length < 2)
                        continue;
                    try {
                        electricLoad.add(Double.parseDouble(parts[1].trim()));
                    } catch (NumberFormatException ignored) {
                    }
                }
                System.out.println("DATASET Rows loaded: " + electricLoad.size());
            } catch (FileNotFoundException e) {
                System.out.println("DATASET Read error using synthetic data.");
            }
        } else {
            System.out.println("DATASET CSV not found using synthetic 10-sample dataset.");
        }

        int totalRows = 0;
        // Total fitness cases = 201604 - n
        // int totalRows = electricLoad.size() - n;
        if (numForDataset == 0) {
            totalRows = 101604 - n;

        } else if (numForDataset == 1) {
            totalRows = 10000 - n;
        } else {

            totalRows = electricLoad.size() - n;
        }
        double[][] X = new double[totalRows][n];
        double[] y = new double[totalRows];

        for (int i = 0; i < totalRows; i++) {
            // Fill load columns
            for (int load = 0; load < n; load++) {

                X[i][load] = electricLoad.get(i + load);
            }
            // value right after the loads
            y[i] = electricLoad.get(i + n);
        }
        // Train/Test split (80/20)
        int splitIndex = (int) (totalRows * 0.8);

        // Training data
        X_train = new double[splitIndex][n];
        y_train = new double[splitIndex];

        // Test data
        X_test = new double[totalRows - splitIndex][n];
        y_test = new double[totalRows - splitIndex];

        for (int i = 0; i < splitIndex; i++) {
            X_train[i] = X[i];
            y_train[i] = y[i];
        }

        for (int i = splitIndex; i < totalRows; i++) {
            X_test[i - splitIndex] = X[i];
            y_test[i - splitIndex] = y[i];
        }

    }

    // POPULATION INITIALISATION
    // Creates POP_SIZE individuals, each with a random codon array
    static List<GEIndividual> initPopulation() {

        List<GEIndividual> pop = new ArrayList<>();
        for (int i = 0; i < POP_SIZE; i++) {
            int[] codons = new int[CODON_LENGTH];
            for (int j = 0; j < CODON_LENGTH; j++)
                codons[j] = rng.nextInt(CODON_MAX + 1);
            GEIndividual ind = new GEIndividual(codons);

            pop.add(ind);
        }
        return pop;
    }

    // GENOTYPE to PHENOTYPE MAPPING
    // Maps every individual's codon array to a phenotype string using the BNF
    // grammar, then evaluates fitness.
    static void mapAndEvaluateAll(List<GEIndividual> pop) {
        for (int i = 0; i < pop.size(); i++) {

            BNFGrammar.mapGenotype(pop.get(i));
            evaluateFitness(pop.get(i));

        }
    }

    // FITNESS EVALUATION
    // Evaluates one individual's fitness, MSE on training data.
    static void evaluateFitness(GEIndividual ind) {
        if (!ind.valid || ind.phenotype == null) {
            ind.fitness = 1e9;
            return;
        }
        ind.fitness = computeMSE(ind.phenotype, y_train, X_train);
    }

    // MSE = (1/n) * sum( (predicted_i - actual_i)^2 )
    static double computeMSE(String phenotype, double[] targets, double[][] inputs) {
        double sumSq = 0.0;
        for (int i = 0; i < targets.length; i++) {
            double pred = evaluator.evaluate(phenotype, inputs[i]);
            double error = pred - targets[i];
            sumSq += error * error;
        }
        return sumSq / targets.length;
    }

    // PHENOTYPE SIMILARITY
    // Computes the phenotype similarity.
    static int[] computeSimilarityMatrix(List<GEIndividual> pop) {
        int n = pop.size();
        int[] simCounts = new int[n];

        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                double sim = phenotypeSimilarity(pop.get(i), pop.get(j));
                if (sim >= SIM_THRESHOLD) {
                    simCounts[i]++;
                    simCounts[j]++;
                }
            }
        }

        return simCounts;
    }

    // Character-level phenotype similarity between two individuals.
    static double phenotypeSimilarity(GEIndividual a, GEIndividual b) {
        if (!a.valid || !b.valid)
            return 0.0;
        String sA = a.phenotype;
        String sB = b.phenotype;
        int total = Math.max(sA.length(), sB.length());
        if (total == 0)
            return 1.0;
        int matches = 0;
        int compare = Math.min(sA.length(), sB.length());
        for (int k = 0; k < compare; k++)
            if (sA.charAt(k) == sB.charAt(k))
                matches++;
        return (double) matches / total;
    }

    // Average similarity between one individual and a group.
    static double avgSimilarityToGroup(GEIndividual ind, List<GEIndividual> group) {
        if (group.isEmpty())
            return 0.0;
        double total = 0.0;
        for (GEIndividual m : group)
            total += phenotypeSimilarity(ind, m);
        return total / group.size();
    }

    // Computes the combined selection score for one candidate
    // Score = ALPHA * StructureDiversity + (1 - ALPHA) * FitnessScore
    static double combinedScore(GEIndividual candidate, List<GEIndividual> newPop) {
        double structDiv = 1.0 - avgSimilarityToGroup(candidate, newPop);
        double fitnessScore = 1.0 / (1.0 + candidate.fitness);
        double score = ALPHA * structDiv + (1.0 - ALPHA) * fitnessScore;

        return score;
    }

    // Combined-score tournament selection.
    static GEIndividual selectCombined(List<GEIndividual> pop,
            List<GEIndividual> newPop) {

        int actualK = Math.min(TOURNAMENT_SIZE, pop.size());

        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < pop.size(); i++)
            indices.add(i);
        Collections.shuffle(indices, rng);

        List<Integer> tournament = new ArrayList<>();
        for (int k = 0; k < actualK; k++)
            tournament.add(indices.get(k));

        int bestIdx = tournament.get(0);
        double bestScore = -1.0;

        for (int idx : tournament) {
            double score = combinedScore(pop.get(idx), newPop);
            if (score > bestScore) {
                bestScore = score;
                bestIdx = idx;
            }
        }

        return pop.get(bestIdx);
    }

    // CROSSOVER, one-point on the INTEGER CODON ARRAY
    // Single-point crossover on the codon chromosome.
    static GEIndividual[] crossover(GEIndividual parentA, GEIndividual parentB) {
        int cut = 1 + rng.nextInt(CODON_LENGTH - 1);

        int[] codonsA = new int[CODON_LENGTH];
        int[] codonsB = new int[CODON_LENGTH];

        for (int i = 0; i < CODON_LENGTH; i++) {
            if (i <= cut) {
                codonsA[i] = parentA.codons[i];
                codonsB[i] = parentB.codons[i];
            } else {
                codonsA[i] = parentB.codons[i];
                codonsB[i] = parentA.codons[i];
            }
        }

        GEIndividual childA = new GEIndividual(codonsA);
        GEIndividual childB = new GEIndividual(codonsB);

        return new GEIndividual[] { childA, childB };
    }

    // MUTATION, per-codon random replacement on the INTEGER CODON ARRAY
    // Per-codon mutation on the codon chromosome.
    static void mutateCodons(GEIndividual ind, double mutProb) {
        for (int i = 0; i < CODON_LENGTH; i++) {
            if (rng.nextDouble() < mutProb) {
                ind.codons[i] = rng.nextInt(CODON_MAX + 1);
            }
        }

    }

    // COMPUTATIONAL EFFORT
    static void computeEffort(List<Integer> succGen, int popSize) {
        int targetGen = MAX_GEN - 1;
        int successes = 0;
        for (int g : succGen)
            if (g >= 0 && g <= targetGen)
                successes++;
        double p = (double) successes / succGen.size();

        System.out.println("\n" + sep('=', 70));
        System.out.println("  COMPUTATIONAL EFFORT  [I]");
        System.out.println(sep('=', 70));
        System.out.printf("  Successes by gen %d  : %d / %d%n",
                targetGen, successes, succGen.size());
        System.out.printf("  Success probability p = %.4f%n", p);

        if (p > 0.0) {
            int runsNeeded = (int) Math.ceil(Math.log(0.01) / Math.log(1.0 - p));
            long effort = (long) runsNeeded * popSize * (targetGen + 1);
            System.out.printf("  Runs needed (99%% confidence) : %d%n", runsNeeded);
            System.out.printf("  Computational effort          : %,d evaluations%n", effort);
        } else {
            System.out.println("  No successes recorded - cannot estimate effort.");
        }
    }

    // REPORT
    static void printReport(List<Double> trains, List<Double> tests,
            List<Double> avgs, List<Long> times,
            List<String> eqs) {
        int n = trains.size();
        double sumTrain = 0, sumTest = 0, sumAvg = 0;
        long sumTime = 0;
        for (int i = 0; i < n; i++) {
            sumTrain += trains.get(i);
            sumTest += tests.get(i);
            sumAvg += avgs.get(i);
            sumTime += times.get(i);
        }
        double meanTrain = sumTrain / n;
        double meanTest = sumTest / n;
        double meanAvg = sumAvg / n;

        double varTrain = 0, varTest = 0;
        for (int i = 0; i < n; i++) {
            varTrain += Math.pow(trains.get(i) - meanTrain, 2);
            varTest += Math.pow(tests.get(i) - meanTest, 2);
        }
        double stdTrain = Math.sqrt(varTrain / n);
        double stdTest = Math.sqrt(varTest / n);

        int bestRun = findBestIdx(trains);

        System.out.println("\n" + sep('=', 70));
        System.out.printf("  AGGREGATE RESULTS  (%d runs)%n", n);
        System.out.println(sep('=', 70));
        System.out.printf("  Train MSE   mean=%.6f   std=%.6f%n", meanTrain, stdTrain);
        System.out.printf("  Test  MSE   mean=%.6f   std=%.6f%n", meanTest, stdTest);
        System.out.printf("  Avg best-gen MSE          : %.6f%n", meanAvg);
        System.out.printf("  Best run                  : %d  (train MSE=%.6f)%n",
                bestRun, trains.get(bestRun));
        System.out.println("  Best expression           : " + eqs.get(bestRun));
        System.out.printf("  Total wall time (ms)      : %d%n", endTime - startTime);
        System.out.printf("  Mean run time   (ms)      : %.1f%n", (double) sumTime / n);
    }

    // PREDICTION
    static void printPrediction(GEIndividual bestInd) {
        System.out.println("\n" + sep('=', 70));
        System.out.println("  PREDICTION - overall best model on test set");
        System.out.println(sep('=', 70));
        System.out.println("  Phenotype : " + bestInd.phenotype);
        System.out.println("  Codons (first 20): "
                + Arrays.toString(Arrays.copyOf(bestInd.codons, 20)));
        System.out.printf("  Train MSE : %.6f%n%n",
                computeMSE(bestInd.phenotype, y_train, X_train));

        System.out.printf("  %-6s  %-10s  %-12s  %-10s%n",
                "Row", "Actual", "Predicted", "AbsErr");
        System.out.println("  " + sep('-', 44));

        double sumAbsErr = 0.0;
        for (int i = 0; i < y_test.length; i++) {
            double pred = evaluator.evaluate(bestInd.phenotype, X_test[i]);
            double absErr = Math.abs(pred - y_test[i]);
            sumAbsErr += absErr;
        }
        System.out.println("  " + sep('-', 44));
        System.out.printf("  Mean Absolute Error : %.4f%n", sumAbsErr / y_test.length);
        System.out.printf("  Test MSE            : %.6f%n",
                computeMSE(bestInd.phenotype, y_test, X_test));
    }

    // POPULATION UTILITIES
    static GEIndividual getBest(List<GEIndividual> pop) {
        GEIndividual best = pop.get(0);
        for (GEIndividual ind : pop)
            if (ind.fitness < best.fitness)
                best = ind;
        return best;
    }

    static double getBestFitness(List<GEIndividual> pop) {
        return getBest(pop).fitness;
    }

    static double getAvgFitness(List<GEIndividual> pop) {
        double sum = 0;
        int cnt = 0;
        for (GEIndividual ind : pop) {
            if (ind.fitness < 1e8) {
                sum += ind.fitness;
                cnt++;
            }
        }
        return cnt > 0 ? sum / cnt : 0.0;
    }

    static int countValid(List<GEIndividual> pop) {
        int cnt = 0;
        for (GEIndividual ind : pop)
            if (ind.valid)
                cnt++;
        return cnt;
    }

    // GENERAL UTILITIES
    static int findBestIdx(List<Double> vals) {
        int best = 0;
        for (int i = 1; i < vals.size(); i++)
            if (vals.get(i) < vals.get(best))
                best = i;
        return best;
    }

    static String sep(char c, int n) {
        char[] arr = new char[n];
        Arrays.fill(arr, c);
        return new String(arr);
    }

    static void tick() {
        startTime = System.currentTimeMillis();
    }

    static void tock() {
        endTime = System.currentTimeMillis();
    }

    static void printBanner() {

        System.out.println("  Structure-Based GRAMMATICAL EVOLUTION");
        System.out.println("  Assignment 3 - Electricity Load Prediction");
        System.out.println(sep('-', 70));
    }

    static void printConfig() {
        System.out.println("\n" + sep('-', 70));
        System.out.println("  GE HYPERPARAMETERS");
        System.out.println(sep('-', 70));
        System.out.printf("  Codon length       : %d integers%n", CODON_LENGTH);
        System.out.printf("  Codon range        : [0, %d]%n", CODON_MAX);
        System.out.printf("  Max wraps          : %d%n", MAX_WRAPS);
        System.out.printf("  Mapping max depth  : %d%n", MAPPING_MAX_DEPTH);
        System.out.printf("  Population size    : %d%n", POP_SIZE);
        System.out.printf("  Generations        : %d%n", MAX_GEN);
        System.out.printf("  Tournament size    : %d%n", TOURNAMENT_SIZE);
        // System.out.printf(" Crossover rate : %.2f%n", CROSSOVER_RATE);
        System.out.printf("  Mutation prob/codon: %.4f (LOW=%.4f  HIGH=%.4f)%n",
                MUTATION_PROB, MUTATION_LOW, MUTATION_HIGH);
        System.out.printf("  Sim threshold      : %.2f%n", SIM_THRESHOLD);
        System.out.printf("  Cvg fraction       : %.2f%n", CVG_FRACTION);
        System.out.printf("  Independent runs   : %d%n", NO_OF_RUNS);
        System.out.printf("  Base seed          : %d%n", SEED);
        System.out.printf("  Success threshold  : %.6f%n", SUCCESS_THRESHOLD);
        System.out.printf("  Train samples      : %d%n", y_train.length);
        System.out.printf("  Test  samples      : %d%n", y_test.length);
    }
}

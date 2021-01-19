package util;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class Helper {
    /**
     * getColLength
     * Returns an integer 1.2 times greater than the longest sequence
     * in the given sequences.
     *
     * @param seqs MSA util.Sequences
     * @return column length
     */
    public static int getColLength(String[] seqs) {
        int length = 0;
        for (String seq : seqs) {
            if (seq.length() > length)
                length = seq.length();
        }
        return (int) Math.ceil(length * 1.2);
    }

    /**
     * An aggregated function.
     *
     * @param bitmatrix
     * @param seq
     * @param w1
     * @param w2
     * @param checkInfeasibility
     * @param ops
     * @return
     */
    public static double aggregatedFunction(int[][] bitmatrix,
                                            String[] seq,
                                            double w1,
                                            double w2,
                                            boolean checkInfeasibility,
                                            Operator[] ops) {
        if (checkInfeasibility && ops == null)
            throw new IllegalArgumentException("Checking for infeasibility but ops are null!");

        if (checkInfeasibility && infeasible(bitmatrix, seq, ops))
            return Double.MIN_VALUE;

        String[] strings = bitsToStrings(bitmatrix, seq);

        int nMax = maxNumOfIndels(bitmatrix, seq); // total number of indels

        int nI = numOfInsertedIndels(bitmatrix, seq); // number of indels before last char

        return (w1 * numOfAlignedChars(strings)) + (w2 * (nMax - nI));
    }

    /**
     * Counts the number of indels that are found
     * before the last character in a sequence.
     */
    public static int numOfInsertedIndels(int[][] bitmatrix,
                                          String[] seq) {
        /*
        Remember: bit 0 means character from sequence
                     bit 1 means indel
         */

        int count = 0;

        for (int i = 0; i < seq.length; i++) { // Go through each bitstring in bitmatrix
            int tmp = 0;
            boolean hitLastChar = false;
            for (int bit : bitmatrix[i]) { // Go through each bit in bitstring
                if (bit == 0) { // we hit a character
                    tmp += 1;
                    if (tmp == seq[i].length())
                        hitLastChar = true;
                } else // we hit an indel
                    count += 1;

                if (hitLastChar)
                    break;
            }
        }

        return count;
    }

    /**
     * Counts the total number of indels according to the difference between
     * bitmatrix row length and sequence length.
     * <p>
     * Does not consider individual bits.
     */
    public static int maxNumOfIndels(int[][] bitmatrix,
                                     String[] seq) {
        int count = 0;
        for (int i = 0; i < bitmatrix.length; i++)
            count += bitmatrix[i].length - seq[i].length();
        return count;
    }

    /**
     * Determines if a bitmatrix and sequence combination is infeasible.
     * <p>
     * Counts the number of 0 bits in each row of the bitmatrix, and if it
     * doesn't match up with an operator comparison between the 0 bit count
     * and the length of the sequence, returns True for infeasibility.
     * <p>
     * The operator can be anything within the parameter format, however
     * usually it is >, <, =. For easy access to these use the operator
     * module.
     * <p>
     * Ex. Say ops is [>, <], then bitmatrix is infeasible when # of 0's > # of chars,
     * and when # of 0's < # of chars. Thus it will only be feasible when # of 0's = # of chars.
     */
    public static boolean infeasible(int[][] bitmatrix,
                                     String[] seq,
                                     Operator[] ops) {
        for (int i = 0; i < seq.length; i++) {
            int count = 0;
            for (int bit : bitmatrix[i])
                if (bit == 0)
                    count += 1;

            for (Operator op : ops)
                if (op.operation(count, seq[i].length()))
                    return true;
        }

        return false;
    }

    /**
     * Converts a list of sequences into a list of strings with indels, according to the bitmatrix provided.
     * <p>
     * Note: A bit of 0 means a character, and a bit of 1 means indel.
     * However, if the number of 0 bits exceeds the length of the sequence,
     * indels will be inserted instead.
     */
    public static String[] bitsToStrings(int[][] bitmatrix,
                                         String[] seq) {
        String[] result = new String[seq.length];
        int i = 0;
        for (int[] bitlist : bitmatrix) {
            int j = 0;
            StringBuilder tmp = new StringBuilder();
            for (int bit : bitlist) {
                if (bit == 0 && j < seq[i].length()) {
                    tmp.append(seq[i].charAt(j));
                    j++;
                } else
                    tmp.append("-");
            }
            result[i] = tmp.toString();
            i++;
        }
        return result;
    }

    /**
     * Counts the number of aligned characters in a list of strings.
     * <p>
     * Important: Assumes that each string is of the same length.
     */
    public static int numOfAlignedChars(String[] strings) {
        if (strings.length < 1)
            throw new IllegalArgumentException("List of strings is empty");
        else if (strings.length == 1) {
            System.out.println("Warning: only 1 string in numOfAlignedChars function!");
            return 0;
        }

        int result = 0;
        HashMap<Character, Integer> charList = new HashMap<>();

        /*
            The charList is a mapping, of each character to
            the amount of matchings that character has.
         */
        for (int i = 0; i < strings[0].length(); i++) {
            HashMap<Character, Integer> tmpMapping = new HashMap<>(); // Temp before adding to main charList
            // Fill out the tmpMapping for all characters in column i
            for (String str : strings) {
                Character c = str.charAt(i);
                if (c != '-') {
                    if (tmpMapping.containsKey(c))
                        tmpMapping.put(c, tmpMapping.get(c) + 1);
                    else
                        tmpMapping.put(c, 1);
                }
            }

            // Now add the tmpMapping into main charList
            for (Character c : tmpMapping.keySet()) {
                // if there's more than one of that character in the col, there's a matching
                if (tmpMapping.get(c) > 1) {
                    if (charList.containsKey(c))
                        charList.put(c, charList.get(c) + tmpMapping.get(c));
                    else
                        charList.put(c, tmpMapping.get(c));
                }
            }
        }

        for (Integer i : charList.values())
            result += i;

        return result;
    }

    /**
     * The classic sigmoid function.
     */
    public static double Sigmoid(double x) {
        return 1 / (1 + Math.exp(-x));
    }

    /**
     * Generating a bitMatrix from the
     * angular modulated generation function.
     * <p>
     * Values are sampled from 0, 1, ..., n-1
     *
     * @param pos
     * @param seq
     * @param colLength
     * @return
     */
    public static int[][] genBitMatrix(double[] pos,
                                       String[] seq,
                                       int colLength) {
        int[][] bitmatrix = new int[seq.length][colLength];

        for (int i = 0; i < seq.length; i++) {
            for (int j = 0; j < colLength; j++) {
                double value = gen((i * colLength) + j, pos[0], pos[1], pos[2], pos[3]);
                if (value > 0)
                    bitmatrix[i][j] = 1;
                else
                    bitmatrix[i][j] = 0;
            }
        }

        return bitmatrix;
    }

    /**
     * Angular Modulated Generation Function
     *
     * @param x
     * @param a
     * @param b
     * @param c
     * @param d
     * @return float
     */
    private static double gen(double x,
                              double a,
                              double b,
                              double c,
                              double d) {
        return Math.sin(2 * Math.PI * (x - a) * b * Math.cos(2 * Math.PI * c * (x - a))) + d;
    }

    /**
     * Checks if one bitmatrix dominates another.
     * <p>
     * This means that for each FitnessFunction,
     * a bitmatrix must be at least as good as
     * the other in all functions, and in at least
     * one function it must be better.
     *
     * @param seq
     * @param bm1
     * @param bm2
     * @return
     */
    private static boolean dominates(String[] seq,
                                     int[][] bm1,
                                     int[][] bm2,
                                     FitnessFunction[] f) {
        boolean betterInAtLeastOne = false;

        for (FitnessFunction ff : f) {
            double res1 = ff.calculate(bm1, seq);
            double res2 = ff.calculate(bm2, seq);

            if (res1 < res2)
                return false;
            else if (res1 > res2)
                betterInAtLeastOne = true;
        }

        return betterInAtLeastOne;
    }

    /**
     * Checks that two bitmatrixes are the same.
     *
     * @param x
     * @param y
     * @return
     */
    private static boolean theSame(int[][] x,
                                   int[][] y) {
        for (int i = 0; i < x.length; i++) {
            if (x[i].length != y[i].length)
                return false;

            for (int j = 0; j < x[i].length; j++) {
                if (x[i][j] != y[i][j])
                    return false;
            }
        }

        return true;
    }

    /**
     * Solution x is attempted to be added into the
     * archive if it is not dominated.
     * <p>
     * If after addition to the archive, and previous
     * solutions are now dominated by this new addition,
     * then those dominated solutions are removed.
     * <p>
     * Also, if the archive is full then the most crowded
     * solution is removed.
     *
     * @param seq
     * @param sArchive
     * @param x
     * @param n
     * @return
     */
    public static void addToArchive(String[] seq,
                                    ArrayList<int[][]> sArchive,
                                    int[][] x,
                                    FitnessFunction[] f,
                                    int n) {
        ArrayList<int[][]> aDominated = new ArrayList<>();

        for (int[][] s : sArchive) {
            // archive entry dominates x
            if (dominates(seq, s, x, f))
                return;
            else if (dominates(seq, x, s, f)) // x dominates archive entry
                aDominated.add(s);

            // x and archive entry are the exact same
            if (theSame(x, s))
                return;
        }

        /*
         When function enters here, x is not dominated by any
         entry in the archive, but may have dominated some
         entries. So add it to archive, but remove all the
         previous entries that are dominated by it.

         If the limit of the archive has been reached, then
         calculate crowding distances and remove the most
         crowded solution.
         */

        sArchive.removeIf(aDominated::contains);

        /*
             If this is true then the archive is completely full so
             remove most crowded solution.
         */
        if (sArchive.size() == n) {
            // Updated crowding distance: collect unique fitnesses first and perform calculations on these
            ArrayList<Double> uniqueFitnesses = getUniqueFitnesses(seq, f, sArchive);
            double mostCrowdedFitness = getCrowdingDistances(uniqueFitnesses, f).get(0).getFirst();

            /*
                Distances arraylist will be sorted in asc. order by distance

                However, we must match the fitness to an archive entry in
                order to remove it, so this is done here, and after the
                most crowded fitness is removed, x is added and we exit
                this function.
             */
            for (FitnessFunction ff : f) {
                for (int j = 0; j < sArchive.size(); j++) {
                    if (ff.calculate(sArchive.get(j), seq) == mostCrowdedFitness) {
                        sArchive.remove(j);
                        sArchive.add(x);
                        return;
                    }
                }
            }
        }

        sArchive.add(x);
    }

    private static ArrayList<Double> getUniqueFitnesses(String[] seq,
                                                        FitnessFunction[] f,
                                                        ArrayList<int[][]> sArchive) {
        ArrayList<Double> uniqueFitnesses = new ArrayList<>();

        for (FitnessFunction ff : f) {
            for (int[][] s : sArchive) {
                double fitness = ff.calculate(s, seq);
                if (!uniqueFitnesses.contains(fitness))
                    uniqueFitnesses.add(fitness);
            }
        }

        return uniqueFitnesses;
    }

    /**
     * Get a guide from the archive.
     *
     * @param seq      sequences
     * @param sArchive the archive
     * @param f        list of functions
     * @param k        number of elements for tournament selection
     * @return archive value as guide
     */
    public static int[][] archiveGuide(String[] seq,
                                       ArrayList<int[][]> sArchive,
                                       FitnessFunction[] f,
                                       int k) {
        ArrayList<Pair<Double, Double>> distances = getCrowdingDistances(getUniqueFitnesses(seq, f, sArchive), f);

        int index = ThreadLocalRandom.current().nextInt(0, distances.size());

        for (int i = 0; i < k; i++) {
            int tmp = ThreadLocalRandom.current().nextInt(0, distances.size());
            /*
                This works because in the getCrowdingDistances
                function, all elements are sorted in asc. order
                by crowding distance value.
             */
            if (index > tmp)
                index = tmp;
        }

        double fitness = distances.get(index).getFirst();

        for (FitnessFunction ff : f)
            for (int[][] s : sArchive)
                if (ff.calculate(s, seq) == fitness)
                    return s;

        throw new RuntimeException("archive guide empty");
    }

    /**
     * Crowding distance calculation for archive
     * entry removal.
     * <p>
     * Uses the updated crowding distance calculation from:
     * Fortin, F.-A., & Parizeau, M. (2013). Revisiting the NSGA-II crowding-distance computation. Proceeding of the Fifteenth Annual Conference on Genetic and Evolutionary Computation Conference - GECCO  ’13. doi:10.1145/2463372.2463456
     *
     * @param uniqueFitnesses
     * @param f
     */
    private static ArrayList<Pair<Double, Double>> getCrowdingDistances(ArrayList<Double> uniqueFitnesses,
                                                                        FitnessFunction[] f) {
        /*
            TODO: double check this function works properly with the newly
            suggested crowding distance on unique fitness values.
         */
        /*
            Unique fitnesses are used for distance calculations,
            but the structure must be updated to contain a pair
            so we can match a fitness with its corresponding distance.
         */
        ArrayList<Pair<Double, Double>> uniqueDistancedFitnesses = new ArrayList<>();
        for (double d : uniqueFitnesses) {
            uniqueDistancedFitnesses.add(new Pair<>(d, 0.0));
        }

        // Now the uniqueDistancedFitnesses structure is set up and ready to be used

        for (FitnessFunction ff : f) {
            // Sort uniqueDistancedFitnesses in ascending order
            Comparator<Pair<Double, Double>> c = (o1, o2) -> o1.getFirst() > o2.getFirst() ? 1 : 0;
            uniqueDistancedFitnesses.sort(c);

            // Set the boundary distances to infinite (max value)
            uniqueDistancedFitnesses.get(0).setSecond(Double.MAX_VALUE);
            uniqueDistancedFitnesses.get(uniqueDistancedFitnesses.size() - 1).setSecond(Double.MAX_VALUE);

            for (int i = 1; i < uniqueDistancedFitnesses.size() - 1; i++) {
                double curr = uniqueDistancedFitnesses.get(i).getSecond();
                double next = uniqueDistancedFitnesses.get(i + 1).getFirst();
                double prev = uniqueDistancedFitnesses.get(i - 1).getFirst();
                uniqueDistancedFitnesses.get(i).setSecond(curr + next - prev);
            }
        }

        // This will sort uniqueDistancedFitnesses based on distance (ascending)
        Comparator<Pair<Double, Double>> c = (o1, o2) -> o1.getSecond() > o2.getSecond() ? 1 : 0;
        uniqueDistancedFitnesses.sort(c);

        return uniqueDistancedFitnesses;
    }
}

import java.util.*;

public class ThresholdSweepAnalyzer {

    public static class SweepResult {
        public double varH, varD, varG, varC, varR, varP;
        @Override
        public String toString() {
            return String.format(
                    "Var(H)=%.6e, Var(D)=%.6e, Var(G)=%.6e, Var(C)=%.6e, Var(R)=%.6e, Var(P)=%.6e",
                    varH, varD, varG, varC, varR, varP
            );
        }
    }

    /** Main method: sweep thresholds and calculate the variance */
    public static SweepResult sweep(SandboxCity city, double resolution) {
        List<Double> Hs = new ArrayList<>();
        List<Double> Ds = new ArrayList<>();
        List<Double> Gs = new ArrayList<>();
        List<Double> Cs = new ArrayList<>();
        List<Double> Rs = new ArrayList<>();
        List<Double> Ps = new ArrayList<>();

        List<double[]> combos = generateThresholdCombos(resolution, city.transit.levels.length);

        // Saved original transitValue
        double[][] originalTransitValues = new double[city.getSize()][city.getSize()];
        for (int i = 0; i < city.getSize(); i++) {
            for (int j = 0; j < city.getSize(); j++) {
                originalTransitValues[i][j] = city.grid[i][j].transitValue;
            }
        }

        for (double[] th : combos) {
            // Using threshold to generate new TransitLevels
            TransitLevels tl = fromThresholds(th);

            // Updating level
            for (int i = 0; i < city.getSize(); i++) {
                for (int j = 0; j < city.getSize(); j++) {
                    double val = originalTransitValues[i][j];
                    city.grid[i][j].level = tl.indexOf(val);
                }
            }

            // Calculating index
            double H = city.computeEntropyRatio();
            double D = city.computeDissimilarity();
            double G = city.computeGini();
            double C = city.computeVariation();
            double R = city.computeRelativeDiversity();
            double P = city.computeExposure();

            Hs.add(H); Ds.add(D); Gs.add(G); Cs.add(C); Rs.add(R); Ps.add(P);
        }

        SweepResult result = new SweepResult();
        result.varH = variance(Hs);
        result.varD = variance(Ds);
        result.varG = variance(Gs);
        result.varC = variance(Cs);
        result.varR = variance(Rs);
        result.varP = variance(Ps);

        return result;
    }

    /** Using resolution to generate all possible combination of decreasing thresholds*/
    private static List<double[]> generateThresholdCombos(double resolution, int levels) {
        List<Double> values = new ArrayList<>();
        for (double v = resolution; v <= 1.0; v += resolution) values.add(v);
        Collections.reverse(values); // reverse to order

        List<double[]> combos = new ArrayList<>();
        int needed = levels - 1; // n cats need n-1 thresholds
        generateCombinations(values, needed, new double[needed], 0, combos);
        return combos;
    }

    private static void generateCombinations(List<Double> values, int needed, double[] current, int idx, List<double[]> out) {
        if (idx == needed) {
            // Checking for monotonically decreasing
            for (int i = 1; i < current.length; i++) {
                if (!(current[i - 1] > current[i])) return;
            }
            out.add(Arrays.copyOf(current, current.length));
            return;
        }
        for (double v : values) {
            current[idx] = v;
            generateCombinations(values, needed, current, idx + 1, out);
        }
    }

    /** Using Thresholds to contruct TransitLevels */
    private static TransitLevels fromThresholds(double[] th) {
        int n = th.length + 1;
        TransitLevels.Level[] levels = new TransitLevels.Level[n];

        double prev = Double.POSITIVE_INFINITY;
        for (int idx = 0; idx < n; idx++) {
            double cur = (idx == n - 1) ? Double.NEGATIVE_INFINITY : th[idx];
            String name = "T" + (idx + 1);
            levels[idx] = new TransitLevels.Level(name, cur, prev);
            prev = cur;
        }
        return new TransitLevels(levels);
    }

    /** Calculating Variance */
    private static double variance(List<Double> values) {
        double[] arr = values.stream()
                .filter(d -> d != null && !d.isNaN() && !d.isInfinite())
                .mapToDouble(Double::doubleValue).toArray();
        if (arr.length == 0) return Double.NaN;
        double mean = Arrays.stream(arr).average().orElse(Double.NaN);
        return Arrays.stream(arr).map(v -> (v - mean) * (v - mean)).average().orElse(Double.NaN);
    }

}

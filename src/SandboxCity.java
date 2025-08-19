import java.util.Random;

public class SandboxCity {


    /* ---------- Public Generation Mode---------- */
    public enum GenerationMode {
        HIERARCHICAL,       // π ~ Dir(β)，每格 p ~ Dir(κ·π)
        DIRICHLET1,         // 每格 p ~ Dir(1,...,1)
        UNIFORM_NORMALIZE   // 旧版：先取 U(0,1) 后整体归一
    }


    /* ---------- Basic Status ---------- */
    private final int size;
    public final Neighborhood[][] grid;
    public final String[] races;
    public final String[] transitCategories;
    public final TransitLevels transit;

    /* ---------- Random and Parameters ---------- */
    private double beta = 0.3;         // 抽 π 的浓度（越小越偏）
    private double kappa = 0.6;        // 单格浓度（越小越尖）
    private final Random rand;
    private long seed;                 // 当前种子（可选）
    private GenerationMode mode = GenerationMode.HIERARCHICAL;

    /* ---------- Debugger ---------- */
    private static final boolean DEBUG = false;
    private double[] lastPi = null;    // 最近一次 HIERARCHICAL/DIRICHLET1 时的 π（UNIFORM_NORMALIZE 下为 null）


//    public SandboxCity(int size, int raceCount, int transitCount) {
//        this.size = size;
//        this.grid = new Neighborhood[size][size];
//
//        races = new String[raceCount];
//        for (int i = 0; i < raceCount; i++) {
//            races[i] = "Race" + (i+1);
//        }
//        transitCategories = new String[transitCount];
//        for (int j = 0; j < transitCount; j++) {
//            transitCategories[j] = "T" + (j+1);
//        }
//
//        transit = new TransitLevels(transitCount, rand);
//    }
    /* ---------- Constructor ---------- */
    public SandboxCity(int size, int raceCount, int transitCount) {
        this(size, raceCount, transitCount, 0.3, 0.6, System.nanoTime());
    }

    public SandboxCity(int size, int raceCount, int transitCount, double beta, double kappa, long seed) {
        this.size = size;
        this.grid = new Neighborhood[size][size];
        this.beta = Math.max(1e-6, beta);
        this.kappa = Math.max(1e-6, kappa);
        this.seed = seed;

        this.rand = new Random(seed);

        races = new String[raceCount];
        for (int i = 0; i < raceCount; i++) races[i] = "Race" + (i + 1);

        transitCategories = new String[transitCount];
        for (int j = 0; j < transitCount; j++) transitCategories[j] = "T" + (j + 1);

        transit = new TransitLevels(transitCount, rand);
    }

    /* ---------- Accessors ---------- */
    public int getSize() { return size; }
    public double getBeta()  { return beta;  }
    public double getKappa() { return kappa; }
    public void setParams(double beta, double kappa) {
        this.beta  = Math.max(1e-6, beta);
        this.kappa = Math.max(1e-6, kappa);
    }
    public GenerationMode getMode() { return mode; }
    public void setMode(GenerationMode m) { this.mode = (m == null ? GenerationMode.HIERARCHICAL : m); }
    public long getSeed() { return seed; }
    public void reseed(long newSeed) { this.seed = newSeed; rand.setSeed(newSeed); }
    public double[] getLastPi() {
        return (lastPi == null) ? null : java.util.Arrays.copyOf(lastPi, lastPi.length);
    }

    /* ---------- 采样与初始化 ---------- */
    public void initializeRandom() {
        int K = races.length;
        double[] pi;
        double kappaLocal;

        switch (mode) {
            case DIRICHLET1:
                // π 均匀，κ=K => α_i = 1
                pi = new double[K];
                java.util.Arrays.fill(pi, 1.0 / K);
                kappaLocal = K;
                this.lastPi = java.util.Arrays.copyOf(pi, K);
                break;
            case UNIFORM_NORMALIZE:
                // 本模式不使用 π；这里把 lastPi 置空仅用于 snapshot 展示
                pi = null;
                kappaLocal = Double.NaN; // 无意义
                this.lastPi = null;
                break;
            default: // HIERARCHICAL
                pi = dirichletSymmetric(K, beta, rand); // 全局配比
                kappaLocal = kappa;
                this.lastPi = java.util.Arrays.copyOf(pi, K);
                break;
        }

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                double[] rp;
                if (mode == GenerationMode.UNIFORM_NORMALIZE) {
                    rp = randomProbabilitiesUniform(K); // 旧版：归一化均匀数
                } else {
                    double[] alphaCell = new double[K];
                    for (int k = 0; k < K; k++)
                        alphaCell[k] = Math.max(pi[k] * kappaLocal, 1e-8); // 数值下限更稳
                    rp = dirichlet(alphaCell, rand);
                }

                double tVal = rand.nextDouble();
                int lvl = transit.indexOf(tVal);
                grid[i][j] = new Neighborhood(rp, lvl, tVal);
            }
        }
        if (DEBUG) debugSnapshot(lastPi);
    }

    /* ---------- 快照（给 GUI/控制台） ---------- */
    public String snapshotString() {
        int K = races.length, m = size * size;
        double[] avg = new double[K];
        double maxCellShare = 0, minCellShare = 1;

        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) {
                double[] rp = grid[i][j].raceProps;
                double mx = 0;
                for (int k = 0; k < K; k++) {
                    avg[k] += rp[k];
                    if (rp[k] > mx) mx = rp[k];
                }
                if (mx > maxCellShare) maxCellShare = mx;
                if (mx < minCellShare) minCellShare = mx;
            }
        for (int k = 0; k < K; k++) avg[k] /= m;

        StringBuilder sb = new StringBuilder();
        sb.append("mode=").append(mode)
                .append(", β=").append(beta)
                .append(", κ=").append(kappa)
                .append(", size=").append(size).append(" (m=").append(m).append(")\n");
        sb.append("π     = ")
                .append(lastPi == null ? "N/A (UNIFORM_NORMALIZE)" : java.util.Arrays.toString(lastPi))
                .append('\n');
        sb.append("avg   = ").append(java.util.Arrays.toString(avg)).append('\n');
        sb.append(String.format("cell-max share range = [%.3f, %.3f]%n", minCellShare, maxCellShare));
        return sb.toString();
    }
    public void printSnapshot() { System.out.print(snapshotString()); }


    /* ---------- Dirichlet / Gamma / legacy uniform  ---------- */
    private double gamma(double a, Random r){
        if (a < 1.0) return gamma(a + 1.0, r) * Math.pow(r.nextDouble(), 1.0/a);
        double d = a - 1.0/3.0, c = 1.0/Math.sqrt(9.0*d);
        while(true){
            double x = r.nextGaussian();
            double v = 1 + c*x; if (v <= 0) continue; v = v*v*v;
            double u = r.nextDouble();
            if (u < 1 - 0.0331*x*x*x*x) return d*v;
            if (Math.log(u) < 0.5*x*x + d*(1 - v + Math.log(v))) return d*v;
        }
    }
    private double[] dirichlet(double[] alpha, Random r){
        int K = alpha.length; double[] x = new double[K]; double s=0;
        for(int i=0;i<K;i++){ x[i] = gamma(alpha[i], r); s += x[i]; }
        for(int i=0;i<K;i++) x[i] /= s; return x;
    }
    private double[] dirichletSymmetric(int K, double a, Random r){
        double[] alpha = new double[K]; java.util.Arrays.fill(alpha, a);
        return dirichlet(alpha, r);
    }
    private double[] randomProbabilitiesUniform(int n){
        double[] a = new double[n];
        double s = 0;
        for (int i=0;i<n;i++){ a[i] = rand.nextDouble(); s += a[i]; }
        for (int i=0;i<n;i++) a[i] /= s;
        return a;
    }

    /* ---------- 可选：控制台 debug ---------- */
    private void debugSnapshot(double[] pi){
        int K = races.length, m = size*size;
        double[] avg = new double[K];
        double maxCellShare = 0, minCellShare = 1;

        for (int i=0;i<size;i++)
            for (int j=0;j<size;j++){
                double[] rp = grid[i][j].raceProps;
                double mx = 0;
                for (int k=0;k<K;k++) {
                    avg[k] += rp[k];
                    if (rp[k] > mx) mx = rp[k];
                }
                if (mx>maxCellShare) maxCellShare=mx;
                if (mx<minCellShare) minCellShare=mx;
            }
        for (int k=0;k<K;k++) avg[k] /= m;

        System.out.println("mode=" + mode + ", β=" + beta + ", κ=" + kappa);
        System.out.println("π     = " + (pi==null? "N/A" : java.util.Arrays.toString(pi)));
        System.out.println("avg   = " + java.util.Arrays.toString(avg));
        System.out.printf ("cell-max share range = [%.3f, %.3f]%n", minCellShare, maxCellShare);
    }



//    private double gamma(double a, Random r){
//        if (a < 1.0) return gamma(a + 1.0, r) * Math.pow(r.nextDouble(), 1.0/a);
//        double d = a - 1.0/3.0, c = 1.0/Math.sqrt(9.0*d);
//        while(true){
//            double x = r.nextGaussian();
//            double v = 1 + c*x; if (v <= 0) continue; v = v*v*v;
//            double u = r.nextDouble();
//            if (u < 1 - 0.0331*x*x*x*x) return d*v;
//            if (Math.log(u) < 0.5*x*x + d*(1 - v + Math.log(v))) return d*v;
//        }
//    }
//    private double[] dirichlet(double[] alpha, Random r){
//        int K = alpha.length; double[] x = new double[K]; double s=0;
//        for(int i=0;i<K;i++){ x[i] = gamma(alpha[i], r); s += x[i]; }
//        for(int i=0;i<K;i++) x[i] /= s; return x;
//    }
//    private double[] dirichletSymmetric(int K, double a, Random r){
//        double[] alpha = new double[K]; java.util.Arrays.fill(alpha, a);
//        return dirichlet(alpha, r);
//    }

//    private static final boolean DEBUG = true;

//    private void debugSnapshot(double[] pi){
//        int K = races.length, m = size*size;
//        double[] avg = new double[K];
//        double maxCellShare = 0, minCellShare = 1;
//
//        for (int i=0;i<size;i++)
//            for (int j=0;j<size;j++){
//                double[] rp = grid[i][j].raceProps;
//                for (int k=0;k<K;k++) avg[k] += rp[k];
//                double cellMax = 0;
//                for (double v: rp) if (v>cellMax) cellMax = v;
//                if (cellMax>maxCellShare) maxCellShare=cellMax;
//                if (cellMax<minCellShare) minCellShare=cellMax;
//            }
//        for (int k=0;k<K;k++) avg[k] /= m;
//
//        System.out.println("π     = " + java.util.Arrays.toString(pi));
//        System.out.println("avg   = " + java.util.Arrays.toString(avg));
//        System.out.printf ("cell-max share range = [%.3f, %.3f]%n", minCellShare, maxCellShare);
//    }


//    public void initializeRandom() {
//        int K = races.length;
//        double[] pi = dirichletSymmetric(K, 0.3, rand); // 全局配比，可能很偏
//        double kappa = 0.6;                               // 浓度，小=每格更散
//
//        for(int i=0;i<size;i++){
//            for(int j=0;j<size;j++){
//                double[] alphaCell = new double[K];
//                for(int k=0;k<K;k++) alphaCell[k] = pi[k]*kappa;
//                double[] rp = dirichlet(alphaCell, rand);     // 围着 π 抽
//                double tVal = rand.nextDouble();
//                int lvl = transit.indexOf(tVal);
//                grid[i][j] = new Neighborhood(rp, lvl, tVal);
//            }
//        }
//
//        if (DEBUG) debugSnapshot(pi);
//
//
////        Random rand = new Random();
////        for(int i=0;i<size;i++){
////            for(int j=0;j<size;j++){
////                double[] rp = randomProbabilities(races.length);
////                double tVal = rand.nextDouble();
////                int lvl = transit.indexOf(tVal);
////                grid[i][j] = new Neighborhood(rp, lvl, tVal);
////            }
////        }
//    }



//    public void initializeRandom() {
//        int K = races.length;
//        double[] pi = dirichletSymmetric(K, beta, rand); // 全局配比
//        this.lastPi = java.util.Arrays.copyOf(pi, K);    // 记录下来做 snapshot
//
//        for (int i = 0; i < size; i++) {
//            for (int j = 0; j < size; j++) {
//                double[] alphaCell = new double[K];
//                for (int k = 0; k < K; k++)
//                    alphaCell[k] = Math.max(pi[k] * kappa, 1e-8); // 数值下限更稳
//                double[] rp = dirichlet(alphaCell, rand);         // 围着 π 抽
//                double tVal = rand.nextDouble();
//                int lvl = transit.indexOf(tVal);
//                grid[i][j] = new Neighborhood(rp, lvl, tVal);
//            }
//        }
//        if (DEBUG) debugSnapshot(pi);
//    }

//    public String snapshotString() {
//        int K = races.length, m = size * size;
//        double[] avg = new double[K];
//        double maxCellShare = 0, minCellShare = 1;
//
//        for (int i = 0; i < size; i++)
//            for (int j = 0; j < size; j++) {
//                double[] rp = grid[i][j].raceProps;
//                double mx = 0;
//                for (int k = 0; k < K; k++) {
//                    avg[k] += rp[k];
//                    if (rp[k] > mx) mx = rp[k];
//                }
//                if (mx > maxCellShare) maxCellShare = mx;
//                if (mx < minCellShare) minCellShare = mx;
//            }
//        for (int k = 0; k < K; k++) avg[k] /= m;
//
//        StringBuilder sb = new StringBuilder();
//        sb.append("β=").append(beta).append(", κ=").append(kappa)
//                .append(", size=").append(size).append(" (m=").append(m).append(")\n");
//        sb.append("π     = ").append(java.util.Arrays.toString(lastPi)).append('\n');
//        sb.append("avg   = ").append(java.util.Arrays.toString(avg)).append('\n');
//        sb.append(String.format("cell-max share range = [%.3f, %.3f]%n", minCellShare, maxCellShare));
//        return sb.toString();
//    }
//    public void printSnapshot() { System.out.print(snapshotString()); }

    // Dirichlet(1)
//    private double[] randomProbabilities(int n){
//        double[] a = new double[n];
//        double sum = 0;
//        for(int k=0; k<n; k++){
//            double u = rand.nextDouble();
//            a[k] = -Math.log(1 - u);   // Exp(1)；Gamma(1,1)
//            sum += a[k];
//        }
//        for(int k=0; k<n; k++) a[k] /= sum;
//        return a;
//    }

//    private double[] randomProbabilities(int n){
//        double[] a = new double[n];
//        double sum = 0;
//        for(int k=0; k<n; k++){
//            a[k] = rand.nextDouble();
//            sum += a[k];
//        }
//        for(int k=0; k<n; k++) a[k] /= sum;
//        return a;
//    }

    public double computeTotalEntropy(){
        int n = races.length;
        double total = size*size, HX = 0;
        double[] pX = new double[n];
        for(int i=0;i<size;i++)
            for(int j=0;j<size;j++)
                for(int x=0;x<n;x++)
                    pX[x] += grid[i][j].raceProps[x];
        for(int x=0;x<n;x++){
            pX[x] /= total;
            if(pX[x]>0) HX -= pX[x]*log2(pX[x]);
        }
        return HX;
    }

    public double computeConditionalEntropy() {
        int M = transit.levels.length;   // 档数
        int R = races.length;
        double cells = size * size;

        /* ① 聚合到层级：sumY[y][x] = 该层中 Race x 的人口占比 */
        double[][] sumY = new double[M][R];
        double[]   pY   = new double[M];           // 层级总体比例

        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) {
                Neighborhood nb = grid[i][j];
                int y = nb.level;
                for (int x = 0; x < R; x++)
                    sumY[y][x] += nb.raceProps[x] / cells;
                pY[y] += 1.0 / cells;
            }

        /* ② 按公式：H(X|Y) = Σ_y p(y) * H(X | y) */
        double HXgY = 0;
        for (int y = 0; y < M; y++) {
            if (pY[y] == 0) continue;
            double h = 0;
            for (int x = 0; x < R; x++) {
                double p_x_given_y = sumY[y][x] / pY[y];   // π_jm
                if (p_x_given_y > 0)
                    h -= p_x_given_y * log2(p_x_given_y);
            }
            HXgY += pY[y] * h;
        }
        return HXgY;
    }

    public double computeEntropyRatio(){
        double hxy = computeConditionalEntropy();
        double hx  = computeTotalEntropy();
        return hx>0? 1 - hxy/hx : Double.NaN;
    }

    /* ───────────────────────────
     *  Helper – total pXY Table
     * ─────────────────────────── */
    private double[][] pXY() {                 // [y][x]
        int M = transit.levels.length;
        int R = races.length;
        double[][] out = new double[M][R];
        double cells = size * size;
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) {
                Neighborhood nb = grid[i][j];
                for (int x = 0; x < R; x++)
                    out[nb.level][x] += nb.raceProps[x] / cells;
            }
        return out;
    }
    private double[] pY(double[][] pXY) {
        int M = pXY.length;
        double[] out = new double[M];
        for (int y = 0; y < M; y++)
            for (double v : pXY[y]) out[y] += v;
        return out;
    }
    private double[] pX(double[][] pXY) {
        int R = pXY[0].length;
        double[] out = new double[R];
        for (double[] row : pXY)
            for (int x = 0; x < R; x++) out[x] += row[x];
        return out;
    }
    /* ───────────────────────────
     * 1. Duncan Dissimilarity
     * ─────────────────────────── */
    public double computeDissimilarity() {
        double[][] pXY = pXY();
        double[] pY = pY(pXY);
        double[] pX = pX(pXY);
        int M = pY.length, R = pX.length;

        double raw = 0;
        for (int y = 0; y < M; y++) if (pY[y] > 0)
            for (int x = 0; x < R; x++)
                raw += pY[y] * Math.abs(pXY[y][x]/pY[y] - pX[x]);

        double maxD = 0;
        for (double px : pX) maxD += 2 * px * (1 - px);
        return maxD > 0 ? raw / maxD : Double.NaN;
    }
    /* ───────────────────────────
     * 2. Gini (multi-group)
     * ─────────────────────────── */
    public double computeGini() {
        double[][] pXY = pXY();
        double[] pY = pY(pXY);
        double[] pX = pX(pXY);
        int M = pY.length, R = pX.length;

        double G = 0;
        for (int m = 0; m < R; m++) {
            for (int j = 0; j < M; j++) {
                if (pY[j] == 0) continue;               // ✅ 跳过空等级
                double pij = pXY[j][m] / pY[j];
                for (int j2 = 0; j2 < M; j2++) {
                    if (pY[j2] == 0) continue;          // ✅ 跳过空等级
                    double pij2 = pXY[j2][m] / pY[j2];
                    G += pY[j] * pY[j2] * Math.abs(pij - pij2) / 2.0;
                }
            }
        }

        double I = 0; for (double px : pX) I += px * (1 - px);
        return I > 0 ? G / I : Double.NaN;
    }


    /* ───────────────────────────
     * 3. Variation C
     * ─────────────────────────── */
    public double computeVariation() {
        double[][] pXY = pXY();
        double[] pY = pY(pXY);
        double[] pX = pX(pXY);
        int J = pY.length;   // transit bins
        int M = pX.length;   // races

        double C = 0;
        for (int x = 0; x < M; x++) {
            double denom = (M - 1) * pX[x];
            if (denom == 0) continue;
            for (int y = 0; y < J; y++) {
                double pij = (pY[y] > 0) ? pXY[y][x] / pY[y] : 0;
                C += pY[y] * Math.pow(pij - pX[x], 2) / denom;
            }
        }
        return C;
    }

    /* ───────────────────────────
     * 4. Relative Diversity R
     * ─────────────────────────── */
    public double computeRelativeDiversity() {
        double[][] pXY = pXY();
        double[] pY = pY(pXY);
        double[] pX = pX(pXY);
        int M = pY.length, R = pX.length;

        double I = 0; for (double px : pX) I += px*(1-px);
        double Rv = 0;
        for (int y = 0; y < M; y++) if (pY[y]>0)
            for (int x = 0; x < R; x++) {
                double pij = pXY[y][x]/pY[y];
                Rv += pY[y]*Math.pow(pij - pX[x],2)/I;
            }
        return Rv;
    }
    /* ───────────────────────────
     * 5. Exposure P
     * ─────────────────────────── */
    public double computeExposure() {
        double[][] pXY = pXY();
        double[] pY = pY(pXY);
        double[] pX = pX(pXY);
        int M = pY.length, R = pX.length;

        double P = 0;
        for (int y = 0; y < M; y++) if (pY[y]>0)
            for (int x = 0; x < R; x++) {
                double pij = pXY[y][x]/pY[y];
                double denom = 1 - pX[x];
                if (denom>0) P += pY[y]*Math.pow(pij - pX[x],2)/denom;
            }
        return P;
    }





    private double log2(double v){ return Math.log(v)/Math.log(2); }


    public static class Neighborhood {
            double[] raceProps;
            int level;
            double transitValue;

            Neighborhood(double[] props, int level, double tValue){
                this.raceProps     = props;
                this.level  = level;
                this.transitValue  = tValue;
            }
        }
}

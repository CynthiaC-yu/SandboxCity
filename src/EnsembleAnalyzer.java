import java.util.*;
import static java.lang.Math.*;

/** Massive-scale ensemble engine (no UI, streaming accumulators). */
public final class EnsembleAnalyzer {

    /* 直方图分箱模式 */
    public enum BinMode { FIXED, AUTO }

    /* 生成模式（映射到你已有的 SandboxCity.GenerationMode 与四种快捷档） */
    public enum GenPreset { LEGACY_UNIFORM, DIRICHLET1, LOW_SEG, HIGH_SEG, CUSTOM_HIER }

    /* 运行参数 */
    public static final class Config {
        public int runs = 100;              // 多少组城市
        public long cells = 1_000_000L;     // 每组城市的格子数（可上百万）
        public int races = 4;               // K
        public int levels = 5;              // T
        public GenPreset preset = GenPreset.CUSTOM_HIER;
        public double beta = 0.3;           // 仅 CUSTOM_HIER 用
        public double kappa = 0.6;          // 仅 CUSTOM_HIER 用
        public long seed = System.nanoTime();

        /* 阈值扫动（Monte-Carlo） */
        public boolean computeSweepVariance = true;
        public int sweepSamples = 50;       // 每次抽多少组阈值
        public double sweepResolution = 0.05; // 分位桶宽 Q = 1/res（也用于“等分分位”）
        public boolean useEqualQuantiles = true; // “每个等级占 1/T” 的约束
        public double minGap = 0.02;        // 随机阈值时的最小间距（绝对宽度）

        /* 直方图控制 */
        public BinMode binMode = BinMode.AUTO;
        public int fixedBins = 30;
    }

    /* 结果：六个指标的两组数组 */
    public static final class Result {
        public double[] H, D, G, C, R, P;        // 分割指数（每 run 一条）
        public double[] vH, vD, vG, vC, vR, vP;  // 阈值扫动方差（每 run 一条）
        public Config cfg;                       // 备查
    }

    /* 进度回调（可选） */
    public interface Progress {
        void onProgress(int done, int total);
    }

    /* ======= 主入口：运行 ensemble ======= */
    public static Result run(Config cfg, Progress cb) {
        Random rnd = new Random(cfg.seed);
        int K = cfg.races, T = cfg.levels;
        int Q = max(1, (int) round(1.0 / cfg.sweepResolution));  // transit 分位桶数

        double[] H = new double[cfg.runs], D = new double[cfg.runs], G = new double[cfg.runs],
                C = new double[cfg.runs], R = new double[cfg.runs], P = new double[cfg.runs];
        double[] vH = new double[cfg.runs], vD = new double[cfg.runs], vG = new double[cfg.runs],
                vC = new double[cfg.runs], vR = new double[cfg.runs], vP = new double[cfg.runs];

        for (int run = 0; run < cfg.runs; run++) {
            // ① 选择 π 与 κ（按预设）
            double[] pi; double kappaLocal;
            switch (cfg.preset) {
                case LEGACY_UNIFORM:
                    pi = null; kappaLocal = Double.NaN;
                    break;
                case DIRICHLET1:
                    pi = new double[K]; Arrays.fill(pi, 1.0 / K); kappaLocal = K;
                    break;
                case LOW_SEG:    // 城内很糊、中心近均匀
                    pi = new double[K]; Arrays.fill(pi, 1.0 / K);
                    kappaLocal = 5.0 * K; // α_i ≈5：非常混合
                    break;
                case HIGH_SEG:   // 城内很尖、中心可偏
                    pi = dirichletSym(K, 0.3, rnd);
                    kappaLocal = 0.6;
                    break;
                default:         // CUSTOM_HIER
                    pi = dirichletSym(K, cfg.beta, rnd);
                    kappaLocal = cfg.kappa;
            }

            // ② 累计到 transit 分位桶（Q×K），不创建任何格子对象
            double[][] binSums = new double[Q][K];  // 每桶各人群的“份额”总和
            for (long c = 0; c < cfg.cells; c++) {
                double[] rp;
                if (cfg.preset == GenPreset.LEGACY_UNIFORM) {
                    rp = normalizeUniform(K, rnd);
                } else {
                    rp = dirichlet(alphaFrom(pi, kappaLocal), rnd);
                }
                double t = rnd.nextDouble();
                int b = (int) min(Q - 1, floor(t * Q));
                double[] row = binSums[b];
                for (int k = 0; k < K; k++) row[k] += rp[k];
            }

            // ③ “等分分位”作为默认等级（T1 为最高分位）
            int[] cutsEq = new int[T - 1];
            int per = Q / T, rem = Q % T, idx = Q;
            for (int j = 0; j < T - 1; j++) {
                int width = per + (j < rem ? 1 : 0);
                idx -= width;
                cutsEq[j] = idx;        // 以“高→低”的边界索引
            }
            Metrics base = metricsFromBins(binSums, cutsEq, cfg.cells);
            H[run] = base.H; D[run] = base.D; G[run] = base.G; C[run] = base.C; R[run] = base.R; P[run] = base.P;

            // ④ 阈值扫动（Monte-Carlo）：抽 S 组 cuts，计算指标方差
            if (cfg.computeSweepVariance) {
                // 预构前缀和，便于 O(1) 区间求和
                double[][] pref = prefixFromBins(binSums);  // [Q+1][K]
                double[] sh = new double[cfg.sweepSamples],
                        sd = new double[cfg.sweepSamples],
                        sg = new double[cfg.sweepSamples],
                        sc = new double[cfg.sweepSamples],
                        sr = new double[cfg.sweepSamples],
                        sp = new double[cfg.sweepSamples];

                for (int s = 0; s < cfg.sweepSamples; s++) {
                    int[] cuts = (cfg.useEqualQuantiles)
                            ? cutsEq
                            : randomCutsByQuota(T, Q, (int) round(cfg.minGap * Q), rnd);
                    Metrics m = metricsFromBinsFast(pref, cuts, cfg.cells);
                    sh[s] = m.H; sd[s] = m.D; sg[s] = m.G; sc[s] = m.C; sr[s] = m.R; sp[s] = m.P;
                }
                vH[run] = variance(sh); vD[run] = variance(sd); vG[run] = variance(sg);
                vC[run] = variance(sc); vR[run] = variance(sr); vP[run] = variance(sp);
            } else {
                vH[run] = vD[run] = vG[run] = vC[run] = vR[run] = vP[run] = Double.NaN;
            }

            if (cb != null) cb.onProgress(run + 1, cfg.runs);
        }

        Result out = new Result();
        out.H = H; out.D = D; out.G = G; out.C = C; out.R = R; out.P = P;
        out.vH = vH; out.vD = vD; out.vG = vG; out.vC = vC; out.vR = vR; out.vP = vP;
        out.cfg = cfg;
        return out;
    }

    /* ======= ── 下方是算法小工具 ── ======= */

    private static double[] alphaFrom(double[] pi, double kappa) {
        if (pi == null) return null;
        double[] a = new double[pi.length];
        for (int i = 0; i < a.length; i++) a[i] = max(1e-8, pi[i] * kappa);
        return a;
    }

    private static double[][] prefixFromBins(double[][] bins){ // [Q+1][K]，pref[q] = sum_{0..q-1}
        int Q = bins.length, K = bins[0].length;
        double[][] pref = new double[Q + 1][K];
        for (int q = 1; q <= Q; q++) {
            for (int k = 0; k < K; k++)
                pref[q][k] = pref[q-1][k] + bins[q-1][k];
        }
        return pref;
    }

    /** cuts: 长度 T-1，按“从高到低”的 bin 边界索引（0..Q） */
    private static Metrics metricsFromBins(double[][] bins, int[] cuts, long cells) {
        return metricsFromBinsFast(prefixFromBins(bins), cuts, cells);
    }

    private static Metrics metricsFromBinsFast(double[][] pref, int[] cuts, long cells) {
        int Q = pref.length - 1, K = pref[0].length, T = cuts.length + 1;
        // 先把“高→低”的 cuts 转成每个 level 的 [lo, hi) 区间（lo/hi 均为 0..Q 的整数）
        int[] lo = new int[T], hi = new int[T];
        int curHi = Q;  // 最高端
        for (int t = 0; t < T; t++) {
            int curLo = (t == T - 1) ? 0 : cuts[t];
            lo[t] = curLo; hi[t] = curHi;
            curHi = curLo;
        }
        // pXY[y][x] = sum_{q in [lo_y, hi_y)} binSums[q][x] / cells
        double[][] pXY = new double[T][K];
        for (int y = 0; y < T; y++) {
            for (int k = 0; k < K; k++) {
                pXY[y][k] = (pref[hi[y]][k] - pref[lo[y]][k]) / (double) cells;
            }
        }
        return metricsFrom_pXY(pXY);
    }

    /** 完整指标（复用了你在 SandboxCity 里的公式） */
    private static Metrics metricsFrom_pXY(double[][] pXY) {
        int T = pXY.length, K = pXY[0].length;
        double[] pY = new double[T], pX = new double[K];
        for (int y = 0; y < T; y++) for (int k = 0; k < K; k++) { pY[y] += pXY[y][k]; pX[k] += pXY[y][k]; }

        // H ratio
        double HX = 0, HXgY = 0;
        for (int k = 0; k < K; k++) if (pX[k] > 0) HX -= pX[k] * log2(pX[k]);
        for (int y = 0; y < T; y++) if (pY[y] > 0) {
            double h = 0;
            for (int k = 0; k < K; k++) {
                double p = (pY[y] > 0) ? pXY[y][k] / pY[y] : 0;
                if (p > 0) h -= p * log2(p);
            }
            HXgY += pY[y] * h;
        }
        double Hratio = (HX > 0) ? 1 - HXgY / HX : Double.NaN;

        // Duncan D
        double raw = 0; for (int y = 0; y < T; y++) if (pY[y] > 0)
            for (int k = 0; k < K; k++) raw += pY[y] * Math.abs(pXY[y][k]/pY[y] - pX[k]);
        double maxD = 0; for (double px : pX) maxD += 2 * px * (1 - px);
        double D = maxD > 0 ? raw / maxD : Double.NaN;

        // Gini (multi-group)
        double G = 0;
        for (int m = 0; m < K; m++) {
            for (int j = 0; j < T; j++) if (pY[j] > 0) {
                double pij = pXY[j][m] / pY[j];
                for (int j2 = 0; j2 < T; j2++) if (pY[j2] > 0) {
                    double pij2 = pXY[j2][m] / pY[j2];
                    G += pY[j] * pY[j2] * Math.abs(pij - pij2) / 2.0;
                }
            }
        }
        double I = 0; for (double px : pX) I += px * (1 - px);
        G = I > 0 ? G / I : Double.NaN;

        // Variation C
        double C = 0;
        for (int m = 0; m < K; m++) {
            double denom = (K - 1) * pX[m];
            if (denom == 0) continue;
            for (int j = 0; j < T; j++) {
                double pij = (pY[j] > 0) ? pXY[j][m] / pY[j] : 0;
                C += pY[j] * (pij - pX[m]) * (pij - pX[m]) / denom;
            }
        }

        // Relative Diversity R
        double Rv = 0;
        for (int y = 0; y < T; y++) if (pY[y]>0)
            for (int x = 0; x < K; x++) {
                double pij = pXY[y][x]/pY[y];
                Rv += pY[y] * (pij - pX[x]) * (pij - pX[x]) / I;
            }

        // Exposure P
        double P = 0;
        for (int y = 0; y < T; y++) if (pY[y]>0)
            for (int x = 0; x < K; x++) {
                double pij = pXY[y][x]/pY[y];
                double denom = 1 - pX[x];
                if (denom>0) P += pY[y] * (pij - pX[x]) * (pij - pX[x]) / denom;
            }

        Metrics M = new Metrics();
        M.H = Hratio; M.D = D; M.G = G; M.C = C; M.R = Rv; M.P = P;
        return M;
    }

    private static int[] randomCutsByQuota(int T, int Q, int minGapBins, Random rnd) {
        // Q 个分位桶被分成 T 段：每段目标宽度 ≈ Q/T（处理余数）
        int[] cuts = new int[T - 1];          // 高→低的边界索引
        int per = Q / T, rem = Q % T;
        int hi = Q;                           // 从高端往低端走
        for (int j = 0; j < T - 1; j++) {
            int width = per + (j < rem ? 1 : 0);
            int target = hi - width;          // 等分边界
            // 在这一段里给边界一点“抖动”，并保证与下一段至少有 minGapBins
            int jitter = Math.max(0, width / 3);                // 抖动尺度（可调）
            int loBound = Math.max(hi - width + minGapBins, target - jitter);
            int hiBound = Math.min(hi - minGapBins,             target + jitter);
            if (hiBound < loBound) hiBound = loBound;           // 兜底
            cuts[j] = loBound + rnd.nextInt(Math.max(1, hiBound - loBound + 1));
            hi = cuts[j];
        }
        return cuts;
    }

    /* 指标容器 */
    private static final class Metrics { double H, D, G, C, R, P; }

    /* ===== Gamma/Dirichlet & 旧 Uniform 采样 —— 与 SandboxCity 一致 ===== */
    private static double[] dirichletSym(int K, double a, Random r){
        double[] alpha = new double[K]; Arrays.fill(alpha, a);
        return dirichlet(alpha, r);
    }
    private static double[] dirichlet(double[] alpha, Random r){
        int K = alpha.length; double[] x = new double[K]; double s = 0;
        for (int i = 0; i < K; i++) { x[i] = gamma(max(1e-8, alpha[i]), r); s += x[i]; }
        for (int i = 0; i < K; i++) x[i] /= s; return x;
    }
    private static double gamma(double a, Random r){
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
    private static double[] normalizeUniform(int n, Random r){
        double[] a = new double[n]; double s = 0;
        for (int i = 0; i < n; i++) { a[i] = r.nextDouble(); s += a[i]; }
        for (int i = 0; i < n; i++) a[i] /= s;
        return a;
    }
    private static double variance(double[] v){
        int n = 0; double s=0, s2=0;
        for (double x: v) if (!Double.isNaN(x) && !Double.isInfinite(x)) { n++; s+=x; s2+=x*x; }
        if (n <= 1) return Double.NaN;
        double m = s/n;
        return (s2/n - m*m);
    }
    private static double log2(double v){ return Math.log(v)/Math.log(2); }
}

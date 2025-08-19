import javax.swing.*;
import java.awt.*;
import java.util.concurrent.ExecutionException;

public class EnsembleDialog extends JDialog {

    private final SandboxCityGUI owner;
    // controls
    private JSpinner spRuns, spCells, spRaces, spLevels, spBeta, spKappa, spSweepS, spRes, spMinGap, spFixedBins;
    private JComboBox<String> cmbPreset, cmbBinMode;
    private JCheckBox ckSweep, ckEqualQ;
    private JButton btnRun;
    private JPanel histGrid;

    public EnsembleDialog(SandboxCityGUI owner) {
        super(owner, "Ensemble Analysis", false);
        this.owner = owner;
        setLayout(new BorderLayout(8,8));

        add(buildControls(), BorderLayout.NORTH);
        add(buildHistGrid(), BorderLayout.CENTER);

        setSize(1100, 820);
        setLocationRelativeTo(owner);
    }

    private JPanel buildControls() {
        JPanel p = new JPanel(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.insets = new Insets(4,6,4,6); c.anchor = GridBagConstraints.WEST;

        int row = 0;
        // row0: runs, cells
        c.gridx=0; c.gridy=row; p.add(new JLabel("Runs:"), c);
        spRuns = new JSpinner(new SpinnerNumberModel(200, 10, 200000, 10)); p.add(spRuns, next(c));
        p.add(new JLabel("Cells per run:"), next(c));
        spCells = new JSpinner(new SpinnerNumberModel(1_000_000, 10_000, 50_000_000, 50_000)); p.add(spCells, next(c));

        // row1: races, levels
        row++; c.gridx=0; c.gridy=row; p.add(new JLabel("Races (K):"), c);
        int curK = (owner.getCity()!=null)? owner.getCity().races.length : 4;
        spRaces = new JSpinner(new SpinnerNumberModel(curK, 2, 12, 1)); p.add(spRaces, next(c));
        p.add(new JLabel("Transit levels (T):"), next(c));
        int curT = (owner.getCity()!=null)? owner.getCity().transit.levels.length : 5;
        spLevels = new JSpinner(new SpinnerNumberModel(curT, 2, 12, 1)); p.add(spLevels, next(c));

        // row2: preset + beta/kappa
        row++; c.gridx=0; c.gridy=row; p.add(new JLabel("Mode:"), c);
        cmbPreset = new JComboBox<>(new String[]{
                "Custom (β, κ)", "Legacy uniform", "Dirichlet(1)", "Low Segregation", "High Segregation"
        });
        p.add(cmbPreset, next(c));
        p.add(new JLabel("β:"), next(c));
        spBeta  = new JSpinner(new SpinnerNumberModel(0.3, 0.01, 50.0, 0.05)); p.add(spBeta, next(c));
        p.add(new JLabel("κ:"), next(c));
        spKappa = new JSpinner(new SpinnerNumberModel(0.6, 0.01, 200.0, 0.05)); p.add(spKappa, next(c));

        cmbPreset.addActionListener(e -> {
            boolean custom = cmbPreset.getSelectedIndex()==0;
            spBeta.setEnabled(custom); spKappa.setEnabled(custom);
        });
        cmbPreset.setSelectedIndex(0);

        // row3: sweep
        row++; c.gridx=0; c.gridy=row; p.add(new JLabel("Sweep variance:"), c);
        ckSweep = new JCheckBox("enable", true); p.add(ckSweep, next(c));
        p.add(new JLabel("samples:"), next(c));
        spSweepS = new JSpinner(new SpinnerNumberModel(50, 5, 5000, 5)); p.add(spSweepS, next(c));
        p.add(new JLabel("resolution (1/Q):"), next(c));
        spRes = new JSpinner(new SpinnerNumberModel(0.05, 0.005, 0.5, 0.005)); p.add(spRes, next(c));

        // row4: transit thresholds restraint
        row++; c.gridx=0; c.gridy=row; p.add(new JLabel("Transit levels:"), c);
        ckEqualQ = new JCheckBox("equal quantiles (1/T each)", true); p.add(ckEqualQ, next(c));
        p.add(new JLabel("min-gap (if random):"), next(c));
        spMinGap = new JSpinner(new SpinnerNumberModel(0.02, 0.0, 0.5, 0.005)); p.add(spMinGap, next(c));

        // row5: histogram binning
        row++; c.gridx=0; c.gridy=row; p.add(new JLabel("Histogram binning:"), c);
        cmbBinMode = new JComboBox<>(new String[]{"Auto Adjusted","Fixed"});
        p.add(cmbBinMode, next(c));
        p.add(new JLabel("bins (Fixed):"), next(c));
        spFixedBins = new JSpinner(new SpinnerNumberModel(30, 5, 200, 1)); p.add(spFixedBins, next(c));
        cmbBinMode.addActionListener(e -> spFixedBins.setEnabled(cmbBinMode.getSelectedIndex()==1));
        spFixedBins.setEnabled(false);

        // row6: run
        row++; c.gridx=0; c.gridy=row; btnRun = new JButton("Run Ensemble");
        btnRun.addActionListener(e -> runAsync()); p.add(btnRun, c);

        return p;
    }

    private JPanel buildHistGrid() {
        histGrid = new JPanel(new GridLayout(3, 4, 8, 8)); // 12 charts
        for (int i=0;i<12;i++) histGrid.add(new HistogramPanel());
        JScrollPane sp = new JScrollPane(histGrid);
        sp.setBorder(BorderFactory.createTitledBorder("Distributions"));
        JPanel wrap = new JPanel(new BorderLayout());
        wrap.add(sp, BorderLayout.CENTER);
        return wrap;
    }

    private void runAsync() {
        EnsembleAnalyzer.Config cfg = new EnsembleAnalyzer.Config();
        cfg.runs   = (Integer) spRuns.getValue();
        cfg.cells  = ((Number) spCells.getValue()).longValue();
        cfg.races  = (Integer) spRaces.getValue();
        cfg.levels = (Integer) spLevels.getValue();

        switch (cmbPreset.getSelectedIndex()) {
            case 1 -> cfg.preset = EnsembleAnalyzer.GenPreset.LEGACY_UNIFORM;
            case 2 -> cfg.preset = EnsembleAnalyzer.GenPreset.DIRICHLET1;
            case 3 -> cfg.preset = EnsembleAnalyzer.GenPreset.LOW_SEG;
            case 4 -> cfg.preset = EnsembleAnalyzer.GenPreset.HIGH_SEG;
            default -> cfg.preset = EnsembleAnalyzer.GenPreset.CUSTOM_HIER;
        }
        cfg.beta  = ((Number) spBeta.getValue()).doubleValue();
        cfg.kappa = ((Number) spKappa.getValue()).doubleValue();

        cfg.computeSweepVariance = ckSweep.isSelected();
        cfg.sweepSamples = (Integer) spSweepS.getValue();
        cfg.sweepResolution = ((Number) spRes.getValue()).doubleValue();
        cfg.useEqualQuantiles = ckEqualQ.isSelected();
        cfg.minGap = ((Number) spMinGap.getValue()).doubleValue();

        cfg.binMode = (cmbBinMode.getSelectedIndex()==0)? EnsembleAnalyzer.BinMode.AUTO : EnsembleAnalyzer.BinMode.FIXED;
        cfg.fixedBins = (Integer) spFixedBins.getValue();

        btnRun.setEnabled(false);

        SwingWorker<EnsembleAnalyzer.Result, Void> wk = new SwingWorker<>() {
            @Override protected EnsembleAnalyzer.Result doInBackground() {
                return EnsembleAnalyzer.run(cfg, (done,total) -> setTitle("Ensemble Analysis  ("+done+"/"+total+")"));
            }
            @Override protected void done() {
                try {
                    EnsembleAnalyzer.Result r = get();
                    updateHists(r);
                } catch (InterruptedException | ExecutionException ex) {
                    JOptionPane.showMessageDialog(EnsembleDialog.this, "Failed: "+ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
                } finally {
                    btnRun.setEnabled(true);
                    setTitle("Ensemble Analysis");
                }
            }
        };
        wk.execute();
    }

    private void updateHists(EnsembleAnalyzer.Result r) {
        boolean auto = r.cfg.binMode == EnsembleAnalyzer.BinMode.AUTO;
        int bins = r.cfg.fixedBins;

        String[] names = {"Theil H","Dissimilarity D","Gini","Variation C","Rel. Diversity R","Exposure P"};
        double[][] A = { r.H, r.D, r.G, r.C, r.R, r.P };
        double[][] B = { r.vH, r.vD, r.vG, r.vC, r.vR, r.vP };

        for (int i=0;i<6;i++) {
            ((HistogramPanel)histGrid.getComponent(i)).setData("Score — " + names[i], A[i], auto, bins);
        }
        for (int i=0;i<6;i++) {
            ((HistogramPanel)histGrid.getComponent(6+i)).setData("Sweep Var — " + names[i], B[i], auto, bins);
        }
        histGrid.revalidate(); histGrid.repaint();
    }

    // 供外部调用
    public void open() { setVisible(true); }

    // 让外部能拿到 city
    public SandboxCityGUI getOwnerGUI(){ return owner; }

    // GridBag 小助手：把 c.gridx 往右移一格并返回同一个对象。
    private GridBagConstraints next(GridBagConstraints c) {
        c.gridx++;
        return c;
    }

}

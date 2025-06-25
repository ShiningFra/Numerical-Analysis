package FiniteVolume1;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Line2D;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

public class ODEFiniteVolume extends JFrame {
    
    // Interfaces fonctionnelles pour les fonctions source
    public interface SourceFunction {
        double apply(double x);
        String getName();
    }
    
    // Fonctions source disponibles
    public static class SinPiX implements SourceFunction {
        public double apply(double x) { return Math.sin(Math.PI * x); }
        public String getName() { return "sin(πx)"; }
    }
    
    public static class XCubed implements SourceFunction {
        public double apply(double x) { return x * x * x; }
        public String getName() { return "x³"; }
    }
    
    // Classe abstraite pour les équations différentielles
    public abstract static class ODESolver {
        protected int n;
        protected double h;
        protected SourceFunction f;
        protected double[] x;
        protected double[] u_numerical;
        protected double[] u_exact;
        
        public ODESolver(int n, SourceFunction f) {
            this.n = n;
            this.h = 1.0 / (n + 1);
            this.f = f;
            this.x = new double[n + 2];
            this.u_numerical = new double[n + 2];
            this.u_exact = new double[n + 2];
            
            // Initialisation du maillage
            for (int i = 0; i <= n + 1; i++) {
                x[i] = i * h;
            }
        }
        
        public abstract void solve();
        public abstract void computeExactSolution();
        public abstract String getEquationType();
        
        public double computeError() {
            double error = 0.0;
            for (int i = 1; i <= n; i++) {
                error += Math.pow(u_numerical[i] - u_exact[i], 2);
            }
            return Math.sqrt(error * h);
        }
        
        public double[] getX() { return x; }
        public double[] getNumerical() { return u_numerical; }
        public double[] getExact() { return u_exact; }
    }
    
    // Solver pour -u'' + u = f
    public static class Solver1 extends ODESolver {
        public Solver1(int n, SourceFunction f) {
            super(n, f);
        }
        
        @Override
        public void solve() {
            // Matrice tridiagonale pour -u'' + u = f
            double[] a = new double[n];  // sous-diagonale
            double[] b = new double[n];  // diagonale
            double[] c = new double[n];  // sur-diagonale
            double[] d = new double[n];  // second membre
            
            for (int i = 0; i < n; i++) {
                a[i] = (i > 0) ? -1.0 / (h * h) : 0.0;
                b[i] = 2.0 / (h * h) + 1.0;
                c[i] = (i < n - 1) ? -1.0 / (h * h) : 0.0;
                d[i] = f.apply(x[i + 1]);
            }
            
            solveTridiagonal(a, b, c, d);
        }
        
        @Override
        public void computeExactSolution() {
            if (f instanceof SinPiX) {
                for (int i = 0; i <= n + 1; i++) {
                    u_exact[i] = Math.sin(Math.PI * x[i]) / (Math.PI * Math.PI + 1);
                }
            } else if (f instanceof XCubed) {
                for (int i = 0; i <= n + 1; i++) {
                    double xi = x[i];
                    u_exact[i] = xi*xi*xi - 6*xi;
                }
            }
        }
        
        @Override
        public String getEquationType() {
            return "-u'' + u = " + f.getName();
        }
        
        private void solveTridiagonal(double[] a, double[] b, double[] c, double[] d) {
            // Thomas algorithm
            for (int i = 1; i < n; i++) {
                double m = a[i] / b[i - 1];
                b[i] = b[i] - m * c[i - 1];
                d[i] = d[i] - m * d[i - 1];
            }
            
            u_numerical[n] = d[n - 1] / b[n - 1];
            for (int i = n - 1; i >= 1; i--) {
                u_numerical[i] = (d[i - 1] - c[i - 1] * u_numerical[i + 1]) / b[i - 1];
            }
        }
    }
    
    // Solver pour -u'' + u' = f
    public static class Solver2 extends ODESolver {
        public Solver2(int n, SourceFunction f) {
            super(n, f);
        }
        
        @Override
        public void solve() {
            double[] a = new double[n];
            double[] b = new double[n];
            double[] c = new double[n];
            double[] d = new double[n];
            
            for (int i = 0; i < n; i++) {
                a[i] = (i > 0) ? (-1.0 / (h * h) - 1.0 / (2 * h)) : 0.0;
                b[i] = 2.0 / (h * h);
                c[i] = (i < n - 1) ? (-1.0 / (h * h) + 1.0 / (2 * h)) : 0.0;
                d[i] = f.apply(x[i + 1]);
            }
            
            solveTridiagonal(a, b, c, d);
        }
        
        @Override
        public void computeExactSolution() {
            if (f instanceof SinPiX) {
                for (int i = 0; i <= n + 1; i++) {
                    double xi = x[i];
                    u_exact[i] = (Math.PI * Math.sin(Math.PI * xi) - Math.cos(Math.PI * xi) + 1) / (Math.PI * Math.PI + Math.PI);
                }
            } else if (f instanceof XCubed) {
                for (int i = 0; i <= n + 1; i++) {
                    double xi = x[i];
                    u_exact[i] = xi*xi*xi*xi/4 - 6*xi*xi + 12*xi - 12;
                }
            }
        }
        
        @Override
        public String getEquationType() {
            return "-u'' + u' = " + f.getName();
        }
        
        private void solveTridiagonal(double[] a, double[] b, double[] c, double[] d) {
            for (int i = 1; i < n; i++) {
                double m = a[i] / b[i - 1];
                b[i] = b[i] - m * c[i - 1];
                d[i] = d[i] - m * d[i - 1];
            }
            
            u_numerical[n] = d[n - 1] / b[n - 1];
            for (int i = n - 1; i >= 1; i--) {
                u_numerical[i] = (d[i - 1] - c[i - 1] * u_numerical[i + 1]) / b[i - 1];
            }
        }
    }
    
    // Solver pour -u'' = f
    public static class Solver3 extends ODESolver {
        public Solver3(int n, SourceFunction f) {
            super(n, f);
        }
        
        @Override
        public void solve() {
            double[] a = new double[n];
            double[] b = new double[n];
            double[] c = new double[n];
            double[] d = new double[n];
            
            for (int i = 0; i < n; i++) {
                a[i] = (i > 0) ? -1.0 / (h * h) : 0.0;
                b[i] = 2.0 / (h * h);
                c[i] = (i < n - 1) ? -1.0 / (h * h) : 0.0;
                d[i] = f.apply(x[i + 1]);
            }
            
            solveTridiagonal(a, b, c, d);
        }
        
        @Override
        public void computeExactSolution() {
            if (f instanceof SinPiX) {
                for (int i = 0; i <= n + 1; i++) {
                    u_exact[i] = -Math.sin(Math.PI * x[i]) / (Math.PI * Math.PI);
                }
            } else if (f instanceof XCubed) {
                for (int i = 0; i <= n + 1; i++) {
                    double xi = x[i];
                    u_exact[i] = xi*xi*xi*xi*xi/20 - xi/6;
                }
            }
        }
        
        @Override
        public String getEquationType() {
            return "-u'' = " + f.getName();
        }
        
        private void solveTridiagonal(double[] a, double[] b, double[] c, double[] d) {
            for (int i = 1; i < n; i++) {
                double m = a[i] / b[i - 1];
                b[i] = b[i] - m * c[i - 1];
                d[i] = d[i] - m * d[i - 1];
            }
            
            u_numerical[n] = d[n - 1] / b[n - 1];
            for (int i = n - 1; i >= 1; i--) {
                u_numerical[i] = (d[i - 1] - c[i - 1] * u_numerical[i + 1]) / b[i - 1];
            }
        }
    }
    
    // Classe pour les graphiques
    public static class PlotPanel extends JPanel {
        private List<PlotData> plots;
        private String title;
        private String xLabel;
        private String yLabel;
        
        public static class PlotData {
            double[] x, y;
            Color color;
            String label;
            boolean dashed;
            
            public PlotData(double[] x, double[] y, Color color, String label, boolean dashed) {
                this.x = x;
                this.y = y;
                this.color = color;
                this.label = label;
                this.dashed = dashed;
            }
        }
        
        public PlotPanel(String title, String xLabel, String yLabel) {
            this.plots = new ArrayList<>();
            this.title = title;
            this.xLabel = xLabel;
            this.yLabel = yLabel;
            setPreferredSize(new Dimension(800, 600));
            setBackground(Color.WHITE);
        }
        
        public void addPlot(double[] x, double[] y, Color color, String label, boolean dashed) {
            plots.add(new PlotData(x, y, color, label, dashed));
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            if (plots.isEmpty()) return;
            
            int margin = 80;
            int plotWidth = getWidth() - 2 * margin;
            int plotHeight = getHeight() - 2 * margin;
            
            // Trouver les limites
            double xMin = Double.MAX_VALUE, xMax = -Double.MAX_VALUE;
            double yMin = Double.MAX_VALUE, yMax = -Double.MAX_VALUE;
            
            for (PlotData plot : plots) {
                for (int i = 0; i < plot.x.length; i++) {
                    xMin = Math.min(xMin, plot.x[i]);
                    xMax = Math.max(xMax, plot.x[i]);
                    yMin = Math.min(yMin, plot.y[i]);
                    yMax = Math.max(yMax, plot.y[i]);
                }
            }
            
            // Ajouter une marge
            double xRange = xMax - xMin;
            double yRange = yMax - yMin;
            xMin -= 0.1 * xRange;
            xMax += 0.1 * xRange;
            yMin -= 0.1 * yRange;
            yMax += 0.1 * yRange;
            
            // Dessiner les axes
            g2d.setColor(Color.BLACK);
            g2d.setStroke(new BasicStroke(2));
            g2d.drawLine(margin, margin, margin, margin + plotHeight);
            g2d.drawLine(margin, margin + plotHeight, margin + plotWidth, margin + plotHeight);
            
            // Dessiner les graduations
            for (int i = 0; i <= 10; i++) {
                int x = margin + i * plotWidth / 10;
                int y = margin + plotHeight + 5;
                g2d.drawLine(x, margin + plotHeight, x, y);
                
                double xVal = xMin + i * (xMax - xMin) / 10;
                g2d.drawString(String.format("%.2f", xVal), x - 15, y + 15);
            }
            
            for (int i = 0; i <= 10; i++) {
                int x = margin - 5;
                int y = margin + plotHeight - i * plotHeight / 10;
                g2d.drawLine(x, y, margin, y);
                
                double yVal = yMin + i * (yRange) / 10;
                g2d.drawString(String.format("%.3f", yVal), x - 40, y + 5);
            }
            
            // Dessiner les courbes
            for (PlotData plot : plots) {
                g2d.setColor(plot.color);
                if (plot.dashed) {
                    g2d.setStroke(new BasicStroke(2, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{5}, 0));
                } else {
                    g2d.setStroke(new BasicStroke(2));
                }
                
                for (int i = 0; i < plot.x.length - 1; i++) {
                    int x1 = margin + (int) ((plot.x[i] - xMin) / (xMax - xMin) * plotWidth);
                    int y1 = margin + plotHeight - (int) ((plot.y[i] - yMin) / (yMax - yMin) * plotHeight);
                    int x2 = margin + (int) ((plot.x[i + 1] - xMin) / (xMax - xMin) * plotWidth);
                    int y2 = margin + plotHeight - (int) ((plot.y[i + 1] - yMin) / (yMax - yMin) * plotHeight);
                    
                    g2d.drawLine(x1, y1, x2, y2);
                }
            }
            
            // Titre et labels
            g2d.setColor(Color.BLACK);
            g2d.setFont(new Font("Arial", Font.BOLD, 16));
            FontMetrics fm = g2d.getFontMetrics();
            int titleWidth = fm.stringWidth(title);
            g2d.drawString(title, (getWidth() - titleWidth) / 2, 30);
            
            g2d.setFont(new Font("Arial", Font.PLAIN, 12));
            g2d.drawString(xLabel, getWidth() / 2 - 20, getHeight() - 20);
            
            // Y label (vertical)
            Graphics2D g2dCopy = (Graphics2D) g2d.create();
            g2dCopy.rotate(-Math.PI / 2);
            g2dCopy.drawString(yLabel, -getHeight() / 2 - 20, 20);
            g2dCopy.dispose();
            
            // Légende
            int legendY = margin + 20;
            for (PlotData plot : plots) {
                g2d.setColor(plot.color);
                if (plot.dashed) {
                    g2d.setStroke(new BasicStroke(2, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{5}, 0));
                } else {
                    g2d.setStroke(new BasicStroke(2));
                }
                g2d.drawLine(margin + plotWidth - 150, legendY, margin + plotWidth - 120, legendY);
                g2d.setColor(Color.BLACK);
                g2d.drawString(plot.label, margin + plotWidth - 115, legendY + 5);
                legendY += 20;
            }
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            // Choisir la fonction source
            SourceFunction[] functions = {new SinPiX(), new XCubed()};
            String[] functionNames = {"sin(πx)", "x³"};
            
            int choice = JOptionPane.showOptionDialog(null,
                "Choisissez la fonction source f(x):",
                "Sélection de fonction",
                JOptionPane.DEFAULT_OPTION,
                JOptionPane.QUESTION_MESSAGE,
                null,
                functionNames,
                functionNames[0]);
            
            if (choice == -1) return;
            
            SourceFunction selectedFunction = functions[choice];
            
            // Tester les trois équations
            Class<?>[] solverClasses = {/*Solver1.class, Solver2.class, */Solver3.class};
            int[] meshSizes = {10, 20, 40, 80, 160, 320, 640};
            
            for (Class<?> solverClass : solverClasses) {
                try {
                    // Calcul des erreurs et ordre de convergence
                    System.out.println("\n=== " + solverClass.getSimpleName() + " ===");
                    
                    List<Double> errors = new ArrayList<>();
                    List<Integer> sizes = new ArrayList<>();
                    
                    for (int size : meshSizes) {
                        ODESolver solver = (ODESolver) solverClass.getConstructor(int.class, SourceFunction.class)
                                                                .newInstance(size, selectedFunction);
                        solver.solve();
                        solver.computeExactSolution();
                        
                        double error = solver.computeError();
                        errors.add(error);
                        sizes.add(size);
                        
                        System.out.printf("Mailles: %d, Erreur: %.6e\n", size, error);
                    }
                    
                    // Calcul de l'ordre de convergence
                    if (errors.size() >= 2) {
                        double order = Math.log(errors.get(1) / errors.get(0)) / Math.log((double) sizes.get(0) / sizes.get(1));
                        System.out.printf("Ordre de convergence: %.2f\n", order);
                    }
                    
                    // Graphique des solutions pour 50 mailles
                    ODESolver solver50 = (ODESolver) solverClass.getConstructor(int.class, SourceFunction.class)
                                                               .newInstance(640, selectedFunction);
                    solver50.solve();
                    solver50.computeExactSolution();
                    
                    JFrame frame = new JFrame("Solutions - " + solver50.getEquationType());
                    PlotPanel plotPanel = new PlotPanel(
                        solver50.getEquationType() + " (640 mailles)",
                        "x", "u(x)"
                    );
                    
                    plotPanel.addPlot(solver50.getX(), solver50.getExact(), Color.BLUE, "Solution exacte", false);
                    plotPanel.addPlot(solver50.getX(), solver50.getNumerical(), Color.RED, "Solution numérique", true);
                    
                    frame.add(plotPanel);
                    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                    frame.pack();
                    frame.setVisible(true);
                    
                    // Graphique de convergence
                    JFrame convergenceFrame = new JFrame("Convergence - " + solver50.getEquationType());
                    PlotPanel convergencePlot = new PlotPanel(
                        "Évolution de l'erreur - " + solver50.getEquationType(),
                        "Nombre de mailles", "Erreur L2"
                    );
                    
                    double[] meshSizesDouble = sizes.stream().mapToDouble(Integer::doubleValue).toArray();
                    double[] errorsArray = errors.stream().mapToDouble(Double::doubleValue).toArray();
                    
                    convergencePlot.addPlot(meshSizesDouble, errorsArray, Color.GREEN, "Erreur L2", false);
                    
                    convergenceFrame.add(convergencePlot);
                    convergenceFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                    convergenceFrame.pack();
                    convergenceFrame.setLocationRelativeTo(null);
                    convergenceFrame.setVisible(true);
                    
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }
}
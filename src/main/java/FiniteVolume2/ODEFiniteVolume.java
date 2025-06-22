package FiniteVolume2;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;

public class ODEFiniteVolume extends JFrame {
    private static final double PI = Math.PI;
    private static final int PANEL_SIZE = 600;
    
    // Classe abstraite pour les différents types d'équations
    abstract static class ODESolver {
        protected int n;
        protected double h;
        protected double[][] u;
        protected double[][] f;
        protected double[][] exactSolution;
        
        public ODESolver(int n) {
            this.n = n;
            this.h = 1.0 / (n - 1);
            this.u = new double[n][n];
            this.f = new double[n][n];
            this.exactSolution = new double[n][n];
            initializeSourceTerm();
            computeExactSolution();
        }
        
        private void initializeSourceTerm() {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double x = i * h;
                    double y = j * h;
                    f[i][j] = Math.sin(PI * x) + Math.pow(y, 3);
                }
            }
        }
        
        protected abstract void computeExactSolution();
        protected abstract void solve();
        protected abstract String getName();
        
        public double computeError() {
            double error = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    error += Math.pow(u[i][j] - exactSolution[i][j], 2);
                }
            }
            return Math.sqrt(error * h * h);
        }
        
        public double[][] getSolution() { return u; }
        public double[][] getExactSolution() { return exactSolution; }
        public int getN() { return n; }
        public double getH() { return h; }
    }
    
    // Classe pour -u'' + u = f
    static class PoissonHelmholtz extends ODESolver {
        public PoissonHelmholtz(int n) {
            super(n);
        }
        
        @Override
        protected void computeExactSolution() {
            // Solution exacte approximée pour -u'' + u = sin(πx) + y³
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double x = i * h;
                    double y = j * h;
                    // Solution approximée basée sur une série
                    exactSolution[i][j] = Math.sin(PI * x) / (PI * PI + 1) + 
                                        Math.pow(y, 3) / (1 + 6);
                }
            }
        }
        
        @Override
        protected void solve() {
            // Méthode des volumes finis pour -u'' + u = f
            // Système linéaire Ax = b
            int size = n * n;
            double[][] A = new double[size][size];
            double[] b = new double[size];
            
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int idx = i * n + j;
                    
                    if (i == 0 || i == n-1 || j == 0 || j == n-1) {
                        // Conditions aux limites (Dirichlet u = 0)
                        A[idx][idx] = 1.0;
                        b[idx] = 0.0;
                    } else {
                        // Points intérieurs
                        double coeff = 4.0/(h*h) + 1.0;
                        A[idx][idx] = coeff;
                        A[idx][idx-n] = -1.0/(h*h); // u[i-1][j]
                        A[idx][idx+n] = -1.0/(h*h); // u[i+1][j]
                        A[idx][idx-1] = -1.0/(h*h); // u[i][j-1]
                        A[idx][idx+1] = -1.0/(h*h); // u[i][j+1]
                        b[idx] = f[i][j];
                    }
                }
            }
            
            // Résolution du système linéaire par Gauss-Seidel
            double[] x = new double[size];
            gaussSeidel(A, b, x, 1000, 1e-10);
            
            // Récupération de la solution
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    u[i][j] = x[i * n + j];
                }
            }
        }
        
        @Override
        protected String getName() {
            return "-u'' + u = f";
        }
    }
    
    // Classe pour -u'' + u' = f
    static class ConvectionDiffusion extends ODESolver {
        public ConvectionDiffusion(int n) {
            super(n);
        }
        
        @Override
        protected void computeExactSolution() {
            // Solution exacte approximée pour -u'' + u' = sin(πx) + y³
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double x = i * h;
                    double y = j * h;
                    exactSolution[i][j] = Math.sin(PI * x) / (PI * PI - PI) + 
                                        Math.pow(y, 3) / 6;
                }
            }
        }
        
        @Override
        protected void solve() {
            int size = n * n;
            double[][] A = new double[size][size];
            double[] b = new double[size];
            
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int idx = i * n + j;
                    
                    if (i == 0 || i == n-1 || j == 0 || j == n-1) {
                        A[idx][idx] = 1.0;
                        b[idx] = 0.0;
                    } else {
                        double coeff = 4.0/(h*h);
                        A[idx][idx] = coeff;
                        A[idx][idx-n] = -1.0/(h*h) - 0.5/h; // u[i-1][j]
                        A[idx][idx+n] = -1.0/(h*h) + 0.5/h; // u[i+1][j]
                        A[idx][idx-1] = -1.0/(h*h) - 0.5/h; // u[i][j-1]
                        A[idx][idx+1] = -1.0/(h*h) + 0.5/h; // u[i][j+1]
                        b[idx] = f[i][j];
                    }
                }
            }
            
            double[] x = new double[size];
            gaussSeidel(A, b, x, 1000, 1e-10);
            
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    u[i][j] = x[i * n + j];
                }
            }
        }
        
        @Override
        protected String getName() {
            return "-u'' + u' = f";
        }
    }
    
    // Classe pour -u'' = f (Poisson pur)
    static class Poisson extends ODESolver {
        public Poisson(int n) {
            super(n);
        }
        
        @Override
        protected void computeExactSolution() {
            // Solution exacte pour -u'' = sin(πx) + y³
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double x = i * h;
                    double y = j * h;
                    exactSolution[i][j] = Math.sin(PI * x) / (PI * PI) + 
                                        Math.pow(y, 3) / 6;
                }
            }
        }
        
        @Override
        protected void solve() {
            int size = n * n;
            double[][] A = new double[size][size];
            double[] b = new double[size];
            
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int idx = i * n + j;
                    
                    if (i == 0 || i == n-1 || j == 0 || j == n-1) {
                        A[idx][idx] = 1.0;
                        b[idx] = 0.0;
                    } else {
                        A[idx][idx] = 4.0/(h*h);
                        A[idx][idx-n] = -1.0/(h*h);
                        A[idx][idx+n] = -1.0/(h*h);
                        A[idx][idx-1] = -1.0/(h*h);
                        A[idx][idx+1] = -1.0/(h*h);
                        b[idx] = f[i][j];
                    }
                }
            }
            
            double[] x = new double[size];
            gaussSeidel(A, b, x, 1000, 1e-10);
            
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    u[i][j] = x[i * n + j];
                }
            }
        }
        
        @Override
        protected String getName() {
            return "-u'' = f";
        }
    }
    
    // Méthode de Gauss-Seidel pour résoudre le système linéaire
    static void gaussSeidel(double[][] A, double[] b, double[] x, int maxIter, double tol) {
        int n = A.length;
        double[] xOld = new double[n];
        
        for (int iter = 0; iter < maxIter; iter++) {
            System.arraycopy(x, 0, xOld, 0, n);
            
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        sum += A[i][j] * x[j];
                    }
                }
                x[i] = (b[i] - sum) / A[i][i];
            }
            
            // Vérification de la convergence
            double error = 0.0;
            for (int i = 0; i < n; i++) {
                error += Math.pow(x[i] - xOld[i], 2);
            }
            if (Math.sqrt(error) < tol) break;
        }
    }
    
    // Panneau de visualisation avec colormap et info sur survol
    static class VisualizationPanel extends JPanel implements MouseMotionListener {
        private double[][] data;
        private String title;
        private BufferedImage image;
        private int mouseX = -1, mouseY = -1;
        private double mouseValue = 0.0;
        
        public VisualizationPanel(double[][] data, String title) {
            this.data = data;
            this.title = title;
            setPreferredSize(new Dimension(PANEL_SIZE, PANEL_SIZE));
            addMouseMotionListener(this);
            createImage();
        }
        
        private void createImage() {
            int n = data.length;
            image = new BufferedImage(PANEL_SIZE, PANEL_SIZE, BufferedImage.TYPE_INT_RGB);
            Graphics2D g2d = image.createGraphics();
            
            // Trouver min et max pour la normalisation
            double min = Double.MAX_VALUE, max = -Double.MAX_VALUE;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    min = Math.min(min, data[i][j]);
                    max = Math.max(max, data[i][j]);
                }
            }
            
            double range = max - min;
            if (range == 0) range = 1.0;
            
            // Dessiner la heatmap
            int pixelSize = PANEL_SIZE / n;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double normalized = (data[i][j] - min) / range;
                    Color color = getHeatmapColor(normalized);
                    g2d.setColor(color);
                    g2d.fillRect(i * pixelSize, j * pixelSize, pixelSize, pixelSize);
                }
            }
            
            g2d.dispose();
        }
        
        private Color getHeatmapColor(double value) {
            // Colormap bleu -> cyan -> vert -> jaune -> rouge
            if (value < 0.25) {
                float ratio = (float)(value / 0.25);
                return new Color(0, (int)(255 * ratio), 255);
            } else if (value < 0.5) {
                float ratio = (float)((value - 0.25) / 0.25);
                return new Color(0, 255, (int)(255 * (1 - ratio)));
            } else if (value < 0.75) {
                float ratio = (float)((value - 0.5) / 0.25);
                return new Color((int)(255 * ratio), 255, 0);
            } else {
                float ratio = (float)((value - 0.75) / 0.25);
                return new Color(255, (int)(255 * (1 - ratio)), 0);
            }
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            if (image != null) {
                g.drawImage(image, 0, 0, null);
            }
            
            // Titre
            g.setColor(Color.WHITE);
            g.setFont(new Font("Arial", Font.BOLD, 14));
            g.drawString(title, 10, 20);
            
            // Affichage de la valeur sous le curseur
            if (mouseX >= 0 && mouseY >= 0) {
                String valueStr = String.format("Valeur: %.4f", mouseValue);
                FontMetrics fm = g.getFontMetrics();
                int strWidth = fm.stringWidth(valueStr);
                int x = mouseX + 10;
                if (x + strWidth > getWidth()) x = mouseX - strWidth - 10;
                
                g.setColor(Color.BLACK);
                g.fillRect(x-2, mouseY-17, strWidth+4, 20);
                g.setColor(Color.WHITE);
                g.drawString(valueStr, x, mouseY);
            }
        }
        
        @Override
        public void mouseMoved(MouseEvent e) {
            mouseX = e.getX();
            mouseY = e.getY();
            
            int n = data.length;
            int pixelSize = PANEL_SIZE / n;
            int i = mouseX / pixelSize;
            int j = mouseY / pixelSize;
            
            if (i >= 0 && i < n && j >= 0 && j < n) {
                mouseValue = data[i][j];
                repaint();
            }
        }
        
        @Override
        public void mouseDragged(MouseEvent e) {
            mouseMoved(e);
        }
    }
    
    // Panneau pour les graphiques d'erreur
    static class ErrorPanel extends JPanel {
        private List<Integer> meshSizes;
        private List<Double> errors;
        private String title;
        
        public ErrorPanel(List<Integer> meshSizes, List<Double> errors, String title) {
            this.meshSizes = meshSizes;
            this.errors = errors;
            this.title = title;
            setPreferredSize(new Dimension(400, 300));
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int margin = 50;
            int width = getWidth() - 2 * margin;
            int height = getHeight() - 2 * margin;
            
            // Axes
            g2d.setColor(Color.BLACK);
            g2d.drawLine(margin, getHeight() - margin, getWidth() - margin, getHeight() - margin);
            g2d.drawLine(margin, margin, margin, getHeight() - margin);
            
            // Titre
            g2d.setFont(new Font("Arial", Font.BOLD, 14));
            FontMetrics fm = g2d.getFontMetrics();
            int titleWidth = fm.stringWidth(title);
            g2d.drawString(title, (getWidth() - titleWidth) / 2, 20);
            
            if (meshSizes.size() < 2) return;
            
            // Échelles logarithmiques
            double minMesh = Math.log10(meshSizes.get(0));
            double maxMesh = Math.log10(meshSizes.get(meshSizes.size() - 1));
            double minError = Math.log10(errors.stream().mapToDouble(Double::doubleValue).min().orElse(1e-10));
            double maxError = Math.log10(errors.stream().mapToDouble(Double::doubleValue).max().orElse(1.0));
            
            // Points et lignes
            g2d.setColor(Color.BLUE);
            g2d.setStroke(new BasicStroke(2));
            
            for (int i = 0; i < meshSizes.size() - 1; i++) {
                double x1 = margin + width * (Math.log10(meshSizes.get(i)) - minMesh) / (maxMesh - minMesh);
                double y1 = getHeight() - margin - height * (Math.log10(errors.get(i)) - minError) / (maxError - minError);
                double x2 = margin + width * (Math.log10(meshSizes.get(i+1)) - minMesh) / (maxMesh - minMesh);
                double y2 = getHeight() - margin - height * (Math.log10(errors.get(i+1)) - minError) / (maxError - minError);
                
                g2d.drawLine((int)x1, (int)y1, (int)x2, (int)y2);
                g2d.fillOval((int)x1-3, (int)y1-3, 6, 6);
            }
            
            // Dernier point
            double xLast = margin + width * (Math.log10(meshSizes.get(meshSizes.size()-1)) - minMesh) / (maxMesh - minMesh);
            double yLast = getHeight() - margin - height * (Math.log10(errors.get(errors.size()-1)) - minError) / (maxError - minError);
            g2d.fillOval((int)xLast-3, (int)yLast-3, 6, 6);
            
            // Labels des axes
            g2d.setColor(Color.BLACK);
            g2d.setFont(new Font("Arial", Font.PLAIN, 10));
            g2d.drawString("Nombre de mailles", getWidth()/2 - 40, getHeight() - 10);
            
            // Rotation pour le label y
            g2d.rotate(-Math.PI/2);
            g2d.drawString("Erreur L2", -getHeight()/2 - 20, 15);
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            // Création des solveurs
            ODESolver[] solvers = {
                new PoissonHelmholtz(50),
                new ConvectionDiffusion(50),
                new Poisson(50)
            };
            
            // Calcul des solutions et erreurs pour différentes tailles de mailles
            int[] meshSizes = {20, 50, 100, 200};
            
            for (ODESolver solver : solvers) {
                System.out.println("\n=== " + solver.getName() + " ===");
                
                List<Integer> sizes = new ArrayList<>();
                List<Double> errors = new ArrayList<>();
                
                // Test de convergence
                for (int n : meshSizes) {
                    ODESolver testSolver;
                    if (solver instanceof PoissonHelmholtz) {
                        testSolver = new PoissonHelmholtz(n);
                    } else if (solver instanceof ConvectionDiffusion) {
                        testSolver = new ConvectionDiffusion(n);
                    } else {
                        testSolver = new Poisson(n);
                    }
                    
                    testSolver.solve();
                    double error = testSolver.computeError();
                    
                    sizes.add(n);
                    errors.add(error);
                    
                    System.out.printf("Mailles: %d, Erreur L2: %.6e\n", n, error);
                }
                
                // Calcul de l'ordre de convergence
                if (errors.size() >= 2) {
                    double order = Math.log(errors.get(errors.size()-1) / errors.get(errors.size()-2)) / 
                                  Math.log((double)sizes.get(sizes.size()-2) / sizes.get(sizes.size()-1));
                    System.out.printf("Ordre de convergence estimé: %.2f\n", order);
                }
                
                // Résolution pour n=50 pour visualisation
                solver.solve();
                
                // Création de l'interface graphique
                JFrame frame = new JFrame(solver.getName());
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                frame.setLayout(new GridLayout(2, 2));
                
                // Panneaux de visualisation
                frame.add(new VisualizationPanel(solver.getSolution(), "Solution numérique"));
                frame.add(new VisualizationPanel(solver.getExactSolution(), "Solution exacte"));
                
                // Calcul de l'erreur point par point
                double[][] errorMap = new double[solver.getN()][solver.getN()];
                for (int i = 0; i < solver.getN(); i++) {
                    for (int j = 0; j < solver.getN(); j++) {
                        errorMap[i][j] = Math.abs(solver.getSolution()[i][j] - solver.getExactSolution()[i][j]);
                    }
                }
                frame.add(new VisualizationPanel(errorMap, "Erreur absolue"));
                
                // Graphique d'évolution de l'erreur
                frame.add(new ErrorPanel(sizes, errors, "Évolution de l'erreur"));
                
                frame.pack();
                frame.setVisible(true);
                frame.setLocationRelativeTo(null);
            }
        });
    }
}
package FiniteDifference2;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class ODEFiniteDifference {
    private static final double PI = Math.PI;
    
    // Classe pour la méthode directe
    static class DirectMethod {
        public static double[][] solve(int n, double h) {
            double[][] u = new double[n+1][n+1];
            double[][] f = new double[n+1][n+1];
            
            // Calcul de f(x,y) = sin(π*x) + y³
            for (int i = 0; i <= n; i++) {
                for (int j = 0; j <= n; j++) {
                    double x = i * h;
                    double y = j * h;
                    f[i][j] = Math.sin(PI * x) + Math.pow(y, 3);
                }
            }
            
            // Résolution par méthode de Gauss (élimination)
            int size = (n-1) * (n-1);
            double[][] A = new double[size][size];
            double[] b = new double[size];
            
            // Construction de la matrice A et du vecteur b
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    int idx = (i-1) * (n-1) + (j-1);
                    A[idx][idx] = -4.0 / (h*h);
                    b[idx] = -f[i][j];
                    
                    if (i > 1) A[idx][idx - (n-1)] = 1.0 / (h*h);
                    if (i < n-1) A[idx][idx + (n-1)] = 1.0 / (h*h);
                    if (j > 1) A[idx][idx - 1] = 1.0 / (h*h);
                    if (j < n-1) A[idx][idx + 1] = 1.0 / (h*h);
                }
            }
            
            // Résolution du système Ax = b
            double[] x = gaussElimination(A, b);
            
            // Reconstruction de la solution
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    int idx = (i-1) * (n-1) + (j-1);
                    u[i][j] = x[idx];
                }
            }
            
            return u;
        }
        
        private static double[] gaussElimination(double[][] A, double[] b) {
            int n = A.length;
            double[] x = new double[n];
            
            // Élimination
            for (int k = 0; k < n-1; k++) {
                for (int i = k+1; i < n; i++) {
                    if (Math.abs(A[k][k]) > 1e-10) {
                        double factor = A[i][k] / A[k][k];
                        for (int j = k; j < n; j++) {
                            A[i][j] -= factor * A[k][j];
                        }
                        b[i] -= factor * b[k];
                    }
                }
            }
            
            // Substitution arrière
            for (int i = n-1; i >= 0; i--) {
                x[i] = b[i];
                for (int j = i+1; j < n; j++) {
                    x[i] -= A[i][j] * x[j];
                }
                if (Math.abs(A[i][i]) > 1e-10) {
                    x[i] /= A[i][i];
                }
            }
            
            return x;
        }
    }
    
    // Classe pour Gauss-Seidel avec relaxation
    static class GaussSeidelRelaxation {
        public static double[][] solve(int n, double h, double omega, int maxIter, double tol) {
            double[][] u = new double[n+1][n+1];
            double[][] f = new double[n+1][n+1];
            
            // Calcul de f(x,y)
            for (int i = 0; i <= n; i++) {
                for (int j = 0; j <= n; j++) {
                    double x = i * h;
                    double y = j * h;
                    f[i][j] = Math.sin(PI * x) + Math.pow(y, 3);
                }
            }
            
            // Itérations Gauss-Seidel avec relaxation
            for (int iter = 0; iter < maxIter; iter++) {
                double maxChange = 0.0;
                
                for (int i = 1; i < n; i++) {
                    for (int j = 1; j < n; j++) {
                        double oldU = u[i][j];
                        double newU = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] + h*h*f[i][j]) / 4.0;
                        u[i][j] = (1.0 - omega) * oldU + omega * newU;
                        maxChange = Math.max(maxChange, Math.abs(u[i][j] - oldU));
                    }
                }
                
                if (maxChange < tol) break;
            }
            
            return u;
        }
    }
    
    // Classe pour Gauss-Seidel parallélisé
    static class GaussSeidelParallel {
        public static double[][] solve(int n, double h, double omega, int maxIter, double tol) {
            double[][] u = new double[n+1][n+1];
            double[][] f = new double[n+1][n+1];
            
            // Calcul de f(x,y)
            for (int i = 0; i <= n; i++) {
                for (int j = 0; j <= n; j++) {
                    double x = i * h;
                    double y = j * h;
                    f[i][j] = Math.sin(PI * x) + Math.pow(y, 3);
                }
            }
            
            ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
            
            // Itérations parallèles
            for (int iter = 0; iter < maxIter; iter++) {
                final int currentIter = iter;
                List<Future<Double>> futures = new ArrayList<>();
                
                // Division en chunks pour parallélisation
                int numThreads = Runtime.getRuntime().availableProcessors();
                int chunkSize = Math.max(1, (n-1) / numThreads);
                
                for (int t = 0; t < numThreads; t++) {
                    final int startI = Math.max(1, 1 + t * chunkSize);
                    final int endI = Math.min(n, 1 + (t + 1) * chunkSize);
                    
                    futures.add(executor.submit(() -> {
                        double maxChange = 0.0;
                        for (int i = startI; i < endI; i++) {
                            for (int j = 1; j < n; j++) {
                                double oldU = u[i][j];
                                double newU = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] + h*h*f[i][j]) / 4.0;
                                u[i][j] = (1.0 - omega) * oldU + omega * newU;
                                maxChange = Math.max(maxChange, Math.abs(u[i][j] - oldU));
                            }
                        }
                        return maxChange;
                    }));
                }
                
                // Attendre tous les threads et calculer le changement maximum
                double maxChange = 0.0;
                try {
                    for (Future<Double> future : futures) {
                        maxChange = Math.max(maxChange, future.get());
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
                
                if (maxChange < tol) break;
            }
            
            executor.shutdown();
            return u;
        }
    }
    
    // Solution analytique approximative pour comparaison
    public static double[][] exactSolution(int n, double h) {
        double[][] exact = new double[n+1][n+1];
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n; j++) {
                double x = i * h;
                double y = j * h;
                // Solution approximative basée sur la forme de f
                exact[i][j] = Math.sin(PI * x) / (PI * PI) + y * y * y * y * y / 20.0;
            }
        }
        return exact;
    }
    
    // Calcul de l'erreur
    public static double calculateError(double[][] numerical, double[][] exact) {
        double error = 0.0;
        int n = numerical.length - 1;
        
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n; j++) {
                error += Math.pow(numerical[i][j] - exact[i][j], 2);
            }
        }
        
        return Math.sqrt(error / ((n+1) * (n+1)));
    }
    
    // Calcul de l'ordre de convergence
    public static double calculateConvergenceOrder(double error1, double error2, double h1, double h2) {
        return Math.log(error1 / error2) / Math.log(h1 / h2);
    }
    
    // Classe pour la visualisation
    static class Visualization extends JPanel implements MouseMotionListener {
        private double[][] data;
        private String title;
        private double minVal, maxVal;
        private int mouseX, mouseY;
        private boolean showMousePos = false;
        
        public Visualization(double[][] data, String title) {
            this.data = data;
            this.title = title;
            this.addMouseMotionListener(this);
            
            // Trouver min et max pour la normalisation des couleurs
            minVal = Double.MAX_VALUE;
            maxVal = Double.MIN_VALUE;
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data[i].length; j++) {
                    minVal = Math.min(minVal, data[i][j]);
                    maxVal = Math.max(maxVal, data[i][j]);
                }
            }
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int width = getWidth() - 100;
            int height = getHeight() - 100;
            int n = data.length - 1;
            
            double cellWidth = (double) width / n;
            double cellHeight = (double) height / n;
            
            // Dessiner la grille colorée
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double value = data[i][j];
                    Color color = getColorForValue(value);
                    g2d.setColor(color);
                    
                    int x = 50 + (int)(i * cellWidth);
                    int y = 50 + (int)((n-1-j) * cellHeight);
                    
                    g2d.fillRect(x, y, (int)cellWidth + 1, (int)cellHeight + 1);
                }
            }
            
            // Dessiner les axes
            g2d.setColor(Color.BLACK);
            g2d.drawRect(50, 50, width, height);
            
            // Titre
            g2d.setFont(new Font("Arial", Font.BOLD, 16));
            FontMetrics fm = g2d.getFontMetrics();
            int titleWidth = fm.stringWidth(title);
            g2d.drawString(title, (getWidth() - titleWidth) / 2, 30);
            
            // Échelle de couleurs
            drawColorScale(g2d, getWidth() - 40, 50, 20, height);
            
            // Afficher la valeur au pointeur
            if (showMousePos) {
                int gridX = (int)((mouseX - 50) / cellWidth);
                int gridY = n - 1 - (int)((mouseY - 50) / cellHeight);
                
                if (gridX >= 0 && gridX < n && gridY >= 0 && gridY < n) {
                    double value = data[gridX][gridY];
                    String valueText = String.format("(%.2f, %.2f): %.4f", 
                        gridX * (1.0/n), gridY * (1.0/n), value);
                    
                    g2d.setColor(Color.WHITE);
                    g2d.fillRect(mouseX + 5, mouseY - 20, 150, 20);
                    g2d.setColor(Color.BLACK);
                    g2d.drawString(valueText, mouseX + 8, mouseY - 5);
                }
            }
        }
        
        private Color getColorForValue(double value) {
            double normalized = (value - minVal) / (maxVal - minVal);
            
            // Palette de couleurs : bleu -> vert -> jaune -> rouge
            if (normalized < 0.25) {
                float ratio = (float)(normalized / 0.25);
                return new Color(0, (int)(255 * ratio), 255);
            } else if (normalized < 0.5) {
                float ratio = (float)((normalized - 0.25) / 0.25);
                return new Color(0, 255, (int)(255 * (1 - ratio)));
            } else if (normalized < 0.75) {
                float ratio = (float)((normalized - 0.5) / 0.25);
                return new Color((int)(255 * ratio), 255, 0);
            } else {
                float ratio = (float)((normalized - 0.75) / 0.25);
                return new Color(255, (int)(255 * (1 - ratio)), 0);
            }
        }
        
        private void drawColorScale(Graphics2D g2d, int x, int y, int width, int height) {
            for (int i = 0; i < height; i++) {
                double value = minVal + (maxVal - minVal) * (1.0 - (double)i / height);
                Color color = getColorForValue(value);
                g2d.setColor(color);
                g2d.drawLine(x, y + i, x + width, y + i);
            }
            
            g2d.setColor(Color.BLACK);
            g2d.drawRect(x, y, width, height);
            
            // Labels
            g2d.setFont(new Font("Arial", Font.PLAIN, 10));
            g2d.drawString(String.format("%.3f", maxVal), x + width + 5, y + 10);
            g2d.drawString(String.format("%.3f", minVal), x + width + 5, y + height);
        }
        
        @Override
        public void mouseMoved(MouseEvent e) {
            mouseX = e.getX();
            mouseY = e.getY();
            showMousePos = true;
            repaint();
        }
        
        @Override
        public void mouseDragged(MouseEvent e) {
            mouseMoved(e);
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            // Paramètres
            int[] meshSizes = {10, 20, 40/*, 80, 160, 320, 640*/};
            double omega = 1.2; // Facteur de relaxation
            int maxIter = 10000;
            double tol = 1e-8;
            
            List<Double> errors = new ArrayList<>();
            List<Double> convergenceOrders = new ArrayList<>();
            
            System.out.println("Résolution de -u'' = sin(πx) + y³");
            System.out.println("=====================================");
            
            JFrame mainFrame = new JFrame("Résolution d'Équation Différentielle 2D");
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new GridLayout(2, 2));
            
            double previousError = 0;
            double previousH = 0;
            
            for (int meshSize : meshSizes) {
                double h = 1.0 / meshSize;
                
                System.out.printf("\nMaillage %dx%d (h = %.6f):\n", meshSize, meshSize, h);
                
                // Résolution avec les trois méthodes
                long startTime = System.currentTimeMillis();
                double[][] directSol = DirectMethod.solve(meshSize, h);
                long directTime = System.currentTimeMillis() - startTime;
                
                startTime = System.currentTimeMillis();
                double[][] gsRelaxSol = GaussSeidelRelaxation.solve(meshSize, h, omega, maxIter, tol);
                long gsRelaxTime = System.currentTimeMillis() - startTime;
                
                startTime = System.currentTimeMillis();
                double[][] gsParallelSol = GaussSeidelParallel.solve(meshSize, h, omega, maxIter, tol);
                long gsParallelTime = System.currentTimeMillis() - startTime;
                
                // Solution exacte
                double[][] exactSol = exactSolution(meshSize, h);
                
                // Calcul des erreurs
                double errorDirect = calculateError(directSol, exactSol);
                double errorGSRelax = calculateError(gsRelaxSol, exactSol);
                double errorGSParallel = calculateError(gsParallelSol, exactSol);
                
                System.out.printf("Méthode directe: Erreur = %.6e, Temps = %d ms\n", errorDirect, directTime);
                System.out.printf("Gauss-Seidel + Relaxation: Erreur = %.6e, Temps = %d ms\n", errorGSRelax, gsRelaxTime);
                System.out.printf("Gauss-Seidel Parallèle: Erreur = %.6e, Temps = %d ms\n", errorGSParallel, gsParallelTime);
                
                errors.add(errorDirect);
                
                // Calcul de l'ordre de convergence
                if (previousError > 0) {
                    double order = calculateConvergenceOrder(previousError, errorDirect, previousH, h);
                    convergenceOrders.add(order);
                    System.out.printf("Ordre de convergence: %.2f\n", order);
                }
                
                previousError = errorDirect;
                previousH = h;
                
                // Visualisation pour certaines tailles de maillage
                if (meshSize == 40) {
                    mainFrame.add(new Visualization(exactSol, "Solution Exacte (40x40)"));
                    mainFrame.add(new Visualization(directSol, "Méthode Directe (40x40)"));
                    mainFrame.add(new Visualization(gsRelaxSol, "Gauss-Seidel Relaxation (40x40)"));
                    mainFrame.add(new Visualization(gsParallelSol, "Gauss-Seidel Parallèle (40x40)"));
                }
            }
            
            // Affichage du graphique de convergence
            System.out.println("\nOrdres de convergence moyens:");
            double avgOrder = convergenceOrders.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
            System.out.printf("Ordre moyen: %.2f\n", avgOrder);
            
            // Fenêtre principale
            mainFrame.setSize(1200, 800);
            mainFrame.setLocationRelativeTo(null);
            mainFrame.setVisible(true);
            
            // Graphique de l'évolution de l'erreur
            JFrame errorFrame = new JFrame("Évolution de l'Erreur");
            errorFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            errorFrame.add(new ErrorPlot(meshSizes, errors));
            errorFrame.setSize(600, 400);
            errorFrame.setLocationRelativeTo(null);
            errorFrame.setVisible(true);
        });
    }
    
    // Classe pour tracer l'évolution de l'erreur
    static class ErrorPlot extends JPanel {
        private int[] meshSizes;
        private List<Double> errors;
        
        public ErrorPlot(int[] meshSizes, List<Double> errors) {
            this.meshSizes = meshSizes;
            this.errors = errors;
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int width = getWidth() - 100;
            int height = getHeight() - 100;
            
            // Axes
            g2d.setColor(Color.BLACK);
            g2d.drawLine(50, height + 50, width + 50, height + 50); // X
            g2d.drawLine(50, 50, 50, height + 50); // Y
            
            // Titre
            g2d.setFont(new Font("Arial", Font.BOLD, 16));
            FontMetrics fm = g2d.getFontMetrics();
            String title = "Évolution de l'Erreur en fonction du Nombre de Mailles";
            int titleWidth = fm.stringWidth(title);
            g2d.drawString(title, (getWidth() - titleWidth) / 2, 30);
            
            // Labels des axes
            g2d.setFont(new Font("Arial", Font.PLAIN, 12));
            g2d.drawString("Nombre de mailles", width/2, height + 80);
            
            g2d.rotate(-Math.PI/2);
            g2d.drawString("Erreur (log)", -height/2 - 50, 20);
            g2d.rotate(Math.PI/2);
            
            if (errors.isEmpty()) return;
            
            // Échelle logarithmique
            double minError = errors.stream().mapToDouble(Double::doubleValue).min().orElse(1e-10);
            double maxError = errors.stream().mapToDouble(Double::doubleValue).max().orElse(1.0);
            
            // Tracer les points et les lignes
            g2d.setColor(Color.BLUE);
            g2d.setStroke(new BasicStroke(2));
            
            for (int i = 0; i < Math.min(meshSizes.length, errors.size()); i++) {
                double x = 50 + (double)width * i / (meshSizes.length - 1);
                double logError = Math.log10(errors.get(i));
                double logMin = Math.log10(minError);
                double logMax = Math.log10(maxError);
                double y = 50 + height - (double)height * (logError - logMin) / (logMax - logMin);
                
                // Point
                g2d.fillOval((int)x - 3, (int)y - 3, 6, 6);
                
                // Ligne vers le point suivant
                if (i < Math.min(meshSizes.length, errors.size()) - 1) {
                    double nextX = 50 + (double)width * (i + 1) / (meshSizes.length - 1);
                    double nextLogError = Math.log10(errors.get(i + 1));
                    double nextY = 50 + height - (double)height * (nextLogError - logMin) / (logMax - logMin);
                    g2d.drawLine((int)x, (int)y, (int)nextX, (int)nextY);
                }
                
                // Label
                g2d.setColor(Color.BLACK);
                g2d.drawString(String.valueOf(meshSizes[i]), (int)x - 10, height + 65);
                g2d.setColor(Color.BLUE);
            }
        }
    }
}
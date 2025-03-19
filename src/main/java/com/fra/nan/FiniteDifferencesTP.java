/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.fra.nan;

/**
 *
 * @author Roddier
 */
import java.awt.BorderLayout;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Programme complet pour résoudre -u'' = f sur (0,1) par différences finies,
 * calculer la vitesse de convergence en norme L∞ et afficher les graphiques.
 */
public class FiniteDifferencesTP {

    /**
     * Résout le système discret obtenu avec la méthode des différences finies.
     * Le système est: u_{i-1} - 2u_i + u_{i+1} = h^2 * u''(x_i) pour i = 1,...,N.
     *
     * @param N nombre de points intérieurs (maillage de taille N+2)
     * @param u0 condition limite en 0
     * @param u1 condition limite en 1
     * @param secondDeriv fonction donnant u''(x)
     * @return solution approximative aux points intérieurs (u1,...,u_N)
     */
    public static double[] solveFD(int N, double u0, double u1, Function secondDeriv) {
        double h = 1.0 / (N + 1);
        double[] x = new double[N];
        double[] b = new double[N];
        double[][] A = new double[N][N];

        // Construction du second membre : b[i] = h^2 * u''(x_{i+1})
        for (int i = 0; i < N; i++) {
            x[i] = (i + 1) * h;
            b[i] = secondDeriv.eval(x[i]) * h * h;
        }
        // Prise en compte des conditions limites
        b[0]    += u0;  // pour i = 1, la contribution de u(0)
        b[N - 1] += u1;  // pour i = N, la contribution de u(1)

        // Construction de la matrice tridiagonale
        for (int i = 0; i < N; i++) {
            A[i][i] = -2;
            if (i > 0) {
                A[i][i - 1] = 1;
            }
            if (i < N - 1) {
                A[i][i + 1] = 1;
            }
        }
        // Résolution par élimination de Gauss
        return gaussElimination(A, b);
    }

    /**
     * Méthode d'élimination de Gauss pour résoudre Ax = b.
     */
    public static double[] gaussElimination(double[][] A, double[] b) {
        int N = b.length;
        // Phase d'élimination
        for (int i = 0; i < N; i++) {
            for (int k = i + 1; k < N; k++) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < N; j++) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
        // Phase de substitution arrière
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            x[i] = b[i];
            for (int j = i + 1; j < N; j++) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
        }
        return x;
    }

    /**
     * Calcule et affiche l'erreur en norme L∞ pour différents maillages
     * et estime l'ordre de convergence.
     *
     * @param meshSizes tableau des tailles de maillage (nombre de points intérieurs)
     * @param uExact fonction de la solution exacte u(x)
     * @param uSecond fonction donnant u''(x)
     */
    public static void computeConvergence(int[] meshSizes, Function uExact, Function uSecond) {
        double prevError = -1;
        System.out.println("Convergence pour u(x) = " + uExact.getDescription());
        for (int N : meshSizes) {
            double h = 1.0 / (N + 1);
            double[] x = new double[N];
            for (int i = 0; i < N; i++) {
                x[i] = (i + 1) * h;
            }
            double[] uNum = solveFD(N, uExact.eval(0), uExact.eval(1), uSecond);
            double maxError = 0;
            for (int i = 0; i < N; i++) {
                double err = Math.abs(uNum[i] - uExact.eval(x[i]));
                if (err > maxError) {
                    maxError = err;
                }
            }
            System.out.printf("N = %d, h = %.5e, erreur = %.5e", N, h, maxError);
            if (prevError > 0) {
                double order = Math.log(prevError / maxError) / Math.log(2);
                System.out.printf("  --> ordre ~ %.2f", order);
            }
            System.out.println();
            prevError = maxError;
        }
        System.out.println();
    }

    public static void main(String[] args) {
        // Maillages à tester
        int[] meshSizes = {10, 20, 40, 80, 160, 320};

        // Cas 1 : u(x) = sin(pi*x)
        // u(0) = 0, u(1) = sin(pi) = 0, et u''(x) = -pi^2 sin(pi*x)
        Function uExactSin = new Function() {
            public double eval(double x) {
                return Math.sin(Math.PI * x);
            }
            public String getDescription() {
                return "sin(pi*x)";
            }
        };
        Function uSecondSin = new Function() {
            public double eval(double x) {
                return -Math.PI * Math.PI * Math.sin(Math.PI * x);
            }
            public String getDescription() {
                return "u''(x) de sin(pi*x)";
            }
        };

        // Cas 2 : u(x) = x^3
        // u(0) = 0, u(1) = 1, et u''(x) = 6x
        Function uExactCubic = new Function() {
            public double eval(double x) {
                return x * x * x;
            }
            public String getDescription() {
                return "x^3";
            }
        };
        Function uSecondCubic = new Function() {
            public double eval(double x) {
                return 6 * x;
            }
            public String getDescription() {
                return "u''(x) de x^3";
            }
        };

        // Calcul de l'ordre de convergence pour les deux cas
        computeConvergence(meshSizes, uExactSin, uSecondSin);
        computeConvergence(meshSizes, uExactCubic, uSecondCubic);

        // Représentation graphique (pour le cas sin(pi*x) sur le maillage le plus fin)
        int N = 320;
        double h = 1.0 / (N + 1);
        int Ntot = N + 2; // incluant les bornes
        double[] x = new double[Ntot];
        double[] uExactVals = new double[Ntot];
        double[] uNumVals = new double[Ntot];

        // Calcul des points (avec conditions aux limites)
        for (int i = 0; i < Ntot; i++) {
            x[i] = i * h;
            uExactVals[i] = uExactSin.eval(x[i]);
        }
        // Solution numérique (pour les points intérieurs)
        double[] uInterior = solveFD(N, uExactSin.eval(0), uExactSin.eval(1), uSecondSin);
        uNumVals[0] = uExactSin.eval(0);
        uNumVals[Ntot - 1] = uExactSin.eval(1);
        for (int i = 0; i < N; i++) {
            uNumVals[i + 1] = uInterior[i];
        }
        // Affichage des courbes (solution et erreur)
        Plotter.plotTwoSeries(x, uExactVals, uNumVals, "Solution exacte et approchée pour sin(pi*x)",
                "x", "u(x)");
        Plotter.plotSingleSeries(x, uNumVals, "Solution approchée pour sin(pi*x)",
                "x", "u(x)");
        Plotter.plotSingleSeries(x, uExactVals, "Solution exacte pour sin(pi*x)",
                "x", "u(x)");
        // Courbe d'erreur
        double[] error = new double[Ntot];
        for (int i = 0; i < Ntot; i++) {
            error[i] = Math.abs(uExactVals[i] - uNumVals[i]);
        }
        Plotter.plotSingleSeries(x, error, "Erreur |u_exact - u_num|", "x", "Erreur");
    }
}

/**
 * Interface pour définir une fonction réelle.
 */
interface Function {
    double eval(double x);
    default String getDescription() {
        return "";
    }
}

/**
 * Classe utilitaire pour afficher des graphiques à l'aide de JFreeChart.
 */
class Plotter {
    /**
     * Affiche deux séries (par exemple, la solution exacte et l'approchée) dans un même graphique.
     */
    public static void plotTwoSeries(double[] x, double[] series1, double[] series2,
                                     String title, String xLabel, String yLabel) {
        XYSeries s1 = new XYSeries("Exacte");
        XYSeries s2 = new XYSeries("Approchée");
        for (int i = 0; i < x.length; i++) {
            s1.add(x[i], series1[i]);
            s2.add(x[i], series2[i]);
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(s1);
        dataset.addSeries(s2);
        JFreeChart chart = ChartFactory.createXYLineChart(title, xLabel, yLabel,
                dataset, PlotOrientation.VERTICAL, true, true, false);
        showChart(chart, title);
    }

    /**
     * Affiche une seule série (par exemple, la courbe d'erreur).
     */
    public static void plotSingleSeries(double[] x, double[] series, String title,
                                        String xLabel, String yLabel) {
        XYSeries s = new XYSeries(title);
        for (int i = 0; i < x.length; i++) {
            s.add(x[i], series[i]);
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(s);
        JFreeChart chart = ChartFactory.createXYLineChart(title, xLabel, yLabel,
                dataset, PlotOrientation.VERTICAL, true, true, false);
        showChart(chart, title);
    }

    /**
     * Affiche le graphique dans une fenêtre.
     */
    private static void showChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        ChartPanel chartPanel = new ChartPanel(chart);
        frame.getContentPane().add(chartPanel, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);
    }
}

